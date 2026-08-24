/*******************************************************************************
Copyright (c) 2015, Jonathan Hiller
To cite academic use of Voxelyze: Jonathan Hiller and Hod Lipson "Dynamic Simulation of Soft Multimaterial 3D-Printed Objects" Soft Robotics. March 2014, 1(1): 88-101.
Available at http://online.liebertpub.com/doi/pdfplus/10.1089/soro.2013.0010

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

/*******************************************************************************
Copyright (c) 2026, YoonSik Shim (Based on Voxelyze by Jonathan Hiller)

[Modification Notice]
- Newly created file to implement nested parallelism architecture (doTimeStep_Nested, updateCollisions_Nested).
- Added Windows Processor Group Thread Binding (PinToGroupIfNeeded / SetThreadGroupAffinity) to optimize performance for large-scale multi-core systems (e.g., AMD Threadripper 3990X).
- Resolved deadlocks in collision processing using copyprivate and team-wide execution blocks.
*******************************************************************************/


#ifdef _WIN32
// GetThreadGroupAffinity / SetThreadGroupAffinity 프로토타입을 보려면
// 프로젝트에서 _WIN32_WINNT >= 0x0601(Windows7)로 잡혀 있어야 합니다.
// (MSVC 설정/전처리 정의에 넣거나, include 전에 define)
#ifndef _WIN32_WINNT
	#define _WIN32_WINNT 0x0601
#endif

#include <windows.h>
#include <processtopologyapi.h>

static inline KAFFINITY BuildFullMaskForGroup(WORD group)
{
	DWORD n = GetActiveProcessorCount(group); // group 내 활성 논리프로세서 수(<=64)
	if (n >= 64) return ~KAFFINITY(0);
	if (n == 0)  return 0;
	return (KAFFINITY(1) << n) - 1;
}

static inline WORD GetMyPrimaryGroup()
{
	GROUP_AFFINITY ga{};
	GetThreadGroupAffinity(GetCurrentThread(), &ga);
	return ga.Group; // 현재 스레드의 primary group (Win11에서도 primary group 기준) :contentReference[oaicite:0]{index=0}
}

static inline void PinToGroupIfNeeded(WORD group)
{
	// 같은 OS 스레드에서 doTimeStep이 수천 번 호출되므로,
	// 매번 SetThreadGroupAffinity를 때리지 않도록 thread_local 캐시.
	static thread_local int pinned_group = -1;
	if (pinned_group == (int)group) return;

	GROUP_AFFINITY ga{};
	ga.Group = group;
	ga.Mask  = BuildFullMaskForGroup(group);

	if (SetThreadGroupAffinity(GetCurrentThread(), &ga, nullptr))
	{
		pinned_group = (int)group;
	}
	else
	{
		// 디버깅 필요 시:
		// printf("SetThreadGroupAffinity failed err=%lu\n", GetLastError());
	}
}
#endif



#include <omp.h>
#include "Voxelyze.h"


#define VER_02_19

#ifdef VER_02_19


bool CVoxelyze::doTimeStep_Nested(double dt)
{
    if (dt < 0) dt = recommendedTimeStep();

    int linkCount = (int)linksList.size();
    int voxCount  = (int)voxelsList.size();

    // 내부 병렬을 쓸지 말지 (num_thread=4, is_thread=true일 때만)
    const bool use_inner_mt = (is_thread && num_thread > 1);

    // "부모(outer) 스레드가 속한 그룹을 읽어둔다.
    // outer에서 이미 SetThreadGroupAffinity로 group0/1로 박아놨다면 여기 값이 그대로 나옵니다.
    WORD parent_group = 0;
#ifdef _WIN32
    if (use_inner_mt) parent_group = GetMyPrimaryGroup();
#endif

    // OpenMP 2.x에서 bool atomic write 같은 걸 피하려고 int 플래그 사용
    int diverged_flag = 0;

    #pragma omp parallel num_threads(num_thread) if(use_inner_mt) firstprivate(parent_group) shared(diverged_flag, linkCount, voxCount)
    {
        // 여기(병렬영역 시작 직후)가“inner 스레드 affinity 설정 위치입니다.
#ifdef _WIN32
        if (use_inner_mt) PinToGroupIfNeeded(parent_group);
#endif

        // [Phase 1] 링크 기하학 정보 갱신
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)
            linksList[i]->preUpdateGeometry();

        // [Phase 2] 복셀 캐시 갱신
        #pragma omp for schedule(static)
        for (int i = 0; i < voxCount; i++)
            voxelsList[i]->poissonsStrain();

        // [Phase 3] 힘 계산 + 발산 체크
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)
        {
            linksList[i]->finalUpdateForces();

            if (linksList[i]->axialStrain() > 100) {
                // OpenMP 2.x에서도 안전하게: atomic update
                #pragma omp atomic
                diverged_flag |= 1;
            }
        }


        // updateCollisions_Nested()는 single이 아니라 "팀 전체"가 호출해야
        //    updateCollisions_Nested() 내부의 omp for가 정상 동작(교착 방지)
        int do_coll = 0;
        #pragma omp single copyprivate(do_coll)
        {
            do_coll = (!diverged_flag && collisions) ? 1 : 0;
        }

        if (do_coll)
        {
            // 팀 전체 호출
            updateCollisions_Nested();
        }



        // [Phase 4] 복셀 timestep (발산이면 모두 스킵)
        if (!diverged_flag)
        {
            // 주의: 이 omp for는 “팀 전체가 동일 조건으로” 만나야 합니다.
            // diverged_flag는 위 omp for 끝의 implicit barrier 이후 확정되므로,
            // 여기서 모든 스레드가 같은 값(0/1)을 보게 됩니다.
            #pragma omp for schedule(static)
            for (int i = 0; i < voxCount; i++)
                voxelsList[i]->timeStep(dt);
        }
    } // end parallel

    if (diverged_flag) return false;

    currentTime += dt;
    return true;
}




void CVoxelyze::updateCollisions_Nested()
{
    const double watchRadiusVx  = 2.0 * boundingRadius + watchDistance; // outer radius to track all voxels within
    const double watchRadiusMm  = (double)(voxSize * watchRadiusVx);    // outer radius to track all voxels within (mm)
    const double watchRadiusMm2 = watchRadiusMm * watchRadiusMm;
    const double recalcDist     = (double)(voxSize * watchDistance * 0.5);
    const double recalcDist2    = recalcDist * recalcDist;

    const int voxCount = (int)voxelsList.size();

    // ------------------------------------------------------------
    // 안전장치: 만약 parallel 밖에서 호출되는 경로가 있으면 직렬로 처리
    // (doTimeStep_Nested 내부에서 호출될 때는 omp_in_parallel()==1)
    // ------------------------------------------------------------
    if (!omp_in_parallel())
    {
        // ---- 원래 직렬 로직 ----
        if (nearbyStale) {
            for (int i = 0; i < voxCount; i++)
                voxelsList[i]->generateNearby(watchRadiusVx * 2.0, false);

            nearbyStale = false;
            collisionsStale = true;
        }

        for (int i = 0; i < voxCount; i++)
        {
            CVX_Voxel* pV = voxelsList[i];
            if (pV->isSurface() && (pV->pos - *pV->lastColWatchPosition).Length2() > recalcDist2)
            {
                collisionsStale = true;
            }
        }

        if (collisionsStale) regenerateCollisions(watchRadiusMm2);

        const int colCount = (int)collisionsList.size();
        for (int i = 0; i < colCount; i++)
            collisionsList[i]->updateContactForce();

        return;
    }

    // ------------------------------------------------------------
    // 여기부터는: doTimeStep_Nested()의 inner parallel 팀 내부에서 실행
    // => 새 parallel 생성 금지, omp for / single만 사용
    // ------------------------------------------------------------

    // 1) nearbyStale 처리: 자료구조/상태 변경 가능 -> single
    #pragma omp single
    {
        if (nearbyStale)
        {
            for (int i = 0; i < voxCount; i++)
                voxelsList[i]->generateNearby(watchRadiusVx * 2.0, false);

            nearbyStale = false;
            collisionsStale = true;
        }
    }
    // single 끝 implicit barrier: 팀 동기화 완료

    // 2) stale 체크 루프: 병렬화(검증된 구간)
    int local_hit = 0;

    #pragma omp for schedule(static)
    for (int i = 0; i < voxCount; i++)
    {
        CVX_Voxel* pV = voxelsList[i];
        if (pV->isSurface() &&
            (pV->pos - *pV->lastColWatchPosition).Length2() > recalcDist2)
        {
            local_hit = 1;
        }
    }

    // collisionsStale은 shared 멤버이므로 동시 write 레이스 방지
    if (local_hit)
    {
        #pragma omp critical(COLLISIONS_STALE_SET)
        {
            collisionsStale = true;
        }
    }

    // regenerateCollisions 전에 모든 스레드의 stale 반영이 끝나야 함
    #pragma omp barrier

    // 3) regenerateCollisions: collisionsList 재구성 가능 -> single
    int colCount = 0;
    #pragma omp single copyprivate(colCount)
    {
        if (collisionsStale)
            regenerateCollisions(watchRadiusMm2);

        // colCount는 모든 스레드가 동일 값 사용해야 하므로 copyprivate로 브로드캐스트
        colCount = (int)collisionsList.size();
    }
    // single 끝 implicit barrier: collisionsList 안정화 + colCount 공유 완료

    // 4) contact force 업데이트 루프: 병렬화(검증된 구간 + updateContactForce 안전)
    #pragma omp for schedule(static)
    for (int i = 0; i < colCount; i++)
    {
        collisionsList[i]->updateContactForce();
    }
    // omp for implicit barrier: 다음 단계에서 force를 읽어도 안전
}



#else



bool CVoxelyze::doTimeStep_Nested(double dt)
{
    if (dt < 0) dt = recommendedTimeStep();

    int linkCount = (int)linksList.size();
    int voxCount  = (int)voxelsList.size();

    // 내부 병렬을 쓸지 말지 (num_thread=4, is_thread=true일 때만)
    const bool use_inner_mt = (is_thread && num_thread > 1);

    // "부모(outer) 스레드가 속한 그룹을 읽어둔다.
    // outer에서 이미 SetThreadGroupAffinity로 group0/1로 박아놨다면 여기 값이 그대로 나옵니다.
    WORD parent_group = 0;
#ifdef _WIN32
    if (use_inner_mt) parent_group = GetMyPrimaryGroup();
#endif

    // OpenMP 2.x에서 bool atomic write 같은 걸 피하려고 int 플래그 사용
    int diverged_flag = 0;

    #pragma omp parallel num_threads(num_thread) if(use_inner_mt) firstprivate(parent_group) shared(diverged_flag, linkCount, voxCount)
    {
        // 여기(병렬영역 시작 직후)가“inner 스레드 affinity 설정 위치입니다.
#ifdef _WIN32
        if (use_inner_mt) PinToGroupIfNeeded(parent_group);
#endif

        // [Phase 1] 링크 기하학 정보 갱신
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)
            linksList[i]->preUpdateGeometry();

        // [Phase 2] 복셀 캐시 갱신
        #pragma omp for schedule(static)
        for (int i = 0; i < voxCount; i++)
            voxelsList[i]->poissonsStrain();

        // [Phase 3] 힘 계산 + 발산 체크
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)
        {
            linksList[i]->finalUpdateForces();

            if (linksList[i]->axialStrain() > 100) {
                // OpenMP 2.x에서도 안전하게: atomic update
                #pragma omp atomic
                diverged_flag |= 1;
            }
        }

        // updateCollisions는 1개 스레드만 실행
        #pragma omp single
        {
            if (!diverged_flag && collisions)
                updateCollisions_Nested();
        }

        // [Phase 4] 복셀 timestep (발산이면 모두 스킵)
        if (!diverged_flag)
        {
            // 주의: 이 omp for는 “팀 전체가 동일 조건으로” 만나야 합니다.
            // diverged_flag는 위 omp for 끝의 implicit barrier 이후 확정되므로,
            // 여기서 모든 스레드가 같은 값(0/1)을 보게 됩니다.
            #pragma omp for schedule(static)
            for (int i = 0; i < voxCount; i++)
                voxelsList[i]->timeStep(dt);
        }
    } // end parallel

    if (diverged_flag) return false;

    currentTime += dt;
    return true;
}


void CVoxelyze::updateCollisions_Nested()
{

	double watchRadiusVx = 2.0*boundingRadius + watchDistance; //outer radius to track all voxels within
	double watchRadiusMm = (double)(voxSize*watchRadiusVx); //outer radius to track all voxels within
	double recalcDist = (double)(voxSize*watchDistance*0.5); //if the voxel moves further than this radius, recalc! //1/2 the allowabl, accounting for 0.5x radius of the voxel iself

	//if voxels have been added/removed, regenerate everybody's nearby list
	if (nearbyStale){
		for (std::vector<CVX_Voxel*>::iterator it=voxelsList.begin(); it != voxelsList.end(); it++)
		{
			(*it)->generateNearby(watchRadiusVx*2.0, false);
		}
		nearbyStale = false;
		collisionsStale = true;
	}

	//check if any voxels have moved far enough to make collisions stale
	int voxCount = (int)voxelsList.size();

	
//OMP_PRAG
	for (int i=0; i<voxCount; i++)
	{
		CVX_Voxel* pV = voxelsList[i]; //(*it);
		if (pV->isSurface() && (pV->pos - *pV->lastColWatchPosition).Length2() > recalcDist*recalcDist)
		{
			collisionsStale = true;	
		}
	}



	if (collisionsStale) regenerateCollisions(watchRadiusMm*watchRadiusMm);

	//update the forces!

	int colCount = (int)collisionsList.size();


//OMP_PRAG
	for (int i=0; i<colCount; i++)
	{
		collisionsList[i]->updateContactForce();
	}

}


#endif


