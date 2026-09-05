// ============================================================================
//  [PATCH 2026-09-05] inner thread >= 2 에서만 나타나는 메모리 증가 대응
// ----------------------------------------------------------------------------
//  [배경]
//   * 관측: inner thread 를 모두 1 로 두면 메모리 증가 없음.
//           516복셀 로봇에 2개 이상 주면 빠르게 증가.
//           2복셀 로봇에만 2개 주면 아주 천천히 증가.
//   * 소거됨: vcomp 중첩 fork 자체는 무죄 (fork/barrier 만 반복하는 최소
//             재현 프로그램에서 증가량 0). 배리어 구현도 무죄.
//   * 소거됨: CVX_Collision::updateContactForce() 는 this->force 만 쓰고
//             복셀은 읽기만 하므로 omp for 로 병렬화해도 레이스 없음.
//             (pull 모델: CVX_Voxel::force() 가 contactForce() 로 당겨 읽음)
//
//  [남은 원인 - 이 패치가 겨냥하는 것]
//   #pragma omp single 은 "먼저 도착한 스레드"가 실행한다. 즉 스텝마다 실행자가
//   바뀐다. 그런데 single 안에서 regenerateCollisions() / generateNearby() 가
//   힙을 대량으로 건드린다:
//       clearCollisions();                                 // N x delete
//       collisionsList.push_back(new CVX_Collision(...));  // N x new
//       pV1->colWatch->push_back(...);                     // 2N x push_back
//   결과적으로 new 는 A 스레드, 다음 스텝의 delete 는 B 스레드에서 일어난다.
//   Windows LFH 는 교차 해제된 블록을 할당 슬롯으로 즉시 회수하지 못하고
//   서브세그먼트를 계속 새로 확보한다 -> Private 메모리 증가.
//   팀 크기가 1 이면 실행자가 항상 동일하므로 이 현상이 없다. (관측과 일치)
//
//  [수정 내용]
//   1) updateCollisions_Nested()
//      - 힙을 건드리는 블록 전부를 #pragma omp single -> #pragma omp master 로
//        변경. master 는 이 parallel 을 만든 outer 스레드이고, 그 스레드는
//        Loop_Voxel_Unity() 의 schedule(static) 덕에 로봇별로 세션 내내 고정.
//        -> new/delete 가 항상 같은 스레드/같은 힙 슬롯에서 일어난다.
//      - master 는 암묵적 배리어가 없으므로 #pragma omp barrier 를 명시.
//      - copyprivate 는 single 전용이므로 제거. 배리어 뒤에 collisionsList.size()
//        를 그냥 읽으면 전 스레드가 같은 값을 본다.
//      - 스테일 체크(거리 비교 voxCount 회, 1us 미만)도 master 블록에 흡수.
//        병렬화 이득보다 omp for 배리어 + critical 비용이 크다. break 로 조기 탈출.
//      - critical(COLLISIONS_STALE_SET) 제거 (master 단독 실행이므로 불필요).
//      - contact force 루프는 병렬 유지 (레이스 없음이 확인됨).
//      => 이 함수의 동기화 지점 4개 -> 2개
//
//   2) doTimeStep_Nested()
//      - single copyprivate(do_coll) 제거. Phase 3 의 암묵적 배리어 뒤에서
//        diverged_flag 를 지역 const 로 스냅샷하면 전 스레드가 동일 값을 본다.
//        collisions 는 읽기 전용 멤버라 원래부터 동일 값이다.
//        -> 배리어 1개 감소 + copyprivate 브로드캐스트 의존성 제거
//      => 스텝당 배리어 10개 -> 6개
//
//   3) 진단용 계측 추가
//      - VX_PIN_INNER_THREADS : 1(기본, 기존 동작 유지) / 0(핀 호출 제거)
//        PinToGroupIfNeeded() 는 use_inner_mt 로 게이트되는 유일한 코드이며
//        최소 재현 프로그램에 없던 변수다. 위 수정 후에도 증가가 남으면
//        이 값을 0 으로 바꿔 재측정할 것.
//        (부수 효과: 지금 이 핀은 Optimize_OpenMP_ProcessGroups() 의 1:1 핀을
//         그룹 전체 마스크로 덮어쓰고, outer 스레드 0/1 이 모두 Group 0 이라
//         inner 스레드 전부가 Group 0 에 몰린다. 끄면 그 문제도 사라진다.)
//      - Get_Vx_Pin_FirstTouch_Count() : "처음 보는 스레드" 누적 카운트.
//        1초마다 찍어서 선형 증가하면 vcomp 가 중첩 팀 스레드를 재생성하는 것.
//        한 자릿수에서 멈추면 스레드는 풀링되고 있다는 뜻.
//
//  [손대지 않은 것]
//   * omp_in_parallel() == false 인 직렬 경로 (원본 그대로)
//   * Phase 1~4 의 #pragma omp for schedule(static) 4개. 링크->복셀->링크->복셀
//     의존성을 암묵적 배리어가 지켜주고 있어 하나도 제거할 수 없다.
//     schedule(static) 도 캐시 지역성 때문에 올바른 선택이므로 유지.
//   * 파일 하단의 주석처리된 구버전 블록
//
//  [검증 절차]
//   0. VoxelEngineCore.cs 의 ConnectConsoleOutput() 을 먼저 주석 처리해
//      배경 누수를 제거한다 (미검증 최우선 용의자).
//   1. inner=1 로 30분 이상 돌려 기준선이 진짜 0 인지 확인.
//   2. 이 패치 적용 후 inner=3 으로 30분. 기준선과 같으면 해결.
//   3. 그래도 증가하면 VX_PIN_INNER_THREADS 를 0 으로 두고 재측정.
//   4. VMMap Compare 로 Heap / Private Data / Thread Stack 중 어느 행이
//      늘어나는지 확인.
// ============================================================================

/*******************************************************************************
Copyright (c) 2015, Jonathan Hiller
To cite academic use of Voxelyze: Jonathan Hiller and Hod Lipson "Dynamic Simulation of Soft Multimaterial 3D-Printed Objects" Soft Robotics. March 2014, 1(1): 88-101.
Available at http://online.liebertpub.com/doi/pdfplus/10.1089/soro.2013.0010

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
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
#include <atomic>

// ----------------------------------------------------------------------------
// [PATCH] inner 스레드 프로세서 그룹 핀 사용 여부.
//   1 : 기존 동작 유지 (기본)
//   0 : PinToGroupIfNeeded() 호출 제거.
//       use_inner_mt 로 게이트되는 유일한 코드 경로를 없애는 실험용 스위치.
// ----------------------------------------------------------------------------
#define VX_PIN_INNER_THREADS 1

// [PATCH] 진단용: "처음 보는 스레드" 누적 카운트.
//   선형 증가       -> vcomp 가 중첩 팀 스레드를 매 fork 마다 재생성하고 있음
//   한 자릿수 정지  -> 스레드는 정상적으로 풀링되고 있음
static std::atomic<long long> g_vxPinFirstTouch(0);

long long Get_Vx_Pin_FirstTouch_Count(void)
{
	return g_vxPinFirstTouch.load(std::memory_order_relaxed);
}

#if VX_PIN_INNER_THREADS

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
	return ga.Group; // 현재 스레드의 primary group
}

static inline void PinToGroupIfNeeded(WORD group)
{
	// 같은 OS 스레드에서 doTimeStep이 수천 번 호출되므로,
	// 매번 SetThreadGroupAffinity를 때리지 않도록 thread_local 캐시.
	static thread_local int pinned_group = -1;
	if (pinned_group == (int)group) return;

	// [PATCH] 여기 도달했다 = 이 스레드는 처음 보는 스레드다.
	g_vxPinFirstTouch.fetch_add(1, std::memory_order_relaxed);

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

#endif // VX_PIN_INNER_THREADS
#endif // _WIN32



#include <omp.h>
#include "Voxelyze.h"


//#define VER_02_19
//#ifdef VER_02_19


bool CVoxelyze::doTimeStep_Nested(double dt)
{
    if (dt < 0) dt = recommendedTimeStep();

    const int linkCount = (int)linksList.size();
    const int voxCount  = (int)voxelsList.size();

    // 내부 병렬을 쓸지 말지 (num_thread>1 이고 is_thread=true일 때만)
    const bool use_inner_mt = (is_thread && num_thread > 1);

    // 부모(outer) 스레드가 속한 그룹을 읽어둔다.
    // outer에서 이미 SetThreadGroupAffinity로 group0/1로 박아놨다면 그 값이 그대로 나온다.
    WORD parent_group = 0;
#if defined(_WIN32) && VX_PIN_INNER_THREADS
    if (use_inner_mt) parent_group = GetMyPrimaryGroup();
#endif

    // OpenMP 2.x에서 bool atomic write 같은 걸 피하려고 int 플래그 사용
    int diverged_flag = 0;

    #pragma omp parallel num_threads(num_thread) if(use_inner_mt) firstprivate(parent_group) shared(diverged_flag)
    {
#if defined(_WIN32) && VX_PIN_INNER_THREADS
        // inner 스레드 affinity 설정 위치 (병렬영역 시작 직후)
        if (use_inner_mt) PinToGroupIfNeeded(parent_group);
#endif

        // [Phase 1] 링크 기하학 정보 갱신
        //   쓰기: 링크 자기 멤버 / 읽기: 복셀의 pos, orient
        //   링크 하나 = 스레드 하나 소유 -> 레이스 없음
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)	linksList[i]->preUpdateGeometry();
        // 암묵적 배리어: 링크 기하가 확정된다 (Phase 2 가 이것을 읽는다)

        // [Phase 2] 복셀 캐시 갱신  ★★ 삭제·병합 금지 ★★
        //   Phase 3 의 finalUpdateForces() → updateTransverseInfo()
        //   → pV->transverseArea()/transverseStrainSum() → poissonsStrain() 경로는
        //   poissonsStrainInvalid == true 이면 pStrain(24바이트)에 '쓰기'를 한다.
        //   복셀 하나는 최대 6개 링크에 속하고 그 링크들은 서로 다른 스레드에 있으므로,
        //   Phase 2 가 미리 캐시를 채워 flag 를 false 로 만들어두는 것이 유일한 방어다.
        #pragma omp for schedule(static)
        for (int i = 0; i < voxCount; i++)
            voxelsList[i]->poissonsStrain();
        // 암묵적 배리어: pStrain 이 확정된다 (Phase 3 가 이것을 읽는다)

        // [Phase 3] 힘 계산 + 발산 체크
        //   쓰기: 링크의 force/moment / 읽기: 양쪽 복셀의 pStrain
        #pragma omp for schedule(static)
        for (int i = 0; i < linkCount; i++)
        {
            linksList[i]->finalUpdateForces();

			const double as = linksList[i]->axialStrain();

			//if (linksList[i]->axialStrain() > 100) {
			if (!(as == as) || as > 100) 
			{            
                // OpenMP 2.x에서도 안전하게: atomic update
                #pragma omp atomic
                diverged_flag |= 1;
            }
        }
        // 암묵적 배리어: 링크 힘 + diverged_flag 가 확정된다

        // --------------------------------------------------------------------
        // [PATCH] single copyprivate(do_coll) 제거.
        //   diverged_flag 는 바로 위 omp for 의 암묵적 배리어에서 확정되었고
        //   collisions 는 읽기 전용 멤버이므로, 전 스레드가 이미 동일한 값을 본다.
        //   따라서 브로드캐스트가 필요 없다. 지역 const 로 스냅샷해서
        //   "팀 전체가 동일 조건으로 워크셰어링을 만난다"는 것을 코드로 못 박는다.
        //   (배리어 1개 감소 + copyprivate 의존성 제거)
        // --------------------------------------------------------------------
        const int diverged = diverged_flag;

        if (!diverged && collisions)
        {
            // 팀 전체가 호출해야 한다 (내부에 omp master/barrier/for 가 있다)
            updateCollisions_Nested();
        }

        // [Phase 4] 복셀 timestep (발산이면 모두 스킵)
        //   쓰기: 복셀 자기 pos/orient/모멘텀
        //   읽기: 6개 링크의 force(Phase 3) + colWatch 의 contactForce(위 단계)
        //         둘 다 배리어로 확정된 뒤이므로 안전
        if (!diverged)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < voxCount; i++)	voxelsList[i]->timeStep(dt);
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
    //  [PATCH] 이 경로는 원본 그대로 유지.
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
    // => 새 parallel 생성 금지, omp master / barrier / for 만 사용
    // ------------------------------------------------------------

    // ============================================================
    // [PATCH] 1) + 2) + 3) 을 하나의 master 블록으로 통합
    // ------------------------------------------------------------
    //  왜 single 이 아니라 master 인가:
    //    single 은 "먼저 도착한 스레드"가 실행하므로 스텝마다 실행자가 바뀐다.
    //    이 블록 안의 generateNearby() / regenerateCollisions() 는
    //    N x new + N x delete + 2N x push_back 규모로 힙을 대량 사용한다.
    //    실행자가 바뀌면 new 는 A 스레드, 다음 스텝 delete 는 B 스레드에서
    //    일어나고, Windows LFH 는 교차 해제된 블록을 할당 슬롯으로 회수하지
    //    못해 서브세그먼트를 계속 새로 확보한다 (= 관측된 메모리 증가).
    //
    //    master 는 이 parallel 을 만든 outer 스레드이고, Loop_Voxel_Unity() 의
    //    #pragma omp parallel for schedule(static) 덕에 로봇 T 는 항상 outer
    //    스레드 T 가 처리한다. 즉 master 는 로봇별로 세션 내내 고정된 OS
    //    스레드다. -> new/delete 가 항상 같은 힙 슬롯에서 일어난다.
    //
    //  스테일 체크도 여기 흡수한 이유:
    //    거리 비교 voxCount 회는 1us 미만이다. 병렬화해서 얻는 이득보다
    //    omp for 의 배리어 + critical 비용이 크다. 게다가 결과가 bool 하나라
    //    첫 히트에서 break 로 빠져나올 수 있다.
    // ============================================================
    #pragma omp master
    {
        // 1) 복셀 추가/삭제가 있었으면 nearby 리스트 전체 재생성
        if (nearbyStale)
        {
            for (int i = 0; i < voxCount; i++)
                voxelsList[i]->generateNearby(watchRadiusVx * 2.0, false);

            nearbyStale = false;
            collisionsStale = true;
        }

        // 2) 복셀이 충분히 움직였으면 충돌 리스트를 stale 로 표시
        if (!collisionsStale)
        {
            for (int i = 0; i < voxCount; i++)
            {
                CVX_Voxel* pV = voxelsList[i];
                if (pV->isSurface() &&
                    (pV->pos - *pV->lastColWatchPosition).Length2() > recalcDist2)
                {
                    collisionsStale = true;
                    break;                  // 하나만 찾으면 충분
                }
            }
        }

        // 3) collisionsList 재구성 (여기가 힙을 가장 많이 쓰는 곳)
        if (collisionsStale) regenerateCollisions(watchRadiusMm2);
    }
    #pragma omp barrier
    // [PATCH] master 는 암묵적 배리어가 없다. 이 배리어가 없으면
    //         아래에서 collisionsList 를 읽는 것이 즉시 데이터 레이스가 된다.

    // [PATCH] 배리어 뒤이므로 shared 멤버를 그냥 읽어도 전 스레드가 같은 값을 본다.
    //         copyprivate 불필요 (master 에서는 애초에 사용할 수도 없다).
    const int colCount = (int)collisionsList.size();

    // ============================================================
    // 4) contact force 업데이트: 병렬 유지
    // ------------------------------------------------------------
    //   CVX_Collision::updateContactForce() 는 this->force 에만 쓰고
    //   pV1/pV2 는 position()/velocity()/baseSizeAverage() 로 읽기만 한다.
    //   omp for 가 collisionsList 를 인덱스로 가르므로 충돌 객체 1개 =
    //   스레드 1개 소유가 되어 쓰기 대상이 겹치지 않는다 -> 레이스 없음.
    //   (복셀에 힘을 더하는 것은 push 가 아니라 pull 이다:
    //    CVX_Voxel::force() 가 colWatch 를 돌며 contactForce(this) 를 읽는다)
    // ============================================================
    #pragma omp for schedule(static)
    for (int i = 0; i < colCount; i++)
    {
        collisionsList[i]->updateContactForce();
    }
    // omp for 의 암묵적 배리어: Phase 4 에서 contactForce() 를 읽어도 안전
}






/*

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


*/