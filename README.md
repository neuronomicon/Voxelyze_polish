# [Voxelyze Engine]: OpenMP Inconsistency Fix & Large-Thread Optimization

This repository is a forked and modified version of the original [Voxelyze Engine](https://github.com/jonhiller/Voxelyze)
All original credits and copyrights belong to Jonathan Hiller and contributors.

**Modified and maintained by [Y.S.Shim](https://github.com/neuronomicon)**

### Key Modifications in this Fork
* **OpenMP Inconsistency Fix:** Resolved non-deterministic behavior when using OpenMP by implementing a 3-Phase Update architecture.
* **Precision Enhancement:** Hardcoded floating-point variables from `float` to `double` for mathematical stability.
* **Large Thread Group Fix:** Added processor group affinity binding for Windows environments to support CPUs with 64+ threads (e.g., AMD Threadripper 3990X).

**[Disclaimer]**
Please note that this fixed version is an early patch modified for a personal project. Floating-point types have been hardcoded to `double`, and some parts of the code may appear somewhat unpolished. However, the engine functions perfectly and can be used as-is without any operational issues.

### Intro
The Fixed version represents a major overhaul of the original engine, focusing on three core improvements: **1) resolving non-deterministic arithmetic issues (Inconsistency)**, **2) enhancing the physical precision of the engine**, and **3) optimizing performance and eliminating cache bottlenecks in multi-core/multi-threaded environments**. Additionally, it includes **critical fixes for thread group management when using OpenMP with large thread counts (64 or more threads, e.g., AMD Threadripper 3990X and above) in a Windows environment.**

### $\color{orange}{\text{▶ Inconsistency Problem in the Original OpenMP Multi-threaded Simulation}}$
Many researchers widely use this engine for studies related to the evolution of voxel-based robots. For these applications, it is an absolute requirement that identical voxel parameters starting from the same initial conditions produce the exact same simulation progress. Although researchers might have already been aware of this problem, to the best of my knowledge, I have not seen this issue of OpenMP inconsistency raised anywhere online or in existing literature up to this point. However, in the original version of the Voxelyze engine, enabling OpenMP introduced non-deterministic behavior, causing simulation values to fluctuate across runs despite identical starting conditions (using \#define USE_OMP). This fixed version successfully identifies the root cause of this inconsistency and restructures the underlying code. As a result, the engine now guarantees perfectly deterministic and identical simulation outcomes, even with OpenMP fully active.

## Change List

### 1. Enhancing Numerical Precision and Mathematical Stability (Precision & Math)

To reduce accumulating rounding errors during physical computations and improve overall accuracy, the fundamental data types and arithmetic logic across the system have been upgraded.

*   **Change in Default Floating-Point Type:** The default template types for core data structures like `Vec3D.h`, `Quat3D.h`, and `Array3D.h` were changed from `float` to `double`. Consequently, variables handling mass, density, stiffness, stress, and strain in `CVX_Material` and its derived classes were uniformly upgraded to `double`.
*   **Epsilon-Based Comparisons:** To prevent errors caused by direct `== 0.0` floating-point comparisons, a margin-based approach (Epsilon) was applied.
    *   In `VX_Voxel.cpp`, floor collision detection now uses `floorPenetration() >= 1E-14` instead of `floorPenetration() >= 0`.
    *   A macro `#define DP_EPSILON 1.192092896e-07` was added to `VX_MaterialLink.cpp` and `VX_Material.cpp` to perform precise floating-point comparisons, such as `fabs(f1+1.0) < DP_EPSILON`.
*   **Vector/Quaternion Operation Optimization:** 
    *   In `Vec3D.h`, the `RotX`, `RotY`, and `RotZ` methods were optimized to cache `sin` and `cos` values as local variables to avoid redundant calculations.
    *   In `Quat3D.h`, the `RotateVec3D` method was rewritten using an optimized Rodrigues' rotation formula leveraging vector cross products (`v + (t * w) + q_vec.Cross(t)`) instead of heavy matrix multiplications, significantly improving computation speed.

### 2. Preventing Memory False Sharing

Memory alignment techniques were introduced to prevent CPU cache line invalidation (false sharing), which occurs when multiple threads update adjacent objects in memory simultaneously in a multi-threaded environment.

*   **Forced Class Memory Alignment:** The `alignas(64)` attribute was applied to the declarations of core physical simulation objects, including `CVoxelyze`, `CVX_Collision`, `CVX_External`, `CVX_Link`, `CVX_Material`, and `CVX_Voxel`.
*   **Dynamic Allocation Optimization:** The `new`, `delete`, `new[]`, and `delete[]` operators were overloaded within these classes. By using `_aligned_malloc(size, 64)` and `_aligned_free` instead of the standard `malloc`, the starting memory address of instantiated objects is forced to be a multiple of the L1 cache line size (64 bytes), thereby blocking memory interference between threads.

### 3. Resolving Non-Deterministic Operations (3-Phase Update)

The critical error where simulation results varied every time due to read-write data races during OpenMP multi-threading execution has been structurally resolved.

*   **Function Splitting and State Preservation:** The massive `updateForces()` function in the original `VX_Link.cpp` was split into two separate functions: `preUpdateGeometry()` and `finalUpdateForces()`. Furthermore, to prevent the loss of velocity information required for damping calculations, member variables `Vec3D<> dPos2, dAngle1, dAngle2;` were added to `VX_Link.h` to cache intermediate values.
*   **Applying Implicit Barriers via Loop Separation:** In `Voxelyze.cpp`, the operations handled by a single parallel loop in `doTimeStep()` were separated into three distinct phases:
    1.  `linksList[i]->preUpdateGeometry()`: Updates the geometric information of each link (Write-only).
    2.  `voxelsList[i]->poissonsStrain()`: Synchronizes the Poisson's strain cache for voxels (No collision).
    3.  `linksList[i]->finalUpdateForces()`: Calculates the final force and stress based on synchronized data (Read-only).
*   This clear separation guarantees 100% identical and deterministic simulation results regardless of thread scheduling environments.

### 4. Introduction of Nested Parallelism and Core Allocation (`Voxelyze_Nested.cpp`)

To simultaneously evaluate dozens of voxel robots in algorithms like neuroevolution, a new nested parallelism architecture was added, parallelizing both the outer loop (number of robots) and the inner loop (physics engine updates).

*   **Addition of `Voxelyze_Nested.cpp`:** This newly created file implements the `doTimeStep_Nested` and `updateCollisions_Nested` functions which did not exist in the original version.
*   **Deadlock Prevention (Team-wide Execution):** Previously, collision processing used `#pragma omp single` so that only the master thread executed it, which could cause deadlocks in nested environments. In the new code, the condition is evaluated using `#pragma omp single copyprivate`, and then all threads in the inner parallel team simultaneously enter `updateCollisions_Nested()` to share the workload safely.
*   **Windows Processor Group Thread Binding (Thread Affinity):** For systems with 64 or more threads (e.g., AMD Threadripper 3990X), an issue where cores remain underutilized due to Windows OS processor group divisions was addressed. The `PinToGroupIfNeeded()` function, which calls the `SetThreadGroupAffinity` API, was implemented to explicitly bind threads to specific CPU groups, maximizing hardware resource utilization.

### 5. Data Access Optimization and Feature Expansion

*   **Removing `Array3D.h` Access Overhead:** The structure that constantly checked boundaries for safety inside `operator[]` and `operator()` was modified. For non-debug builds (when `#ifdef _DEBUG` is false), it now calls `getIndexFast()` to access memory directly, eliminating severe data retrieval bottlenecks.
*   **Permanent Preservation of Voxel Position Data (`VX_Voxel.h/cpp`):** Variables `originPos` and `originOrient`, along with methods `Set_OriginPos()` and `original_Orient()`, were added to remember the voxel's initial position and orientation. Previously, `originalPosition()` calculated the position mathematically every time it was called, but it now simply returns the cached `originPos`. Additionally, a backdoor method `Set_Pos_Direct()` was introduced to manually control voxel positions from external wrappers.
*   **Multi-thread Control Flags:** Variables `num_thread`, `is_thread`, and `is_nested` were added inside `Voxelyze.h`, providing a clean interface to dynamically control thread usage and thread count at runtime from external modules (like a DLL wrapper).

▶ Modified Dec 2025, Published Aug 2026 by Yoonsik Shim (Y.S.Shim, NeuronomicoN)
    Department of Software Engineering, Pai Chai University, Daejeon, South Korea.
