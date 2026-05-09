/* ---------------------------------------------------------------------------
 * RADAU5 — Rootfinding (event detection) internal header
 *
 * Declares internal functions for sign-change detection and Illinois method
 * root refinement using the RADAU5 continuous output polynomial.
 * ---------------------------------------------------------------------------*/

#ifndef RADAU5_ROOT_H_
#define RADAU5_ROOT_H_

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declaration */
typedef struct Radau5Mem_* Radau5Mem;

/* ---------------------------------------------------------------------------
 * Internal rootfinding functions (called from radau5.c)
 * ---------------------------------------------------------------------------*/

/* Initialize root function values at t0 (called once before first step) */
int radau5_root_Check1(Radau5Mem rmem);

/* Post-root re-entry: re-establish baseline after RADAU5_ROOT_RETURN */
int radau5_root_Check2(Radau5Mem rmem);

/* Post-step check: detect sign changes and refine root location */
int radau5_root_Check3(Radau5Mem rmem);

/* Free all rootfinding memory */
void radau5_root_Free(Radau5Mem rmem);

#ifdef __cplusplus
}
#endif

#endif /* RADAU5_ROOT_H_ */
