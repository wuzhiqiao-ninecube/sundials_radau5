/* ---------------------------------------------------------------------------
 * radau5_colgroup.h — Column grouping for sparse DQ Jacobian
 *
 * Implements the Curtis-Powell-Reid (CPR) column grouping technique:
 * given a CSC sparsity pattern, groups structurally independent columns
 * so that multiple Jacobian columns can be computed per RHS evaluation.
 * ---------------------------------------------------------------------------*/

#ifndef RADAU5_COLGROUP_H_
#define RADAU5_COLGROUP_H_

#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute column grouping from a CSC sparsity pattern.
 *
 * Tries two orderings (natural + reverse-degree) with first-fit greedy
 * coloring and returns whichever produces fewer groups.
 *
 * Input:
 *   colptrs[n+1], rowinds[nnz] — CSC sparsity pattern
 *   n    — number of columns (= number of rows, square matrix)
 *   nnz  — number of nonzeros
 *
 * Output:
 *   col_group_out — allocated array of length n. col_group[j] = group
 *                   number (0-based) for column j, or -1 if column j
 *                   has no nonzero entries. Caller must free().
 *   ngroups_out   — total number of groups
 *
 * Returns 0 on success, -1 on memory allocation failure.
 */
int radau5_ComputeColGroup(const sunindextype* colptrs,
                           const sunindextype* rowinds,
                           sunindextype n, sunindextype nnz,
                           sunindextype** col_group_out,
                           sunindextype* ngroups_out);

/* Build group lookup arrays from a column grouping.
 *
 * Input:
 *   col_group[n] — group assignment (from radau5_ComputeColGroup)
 *   n            — number of columns
 *   ngroups      — number of groups
 *
 * Output:
 *   group_offsets_out — allocated array of length ngroups+1 (CSC-style).
 *                       group_offsets[g] = start index in group_cols for
 *                       group g. Caller must free().
 *   group_cols_out    — allocated array of length n. Column indices
 *                       ordered by group. Caller must free().
 *
 * Returns 0 on success, -1 on memory allocation failure.
 */
int radau5_BuildGroupLookup(const sunindextype* col_group,
                            sunindextype n, sunindextype ngroups,
                            sunindextype** group_offsets_out,
                            sunindextype** group_cols_out);

#ifdef __cplusplus
}
#endif

#endif /* RADAU5_COLGROUP_H_ */
