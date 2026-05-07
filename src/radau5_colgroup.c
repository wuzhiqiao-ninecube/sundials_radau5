/* ---------------------------------------------------------------------------
 * radau5_colgroup.c — Column grouping (graph coloring) for sparse DQ Jacobian
 *
 * Implements the Curtis-Powell-Reid technique using first-fit greedy
 * coloring, following MATLAB's colgroup.m approach. Two column orderings
 * are tried (natural order and reverse-degree order) and the one producing
 * fewer groups is returned.
 *
 * References:
 *   T.F. Coleman, B.S. Garbow, J.J. More, "Software for estimating
 *   sparse Jacobian matrices", ACM TOMS 11 (1984) 329-345.
 *
 *   M.W. Reichelt, L.F. Shampine, "The MATLAB ODE Suite",
 *   SIAM J. Sci. Comput. 18-1 (1997).
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include "radau5_colgroup.h"

/* ---------------------------------------------------------------------------
 * Helper: first-fit greedy coloring
 *
 * Colors columns in the order given by perm[0..n-1].
 * Uses row_offsets/row_cols (CSR-style row-to-column mapping) to detect
 * conflicts: two columns conflict if they share a nonzero row.
 *
 * group[j] is set to the 0-based group number for column j.
 * Returns the number of groups used.
 * ---------------------------------------------------------------------------*/
static sunindextype firstfit_color(
    sunindextype n,
    const sunindextype* colptrs,
    const sunindextype* rowinds,
    const sunindextype* row_offsets,
    const sunindextype* row_cols,
    const sunindextype* perm,
    sunindextype* group)
{
  sunindextype ngroups = 0;

  /* forbidden[g] = current column index if group g is forbidden for this col */
  sunindextype* forbidden = (sunindextype*)calloc((size_t)n, sizeof(sunindextype));
  if (!forbidden) return -1;

  /* Initialize: -1 means "not forbidden by anyone" */
  for (sunindextype i = 0; i < n; i++)
    forbidden[i] = -1;

  for (sunindextype idx = 0; idx < n; idx++)
  {
    sunindextype j = perm[idx]; /* column to color */

    /* Skip all-zero columns */
    if (colptrs[j + 1] == colptrs[j])
    {
      group[j] = -1;
      continue;
    }

    /* Mark forbidden groups: for each row i in column j's pattern,
     * look at all other columns that also have row i. Their groups
     * are forbidden for column j. */
    for (sunindextype p = colptrs[j]; p < colptrs[j + 1]; p++)
    {
      sunindextype row = rowinds[p];
      for (sunindextype q = row_offsets[row]; q < row_offsets[row + 1]; q++)
      {
        sunindextype other_col = row_cols[q];
        if (group[other_col] >= 0) /* already colored */
          forbidden[group[other_col]] = j; /* mark as forbidden for col j */
      }
    }

    /* Find the smallest non-forbidden group */
    sunindextype g;
    for (g = 0; g < ngroups; g++)
    {
      if (forbidden[g] != j) break;
    }
    if (g == ngroups) ngroups++;
    group[j] = g;
  }

  free(forbidden);
  return ngroups;
}

/* ---------------------------------------------------------------------------
 * Comparison function for qsort: sort by decreasing degree
 * ---------------------------------------------------------------------------*/
typedef struct { sunindextype col; sunindextype degree; } col_degree_t;

static int cmp_degree_desc(const void* a, const void* b)
{
  const col_degree_t* ca = (const col_degree_t*)a;
  const col_degree_t* cb = (const col_degree_t*)b;
  if (cb->degree != ca->degree)
    return (cb->degree > ca->degree) ? 1 : -1;
  return (ca->col > cb->col) ? 1 : ((ca->col < cb->col) ? -1 : 0);
}

/* ---------------------------------------------------------------------------
 * radau5_ComputeColGroup
 * ---------------------------------------------------------------------------*/
int radau5_ComputeColGroup(const sunindextype* colptrs,
                           const sunindextype* rowinds,
                           sunindextype n, sunindextype nnz,
                           sunindextype** col_group_out,
                           sunindextype* ngroups_out)
{
  if (n == 0) {
    *col_group_out = NULL;
    *ngroups_out = 0;
    return 0;
  }

  /* --- Step 1: Build row-to-columns reverse mapping (CSR-style) --- */
  sunindextype* row_offsets = (sunindextype*)calloc((size_t)(n + 1), sizeof(sunindextype));
  sunindextype* row_cols    = (sunindextype*)malloc((size_t)nnz * sizeof(sunindextype));
  if (!row_offsets || !row_cols) goto fail_early;

  /* Count entries per row */
  for (sunindextype k = 0; k < nnz; k++)
    row_offsets[rowinds[k] + 1]++;

  /* Prefix sum */
  for (sunindextype i = 0; i < n; i++)
    row_offsets[i + 1] += row_offsets[i];

  /* Fill row_cols */
  {
    sunindextype* tmp_offsets = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
    if (!tmp_offsets) goto fail_early;
    memcpy(tmp_offsets, row_offsets, (size_t)n * sizeof(sunindextype));

    for (sunindextype j = 0; j < n; j++)
      for (sunindextype p = colptrs[j]; p < colptrs[j + 1]; p++)
      {
        sunindextype row = rowinds[p];
        row_cols[tmp_offsets[row]++] = j;
      }
    free(tmp_offsets);
  }

  /* --- Step 2: First-fit coloring with natural order --- */
  sunindextype* perm_natural = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  sunindextype* g1 = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  if (!perm_natural || !g1) goto fail;

  for (sunindextype j = 0; j < n; j++) {
    perm_natural[j] = j;
    g1[j] = -1;
  }

  sunindextype ngroups1 = firstfit_color(n, colptrs, rowinds,
                                         row_offsets, row_cols,
                                         perm_natural, g1);
  free(perm_natural);
  if (ngroups1 < 0) goto fail;

  /* --- Step 3: First-fit coloring with reverse-degree order --- */
  col_degree_t* cd = (col_degree_t*)malloc((size_t)n * sizeof(col_degree_t));
  sunindextype* perm_degree = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  sunindextype* g2 = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  if (!cd || !perm_degree || !g2) {
    free(cd); free(perm_degree); free(g2);
    goto fail;
  }

  for (sunindextype j = 0; j < n; j++) {
    cd[j].col = j;
    cd[j].degree = colptrs[j + 1] - colptrs[j];
    g2[j] = -1;
  }
  qsort(cd, (size_t)n, sizeof(col_degree_t), cmp_degree_desc);
  for (sunindextype i = 0; i < n; i++)
    perm_degree[i] = cd[i].col;
  free(cd);

  sunindextype ngroups2 = firstfit_color(n, colptrs, rowinds,
                                         row_offsets, row_cols,
                                         perm_degree, g2);
  free(perm_degree);
  if (ngroups2 < 0) { free(g2); goto fail; }

  /* --- Step 4: Select the better result --- */
  if (ngroups1 <= ngroups2) {
    free(g2);
    *col_group_out = g1;
    *ngroups_out = ngroups1;
  } else {
    free(g1);
    *col_group_out = g2;
    *ngroups_out = ngroups2;
  }

  free(row_offsets);
  free(row_cols);
  return 0;

fail:
  free(g1);
fail_early:
  free(row_offsets);
  free(row_cols);
  return -1;
}

/* ---------------------------------------------------------------------------
 * radau5_BuildGroupLookup
 * ---------------------------------------------------------------------------*/
int radau5_BuildGroupLookup(const sunindextype* col_group,
                            sunindextype n, sunindextype ngroups,
                            sunindextype** group_offsets_out,
                            sunindextype** group_cols_out)
{
  sunindextype* offsets = (sunindextype*)calloc((size_t)(ngroups + 1),
                                                sizeof(sunindextype));
  sunindextype* cols = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  if (!offsets || !cols) { free(offsets); free(cols); return -1; }

  /* Count columns per group */
  for (sunindextype j = 0; j < n; j++)
    if (col_group[j] >= 0)
      offsets[col_group[j] + 1]++;

  /* Prefix sum */
  for (sunindextype g = 0; g < ngroups; g++)
    offsets[g + 1] += offsets[g];

  /* Scatter columns into group order */
  sunindextype* tmp = (sunindextype*)malloc((size_t)ngroups * sizeof(sunindextype));
  if (!tmp) { free(offsets); free(cols); return -1; }
  memcpy(tmp, offsets, (size_t)ngroups * sizeof(sunindextype));

  for (sunindextype j = 0; j < n; j++)
    if (col_group[j] >= 0)
      cols[tmp[col_group[j]]++] = j;

  free(tmp);
  *group_offsets_out = offsets;
  *group_cols_out = cols;
  return 0;
}
