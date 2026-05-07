/* ---------------------------------------------------------------------------
 * test_radau5_colgroup.c — Unit tests for radau5_colgroup
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "radau5_colgroup.h"
#include <sundials/sundials_types.h>

static int nfail = 0;

#define CHECK(cond, msg)                                          \
  do {                                                            \
    if (!(cond)) {                                                \
      fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__, msg); \
      nfail++;                                                    \
    }                                                             \
  } while (0)

/* ---------------------------------------------------------------------------
 * verify_no_conflicts: for each group, check no two columns share a row
 * ---------------------------------------------------------------------------*/
static int verify_no_conflicts(const sunindextype* col_group,
                                sunindextype n,
                                const sunindextype* colptrs,
                                const sunindextype* rowinds,
                                sunindextype ngroups)
{
  /* For each row, track which group last claimed it */
  sunindextype* row_group = (sunindextype*)malloc((size_t)n * sizeof(sunindextype));
  if (!row_group) return -1;
  for (sunindextype i = 0; i < n; i++) row_group[i] = -1;

  int ok = 1;
  for (sunindextype j = 0; j < n; j++)
  {
    sunindextype g = col_group[j];
    if (g < 0) continue; /* empty column */
    for (sunindextype p = colptrs[j]; p < colptrs[j + 1]; p++)
    {
      sunindextype row = rowinds[p];
      if (row_group[row] == g)
      {
        fprintf(stderr, "  conflict: cols %lld and %lld both in group %lld share row %lld\n",
                (long long)j, (long long)row_group[row], (long long)g, (long long)row);
        ok = 0;
      }
      row_group[row] = g;
    }
  }
  free(row_group);
  return ok;
}

/* ---------------------------------------------------------------------------
 * verify_group_lookup: check group_offsets/group_cols are consistent
 * ---------------------------------------------------------------------------*/
static int verify_group_lookup(const sunindextype* col_group,
                                sunindextype n,
                                sunindextype ngroups,
                                const sunindextype* group_offsets,
                                const sunindextype* group_cols)
{
  /* Count expected columns per group */
  sunindextype* cnt = (sunindextype*)calloc((size_t)ngroups, sizeof(sunindextype));
  if (!cnt) return -1;
  for (sunindextype j = 0; j < n; j++)
    if (col_group[j] >= 0) cnt[col_group[j]]++;

  int ok = 1;
  for (sunindextype g = 0; g < ngroups; g++)
  {
    sunindextype span = group_offsets[g + 1] - group_offsets[g];
    if (span != cnt[g])
    {
      fprintf(stderr, "  group %lld: expected %lld cols, got %lld\n",
              (long long)g, (long long)cnt[g], (long long)span);
      ok = 0;
    }
    /* Each listed column must map back to this group */
    for (sunindextype k = group_offsets[g]; k < group_offsets[g + 1]; k++)
    {
      sunindextype c = group_cols[k];
      if (col_group[c] != g)
      {
        fprintf(stderr, "  group_cols[%lld]=%lld but col_group[%lld]=%lld\n",
                (long long)k, (long long)c, (long long)c, (long long)col_group[c]);
        ok = 0;
      }
    }
  }
  free(cnt);
  return ok;
}

/* ---------------------------------------------------------------------------
 * Test 1: Diagonal matrix (n=10) → ngroups = 1
 * ---------------------------------------------------------------------------*/
static void test_diagonal(void)
{
  const sunindextype n = 10;
  sunindextype colptrs[11], rowinds[10];
  for (sunindextype j = 0; j < n; j++) {
    colptrs[j] = j;
    rowinds[j] = j;
  }
  colptrs[n] = n;

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(colptrs, rowinds, n, n, &col_group, &ngroups);
  CHECK(ret == 0, "diagonal: ComputeColGroup returned error");
  CHECK(ngroups == 1, "diagonal: expected ngroups=1");
  CHECK(verify_no_conflicts(col_group, n, colptrs, rowinds, ngroups) == 1,
        "diagonal: conflict detected");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "diagonal: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "diagonal: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * Test 2: Full matrix (n=5) → ngroups = 5
 * ---------------------------------------------------------------------------*/
static void test_full(void)
{
  const sunindextype n = 5;
  const sunindextype nnz = n * n;
  sunindextype colptrs[6];
  sunindextype rowinds[25];
  for (sunindextype j = 0; j < n; j++) {
    colptrs[j] = j * n;
    for (sunindextype i = 0; i < n; i++)
      rowinds[j * n + i] = i;
  }
  colptrs[n] = nnz;

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(colptrs, rowinds, n, nnz, &col_group, &ngroups);
  CHECK(ret == 0, "full: ComputeColGroup returned error");
  CHECK(ngroups == 5, "full: expected ngroups=5");
  CHECK(verify_no_conflicts(col_group, n, colptrs, rowinds, ngroups) == 1,
        "full: conflict detected");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "full: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "full: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * Test 3: Tridiagonal (n=10) → ngroups = 3
 * ---------------------------------------------------------------------------*/
static void test_tridiagonal(void)
{
  const sunindextype n = 10;
  /* Build CSC tridiagonal: col j has rows max(0,j-1)..min(n-1,j+1) */
  sunindextype colptrs[11];
  sunindextype rowinds[28]; /* at most 3 per col, 2 for corners */
  sunindextype nnz = 0;
  colptrs[0] = 0;
  for (sunindextype j = 0; j < n; j++) {
    if (j > 0)     rowinds[nnz++] = j - 1;
    rowinds[nnz++] = j;
    if (j < n - 1) rowinds[nnz++] = j + 1;
    colptrs[j + 1] = nnz;
  }

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(colptrs, rowinds, n, nnz, &col_group, &ngroups);
  CHECK(ret == 0, "tridiag: ComputeColGroup returned error");
  CHECK(ngroups == 3, "tridiag: expected ngroups=3");
  CHECK(verify_no_conflicts(col_group, n, colptrs, rowinds, ngroups) == 1,
        "tridiag: conflict detected");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "tridiag: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "tridiag: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * Test 4: Pollu sparsity pattern (n=20, nnz=86)
 * ---------------------------------------------------------------------------*/
static const sunindextype pollu_colptrs[21] = {
  0, 10, 19, 21, 27, 31, 41, 45, 46, 52, 56, 62, 63, 66, 69, 70, 73, 77, 78, 83, 86
};
static const sunindextype pollu_rowinds[86] = {
   0,  1,  2,  3,  5, 10, 12, 14, 18, 19,
   0,  1,  3,  4,  5,  9, 10, 11, 13,
   2,  3,
   0,  1,  2,  3, 15, 18,
   0,  1,  4,  5,
   0,  4,  5,  6,  7,  8, 10, 14, 16, 17,
   4,  5,  6,  7,
   7,
   4,  5,  7,  8,  9, 10,
   0,  1,  9, 13,
   0,  1,  9, 10, 11, 12,
  11,
   0, 10, 12,
   4,  6, 13,
  14,
   2,  5, 15,
   4,  5, 16, 17,
  17,
   0,  1,  2, 18, 19,
   0, 18, 19
};

static void test_pollu(void)
{
  const sunindextype n = 20;
  const sunindextype nnz = 86;

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(pollu_colptrs, pollu_rowinds, n, nnz,
                                   &col_group, &ngroups);
  CHECK(ret == 0, "pollu: ComputeColGroup returned error");
  CHECK(ngroups < n, "pollu: expected ngroups < 20");
  CHECK(verify_no_conflicts(col_group, n, pollu_colptrs, pollu_rowinds, ngroups) == 1,
        "pollu: conflict detected");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "pollu: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "pollu: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * Test 5: Matrix with all-zero columns (cols 1,3 empty in 5x5)
 * ---------------------------------------------------------------------------*/
static void test_empty_cols(void)
{
  /* 5x5, cols 1 and 3 are empty */
  const sunindextype n = 5;
  /* col 0: row 0; col 1: empty; col 2: row 2; col 3: empty; col 4: row 4 */
  sunindextype colptrs[6] = {0, 1, 1, 2, 2, 3};
  sunindextype rowinds[3] = {0, 2, 4};
  const sunindextype nnz = 3;

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(colptrs, rowinds, n, nnz, &col_group, &ngroups);
  CHECK(ret == 0, "empty_cols: ComputeColGroup returned error");
  CHECK(col_group[1] == -1, "empty_cols: col 1 should have group=-1");
  CHECK(col_group[3] == -1, "empty_cols: col 3 should have group=-1");
  CHECK(col_group[0] >= 0, "empty_cols: col 0 should have valid group");
  CHECK(col_group[2] >= 0, "empty_cols: col 2 should have valid group");
  CHECK(col_group[4] >= 0, "empty_cols: col 4 should have valid group");
  CHECK(verify_no_conflicts(col_group, n, colptrs, rowinds, ngroups) == 1,
        "empty_cols: conflict detected");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "empty_cols: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "empty_cols: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * Test 6: Single column (n=1) → ngroups = 1
 * ---------------------------------------------------------------------------*/
static void test_single(void)
{
  const sunindextype n = 1;
  sunindextype colptrs[2] = {0, 1};
  sunindextype rowinds[1] = {0};

  sunindextype* col_group = NULL;
  sunindextype ngroups = 0;
  int ret = radau5_ComputeColGroup(colptrs, rowinds, n, 1, &col_group, &ngroups);
  CHECK(ret == 0, "single: ComputeColGroup returned error");
  CHECK(ngroups == 1, "single: expected ngroups=1");
  CHECK(col_group[0] == 0, "single: col 0 should be in group 0");

  sunindextype* go = NULL; sunindextype* gc = NULL;
  ret = radau5_BuildGroupLookup(col_group, n, ngroups, &go, &gc);
  CHECK(ret == 0, "single: BuildGroupLookup returned error");
  CHECK(verify_group_lookup(col_group, n, ngroups, go, gc) == 1,
        "single: group lookup inconsistent");

  free(col_group); free(go); free(gc);
}

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(void)
{
  test_diagonal();
  test_full();
  test_tridiagonal();
  test_pollu();
  test_empty_cols();
  test_single();

  if (nfail == 0)
    printf("test_radau5_colgroup: PASSED\n");
  else
    printf("test_radau5_colgroup: %d FAILURES\n", nfail);

  return nfail;
}
