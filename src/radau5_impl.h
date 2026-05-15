/* ---------------------------------------------------------------------------
 * RADAU5 internal data structures
 * ---------------------------------------------------------------------------*/

#ifndef RADAU5_IMPL_H_
#define RADAU5_IMPL_H_

#include "radau5.h"
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of stages supported */
#define RADAU5_NS_MAX 7
#define RADAU5_NPAIRS_MAX 3  /* (RADAU5_NS_MAX-1)/2 */

/* ---------------------------------------------------------------------------
 * Solver memory structure
 * ---------------------------------------------------------------------------*/
typedef struct Radau5Mem_
{
  SUNContext sunctx;
  sunindextype n; /* problem dimension */

  /* Variable-order support */
  int ns;      /* current number of stages (3, 5, or 7) */
  int npairs;  /* (ns-1)/2: number of complex eigenvalue pairs */

  /* User functions */
  Radau5RhsFn rhs;
  Radau5JacFn jac;       /* NULL → DQ Jacobian */
  Radau5MassFn mas;      /* NULL → identity mass */
  Radau5SolOutFn solout; /* NULL → no output callback */
  void* user_data;

  /* Current time and state */
  sunrealtype tn; /* current time */
  N_Vector ycur;  /* current solution (owned) */
  N_Vector fn;    /* f(tn, ycur) */

  /* SUNDIALS objects — Jacobian and mass */
  SUNMatrix J;       /* Jacobian template (cloned from user) */
  SUNMatrix Jsaved;  /* saved Jacobian for reuse */
  SUNMatrix M;       /* mass matrix (NULL if identity) */

  /* Real system E1 (n×n): E1 = fac1*M - J */
  SUNMatrix E1;
  SUNLinearSolver LS_E1;

  /* Realified complex systems E2 (2n×2n) — one per complex eigenvalue pair */
  SUNMatrix E2[RADAU5_NPAIRS_MAX];
  SUNLinearSolver LS_E2[RADAU5_NPAIRS_MAX];

  /* Stage increments Z_i = Y_i - y
   * Ordering: z[0]=real eigenvalue, z[1],z[2]=complex pair 0, etc. */
  N_Vector z[RADAU5_NS_MAX];
  /* Transformed stage increments in TI-space */
  N_Vector f[RADAU5_NS_MAX];

  /* Continuous output: cont[0]=ycur copy, cont[1..ns]=divided-diff coeffs */
  N_Vector cont[RADAU5_NS_MAX + 1];

  /* Error weight vector: scal_i = atol_i + rtol_i * |y_i| */
  N_Vector scal;
  N_Vector ewt; /* alias or separate for DQ Jacobian */

  /* Scratch vectors */
  N_Vector tmp1, tmp2, tmp3;

  /* For realified complex solve (length 2n) — one per pair */
  N_Vector rhs2[RADAU5_NPAIRS_MAX];
  N_Vector sol2[RADAU5_NPAIRS_MAX];
  N_Vector y2n[RADAU5_NPAIRS_MAX];

  /* DAE support */
  N_Vector id; /* differential/algebraic indicator (NULL if ODE) */

  /* Tolerances */
  sunrealtype rtol_s; /* scalar rtol */
  sunrealtype atol_s; /* scalar atol */
  N_Vector rtol_v;    /* vector rtol (NULL if scalar) */
  N_Vector atol_v;    /* vector atol (NULL if scalar) */
  int itol;           /* 0=scalar, 1=vector */

  /* Method constants (variable-order Radau IIA, ns=3,5,7) */
  sunrealtype c[RADAU5_NS_MAX];         /* collocation nodes [ns] */
  sunrealtype dd[RADAU5_NS_MAX];        /* error estimation coefficients [ns] */
  sunrealtype u1;                        /* real eigenvalue of A^{-1} */
  sunrealtype alph[RADAU5_NPAIRS_MAX];  /* complex eigenvalue real parts */
  sunrealtype beta_eig[RADAU5_NPAIRS_MAX]; /* complex eigenvalue imag parts */
  sunrealtype T_mat[RADAU5_NS_MAX * RADAU5_NS_MAX];   /* T[ns*ns] row-major */
  sunrealtype TI_mat[RADAU5_NS_MAX * RADAU5_NS_MAX];  /* TI[ns*ns] row-major */

  /* Schur decomposition option */
  int use_schur;
  sunrealtype US_mat[RADAU5_NS_MAX * RADAU5_NS_MAX];  /* Schur vectors */
  sunrealtype TS_mat[RADAU5_NS_MAX * RADAU5_NS_MAX];  /* Schur form */

  /* Step size control */
  sunrealtype h;
  sunrealtype hold;
  sunrealtype hmax;
  sunrealtype hopt;
  sunrealtype safe;         /* safety factor (default 0.9) */
  sunrealtype facl;         /* max step ratio (default 8.0) */
  sunrealtype facr;         /* min step ratio (default 0.2) */
  sunrealtype fnewt;        /* Newton convergence threshold */
  sunrealtype thet;         /* Jacobian reuse threshold (default 0.001) */
  sunrealtype quot1, quot2; /* step ratio bounds for Jac reuse (1,6) */
  int pred;                 /* 1=Gustafsson, 2=classical */
  sunrealtype hacc, erracc; /* Gustafsson state */
  sunrealtype faccon;       /* Newton convergence factor */
  sunrealtype theta;        /* Newton contraction rate */
  sunrealtype hhfac;        /* step-size reduction factor predicted by Newton */

  /* DAE index support */
  sunindextype nind1, nind2, nind3;

  /* Newton iteration */
  int nit;    /* max Newton iterations (default 7) */
  int startn; /* 0=extrapolate, 1=zero starting values */

  /* Variable-order selection */
  int nsmin, nsmax;        /* order bounds (default both 3) */
  int variab;              /* 1 if variable order enabled (nsmin < nsmax) */
  int ichan;               /* steps since last order change */
  sunrealtype thetat;      /* smoothed Newton contraction rate */
  int nsnew;               /* proposed new ns */
  int change;              /* 1 if order change pending this step */
  int unexp;               /* 1 if unexpected step rejection occurred */
  int unexn;               /* 1 if Newton predicted hhfac <= 0.5 */
  int ikeep;               /* 1 if skipdecomp path wants to retry with higher order */
  int newt_prev;           /* Newton iteration count from previous accepted step */
  sunrealtype vitu;        /* theta threshold to increase order (default 0.002) */
  sunrealtype vitd;        /* theta threshold to decrease order (default 0.8) */
  sunrealtype hhou;        /* step ratio upper bound for order increase (default 1.2) */
  sunrealtype hhod;        /* step ratio lower bound for order increase (default 0.8) */

  /* Counters */
  long int nstep, naccpt, nrejct;
  long int nfcn, njac, ndec, nsol, nnewt;
  int nsing; /* consecutive singular matrix count */

  /* State flags */
  int first;
  int reject;
  int last;
  int caljac;
  int skipdecomp;      /* 1 = skip both Jacobian and decomposition */
  sunrealtype hnew_save; /* saved hnew when skipdecomp=1 (h kept unchanged) */
  int tol_transformed;  /* 1 = tolerances already transformed */
  int setup_done;      /* Radau5Init completed */
  int ls_done;         /* Radau5SetLinearSolver completed */
  int mass_evaluated;  /* 1 if mass matrix has been evaluated */
  long int mxstep;     /* max steps (default 100000) */

  /* Continuous output state */
  sunrealtype xsol, xold, hsol;

  /* Matrix type info */
  SUNMatrix_ID mat_id;
  sunindextype mu, ml; /* band widths (if band matrix) */

  /* Column grouping for sparse DQ Jacobian */
  sunindextype* col_group;       /* length n: group number for each column, -1 = all-zero */
  sunindextype  ngroups;         /* total number of groups */
  sunindextype* group_offsets;   /* length ngroups+1: CSC-style offsets into group_cols */
  sunindextype* group_cols;      /* length n: column indices ordered by group */
  sunindextype* sp_colptrs;      /* length n+1: copy of sparsity pattern column pointers */
  sunindextype* sp_rowinds;      /* length sp_nnz: copy of sparsity pattern row indices */
  sunindextype  sp_nnz;          /* number of nonzeros in sparsity pattern */

  /* Rootfinding (event detection) */
  Radau5RootFn gfun;          /* user root function (NULL = inactive) */
  int          nrtfn;         /* number of root functions */
  sunrealtype* glo;           /* g values at step start [nrtfn] */
  sunrealtype* ghi;           /* g values at step end [nrtfn] */
  sunrealtype* grout;         /* g values at located root [nrtfn] */
  int*         iroots;        /* output: +1 rising, -1 falling, 0 none [nrtfn] */
  int*         rootdir;       /* direction filter: +1, -1, 0 [nrtfn] */
  int*         gactive;       /* 1=active, 0=masked [nrtfn] */
  sunrealtype  troot;         /* located root time */
  N_Vector     y_root;        /* interpolated y at troot (scratch) */
  long int     nge;           /* total g evaluations */
  int          root_active;   /* 1 after Radau5RootInit with nrtfn>0 */
  int          root_init_done;/* 1 after radau5_root_Check1 has run */
  int          irfnd;         /* 1 = last return was ROOT_RETURN */

} * Radau5Mem;

/* ---------------------------------------------------------------------------
 * Access macro
 * ---------------------------------------------------------------------------*/
#define RADAU5_MEM(mem) ((Radau5Mem)(mem))

/* ---------------------------------------------------------------------------
 * Internal function prototypes
 * ---------------------------------------------------------------------------*/

/* radau5_linsys.c */
int radau5_InitConstants(Radau5Mem rmem);
int radau5_ChangeOrder(Radau5Mem rmem, int nsnew);
int radau5_DQJacDense(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy);
int radau5_DQJacBand(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy);
int radau5_DQJacSparse(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy);
int radau5_BuildE1(Radau5Mem rmem, sunrealtype fac1);
int radau5_BuildE2(Radau5Mem rmem, int pair_idx, sunrealtype alphn, sunrealtype betan);
int radau5_DecompE1(Radau5Mem rmem);
int radau5_DecompE2(Radau5Mem rmem, int pair_idx);
int radau5_MassMult(Radau5Mem rmem, N_Vector x, N_Vector result);
SUNMatrix radau5_SparseUnion(SUNMatrix A, SUNMatrix B, SUNContext sunctx);

/* radau5_newt.c */
int radau5_Newton(Radau5Mem rmem, int* newt_out);

/* radau5_estrad.c */
int radau5_ErrorEstimate(Radau5Mem rmem, sunrealtype* err);

/* radau5_step.c */
int radau5_Step(Radau5Mem rmem);

/* radau5_contr.c */
void radau5_UpdateContinuousOutput(Radau5Mem rmem);

/* radau5_ic.c */
int radau5_CalcIC(Radau5Mem rmem, N_Vector id);

/* radau5_root.c */
int radau5_root_Check1(Radau5Mem rmem);
int radau5_root_Check2(Radau5Mem rmem);
int radau5_root_Check3(Radau5Mem rmem);
void radau5_root_Free(Radau5Mem rmem);

/* Utility */
void radau5_ComputeScal(Radau5Mem rmem, N_Vector y);

#ifdef __cplusplus
}
#endif

#endif /* RADAU5_IMPL_H_ */
