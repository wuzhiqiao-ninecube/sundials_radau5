/* ---------------------------------------------------------------------------
 * radau5_bruss.c — 1D Brusselator problem (IVPtestset)
 *
 * Reaction-diffusion system on [0,1] with N=500 interior points:
 *   u_i' = 1 + u_i^2*v_i - 4*u_i + gamma*(u_{i-1} - 2*u_i + u_{i+1})
 *   v_i' = 3*u_i - u_i^2*v_i + gamma*(v_{i-1} - 2*v_i + v_{i+1})
 *
 * Boundary conditions: u(0)=u(1)=1, v(0)=v(1)=3
 * Initial conditions: u_i = 1 + 0.5*sin(2*pi*x_i), v_i = 3
 * Parameters: N=500, gamma = 0.02*(N+1)^2, t in [0, 10]
 *
 * System dimension: ND = 2*N = 1000 (interleaved u,v pairs)
 * Jacobian: banded with ml=2, mu=2
 *
 * Reference: IVPtestset 2.4, bruss problem
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include "radau5.h"

#define N_GRID  500
#define NEQ     (2 * N_GRID)  /* 1000 */
#define ML      2
#define MU      2
#define PI_VAL  3.14159265358979324

static sunrealtype GAMMA;
static sunrealtype GAMMA2;

/* ---------------------------------------------------------------------------
 * RHS: translates Fortran FBRUS exactly
 * ---------------------------------------------------------------------------*/
static int rhs_bruss(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  int n = N_GRID;

  for (int i = 1; i <= n; i++) {
    int iu = 2 * i - 2;  /* 0-based index for u_i */
    int iv = 2 * i - 1;  /* 0-based index for v_i */

    sunrealtype ui = yv[iu];
    sunrealtype vi = yv[iv];

    /* Boundary conditions: u(0)=1, v(0)=3, u(N+1)=1, v(N+1)=3 */
    sunrealtype uim, vim, uip, vip;
    if (i == 1) {
      uim = 1.0; vim = 3.0;
    } else {
      uim = yv[iu - 2]; vim = yv[iv - 2];
    }
    if (i == n) {
      uip = 1.0; vip = 3.0;
    } else {
      uip = yv[iu + 2]; vip = yv[iv + 2];
    }

    sunrealtype prod = ui * ui * vi;
    f[iu] = 1.0 + prod - 4.0 * ui + GAMMA * (uim - 2.0 * ui + uip);
    f[iv] = 3.0 * ui - prod + GAMMA * (vim - 2.0 * vi + vip);
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian: banded (ml=2, mu=2), translates Fortran JBRUS
 *
 * For each grid point i (1-based), the 2x2 diagonal block is:
 *   dF(iu)/dy(iu) = 2*u*v - 4 - gamma2    (band row 3, col iu)
 *   dF(iu)/dy(iv) = u^2                    (band row 2, col iv)
 *   dF(iv)/dy(iu) = 3 - 2*u*v              (band row 4, col iu)
 *   dF(iv)/dy(iv) = -u^2 - gamma2          (band row 3, col iv)
 * Off-diagonal coupling (diffusion):
 *   dF(iu)/dy(iu-2) = gamma                (band row 5, col iu-2) -> (band row 1, col iu)
 *   dF(iv)/dy(iv-2) = gamma                (band row 5, col iv-2) -> (band row 1, col iv)
 *   dF(iu)/dy(iu+2) = gamma                (band row 1, col iu+2)
 *   dF(iv)/dy(iv+2) = gamma                (band row 1, col iv+2)
 * ---------------------------------------------------------------------------*/
static int jac_bruss(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);

  int n = N_GRID;

  for (int i = 1; i <= n; i++) {
    int iu = 2 * i - 2;
    int iv = 2 * i - 1;
    sunrealtype ui = yv[iu];
    sunrealtype vi = yv[iv];
    sunrealtype uivi = ui * vi;
    sunrealtype ui2 = ui * ui;

    /* Diagonal block */
    SM_ELEMENT_B(J, iu, iu) = 2.0 * uivi - 4.0 - GAMMA2;
    SM_ELEMENT_B(J, iu, iv) = ui2;
    SM_ELEMENT_B(J, iv, iu) = 3.0 - 2.0 * uivi;
    SM_ELEMENT_B(J, iv, iv) = -ui2 - GAMMA2;

    /* Off-diagonal: zero the cross-species coupling positions */
    SM_ELEMENT_B(J, iv, iv - 1) = 0.0;  /* dF(iv)/dy(iu) already set above */
    SM_ELEMENT_B(J, iu, iu + 1) = 0.0;  /* dF(iu)/dy(iv) already set above */
  }

  /* Diffusion coupling: sub/super-diagonal at distance 2 */
  for (int i = 0; i < NEQ - 2; i++) {
    SM_ELEMENT_B(J, i + 2, i) = GAMMA;  /* lower: row i+2, col i */
    SM_ELEMENT_B(J, i, i + 2) = GAMMA;  /* upper: row i, col i+2 */
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Reference solution at t=10 (from res_exact_pic.txt)
 * 143 values: y[0..994 step 7] (Fortran: y(1..995 step 7))
 * ---------------------------------------------------------------------------*/
static const sunrealtype yref[] = {
  0.9949197002317599, 3.0213845767604077,
  0.9594350193986054, 3.0585989778165419,
  0.9243010095428502, 3.0952478919989637,
  0.8897959106772672, 3.1310118289054731,
  0.8561653620284367, 3.1656101198770159,
  0.8236197147449046, 3.1988043370624344,
  0.7923328094811884, 3.2303999530641514,
  0.7624421042573115, 3.2602463873623941,
  0.7340499750795348, 3.2882356529108807,
  0.7072259700779899, 3.3142998590079271,
  0.6820097782458483, 3.3384078449410937,
  0.6584146743834650, 3.3605612157873943,
  0.6364312187752559, 3.3807900316323134,
  0.6160310186921587, 3.3991483695914764,
  0.5971703941198909, 3.4157099395342736,
  0.5797938277687891, 3.4305638938070224,
  0.5638371159206763, 3.4438109320334580,
  0.5492301695479158, 3.4555597666485198,
  0.5358994429426996, 3.4659239846027008,
  0.5237699892215797, 3.4750193162238476,
  0.5127671585747183, 3.4829613034792271,
  0.5028179665048467, 3.4898633463634923,
  0.4938521662914935, 3.4958350971335204,
  0.4858030633656755, 3.5009811668111510,
  0.4786081100251151, 3.5054001059792705,
  0.4722093177200750, 3.5091836216744015,
  0.4665535216425440, 3.5124159935026285,
  0.4615925290790646, 3.5151736544621075,
  0.4572831793403656, 3.5175249049438184,
  0.4535873393501199, 3.5195297317024448,
  0.4504718553589467, 3.5212397070273984,
  0.4479084778719241, 3.5226979467564341,
  0.4458737738041973, 3.5239391090719634,
  0.4443490371324889, 3.5249894191569453,
  0.4433202068820853, 3.5258667077466495,
  0.4427777991494095, 3.5265804544017270,
  0.4427168579654424, 3.5271318289682063,
  0.4431369281018266, 3.5275137272135266,
  0.4440420513508381, 3.5277107990730161,
  0.4454407863109616, 3.5276994703501980,
  0.4473462502188303, 3.5274479611304068,
  0.4497761798232572, 3.5269163066324394,
  0.4527530066369863, 3.5260563887768472,
  0.4563039400688689, 3.5248119894251024,
  0.4604610498812091, 3.5231188790654930,
  0.4652613370907894, 3.5209049576992761,
  0.4707467798082714, 3.5180904678044698,
  0.4769643375804777, 3.5145883024867057,
  0.4839658945842979, 3.5103044351908528,
  0.4918081185812277, 3.5051385005173827,
  0.5005522089940899, 3.4989845585737802,
  0.5102635039989190, 3.4917320776245013,
  0.5210109134090777, 3.4832671712209993,
  0.5328661417420937, 3.4734741260299615,
  0.5459026646938675, 3.4622372546582585,
  0.5601944229089820, 3.4494431032230182,
  0.5758142001453760, 3.4349830354873361,
  0.5928316594749734, 3.4187562033108012,
  0.6113110218368440, 3.4006728962523969,
  0.6313083867734524, 3.3806582409098729,
  0.6528687160104193, 3.3586561928427350,
  0.6760225267555723, 3.3346337311179157,
  0.7007823726569721, 3.3085851288057930,
  0.7271392249346637, 3.2805361342349380,
  0.7550589020044152, 3.2505478606008622,
  0.7844787296769868, 3.2187201496972175,
  0.8153046416214843, 3.1851941538893653,
  0.8474089465959840, 3.1501538739882800,
  0.8806289904192589, 3.1138264039027113,
  0.9147669230929857, 3.0764806689389470,
  0.9495907429372025, 3.0384245041548366,
  0.9848367306701233
};
#define NREF 143

/* Reference indices: Fortran outputs y(1..995 step 7), i.e. 0-based indices 0..994 step 7 */
static void get_ref_indices(int* indices)
{
  int cnt = 0;
  for (int i = 0; i < 995; i += 7)
    indices[cnt++] = i;
}

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 13;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  /* Compute derived constants */
  double anp1 = (double)(N_GRID + 1);
  double usdelq = anp1 * anp1;
  GAMMA = 0.02 * usdelq;
  GAMMA2 = 2.0 * GAMMA;

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Initial conditions: u_i = 1 + 0.5*sin(2*pi*x_i), v_i = 3 */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  for (int i = 1; i <= N_GRID; i++) {
    double xi = (double)i / anp1;
    y0v[2 * i - 2] = 1.0 + 0.5 * sin(2.0 * PI_VAL * xi);
    y0v[2 * i - 1] = 3.0;
  }

  /* Solver setup */
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_bruss, 0.0, y0);

  /* Band Jacobian: ml=2, mu=2 */
  SUNMatrix Jt = SUNBandMatrix(NEQ, ML, MU, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_bruss);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Solve to t=10 */
  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 10.0, yout, &tret);

  /* Compare with reference at selected indices */
  sunrealtype* yd = N_VGetArrayPointer(yout);
  int ref_indices[NREF];
  get_ref_indices(ref_indices);

  sunrealtype maxerr = 0.0;
  for (int r = 0; r < NREF; r++) {
    int idx = ref_indices[r];
    sunrealtype err = fabs(yd[idx] - yref[r]);
    if (err > maxerr) maxerr = err;
  }

  printf("=== Brusselator 1D (n=%d, rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         NEQ, rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.10e\n", ret, tret);
  printf("max |y - yref| at ref points = %.6e\n", maxerr);

  /* Print first few u,v pairs */
  printf("\nFirst 5 grid points (u, v):\n");
  for (int i = 0; i < 5; i++)
    printf("  u[%d] = %18.14e  v[%d] = %18.14e\n",
           i + 1, yd[2 * i], i + 1, yd[2 * i + 1]);

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("\nnstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  int pass = (ret == 0 && maxerr < 1.0e-3);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
