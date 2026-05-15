/* ---------------------------------------------------------------------------
 * radau5_linsys.c — Linear system infrastructure for RADAU5
 *
 * Implements:
 *   radau5_InitConstants  — Radau IIA method constants (from Fortran radau5.f)
 *   radau5_DQJacDense     — Dense finite-difference Jacobian (CVODE-style)
 *   radau5_BuildE1        — Assemble E1 = fac1*M - J  (n×n real system)
 *   radau5_BuildE2        — Assemble E2 realified 2n×2n complex system
 *   radau5_DecompE1       — Factor E1 via SUNLinSolSetup
 *   radau5_DecompE2       — Factor E2 via SUNLinSolSetup
 *   radau5_ComputeScal    — Error weight vector scal[i] = atol[i] + rtol[i]*|y[i]|
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "radau5_impl.h"

/* ===========================================================================
 * Precomputed method constants for ns=5 and ns=7
 * Source: Fortran radau.f (Hairer & Wanner), subroutines COERTV and COERCV
 *
 * Stage ordering convention (matching Fortran):
 *   Row/index 0: real eigenvalue component (solved with E1)
 *   Rows 1,2: complex pair 0 (solved with E2)
 *   Rows 3,4: complex pair 1 (for ns=5)
 *   Rows 5,6: complex pair 2 (for ns=7)
 * ===========================================================================*/

/* --- NS=5 constants --- */
static const sunrealtype C5[5] = {
  0.5710419611451768219312e-01, 0.2768430136381238276800e+00,
  0.5835904323689168200567e+00, 0.8602401356562194478479e+00, 1.0
};
static const sunrealtype DD5[5] = {
  -0.2778093394406463730479e+02,  0.3641478498049213152712e+01,
  -0.1252547721169118720491e+01,  0.5920031671845428725662e+00,
  -0.2000000000000000000000e+00
};
static const sunrealtype U1_5  = 0.6286704751729276645173e+01;
static const sunrealtype ALPH5[2] = { 0.3655694325463572258243e+01, 0.5700953298671789419170e+01 };
static const sunrealtype BETA5[2] = { 0.6543736899360077294021e+01, 0.3210265600308549888425e+01 };

/* T5: 5x5 row-major. Row 4 has implicit entries from Fortran back-transform:
 *   ZZ(N4) = T551*Z1 + 1*Z2 + 0*Z3 + 1*Z4 + 0*Z5 */
static const sunrealtype T5[25] = {
  -0.1251758622050104589014e-01, -0.1024204781790882707009e-01,  0.4767387729029572386318e-01, -0.1147851525522951470794e-01, -0.1401985889287541028108e-01,
  -0.1491670151895382429004e-02,  0.5017286451737105816299e-01, -0.9433181918161143698066e-01, -0.7668830749180162885157e-02,  0.2470857842651852681253e-01,
  -0.7298187638808714862266e-01, -0.2305395340434179467214e+00,  0.1027030453801258997922e+00,  0.1939846399882895091122e-01,  0.8180035370375117083639e-01,
  -0.3800914400035681041264e+00,  0.3778939022488612495439e+00,  0.4667441303324943592896e+00,  0.4076011712801990666217e+00,  0.1996824278868025259365e+00,
  -0.9219789736812104884883e+00,  1.0, 0.0, 1.0, 0.0
};
static const sunrealtype TI5[25] = {
  -0.3004156772154440162771e+02, -0.1386510785627141316518e+02, -0.3480002774795185561828e+01,  0.1032008797825263422771e+01, -0.8043030450739899174753e+00,
   0.5344186437834911598895e+01,  0.4593615567759161004454e+01, -0.3036360323459424298646e+01,  0.1050660190231458863860e+01, -0.2727786118642962705386e+00,
   0.3748059807439804860051e+01, -0.3984965736343884667252e+01, -0.1044415641608018792942e+01,  0.1184098568137948487231e+01, -0.4499177701567803688988e+00,
  -0.3304188021351900000806e+02, -0.1737695347906356701945e+02, -0.1721290632540055611515e+00, -0.9916977798254264258817e-01,  0.5312281158383066671849e+00,
  -0.8611443979875291977700e+01,  0.9699991409528808231336e+01,  0.1914728639696874284851e+01,  0.2418692006084940026427e+01, -0.1047463487935337418694e+01
};

/* Schur decomposition for ns=5 (user-provided) */
static const sunrealtype US5[25] = {
   104.376011865108e-003,    2.75321569933667e-003,    74.5464030994560e-003,   -411.612705588829e-003,   -902.283703905388e-003,
  -216.061922297273e-003,   -32.2153757820648e-003,    103.014361885950e-003,    876.499915903360e-003,   -416.431691833195e-003,
   291.316048441817e-003,    192.522319893028e-003,   -921.003086216347e-003,    137.473507673671e-003,   -104.520170929409e-003,
   869.441565532767e-003,   -401.017866194355e-003,    215.864391849722e-003,    188.938760166232e-003,    30.9958764201555e-003,
  -318.793378106442e-003,   -895.027606670579e-003,   -298.306477538585e-003,   -87.8979955060830e-003,   -24.1568461838729e-003
};
static const sunrealtype TS5[25] = {
   3.65569432546356e+000,   -2.76149775889893e+000,    9.42380616684492e+000,    4.78682213405259e+000,    1.66044337303607e+000,
   15.5062565124514e+000,    3.65569432546356e+000,    12.7085144861966e+000,    8.98581391232333e+000,    4.21459604790153e+000,
   0.0,    0.0,    5.70095329867178e+000,    9.65721699665011e+000,    4.07657545572918e+000,
   0.0,    0.0,   -1.06716098728022e+000,    5.70095329867178e+000,    9.26586374912508e+000,
   0.0,    0.0,    0.0,    0.0,    6.28670475172929e+000
};

/* --- NS=7 constants --- */
static const sunrealtype C7[7] = {
  0.2931642715978489197205e-01, 0.1480785996684842918500e+00,
  0.3369846902811542990971e+00, 0.5586715187715501320814e+00,
  0.7692338620300545009169e+00, 0.9269456713197411148519e+00, 1.0
};
static const sunrealtype DD7[7] = {
  -0.5437443689412861451458e+02,  0.7000024004259186512041e+01,
  -0.2355661091987557192256e+01,  0.1132289066106134386384e+01,
  -0.6468913267673587118673e+00,  0.3875333853753523774248e+00,
  -0.1428571428571428571429e+00
};
static const sunrealtype U1_7  = 0.8936832788405216337302e+01;
static const sunrealtype ALPH7[3] = { 0.4378693561506806002523e+01, 0.7141055219187640105775e+01, 0.8511834825102945723051e+01 };
static const sunrealtype BETA7[3] = { 0.1016969328379501162732e+02, 0.6623045922639275970621e+01, 0.3281013624325058830036e+01 };

/* T7: 7x7 row-major. Row 6 has implicit entries:
 *   ZZ(N6) = T771*Z1 + 1*Z2 + 0*Z3 + 1*Z4 + 0*Z5 + 1*Z6 + 0*Z7 */
static const sunrealtype T7[49] = {
  -0.2153754627310526422828e-02,  0.2156755135132077338691e-01,  0.8783567925144144407326e-02, -0.4055161452331023898198e-02,  0.4427232753268285479678e-02, -0.1238646187952874056377e-02, -0.2760617480543852499548e-02,
   0.1600025077880428526831e-02, -0.3813164813441154669442e-01, -0.2152556059400687552385e-01,  0.8415568276559589237177e-02, -0.4031949570224549492304e-02, -0.6666635339396338181761e-04,  0.3185474825166209848748e-02,
  -0.4059107301947683091650e-02,  0.5739650893938171539757e-01,  0.5885052920842679105612e-01, -0.8560431061603432060177e-02, -0.6923212665023908924141e-02, -0.2352180982943338340535e-02,  0.4169077725297562691409e-03,
  -0.1575048807937684420346e-01, -0.3821469359696835048464e-01, -0.1657368112729438512412e+00, -0.3737124230238445741907e-01,  0.8239007298507719404499e-02,  0.3115071152346175252726e-02,  0.2511660491343882192836e-01,
  -0.1129776610242208076086e+00, -0.2491742124652636863308e+00,  0.2735633057986623212132e+00,  0.5366761379181770094279e-02,  0.1932111161012620144312e+00,  0.1017177324817151468081e+00,  0.9504502035604622821039e-01,
  -0.4583810431839315010281e+00,  0.5315846490836284292051e+00,  0.4863228366175728940567e+00,  0.5265742264584492629141e+00,  0.2755343949896258141929e+00,  0.5217519452747652852946e+00,  0.1280719446355438944141e+00,
  -0.8813915783538183763135e+00,  1.0, 0.0, 1.0, 0.0, 1.0, 0.0
};
static const sunrealtype TI7[49] = {
  -0.2581319263199822292761e+03, -0.1890737630813985089520e+03, -0.4908731481793013119445e+02, -0.4110647469661428418112e+01, -0.4053447889315563304175e+01,  0.3112755366607346076554e+01, -0.1646774913558444650169e+01,
  -0.3007390169451292131731e+01, -0.1101586607876577132911e+02,  0.1487799456131656281486e+01,  0.2130388159559282459432e+01, -0.1816141086817565624822e+01,  0.1134325587895161100083e+01, -0.4146990459433035319930e+00,
  -0.8441963188321084681757e+01, -0.6505252740575150028169e+00,  0.6940670730369876478804e+01, -0.3205047525597898431565e+01,  0.1071280943546478589783e+01, -0.3548507491216221879730e+00,  0.9198549132786554154409e-01,
   0.7467833223502269977153e+02,  0.8740858897990081640204e+02,  0.4024158737379997877014e+01, -0.3714806315158364186639e+01, -0.3430093985982317350741e+01,  0.2696604809765312378853e+01, -0.9386927436075461933568e+00,
   0.5835652885190657724237e+02, -0.1006877395780018096325e+02, -0.3036638884256667120811e+02, -0.1020020865184865985027e+01, -0.1124175003784249621267e+00,  0.1890640831000377622800e+01, -0.9716486393831482282172e+00,
  -0.2991862480282520966786e+03, -0.2430407453687447911819e+03, -0.4877710407803786921219e+02, -0.2038671905741934405280e+01,  0.1673560239861084944268e+01, -0.1087374032057106164456e+01,  0.9019382492960993738427e+00,
  -0.9307650289743530591157e+02,  0.2388163105628114427703e+02,  0.3927888073081384382710e+02,  0.1438891568549108006988e+02, -0.3510438399399361221087e+01,  0.4863284885566180701215e+01, -0.2246482729591239916400e+01
};

/* Schur decomposition for ns=7 (user-provided) */
static const sunrealtype US7[49] = {
   8.26531181102153e-003,    19.3703898664289e-003,   -104.315969136059e-003,   -81.4241406430193e-003,   -181.772743196840e-003,   -559.840169567410e-003,   -797.234228373487e-003,
  -25.2895100166110e-003,   -35.2061557954432e-003,    204.229798631992e-003,    81.2060167635775e-003,    458.758681238129e-003,    631.157581991491e-003,   -583.949756873537e-003,
   85.1613253344128e-003,    57.2237885650667e-003,   -351.097510657236e-003,    147.755331790228e-003,   -751.240218165059e-003,    506.977266649751e-003,   -151.604987843608e-003,
  -282.035421030503e-003,   -58.3426730538428e-003,    272.480417751438e-003,   -875.341487409715e-003,   -224.492116364809e-003,    161.325654110381e-003,   -12.6956355624324e-003,
   572.286856307417e-003,   -163.800106226845e-003,    737.066906580986e-003,    123.716054659517e-003,   -292.125545934809e-003,   -39.8739638324410e-003,   -12.5189760381881e-003,
   684.163402705269e-003,    520.591065432539e-003,   -268.339696490930e-003,   -380.719854826641e-003,    202.315726412850e-003,    54.1064193617782e-003,    9.61366987103237e-003,
  -341.912105301555e-003,    832.983404940105e-003,    366.905315301544e-003,    196.132772552449e-003,   -124.384801238100e-003,   -25.4899598220441e-003,   -5.08602460080364e-003
};
static const sunrealtype TS7[49] = {
   4.37869356150681e+000,    4.75764227602788e+000,   -12.9046084298547e+000,   -11.2528647710127e+000,    9.09988536341528e+000,    3.62412061623971e+000,    1.47489567635798e+000,
  -21.7382172694186e+000,    4.37869356150681e+000,    17.2961366873499e+000,    20.3521232170655e+000,   -17.3116977259295e+000,   -8.50794959299232e+000,   -4.23541068770328e+000,
   0.0,    0.0,    7.14105521918778e+000,    16.5371636507187e+000,   -9.26428938178073e+000,   -2.24488757181613e+000,    158.402801655936e-003,
   0.0,    0.0,   -2.65249460063750e+000,    7.14105521918778e+000,   -13.0179588276324e+000,   -5.88900159332369e+000,   -2.51500945088783e+000,
   0.0,    0.0,    0.0,    0.0,    8.51183482510229e+000,    14.1545546297742e+000,    8.90382230904172e+000,
   0.0,    0.0,    0.0,    0.0,   -760.536144341915e-003,    8.51183482510229e+000,    14.1567371728117e+000,
   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.93683278840622e+000
};

/* ===========================================================================
 * radau5_InitConstants
 *
 * Initialize all Radau IIA method constants for the current ns.
 * ===========================================================================*/

/* ---------------------------------------------------------------------------
 * radau5_InitConstants
 *
 * Initialize all Radau IIA 3-stage order-5 method constants.
 * Transcribed exactly from Fortran radau5.f lines 706-737.
 * ---------------------------------------------------------------------------*/
int radau5_InitConstants(Radau5Mem rmem)
{
  int ns = rmem->ns;
  int npairs = rmem->npairs;

  if (ns == 3)
  {
    /* --- NS=3 (order 5) constants --- */
    sunrealtype SQ6 = sqrt(SUN_RCONST(6.0));
    sunrealtype c1  = (SUN_RCONST(4.0) - SQ6) / SUN_RCONST(10.0);
    sunrealtype c2  = (SUN_RCONST(4.0) + SQ6) / SUN_RCONST(10.0);

    rmem->c[0] = c1;
    rmem->c[1] = c2;
    rmem->c[2] = SUN_RCONST(1.0);

    rmem->dd[0] = -(SUN_RCONST(13.0) + SUN_RCONST(7.0) * SQ6) / SUN_RCONST(3.0);
    rmem->dd[1] = (-SUN_RCONST(13.0) + SUN_RCONST(7.0) * SQ6) / SUN_RCONST(3.0);
    rmem->dd[2] = -SUN_RCONST(1.0) / SUN_RCONST(3.0);

    sunrealtype u1 = (SUN_RCONST(6.0) + pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
                                        - pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
                     / SUN_RCONST(30.0);
    sunrealtype alph = (SUN_RCONST(12.0) - pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
                                          + pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
                       / SUN_RCONST(60.0);
    sunrealtype beta = (pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
                       + pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
                       * sqrt(SUN_RCONST(3.0)) / SUN_RCONST(60.0);
    sunrealtype cno = alph * alph + beta * beta;
    u1   = SUN_RCONST(1.0) / u1;
    alph = alph / cno;
    beta = beta / cno;

    rmem->u1 = u1;
    rmem->alph[0] = alph;
    rmem->beta_eig[0] = beta;

    /* Eigenvector matrix T (3x3 row-major) */
    rmem->T_mat[0] =  SUN_RCONST(9.1232394870892942792e-02);
    rmem->T_mat[1] = -SUN_RCONST(0.14125529502095420843);
    rmem->T_mat[2] = -SUN_RCONST(3.0029194105147424492e-02);
    rmem->T_mat[3] =  SUN_RCONST(0.24171793270710701896);
    rmem->T_mat[4] =  SUN_RCONST(0.20412935229379993199);
    rmem->T_mat[5] =  SUN_RCONST(0.38294211275726193779);
    rmem->T_mat[6] =  SUN_RCONST(0.96604818261509293619);
    rmem->T_mat[7] =  SUN_RCONST(1.0);  /* T32 */
    rmem->T_mat[8] =  SUN_RCONST(0.0);  /* T33 */

    /* Inverse eigenvector matrix TI (3x3 row-major) */
    rmem->TI_mat[0] =  SUN_RCONST(4.3255798900631553510);
    rmem->TI_mat[1] =  SUN_RCONST(0.33919925181580986954);
    rmem->TI_mat[2] =  SUN_RCONST(0.54177053993587487119);
    rmem->TI_mat[3] = -SUN_RCONST(4.1787185915519047273);
    rmem->TI_mat[4] = -SUN_RCONST(0.32768282076106238708);
    rmem->TI_mat[5] =  SUN_RCONST(0.47662355450055045196);
    rmem->TI_mat[6] = -SUN_RCONST(0.50287263494578687595);
    rmem->TI_mat[7] =  SUN_RCONST(2.5719269498556054292);
    rmem->TI_mat[8] = -SUN_RCONST(0.59603920482822492497);

    /* Schur constants for ns=3 */
    if (rmem->use_schur)
    {
      rmem->US_mat[0] =  SUN_RCONST( 0.138665108751908);
      rmem->US_mat[1] =  SUN_RCONST( 0.046278149309488);
      rmem->US_mat[2] =  SUN_RCONST( 0.989257459163847);
      rmem->US_mat[3] =  SUN_RCONST(-0.229641242351741);
      rmem->US_mat[4] =  SUN_RCONST(-0.970178886551833);
      rmem->US_mat[5] =  SUN_RCONST( 0.077574660168092);
      rmem->US_mat[6] =  SUN_RCONST(-0.963346711950568);
      rmem->US_mat[7] =  SUN_RCONST( 0.237931210616713);
      rmem->US_mat[8] =  SUN_RCONST( 0.123902589111344);

      rmem->TS_mat[0] =  SUN_RCONST( 2.68108287362775);
      rmem->TS_mat[1] =  SUN_RCONST(-8.42387579808538);
      rmem->TS_mat[2] =  SUN_RCONST(-4.08857781305964);
      rmem->TS_mat[3] =  SUN_RCONST( 1.10461319985220);
      rmem->TS_mat[4] =  SUN_RCONST( 2.68108287362775);
      rmem->TS_mat[5] =  SUN_RCONST( 4.70015200634394);
      rmem->TS_mat[6] =  SUN_RCONST( 0.0);
      rmem->TS_mat[7] =  SUN_RCONST( 0.0);
      rmem->TS_mat[8] =  SUN_RCONST( 3.63783425274450);

      rmem->u1 = rmem->TS_mat[8]; /* real eigenvalue from 1x1 block */

      /* T = US, TI = US^T */
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          rmem->T_mat[i*3+j] = rmem->US_mat[i*3+j];
          rmem->TI_mat[i*3+j] = rmem->US_mat[j*3+i]; /* transpose */
        }
    }
  }
  else if (ns == 5)
  {
    /* --- NS=5 (order 9) constants from precomputed tables --- */
    memcpy(rmem->c, C5, 5 * sizeof(sunrealtype));
    memcpy(rmem->dd, DD5, 5 * sizeof(sunrealtype));
    rmem->u1 = U1_5;
    memcpy(rmem->alph, ALPH5, 2 * sizeof(sunrealtype));
    memcpy(rmem->beta_eig, BETA5, 2 * sizeof(sunrealtype));
    memcpy(rmem->T_mat, T5, 25 * sizeof(sunrealtype));
    memcpy(rmem->TI_mat, TI5, 25 * sizeof(sunrealtype));

    if (rmem->use_schur)
    {
      memcpy(rmem->US_mat, US5, 25 * sizeof(sunrealtype));
      memcpy(rmem->TS_mat, TS5, 25 * sizeof(sunrealtype));
      rmem->u1 = rmem->TS_mat[4*5+4]; /* TS[4][4] = real eigenvalue */
      /* T = US, TI = US^T */
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) {
          rmem->T_mat[i*5+j] = rmem->US_mat[i*5+j];
          rmem->TI_mat[i*5+j] = rmem->US_mat[j*5+i];
        }
    }
  }
  else if (ns == 7)
  {
    /* --- NS=7 (order 13) constants from precomputed tables --- */
    memcpy(rmem->c, C7, 7 * sizeof(sunrealtype));
    memcpy(rmem->dd, DD7, 7 * sizeof(sunrealtype));
    rmem->u1 = U1_7;
    memcpy(rmem->alph, ALPH7, 3 * sizeof(sunrealtype));
    memcpy(rmem->beta_eig, BETA7, 3 * sizeof(sunrealtype));
    memcpy(rmem->T_mat, T7, 49 * sizeof(sunrealtype));
    memcpy(rmem->TI_mat, TI7, 49 * sizeof(sunrealtype));

    if (rmem->use_schur)
    {
      memcpy(rmem->US_mat, US7, 49 * sizeof(sunrealtype));
      memcpy(rmem->TS_mat, TS7, 49 * sizeof(sunrealtype));
      rmem->u1 = rmem->TS_mat[6*7+6]; /* TS[6][6] = real eigenvalue */
      /* T = US, TI = US^T */
      for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++) {
          rmem->T_mat[i*7+j] = rmem->US_mat[i*7+j];
          rmem->TI_mat[i*7+j] = rmem->US_mat[j*7+i];
        }
    }
  }
  else
  {
    return RADAU5_ILL_INPUT;
  }

  /* fnewt will be computed after tolerance transformation in Radau5Solve */
  rmem->fnewt = SUN_RCONST(0.0);

  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * radau5_ChangeOrder
 *
 * Runtime order change: updates ns, npairs, reloads method constants,
 * rebuilds E2 matrices and linear solvers for the new number of pairs.
 * N_Vectors (z, f, cont) are pre-allocated to RADAU5_NS_MAX so no
 * reallocation is needed.
 *
 * For dense/band matrices, E2 matrices are rebuilt from the J template.
 * For sparse matrices with variable order, this is not yet supported
 * (would require rebuilding CSC patterns).
 * ===========================================================================*/
int radau5_ChangeOrder(Radau5Mem rmem, int nsnew)
{
  int old_npairs = rmem->npairs;
  int new_npairs = (nsnew - 1) / 2;

  /* Update ns and npairs */
  rmem->ns = nsnew;
  rmem->npairs = new_npairs;

  /* Reload method constants for new ns */
  int ret = radau5_InitConstants(rmem);
  if (ret != RADAU5_SUCCESS) return ret;

  /* Update max Newton iterations */
  rmem->nit = 7 + (nsnew - 3) * 5 / 2;

  /* Recompute fnewt for new ns (Fortran: EXPMI=1/EXPMNS, FNEWT=...) */
  {
    sunrealtype uround = SUN_UNIT_ROUNDOFF;
    sunrealtype expmns = (sunrealtype)(nsnew + 1) / (SUN_RCONST(2.0) * (sunrealtype)nsnew);
    sunrealtype expmi  = SUN_RCONST(1.0) / expmns;

    sunrealtype rtol_orig = (rmem->itol == 0) ? rmem->rtol_s
                                               : N_VMin(rmem->rtol_v);
    if (rtol_orig <= SUN_RCONST(0.0)) rtol_orig = SUN_RCONST(1.0e-6);
    rmem->fnewt = SUNMAX(SUN_RCONST(10.0) * uround / rtol_orig,
                         SUNMIN(SUN_RCONST(0.03),
                                SUNRpowerR(rtol_orig, expmi - SUN_RCONST(1.0))));
  }

  /* E2/LS_E2 are pre-allocated to RADAU5_NPAIRS_MAX in SetLinearSolver,
   * so no reallocation is needed here. */

  /* Recompute error weights */
  radau5_ComputeScal(rmem, rmem->ycur);

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_SparseLookup
 *
 * Binary search for value at (row, col) in a CSC sparse matrix.
 * Returns 0.0 if the entry is not in the sparsity pattern.
 * ---------------------------------------------------------------------------*/
static sunrealtype radau5_SparseLookup(SUNMatrix A, sunindextype row,
                                       sunindextype col)
{
  sunindextype *Ap = SM_INDEXPTRS_S(A);
  sunindextype *Ai = SM_INDEXVALS_S(A);
  sunrealtype  *Ad = SM_DATA_S(A);
  sunindextype lo = Ap[col], hi = Ap[col + 1];
  while (lo < hi)
  {
    sunindextype mid = lo + (hi - lo) / 2;
    if (Ai[mid] == row) return Ad[mid];
    else if (Ai[mid] < row) lo = mid + 1;
    else hi = mid;
  }
  return SUN_RCONST(0.0);
}

/* ---------------------------------------------------------------------------
 * radau5_SparseUnion
 *
 * Create a new n×n CSC sparse matrix whose sparsity pattern is the union
 * of A's and B's patterns. Values are zeroed. Row indices are sorted.
 * ---------------------------------------------------------------------------*/
SUNMatrix radau5_SparseUnion(SUNMatrix A, SUNMatrix B, SUNContext sunctx)
{
  sunindextype n = SM_COLUMNS_S(A);
  sunindextype *Ap = SM_INDEXPTRS_S(A);
  sunindextype *Ai = SM_INDEXVALS_S(A);
  sunindextype *Bp = SM_INDEXPTRS_S(B);
  sunindextype *Bi = SM_INDEXVALS_S(B);

  /* First pass: count union NNZ per column */
  sunindextype total_nnz = 0;
  sunindextype j;
  for (j = 0; j < n; j++)
  {
    sunindextype ka = Ap[j], ka_end = Ap[j + 1];
    sunindextype kb = Bp[j], kb_end = Bp[j + 1];
    while (ka < ka_end && kb < kb_end)
    {
      if (Ai[ka] < Bi[kb])      { total_nnz++; ka++; }
      else if (Ai[ka] > Bi[kb]) { total_nnz++; kb++; }
      else                       { total_nnz++; ka++; kb++; }
    }
    total_nnz += (ka_end - ka) + (kb_end - kb);
  }

  SUNMatrix C = SUNSparseMatrix(n, n, total_nnz, CSC_MAT, sunctx);
  if (!C) return NULL;

  sunindextype *Cp = SM_INDEXPTRS_S(C);
  sunindextype *Ci = SM_INDEXVALS_S(C);

  /* Second pass: fill pattern */
  sunindextype pos = 0;
  for (j = 0; j < n; j++)
  {
    Cp[j] = pos;
    sunindextype ka = Ap[j], ka_end = Ap[j + 1];
    sunindextype kb = Bp[j], kb_end = Bp[j + 1];
    while (ka < ka_end && kb < kb_end)
    {
      if (Ai[ka] < Bi[kb])      { Ci[pos++] = Ai[ka++]; }
      else if (Ai[ka] > Bi[kb]) { Ci[pos++] = Bi[kb++]; }
      else                       { Ci[pos++] = Ai[ka++]; kb++; }
    }
    while (ka < ka_end) Ci[pos++] = Ai[ka++];
    while (kb < kb_end) Ci[pos++] = Bi[kb++];
  }
  Cp[n] = pos;

  return C;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacDense
 *
 * Column-by-column forward-difference Jacobian for dense matrices.
 * Mirrors CVODE's cvLsDenseDQJac approach.
 * ---------------------------------------------------------------------------*/
int radau5_DQJacDense(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype j, i, n;
  sunrealtype  srur, inc, inc_inv, ysave;
  sunrealtype *y_data, *fy_data, *tmp1_data;
  SUNMatrix    J;

  n         = rmem->n;
  J         = rmem->J;
  srur      = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fy_data   = N_VGetArrayPointer(fy);

  for (j = 0; j < n; j++)
  {
    y_data    = N_VGetArrayPointer(y);
    ysave     = y_data[j];

    /* Fortran radau5.f: DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE))) */
    inc = SUNRsqrt(SUN_UNIT_ROUNDOFF * SUNMAX(SUN_RCONST(1.0e-5), fabs(ysave)));

    y_data[j] = ysave + inc;
    {
      int rhsret = rmem->rhs(t, y, rmem->tmp1, rmem->user_data);
      if (rhsret < 0) { y_data[j] = ysave; return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) { y_data[j] = ysave; return RADAU5_RHSFUNC_RECVR; }
    }
    y_data[j] = ysave;

    inc_inv   = SUN_RCONST(1.0) / inc;
    tmp1_data = N_VGetArrayPointer(rmem->tmp1);

    for (i = 0; i < n; i++)
      SM_ELEMENT_D(J, i, j) = (tmp1_data[i] - fy_data[i]) * inc_inv;

    rmem->nfcn++;
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacBand
 *
 * Band finite-difference Jacobian using Curtis-Powell-Reid (CPR) grouped
 * column perturbation.  Mirrors CVODE's cvLsBandDQJac.
 * Total RHS evaluations: min(mu + ml + 1, n).
 * ---------------------------------------------------------------------------*/
int radau5_DQJacBand(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype i, j, n, group, i1, i2;
  sunindextype mu, ml, width, ngroups;
  sunrealtype  srur, fnorm, minInc, inc, inc_inv;
  sunrealtype *y_data, *fy_data, *ytemp_data, *ftemp_data, *scal_data;
  sunrealtype *col_j;
  sunrealtype  inc_arr[1024]; /* stack buffer for per-column inc values */
  sunrealtype *incs;
  SUNMatrix    J;

  n    = rmem->n;
  J    = rmem->J;
  mu   = rmem->mu;
  ml   = rmem->ml;
  srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fnorm  = N_VWrmsNorm(fy, rmem->scal);
  minInc = (fnorm != SUN_RCONST(0.0))
           ? (SUN_RCONST(1000.0) * fabs(rmem->h) * SUN_UNIT_ROUNDOFF
              * (sunrealtype)n * fnorm)
           : SUN_RCONST(1.0);

  y_data    = N_VGetArrayPointer(y);
  fy_data   = N_VGetArrayPointer(fy);
  scal_data = N_VGetArrayPointer(rmem->scal);

  /* Use tmp1 as ytemp (copy of y), tmp2 as ftemp */
  N_VScale(SUN_RCONST(1.0), y, rmem->tmp1);
  ytemp_data = N_VGetArrayPointer(rmem->tmp1);
  ftemp_data = N_VGetArrayPointer(rmem->tmp2);

  /* Allocate inc array if n > stack buffer size */
  incs = (n <= 1024) ? inc_arr : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));

  width   = ml + mu + 1;
  ngroups = SUNMIN(width, n);

  SUNMatZero(J);

  for (group = 0; group < ngroups; group++)
  {
    /* Perturb all columns in this group simultaneously */
    for (j = group; j < n; j += width)
    {
      inc = SUNMAX(srur * fabs(y_data[j]), minInc / scal_data[j]);
      inc = (y_data[j] + inc) - y_data[j];
      if (inc == SUN_RCONST(0.0)) inc = srur;
      incs[j]        = inc;
      ytemp_data[j] += inc;
    }

    /* One RHS evaluation for the whole group */
    {
      int rhsret = rmem->rhs(t, rmem->tmp1, rmem->tmp2, rmem->user_data);
      if (rhsret < 0) { if (n > 1024) free(incs); return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) {
        for (j = group; j < n; j += width)
          ytemp_data[j] = y_data[j];
        if (n > 1024) free(incs);
        return RADAU5_RHSFUNC_RECVR;
      }
    }
    rmem->nfcn++;

    /* Restore ytemp and fill Jacobian columns */
    for (j = group; j < n; j += width)
    {
      ytemp_data[j] = y_data[j];
      inc_inv = SUN_RCONST(1.0) / incs[j];
      col_j   = SUNBandMatrix_Column(J, j);

      i1 = SUNMAX(0, j - mu);
      i2 = SUNMIN(j + ml, n - 1);
      for (i = i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j, i, j) =
          (ftemp_data[i] - fy_data[i]) * inc_inv;
    }
  }

  if (n > 1024) free(incs);

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacSparse
 *
 * Sparse finite-difference Jacobian using column grouping (CPR technique).
 * Columns in the same group are structurally independent (no shared nonzero
 * rows), so they can be perturbed simultaneously in a single RHS evaluation.
 * Total RHS evaluations: ngroups (typically << n for sparse patterns).
 *
 * Requires that Radau5SetSparsityPattern has been called to compute the
 * column grouping (col_group, group_offsets, group_cols).
 * ---------------------------------------------------------------------------*/
int radau5_DQJacSparse(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype j, k, p, n;
  sunrealtype  srur, fnorm, minInc, inc, inc_inv;
  sunrealtype *y_data, *fy_data, *ytemp_data, *ftemp_data, *scal_data;
  sunrealtype  inc_arr[1024];
  sunrealtype *incs;

  n    = rmem->n;
  srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fnorm  = N_VWrmsNorm(fy, rmem->scal);
  minInc = (fnorm != SUN_RCONST(0.0))
           ? (SUN_RCONST(1000.0) * fabs(rmem->h) * SUN_UNIT_ROUNDOFF
              * (sunrealtype)n * fnorm)
           : SUN_RCONST(1.0);

  y_data    = N_VGetArrayPointer(y);
  fy_data   = N_VGetArrayPointer(fy);
  scal_data = N_VGetArrayPointer(rmem->scal);

  /* Use tmp1 as ytemp (copy of y), tmp2 as ftemp */
  N_VScale(SUN_RCONST(1.0), y, rmem->tmp1);
  ytemp_data = N_VGetArrayPointer(rmem->tmp1);
  ftemp_data = N_VGetArrayPointer(rmem->tmp2);

  /* Per-column increment storage */
  incs = (n <= 1024) ? inc_arr : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));

  /* Zero J data array (preserve CSC structure) */
  sunrealtype* J_data = SM_DATA_S(rmem->J);
  sunindextype* Jp = SM_INDEXPTRS_S(rmem->J);
  sunindextype nnz_J = Jp[n];
  for (p = 0; p < nnz_J; p++) J_data[p] = SUN_RCONST(0.0);

  /* Sparsity pattern and group lookup from SetSparsityPattern */
  const sunindextype* sp_colptrs    = rmem->sp_colptrs;
  const sunindextype* sp_rowinds    = rmem->sp_rowinds;
  const sunindextype* group_offsets = rmem->group_offsets;
  const sunindextype* group_cols    = rmem->group_cols;

  for (sunindextype group = 0; group < rmem->ngroups; group++)
  {
    /* Perturb all columns in this group simultaneously */
    for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
    {
      j = group_cols[k];
      inc = SUNMAX(srur * fabs(y_data[j]), minInc / scal_data[j]);
      inc = (y_data[j] + inc) - y_data[j];
      if (inc == SUN_RCONST(0.0)) inc = srur;
      incs[j]        = inc;
      ytemp_data[j] += inc;
    }

    /* One RHS evaluation for the whole group */
    {
      int rhsret = rmem->rhs(t, rmem->tmp1, rmem->tmp2, rmem->user_data);
      if (rhsret < 0) { if (n > 1024) free(incs); return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) {
        for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
          ytemp_data[group_cols[k]] = y_data[group_cols[k]];
        if (n > 1024) free(incs);
        return RADAU5_RHSFUNC_RECVR;
      }
    }
    rmem->nfcn++;

    /* Extract Jacobian columns and restore ytemp */
    for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
    {
      j = group_cols[k];
      ytemp_data[j] = y_data[j];
      inc_inv = SUN_RCONST(1.0) / incs[j];

      for (p = sp_colptrs[j]; p < sp_colptrs[j + 1]; p++)
      {
        sunindextype row = sp_rowinds[p];
        J_data[p] = (ftemp_data[row] - fy_data[row]) * inc_inv;
      }
    }
  }

  if (n > 1024) free(incs);

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_BuildE1
 *
 * Build E1 = fac1*M - J  (or fac1*I - J when M is NULL).
 * Dispatches on matrix type: dense, band, or sparse.
 * ---------------------------------------------------------------------------*/
int radau5_BuildE1(Radau5Mem rmem, sunrealtype fac1)
{
  sunindextype i, j, k, n;
  SUNMatrix    J, M, E1;

  n  = rmem->n;
  J  = rmem->J;
  M  = rmem->M;
  E1 = rmem->E1;

  if (rmem->mat_id == SUNMATRIX_SPARSE)
  {
    /* --- Sparse (CSC) case ---
     * When M is sparse, E1 has the union(J,M) sparsity pattern.
     * When M is NULL/dense/band, E1 has J's pattern (cloned).
     * Fill: E1[i,j] = fac1*M[i,j] - J[i,j]  (or fac1*delta(i,j) - J[i,j]) */
    sunindextype *E1p = SM_INDEXPTRS_S(E1);
    sunindextype *E1i = SM_INDEXVALS_S(E1);
    sunrealtype  *E1d = SM_DATA_S(E1);

    if (M != NULL && SUNMatGetID(M) == SUNMATRIX_SPARSE)
    {
      /* E1 pattern = union(J, M). Look up J and M values by binary search. */
      for (j = 0; j < n; j++)
      {
        for (k = E1p[j]; k < E1p[j+1]; k++)
        {
          i = E1i[k];
          sunrealtype jij = radau5_SparseLookup(J, i, j);
          sunrealtype mij = radau5_SparseLookup(M, i, j);
          E1d[k] = fac1 * mij - jij;
        }
      }
    }
    else
    {
      /* E1 pattern = J's pattern. Iterate over J directly. */
      sunindextype *Jp = SM_INDEXPTRS_S(J);
      sunindextype *Ji = SM_INDEXVALS_S(J);
      sunrealtype  *Jd = SM_DATA_S(J);

      for (j = 0; j < n; j++)
      {
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          E1d[k] = -Jd[k];
          if (M == NULL)
          {
            if (i == j) E1d[k] += fac1;
          }
          else if (SUNMatGetID(M) == SUNMATRIX_DENSE)
          {
            E1d[k] += fac1 * SM_ELEMENT_D(M, i, j);
          }
          else if (SUNMatGetID(M) == SUNMATRIX_BAND)
          {
            sunindextype muM = SM_UBAND_B(M);
            sunindextype mlM = SM_LBAND_B(M);
            if (i >= j - muM && i <= j + mlM)
              E1d[k] += fac1 * SM_ELEMENT_B(M, i, j);
          }
        }
      }
    }
  }
  else if (rmem->mat_id == SUNMATRIX_BAND)
  {
    /* --- Band case --- */
    sunindextype mu = rmem->mu;
    sunindextype ml = rmem->ml;
    sunindextype i1, i2;

    SUNMatZero(E1);

    for (j = 0; j < n; j++)
    {
      i1 = SUNMAX(0, j - mu);
      i2 = SUNMIN(j + ml, n - 1);
      for (i = i1; i <= i2; i++)
        SM_ELEMENT_B(E1, i, j) = -SM_ELEMENT_B(J, i, j);

      if (M == NULL)
      {
        SM_ELEMENT_B(E1, j, j) += fac1;
      }
      else
      {
        /* M may be band or dense — check its type */
        if (SUNMatGetID(M) == SUNMATRIX_BAND)
        {
          sunindextype muM = SM_UBAND_B(M);
          sunindextype mlM = SM_LBAND_B(M);
          sunindextype im1 = SUNMAX(0, j - muM);
          sunindextype im2 = SUNMIN(j + mlM, n - 1);
          for (i = im1; i <= im2; i++)
            SM_ELEMENT_B(E1, i, j) += fac1 * SM_ELEMENT_B(M, i, j);
        }
        else
        {
          /* Dense M with band E1 — only copy within band */
          for (i = i1; i <= i2; i++)
            SM_ELEMENT_B(E1, i, j) += fac1 * SM_ELEMENT_D(M, i, j);
        }
      }
    }
  }
  else
  {
    /* --- Dense case --- */
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < n; i++)
        SM_ELEMENT_D(E1, i, j) = -SM_ELEMENT_D(J, i, j);

      if (M == NULL)
      {
        SM_ELEMENT_D(E1, j, j) += fac1;
      }
      else
      {
        for (i = 0; i < n; i++)
          SM_ELEMENT_D(E1, i, j) += fac1 * SM_ELEMENT_D(M, i, j);
      }
    }
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_BuildE2
 *
 * Build the 2n×2n system for the coupled block.
 *
 * Eigenvalue mode (use_schur==0):
 *   E2 = [ alphn*M - J,   -betan*M ]
 *        [ betan*M,     alphn*M - J ]
 *
 * Schur mode (use_schur==1):
 *   E2 = [ TS[0][0]/h*M - J,   TS[0][1]/h*M       ]
 *        [ TS[1][0]/h*M,        TS[1][1]/h*M - J   ]
 *
 * (or with I in place of M when M is NULL)
 * J may be dense, band, or sparse. E2 is dense for dense/band J,
 * sparse (CSC) for sparse J.
 * ---------------------------------------------------------------------------*/
int radau5_BuildE2(Radau5Mem rmem, int pair_idx, sunrealtype alphn, sunrealtype betan)
{
  sunindextype i, j, k, n;
  sunrealtype  jij, mij;
  SUNMatrix    J, M, E2;

  /* For Schur mode, the 2×2 block has 4 independent coefficients.
   * For eigenvalue mode: a00 = a11 = alphn, a01 = -betan, a10 = betan. */
  sunrealtype a00, a01, a10, a11;

  n  = rmem->n;
  J  = rmem->J;
  M  = rmem->M;
  E2 = rmem->E2[pair_idx];

  if (rmem->use_schur)
  {
    sunrealtype h = rmem->h;
    int ns = rmem->ns;
    int r0 = 2 * pair_idx;  /* row/col offset of this 2x2 diagonal block */
    a00 = rmem->TS_mat[r0 * ns + r0]       / h;
    a01 = rmem->TS_mat[r0 * ns + (r0 + 1)] / h;
    a10 = rmem->TS_mat[(r0 + 1) * ns + r0] / h;
    a11 = rmem->TS_mat[(r0 + 1) * ns + (r0 + 1)] / h;
  }
  else
  {
    a00 = alphn;
    a01 = -betan;
    a10 = betan;
    a11 = alphn;
  }

  if (rmem->mat_id == SUNMATRIX_SPARSE)
  {
    /* --- Sparse (CSC) case ---
     * Build the 2n×2n CSC matrix.
     * Column layout of E2 (2n columns):
     *   col j     (j=0..n-1): top-left block + bottom-left block
     *   col j+n   (j=0..n-1): top-right block + bottom-right block
     *
     * When M is sparse, diagonal blocks use union(J,M) pattern and
     * off-diagonal blocks use M's pattern. When M is NULL or dense,
     * the original approach is used.
     */
    sunindextype *Jp = SM_INDEXPTRS_S(J);
    sunindextype *Ji = SM_INDEXVALS_S(J);
    sunrealtype  *Jd = SM_DATA_S(J);

    sunindextype *E2p = SM_INDEXPTRS_S(E2);
    sunindextype *E2i = SM_INDEXVALS_S(E2);
    sunrealtype  *E2d = SM_DATA_S(E2);

    sunindextype pos = 0;
    int m_is_sparse = (M != NULL && SUNMatGetID(M) == SUNMATRIX_SPARSE);

    if (m_is_sparse)
    {
      /* E1 already holds the union(J,M) pattern from SetLinearSolver */
      SUNMatrix E1 = rmem->E1;
      sunindextype *Up = SM_INDEXPTRS_S(E1);  /* union pattern */
      sunindextype *Ui = SM_INDEXVALS_S(E1);
      sunindextype *Mp = SM_INDEXPTRS_S(M);
      sunindextype *Mi = SM_INDEXVALS_S(M);

      for (j = 0; j < n; j++)
      {
        E2p[j] = pos;

        /* Top-left block: a00*M[i,j] - J[i,j] over union pattern */
        for (k = Up[j]; k < Up[j+1]; k++)
        {
          i = Ui[k];
          E2i[pos] = i;
          E2d[pos] = a00 * radau5_SparseLookup(M, i, j)
                   - radau5_SparseLookup(J, i, j);
          pos++;
        }

        /* Bottom-left block: a10*M[i,j] over M's pattern */
        for (k = Mp[j]; k < Mp[j+1]; k++)
        {
          E2i[pos] = Mi[k] + n;
          E2d[pos] = a10 * SM_DATA_S(M)[k];
          pos++;
        }
      }

      for (j = 0; j < n; j++)
      {
        E2p[j + n] = pos;

        /* Top-right block: a01*M[i,j] over M's pattern */
        for (k = Mp[j]; k < Mp[j+1]; k++)
        {
          E2i[pos] = Mi[k];
          E2d[pos] = a01 * SM_DATA_S(M)[k];
          pos++;
        }

        /* Bottom-right block: a11*M[i,j] - J[i,j] over union pattern */
        for (k = Up[j]; k < Up[j+1]; k++)
        {
          i = Ui[k];
          E2i[pos] = i + n;
          E2d[pos] = a11 * radau5_SparseLookup(M, i, j)
                   - radau5_SparseLookup(J, i, j);
          pos++;
        }
      }

      E2p[2 * n] = pos;
    }
    else
    {
      /* M is NULL or dense — original approach */
      for (j = 0; j < n; j++)
      {
        E2p[j] = pos;

        /* Top-left block: a00*M[i,j] - J[i,j] for each (i,j) in J's pattern */
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          mij = (M == NULL) ? ((i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0))
                            : SM_ELEMENT_D(M, i, j);
          E2i[pos] = i;
          E2d[pos] = a00 * mij - Jd[k];
          pos++;
        }

        /* Bottom-left block: a10*M[i,j] */
        if (M == NULL)
        {
          E2i[pos] = j + n;
          E2d[pos] = a10;
          pos++;
        }
        else
        {
          for (i = 0; i < n; i++)
          {
            mij = SM_ELEMENT_D(M, i, j);
            if (mij != SUN_RCONST(0.0))
            {
              E2i[pos] = i + n;
              E2d[pos] = a10 * mij;
              pos++;
            }
          }
        }
      }

      for (j = 0; j < n; j++)
      {
        E2p[j + n] = pos;

        /* Top-right block: a01*M[i,j] */
        if (M == NULL)
        {
          E2i[pos] = j;
          E2d[pos] = a01;
          pos++;
        }
        else
        {
          for (i = 0; i < n; i++)
          {
            mij = SM_ELEMENT_D(M, i, j);
            if (mij != SUN_RCONST(0.0))
            {
              E2i[pos] = i;
              E2d[pos] = a01 * mij;
              pos++;
            }
          }
        }

        /* Bottom-right block: a11*M[i,j] - J[i,j] */
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          mij = (M == NULL) ? ((i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0))
                            : SM_ELEMENT_D(M, i, j);
          E2i[pos] = i + n;
          E2d[pos] = a11 * mij - Jd[k];
          pos++;
        }
      }

      E2p[2 * n] = pos;
    }
  }
  else
  {
    /* --- Dense / Band case (E2 is always dense) --- */
    int j_is_band = (rmem->mat_id == SUNMATRIX_BAND);
    sunindextype mu = rmem->mu;
    sunindextype ml = rmem->ml;

    for (j = 0; j < n; j++)
    {
      for (i = 0; i < n; i++)
      {
        /* Read J(i,j) */
        if (j_is_band)
          jij = (i >= j - mu && i <= j + ml) ? SM_ELEMENT_B(J, i, j)
                                              : SUN_RCONST(0.0);
        else
          jij = SM_ELEMENT_D(J, i, j);

        /* Read M(i,j) or identity */
        if (M == NULL)
        {
          mij = (i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
        }
        else if (SUNMatGetID(M) == SUNMATRIX_BAND)
        {
          sunindextype muM = SM_UBAND_B(M);
          sunindextype mlM = SM_LBAND_B(M);
          mij = (i >= j - muM && i <= j + mlM) ? SM_ELEMENT_B(M, i, j)
                                                 : SUN_RCONST(0.0);
        }
        else
        {
          mij = SM_ELEMENT_D(M, i, j);
        }

        SM_ELEMENT_D(E2, i,     j    ) = a00 * mij - jij;
        SM_ELEMENT_D(E2, i,     j + n) = a01 * mij;
        SM_ELEMENT_D(E2, i + n, j    ) = a10 * mij;
        SM_ELEMENT_D(E2, i + n, j + n) = a11 * mij - jij;
      }
    }
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DecompE1
 *
 * Factor E1 via SUNLinSolSetup. Increments ndec.
 * ---------------------------------------------------------------------------*/
int radau5_DecompE1(Radau5Mem rmem)
{
  int retval;

  retval = SUNLinSolSetup(rmem->LS_E1, rmem->E1);
  rmem->ndec++;

  return (retval == 0) ? RADAU5_SUCCESS : RADAU5_SINGULAR_MATRIX;
}

/* ---------------------------------------------------------------------------
 * radau5_DecompE2
 *
 * Factor E2 via SUNLinSolSetup. Increments ndec.
 * ---------------------------------------------------------------------------*/
int radau5_DecompE2(Radau5Mem rmem, int pair_idx)
{
  int retval;

  retval = SUNLinSolSetup(rmem->LS_E2[pair_idx], rmem->E2[pair_idx]);
  rmem->ndec++;

  return (retval == 0) ? RADAU5_SUCCESS : RADAU5_SINGULAR_MATRIX;
}

/* ---------------------------------------------------------------------------
 * radau5_ComputeScal
 *
 * Compute scal[i] = atol[i] + rtol[i] * |y[i]|
 * Matches Fortran radau5.f lines 780-788.
 * ---------------------------------------------------------------------------*/
void radau5_ComputeScal(Radau5Mem rmem, N_Vector y)
{
  sunindextype  i, n;
  sunrealtype  *y_data, *scal_data;
  sunrealtype  *rtol_v_data, *atol_v_data;

  n         = rmem->n;
  y_data    = N_VGetArrayPointer(y);
  scal_data = N_VGetArrayPointer(rmem->scal);

  if (rmem->itol == 0)
  {
    /* Scalar tolerances */
    for (i = 0; i < n; i++)
      scal_data[i] = rmem->atol_s + rmem->rtol_s * fabs(y_data[i]);
  }
  else
  {
    /* Vector tolerances */
    rtol_v_data = N_VGetArrayPointer(rmem->rtol_v);
    atol_v_data = N_VGetArrayPointer(rmem->atol_v);
    for (i = 0; i < n; i++)
      scal_data[i] = atol_v_data[i] + rtol_v_data[i] * fabs(y_data[i]);
  }
}

/* ---------------------------------------------------------------------------
 * radau5_MassMult
 *
 * Compute result = M * x, dispatching on M's matrix type (dense or band).
 * ---------------------------------------------------------------------------*/
int radau5_MassMult(Radau5Mem rmem, N_Vector x, N_Vector result)
{
  sunindextype i, j, n = rmem->n;
  SUNMatrix M = rmem->M;
  sunrealtype *xd = N_VGetArrayPointer(x);
  sunrealtype *rd = N_VGetArrayPointer(result);

  if (M == NULL)
  {
    /* Identity mass — just copy */
    N_VScale(SUN_RCONST(1.0), x, result);
    return RADAU5_SUCCESS;
  }

  SUNMatrix_ID mid = SUNMatGetID(M);

  if (mid == SUNMATRIX_DENSE)
  {
    for (i = 0; i < n; i++)
    {
      sunrealtype sum = SUN_RCONST(0.0);
      for (j = 0; j < n; j++)
        sum += SM_ELEMENT_D(M, i, j) * xd[j];
      rd[i] = sum;
    }
  }
  else if (mid == SUNMATRIX_BAND)
  {
    sunindextype muM = SM_UBAND_B(M);
    sunindextype mlM = SM_LBAND_B(M);
    for (i = 0; i < n; i++)
    {
      sunrealtype sum = SUN_RCONST(0.0);
      sunindextype j1 = SUNMAX(0, i - mlM);
      sunindextype j2 = SUNMIN(n - 1, i + muM);
      for (j = j1; j <= j2; j++)
        sum += SM_ELEMENT_B(M, i, j) * xd[j];
      rd[i] = sum;
    }
  }
  else
  {
    /* Fallback: use SUNMatMatvec if available */
    SUNErrCode err = SUNMatMatvec(M, x, result);
    if (err != SUN_SUCCESS) return RADAU5_LSOLVE_FAIL;
  }

  return RADAU5_SUCCESS;
}
