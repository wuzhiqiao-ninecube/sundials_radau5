function tba_ode15s
%TBA_ODE15S Two Bit Adding Unit (DAE index-1, n=350) solved with ode15s.
%
%  MOSFET circuit simulation of a two-bit adder computing:
%    A1*2+A0 + B1*2+B0 + CIN = C*4 + S1*2 + S0
%
%  State vector (1-based):
%    y(1:175)   = charges g(U)
%    y(176:350) = node voltages U
%
%  M*y' = f(t,y) where M = diag(1,...,1,0,...,0) (175 ones, 175 zeros)
%  t in [0, 320] with 63 discontinuities at t=5*k, k=1..63
%
%  Reference solution at t=315:
%    y(224) ~ 0.2040  (S0: node 49)
%    y(305) ~ 4.9972  (S1: node 130)
%    y(323) ~ 0.2039  (C:  node 148)
%
%  Based on IVPtestset 2.4 tba.f and radau5_tba.c

clc; close all;
fprintf('Two Bit Adding Unit (n=350, DAE-1) with ode15s\n');

%% Parameters
p = init_params();

%% Mass matrix (constant, singular): 175 ones then 175 zeros
M = diag([ones(175,1); zeros(175,1)]);

%% Initial conditions
U0 = init_voltages();
G0 = GCN(U0, p);
y0 = [G0; U0];

%% Tolerances
Tol = 1e-5;
rtol = Tol;
atol = Tol;

%% Discontinuity times
disc = 0:5:320;

%% ODE options
opts = odeset('Mass', M, ...
              'MassSingular', 'yes', ...
              'RelTol', rtol, ...
              'AbsTol', atol, ...
              'InitialStep', 4e-5, ...
              'MaxStep', 0.5, ...
              'Jacobian', @(t,y) jac_tba(t, y, p));

%% Segmented integration
y = y0;
total_steps = 0;
tout_all = [];
yout_all = [];

tic;
for seg = 1:length(disc)-1
    tspan = [disc(seg), disc(seg+1)];
    if tspan(2) <= tspan(1), continue; end

    [tsol, ysol] = ode15s(@(t,y) rhs_tba(t, y, p), tspan, y, opts);

    y = ysol(end,:)';
    total_steps = total_steps + length(tsol) - 1;
    tout_all = [tout_all; tsol];
    yout_all = [yout_all; ysol];
end
elapsed = toc;

%% Reference solution at t=315
yref_224 = 0.2040419147264534;   % S0: node 49
yref_305 = 0.4997238455712048e1; % S1: node 130
yref_323 = 0.2038985905095614;   % C:  node 148

fprintf('\n=== Results at t=%.1f (elapsed %.2f s, steps=%d) ===\n', ...
        tout_all(end), elapsed, total_steps);
fprintf('y(224) (S0) = %16.10e  ref = %16.10e  rel_err = %.3e\n', ...
        y(224), yref_224, abs(y(224)-yref_224)/abs(yref_224));
fprintf('y(305) (S1) = %16.10e  ref = %16.10e  rel_err = %.3e\n', ...
        y(305), yref_305, abs(y(305)-yref_305)/abs(yref_305));
fprintf('y(323) (C)  = %16.10e  ref = %16.10e  rel_err = %.3e\n', ...
        y(323), yref_323, abs(y(323)-yref_323)/abs(yref_323));

maxerr = max([abs(y(224)-yref_224)/abs(yref_224), ...
              abs(y(305)-yref_305)/abs(yref_305), ...
              abs(y(323)-yref_323)/abs(yref_323)]);
fprintf('max_rel_err = %.3e\n', maxerr);
if maxerr < 0.1
    fprintf('PASSED\n');
else
    fprintf('FAILED\n');
end

%% Plot output signals
figure(1);
subplot(3,1,1); plot(tout_all, yout_all(:,224)); title('S0 (node 49)'); ylabel('V');
subplot(3,1,2); plot(tout_all, yout_all(:,305)); title('S1 (node 130)'); ylabel('V');
subplot(3,1,3); plot(tout_all, yout_all(:,323)); title('C (node 148)'); ylabel('V'); xlabel('t');

end % tba_ode15s

%% =========================================================================
%  Local functions
%  =========================================================================

%% --- RHS function --------------------------------------------------------
function f = rhs_tba(t, y, p)
    N = 175;
    U = y(N+1:2*N);  % voltages

    % Circuit current equations
    res = FCN(N, t, U, p);

    % Charge function
    G = GCN(U, p);

    f = zeros(2*N, 1);
    f(1:N) = res;              % M*y'(1:175) = FCN(t, U)
    f(N+1:2*N) = y(1:N) - G;  % 0 = y(1:175) - GCN(U)
end

%% --- Numerical Jacobian (sparse, via finite differences) -----------------
function J = jac_tba(t, y, p)
    n = length(y);
    f0 = rhs_tba(t, y, p);
    J = zeros(n, n);
    sqrteps = sqrt(eps);
    for j = 1:n
        yp = y;
        hj = sqrteps * max(abs(y(j)), 1e-5);
        yp(j) = yp(j) + hj;
        fp = rhs_tba(t, yp, p);
        J(:,j) = (fp - f0) / hj;
    end
    J = sparse(J);
end

%% --- MOSFET model functions ----------------------------------------------
function c = CBDBS(V, p)
    PHIB = 0.87;
    if V <= 0
        c = p.CBD / sqrt(1 - V/PHIB);
    else
        c = p.CBD * (1 + V/(2*PHIB));
    end
end

function i = fn_IBS(VBS, p)
    if VBS <= 0
        i = -p.CURIS * (exp(VBS/p.VTH) - 1);
    else
        i = 0;
    end
end

function i = fn_IBD(VBD, p)
    if VBD <= 0
        i = -p.CURIS * (exp(VBD/p.VTH) - 1);
    else
        i = 0;
    end
end

function val = GDSP(NED, VDS, VGS, VBS, p)
    [VT0, CGAMMA, PHI, BETA] = mosfet_params(NED, p);
    if PHI - VBS < 0 || PHI < 0
        error('GDSP: illegal input');
    end
    VTE = VT0 + CGAMMA*(sqrt(PHI-VBS) - sqrt(PHI));
    if VGS - VTE <= 0
        val = 0;
    elseif VGS - VTE <= VDS
        val = -BETA*(VGS-VTE)^2*(1 + p.DELTA*VDS);
    else
        val = -BETA*VDS*(2*(VGS-VTE)-VDS)*(1 + p.DELTA*VDS);
    end
end

function val = GDSM(NED, VDS, VGD, VBD, p)
    [VT0, CGAMMA, PHI, BETA] = mosfet_params(NED, p);
    if PHI - VBD < 0 || PHI < 0
        error('GDSM: illegal input');
    end
    VTE = VT0 + CGAMMA*(sqrt(PHI-VBD) - sqrt(PHI));
    if VGD - VTE <= 0
        val = 0;
    elseif VGD - VTE <= -VDS
        val = BETA*(VGD-VTE)^2*(1 - p.DELTA*VDS);
    else
        val = -BETA*VDS*(2*(VGD-VTE)+VDS)*(1 - p.DELTA*VDS);
    end
end

function val = fn_IDS(NED, VDS, VGS, VBS, VGD, VBD, p)
    if VDS > 0
        val = GDSP(NED, VDS, VGS, VBS, p);
    elseif VDS == 0
        val = 0;
    else
        val = GDSM(NED, VDS, VGD, VBD, p);
    end
end

function [VT0, CGAMMA, PHI, BETA] = mosfet_params(NED, p)
    if NED == 0
        VT0 = -2.43; CGAMMA = 0.2; PHI = 1.28;
        BETA = 53.5e-6 * p.CTIME * p.STIFF;
    elseif NED == 1
        VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
        BETA = 4*43.7e-6 * p.CTIME * p.STIFF;
    elseif NED == 2
        VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
        BETA = 8*43.7e-6 * p.CTIME * p.STIFF;
    else
        VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
        BETA = 12*43.7e-6 * p.CTIME * p.STIFF;
    end
end

%% --- PULSE input signal --------------------------------------------------
function [VIN, VIND] = PULSE(X, LOW, HIGH, DELAY, T1, T2, T3, PERIOD)
    TIME = mod(X, PERIOD);
    if TIME > (DELAY + T1 + T2 + T3)
        VIN = LOW; VIND = 0;
    elseif TIME > (DELAY + T1 + T2)
        VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW;
        VIND = -((HIGH-LOW)/T3);
    elseif TIME > (DELAY + T1)
        VIN = HIGH; VIND = 0;
    elseif TIME > DELAY
        VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW;
        VIND = ((HIGH-LOW)/T1);
    else
        VIN = LOW; VIND = 0;
    end
end

%% --- Gate subroutines (1-based indexing) ---------------------------------
function F = gate_NOR(N, I, U1, U2, U1D, U2D, Y, F, p)
    F(I) = -(Y(I)-Y(I+4))/p.RGS ...
           - fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+1) = -(Y(I+1)-p.VDD)/p.RGD ...
             + fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+2) = -(Y(I+2)-p.VBB)/p.RBS + fn_IBS(Y(I+2)-Y(I+4), p);
    F(I+3) = -(Y(I+3)-p.VBB)/p.RBD + fn_IBD(Y(I+3)-p.VDD, p);

    F(I+4) = -(Y(I+4)-Y(I))/p.RGS - fn_IBS(Y(I+2)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+6))/p.RGD - fn_IBD(Y(I+8)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+10))/p.RGD - fn_IBD(Y(I+12)-Y(I+4), p);

    F(I+5) = p.CGS*U1D - Y(I+5)/p.RGS ...
             - fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+6) = p.CGD*U1D - (Y(I+6)-Y(I+4))/p.RGD ...
             + fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+7) = -(Y(I+7)-p.VBB)/p.RBS + fn_IBS(Y(I+7), p);
    F(I+8) = -(Y(I+8)-p.VBB)/p.RBD + fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+9) = p.CGS*U2D - Y(I+9)/p.RGS ...
             - fn_IDS(1, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+10) = p.CGD*U2D - (Y(I+10)-Y(I+4))/p.RGD ...
              + fn_IDS(1, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+11) = -(Y(I+11)-p.VBB)/p.RBS + fn_IBS(Y(I+11), p);
    F(I+12) = -(Y(I+12)-p.VBB)/p.RBD + fn_IBD(Y(I+12)-Y(I+4), p);
end

function F = gate_ANDOI(N, I, U1, U2, U3, U1D, U2D, U3D, Y, F, p)
    F(I) = -(Y(I)-Y(I+4))/p.RGS ...
           - fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+1) = -(Y(I+1)-p.VDD)/p.RGD ...
             + fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+2) = -(Y(I+2)-p.VBB)/p.RBS + fn_IBS(Y(I+2)-Y(I+4), p);
    F(I+3) = -(Y(I+3)-p.VBB)/p.RBD + fn_IBD(Y(I+3)-p.VDD, p);

    F(I+4) = -(Y(I+4)-Y(I))/p.RGS - fn_IBS(Y(I+2)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+6))/p.RGD - fn_IBD(Y(I+8)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+10))/p.RGD - fn_IBD(Y(I+12)-Y(I+4), p);

    F(I+5) = p.CGS*U1D - Y(I+5)/p.RGS ...
             - fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+6) = p.CGD*U1D - (Y(I+6)-Y(I+4))/p.RGD ...
             + fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+7) = -(Y(I+7)-p.VBB)/p.RBS + fn_IBS(Y(I+7), p);
    F(I+8) = -(Y(I+8)-p.VBB)/p.RBD + fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+9) = p.CGS*U2D - (Y(I+9)-Y(I+13))/p.RGS ...
             - fn_IDS(2, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11)-Y(I+13), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+10) = p.CGD*U2D - (Y(I+10)-Y(I+4))/p.RGD ...
              + fn_IDS(2, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11)-Y(I+13), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+11) = -(Y(I+11)-p.VBB)/p.RBS + fn_IBS(Y(I+11)-Y(I+13), p);
    F(I+12) = -(Y(I+12)-p.VBB)/p.RBD + fn_IBD(Y(I+12)-Y(I+4), p);

    F(I+13) = -(Y(I+13)-Y(I+9))/p.RGS - fn_IBS(Y(I+11)-Y(I+13), p) ...
              -(Y(I+13)-Y(I+15))/p.RGD - fn_IBD(Y(I+17)-Y(I+13), p);

    F(I+14) = p.CGS*U3D - Y(I+14)/p.RGS ...
              - fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+13), p);
    F(I+15) = p.CGD*U3D - (Y(I+15)-Y(I+13))/p.RGD ...
              + fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+13), p);
    F(I+16) = -(Y(I+16)-p.VBB)/p.RBS + fn_IBS(Y(I+16), p);
    F(I+17) = -(Y(I+17)-p.VBB)/p.RBD + fn_IBD(Y(I+17)-Y(I+13), p);
end

function F = gate_NAND(N, I, U1, U2, U1D, U2D, Y, F, p)
    F(I) = -(Y(I)-Y(I+4))/p.RGS ...
           - fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+1) = -(Y(I+1)-p.VDD)/p.RGD ...
             + fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+2) = -(Y(I+2)-p.VBB)/p.RBS + fn_IBS(Y(I+2)-Y(I+4), p);
    F(I+3) = -(Y(I+3)-p.VBB)/p.RBD + fn_IBD(Y(I+3)-p.VDD, p);

    F(I+4) = -(Y(I+4)-Y(I))/p.RGS - fn_IBS(Y(I+2)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+6))/p.RGD - fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+5) = p.CGS*U1D - (Y(I+5)-Y(I+9))/p.RGS ...
             - fn_IDS(2, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7)-Y(I+9), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+6) = p.CGD*U1D - (Y(I+6)-Y(I+4))/p.RGD ...
             + fn_IDS(2, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7)-Y(I+9), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+7) = -(Y(I+7)-p.VBB)/p.RBS + fn_IBS(Y(I+7)-Y(I+9), p);
    F(I+8) = -(Y(I+8)-p.VBB)/p.RBD + fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+9) = -(Y(I+9)-Y(I+5))/p.RGS - fn_IBS(Y(I+7)-Y(I+9), p) ...
             -(Y(I+9)-Y(I+11))/p.RGD - fn_IBD(Y(I+13)-Y(I+9), p);

    F(I+10) = p.CGS*U2D - Y(I+10)/p.RGS ...
              - fn_IDS(2, Y(I+11)-Y(I+10), U2-Y(I+10), Y(I+12), U2-Y(I+11), Y(I+13)-Y(I+9), p);
    F(I+11) = p.CGD*U2D - (Y(I+11)-Y(I+9))/p.RGD ...
              + fn_IDS(2, Y(I+11)-Y(I+10), U2-Y(I+10), Y(I+12), U2-Y(I+11), Y(I+13)-Y(I+9), p);
    F(I+12) = -(Y(I+12)-p.VBB)/p.RBS + fn_IBS(Y(I+12), p);
    F(I+13) = -(Y(I+13)-p.VBB)/p.RBD + fn_IBD(Y(I+13)-Y(I+9), p);
end

function F = gate_ORANI(N, I, U1, U2, U3, U1D, U2D, U3D, Y, F, p)
    F(I) = -(Y(I)-Y(I+4))/p.RGS ...
           - fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+1) = -(Y(I+1)-p.VDD)/p.RGD ...
             + fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+2) = -(Y(I+2)-p.VBB)/p.RBS + fn_IBS(Y(I+2)-Y(I+4), p);
    F(I+3) = -(Y(I+3)-p.VBB)/p.RBD + fn_IBD(Y(I+3)-p.VDD, p);

    F(I+4) = -(Y(I+4)-Y(I))/p.RGS - fn_IBS(Y(I+2)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+6))/p.RGD - fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+5) = p.CGS*U1D - (Y(I+5)-Y(I+9))/p.RGS ...
             - fn_IDS(2, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7)-Y(I+9), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+6) = p.CGD*U1D - (Y(I+6)-Y(I+4))/p.RGD ...
             + fn_IDS(2, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7)-Y(I+9), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+7) = -(Y(I+7)-p.VBB)/p.RBS + fn_IBS(Y(I+7)-Y(I+9), p);
    F(I+8) = -(Y(I+8)-p.VBB)/p.RBD + fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+9) = -(Y(I+9)-Y(I+5))/p.RGS - fn_IBS(Y(I+7)-Y(I+9), p) ...
             -(Y(I+9)-Y(I+11))/p.RGD - fn_IBD(Y(I+13)-Y(I+9), p) ...
             -(Y(I+9)-Y(I+15))/p.RGD - fn_IBD(Y(I+17)-Y(I+9), p);

    F(I+10) = p.CGS*U2D - Y(I+10)/p.RGS ...
              - fn_IDS(2, Y(I+11)-Y(I+10), U2-Y(I+10), Y(I+12), U2-Y(I+11), Y(I+13)-Y(I+9), p);
    F(I+11) = p.CGD*U2D - (Y(I+11)-Y(I+9))/p.RGD ...
              + fn_IDS(2, Y(I+11)-Y(I+10), U2-Y(I+10), Y(I+12), U2-Y(I+11), Y(I+13)-Y(I+9), p);
    F(I+12) = -(Y(I+12)-p.VBB)/p.RBS + fn_IBS(Y(I+12), p);
    F(I+13) = -(Y(I+13)-p.VBB)/p.RBD + fn_IBD(Y(I+13)-Y(I+9), p);

    F(I+14) = p.CGS*U3D - Y(I+14)/p.RGS ...
              - fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+9), p);
    F(I+15) = p.CGD*U3D - (Y(I+15)-Y(I+9))/p.RGD ...
              + fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+9), p);
    F(I+16) = -(Y(I+16)-p.VBB)/p.RBS + fn_IBS(Y(I+16), p);
    F(I+17) = -(Y(I+17)-p.VBB)/p.RBD + fn_IBD(Y(I+17)-Y(I+9), p);
end

function F = gate_ANDOIP(N, I, U1, U2, U3, U1D, U2D, U3D, Y, F, p)
    % Same as ANDOI but result node I+4 has extra coupling to nodes 163,165 (1-based)
    F(I) = -(Y(I)-Y(I+4))/p.RGS ...
           - fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+1) = -(Y(I+1)-p.VDD)/p.RGD ...
             + fn_IDS(0, Y(I+1)-Y(I), Y(I+4)-Y(I), Y(I+2)-Y(I+4), Y(I+4)-Y(I+1), Y(I+3)-p.VDD, p);
    F(I+2) = -(Y(I+2)-p.VBB)/p.RBS + fn_IBS(Y(I+2)-Y(I+4), p);
    F(I+3) = -(Y(I+3)-p.VBB)/p.RBD + fn_IBD(Y(I+3)-p.VDD, p);

    % Result node I+4 with extra coupling to Y(163), Y(165)
    F(I+4) = -(Y(I+4)-Y(I))/p.RGS - fn_IBS(Y(I+2)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+6))/p.RGD - fn_IBD(Y(I+8)-Y(I+4), p) ...
             -(Y(I+4)-Y(I+10))/p.RGD - fn_IBD(Y(I+12)-Y(I+4), p) ...
             -(Y(I+4)-Y(163))/p.RGD - fn_IBD(Y(165)-Y(I+4), p);

    F(I+5) = p.CGS*U1D - Y(I+5)/p.RGS ...
             - fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+6) = p.CGD*U1D - (Y(I+6)-Y(I+4))/p.RGD ...
             + fn_IDS(1, Y(I+6)-Y(I+5), U1-Y(I+5), Y(I+7), U1-Y(I+6), Y(I+8)-Y(I+4), p);
    F(I+7) = -(Y(I+7)-p.VBB)/p.RBS + fn_IBS(Y(I+7), p);
    F(I+8) = -(Y(I+8)-p.VBB)/p.RBD + fn_IBD(Y(I+8)-Y(I+4), p);

    F(I+9) = p.CGS*U2D - (Y(I+9)-Y(I+13))/p.RGS ...
             - fn_IDS(2, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11)-Y(I+13), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+10) = p.CGD*U2D - (Y(I+10)-Y(I+4))/p.RGD ...
              + fn_IDS(2, Y(I+10)-Y(I+9), U2-Y(I+9), Y(I+11)-Y(I+13), U2-Y(I+10), Y(I+12)-Y(I+4), p);
    F(I+11) = -(Y(I+11)-p.VBB)/p.RBS + fn_IBS(Y(I+11)-Y(I+13), p);
    F(I+12) = -(Y(I+12)-p.VBB)/p.RBD + fn_IBD(Y(I+12)-Y(I+4), p);

    F(I+13) = -(Y(I+13)-Y(I+9))/p.RGS - fn_IBS(Y(I+11)-Y(I+13), p) ...
              -(Y(I+13)-Y(I+15))/p.RGD - fn_IBD(Y(I+17)-Y(I+13), p);

    F(I+14) = p.CGS*U3D - Y(I+14)/p.RGS ...
              - fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+13), p);
    F(I+15) = p.CGD*U3D - (Y(I+15)-Y(I+13))/p.RGD ...
              + fn_IDS(2, Y(I+15)-Y(I+14), U3-Y(I+14), Y(I+16), U3-Y(I+15), Y(I+17)-Y(I+13), p);
    F(I+16) = -(Y(I+16)-p.VBB)/p.RBS + fn_IBS(Y(I+16), p);
    F(I+17) = -(Y(I+17)-p.VBB)/p.RBD + fn_IBD(Y(I+17)-Y(I+13), p);
end

%% --- FCN: Circuit equations for 175 voltage nodes (1-based) --------------
function F = FCN(N, X, Y, p)
    F = zeros(N, 1);

    % Input signals
    [V1, V1D] = PULSE(X, 0, 5, 0, 5, 5, 5, 20);
    [V2, V2D] = PULSE(X, 0, 5, 10, 5, 15, 5, 40);
    [V3, V3D] = PULSE(X, 0, 5, 30, 5, 35, 5, 80);
    [V4, V4D] = PULSE(X, 0, 5, 70, 5, 75, 5, 160);
    [CIN, CIND] = PULSE(X, 0, 5, 150, 5, 155, 5, 320);

    % NOR-gate 1: nodes 1--13
    F = gate_NOR(N, 1, V1, V2, V1D, V2D, Y, F, p);

    % ANDOI-gate 1: nodes 14--31
    F = gate_ANDOI(N, 14, Y(5), V2, V1, 0, V2D, V1D, Y, F, p);

    % NOR-gate 2: nodes 32--44
    F = gate_NOR(N, 32, Y(18), CIN, 0, CIND, Y, F, p);

    % ANDOI-gate 2: nodes 45--62
    F = gate_ANDOI(N, 45, Y(36), CIN, Y(18), 0, CIND, 0, Y, F, p);

    % ANDOI-gate 3: nodes 63--80
    F = gate_ANDOI(N, 63, Y(5), CIN, Y(18), 0, CIND, 0, Y, F, p);

    % NOR-gate 3: nodes 81--93
    F = gate_NOR(N, 81, V3, V4, V3D, V4D, Y, F, p);

    % ANDOI-gate 4: nodes 94--111
    F = gate_ANDOI(N, 94, Y(85), V4, V3, 0, V4D, V3D, Y, F, p);

    % NAND-gate: nodes 112--125
    F = gate_NAND(N, 112, Y(67), Y(98), 0, 0, Y, F, p);

    % ORANI-gate 1: nodes 126--143
    F = gate_ORANI(N, 126, Y(116), Y(67), Y(98), 0, 0, 0, Y, F, p);

    % ANDOIP-gate 5: nodes 144--161
    F = gate_ANDOIP(N, 144, Y(85), Y(5), Y(98), 0, 0, 0, Y, F, p);

    % Three additional enhancement transistors: nodes 162--175
    F(162) = -(Y(162)-Y(166))/p.RGS ...
             - fn_IDS(3, Y(163)-Y(162), Y(98)-Y(162), Y(164)-Y(166), Y(98)-Y(163), Y(165)-Y(148), p);
    F(163) = -(Y(163)-Y(148))/p.RGD ...
             + fn_IDS(3, Y(163)-Y(162), Y(98)-Y(162), Y(164)-Y(166), Y(98)-Y(163), Y(165)-Y(148), p);
    F(164) = -(Y(164)-p.VBB)/p.RBS + fn_IBS(Y(164)-Y(166), p);
    F(165) = -(Y(165)-p.VBB)/p.RBD + fn_IBD(Y(165)-Y(148), p);

    F(166) = -fn_IBS(Y(164)-Y(166), p) - (Y(166)-Y(162))/p.RGS ...
             -fn_IBD(Y(170)-Y(166), p) - (Y(166)-Y(168))/p.RGD;

    F(167) = -(Y(167)-Y(171))/p.RGS ...
             - fn_IDS(3, Y(168)-Y(167), Y(18)-Y(167), Y(169)-Y(171), Y(18)-Y(168), Y(170)-Y(166), p);
    F(168) = -(Y(168)-Y(166))/p.RGD ...
             + fn_IDS(3, Y(168)-Y(167), Y(18)-Y(167), Y(169)-Y(171), Y(18)-Y(168), Y(170)-Y(166), p);
    F(169) = -(Y(169)-p.VBB)/p.RBS + fn_IBS(Y(169)-Y(171), p);
    F(170) = -(Y(170)-p.VBB)/p.RBD + fn_IBD(Y(170)-Y(166), p);

    F(171) = -fn_IBS(Y(169)-Y(171), p) - (Y(171)-Y(167))/p.RGS ...
             -fn_IBD(Y(175)-Y(171), p) - (Y(171)-Y(173))/p.RGD;

    F(172) = p.CGS*CIND - Y(172)/p.RGS ...
             - fn_IDS(3, Y(173)-Y(172), CIN-Y(172), Y(174), CIN-Y(173), Y(175)-Y(171), p);
    F(173) = p.CGD*CIND - (Y(173)-Y(171))/p.RGD ...
             + fn_IDS(3, Y(173)-Y(172), CIN-Y(172), Y(174), CIN-Y(173), Y(175)-Y(171), p);
    F(174) = -(Y(174)-p.VBB)/p.RBS + fn_IBS(Y(174), p);
    F(175) = -(Y(175)-p.VBB)/p.RBD + fn_IBD(Y(175)-Y(171), p);
end

%% --- GCN: Charge function for 175 voltage nodes (1-based) ---------------
function G = GCN(U, p)
    N = 175;
    G = zeros(N, 1);

    % Ten logical subcircuits
    G = charge_NOR(N, U, 1, G, p);
    G = charge_ANDOI(N, U, 14, G, p);
    G = charge_NOR(N, U, 32, G, p);
    G = charge_ANDOI(N, U, 45, G, p);
    G = charge_ANDOI(N, U, 63, G, p);
    G = charge_NOR(N, U, 81, G, p);
    G = charge_ANDOI(N, U, 94, G, p);
    G = charge_NAND(N, U, 112, G, p);
    G = charge_ORANI(N, U, 126, G, p);
    G = charge_ANDOI(N, U, 144, G, p);

    % Capacitive coupling: result node NOR-gate 1 (node 5)
    G(5) = G(5) + p.CGS*(U(5)-U(19)) + p.CGD*(U(5)-U(20)) ...
                + p.CGS*(U(5)-U(68)) + p.CGD*(U(5)-U(69)) ...
                + p.CGS*(U(5)-U(153)) + p.CGD*(U(5)-U(154));
    G(19)  = G(19)  - p.CGS*U(5);
    G(20)  = G(20)  - p.CGD*U(5);
    G(68)  = G(68)  - p.CGS*U(5);
    G(69)  = G(69)  - p.CGD*U(5);
    G(153) = G(153) - p.CGS*U(5);
    G(154) = G(154) - p.CGD*U(5);

    % Capacitive coupling: result node ANDOI-gate 1 (node 18)
    G(18) = G(18) + p.CGS*(U(18)-U(37)) + p.CGD*(U(18)-U(38)) ...
                  + p.CGS*(U(18)-U(59)) + p.CGD*(U(18)-U(60)) ...
                  + p.CGS*(U(18)-U(77)) + p.CGD*(U(18)-U(78)) ...
                  + p.CGS*(U(18)-U(167)) + p.CGD*(U(18)-U(168));
    G(37) = G(37) - p.CGS*U(18);
    G(38) = G(38) - p.CGD*U(18);
    G(59) = G(59) - p.CGS*U(18);
    G(60) = G(60) - p.CGD*U(18);
    G(77) = G(77) - p.CGS*U(18);
    G(78) = G(78) - p.CGD*U(18);

    % Capacitive coupling: result node NOR-gate 2 (node 36)
    G(36) = G(36) + p.CGS*(U(36)-U(50)) + p.CGD*(U(36)-U(51));
    G(50) = G(50) - p.CGS*U(36);
    G(51) = G(51) - p.CGD*U(36);

    % Capacitive coupling: result node ANDOI-gate 2 = S0 (node 49)
    G(49) = G(49) + p.COUT*U(49);

    % Capacitive coupling: result node ANDOI-gate 3 (node 67)
    G(67) = G(67) + p.CGS*(U(67)-U(117)) + p.CGD*(U(67)-U(118)) ...
                  + p.CGS*(U(67)-U(136)) + p.CGD*(U(67)-U(137));
    G(117) = G(117) - p.CGS*U(67);
    G(118) = G(118) - p.CGD*U(67);
    G(136) = G(136) - p.CGS*U(67);
    G(137) = G(137) - p.CGD*U(67);

    % Capacitive coupling: result node NOR-gate 3 (node 85)
    G(85) = G(85) + p.CGS*(U(85)-U(99)) + p.CGD*(U(85)-U(100)) ...
                  + p.CGS*(U(85)-U(149)) + p.CGD*(U(85)-U(150));
    G(99)  = G(99)  - p.CGS*U(85);
    G(100) = G(100) - p.CGD*U(85);
    G(149) = G(149) - p.CGS*U(85);
    G(150) = G(150) - p.CGD*U(85);

    % Capacitive coupling: result node ANDOI-gate 4 (node 98)
    G(98) = G(98) + p.CGS*(U(98)-U(122)) + p.CGD*(U(98)-U(123)) ...
                  + p.CGS*(U(98)-U(140)) + p.CGD*(U(98)-U(141)) ...
                  + p.CGS*(U(98)-U(158)) + p.CGD*(U(98)-U(159)) ...
                  + p.CGS*(U(98)-U(162)) + p.CGD*(U(98)-U(163));
    G(122) = G(122) - p.CGS*U(98);
    G(123) = G(123) - p.CGD*U(98);
    G(140) = G(140) - p.CGS*U(98);
    G(141) = G(141) - p.CGD*U(98);
    G(158) = G(158) - p.CGS*U(98);
    G(159) = G(159) - p.CGD*U(98);

    % Capacitive coupling: result NAND-gate (node 116)
    G(116) = G(116) + p.CGS*(U(116)-U(131)) + p.CGD*(U(116)-U(132));
    G(131) = G(131) - p.CGS*U(116);
    G(132) = G(132) - p.CGD*U(116);

    % Capacitive coupling: result node ORANI-gate = S1 (node 130)
    G(130) = G(130) + p.COUT*U(130);

    % Capacitive coupling: result ANDOI-gate 5 = Cinvers (node 148)
    G(148) = G(148) + CBDBS(U(165)-U(148), p)*(U(148)-U(165)) + p.COUT*U(148);

    % Charge function of three additional transistors: nodes 162--175
    G(162) = G(162) + p.CGS*(U(162)-U(98));
    G(163) = G(163) + p.CGD*(U(163)-U(98));
    G(164) = G(164) + CBDBS(U(164)-U(166), p)*(U(164)-U(166));
    G(165) = G(165) + CBDBS(U(165)-U(148), p)*(U(165)-U(148));
    G(166) = G(166) + CBDBS(U(164)-U(166), p)*(U(166)-U(164)) ...
                    + CBDBS(U(170)-U(166), p)*(U(166)-U(170)) ...
                    + p.CLOAD*U(166);
    G(167) = G(167) + p.CGS*(U(167)-U(18));
    G(168) = G(168) + p.CGD*(U(168)-U(18));
    G(169) = G(169) + CBDBS(U(169)-U(171), p)*(U(169)-U(171));
    G(170) = G(170) + CBDBS(U(170)-U(166), p)*(U(170)-U(166));
    G(171) = G(171) + CBDBS(U(169)-U(171), p)*(U(171)-U(169)) ...
                    + CBDBS(U(175)-U(171), p)*(U(171)-U(175)) ...
                    + p.CLOAD*U(171);
    G(172) = G(172) + p.CGS*U(172);
    G(173) = G(173) + p.CGD*U(173);
    G(174) = G(174) + CBDBS(U(174), p)*U(174);
    G(175) = G(175) + CBDBS(U(175)-U(171), p)*(U(175)-U(171));
end

%% --- Charge gate subroutines (1-based) -----------------------------------
function G = charge_NOR(N, U, I, G, p)
    G(I)   = G(I)   + p.CGS*(U(I)-U(I+4));
    G(I+1) = G(I+1) + p.CGD*(U(I+1)-U(I+4));
    G(I+2) = G(I+2) + CBDBS(U(I+2)-U(I+4), p)*(U(I+2)-U(I+4));
    G(I+3) = G(I+3) + CBDBS(U(I+3)-p.VDD, p)*U(I+3);
    G(I+4) = G(I+4) + p.CGS*(U(I+4)-U(I)) + p.CGD*(U(I+4)-U(I+1)) ...
           + CBDBS(U(I+2)-U(I+4), p)*(U(I+4)-U(I+2)) ...
           + CBDBS(U(I+8)-U(I+4), p)*(U(I+4)-U(I+8)) ...
           + CBDBS(U(I+12)-U(I+4), p)*(U(I+4)-U(I+12)) ...
           + p.CLOAD*U(I+4);
    G(I+5)  = G(I+5)  + p.CGS*U(I+5);
    G(I+6)  = G(I+6)  + p.CGD*U(I+6);
    G(I+7)  = G(I+7)  + CBDBS(U(I+7), p)*U(I+7);
    G(I+8)  = G(I+8)  + CBDBS(U(I+8)-U(I+4), p)*(U(I+8)-U(I+4));
    G(I+9)  = G(I+9)  + p.CGS*U(I+9);
    G(I+10) = G(I+10) + p.CGD*U(I+10);
    G(I+11) = G(I+11) + CBDBS(U(I+11), p)*U(I+11);
    G(I+12) = G(I+12) + CBDBS(U(I+12)-U(I+4), p)*(U(I+12)-U(I+14));
end

function G = charge_ANDOI(N, U, I, G, p)
    G(I)   = G(I)   + p.CGS*(U(I)-U(I+4));
    G(I+1) = G(I+1) + p.CGD*(U(I+1)-U(I+4));
    G(I+2) = G(I+2) + CBDBS(U(I+2)-U(I+4), p)*(U(I+2)-U(I+4));
    G(I+3) = G(I+3) + CBDBS(U(I+3)-p.VDD, p)*U(I+3);
    G(I+4) = G(I+4) + p.CGS*(U(I+4)-U(I)) + p.CGD*(U(I+4)-U(I+1)) ...
           + CBDBS(U(I+2)-U(I+4), p)*(U(I+4)-U(I+2)) ...
           + CBDBS(U(I+8)-U(I+4), p)*(U(I+4)-U(I+8)) ...
           + CBDBS(U(I+12)-U(I+4), p)*(U(I+4)-U(I+12)) ...
           + p.CLOAD*U(I+4);
    G(I+5)  = G(I+5)  + p.CGS*U(I+5);
    G(I+6)  = G(I+6)  + p.CGD*U(I+6);
    G(I+7)  = G(I+7)  + CBDBS(U(I+7), p)*U(I+7);
    G(I+8)  = G(I+8)  + CBDBS(U(I+8)-U(I+4), p)*(U(I+8)-U(I+4));
    G(I+9)  = G(I+9)  + p.CGS*U(I+9);
    G(I+10) = G(I+10) + p.CGD*U(I+10);
    G(I+11) = G(I+11) + CBDBS(U(I+11)-U(I+13), p)*(U(I+11)-U(I+13));
    G(I+12) = G(I+12) + CBDBS(U(I+12)-U(I+4), p)*(U(I+12)-U(I+4));
    G(I+13) = G(I+13) + CBDBS(U(I+11)-U(I+13), p)*(U(I+13)-U(I+11)) ...
            + CBDBS(U(I+17)-U(I+13), p)*(U(I+13)-U(I+17)) ...
            + p.CLOAD*U(I+13);
    G(I+14) = G(I+14) + p.CGS*U(I+14);
    G(I+15) = G(I+15) + p.CGD*U(I+15);
    G(I+16) = G(I+16) + CBDBS(U(I+16), p)*U(I+16);
    G(I+17) = G(I+17) + CBDBS(U(I+17)-U(I+13), p)*(U(I+17)-U(I+13));
end

function G = charge_NAND(N, U, I, G, p)
    G(I)   = G(I)   + p.CGS*(U(I)-U(I+4));
    G(I+1) = G(I+1) + p.CGD*(U(I+1)-U(I+4));
    G(I+2) = G(I+2) + CBDBS(U(I+2)-U(I+4), p)*(U(I+2)-U(I+4));
    G(I+3) = G(I+3) + CBDBS(U(I+3)-p.VDD, p)*U(I+3);
    G(I+4) = G(I+4) + p.CGS*(U(I+4)-U(I)) + p.CGD*(U(I+4)-U(I+1)) ...
           + CBDBS(U(I+2)-U(I+4), p)*(U(I+4)-U(I+2)) ...
           + CBDBS(U(I+8)-U(I+4), p)*(U(I+4)-U(I+8)) ...
           + p.CLOAD*U(I+4);
    G(I+5)  = G(I+5)  + p.CGS*U(I+5);
    G(I+6)  = G(I+6)  + p.CGD*U(I+6);
    G(I+7)  = G(I+7)  + CBDBS(U(I+7)-U(I+9), p)*(U(I+7)-U(I+9));
    G(I+8)  = G(I+8)  + CBDBS(U(I+8)-U(I+4), p)*(U(I+8)-U(I+4));
    G(I+9)  = G(I+9)  + CBDBS(U(I+7)-U(I+9), p)*(U(I+9)-U(I+7)) ...
            + CBDBS(U(I+13)-U(I+9), p)*(U(I+9)-U(I+13)) ...
            + p.CLOAD*U(I+9);
    G(I+10) = G(I+10) + p.CGS*U(I+10);
    G(I+11) = G(I+11) + p.CGD*U(I+11);
    G(I+12) = G(I+12) + CBDBS(U(I+12), p)*U(I+12);
    G(I+13) = G(I+13) + CBDBS(U(I+13)-U(I+9), p)*(U(I+13)-U(I+9));
end

function G = charge_ORANI(N, U, I, G, p)
    G(I)   = G(I)   + p.CGS*(U(I)-U(I+4));
    G(I+1) = G(I+1) + p.CGD*(U(I+1)-U(I+4));
    G(I+2) = G(I+2) + CBDBS(U(I+2)-U(I+4), p)*(U(I+2)-U(I+4));
    G(I+3) = G(I+3) + CBDBS(U(I+3)-p.VDD, p)*U(I+3);
    G(I+4) = G(I+4) + p.CGS*(U(I+4)-U(I)) + p.CGD*(U(I+4)-U(I+1)) ...
           + CBDBS(U(I+2)-U(I+4), p)*(U(I+4)-U(I+2)) ...
           + CBDBS(U(I+8)-U(I+4), p)*(U(I+4)-U(I+8)) ...
           + p.CLOAD*U(I+4);
    G(I+5)  = G(I+5)  + p.CGS*U(I+5);
    G(I+6)  = G(I+6)  + p.CGD*U(I+6);
    G(I+7)  = G(I+7)  + CBDBS(U(I+7)-U(I+9), p)*(U(I+7)-U(I+9));
    G(I+8)  = G(I+8)  + CBDBS(U(I+8)-U(I+4), p)*(U(I+8)-U(I+4));
    G(I+9)  = G(I+9)  + CBDBS(U(I+7)-U(I+9), p)*(U(I+9)-U(I+7)) ...
            + CBDBS(U(I+13)-U(I+9), p)*(U(I+9)-U(I+13)) ...
            + CBDBS(U(I+17)-U(I+9), p)*(U(I+9)-U(I+17)) ...
            + p.CLOAD*U(I+9);
    G(I+10) = G(I+10) + p.CGS*U(I+10);
    G(I+11) = G(I+11) + p.CGD*U(I+11);
    G(I+12) = G(I+12) + CBDBS(U(I+12), p)*U(I+12);
    G(I+13) = G(I+13) + CBDBS(U(I+13)-U(I+9), p)*(U(I+13)-U(I+9));
    G(I+14) = G(I+14) + p.CGS*U(I+14);
    G(I+15) = G(I+15) + p.CGD*U(I+15);
    G(I+16) = G(I+16) + CBDBS(U(I+16), p)*U(I+16);
    G(I+17) = G(I+17) + CBDBS(U(I+17)-U(I+9), p)*(U(I+17)-U(I+9));
end

%% --- Parameter initialization --------------------------------------------
function p = init_params()
    p.CTIME = 1e4;
    p.STIFF = 5.0;
    p.RGS   = 40 / (p.CTIME * p.STIFF);
    p.RGD   = 40 / (p.CTIME * p.STIFF);
    p.RBS   = 100 / (p.CTIME * p.STIFF);
    p.RBD   = 100 / (p.CTIME * p.STIFF);
    p.CGS   = 0.6e-4 * p.CTIME;
    p.CGD   = 0.6e-4 * p.CTIME;
    p.CBD   = 2.4e-5 * p.CTIME;
    p.CBS   = 2.4e-5 * p.CTIME;
    p.DELTA = 0.02;
    p.CURIS = 1e-15 * p.CTIME * p.STIFF;
    p.VTH   = 25.85;
    p.VDD   = 5.0;
    p.VBB   = -2.5;
    p.CLOAD = 0.0;
    p.COUT  = 2e-4 * p.CTIME - p.CLOAD;
end

%% --- Initial voltages (175 values) --------------------------------------
function U = init_voltages()
    U = zeros(175, 1);
    U(  1) =  4.999999999996544e+00;
    U(  2) =  4.999999999999970e+00;
    U(  3) = -2.499999999999975e+00;
    U(  4) = -2.499999999999975e+00;
    U(  5) =  4.999999999996514e+00;
    U(  6) =  0.000000000000000e+00;
    U(  7) =  4.999999999996514e+00;
    U(  8) = -2.499999999999991e+00;
    U(  9) = -2.499999999999975e+00;
    U( 10) =  0.000000000000000e+00;
    U( 11) =  4.999999999996514e+00;
    U( 12) = -2.499999999999991e+00;
    U( 13) = -2.499999999999975e+00;
    U( 14) =  0.215858486765796e+00;
    U( 15) =  4.988182208251953e+00;
    U( 16) = -2.499999999999990e+00;
    U( 17) = -2.499999999999975e+00;
    U( 18) =  0.204040695017748e+00;
    U( 19) =  0.011817791748026e+00;
    U( 20) =  0.192222903269723e+00;
    U( 21) = -2.499999999999991e+00;
    U( 22) = -2.499999999999990e+00;
    U( 23) = -0.228160951881239e+00;
    U( 24) =  0.204040695017748e+00;
    U( 25) = -2.499999999999992e+00;
    U( 26) = -2.499999999999990e+00;
    U( 27) = -0.228160951881241e+00;
    U( 28) =  0.000000000000000e+00;
    U( 29) = -0.228160951881239e+00;
    U( 30) = -2.499999999999991e+00;
    U( 31) = -2.499999999999992e+00;
    U( 32) =  4.999999999996547e+00;
    U( 33) =  4.999999999999970e+00;
    U( 34) = -2.499999999999975e+00;
    U( 35) = -2.499999999999975e+00;
    U( 36) =  4.999999999996517e+00;
    U( 37) =  0.000000000000000e+00;
    U( 38) =  4.999999999996517e+00;
    U( 39) = -2.499999999999991e+00;
    U( 40) = -2.499999999999975e+00;
    U( 41) =  0.000000000000000e+00;
    U( 42) =  4.999999999996517e+00;
    U( 43) = -2.499999999999991e+00;
    U( 44) = -2.499999999999975e+00;
    U( 45) =  0.215858484247529e+00;
    U( 46) =  4.988182208251953e+00;
    U( 47) = -2.499999999999990e+00;
    U( 48) = -2.499999999999975e+00;
    U( 49) =  0.204040692499482e+00;
    U( 50) =  0.011817791748035e+00;
    U( 51) =  0.192222900751447e+00;
    U( 52) = -2.499999999999991e+00;
    U( 53) = -2.499999999999990e+00;
    U( 54) = -0.026041071738432e+00;
    U( 55) =  0.204040692499482e+00;
    U( 56) = -2.499999999999992e+00;
    U( 57) = -2.499999999999990e+00;
    U( 58) = -0.026041071738434e+00;
    U( 59) =  0.000000000000000e+00;
    U( 60) = -0.026041071738432e+00;
    U( 61) = -2.499999999999991e+00;
    U( 62) = -2.499999999999992e+00;
    U( 63) =  0.215858484880918e+00;
    U( 64) =  4.988182208251953e+00;
    U( 65) = -2.499999999999990e+00;
    U( 66) = -2.499999999999975e+00;
    U( 67) =  0.204040693132870e+00;
    U( 68) =  0.011817791748026e+00;
    U( 69) =  0.192222901384845e+00;
    U( 70) = -2.499999999999991e+00;
    U( 71) = -2.499999999999990e+00;
    U( 72) = -0.026041071737961e+00;
    U( 73) =  0.204040693132870e+00;
    U( 74) = -2.499999999999992e+00;
    U( 75) = -2.499999999999990e+00;
    U( 76) = -0.026041071737963e+00;
    U( 77) =  0.000000000000000e+00;
    U( 78) = -0.026041071737961e+00;
    U( 79) = -2.499999999999991e+00;
    U( 80) = -2.499999999999992e+00;
    U( 81) =  4.999999999996546e+00;
    U( 82) =  4.999999999999970e+00;
    U( 83) = -2.499999999999975e+00;
    U( 84) = -2.499999999999975e+00;
    U( 85) =  4.999999999996516e+00;
    U( 86) =  0.000000000000000e+00;
    U( 87) =  4.999999999996516e+00;
    U( 88) = -2.499999999999991e+00;
    U( 89) = -2.499999999999975e+00;
    U( 90) =  0.000000000000000e+00;
    U( 91) =  4.999999999996516e+00;
    U( 92) = -2.499999999999991e+00;
    U( 93) = -2.499999999999975e+00;
    U( 94) =  0.215858481060569e+00;
    U( 95) =  4.988182208251953e+00;
    U( 96) = -2.499999999999990e+00;
    U( 97) = -2.499999999999975e+00;
    U( 98) =  0.204040689312522e+00;
    U( 99) =  0.011817791748023e+00;
    U(100) =  0.192222897564498e+00;
    U(101) = -2.499999999999991e+00;
    U(102) = -2.499999999999990e+00;
    U(103) =  4.734672533390068e+00;
    U(104) =  0.204040689312522e+00;
    U(105) = -2.499999999999977e+00;
    U(106) = -2.499999999999990e+00;
    U(107) =  4.734672533390062e+00;
    U(108) =  0.000000000000000e+00;
    U(109) =  4.734672533390068e+00;
    U(110) = -2.499999999999991e+00;
    U(111) = -2.499999999999977e+00;
    U(112) =  4.999999999996870e+00;
    U(113) =  4.999999999999972e+00;
    U(114) = -2.499999999999975e+00;
    U(115) = -2.499999999999975e+00;
    U(116) =  4.999999999996843e+00;
    U(117) = -0.025968303070038e+00;
    U(118) =  4.999999999996843e+00;
    U(119) = -2.499999999999992e+00;
    U(120) = -2.499999999999975e+00;
    U(121) = -0.025968303070040e+00;
    U(122) =  0.000000000000000e+00;
    U(123) = -0.025968303070038e+00;
    U(124) = -2.499999999999991e+00;
    U(125) = -2.499999999999992e+00;
    U(126) =  4.999999999997699e+00;
    U(127) =  4.999999999999980e+00;
    U(128) = -2.499999999999975e+00;
    U(129) = -2.499999999999975e+00;
    U(130) =  4.999999999997678e+00;
    U(131) =  4.744923533081106e+00;
    U(132) =  4.999999999997678e+00;
    U(133) = -2.499999999999977e+00;
    U(134) = -2.499999999999975e+00;
    U(135) =  4.744923533081098e+00;
    U(136) =  0.000000000000000e+00;
    U(137) =  4.744923533081106e+00;
    U(138) = -2.499999999999991e+00;
    U(139) = -2.499999999999977e+00;
    U(140) =  0.000000000000000e+00;
    U(141) =  4.744923533081106e+00;
    U(142) = -2.499999999999991e+00;
    U(143) = -2.499999999999977e+00;
    U(144) =  0.215858484844162e+00;
    U(145) =  4.988182208251953e+00;
    U(146) = -2.499999999999990e+00;
    U(147) = -2.499999999999975e+00;
    U(148) =  0.204040693096114e+00;
    U(149) =  0.011817791748023e+00;
    U(150) =  0.192222901348091e+00;
    U(151) = -2.499999999999991e+00;
    U(152) = -2.499999999999990e+00;
    U(153) =  0.204040693096045e+00;
    U(154) =  0.204040693096107e+00;
    U(155) = -2.499999999999990e+00;
    U(156) = -2.499999999999990e+00;
    U(157) =  0.204040693096037e+00;
    U(158) =  0.000000000000000e+00;
    U(159) =  0.204040693096037e+00;
    U(160) = -2.499999999999991e+00;
    U(161) = -2.499999999999990e+00;
    U(162) = -0.026017361873565e+00;
    U(163) =  0.204040693096114e+00;
    U(164) = -2.499999999999992e+00;
    U(165) = -2.499999999999990e+00;
    U(166) = -0.026017361873568e+00;
    U(167) = -0.026017590106916e+00;
    U(168) = -0.026017361873565e+00;
    U(169) = -2.499999999999992e+00;
    U(170) = -2.499999999999992e+00;
    U(171) = -0.026017590106918e+00;
    U(172) =  0.000000000000000e+00;
    U(173) = -0.026017590106916e+00;
    U(174) = -2.499999999999991e+00;
    U(175) = -2.499999999999992e+00;
end
