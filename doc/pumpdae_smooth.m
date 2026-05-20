function pumpdae_smooth
%PUMPDAE_SMOOTH Charge pump DAE with smoothed charge functions.
%  The original charge functions QG, QS, QD have C0-but-not-C1 transitions
%  at the saturation/linear boundary (max(ugd-vte,0)) and the depletion/
%  inversion boundary (ugs=vte). This version replaces those hard switches
%  with softplus/sigmoid smoothing to produce C-infinity right-hand sides.
%
%   See also ODEBIM, PUMPDAE.
%
%   Based on pumpdae.m by Wu Zhiqiao (2007)
%   Smoothed version: 2026-05-20

%% Parameters
clc; close all;
VT0   = 0.20;
GAMMA = 0.035;
PHI   = 1.01;
COX   = 4.0e-12;
CAPD  = 0.40e-12;
CAPS  = 1.60e-12;

VHIGH  = 20.0;
DELTAT = 120.0;
T1_PLS = 50.0;
T2_PLS = 60.0;
T3_PLS = 110.0;

% Smoothing parameter (voltage scale for transition width)
EPS_SMOOTH = 1e-4;

%% Mass matrix (constant, singular)
M = zeros(9,9);
M(1,1) = 1;
M(2,2) = 1;
M(2,3) = 1;
M(3,4) = 1;
M(3,5) = 1;

%% Initial conditions
y0 = zeros(9,1);
y0(1) = fn_qgate_smooth(0, 0, 0, VT0, GAMMA, PHI, COX, EPS_SMOOTH);
y0(3) = fn_qsrc_smooth(0, 0, 0, VT0, GAMMA, PHI, COX, EPS_SMOOTH);
y0(5) = fn_qdrain_smooth(0, 0, 0, VT0, GAMMA, PHI, COX, EPS_SMOOTH);

%% Tolerances
Tol = 1e-9;
rtol = Tol;
atol_v = [Tol*1e-6*ones(5,1); Tol*ones(3,1); Tol];

%% Discontinuity times (for Vin piecewise-linear input)
disc = zeros(41,1);
disc(1) = 0;
disc(2) = 50;
disc(3) = 60;
disc(4) = 110;
for i = 5:40
    disc(i) = disc(mod(i-1,4)+1) + floor((i-1)/4)*120;
end
disc(41) = 1200;

%% ode options
opts = odeset('Mass', M, ...
              'MassSingular', 'yes', ...
              'RelTol', rtol, ...
              'AbsTol', atol_v, ...
              'InitialStep', 1, ...
              'MaxStep', 1, ...
              'Jacobian', @(t,y) jac_pump_smooth(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, EPS_SMOOTH), ...
              'DaeIndex', [8,1,0]);

%% Segmented integration (restart at each discontinuity)
y = y0;
t_cur = 0;
total_steps = 0;
tout = [];
yout = [];

fprintf('Solving charge pump (smoothed, eps=%.1e) with odebim...\n', EPS_SMOOTH);
for seg = 1:40
    tspan = [disc(seg), disc(seg+1)];
    if tspan(2) <= tspan(1), continue; end

    odefun = @(t, y) rhs_pump_smooth(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, ...
                               VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS, EPS_SMOOTH);

    [tsol, ysol] = odebim(odefun, tspan, y, opts);

    y = ysol(end,:)';
    total_steps = total_steps + length(tsol) - 1;
    t_cur = tsol(end);
    tout = [tout; tsol];
    yout = [yout; ysol];
end

%% Plot
figure(1); plot(tout, yout(:,6)); title('U1 (gate voltage)');
figure(2); plot(tout, yout(:,7)); title('U2 (source voltage)');
figure(3); plot(tout, yout(:,8)); title('U3 (drain voltage)');
figure(4); plot(tout, yout(:,9)); title('I (current)');

%% Reference solution at t = 1200
yref = zeros(9,1);
yref(1) = 0.126280042987675933170893e-12;
yref(9) = 0.152255686815577679043511e-03;

%% Results
fprintf('\n=== Charge Pump Smoothed (eps=%.1e, Tol=%.1e) ===\n', EPS_SMOOTH, Tol);
fprintf('Final time: t = %.6e\n', t_cur);
fprintf('Total steps: %d\n\n', total_steps);

maxrel = 0;
for i = 1:9
    err = abs(y(i) - yref(i));
    if abs(yref(i)) > 1e-15
        relerr = err / abs(yref(i));
    else
        relerr = err;
    end
    if relerr > maxrel, maxrel = relerr; end
    fprintf('y(%d) = %22.14e  ref = %22.14e  rel_err = %.3e\n', ...
            i, y(i), yref(i), relerr);
end

rel1 = abs(y(1) - yref(1)) / abs(yref(1));
rel9 = abs(y(9) - yref(9)) / abs(yref(9));
maxrel_key = max(rel1, rel9);
fprintf('\nmax_rel_err (y(1),y(9)) = %.3e\n', maxrel_key);
if maxrel_key < 1e-2
    fprintf('PASSED\n');
else
    fprintf('FAILED\n');
end

%% =========================================================================
%  Local functions
%  =========================================================================

%% --- Smoothing helpers ---------------------------------------------------
function y = softplus(x, eps)
    % softplus(x,eps) = eps * log(1 + exp(x/eps))
    % Numerically stable implementation to avoid overflow
    z = x / eps;
    if z > 30
        y = x;           % exp(z) dominates, log(1+exp(z)) ≈ z
    elseif z < -30
        y = 0;           % exp(z) ≈ 0, log(1+exp(z)) ≈ 0
    else
        y = eps * log(1 + exp(z));
    end
end

function dy = dsoftplus(x, eps)
    % Derivative of softplus: sigmoid(x/eps)
    z = x / eps;
    if z > 30
        dy = 1;
    elseif z < -30
        dy = 0;
    else
        dy = 1 / (1 + exp(-z));
    end
end

function s = sigmoid(x, eps)
    % Smooth step function: 1/(1+exp(-x/eps))
    z = x / eps;
    if z > 30
        s = 1;
    elseif z < -30
        s = 0;
    else
        s = 1 / (1 + exp(-z));
    end
end

function ds = dsigmoid(x, eps)
    % Derivative of sigmoid: sigmoid*(1-sigmoid)/eps
    s = sigmoid(x, eps);
    ds = s * (1 - s) / eps;
end

%% --- RHS function --------------------------------------------------------
function f = rhs_pump_smooth(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, ...
                             VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS, EPS)
    v = y;
    vin_t = vin_func(t, VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS);

    vgb = v(6);
    vgs = v(6) - v(7);
    vgd = v(6) - v(8);

    f = zeros(9,1);
    f(1) = -1e-9 * v(9);
    f(2) = 0;
    f(3) = 0;
    f(4) = 1e-9 * (-v(6) + vin_t);
    f(5) = 1e-9 * (v(1) - fn_qgate_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS));
    f(6) = 1e-9 * (v(2) - CAPS * v(7));
    f(7) = 1e-9 * (v(3) - fn_qsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS));
    f(8) = 1e-9 * (v(4) - CAPD * v(8));
    f(9) = 1e-9 * (v(5) - fn_qdrain_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS));
end

%% --- Vin (unchanged, piecewise linear) ----------------------------------
function v = vin_func(t, VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS)
    dummy = mod(t, DELTAT);
    if dummy < T1_PLS
        v = 0;
    elseif dummy < T2_PLS
        v = (dummy - T1_PLS) * 0.10 * VHIGH;
    elseif dummy < T3_PLS
        v = VHIGH;
    else
        v = (DELTAT - dummy) * 0.10 * VHIGH;
    end
end

%% --- Smoothed charge functions -------------------------------------------
function q = fn_qgate_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    % Smoothed version of fn_qgate
    % Eliminates hard switches at ugs=vte and max(ugd-vte,0)
    if (vgs - vgd) <= 0
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    % Region A: q_A = COX*(ugb - vfb)
    q_A = COX * (ugb - vfb);

    % Region B: q_B = COX*GAMMA*(sqrt((GAMMA/2)^2 + ugb - vfb) - GAMMA/2)
    arg_B = (GAMMA/2)^2 + ugb - vfb;
    % Protect sqrt argument (should be positive when ugb > vfb)
    arg_B = max(arg_B, 1e-30);
    q_B = COX * GAMMA * (sqrt(arg_B) - GAMMA/2);

    % Region C: uses smoothed ugdt
    ugst = softplus(ugs - vte, EPS);  % smooth max(ugs-vte, 0) for robustness
    ugdt = softplus(ugd - vte, EPS);  % smooth max(ugd-vte, 0)
    S = ugst + ugdt;
    if S > 1e-30
        H = ugdt + ugst - (ugdt*ugst)/(ugdt + ugst);
    else
        H = 0;
    end
    q_C = COX * ((2/3)*H + GAMMA*sqrt(PHI - ubs));

    % Smooth blending: Region A/B boundary at ugb = vfb (already C1, keep hard)
    % Region B/C boundary at ugs = vte: use sigmoid blending
    sigma_BC = sigmoid(ugs - vte, EPS);

    if ugb <= vfb
        q = q_A;
    else
        q = (1 - sigma_BC) * q_B + sigma_BC * q_C;
    end
end

function q = fn_qsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    % Smoothed version of fn_qsrc
    if (vgs - vgd) <= 0
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    % Region A & B: q = 0
    % Region C: q = -COX*(1/3)*H
    ugst = softplus(ugs - vte, EPS);  % smooth max(ugs-vte, 0) for robustness
    ugdt = softplus(ugd - vte, EPS);  % smooth max(ugd-vte, 0)
    S = ugst + ugdt;
    if S > 1e-30
        H = ugdt + ugst - (ugdt*ugst)/(ugdt + ugst);
    else
        H = 0;
    end
    q_C = -COX * (1/3) * H;

    % Smooth blending at ugs = vte boundary
    sigma_BC = sigmoid(ugs - vte, EPS);

    % In regions A and B, q=0; blend to q_C in region C
    q = sigma_BC * q_C;
end

function q = fn_qdrain_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    % qdrain is identical to qsrc in this model
    q = fn_qsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS);
end

%% --- Smoothed Jacobian ---------------------------------------------------
function J = jac_pump_smooth(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, EPS)
    J = zeros(9,9);

    vgb = y(6);
    vgs = y(6) - y(7);
    vgd = y(6) - y(8);

    % Row 1: f(1) = -1e-9 * y(9)
    J(1,9) = -1e-9;

    % Row 4: f(4) = 1e-9 * (-y(6) + vin(t))
    J(4,6) = -1e-9;

    % Row 5: f(5) = 1e-9 * (y(1) - qgate(vgb, vgs, vgd))
    J(5,1) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqgate_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS);
    J(5,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(5,7) = -1e-9 * (-dq_dvgs);
    J(5,8) = -1e-9 * (-dq_dvgd);

    % Row 6: f(6) = 1e-9 * (y(2) - CAPS * y(7))
    J(6,2) = 1e-9;
    J(6,7) = -1e-9 * CAPS;

    % Row 7: f(7) = 1e-9 * (y(3) - qsrc(vgb, vgs, vgd))
    J(7,3) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS);
    J(7,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(7,7) = -1e-9 * (-dq_dvgs);
    J(7,8) = -1e-9 * (-dq_dvgd);

    % Row 8: f(8) = 1e-9 * (y(4) - CAPD * y(8))
    J(8,4) = 1e-9;
    J(8,8) = -1e-9 * CAPD;

    % Row 9: f(9) = 1e-9 * (y(5) - qdrain(vgb, vgs, vgd))
    J(9,5) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqdrain_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS);
    J(9,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(9,7) = -1e-9 * (-dq_dvgs);
    J(9,8) = -1e-9 * (-dq_dvgd);
end

%% --- Smoothed Jacobian of qgate -----------------------------------------
function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqgate_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    swapped = (vgs - vgd) <= 0;
    if swapped
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    sqrt_phi_ubs = sqrt(PHI - ubs);
    vte = VT0 + GAMMA*(sqrt_phi_ubs - sqrt(PHI));

    % dvte/dubs = -GAMMA/(2*sqrt(PHI-ubs))
    dvte_dubs = -GAMMA / (2*sqrt_phi_ubs);

    % --- Region A derivatives ---
    dqA_dugb = COX;
    dqA_dugs = 0;
    dqA_dugd = 0;

    % --- Region B derivatives ---
    arg_B = (GAMMA/2)^2 + ugb - vfb;
    arg_B = max(arg_B, 1e-30);
    dqB_dugb = COX * GAMMA / (2*sqrt(arg_B));
    dqB_dugs = 0;
    dqB_dugd = 0;

    % --- Region C derivatives (with smoothed ugst and ugdt) ---
    sp_arg_s = ugs - vte;
    ugst = softplus(sp_arg_s, EPS);
    dsp_s = dsoftplus(sp_arg_s, EPS);  % d(ugst)/d(ugs-vte)
    sp_arg_d = ugd - vte;
    ugdt = softplus(sp_arg_d, EPS);
    dsp_d = dsoftplus(sp_arg_d, EPS);  % d(ugdt)/d(ugd-vte)

    S = ugst + ugdt;
    if S > 1e-30
        H = ugdt + ugst - (ugdt*ugst)/S;
        dH_dugst = 1 - (ugdt^2)/(S^2);
        dH_dugdt = 1 - (ugst^2)/(S^2);
    else
        H = 0;
        dH_dugst = 1;
        dH_dugdt = 0;
    end

    % Chain rules for ugst and ugdt w.r.t. (ugb, ugs, ugd)
    % ugst = softplus(ugs - vte, EPS); d(ugs-vte)/dugs = 1 - dvte_dubs
    % d(ugs-vte)/dugb = dvte_dubs (since dubs/dugb = -1)
    dugst_dugs = dsp_s * (1 - dvte_dubs);
    dugst_dugb = dsp_s * dvte_dubs;
    % ugdt = softplus(ugd - vte, EPS)
    % d(ugd-vte)/dugd = 1, d(ugd-vte)/dugs = -dvte_dubs, d(ugd-vte)/dugb = dvte_dubs
    dugdt_dugd = dsp_d * 1;
    dugdt_dugs = dsp_d * (-dvte_dubs);
    dugdt_dugb = dsp_d * dvte_dubs;

    dsqrt_dugs = -GAMMA / (2*sqrt_phi_ubs);  % d(GAMMA*sqrt(PHI-ubs))/dugs
    dsqrt_dugb =  GAMMA / (2*sqrt_phi_ubs);  % d(GAMMA*sqrt(PHI-ubs))/dugb

    dqC_dugs = COX * ((2/3)*(dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs) + dsqrt_dugs);
    dqC_dugb = COX * ((2/3)*(dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb) + dsqrt_dugb);
    dqC_dugd = COX * (2/3) * dH_dugdt * dugdt_dugd;

    % --- Smooth blending of B/C boundary ---
    sigma_BC = sigmoid(ugs - vte, EPS);
    dsigma_BC = dsigmoid(ugs - vte, EPS);
    % d(sigma_BC)/dugs = dsigma_BC * d(ugs-vte)/dugs = dsigma_BC * (1 - dvte_dubs)
    % d(sigma_BC)/dugb = dsigma_BC * d(ugs-vte)/dugb = dsigma_BC * dvte_dubs
    % d(sigma_BC)/dugd = dsigma_BC * d(ugs-vte)/dugd = 0 (vte doesn't depend on ugd)
    %   Actually vte depends on ubs = ugs - ugb, so no ugd dependence.
    dsigma_dugs = dsigma_BC * (1 - dvte_dubs);
    dsigma_dugb = dsigma_BC * dvte_dubs;
    dsigma_dugd = 0;

    if ugb <= vfb
        % Region A (hard boundary, already C1)
        dq_dugb = dqA_dugb;
        dq_dugs = dqA_dugs;
        dq_dugd = dqA_dugd;
    else
        % Blended: q = (1-sigma)*q_B + sigma*q_C
        q_B = COX * GAMMA * (sqrt(arg_B) - GAMMA/2);
        ugst_val = softplus(ugs - vte, EPS);
        ugdt_val = softplus(ugd - vte, EPS);
        S_val = ugst_val + ugdt_val;
        if S_val > 1e-30
            H_val = ugdt_val + ugst_val - (ugdt_val*ugst_val)/S_val;
        else
            H_val = 0;
        end
        q_C = COX * ((2/3)*H_val + GAMMA*sqrt_phi_ubs);

        % d/dx [(1-sigma)*q_B + sigma*q_C]
        % = -dsigma*q_B + (1-sigma)*dqB + dsigma*q_C + sigma*dqC
        % = dsigma*(q_C - q_B) + (1-sigma)*dqB + sigma*dqC
        dq_dugb = dsigma_dugb*(q_C - q_B) + (1-sigma_BC)*dqB_dugb + sigma_BC*dqC_dugb;
        dq_dugs = dsigma_dugs*(q_C - q_B) + (1-sigma_BC)*dqB_dugs + sigma_BC*dqC_dugs;
        dq_dugd = dsigma_dugd*(q_C - q_B) + (1-sigma_BC)*dqB_dugd + sigma_BC*dqC_dugd;
    end

    % Map back from (ugb, ugs, ugd) to (vgb, vgs, vgd)
    dq_dvgb = dq_dugb;
    if swapped
        dq_dvgs = dq_dugd;
        dq_dvgd = dq_dugs;
    else
        dq_dvgs = dq_dugs;
        dq_dvgd = dq_dugd;
    end
end

%% --- Smoothed Jacobian of qsrc ------------------------------------------
function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    swapped = (vgs - vgd) <= 0;
    if swapped
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    sqrt_phi_ubs = sqrt(PHI - ubs);
    vte = VT0 + GAMMA*(sqrt_phi_ubs - sqrt(PHI));
    dvte_dubs = -GAMMA / (2*sqrt_phi_ubs);

    % Region C: q_C = -COX*(1/3)*H with smoothed ugst and ugdt
    sp_arg_s = ugs - vte;
    ugst = softplus(sp_arg_s, EPS);
    dsp_s = dsoftplus(sp_arg_s, EPS);
    sp_arg_d = ugd - vte;
    ugdt = softplus(sp_arg_d, EPS);
    dsp_d = dsoftplus(sp_arg_d, EPS);

    S = ugst + ugdt;
    if S > 1e-30
        H = ugdt + ugst - (ugdt*ugst)/S;
        dH_dugst = 1 - (ugdt^2)/(S^2);
        dH_dugdt = 1 - (ugst^2)/(S^2);
    else
        H = 0;
        dH_dugst = 1;
        dH_dugdt = 0;
    end

    dugst_dugs = dsp_s * (1 - dvte_dubs);
    dugst_dugb = dsp_s * dvte_dubs;
    dugdt_dugd = dsp_d * 1;
    dugdt_dugs = dsp_d * (-dvte_dubs);
    dugdt_dugb = dsp_d * dvte_dubs;

    dqC_dugs = -COX*(1/3)*(dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs);
    dqC_dugb = -COX*(1/3)*(dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb);
    dqC_dugd = -COX*(1/3)*dH_dugdt*dugdt_dugd;

    q_C = -COX*(1/3)*H;

    % Smooth blending: q = sigma * q_C (since q_AB = 0)
    sigma_BC = sigmoid(ugs - vte, EPS);
    dsigma_BC = dsigmoid(ugs - vte, EPS);
    dsigma_dugs = dsigma_BC * (1 - dvte_dubs);
    dsigma_dugb = dsigma_BC * dvte_dubs;
    dsigma_dugd = 0;

    % d/dx [sigma * q_C] = dsigma * q_C + sigma * dqC
    dq_dugb = dsigma_dugb * q_C + sigma_BC * dqC_dugb;
    dq_dugs = dsigma_dugs * q_C + sigma_BC * dqC_dugs;
    dq_dugd = dsigma_dugd * q_C + sigma_BC * dqC_dugd;

    % Map back
    dq_dvgb = dq_dugb;
    if swapped
        dq_dvgs = dq_dugd;
        dq_dvgd = dq_dugs;
    else
        dq_dvgs = dq_dugs;
        dq_dvgd = dq_dugd;
    end
end

function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqdrain_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS)
    % qdrain is identical to qsrc in this model
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc_smooth(vgb, vgs, vgd, VT0, GAMMA, PHI, COX, EPS);
end


end
