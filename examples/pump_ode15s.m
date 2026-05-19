%% pump_ode15s.m — Charge Pump (DAE index-2, n=9) solved with ode15s
%
% MOS transistor charge pump circuit simulation.
% M*y' = f(t,y) where M is a 9x9 singular mass matrix (5 nonzeros)
%
% State vector:
%   y(1:5) = charges (qgate, caps*Vs, qsrc, capd*Vd, qdrain)
%   y(6:8) = node voltages (Vg, Vs, Vd)
%   y(9)   = current
%
% DAE index: components 1..8 are index-1, component 9 is index-2
%
% t in [0, 1200e-9] with 39 discontinuities (pulsed input vin(t))
% Reference solution from IVPtestset (GAMD quadruple precision)

clear; clc;

%% Parameters
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

%% Mass matrix (constant, singular)
M = zeros(9,9);
M(1,1) = 1;
M(2,2) = 1;
M(2,3) = 1;
M(3,4) = 1;
M(3,5) = 1;

%% Initial conditions
y0 = zeros(9,1);
y0(1) = fn_qgate(0, 0, 0, VT0, GAMMA, PHI, COX);
y0(3) = fn_qsrc(0, 0, 0, VT0, GAMMA, PHI, COX);
y0(5) = fn_qdrain(0, 0, 0, VT0, GAMMA, PHI, COX);

%% Tolerances (component-dependent, matching Fortran)
Tol = 1e-7;
rtol_v = Tol * ones(9,1);
atol_v = [Tol*1e-6*ones(5,1); Tol*ones(3,1); Tol];

%% Discontinuity times
disc = zeros(41,1);
disc(1) = 0;
disc(2) = 50;
disc(3) = 60;
disc(4) = 110;
for i = 5:40
    disc(i) = disc(mod(i-1,4)+1) + floor((i-1)/4)*120;
end
disc(41) = 1200;

%% ode15s options
opts = odeset('Mass', M, ...
              'MassSingular', 'yes', ...
              'RelTol', rtol_v, ...
              'AbsTol', atol_v, ...
              'InitialStep', 1, ...
              'MaxStep', 200, ...
              'Jacobian', @(t,y) jac_pump(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS));

%% Segmented integration (restart at each discontinuity)
y = y0;
t_cur = 0;
total_steps = 0;

fprintf('Solving charge pump with ode15s...\n');
for seg = 1:40
    tspan = [disc(seg), disc(seg+1)];
    if tspan(2) <= tspan(1), continue; end

    odefun = @(t, y) rhs_pump(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, ...
                               VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS);

    [tsol, ysol] = ode15s(odefun, tspan, y, opts);

    y = ysol(end,:)';
    total_steps = total_steps + length(tsol) - 1;
    t_cur = tsol(end);
end

%% Reference solution at t = 1200e-9
yref = zeros(9,1);
yref(1) = 0.126280042987675933170893e-12;
yref(9) = 0.152255686815577679043511e-03;

%% Results
fprintf('\n=== Charge Pump (ode15s, Tol=%.1e) ===\n', Tol);
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

function f = rhs_pump(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS, ...
                      VHIGH, DELTAT, T1_PLS, T2_PLS, T3_PLS)
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
    f(5) = 1e-9 * (v(1) - fn_qgate(vgb, vgs, vgd, VT0, GAMMA, PHI, COX));
    f(6) = 1e-9 * (v(2) - CAPS * v(7));
    f(7) = 1e-9 * (v(3) - fn_qsrc(vgb, vgs, vgd, VT0, GAMMA, PHI, COX));
    f(8) = 1e-9 * (v(4) - CAPD * v(8));
    f(9) = 1e-9 * (v(5) - fn_qdrain(vgb, vgs, vgd, VT0, GAMMA, PHI, COX));
end

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

function q = fn_qgate(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    if (vgs - vgd) <= 0
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    if ugb <= vfb
        q = COX * (ugb - vfb);
    elseif ugs <= vte
        q = COX * GAMMA * (sqrt((GAMMA/2)^2 + ugb - vfb) - GAMMA/2);
    else
        ugst = ugs - vte;
        ugdt = max(ugd - vte, 0);
        q = COX * ((2/3)*(ugdt + ugst - (ugdt*ugst)/(ugdt+ugst)) + ...
                   GAMMA*sqrt(PHI - ubs));
    end
end

function q = fn_qsrc(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    if (vgs - vgd) <= 0
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    if ugb <= vfb
        q = 0;
    elseif ugs <= vte
        q = 0;
    else
        ugst = ugs - vte;
        ugdt = max(ugd - vte, 0);
        q = -COX * (1/3) * (ugdt + ugst - (ugdt*ugst)/(ugdt+ugst));
    end
end

function q = fn_qdrain(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    if (vgs - vgd) <= 0
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    if ugb <= vfb
        q = 0;
    elseif ugs <= vte
        q = 0;
    else
        ugst = ugs - vte;
        ugdt = max(ugd - vte, 0);
        q = -COX * (1/3) * (ugdt + ugst - (ugdt*ugst)/(ugdt+ugst));
    end
end

function J = jac_pump(t, y, VT0, GAMMA, PHI, COX, CAPD, CAPS)
    % Analytic Jacobian: J(i,j) = df(i)/dy(j)
    %
    % Chain rule for charge functions called as Q(y(6), y(6)-y(7), y(6)-y(8)):
    %   dQ/dy(6) = dQ/dvgb + dQ/dvgs + dQ/dvgd
    %   dQ/dy(7) = -dQ/dvgs
    %   dQ/dy(8) = -dQ/dvgd

    J = zeros(9,9);

    vgb = y(6);
    vgs = y(6) - y(7);
    vgd = y(6) - y(8);

    % Row 1: f(1) = -1e-9 * y(9)
    J(1,9) = -1e-9;

    % Row 2: f(2) = 0
    % Row 3: f(3) = 0

    % Row 4: f(4) = 1e-9 * (-y(6) + vin(t))
    J(4,6) = -1e-9;

    % Row 5: f(5) = 1e-9 * (y(1) - qgate(vgb, vgs, vgd))
    J(5,1) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqgate(vgb, vgs, vgd, VT0, GAMMA, PHI, COX);
    J(5,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(5,7) = -1e-9 * (-dq_dvgs);
    J(5,8) = -1e-9 * (-dq_dvgd);

    % Row 6: f(6) = 1e-9 * (y(2) - CAPS * y(7))
    J(6,2) = 1e-9;
    J(6,7) = -1e-9 * CAPS;

    % Row 7: f(7) = 1e-9 * (y(3) - qsrc(vgb, vgs, vgd))
    J(7,3) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc(vgb, vgs, vgd, VT0, GAMMA, PHI, COX);
    J(7,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(7,7) = -1e-9 * (-dq_dvgs);
    J(7,8) = -1e-9 * (-dq_dvgd);

    % Row 8: f(8) = 1e-9 * (y(4) - CAPD * y(8))
    J(8,4) = 1e-9;
    J(8,8) = -1e-9 * CAPD;

    % Row 9: f(9) = 1e-9 * (y(5) - qdrain(vgb, vgs, vgd))
    J(9,5) = 1e-9;
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqdrain(vgb, vgs, vgd, VT0, GAMMA, PHI, COX);
    J(9,6) = -1e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    J(9,7) = -1e-9 * (-dq_dvgs);
    J(9,8) = -1e-9 * (-dq_dvgd);
end

function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqgate(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    % Analytic derivatives of qgate w.r.t. (vgb, vgs, vgd)
    swapped = (vgs - vgd) <= 0;
    if swapped
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    if ugb <= vfb
        % Region A: qgate = COX*(ugb - vfb)
        dq_dugb = COX;
        dq_dugs = 0;
        dq_dugd = 0;
    elseif ugs <= vte
        % Region B: qgate = COX*GAMMA*(sqrt(G^2/4 + ugb - vfb) - G/2)
        arg = (GAMMA/2)^2 + ugb - vfb;
        dq_dugb = COX * GAMMA / (2*sqrt(arg));
        dq_dugs = 0;
        dq_dugd = 0;
    else
        % Region C (inversion)
        sqrt_phi_ubs = sqrt(PHI - ubs);
        dvte_dubs = -GAMMA / (2*sqrt_phi_ubs);

        ugst = ugs - vte;
        ugd_above_vte = (ugd > vte);
        if ugd_above_vte
            ugdt = ugd - vte;
        else
            ugdt = 0;
        end

        S = ugst + ugdt;
        if S > 0
            dH_dugst = 1 - (ugdt^2)/(S^2);
            dH_dugdt = 1 - (ugst^2)/(S^2);
        else
            dH_dugst = 1;
            dH_dugdt = 0;
        end

        dugst_dugs = 1 - dvte_dubs;
        dugst_dugb = dvte_dubs;

        if ugd_above_vte
            dugdt_dugd = 1;
            dugdt_dugs = -dvte_dubs;
            dugdt_dugb = dvte_dubs;
        else
            dugdt_dugd = 0;
            dugdt_dugs = 0;
            dugdt_dugb = 0;
        end

        dsqrt_dugs = -GAMMA / (2*sqrt_phi_ubs);
        dsqrt_dugb =  GAMMA / (2*sqrt_phi_ubs);

        dq_dugs = COX * ((2/3)*(dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs) + dsqrt_dugs);
        dq_dugb = COX * ((2/3)*(dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb) + dsqrt_dugb);
        dq_dugd = COX * (2/3) * dH_dugdt * dugdt_dugd;
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

function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    % Analytic derivatives of qsrc w.r.t. (vgb, vgs, vgd)
    swapped = (vgs - vgd) <= 0;
    if swapped
        ugs = vgd; ugd = vgs;
    else
        ugs = vgs; ugd = vgd;
    end
    ugb = vgb;
    ubs = ugs - ugb;

    vfb = VT0 - GAMMA*sqrt(PHI) - PHI;
    vte = VT0 + GAMMA*(sqrt(PHI - ubs) - sqrt(PHI));

    if ugb <= vfb || ugs <= vte
        % Regions A and B: qsrc = 0
        dq_dvgb = 0; dq_dvgs = 0; dq_dvgd = 0;
        return;
    end

    % Region C: qsrc = -COX*(1/3)*H
    sqrt_phi_ubs = sqrt(PHI - ubs);
    dvte_dubs = -GAMMA / (2*sqrt_phi_ubs);

    ugst = ugs - vte;
    ugd_above_vte = (ugd >= vte);
    if ugd_above_vte
        ugdt = ugd - vte;
    else
        ugdt = 0;
    end

    S = ugst + ugdt;
    if S > 0
        dH_dugst = 1 - (ugdt^2)/(S^2);
        dH_dugdt = 1 - (ugst^2)/(S^2);
    else
        dH_dugst = 1;
        dH_dugdt = 0;
    end

    dugst_dugs = 1 - dvte_dubs;
    dugst_dugb = dvte_dubs;

    if ugd_above_vte
        dugdt_dugd = 1;
        dugdt_dugs = -dvte_dubs;
        dugdt_dugb = dvte_dubs;
    else
        dugdt_dugd = 0;
        dugdt_dugs = 0;
        dugdt_dugb = 0;
    end

    dq_dugs = -COX*(1/3)*(dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs);
    dq_dugb = -COX*(1/3)*(dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb);
    dq_dugd = -COX*(1/3)*dH_dugdt*dugdt_dugd;

    dq_dvgb = dq_dugb;
    if swapped
        dq_dvgs = dq_dugd;
        dq_dvgd = dq_dugs;
    else
        dq_dvgs = dq_dugs;
        dq_dvgd = dq_dugd;
    end
end

function [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqdrain(vgb, vgs, vgd, VT0, GAMMA, PHI, COX)
    % qdrain is identical to qsrc in this model
    [dq_dvgb, dq_dvgs, dq_dvgd] = fn_dqsrc(vgb, vgs, vgd, VT0, GAMMA, PHI, COX);
end
