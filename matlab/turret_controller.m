function [V_m, fire_cmd, dbg] = turret_controller(t, r_meas, theta_meas, ...
                                       meas_valid, theta_turret, theta_dot_turret)
%TURRET_CONTROLLER  Predictive-feedback controller for the tracking platform.
%
%   [V_m, fire_cmd]      = turret_controller(t, r_meas, theta_meas, meas_valid, ...
%                                            theta_turret, theta_dot_turret)
%   [V_m, fire_cmd, dbg] = ...   (optional diagnostics for plotting)
%
%   Signal flow (paper Fig. 2):
%     radar (r, theta) -> EKF state estimate -> open-loop prediction at t+Tp
%       -> look-ahead elevation reference theta_ref -> PI compensator -> V_m
%
%   * The EKF (state_estimator.m) runs only when meas_valid = 1.
%   * target_predictor.m propagates the estimate to the present time and to
%     t + Tp on every call; during a measurement loss it simply keeps
%     extrapolating the last estimate.
%   * Feedback law (LTI compensator, Gc(s) = Kp (s + wi)/s):
%         e = theta_ref - theta_turret
%         V_m = Kp ( e + wi * integral(e) )
%     implemented with a trapezoidal integrator, conditional-integration
%     anti-windup and a supply-voltage limit.
%   * fire_cmd: one-sample 0->1 pulse issued by the gated firing logic below
%     (recharge, estimated miss distance, EKF uncertainty).
%
%   The persistent state is reset whenever t == 0 (start of a run).
%   theta_dot_turret is accepted for interface compatibility; the LTI
%   compensator acts on the angle error only.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

persistent Pc S I e_prev t_prev fire_prev t_fire_last

if isempty(Pc), Pc = parameters(); end
if t <= 0 || isempty(t_prev)                  % ---- start of a new run ----
    S = [];  I = 0;  e_prev = 0;  t_prev = t;  fire_prev = 0;  t_fire_last = -Inf;
end
dt = t - t_prev;
if dt < 0, error('turret_controller:timeReversal', 'Time must not decrease.'); end

% -------------------------------------------------- estimation (on demand)
if meas_valid
    if ~(isfinite(r_meas) && isfinite(theta_meas) && r_meas > 0)
        error('turret_controller:badMeasurement', 'Non-finite or non-positive measurement.');
    end
    S = state_estimator(S, t, r_meas, theta_meas, Pc);
end

% ---------------------------------------------------------- prediction ---
if isempty(S) || ~S.init                       % nothing known yet: hold position
    theta_ref = theta_turret;
    pred = [];
else
    pred = target_predictor(S, t, Pc);
    theta_ref = pred.theta_ref;
end

% ------------------------------------------------------ PI compensator ----
e  = theta_ref - theta_turret;
V_unsat = Pc.ctrl.Kp * (e + Pc.ctrl.wi * I);
V_m = min(max(V_unsat, -Pc.ctrl.V_sat), Pc.ctrl.V_sat);
% conditional integration (anti-windup): freeze I if saturated and e drives further
if V_unsat == V_m || sign(e) ~= sign(V_unsat)
    I = I + 0.5*(e + e_prev)*dt;               % trapezoidal rule
end
e_prev = e;  t_prev = t;

% ---------------------------------------------------------- firing logic --
%   Fire (0->1 edge, one-sample pulse) only when ALL of:
%     (a) recharge elapsed:   t - t_last_pulse >= t_recharge
%     (b) filter armed:       t >= t_arm and the estimator is initialised
%     (c) predicted target is inside the mechanical range (reference not clipped)
%     (d) estimated miss  m_hat = r_p |sin(theta_p - theta_turret)| <= k_gate*R_eng
%     (e) prediction uncertainty r_p*sigma_theta_p <= k_sig*R_eng
%   (e) uses the EKF covariance, so it automatically blocks firing when the
%   estimate has degraded (start-up, long radar loss).
fire_cmd = 0;
if ~isempty(pred) && t >= Pc.fire.t_arm && (t - t_fire_last) >= Pc.t_recharge - 1e-9 ...
        && ~fire_prev && pred.theta_p >= Pc.theta_min && pred.theta_p <= Pc.theta_max
    dth = atan2(sin(pred.theta_p - theta_turret), cos(pred.theta_p - theta_turret));
    m_hat = pred.r_p * abs(sin(dth));
    if abs(dth) < pi/2 && m_hat <= Pc.fire.k_gate * Pc.R_eng
        pc = target_predictor(S, t, Pc, true);
        if pred.r_p * pc.sig_theta_p <= Pc.fire.k_sig * Pc.R_eng
            fire_cmd = 1;  t_fire_last = t;
        end
    end
end
fire_prev = fire_cmd;

if nargout >= 3
    dbg.theta_ref = theta_ref;  dbg.e = e;  dbg.I = I;
    dbg.V_unsat = V_unsat;      dbg.pred = pred;
    dbg.S = S;
end
end
