function M = performance_metrics(res, P)
%PERFORMANCE_METRICS  Quantitative metrics of one simulation result.
%
%   M = performance_metrics(res, P)
%
%   All angles in degrees, positions in metres, speeds in m/s, V in volts.
%   "full" = whole interval 0 <= t <= t_end; "post" = t >= P.sim.t_settle
%   (after the estimator's initial acquisition transient).
%
%   Estimation   : est_pos_rms/max   |p_hat(t) - p(t)|            [m]
%                  est_vel_rms/max   |v_hat(t) - v(t)|            [m/s]
%                  raw_pos_rms       |converted radar sample - p|  [m] (unfiltered)
%   Prediction   : pred_ang_rms/max  theta_p(t) - theta_true(t+Tp) [deg]
%                  pred_pos_rms/max  |p_pred(t) - p_true(t+Tp)|     [m]
%   Control      : trk_rms/max       theta_ref - theta_turret      [deg]
%   Look-ahead   : look_rms/max      theta_true(t+Tp) - theta_turret [deg]
%                  (combined prediction + control pointing error)
%   Actuation    : V_peak, I_peak, om_peak, sat_frac (fraction of time at V_sat)
%   Loss (C only): window t in [t_bo0, t_bo0 + t_bo + Tp]; recovery time =
%                  time after measurement return until the position-estimate
%                  error stays below 3*sigma_r.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

r2d = 180/pi;
t   = res.t;
full = true(size(t));
post = t >= P.sim.t_settle - 1e-12;

ep   = hypot(res.xh - res.tx,  res.yh - res.ty);
ev   = hypot(res.vxh - res.tvx, res.vyh - res.tvy);
% raw (unfiltered) Cartesian conversion error at the instants of valid samples
vs   = res.valid > 0;
xr   = res.rm .* cos(res.thm);  yr = res.rm .* sin(res.thm);
er   = hypot(xr - res.tx, yr - res.ty);
pa   = wrap(res.thp - res.fth) * r2d;
pp   = hypot(res.xp - res.fx, res.yp - res.fy);
tk   = wrap(res.th_ref - res.th) * r2d;
lk   = wrap(res.fth - res.th) * r2d;
in   = wrap(res.tth - res.th) * r2d;

for tag = {'full','post'}
    k = tag{1};
    if strcmp(k,'full'), m = full; else, m = post; end
    M.(k).est_pos_rms = rms_(ep(m));    M.(k).est_pos_max = max(ep(m));
    M.(k).est_vel_rms = rms_(ev(m));    M.(k).est_vel_max = max(ev(m));
    M.(k).raw_pos_rms = rms_(er(m & vs));
    M.(k).pred_ang_rms = rms_(pa(m));   M.(k).pred_ang_max = max(abs(pa(m)));
    M.(k).pred_pos_rms = rms_(pp(m));   M.(k).pred_pos_max = max(pp(m));
    M.(k).trk_rms  = rms_(tk(m));       M.(k).trk_max  = max(abs(tk(m)));
    M.(k).look_rms = rms_(lk(m));       M.(k).look_max = max(abs(lk(m)));
    M.(k).inst_rms = rms_(in(m));       M.(k).inst_max = max(abs(in(m)));
end
% ---------------- laser engagement ---------------------------------------
pl = res.pulses;                                   % [t theta miss_now miss_future ok_future ok_now]
M.eng.n_pulses = size(pl,1);
M.eng.success  = any(pl(:,5));                     % scored at t_fire + Tp (primary)
M.eng.success_now = any(pl(:,6));                  % scored at t_fire (alternative reading)
if M.eng.success, M.eng.t_int = pl(find(pl(:,5),1),1); else, M.eng.t_int = NaN; end
if M.eng.n_pulses > 0
    M.eng.t_first_pulse = pl(1,1);  M.eng.miss_first = pl(1,4);
    if M.eng.success, M.eng.miss_hit = pl(find(pl(:,5),1),4); else, M.eng.miss_hit = NaN; end
    M.eng.hit_frac = mean(pl(:,5));
else
    M.eng.t_first_pulse = NaN;  M.eng.miss_first = NaN;  M.eng.miss_hit = NaN;  M.eng.hit_frac = NaN;
end
M.V_peak  = max(abs(res.V));
M.I_peak  = max(abs(res.I));
M.om_peak = max(abs(res.om)) * r2d;                         % deg/s
M.sat_frac = mean(abs(res.V) >= P.ctrl.V_sat - 1e-9);
M.th_range_deg = [min(res.th) max(res.th)] * r2d;
M.nis_mean = mean(res.nis(~isnan(res.nis)));

% ---------------- measurement-loss window (meaningful for scenario C) -----
cfg = [];
if ischar(res.scenario) && isfield(P.scen, upper(res.scenario)), cfg = P.scen.(upper(res.scenario)); end
if ~isempty(cfg) && isfinite(cfg.t_bo0) && cfg.t_bo > 0
    t1 = cfg.t_bo0;  t2 = cfg.t_bo0 + cfg.t_bo;  t3 = t2 + P.Tp;
    w = t >= t1 & t <= t3;
    M.loss.window = [t1 t3];
    M.loss.pred_ang_max = max(abs(pa(w)));
    M.loss.pred_pos_max = max(pp(w));
    M.loss.est_pos_max  = max(ep(w));
    M.loss.trk_max      = max(abs(tk(w)));
    M.loss.look_max     = max(abs(lk(w)));
    M.loss.sig_thp_max  = max(res.sig_thp(w & ~isnan(res.sig_thp))) * r2d;
    bad = find(ep > 3*P.sigma_r & t > t2, 1, 'last');
    if isempty(bad), M.loss.recovery_time = 0;
    else, M.loss.recovery_time = max(0, t(bad) - t2); end
end
end

function y = rms_(x), x = x(~isnan(x)); y = sqrt(mean(x.^2)); end
function a = wrap(a), a = atan2(sin(a), cos(a)); end
