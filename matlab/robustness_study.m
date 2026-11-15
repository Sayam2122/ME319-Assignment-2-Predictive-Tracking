function R = robustness_study(P, ncase)
%ROBUSTNESS_STUDY  Generalisation test of the FROZEN design (no retuning).
%
%   R = robustness_study(P, ncase)
%
%   Part 1 (full closed loop): ncase random bounded manoeuvres
%       (random_scenario idx 201..200+ncase, even idx contain a sensor loss
%       with onset >= 1 s and duration <= 0.5 s), each simulated with
%       simulation_runner and scored with performance_metrics.
%   Part 2 (estimator + predictor only, cheap): sweep of the loss onset
%       t_bo0 = 1:0.5:9 s (duration 0.5 s) on the scenario-B manoeuvre and
%       compare the worst 0.5 s-ahead angle error inside the loss window
%       with the loss-free run of the same manoeuvre.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if nargin < 2, ncase = 10; end
R.ncase = ncase;
fn = {'look_rms','look_max','pred_ang_rms','pred_ang_max','trk_rms','trk_max','est_pos_rms'};
R.tab = zeros(ncase, numel(fn)+3);  R.eng = zeros(ncase, 7);   % success, t_int, miss_hit, n_pulses, success_now, rmax  R.with_loss = false(ncase,1);
for c = 1:ncase
    idx = 200 + c;  wl = mod(c,2) == 0;
    cfg = random_scenario(idx, P, wl);
    res = simulation_runner(cfg, P, idx);
    M   = performance_metrics(res, P);
    for j = 1:numel(fn), R.tab(c,j) = M.post.(fn{j}); end
    R.tab(c,numel(fn)+1) = M.V_peak;  R.tab(c,numel(fn)+2) = M.I_peak;
    R.tab(c,numel(fn)+3) = M.sat_frac;
    R.with_loss(c) = wl;
    e = M.eng;  R.eng(c,:) = [e.success, e.t_int, e.miss_hit, e.n_pulses, e.success_now, mean(res.rm), e.hit_frac];
    fprintf('     engagement: success=%d  t_int=%.2f s  miss=%.2f m  pulses=%d (scored at t_fire: %d)\n', ...
            e.success, e.t_int, e.miss_hit, e.n_pulses, e.success_now);
    fprintf('  case %2d (loss=%d): look rms %.3f max %.3f deg | pred max %.3f | Vpk %.0f V\n', ...
            c, wl, M.post.look_rms, M.post.look_max, M.post.pred_ang_max, M.V_peak);
end
R.fields = [fn, {'V_peak','I_peak','sat_frac'}];

% ---------------- Part 2: loss-onset sweep (estimator + predictor) --------
Ts = P.Ts_radar;  K = round(P.t_end/Ts);
onset = (1:0.5:9)';  R.onset = onset;
R.sweep_loss = zeros(size(onset));  R.sweep_ref = zeros(size(onset));
cfgB = P.scen.B;
for io = 1:numel(onset)
    cfg = cfgB;  cfg.t_bo0 = onset(io);  cfg.t_bo = P.t_bo;
    e_loss = -inf;  e_ref = -inf;
    for pass = 1:2                                % pass 1: with loss, pass 2: loss-free
        cc = cfg; if pass == 2, cc.t_bo0 = NaN; cc.t_bo = 0; end
        S = [];
        for k = 0:K
            t = k*Ts;
            [rm, thm, v, ~] = target_motion_generator(t, cc, P.seed);
            if v, S = state_estimator(S, t, rm, thm, P); end
            if t >= onset(io) && t <= onset(io) + P.t_bo + P.Tp
                pr = target_predictor(S, t, P);
                [~,~,~,tf] = target_motion_generator(t + P.Tp, cc, P.seed);
                e = abs(atan2(sin(pr.theta_p - tf.theta), cos(pr.theta_p - tf.theta)))*180/pi;
                if pass == 1, e_loss = max(e_loss, e); else, e_ref = max(e_ref, e); end
            end
        end
    end
    R.sweep_loss(io) = e_loss;  R.sweep_ref(io) = e_ref;
end
fprintf('  loss-onset sweep: worst window error %.3f deg (loss) vs %.3f deg (no loss)\n', ...
        max(R.sweep_loss), max(R.sweep_ref));
end
