function T = design_study(P, doPrint)
%DESIGN_STUDY  Select the estimator parameters (q, alpha) on a tuning ensemble.
%
%   T = design_study(P) sweeps the jerk PSD q and the acceleration decay rate
%   alpha of the Singer model over a grid and evaluates the 0.5 s-ahead
%   elevation-prediction error with estimator + predictor only (no turret).
%   The ensemble is 16 random bounded manoeuvres (random_scenario idx 101..116,
%   odd idx include a sensor loss), DISJOINT from scenarios A-C and from the
%   robustness-study cases.  The chosen pair is then frozen in parameters.m.
%
%   T: struct with fields q, alpha (grids), rms (nq x na, deg), p95 (deg).
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if nargin < 2, doPrint = true; end
qg = [15 60 200 600 1500];  ag = [0 0.5 1 2 4];
ncase = 16;  Ts = P.Ts_radar;  K = round(P.t_end/Ts);

% -------- pre-extract measurement sequences and truths -------------------
D = cell(ncase,1);
for c = 1:ncase
    cfg = random_scenario(100 + c, P, mod(c,2) == 1);
    d.t = (0:K)'*Ts; d.rm = zeros(K+1,1); d.thm = d.rm; d.v = false(K+1,1);
    d.thf = d.rm;
    for k = 0:K
        [d.rm(k+1), d.thm(k+1), d.v(k+1)] = target_motion_generator(k*Ts, cfg, 100 + c);
        [~,~,~,tf] = target_motion_generator(k*Ts + P.Tp, cfg, 100 + c);
        d.thf(k+1) = tf.theta;
    end
    D{c} = d;
end

% -------- sweep ----------------------------------------------------------
T.q = qg; T.alpha = ag; T.rms = zeros(numel(qg), numel(ag)); T.p95 = T.rms;
for iq = 1:numel(qg)
  for ia = 1:numel(ag)
    Q = P;  Q.est.q = qg(iq);  Q.est.alpha = ag(ia);
    err = [];
    for c = 1:ncase
        d = D{c};  S = [];
        for k = 0:K
            if d.v(k+1), S = state_estimator(S, d.t(k+1), d.rm(k+1), d.thm(k+1), Q); end
            if k*Ts >= P.sim.t_settle
                pr = target_predictor(S, d.t(k+1), Q);
                err(end+1) = atan2(sin(pr.theta_p - d.thf(k+1)), cos(pr.theta_p - d.thf(k+1)))*180/pi;
            end
        end
    end
    T.rms(iq,ia) = sqrt(mean(err.^2));  ea = sort(abs(err)); T.p95(iq,ia) = ea(ceil(0.95*numel(ea)));
  end
end
if doPrint
    fprintf('\nEstimator design study (16-case ensemble): RMS 0.5 s-ahead elevation error [deg]\n');
    fprintf('%8s', 'q\alpha'); fprintf('%10.2f', ag); fprintf('\n');
    for iq = 1:numel(qg), fprintf('%8g', qg(iq)); fprintf('%10.4f', T.rms(iq,:)); fprintf('\n'); end
    [m, i] = min(T.rms(:)); [iq, ia] = ind2sub(size(T.rms), i);
    fprintf('Minimum %.4f deg at q = %g, alpha = %g\n\n', m, qg(iq), ag(ia));
end
end
