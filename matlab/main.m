function main(varargin)
%MAIN  Run the complete ME319 Assignment 2 project end to end.
%
%   main            full run (checks, design study, scenarios A-C, robustness
%                   study with 10 random cases, figures, report data);
%                   takes about 6-8 minutes in Octave.
%   main('reuse')   re-uses the saved design/gate studies but runs all 10 random cases.
%   main('quick')   skips the estimator design study (re-uses
%                   results/data/design_study.mat) and uses 4 random cases.
%
%   Outputs: results/figures/*.png|pdf, results/data/*.mat|csv,
%            (no LaTeX is written: the paper is written by hand from these files)
%
%   Project : ME319 Control Systems (Autumn 2026) - Assignment 2
%   Authors : Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), IIT Bombay
%   Scope   : non-weaponised predictive tracking-platform simulation.

quick = any(strcmpi(varargin, 'quick'));
reuse = quick || any(strcmpi(varargin, 'reuse'));   % 'reuse': keep saved design/gate studies, still run 10 random cases
here  = fileparts(mfilename('fullpath'));
root  = fileparts(here);
addpath(here);
if exist('OCTAVE_VERSION','builtin') ~= 0, pkg load control; end
figdir = fullfile(root, 'results', 'figures');
datadir = fullfile(root, 'results', 'data');
for d = {figdir, datadir}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end

P = parameters();
fprintf('ME319 Assignment 2 - predictive tracking platform (team: see README)\n');

% 1. automated mathematical checks -----------------------------------------
if ~automated_checks(P), error('main:checks', 'Automated checks failed - aborting.'); end

% 2. controller design / analysis ------------------------------------------
D = controller_design(P);
fprintf('Controller: Gc(s) = %.0f (s + %.3f)/s | PM = %.1f deg, GM = %.1f dB, wgc = %.2f rad/s\n', ...
        P.ctrl.Kp, P.ctrl.wi, D.PM_deg, D.GM_dB, D.w_gc);
fprintf('Step: overshoot %.1f %%, rise %.2f s, settling %.2f s; closed-loop poles:\n', ...
        D.step.overshoot_pct, D.step.rise_time, D.step.settling_time);
disp(D.cl_poles.');

% 3. estimator design study (tuning ensemble, disjoint from A-C) ----------
dsfile = fullfile(datadir, 'design_study.mat');
if reuse && exist(dsfile, 'file')
    tmp = load(dsfile); ds = tmp.ds;
else
    ds = design_study(P);  save('-v7', dsfile, 'ds');
end

% 3b. firing-gate tuning (tuning ensemble only; constants are frozen in parameters.m)
ftfile = fullfile(datadir, 'fire_tuning.mat');
if reuse && exist(ftfile, 'file')
    tmp = load(ftfile); FT = tmp.T;
else
    fprintf('Firing-gate tuning on the tuning ensemble ...\n');
    FT = fire_tuning(P, 16);  T = FT;  save('-v7', ftfile, 'T');
end
[~, ia] = max(FT.hit(:) - 1e-3*FT.J(:));  [ib1, ib2] = ind2sub(size(FT.hit), ia);
fprintf('  best gate on tuning set: k_gate = %.2f, k_sig = %.2f (hit rate %.2f, mean t_int %.2f s); frozen values: %.2f / %.2f\n', ...
        FT.kg(ib1), FT.ks(ib2), FT.hit(ib1,ib2), FT.J(ib1,ib2), P.fire.k_gate, P.fire.k_sig);

% 4. scenarios A-C ----------------------------------------------------------
R = struct();  M = struct();
for s = {'A','B','C'}
    fprintf('Simulating scenario %s ...\n', s{1});
    R.(s{1}) = simulation_runner(s{1}, P);
    M.(s{1}) = performance_metrics(R.(s{1}), P);
    m = M.(s{1}).post;
    fprintf('  %s: EKF pos RMS %.2f m (raw %.2f m) | pred. angle RMS %.3f deg | track RMS %.3f deg | look-ahead RMS %.3f deg, max %.3f deg | Vpk %.0f V\n', ...
            s{1}, m.est_pos_rms, m.raw_pos_rms, m.pred_ang_rms, m.trk_rms, m.look_rms, m.look_max, M.(s{1}).V_peak);
end
L = M.C.loss;
fprintf('  C (loss window %.2f-%.2f s): max pred. error %.3f deg, max look-ahead error %.3f deg, EKF 1-sigma(theta_p) up to %.3f deg\n', ...
        L.window(1), L.window(2), L.pred_ang_max, L.look_max, L.sig_thp_max);

% 5. generalisation study (random bounded manoeuvres) ----------------------
fprintf('Robustness study ...\n');
rob = robustness_study(P, ifelse(quick, 4, 10));

% 6. outputs ----------------------------------------------------------------
fprintf('Writing report data and figures ...\n');
Dsave = rmfield(D, {'Gp','Gc','L','T','S','Gcs'});   % tf objects are not portable in MAT files
save('-v7', fullfile(datadir, 'results_ABC.mat'), 'R', 'M', 'Dsave', 'rob', 'ds', 'P');
write_results_csv(M, rob, datadir);
try
    make_figures(R, D, rob, P, figdir);
catch err
    fprintf('FIGURE ERROR: %s\n', err.message);
    for k = 1:numel(err.stack), fprintf('   at %s line %d\n', err.stack(k).name, err.stack(k).line); end
    rethrow(err);
end
fprintf('Done. Figures: %s   Data (CSV): %s\n', figdir, datadir);

end

function s = ifelse(c, a, b), if c, s = a; else, s = b; end, end
