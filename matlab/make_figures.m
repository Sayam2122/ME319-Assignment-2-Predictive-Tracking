function make_figures(R, D, rob, P, outdir)
%MAKE_FIGURES  Generate all figures from simulation results (nothing hand-drawn).
%
%   make_figures(R, D, rob, P, outdir)
%     R      struct with fields A, B, C : results of simulation_runner
%     D      controller_design output
%     rob    robustness_study output (may be [] to skip that figure)
%     P      parameters
%     outdir figure output folder (PNG, 150 dpi, cropped)
%
%   Files written:  fig_controller, fig_estimation, fig_prediction, fig_turret,
%   fig_robustness, fig_overview_A/B/C.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if exist('OCTAVE_VERSION','builtin') ~= 0
    graphics_toolkit('gnuplot');
end
r2d = 180/pi;
set(0, 'defaultaxesfontsize', 9, 'defaultlinelinewidth', 1.2);
col = struct('k',[0 0 0], 'r',[0.85 0.15 0.15], 'b',[0.1 0.3 0.85], ...
             'g',[0.1 0.6 0.2], 'o',[0.95 0.55 0.1], 'gr',[0.6 0.6 0.6]);
bo = [P.t_bo0, P.t_bo0 + P.t_bo];

% =====================================================================
% Fig. controller: Bode of L, root locus, unit-step response of T
% =====================================================================
f = newfig([1000 330]);
w = logspace(-1, 2.5, 800)';
Lw = squeeze(freqresp(D.L, w));  Pw = squeeze(freqresp(D.Gp, w));
subplot(1,3,1);
semilogx(w, 20*log10(abs(Lw)), 'color', col.b); hold on;
semilogx(w, 20*log10(abs(Pw)), '--', 'color', col.gr);
semilogx([w(1) w(end)], [0 0], 'k:');
semilogx(D.w_gc, 0, 'o', 'color', col.r);
semilogx(D.w_pc, -D.GM_dB, 's', 'color', col.r);
axis([0.1 300 -80 60]); grid on; xlabel('frequency [rad/s]'); ylabel('magnitude [dB]');
title(sprintf('(a) open loop L: wgc=%.1f rad/s, GM=%.1f dB', D.w_gc, D.GM_dB));
legend('L = Gc Gp', 'Gp alone', 'location', 'northeast'); legend('boxoff');

subplot(1,3,2);
ph = unwrap(angle(Lw))*r2d;  ph = ph - 360*round(ph(1)/360);
php = unwrap(angle(Pw))*r2d;
semilogx(w, ph, 'color', col.b); hold on;
semilogx([w(1) w(end)], [-180 -180], 'k:');
semilogx(D.w_gc, -180 + D.PM_deg, 'o', 'color', col.r);
semilogx(D.w_pc, -180, 's', 'color', col.r);
axis([0.1 300 -270 -80]); grid on; xlabel('frequency [rad/s]'); ylabel('phase [deg]');
title(sprintf('(b) phase: PM = %.1f deg', D.PM_deg));

subplot(1,3,3);
plot(D.t_step, D.y_step, 'color', col.b); hold on;
plot([0 8], [1.02 1.02], ':', 'color', col.gr); plot([0 8], [0.98 0.98], ':', 'color', col.gr);
axis([0 6 0 1.4]); grid on; xlabel('time [s]'); ylabel('theta / theta\_ref');
title(sprintf('(c) step: OS=%.1f%%, Ts=%.2f s', D.step.overshoot_pct, D.step.settling_time));
saveboth(f, outdir, 'fig_controller');

% root locus in its own figure (kept separate so axes are readable)
f = newfig([480 360]);
Gu = series(tf([1 P.ctrl.wi],[1 0]), D.Gp);          % locus over Kp
[rl, kk] = rlocus(Gu); %#ok<ASGLU>
plot(real(rl)', imag(rl)', '.', 'color', col.gr, 'markersize', 3); hold on;
plot(real(D.cl_poles), imag(D.cl_poles), 'x', 'color', col.r, 'markersize', 9, 'linewidth', 2);
op = roots(D.Gp.den{1});  plot(real(op), imag(op), 'o', 'color', col.b);
plot(-P.ctrl.wi, 0, 'o', 'color', col.g);
axis([-20 5 -75 75]); grid on; xlabel('real axis [1/s]'); ylabel('imaginary axis [rad/s]');
title(sprintf('root locus vs Kp; x = closed-loop poles at Kp = %.0f V/rad', P.ctrl.Kp));
saveboth(f, outdir, 'fig_rootlocus');

% =====================================================================
% Fig. estimation (scenario B)
% =====================================================================
b = R.B;  t = b.t;  vs = b.valid > 0;
xr = b.rm.*cos(b.thm);  yr = b.rm.*sin(b.thm);
f = newfig([1000 560]);
subplot(2,3,1);
plot(b.tx, b.ty, 'k'); hold on;
plot(xr(vs), yr(vs), '.', 'color', col.gr, 'markersize', 3);
plot(b.xh, b.yh, '--', 'color', col.r);
xlabel('x [m]'); ylabel('y [m]'); grid on; title('(a) trajectory, radar samples, EKF');
legend('true', 'radar (converted)', 'EKF', 'location', 'southeast'); legend('boxoff');
subplot(2,3,2);
plot(t(vs), (b.rm(vs) - hypot(b.tx(vs), b.ty(vs))), '.', 'color', col.b, 'markersize', 4);
hold on; plot([0 10], P.sigma_r*[1 1], 'k:', [0 10], -P.sigma_r*[1 1], 'k:');
xlabel('t [s]'); ylabel('range error [m]'); grid on; title('(b) range measurement noise');
subplot(2,3,3);
plot(t(vs), (b.thm(vs) - b.tth(vs))*r2d, '.', 'color', col.b, 'markersize', 4);
hold on; plot([0 10], P.sigma_theta*r2d*[1 1], 'k:', [0 10], -P.sigma_theta*r2d*[1 1], 'k:');
xlabel('t [s]'); ylabel('elevation error [deg]'); grid on; title('(c) elevation measurement noise');
subplot(2,3,4);
er = hypot(xr - b.tx, yr - b.ty);  ep = hypot(b.xh - b.tx, b.yh - b.ty);
plot(t(vs), er(vs), '.', 'color', col.gr, 'markersize', 4); hold on;
plot(t, ep, 'color', col.r);
xlabel('t [s]'); ylabel('position error [m]'); grid on; title('(d) position error');
legend('raw radar', 'EKF', 'location', 'northeast'); legend('boxoff');
subplot(2,3,5);
plot(t, b.tvx, 'k'); hold on; plot(t, b.tvy, 'k');
plot(t, b.vxh, '--', 'color', col.r); plot(t, b.vyh, '--', 'color', col.b);
xlabel('t [s]'); ylabel('velocity [m/s]'); grid on; title('(e) velocity: true (black) vs EKF');
subplot(2,3,6);
nis = b.nis; k = find(~isnan(nis)); ma = filter(ones(1,20)/20, 1, nis(k));
plot(t(k), nis(k), '.', 'color', col.gr, 'markersize', 4); hold on;
plot(t(k), ma, 'color', col.b); plot([0 10], [2 2], 'k:');
xlabel('t [s]'); ylabel('NIS'); grid on; title('(f) innovation NIS (expected mean 2)');
axis([0 10 0 12]);
saveboth(f, outdir, 'fig_estimation');

% =====================================================================
% Fig. prediction
% =====================================================================
f = newfig([1000 560]);
c = R.C;
subplot(2,2,1);
plot(c.t, c.fth*r2d, 'k'); hold on; plot(c.t, c.thp*r2d, '--', 'color', col.r);
plot(c.t, c.th*r2d, 'color', col.b);
xlabel('t [s]'); ylabel('elevation [deg]'); grid on;
title('(a) scenario C: true future angle, prediction, platform');
legend('true theta(t+Tp)', 'predicted theta_p(t)', 'platform theta(t)', 'location', 'southeast'); legend('boxoff');
subplot(2,2,2);
m = c.t >= 0.5 & c.t <= 3.0;
plot(c.t(m), c.fth(m)*r2d, 'k'); hold on; plot(c.t(m), c.thp(m)*r2d, '--', 'color', col.r);
plot(c.t(m), c.th(m)*r2d, 'color', col.b);
yl = ylim; plot(bo(1)*[1 1], yl, ':', 'color', col.g); plot(bo(2)*[1 1], yl, ':', 'color', col.g);
xlabel('t [s]'); ylabel('elevation [deg]'); grid on; title('(b) zoom on the measurement-loss interval');
subplot(2,2,3);
for q = 1:3
    s = {'A','B','C'}; rr = R.(s{q}); cc = {col.g, col.b, col.r};
    plot(rr.t, wrap(rr.thp - rr.fth)*r2d, 'color', cc{q}); hold on;
end
yl = ylim; plot(bo(1)*[1 1], yl, ':', 'color', col.k); plot(bo(2)*[1 1], yl, ':', 'color', col.k);
xlabel('t [s]'); ylabel('prediction error [deg]'); grid on; title('(c) prediction error theta_p(t) - theta(t+Tp)');
legend('A nominal', 'B manoeuvre', 'C manoeuvre + loss', 'location', 'northeast'); legend('boxoff');
subplot(2,2,4);
k = ~isnan(c.sig_thp);
plot(c.t, abs(wrap(c.thp - c.fth))*r2d, 'color', col.gr); hold on;
plot(c.t(k), c.sig_thp(k)*r2d, 'color', col.r);
yl = ylim; plot(bo(1)*[1 1], yl, ':', 'color', col.k); plot(bo(2)*[1 1], yl, ':', 'color', col.k);
xlabel('t [s]'); ylabel('[deg]'); grid on; title('(d) scenario C: |error| vs EKF 1-sigma');
legend('|error|', 'EKF 1-sigma of theta_p', 'location', 'northeast'); legend('boxoff');
saveboth(f, outdir, 'fig_prediction');

% =====================================================================
% Fig. turret response
% =====================================================================
f = newfig([1000 560]);
subplot(2,2,1);
m = c.t >= 0 & c.t <= 10;
plot(c.t(m), c.th_ref(m)*r2d, 'k'); hold on; plot(c.t(m), c.th(m)*r2d, '--', 'color', col.b);
xlabel('t [s]'); ylabel('elevation [deg]'); grid on; title('(a) scenario C: reference vs platform angle');
legend('reference theta_{ref}', 'platform theta', 'location', 'southeast'); legend('boxoff');
subplot(2,2,2);
for q = 1:3
    s = {'A','B','C'}; rr = R.(s{q}); cc = {col.g, col.b, col.r};
    plot(rr.t, wrap(rr.th_ref - rr.th)*r2d, 'color', cc{q}); hold on;
end
xlabel('t [s]'); ylabel('tracking error [deg]'); grid on; title('(b) tracking error theta_{ref} - theta');
legend('A', 'B', 'C', 'location', 'northeast'); legend('boxoff');
subplot(2,2,3);
for q = 1:3
    s = {'A','B','C'}; rr = R.(s{q}); cc = {col.g, col.b, col.r};
    plot(rr.t, wrap(rr.fth - rr.th)*r2d, 'color', cc{q}); hold on;
end
xlabel('t [s]'); ylabel('look-ahead error [deg]'); grid on; title('(c) theta(t+Tp) - platform angle');
subplot(2,2,4);
plot(c.t, c.V, 'color', col.b); hold on; plot(c.t, 10*c.I, 'color', col.o);
xlabel('t [s]'); ylabel('V_m [V]  and  10 x i [A]'); grid on; title('(d) scenario C: motor voltage and current');
legend('V_m', '10 i', 'location', 'northeast'); legend('boxoff');
saveboth(f, outdir, 'fig_turret');

% =====================================================================
% Fig. robustness
% =====================================================================
if ~isempty(rob)
    f = newfig([1000 340]);
    subplot(1,2,1);
    ix = strcmp(rob.fields, 'look_rms');  iy = strcmp(rob.fields, 'look_max');
    n = rob.ncase;
    bar(1:n, [rob.tab(:,ix) rob.tab(:,iy)], 'grouped'); hold on;
    xl = find(rob.with_loss); plot(xl, 0.01*ones(size(xl)), 'v', 'color', col.r, 'markerfacecolor', col.r);
    xlabel('random case index (red marker: sensor loss)'); ylabel('look-ahead error [deg]'); grid on;
    title('(a) random manoeuvres, RMS and max (t >= 1 s)');
    legend('RMS', 'max', 'location', 'northeast'); legend('boxoff');
    subplot(1,2,2);
    plot(rob.onset, rob.sweep_loss, 'o-', 'color', col.r); hold on;
    plot(rob.onset, rob.sweep_ref, 's--', 'color', col.b);
    xlabel('loss onset [s]'); ylabel('max prediction error [deg]'); grid on;
    title('(b) 0.5 s loss: onset sweep, manoeuvre B');
    legend('with loss', 'no loss (same window)', 'location', 'northeast'); legend('boxoff');
    saveboth(f, outdir, 'fig_robustness');
end

% =====================================================================
% Per-scenario overview figures
% =====================================================================
for s = {'A','B','C'}
    plot_overview(R.(s{1}), s{1}, P, outdir, col);
end
% ---- engagement figure: beam miss distance vs time, pulses marked ----
f = newfig([900 520]);
k = 0;
for s = {'A','B','C'}
    k = k + 1;  r = R.(s{1});
    miss = hypot(r.fx, r.fy) .* abs(sin(wrap(r.fth - r.th)));
    subplot(3,1,k); semilogy(r.t, max(miss, 1e-2), 'color', col.b); hold on;
    plot([0 r.t(end)], P.R_eng*[1 1], '--', 'color', col.r);
    pl = r.pulses;
    if ~isempty(pl)
        ok = pl(:,5) > 0;
        plot(pl(ok,1), max(pl(ok,4),1e-2), 'o', 'color', col.g, 'markerfacecolor', col.g);
        plot(pl(~ok,1), max(pl(~ok,4),1e-2), 'x', 'color', col.r, 'markersize', 9);
    end
    grid on; ylabel('miss [m]'); title(sprintf('Scenario %s: beam-to-target(t+T_p) miss distance, pulses (o hit, x miss), R_{eng}=2 m', s{1}));
    if k == 3, xlabel('t [s]'); end
end
saveboth(f, outdir, 'fig_engagement');
close all;
end

% =====================================================================
function plot_overview(r, name, P, outdir, col)
r2d = 180/pi;  t = r.t;  vs = r.valid > 0;
xr = r.rm.*cos(r.thm);  yr = r.rm.*sin(r.thm);
f = newfig([1100 800]);
subplot(3,3,1); plot(r.tx, r.ty, 'k'); hold on; plot(xr(vs), yr(vs), '.', 'color', col.gr, 'markersize', 3);
plot(r.xh, r.yh, '--', 'color', col.r); grid on; xlabel('x [m]'); ylabel('y [m]'); title('trajectory');
subplot(3,3,2); plot(t, r.xh - r.tx, 'color', col.r); hold on; plot(t, r.yh - r.ty, 'color', col.b);
grid on; xlabel('t [s]'); ylabel('error [m]'); title('EKF position error (x red, y blue)');
subplot(3,3,3); plot(t, r.tvx, 'k', t, r.tvy, 'k'); hold on; plot(t, r.vxh, '--', 'color', col.r); plot(t, r.vyh, '--', 'color', col.b);
grid on; xlabel('t [s]'); ylabel('[m/s]'); title('velocity: true (black) vs estimate');
subplot(3,3,4); plot(t(vs), r.rm(vs) - hypot(r.tx(vs), r.ty(vs)), '.', 'color', col.b, 'markersize', 4);
grid on; xlabel('t [s]'); ylabel('[m]'); title('range measurement noise');
subplot(3,3,5); plot(t(vs), (r.thm(vs) - r.tth(vs))*r2d, '.', 'color', col.b, 'markersize', 4);
grid on; xlabel('t [s]'); ylabel('[deg]'); title('elevation measurement noise');
subplot(3,3,6); plot(t, r.fth*r2d, 'k'); hold on; plot(t, r.thp*r2d, '--', 'color', col.r); plot(t, r.th*r2d, 'color', col.b);
grid on; xlabel('t [s]'); ylabel('[deg]'); title('theta(t+Tp) (black), predicted (red), platform (blue)');
subplot(3,3,7); plot(t, wrap(r.th_ref - r.th)*r2d, 'color', col.b); hold on; plot(t, wrap(r.fth - r.th)*r2d, 'color', col.r);
grid on; xlabel('t [s]'); ylabel('[deg]'); title('tracking (blue) and look-ahead (red) error');
subplot(3,3,8); plot(t, r.V, 'color', col.b); grid on; xlabel('t [s]'); ylabel('V_m [V]'); title('control input');
subplot(3,3,9); k = ~isnan(r.sig_thp); plot(t, abs(wrap(r.thp - r.fth))*r2d, 'color', col.gr); hold on;
plot(t(k), r.sig_thp(k)*r2d, 'color', col.r); grid on; xlabel('t [s]'); ylabel('[deg]'); title('|prediction error| vs EKF 1-sigma');
saveboth(f, outdir, ['fig_overview_' name]);
end

function f = newfig(sz)
f = figure('visible', 'off');
set(f, 'position', [50 50 sz(1) sz(2)], 'paperunits', 'points', ...
       'paperposition', [0 0 sz(1) sz(2)], 'papersize', [sz(1) sz(2)]);
end

function saveboth(f, outdir, name)
%SAVEBOTH  Export PNG (150 dpi) and crop the surrounding white margin.
file = fullfile(outdir, [name '.png']);
print(f, '-dpng', '-r150', file);
try
    A = imread(file);
    g = min(A, [], 3);
    [rr, cc] = find(g < 245);
    pad = 10;
    r1 = max(1, min(rr) - pad);  r2 = min(size(A,1), max(rr) + pad);
    c1 = max(1, min(cc) - pad);  c2 = min(size(A,2), max(cc) + pad);
    imwrite(A(r1:r2, c1:c2, :), file);
catch
    % cropping is cosmetic; keep the uncropped image if imread/imwrite is unavailable
end
close(f);
end

function a = wrap(a), a = atan2(sin(a), cos(a)); end
