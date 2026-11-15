function D = controller_design(P)
%CONTROLLER_DESIGN  Analytic PI design and closed-loop analysis.
%
%   D = controller_design(P) returns the compensator, loop, closed-loop
%   transfer functions and all stability / performance numbers used in the
%   paper (no plotting; figures are made by make_figures.m).
%
%   Design reasoning (paper Sec. VI):
%   1. Plant:  Gp(s) = Kw / ( s (s^2 + 2 zeta wn s + wn^2) ),  Kw = N Kt/(Lm J_eq).
%      Below wn it behaves as an integrator:  Gp ~ Kv/s,  Kv = Gp-velocity DC gain.
%   2. Compensator  Gc(s) = Kp (s + wi)/s  adds a second integrator (type 2)
%      -> zero steady-state error to a ramp reference (target with constant
%      angular rate) and finite error Ka^-1 to constant angular acceleration.
%   3. For wc << wn:  L(jw) ~ Kp Kv (jw + wi)/(jw)^2  =>  |L(j wc)| = 1
%      gives  Kp = wc/Kv  (valid when wc >> wi), phase margin ~ atan(wc/wi).
%   4. wc is bounded above by the lightly damped electromechanical resonance
%      (zeta ~ 0.09 at wn ~ 56 rad/s): |L(j wn)| ~ (wc/wn)/(2 zeta) must stay
%      < 1 for gain-stabilisation; wc = 5 rad/s gives ~6 dB gain margin.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

Mdl = system_model(P);
Kp = P.ctrl.Kp;  wi = P.ctrl.wi;
D.Gp = Mdl.Gp;
D.Gc = tf(Kp*[1 wi], [1 0]);
D.L  = series(D.Gc, D.Gp);                 % open loop  Gc*Gp
D.T  = feedback(D.L, 1);                   % closed loop  theta/theta_ref
D.S  = feedback(1, D.L);                   % sensitivity
D.Gcs = D.Gc * D.S;                        % control sensitivity  V/theta_ref

% characteristic polynomial  1 + L = 0 :  s^2 den_p(s)/... + Kp N Kt (s+wi)
nL = conv(Kp*[1 wi], Mdl.num);             % numerator of L
dL = conv([1 0], Mdl.den);                 % denominator of L  (s * s * quadratic)
D.char_poly = dL + [zeros(1, numel(dL)-numel(nL)), nL];   % dL + nL
D.cl_poles  = roots(D.char_poly);
D.stable    = all(real(D.cl_poles) < 0);

% Routh-Hurwitz (4th-order) necessary and sufficient conditions
a = D.char_poly(:)';   % [a4 a3 a2 a1 a0]
b1 = (a(2)*a(3) - a(1)*a(4))/a(2);
c1 = (b1*a(4) - a(2)*a(5))/b1;
D.routh_first_col = [a(1), a(2), b1, c1, a(5)];
D.routh_ok = all(D.routh_first_col > 0);

% analytic critical gain from the last Routh-array entry (c1 = 0):
%   Kp_crit = a3 (a2 - a3 wi) / (a4 * N Kt),  a_k = coefficients of the plant denominator
a4 = Mdl.den(1); a3 = Mdl.den(2); a2 = Mdl.den(3);
D.Kp_crit = a3*(a2 - a3*wi)/(a4*Mdl.num);
D.GM_routh_dB = 20*log10(D.Kp_crit/Kp);

[Gm, Pm, Wcg, Wcp] = margin(D.L);
D.GM_dB = 20*log10(Gm);  D.PM_deg = Pm;  D.w_pc = Wcg;  D.w_gc = Wcp;
D.sens_peak = max(abs(squeeze(freqresp(D.S, logspace(-2,3,4000)'))));
D.Ms_dB = 20*log10(D.sens_peak);
w  = logspace(-2,3,4000)';
Tm = abs(squeeze(freqresp(D.T, w)));
D.T_peak_dB = 20*log10(max(Tm));
i3 = find(Tm < 1/sqrt(2), 1, 'first');  D.bw_cl = w(i3);        % -3 dB bandwidth [rad/s]
D.Gcs_hf = abs(squeeze(freqresp(D.Gcs, 1e3)));                  % V/rad at high freq ~ Kp

% error constants of the type-2 loop
D.Kv_err = Inf;                      % velocity constant (ramp error = 0)
D.Ka     = Kp*P.Kv*wi;               % acceleration constant [1/s^2] = wc*wi
D.e_ramp = 0;

% step response of the linear closed loop (unit step reference)
t  = (0:1e-3:8)';
y  = step(D.T, t);
D.step = step_metrics(t, y, 1);
D.t_step = t;  D.y_step = y;

% sinusoidal reference tracking: |T(jw)| at representative frequencies
D.T_at = struct('w',[0.5 1 2 5], 'mag', abs(squeeze(freqresp(D.T,[0.5 1 2 5]'))));
end
