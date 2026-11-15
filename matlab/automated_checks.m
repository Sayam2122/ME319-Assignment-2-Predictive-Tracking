function ok = automated_checks(P)
%AUTOMATED_CHECKS  Self-tests of the key mathematical relationships.
%
%   ok = automated_checks(P) runs every check, prints PASS/FAIL, and returns
%   true only if all pass.  main.m aborts if any check fails.
%
%   Checks
%    1  reflected inertia and derived constants
%    2  dimensional homogeneity of the electrical and mechanical equations
%       (unit-exponent bookkeeping in kg, m, s, A; rad is dimensionless)
%    3  transfer function == state-space realisation (poles, DC velocity gain)
%    4  closed-form wn, zeta match the numerical open-loop poles
%    5  energy balance of the simulated plant (power in = dissipation + storage)
%    6  Singer model: semigroup property, white-jerk limit, Q symmetric PSD
%    7  EKF Jacobian vs. finite differences
%    8  generator: scenarios A-C respect range/speed/acceleration limits;
%       noise statistics match sigma_r, sigma_theta; loss mask is correct
%    9  noise-free prediction of a constant-velocity target is exact
%   10  closed loop: Routh-Hurwitz, pole locations, margins, type-2 property
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

results = cell(0,3);  % rows: {name, pass, detail}
mk = @(name, pass, detail) {name, logical(pass), detail};

% ---- 1 -----------------------------------------------------------------
results(end+1,:) = mk('1  J_eq = J_turret + N^2 J_m = 3200 kg m^2', abs(P.J_eq - 3200) < 1e-9, ...
              sprintf('J_eq = %.3f', P.J_eq));

% ---- 2 dimensional analysis ---------------------------------------------
U = @(kg,m,s,A) [kg m s A];                 % exponent vector
V_  = U(1,2,-3,-1);  Nm = U(1,2,-2,0);  A_ = U(0,0,0,1);  s_ = U(0,0,1,0);
Ohm = V_ - A_;  H_ = V_ + s_ - A_;
Kt_u = Nm - A_;                             % N m / A
Ke_u = V_ + s_;                             % V s / rad
b_u  = Nm + s_;                             % N m s / rad
J_u  = Nm + 2*s_;                           % kg m^2  == N m s^2
elec = [H_ + A_ - s_;  Ohm + A_;  Ke_u - s_];         % Lm di/dt, Rm i, Ke w_m  -> all V
mech = [J_u - 2*s_;  Kt_u + A_;  b_u - s_];           % J th'', N Kt i, b th'   -> all N m
ok2 = all(all(elec == repmat(V_,3,1))) && all(all(mech == repmat(Nm,3,1))) ...
      && isequal(J_u, U(1,2,0,0)) && isequal(Kt_u, Ke_u);
results(end+1,:) = mk('2  dimensional homogeneity of armature & mechanical equations; [Kt]=[Ke]', ok2, ...
              'unit exponents (kg,m,s,A) all consistent');

% ---- 3 tf vs ss ----------------------------------------------------------
M = system_model(P);
pt = sort(real(roots(M.den))); 
ps = sort(real(eig(M.A)));
Gss = tf(M.sys);
w0 = 1e-3; g_tf = abs(squeeze(freqresp(M.Gp, 1)));  g_ss = abs(squeeze(freqresp(Gss, 1)));
ok3 = norm(sort(abs(eig(M.A))) - sort(abs(roots(M.den)))) < 1e-6*max(abs(roots(M.den))+1) ...
      && abs(g_tf - g_ss)/g_tf < 1e-6 && abs(dcgain(M.Gw) - P.Kv)/P.Kv < 1e-9;
results(end+1,:) = mk('3  tf == state space (poles, |G(j1)|); velocity DC gain = Kv', ok3, ...
              sprintf('Kv = %.6e rad/s/V', P.Kv));

% ---- 4 closed-form wn, zeta ----------------------------------------------
pc = roots(M.den(1:3));                       % quadratic factor
wn_n = abs(pc(1));  zeta_n = -real(pc(1))/abs(pc(1));
ok4 = abs(wn_n - P.wn) < 1e-6*P.wn && abs(zeta_n - P.zeta) < 1e-6;
results(end+1,:) = mk('4  closed-form wn, zeta == numerical poles', ok4, ...
              sprintf('wn = %.3f rad/s, zeta = %.4f', P.wn, P.zeta));

% ---- 5 energy balance (RK4 on turret_plant, V = 100 V, 3 s) ---------------
dt = 1e-4; n = 30000; x = [0;0;0]; V = 100;
Ein = 0; Ediss = 0;
f = @(x) cell2mat_p(V, x);
for k = 1:n
    Ein   = Ein   + V*x(3)*dt;
    Ediss = Ediss + (P.Rm*x(3)^2 + P.b*x(2)^2)*dt;
    k1=f(x); k2=f(x+0.5*dt*k1); k3=f(x+0.5*dt*k2); k4=f(x+dt*k3);
    x = x + dt/6*(k1+2*k2+2*k3+k4);
end
Estore = 0.5*P.Lm*x(3)^2 + 0.5*P.J_eq*x(2)^2;
resid  = abs(Ein - Ediss - Estore)/Ein;
results(end+1,:) = mk('5  energy balance of turret_plant (rel. residual < 1e-3)', resid < 1e-3, ...
              sprintf('rel. residual = %.2e', resid));

% ---- 6 Singer model ------------------------------------------------------
[F1,Q1] = ca_singer_model(0.05, P.est.alpha, P.est.q);
F2 = ca_singer_model(0.07, P.est.alpha, P.est.q);
F3 = ca_singer_model(0.12, P.est.alpha, P.est.q);
semi = norm(F2*F1 - F3, 'fro') < 1e-12;
[Fz,Qz] = ca_singer_model(0.05, 0, 3);
T_ = 0.05; Qw = 3*[T_^5/20 T_^4/8 T_^3/6; T_^4/8 T_^3/3 T_^2/2; T_^3/6 T_^2/2 T_];
lim0 = norm(Qz(1:3,1:3) - Qw, 'fro') < 1e-9*norm(Qw,'fro') && ...
       norm(Fz(1:3,1:3) - [1 T_ T_^2/2; 0 1 T_; 0 0 1],'fro') < 1e-14;
psd  = norm(Q1 - Q1','fro') < 1e-14 && min(eig(Q1)) > -1e-12*max(eig(Q1));
results(end+1,:) = mk('6  Singer model: semigroup, white-jerk limit, Q symmetric PSD', semi && lim0 && psd, ...
              sprintf('semigroup=%d limit=%d psd=%d', semi, lim0, psd));

% ---- 7 EKF Jacobian ------------------------------------------------------
X0 = [1500; 0;0; 600; 0;0];
h  = @(X) [hypot(X(1),X(4)); atan2(X(4),X(1))];
Hn = zeros(2,6); eps_ = 1e-4;
for c = [1 4]
    dX = zeros(6,1); dX(c) = eps_;
    Hn(:,c) = (h(X0+dX) - h(X0-dX))/(2*eps_);
end
x=X0(1); y=X0(4); r2=x^2+y^2; r=sqrt(r2);
Ha = zeros(2,6); Ha(1,1)=x/r; Ha(1,4)=y/r; Ha(2,1)=-y/r2; Ha(2,4)=x/r2;
results(end+1,:) = mk('7  EKF measurement Jacobian == finite difference', norm(Ha-Hn,'fro') < 1e-8, ...
              sprintf('max err = %.2e', max(abs(Ha(:)-Hn(:)))));

% ---- 8 generator ---------------------------------------------------------
okg = true; det = '';
for sc = {'A','B','C'}
    try
        [~,~,~,~] = target_motion_generator(0, sc{1});
    catch err
        okg = false; det = [det sprintf('%s: %s; ', sc{1}, err.message)];
    end
end
% noise statistics from scenario A over its whole record (K = 221 samples)
Kn = round(P.t_end/P.Ts_radar) + 1; dr = zeros(Kn,1); dth = zeros(Kn,1);
for k = 0:Kn-1
    [rm,thm,~,tr] = target_motion_generator(k*P.Ts_radar, 'A');
    dr(k+1) = rm - tr.r;  dth(k+1) = thm - tr.theta;
end
okn = abs(std(dr)/P.sigma_r - 1) < 0.15 && abs(std(dth)/P.sigma_theta - 1) < 0.15 ...
      && abs(mean(dr)) < 0.3*P.sigma_r && abs(mean(dth)) < 0.3*P.sigma_theta;
% loss mask in scenario C
nv = 0; first_bad = NaN;
for k = 0:round(P.t_end/P.Ts_radar)
    [~,~,v] = target_motion_generator(k*P.Ts_radar, 'C');
    if ~v, nv = nv + 1; if isnan(first_bad), first_bad = k*P.Ts_radar; end, end
end
okm = (nv == round(P.t_bo/P.Ts_radar)) && abs(first_bad - P.t_bo0) < 1e-9;
results(end+1,:) = mk('8  generator: limits enforced, noise std ~ sigma, loss mask correct', okg && okn && okm, ...
              sprintf('%s std(r)=%.2f m, std(th)=%.4f deg, lost samples=%d', det, std(dr), ...
              std(dth)*180/pi, nv));

% ---- 9 noise-free constant-velocity prediction ---------------------------
S = [];  c = P.scen.A;
for k = 0:40
    t = k*P.Ts_radar;  [~,~,~,tr] = target_motion_generator(t,'A');
    S = state_estimator(S, t, tr.r, tr.theta, P);      % noise-free "measurements"
end
tN = 40*P.Ts_radar; pr = target_predictor(S, tN, P);
[~,~,~,trf] = target_motion_generator(tN + P.Tp, 'A');
e9 = abs(pr.theta_p - trf.theta)*180/pi;
results(end+1,:) = mk('9  noise-free constant-velocity prediction error < 0.01 deg after 2 s', e9 < 0.01, ...
              sprintf('error = %.2e deg', e9));

% ---- 10 closed loop ------------------------------------------------------
D = controller_design(P);
ok10 = abs(D.GM_routh_dB - D.GM_dB) < 0.05 && D.stable && D.routh_ok && D.PM_deg > 60 && D.GM_dB > 6 && isinf(D.Kv_err) ...
       && abs(abs(squeeze(freqresp(D.L, D.w_gc))) - 1) < 1e-3;
results(end+1,:) = mk('10 closed loop: Routh-Hurwitz, poles in LHP, PM>60 deg, GM>6 dB (== analytic Routh gain limit), |L(jwgc)|=1', ok10, ...
              sprintf('PM=%.1f deg, GM=%.1f dB, max Re(pole)=%.3f', D.PM_deg, D.GM_dB, max(real(D.cl_poles))));

% ---- report --------------------------------------------------------------
ok = true;
fprintf('\n=== Automated checks ===\n');
for k = 1:size(results,1)
    ok = ok && results{k,2};
    fprintf('[%s] %s  (%s)\n', ifelse(results{k,2},'PASS','FAIL'), results{k,1}, results{k,3});
end
fprintf('=== %s ===\n\n', ifelse(ok, 'ALL CHECKS PASSED', 'SOME CHECKS FAILED'));
end

function xd = cell2mat_p(V, x)
[a,b,c] = turret_plant(0, V, x(1), x(2), x(3));
xd = [a;b;c];
end
function s = ifelse(c, a, b), if c, s = a; else, s = b; end, end
