function S = state_estimator(S, t, r_meas, theta_meas, P)
%STATE_ESTIMATOR  Extended Kalman filter (EKF) for noisy range/elevation data.
%
%   S = state_estimator(S, t, r_meas, theta_meas, P)
%
%   Called ONLY when a new radar measurement is available (meas_valid = 1).
%   Between calls (including a measurement-loss interval) the estimate is
%   propagated open-loop by target_predictor.m.  When the next measurement
%   arrives after a gap, the time update below uses the full elapsed time
%   dt = t - S.t, so the covariance honestly reflects the lost updates.
%
%   State (Cartesian, elevation plane):  X = [x; vx; ax; y; vy; ay]
%   Dynamics: independent Singer acceleration model on each axis
%             (ca_singer_model.m, parameters P.est.q, P.est.alpha).
%   Measurement:  z = [r; theta] = [sqrt(x^2+y^2); atan2(y,x)] + v,
%                 v ~ N(0, R),  R = diag(sigma_r^2, sigma_theta^2).
%   EKF Jacobian (row order r, theta; columns x, vx, ax, y, vy, ay):
%        H = [ x/r, 0, 0, y/r, 0, 0;  -y/r^2, 0, 0, x/r^2, 0, 0 ].
%
%   S fields:  init (logical), t, X (6x1), Pc (6x6), nis (innovation NIS),
%              innov (2x1 innovation of the last update).
%   Pass S = [] to start a new filter; the first measurement initialises it:
%      position  <- polar-to-Cartesian conversion of the first sample,
%      velocity/acceleration <- 0 with prior std from a uniform prior over
%      the admissible set [-vmax,vmax], [-amax,amax] (std = bound/sqrt(3)).
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

persistent c_dt c_F c_Q c_par

R = diag([P.sigma_r^2, P.sigma_theta^2]);

% ------------------------------------------------------------------ init --
if isempty(S) || ~isfield(S,'init') || ~S.init
    x0 = r_meas*cos(theta_meas);  y0 = r_meas*sin(theta_meas);
    J  = [cos(theta_meas), -r_meas*sin(theta_meas); ...
          sin(theta_meas),  r_meas*cos(theta_meas)];
    Pp = J*R*J';                                   % Cartesian position covariance
    Pc = zeros(6);
    Pc([1 4],[1 4]) = Pp;
    Pc(2,2) = P.est.sig_v0^2;  Pc(5,5) = P.est.sig_v0^2;
    Pc(3,3) = P.est.sig_a0^2;  Pc(6,6) = P.est.sig_a0^2;
    S = struct('init', true, 't', t, 'X', [x0;0;0;y0;0;0], 'Pc', Pc, ...
               'nis', NaN, 'innov', [0;0]);
    return
end

% ----------------------------------------------------------- time update --
dt = t - S.t;
if dt < 0, error('state_estimator:timeReversal', 'Measurement time went backwards.'); end
if dt > 0
    par = [P.est.alpha, P.est.q];
    if isempty(c_dt) || abs(dt - c_dt) > 1e-9 || ~isequal(par, c_par)
        [c_F, c_Q] = ca_singer_model(dt, P.est.alpha, P.est.q);
        c_dt = dt;  c_par = par;
    end
    X  = c_F * S.X;
    Pc = c_F * S.Pc * c_F' + c_Q;
else
    X = S.X;  Pc = S.Pc;
end

% -------------------------------------------------------- measurement update
x = X(1);  y = X(4);  rho2 = x^2 + y^2;  rho = sqrt(rho2);
h = [rho; atan2(y, x)];
H = zeros(2,6);
H(1,1) =  x/rho;     H(1,4) = y/rho;
H(2,1) = -y/rho2;    H(2,4) = x/rho2;

nu = [r_meas; theta_meas] - h;
nu(2) = atan2(sin(nu(2)), cos(nu(2)));            % wrap angle innovation
Sm = H*Pc*H' + R;
K  = (Pc*H') / Sm;
X  = X + K*nu;
IKH = eye(6) - K*H;
Pc  = IKH*Pc*IKH' + K*R*K';                       % Joseph form (keeps Pc PSD)
Pc  = (Pc + Pc')/2;

S.t = t;  S.X = X;  S.Pc = Pc;
S.nis = nu' * (Sm \ nu);
S.innov = nu;
end
