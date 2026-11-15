function [F, Q] = ca_singer_model(dt, alpha, q)
%CA_SINGER_MODEL  Discrete transition matrix and process noise of the
%   Singer-type acceleration model used for the target estimator/predictor.
%
%   Per axis, with state s = [p; v; a]:
%       ds/dt = Ac s + G w,   Ac = [0 1 0; 0 0 1; 0 0 -alpha],  G = [0;0;1]
%       E[w(t) w(tau)] = q delta(t - tau)
%   (alpha = 0 gives the nearly-constant-acceleration / white-jerk model.)
%
%   The 6-state vector is ordered  X = [x; vx; ax; y; vy; ay]  so that
%   F = blkdiag(F1, F1) and Q = blkdiag(Q1, Q1).
%
%   [F, Q] = ca_singer_model(dt, alpha, q)
%     F  : 6x6 closed-form transition matrix  exp(Ac*dt) per axis
%     Q  : 6x6 process-noise covariance (Van Loan's method); only computed
%          when nargout >= 2 (it needs a matrix exponential).
%
%   Closed form of exp(Ac*dt):
%       [1  dt  (alpha*dt - 1 + e^{-alpha dt})/alpha^2;
%        0  1   (1 - e^{-alpha dt})/alpha;
%        0  0   e^{-alpha dt}]
%   with the limit [1 dt dt^2/2; 0 1 dt; 0 0 1] as alpha -> 0.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if dt < 0, error('ca_singer_model:negativeDt', 'dt must be >= 0.'); end

ead = exp(-alpha*dt);
if alpha*dt < 1e-4                      % series for tiny alpha*dt (avoids 0/0)
    f13 = dt^2/2 - alpha*dt^3/6;
    f23 = dt - alpha*dt^2/2;
else
    f13 = (alpha*dt - 1 + ead)/alpha^2;
    f23 = (1 - ead)/alpha;
end
F1 = [1 dt f13; 0 1 f23; 0 0 ead];
F  = blkdiag(F1, F1);

if nargout >= 2
    Ac = [0 1 0; 0 0 1; 0 0 -alpha];
    G  = [0; 0; 1];
    n  = 3;
    M  = [-Ac, G*q*G'; zeros(n), Ac'] * dt;       % Van Loan block matrix
    E  = expm(M);
    F1v = E(n+1:2*n, n+1:2*n)';                   % == exp(Ac dt)
    Q1  = F1v * E(1:n, n+1:2*n);                  % integral F G q G' F' dtau
    Q1  = (Q1 + Q1')/2;                           % enforce symmetry
    Q   = blkdiag(Q1, Q1);
end
end
