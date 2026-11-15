function pred = target_predictor(S, t, P, with_cov)
%TARGET_PREDICTOR  Open-loop propagation of the EKF estimate and conversion
%   of the predicted target position into a look-ahead elevation reference.
%
%   pred = target_predictor(S, t, P)
%   pred = target_predictor(S, t, P, true)   % also returns angle uncertainty
%
%   The filter state S.X is valid at time S.t (time of the latest radar
%   update).  With the Singer model, the state at any later time is
%       X(t)      = Phi(t - S.t)        * S.X       (current estimate)
%       X(t + Tp) = Phi(t - S.t + Tp)   * S.X       (prediction, horizon Tp)
%   which needs no new data and therefore keeps working during a
%   measurement-loss interval.
%
%   Output fields
%     X_now, X_pred   6x1 states at t and t+Tp  [x vx ax y vy ay]
%     x_hat, y_hat    estimated position at t                      [m]
%     vx_hat, vy_hat  estimated velocity at t                      [m/s]
%     r_hat, theta_hat  estimated range / elevation at t           [m, rad]
%     x_p, y_p        predicted position at t+Tp                   [m]
%     r_p, theta_p    predicted range / elevation at t+Tp          [m, rad]
%     theta_ref       theta_p limited to the mechanical range      [rad]
%     sig_theta_p     (with_cov) 1-sigma uncertainty of theta_p    [rad]
%     sig_pos_hat     (with_cov) 1-sigma position uncertainty at t [m]
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if nargin < 4, with_cov = false; end
d_now  = t - S.t;
d_pred = d_now + P.Tp;
if d_now < -1e-12, error('target_predictor:timeReversal', 't precedes last update.'); end

F_now  = ca_singer_model(max(d_now,0), P.est.alpha, P.est.q);
F_pred = ca_singer_model(d_pred,       P.est.alpha, P.est.q);
Xn = F_now  * S.X;
Xp = F_pred * S.X;

pred.X_now = Xn;   pred.X_pred = Xp;
pred.x_hat = Xn(1); pred.y_hat = Xn(4);
pred.vx_hat = Xn(2); pred.vy_hat = Xn(5);
pred.r_hat = hypot(Xn(1), Xn(4));  pred.theta_hat = atan2(Xn(4), Xn(1));
pred.x_p = Xp(1);  pred.y_p = Xp(4);
pred.r_p = hypot(Xp(1), Xp(4));    pred.theta_p = atan2(Xp(4), Xp(1));
pred.theta_ref = min(max(pred.theta_p, P.theta_min), P.theta_max);

if with_cov
    [~, Qp] = ca_singer_model(d_pred, P.est.alpha, P.est.q);
    Pp = F_pred * S.Pc * F_pred' + Qp;                    % covariance at t+Tp
    xp = Xp(1); yp = Xp(4); r2 = xp^2 + yp^2;
    gth = zeros(1,6); gth(1) = -yp/r2; gth(4) = xp/r2;    % d theta / d X
    pred.sig_theta_p = sqrt(gth*Pp*gth');
    [~, Qn] = ca_singer_model(max(d_now,1e-9), P.est.alpha, P.est.q);
    Pn = F_now * S.Pc * F_now' + Qn;
    pred.sig_pos_hat = sqrt(Pn(1,1) + Pn(4,4));
end
end
