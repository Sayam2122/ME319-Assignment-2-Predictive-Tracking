function pulses = replay_gate(res, P)
%REPLAY_GATE  Offline replay of the controller's gate on logged signals
%   (evaluated at the 10 ms instants where sig_thp was logged).
pulses = zeros(0,6);  tl = -Inf;
for j = 1:numel(res.t)
    t = res.t(j);
    if isnan(res.sig_thp(j)) || t < P.fire.t_arm || (t - tl) < P.t_recharge - 1e-9, continue, end
    if res.thp(j) < P.theta_min || res.thp(j) > P.theta_max, continue, end
    rp  = hypot(res.xp(j), res.yp(j));
    dth = atan2(sin(res.thp(j) - res.th(j)), cos(res.thp(j) - res.th(j)));
    if abs(dth) >= pi/2 || rp*abs(sin(dth)) > P.fire.k_gate*P.R_eng, continue, end
    if rp*res.sig_thp(j) > P.fire.k_sig*P.R_eng, continue, end
    tl = t;
    mf = bm(res.th(j), res.fx(j), res.fy(j));  mn = bm(res.th(j), res.tx(j), res.ty(j));
    pulses(end+1,:) = [t, res.th(j), mn, mf, mf <= P.R_eng, mn <= P.R_eng]; %#ok<AGROW>
end
end

function d = bm(th,x,y), s = x*cos(th)+y*sin(th); if s<0, d=hypot(x,y); else, d=abs(-x*sin(th)+y*cos(th)); end, end
