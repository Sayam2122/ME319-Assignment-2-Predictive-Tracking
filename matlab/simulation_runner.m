function res = simulation_runner(scenario, P, seed)
%SIMULATION_RUNNER  Closed-loop simulation of one scenario (0 <= t <= t_end).
%
%   res = simulation_runner(scenario, P)
%   res = simulation_runner(scenario, P, seed)
%
%   Loop (every P.sim.dt = 1 ms):
%     1. target_motion_generator -> noisy radar sample (new only every 50 ms)
%     2. turret_controller       -> motor voltage V_m (EKF + prediction + PI)
%     3. turret_plant            -> RK4 step of [theta, theta_dot, i] with V_m
%                                   held constant over the step (ZOH)
%     4. mechanical end-stops    -> theta clipped to [theta_min, theta_max],
%                                   outward velocity set to zero
%   The controller only sees (r_meas, theta_meas, meas_valid, theta, theta_dot);
%   the ground truth is used only for logging/evaluation.
%
%   Output res: struct of column vectors (length n+1) with fields
%     t, V, I (current), th, om, th_ref, e, valid, rm, thm,
%     tx, ty, tvx, tvy, tth, (truth now)     fx, fy, fth (truth at t+Tp),
%     xh, yh, vxh, vyh, thh (estimate now),  xp, yp, thp (prediction at t+Tp),
%     sig_thp (angle-prediction 1-sigma), sig_pos (position 1-sigma),
%     nis (innovation NIS at updates, NaN otherwise), scenario, seed.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

if nargin < 3, seed = P.seed; end
dt = P.sim.dt;
n  = round(P.t_end/dt);
nn = n + 1;

% NOTE: turret_controller resets itself at t == 0 and the generator rebuilds its
% record whenever (scenario, seed) changes, so no 'clear' is needed here.

names = {'t','V','I','th','om','th_ref','e','valid','rm','thm', ...
         'tx','ty','tvx','tvy','tth','fx','fy','fth', ...
         'xh','yh','vxh','vyh','thh','xp','yp','thp','sig_thp','sig_pos','nis'};
for k = 1:numel(names), res.(names{k}) = nan(nn,1); end

% ---- initial condition: platform aligned with initial target elevation ----
[~,~,~,tr0] = target_motion_generator(0, scenario, seed);
if P.sim.align0, th = tr0.theta; else, th = 0; end
om = 0;  cur = 0;
res.fire = zeros(nn,1);  pulses = zeros(0,6);   % [t_fire theta miss_now miss_future ok_future ok_now]
t_last_pulse = -Inf;  fire_prev = 0;

for k = 0:n
    t = k*dt;
    [rm, thm, valid, tr] = target_motion_generator(t, scenario, seed);
    [~,~,~,trf] = target_motion_generator(t + P.Tp, scenario, seed);

    [V, fc, dbg] = turret_controller(t, rm, thm, valid, th, om);

    j = k+1;
    res.t(j)=t; res.V(j)=V; res.I(j)=cur; res.th(j)=th; res.om(j)=om;
    res.th_ref(j)=dbg.theta_ref; res.e(j)=dbg.e; res.valid(j)=valid;
    res.rm(j)=rm; res.thm(j)=thm;
    res.tx(j)=tr.x; res.ty(j)=tr.y; res.tvx(j)=tr.vx; res.tvy(j)=tr.vy; res.tth(j)=tr.theta;
    res.fx(j)=trf.x; res.fy(j)=trf.y; res.fth(j)=trf.theta;
    if ~isempty(dbg.pred)
        pr = dbg.pred;
        res.xh(j)=pr.x_hat; res.yh(j)=pr.y_hat; res.vxh(j)=pr.vx_hat; res.vyh(j)=pr.vy_hat;
        res.thh(j)=pr.theta_hat; res.xp(j)=pr.x_p; res.yp(j)=pr.y_p; res.thp(j)=pr.theta_p;
        if mod(k, 10) == 0                          % uncertainty every 10 ms (costly)
            pc = target_predictor(dbg.S, t, P, true);
            res.sig_thp(j) = pc.sig_theta_p;  res.sig_pos(j) = pc.sig_pos_hat;
        end
        if valid, res.nis(j) = dbg.S.nis; end
    end

    % ---- laser pulse bookkeeping (0->1 edge, recharge rule of the assignment) ----
    res.fire(j) = fc;
    if fc && ~fire_prev && (t - t_last_pulse) >= P.t_recharge - 1e-9
        t_last_pulse = t;
        mn = beam_miss(th, tr.x,  tr.y);
        mf = beam_miss(th, trf.x, trf.y);
        pulses(end+1,:) = [t, th, mn, mf, mf <= P.R_eng, mn <= P.R_eng]; %#ok<AGROW>
    end
    fire_prev = fc;

    if k == n, break; end

    % ---- RK4 step of the plant (V held constant) ----
    x0 = [th; om; cur];
    f  = @(x) cell2mat_plant(t, V, x);
    k1 = f(x0);  k2 = f(x0 + 0.5*dt*k1);
    k3 = f(x0 + 0.5*dt*k2);  k4 = f(x0 + dt*k3);
    x1 = x0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    th = x1(1);  om = x1(2);  cur = x1(3);

    % ---- mechanical end-stops ----
    if th >= P.theta_max, th = P.theta_max; om = min(om, 0); end
    if th <= P.theta_min, th = P.theta_min; om = max(om, 0); end
end
res.pulses = pulses;  res.scenario = scenario;  res.seed = seed;
end

function d = beam_miss(th, x, y)
%BEAM_MISS  Minimum distance between the beam (ray from the origin along th)
%   and the point (x,y).  Behind the origin the distance is that to the origin.
sproj = x*cos(th) + y*sin(th);
if sproj < 0, d = hypot(x, y); else, d = abs(-x*sin(th) + y*cos(th)); end
end

function xd = cell2mat_plant(t, V, x)
[a, b, c] = turret_plant(t, V, x(1), x(2), x(3));
xd = [a; b; c];
end
