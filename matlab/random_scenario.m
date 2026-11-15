function cfg = random_scenario(idx, P, with_loss)
%RANDOM_SCENARIO  Deterministic pseudo-random bounded manoeuvre (scenario struct).
%
%   cfg = random_scenario(idx, P, with_loss)
%
%   Draws start position, velocity and smooth bounded x/y manoeuvres
%   (sinusoidal acceleration windowed by a smoothstep, as in parameters.m)
%   from a Park-Miller stream seeded by idx, and REJECTS draws that violate
%   the range, speed, acceleration or elevation-range limits (by calling
%   target_motion_generator, which enforces them).  Optionally adds a sensor
%   loss (onset in [1, 8] s, duration in [0.1, 0.5] s).
%   Used by design_study.m (estimator tuning ensemble) and
%   robustness_study.m (generalisation test) with disjoint idx ranges.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

m = 2147483647; a = 16807;
s = mod(abs(round(idx))*7919 + 104729, m-1) + 1;
for k = 1:20, s = mod(a*s, m); end
U = zeros(600,1);
for k = 1:600, s = mod(a*s, m); U(k) = s/m; end
uni = @(p, lo, hi) lo + (hi-lo)*p;

for att = 0:36
    u = U(att*16 + (1:16));
    rad   = uni(u(1), 1100, 1900);                 % initial range [m]
    el    = uni(u(2), 0, 55)*pi/180;               % initial elevation
    spd   = uni(u(3), 30, 130);                    % speed [m/s]
    head  = el + pi + uni(u(4), -60, 60)*pi/180;   % broadly inbound heading
    amp   = uni(u(5), 8, 19);                      % total accel amplitude (< a_max)
    spl   = uni(u(6), 0, 2*pi);                    % split between x and y
    cfg = struct('p0',[rad*cos(el); rad*sin(el)], 'v0',spd*[cos(head); sin(head)], ...
        'Ax',amp*cos(spl), 'wx',2*pi/uni(u(7),1.5,8), 'phx',uni(u(9),0,2*pi), ...
        'Ay',amp*sin(spl), 'wy',2*pi/uni(u(8),1.5,8), 'phy',uni(u(10),0,2*pi), ...
        't_man',uni(u(11),1,5), 't_ramp',0.5, 't_bo0',NaN, 't_bo',0);
    if with_loss
        cfg.t_bo0 = uni(u(12), 1, 8);  cfg.t_bo = uni(u(13), 0.1, 0.5);
    end
    try
        target_motion_generator(0, cfg, 1);          % throws if limits violated
        ok = true;
        for tt = 0:0.5:(P.t_end + P.Tp)              % elevation within mechanical range
            [~,~,~,tr] = target_motion_generator(tt, cfg, 1);
            if tr.theta < P.theta_min + 2*pi/180 || tr.theta > P.theta_max - 2*pi/180
                ok = false; break
            end
        end
        if ok, return, end
    catch
    end
end
error('random_scenario:noValidDraw', 'No admissible draw for idx = %d.', idx);
end
