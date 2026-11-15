function P = parameters()
%PARAMETERS  Central parameter definition for the ME319 Assignment 2 project.
%
%   P = parameters() returns a struct with every physical, sensor, estimator,
%   controller and simulation parameter used in the repository.  All values
%   are in SI units (m, s, kg, rad, V, A, ohm, H) unless a comment states
%   otherwise.  Every value is tagged below as
%       [ASSIGNMENT]  -- given in the assignment PDF (Sections 2, 3, 4),
%       [DERIVED]     -- computed from [ASSIGNMENT] values,
%       [DESIGN]      -- a design choice made by the author (justified in the
%                        paper; NOT given by the assignment),
%       [ASSUMPTION]  -- a modelling assumption not stated in the assignment.
%
%   Project : ME319 Control Systems (Autumn 2026) - Assignment 2
%   Authors : Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), IIT Bombay
%   Scope   : non-weaponised predictive tracking-platform simulation.

% -------------------------------------------------------------------------
% 1. Elevation-drive physical parameters  [ASSIGNMENT]
% -------------------------------------------------------------------------
P.N        = 100;      % gear ratio, omega_motor = N*omega_turret          [-]
P.J_turret = 1200;     % platform + cradle inertia about trunnion          [kg m^2]
P.J_m      = 0.20;     % motor rotor inertia                               [kg m^2]
P.b        = 1500;     % platform-side viscous damping                     [N m s/rad]
P.Kt       = 10;       % motor torque constant                             [N m/A]
P.Ke       = 10;       % back-emf constant                                 [V s/rad]
P.Rm       = 1.0;      % armature resistance                               [ohm]
P.Lm       = 0.10;     % armature inductance                               [H]
P.theta_min = -10*pi/180;   % mechanical elevation limit (lower)           [rad]
P.theta_max =  85*pi/180;   % mechanical elevation limit (upper)           [rad]

% -------------------------------------------------------------------------
% 2. Derived plant quantities  [DERIVED]
% -------------------------------------------------------------------------
P.J_eq = P.J_turret + P.N^2 * P.J_m;                  % reflected inertia  [kg m^2]
P.Kv   = P.N*P.Kt / (P.Rm*P.b + P.N^2*P.Kt*P.Ke);     % low-freq. velocity gain [rad/s per V]
P.wn   = sqrt((P.Rm*P.b + P.N^2*P.Kt*P.Ke)/(P.Lm*P.J_eq));    % electro-mech. resonance [rad/s]
P.zeta = (P.Rm*P.J_eq + P.Lm*P.b)/(2*P.Lm*P.J_eq*P.wn);       % its damping ratio [-]

% -------------------------------------------------------------------------
% 3. Sensor / target / scenario parameters
% -------------------------------------------------------------------------
P.Ts_radar    = 0.05;           % sensor update interval (20 Hz)       [s]   [ASSIGNMENT]
P.sigma_r     = 1.0;            % range noise std                      [m]   [ASSIGNMENT]
P.sigma_theta = 0.02*pi/180;    % elevation noise std                  [rad] [ASSIGNMENT]
P.r_min  = 800;                 % operating region                     [m]   [ASSIGNMENT]
P.r_max  = 2000;                % operating region                     [m]   [ASSIGNMENT]
P.v_max  = 150;                 % speed bound                          [m/s] [ASSIGNMENT]
P.a_max  = 20;                  % acceleration bound                   [m/s^2][ASSIGNMENT]
P.t_end  = 10;                  % scenario duration                    [s]   [ASSIGNMENT]
P.Tp     = 0.50;                % prediction (look-ahead) horizon      [s]   [ASSIGNMENT]
P.t_bo0  = 1.00;                % measurement-loss onset (scenario C)  [s]   [ASSIGNMENT]
P.t_bo   = 0.50;                % measurement-loss duration (scenario C)[s]  [ASSIGNMENT]
P.seed   = 24240;               % noise seed (reproducibility)         [-]   [DESIGN]

% -------------------------------------------------------------------------
% 4. Estimator / predictor design parameters  [DESIGN]
%    Model: independent Singer-type acceleration model on each Cartesian axis
%       d/dt [p v a]' = [0 1 0; 0 0 1; 0 0 -alpha][p v a]' + [0 0 1]' w,
%       E[w(t) w(s)] = q delta(t-s)
% -------------------------------------------------------------------------
P.est.q      = 200;   % jerk power spectral density          [m^2/s^5]  (from design_study.m)
P.est.alpha  = 2.0;   % acceleration decay rate (1/tau_a)    [1/s]      (from design_study.m)
                      % => stationary accel. std sqrt(q/(2 alpha)) = 7.1 m/s^2, tau_a = 0.5 s
P.est.sig_v0 = P.v_max/sqrt(3);   % prior velocity std: uniform on [-vmax,vmax]  [m/s]
P.est.sig_a0 = P.a_max/sqrt(3);   % prior accel. std:   uniform on [-amax,amax]  [m/s^2]

% -------------------------------------------------------------------------
% 5. Feedback-controller design parameters  [DESIGN]
%    Gc(s) = Kp (s + wi)/s,  Kp = wc/Kv   (analytic; see paper Sec. VI)
% -------------------------------------------------------------------------
P.ctrl.wc    = 5.0;                      % target gain-crossover frequency  [rad/s]
P.ctrl.ratio = 3.0;                      % wc/wi (PI zero at one third of wc) [-]
P.ctrl.wi    = P.ctrl.wc/P.ctrl.ratio;   % PI zero                          [rad/s]
P.ctrl.Kp    = P.ctrl.wc/P.Kv;           % proportional gain                [V/rad]
P.ctrl.Ki    = P.ctrl.Kp*P.ctrl.wi;      % integral gain                    [V/(rad s)]
P.ctrl.V_sat = 1500;     % supply-voltage limit [V]  [ASSUMPTION: the assignment gives none;
                         % chosen so that it stays inactive in scenarios A-C]

% -------------------------------------------------------------------------
% 6. Simulation settings  [DESIGN]/[ASSUMPTION]
% -------------------------------------------------------------------------
P.sim.dt       = 1e-3;                       % integration and controller step (1 kHz) [s]
P.sim.n_rad    = round(P.Ts_radar/P.sim.dt); % simulation steps per radar sample       [-]
P.sim.align0   = true;   % platform starts aligned with true initial target elevation [ASSUMPTION]
P.sim.t_settle = 1.0;    % acquisition transient excluded from "post-acquisition" metrics [s]

% -------------------------------------------------------------------------
% 7. Scenario definitions (smooth bounded manoeuvres)  [DESIGN]
%    a(t) = A sin(w (t - t_man) + ph) * smoothstep((t - t_man)/t_ramp), t >= t_man
% -------------------------------------------------------------------------
base = struct('p0',[1800;500], 'v0',[-90;0], ...
              'Ax',0,'wx',0,'phx',0,'Ay',0,'wy',0,'phy',0, ...
              't_man',1.0,'t_ramp',0.5,'t_bo0',NaN,'t_bo',0);
P.scen.A = base;                                  % nominal: constant velocity
B = base;
B.Ax = 12; B.wx = 2*pi/6.0; B.phx = 0;            % x acceleration amplitude 12 m/s^2, period 6 s
B.Ay = 14; B.wy = 2*pi/4.5; B.phy = pi/2;         % y acceleration amplitude 14 m/s^2, period 4.5 s
P.scen.B = B;                                     % smooth bounded manoeuvre
C = B;
C.t_bo0 = P.t_bo0; C.t_bo = P.t_bo;               % same manoeuvre + measurement loss
P.scen.C = C;

% -------------------------------------------------------------------------
% 8. Engagement (laser) model and firing logic  [ASSIGNMENT]/[DESIGN]
% -------------------------------------------------------------------------
P.R_eng      = 2.0;      % engagement radius                              [m]  [ASSIGNMENT]
P.t_recharge = 1.0;      % minimum time between pulses                    [s]  [ASSIGNMENT]
P.eng.delay  = P.Tp;     % pulse is scored against the target position at t_fire + delay.
                         % [ASSUMPTION] delay = Tp: the turret is commanded to the
                         % predicted position at t+Tp, i.e. the pulse "arrives" Tp later.
                         % The alternative reading (score at t_fire) is also reported.
P.fire.k_gate  = 0.5;    % fire only if estimated miss <= k_gate*R_eng       [-]  [DESIGN, tuned on idx 101-116]
P.fire.k_sig   = 1.5;    % ... and r_p*sigma_theta_p <= k_sig*R_eng          [-]  [DESIGN, tuned on idx 101-116]
P.fire.t_arm   = 0.5;    % no firing before the filter has run this long     [s]  [DESIGN]
end
