function [theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, ...
                                              theta_turret, theta_dot_turret, i)
%TURRET_PLANT  State derivatives of the geared DC-motor elevation drive.
%
%   [theta_dot, theta_ddot, i_dot] = turret_plant(t, V_m, theta, theta_dot, i)
%
%   States : theta [rad] platform elevation, theta_dot [rad/s], i [A] armature current
%   Input  : V_m [V] armature voltage
%   Outputs: d/dt of the three states (theta_dot is returned as the first
%            output, i.e. the time derivative of theta).
%
%   Equations (derivation in paper Sec. III; N = gear ratio, motor speed
%   = N*theta_dot, J_eq = J_turret + N^2 J_m):
%       Lm di/dt      = V_m - Rm i - Ke N theta_dot          (armature loop)
%       J_eq theta'' = N Kt i - b theta_dot                 (reflected mechanics)
%   Assumptions: rigid gear train (no backlash/compliance), ideal gearing,
%   balanced platform (no gravity torque), linear viscous damping, constant
%   parameters.  The mechanical end-stops are enforced in simulation_runner.m,
%   not here, so that this function remains the pure linear model of the
%   assignment.  Arguments t and theta_turret are unused (LTI model) and kept
%   to match the required interface.
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

persistent Pc
if isempty(Pc), Pc = parameters(); end

theta_dot_turret  = theta_dot_turret;                       % d(theta)/dt
theta_ddot_turret = (Pc.N*Pc.Kt*i - Pc.b*theta_dot_turret) / Pc.J_eq;
i_dot             = (V_m - Pc.Rm*i - Pc.Ke*Pc.N*theta_dot_turret) / Pc.Lm;
end
