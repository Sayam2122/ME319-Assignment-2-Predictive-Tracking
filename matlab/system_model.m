function M = system_model(P)
%SYSTEM_MODEL  Linear model of the geared DC-motor elevation drive.
%
%   M = system_model(P) builds, from the physical parameters in parameters.m,
%       M.Gp      voltage -> platform angle   Theta(s)/V(s)           (tf)
%       M.Gw      voltage -> platform speed   Omega(s)/V(s)           (tf)
%       M.sys     state-space model, x = [theta; omega; i], u = V_m, y = theta
%       M.num/den coefficients of Gp
%       M.poles   open-loop poles
%
%   Derivation (see paper Sec. III).  With motor speed omega_m = N*omega:
%     armature loop :  Lm di/dt + Rm i + Ke N omega = V_m
%     motor torque  :  T_m = Kt i                    (at the motor shaft)
%     gear (ideal)  :  T_load = N T_m ,  omega_m = N omega
%     mechanics     :  J_turret omega' + b omega = N T_m - N^2 J_m omega'
%                      => (J_turret + N^2 J_m) omega' + b omega = N Kt i
%   Laplace transform (zero initial conditions):
%     (Lm s + Rm) I = V - Ke N s Theta ;  (J_eq s + b) s Theta = N Kt I
%   =>  Theta/V = N Kt / ( s [ (Lm s + Rm)(J_eq s + b) + N^2 Kt Ke ] )
%             = N Kt / ( s [ Lm J_eq s^2 + (Rm J_eq + Lm b) s + (Rm b + N^2 Kt Ke) ] ).
%   Valid for small perturbations about any operating angle because the
%   model has no angle-dependent term (balanced platform, linear friction).
%
%   Authors: Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194), ME319 Autumn 2026.

N = P.N; Jeq = P.J_eq;
M.num = N*P.Kt;
M.den = [P.Lm*Jeq, P.Rm*Jeq + P.Lm*P.b, P.Rm*P.b + N^2*P.Kt*P.Ke, 0];
M.Gp  = tf(M.num, M.den);
M.Gw  = tf(M.num, M.den(1:3));            % Omega/V = s * Theta/V

A = [0, 1, 0; ...
     0, -P.b/Jeq,           N*P.Kt/Jeq; ...
     0, -P.Ke*N/P.Lm,      -P.Rm/P.Lm];
B = [0; 0; 1/P.Lm];
C = [1 0 0];  D = 0;
M.A = A; M.B = B; M.C = C; M.D = D;
M.sys   = ss(A, B, C, D);
M.poles = roots(M.den);
end
