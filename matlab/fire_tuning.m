function T = fire_tuning(P, ncase)
%FIRE_TUNING  Choose the firing-gate constants (k_gate, k_sig) on the TUNING
%   ensemble (random_scenario idx 101..100+ncase), disjoint from scenarios A-C
%   and from the robustness/generalisation cases (idx 201..).
%   The gate does not influence the plant, so each case is simulated ONCE with
%   firing disabled and the gate is replayed offline on the logged signals.
%   Score = mean over cases of (first-hit time, or 10 s penalty if no hit)
%           with a 20 s penalty on every pulse that misses are not counted:
%   we minimise  J = mean( t_first_hit | hit, t_end | no hit ) and report the hit rate.
if nargin < 2, ncase = 16; end
Q = P;  Q.fire.k_gate = 0;                         % disable firing in the sims
R = cell(ncase,1);
for c = 1:ncase
    idx = 100 + c;  wl = mod(c,2) == 0;
    cfg = random_scenario(idx, P, wl);
    R{c} = simulation_runner(cfg, Q, idx);
    fprintf('  tuning case %2d simulated\n', c);
end
kg = [0.15 0.25 0.35 0.5 0.7 1.0];  ks = [0.25 0.5 0.75 1.0 1.5 2.5];
T.kg = kg; T.ks = ks; T.hit = zeros(numel(kg),numel(ks)); T.J = T.hit; T.meanpulses = T.hit;
for a = 1:numel(kg)
  for b = 1:numel(ks)
    Pp = P; Pp.fire.k_gate = kg(a); Pp.fire.k_sig = ks(b);
    hits = 0; J = 0; npl = 0;
    for c = 1:ncase
        pl = replay_gate(R{c}, Pp);
        npl = npl + size(pl,1);
        if any(pl(:,5)), hits = hits + 1; J = J + pl(find(pl(:,5),1),1); else, J = J + P.t_end; end
    end
    T.hit(a,b) = hits/ncase;  T.J(a,b) = J/ncase;  T.meanpulses(a,b) = npl/ncase;
  end
end
end

