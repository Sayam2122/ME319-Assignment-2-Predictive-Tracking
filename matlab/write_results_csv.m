function write_results_csv(M, rob, datadir)
%WRITE_RESULTS_CSV  Plain-text export of the simulation metrics.
%   The paper is written by hand; its numbers are copied from these CSV files
%   (results/data/metrics_summary.csv and robustness_cases.csv).
%   This function writes no LaTeX.

S = {'A','B','C'};
fid = fopen(fullfile(datadir, 'metrics_summary.csv'), 'w');
fprintf(fid, 'scenario,window,metric,value\n');
for s = 1:3
    m = M.(S{s});
    for g = {'full','post'}
        fl = fieldnames(m.(g{1}));
        for k = 1:numel(fl), fprintf(fid, '%s,%s,%s,%.6g\n', S{s}, g{1}, fl{k}, m.(g{1}).(fl{k})); end
    end
    fprintf(fid, '%s,eng,success,%d\n%s,eng,t_int,%.6g\n%s,eng,n_pulses,%d\n', S{s}, m.eng.success, S{s}, m.eng.t_int, S{s}, m.eng.n_pulses);
    fprintf(fid, '%s,all,V_peak,%.6g\n%s,all,I_peak,%.6g\n%s,all,om_peak_deg_s,%.6g\n%s,all,sat_frac,%.6g\n', ...
            S{s}, m.V_peak, S{s}, m.I_peak, S{s}, m.om_peak, S{s}, m.sat_frac);
end
if isfield(M.C, 'loss')
    fl = fieldnames(M.C.loss);
    for k = 1:numel(fl)
        v = M.C.loss.(fl{k}); fprintf(fid, 'C,loss,%s,%.6g\n', fl{k}, v(1));
    end
end
fclose(fid);


% robustness (random unseen cases): per-case engagement results
if ~isempty(rob) && isfield(rob, 'eng')
    fid = fopen(fullfile(datadir, 'robustness_cases.csv'), 'w');
    fprintf(fid, 'case,with_loss,success,t_int_s,miss_at_hit_m,n_pulses,success_if_scored_at_tfire,hit_fraction\n');
    for c = 1:size(rob.eng, 1)
        fprintf(fid, '%d,%d,%d,%.4g,%.4g,%d,%d,%.4g\n', c, rob.with_loss(c), rob.eng(c,1), rob.eng(c,2), ...
                rob.eng(c,3), rob.eng(c,4), rob.eng(c,5), rob.eng(c,7));
    end
    fclose(fid);
end
end
