function m = step_metrics(t, y, yss)
%STEP_METRICS  Overshoot, 10-90 % rise time and 2 % settling time of a step response.
%   m = step_metrics(t, y)        final value yss = y(end)
%   m = step_metrics(t, y, yss)   explicit steady-state value
%   (Octave lacks stepinfo; this keeps the code portable.)
if nargin < 3, yss = y(end); end
t = t(:); y = y(:);
m.yss = yss;
m.overshoot_pct = max(0, (max(y) - yss)/abs(yss)*100);
i10 = find(y >= 0.1*yss, 1, 'first');  i90 = find(y >= 0.9*yss, 1, 'first');
m.rise_time = t(i90) - t(i10);
ib = find(abs(y - yss) > 0.02*abs(yss), 1, 'last');
if isempty(ib), m.settling_time = 0;
elseif ib == numel(t), m.settling_time = NaN;
else, m.settling_time = t(ib+1); end
m.peak_time = t(find(y == max(y), 1, 'first'));
end
