% Function to assess end of epidemic stats for a single trajectory
function [tdecs, tdec0, didEnd, z0, z, z1, win, R, tz, idtz, tdec0q, tdecq, tdec1q] = ...
    getSingleEpidFnDisease(Iday, Lam, nday, priors, distvals, Rtrue, k, qs, mu)

% Assumptions and notes
% - quantiles of epidemic declaration time for a given trajectory

% Start checking for epidemic end
tst = find(Iday, 1, 'last');
% Account for epidemic not complete
if tst == length(Iday)
    tst = tst - 200;
end
% Times to interrogate for z
tindex = tst+1:length(Iday); 
lenind = length(tindex);

% Variables for partial epi-curve data
Itemp = cell(1, lenind); Ltemp = Itemp;
z0 = zeros(1, lenind);

% True probability of elimination at t > tst
for i = 1:lenind
    % True data considered
    Itemp{i} = Iday(1:tindex(i)); Ltemp{i} = Lam(1:tindex(i));
    % True elimination probabilities given R
    [z0(i), ~] = getProbEliminTrueDisease(Itemp{i}, Rtrue(tindex(i)), distvals);
end

% Elimination probabilities with estimated R
[~, R, ~, z, z1, tdecs, win] = eliminProbFnDisease(k, Iday, nday, Lam, priors, distvals, lenind, Itemp, mu);
% True declaration time
tdec0 = find(z0 > mu, 1, 'first');
        
% Check epidemic status ensuring it ends
if ~any(tdecs == -1) && tdec0 ~= -1 && max(z) > max(qs)
    % Epidemic ended so flag
    didEnd = 1;
    
    % Estimated declaration quantiles
    quant = @(q) getQuants(q, z); tdecq = arrayfun(quant, qs);
    quant1 = @(q) getQuants(q, z1); tdec1q = arrayfun(quant1, qs);
    % True declaration quantiles
    quant0 = @(q) getQuants(q, z0); tdec0q = arrayfun(quant0, qs);
    
    % Shorten range for plotting
    tz = tindex(1):min(tindex(1)+150, tindex(1)+length(z)-1);
    idtz = 1:length(tz); 
else
    % Default outputs when epidemic ongoing
    tdec0 = -1; tz = []; idtz = []; didEnd = 0;
    tdecq = -ones(size(qs)); tdec1q = tdecq; tdec0q = tdecq;
end

% Simple function to find a quantile of a cdf
function x = getQuants(q, cdfcurve)

% Assumption and notes
% - q in [0, 1] is a quantile probabilites
% - x is the value in cdfcurve corresponding to q

% Quantile is first element crossing desired q
x = find(cdfcurve > q, 1, 'first');

