% Assessing the elimination probability of an epidemic
function [ape, R, RInt, z, z1, tdecs, win] = eliminProbFnDisease(k, Iday, nday, Lam, priors, distvals, lenind, Itemp, mu)

% Assumptions and notes
% - estimated elimin probs on fixed R trajectory
% - computes p(I_{t+h} = 0 | I_{t+1}^{t+h-1} = 0, I_1^t)
% - treats ending sequence of zeros as new pseudo-data

% Statistics of epidemic, Rt estimates and It predictions
[~, ~, prob, R, RInt, win] = getNegBinPriors(k, nday, Iday, Lam, priors);
% APE metric
ape = -sum(log(prob(prob ~= 0)));

% Elimination probabilities and R
z = zeros(1, lenind); z1 = z;

% For t > tst get elimination probability
for i = 1:lenind
    % Estimated elimination probability with pseudo-data    
    [z(i), ~, z1(i), ~, ~] = getProbEliminDisease(k, Itemp{i}, priors, distvals);
end

% Find points of declaration
tdecs = [find(z > mu, 1, 'first'), find(z1 > mu, 1, 'first')];
if any(isempty(tdecs))
    % Epidemic did not end
    tdecs = [-1 -1];
end

