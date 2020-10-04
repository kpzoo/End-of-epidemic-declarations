function pdistr = gammaDistrTypes(tmax, distvals)

% Assumptions and notes
% - Chooses between discrete gamma distributions on days
% - Insert max days and calculate serial probabilities
% - p must be a parameter, tday an array of integers
% - m is mean of each distribution
% - various diseases considered

% Gamma distribution with integer shape (Erlang)
pdistr = gammaDistr(distvals.omega, 1:tmax, distvals.pm);

%% Gamma distribution, p is 1/mean, shape param input
function pr = gammaDistr(m, x, shapePm)

% shapePm is shape parameter of gamma
% scalePm is computed from shape and mean

% Scale parameter based on mean = shapePm*scalePm
scalePm = m/shapePm; ratePm = 1/scalePm;

% Gamma (Erlang) probabilities
pr = -log(gamma(shapePm)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(x) - ratePm*x;
pr = exp(pr);

