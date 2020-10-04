% Probability an epidemic is eliminated
function [z, zseq, z1, zseq1, Rcurr] = getProbElimin(k, I, priors, distvals)

% Assumptions and notes
% - generates
% - append with pseudo-data of 0s to compute
% - only input data up to point want to check elimination
% - posterior for R based on k days look-back

% Gamma prior parameters (a = shape, b = scale)
a = priors.a; b = priors.b; 

% Define a zero look-ahead sequence
Iz = zeros(1, 100); idz = length(I) + 1;

% Append epidemic curve with pseudo-data
Icurr = [I Iz]; ncurr = length(Icurr);

% Serial distribution over all new times
serial = serialDistrTypes(ncurr, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Compute successive total infectiousness for this I
Lcurr = zeros(1, ncurr);
for i = 2:ncurr
    % Relevant part of SI: Pomega(1:i-1))
    Lcurr(i) = sum(Icurr(i-1:-1:1).*Pomega(1:i-1));
end

% Range of index time points
ir = idz:ncurr-1; lenir = length(ir);
% Grouped incidence and infectiousness
B = zeros(1, lenir); A = B;

% At each time (tPts) compute historical R, predict at tPts+1
for i = 1:lenir
    % Look-back window of k (or less)
    ii = ir(i);
    idback = ii:-1:max(ii-k, 1); 
    % Relevant incidence sum (B) and total infectiousness sum (A)
    B(i) = sum(Icurr(idback)); A(i) = sum(Lcurr(idback));
end

% Parameters of posterior gamma on R
alpha = a + B;
beta = 1./(1/b + A);

% Posterior mean of R and its 99% confidence
Rcurr = alpha.*beta;

% Estimated probability of elimination sequence
zseq = (1 +  Lcurr(ir+1).*Rcurr./alpha).^(-alpha);
% Estimated probability of elimination
z = prod(zseq);

% Estimated elimination sequence by substitution
zseq1 = exp(-Rcurr.*Lcurr(ir+1));
% Probability of elimination by substitution
z1 = prod(zseq1);



