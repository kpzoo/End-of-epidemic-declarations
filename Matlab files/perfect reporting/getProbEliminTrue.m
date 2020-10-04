% True probability that epidemic is eliminated
function [z0, zseq0] = getProbEliminTrue(I, Rtrue, distvals)

% Assumptions and notes
% - R is known, assumed constant into future
% - append with pseudo-data of 0s to compute
% - only input data up to point want to check elimination
% - posterior for R based on k days look-back

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
ir = idz:ncurr-1; 

% Probability of elimination sequence
zseq0 = exp(-Rtrue*Lcurr(ir+1));
% Probability of elimination
z0 = prod(zseq0);





