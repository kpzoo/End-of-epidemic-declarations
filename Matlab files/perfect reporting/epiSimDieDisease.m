% Simulate dying epidemic via renewal model
function [Iday, Lam, Rtrue, tday, Iwarn] = epiSimDieDisease(nday, ts, Rs, scenNo, distvals)

% Assumptions and notes
% - simpler method of choosing epidemic than epiSimDie
% - allows selection of serial interval distributions
% - Rs are distinct reprod nums, ts are switch points 

% All scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};

% Variable for true R
Rtrue = zeros(1, nday);
% Change-times in real time
%tch = ceil(nday*ts);
tch = ts;

% Reproduction num profiles
switch(scenNo)
    case 1
        % Rapidly controlled epidemic
        Rtrue(1:tch) = Rs(1);
        Rtrue(tch+1:end) = Rs(2);
    case 2
        % Rapid control that recovers
        Rtrue(1:tch(1)) = Rs(1);
        Rtrue(tch(1)+1:tch(2)) = Rs(2);
        Rtrue(tch(2)+1:end) = Rs(3);
    case 3
        % Two stage control
        Rtrue(1:tch(1)) = Rs(1);
        Rtrue(tch(1)+1:tch(2)) = Rs(2);
        Rtrue(tch(2)+1:end) = Rs(3);
    case 4
        % Exponential rise and fall
        trise = 1:tch; tfall = tch+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tch)); Rmax = Rtrue(tch);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.008*(tfall - tch));
    case 5
        % Second (sine) wave dynamics
        Rtrue = Rs(1) + Rs(2)*sind(ts(1)*(1:nday));
end

% Warning about incidence zeros
Iwarn = 0; % will be set conditionally

% Serial distribution over all tday
Pomega = gammaDistrTypes(nday, distvals);
%Pomega = Pomega/sum(Pomega); % only for exp type

% Daily incidence
Iday = zeros(size(nday)); 
% Infectiousness, Poisson rate 
Lam = Iday; rate = Iday;
% Initialise epidemic
Iday(1) = 10;

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Lam(i) = sum(Iday(i-1:-1:1).*Pomega(1:i-1));
    % Rate for ith day incidence
    rate(i) = Lam(i)*Rtrue(i);
    % Renewal incidence
    Iday(i) = poissrnd(rate(i));
end

% Remove start-up 20 days
idz = 20:nday; tday = idz;
% Adjusted vectors - including tday
Iday = Iday(idz); Rtrue = Rtrue(idz); Lam = Lam(idz);

% Remove small epidemics or very large ones
if sum(Iday) < 500 || max(Iday) > 3000 || max(Iday) < 50
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end

