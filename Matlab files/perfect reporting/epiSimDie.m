% Simulate dying epidemic via renewal model
function [Iday, Lam, Rtrue, tday, Iwarn] = epiSimDie(nday, remGap, ts, Rs, scenNo, distvals)

% Assumptions and notes
% - allows selection of serial interval distributions
% - Rs are distinct reprod nums, ts are switch points 
% - option to remove a startup sequence of zeros
% - warns if sequence of consecutive zero incidences
% - computes incidence and its weekly moving avg

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
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Daily incidence
Iday = zeros(size(nday)); 
% Infectiousness, Poisson rate 
Lam = Iday; rate = Iday;
% Initialise epidemic
Iday(1) = 10;

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    %Pomegat = Pomega(1:i-1);
    % Total infectiousness
    %Lam(i) = Iday(i-1:-1:1)*Pomegat';
    Lam(i) = sum(Iday(i-1:-1:1).*Pomega(1:i-1));
    % Rate for ith day incidence
    rate(i) = Lam(i)*Rtrue(i);
    % Renewal incidence
    Iday(i) = poissrnd(rate(i));
end

% Gaps between non-zero indicence values
zerogaps = diff(find(Iday ~= 0));
% Remove startup sequence of zeros if big
try
    z1 = zerogaps(1);
    if z1 > 5 && remGap
        % Update incidence and related vectors
        idz = z1+1:nday;
        % Flag zero incidence regions after startup
        if max(zerogaps(2:end)) > 5
            %warning('Zero incidences beyond startup');
            Iwarn = 1;
        end
    else
        % Flag any zero incidence region
        if max(zerogaps) > 5
            %warning('Sequences of zero incidence');
            Iwarn = 1;
        end
        % Un-truncated day set
        idz = 2:nday;
    end
    
    % Adjusted vectors - including tday
    clear('tday'); tday = idz;
    Iday = Iday(idz);
    Rtrue = Rtrue(idz); Lam = Lam(idz);
catch
    % Probably only zeros in incidence
    assignin('base', 'Ierr', Iday);
    assignin('base', 'zerogaps', zerogaps);
    % Assign warning to recalculate incidence
    Iwarn = 1; tday = [];
end

