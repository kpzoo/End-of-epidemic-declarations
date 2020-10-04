% Possible SI distributions and R scenarios
function [pm, Rs, ts] = setupSIandR(type, scenNo)

% Hyerparameters of generation time
switch(type)
    case 1
        % Geometric distribution - no parameter
        pm = [];
    case 2
        % Gamma distribution - shape parameter
        pm = 20;
    case 3
        % Delta distribution - odd window around mean
        pm = 7;
    case 4
        % Two Gamma distributions for flare-up
        pm = 45;
end

% Specific scenario parameters
switch(scenNo)
    % Rs are distinct reprod nums, ts are switch points
    case 1
        % Rapidly controlled epidemic
        Rs = [2 0.5];
        ts = 100;
    case 2
        % Rapid control that slightly recovers
        %Rs = [2.5 0.4 0.8];
        Rs = [2.5 0.25 0.75];
        ts = [80 110];
    case 3
        % Two stage control
        Rs = [5 0.8 0.4];
        ts = [40 120];
    case 4
        % Exponential rise and fall
        ts = 40;
        Rs = [];
    case 5
        % Second (sine) wave dynamics
        Rs = [1.4 1.3]; ts = 3;
end