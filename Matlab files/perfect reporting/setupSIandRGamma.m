% Possible gamma SI distributions and R scenarios
function [pm, Rs, ts, omega] = setupSIandRGamma(type, scenNo)

% Assumptions and notes
% - considers only gamma distributions from various diseases

% Shape and mean hyerparameters of generation time
switch(type)
    case 1
        % Marburg distribution from van Kerkhove 2015
        omega = 9; pm = (omega^2)/(5.4^2); 
    case 2
        % MERS distribution from Cauchemez 2016
        omega = 6.8; pm = (omega^2)/(4.1^2);
    case 3
        % Measles distribution from Cori 2013
        omega = 14.9; pm = (omega^2)/(3.9^2);
    case 4
        % COVID-19 distribution from Ferguson 2020
        omega = 6.5; pm = (1/0.65)^2;
end

% Specific scenario parameters
switch(scenNo)
    % Rs are distinct reprod nums, ts are switch points
    case 1
        % Rapidly controlled epidemic
        Rs = [2 0.5];
        %ts = 50; % for MERS, COVID
        ts = 80; % for Marburg, Measles
    case 2
        % Rapid control that slightly recovers
        Rs = [2.5 0.25 0.75];
        %ts = [40 80]; % for MERS, COVID
        ts = [60 100]; % for Marburg, Measles
    case 3
        % Two stage control
        Rs = [5 0.8 0.4];
        ts = [40 120];
    case 4
        % Exponential rise and fall
        %ts = 30; % for MERS, COVID
        ts = 40; % for Marburg, Measles
        Rs = [];
    case 5
        % Second (sine) wave dynamics
        Rs = [1.4 1.3]; ts = 3;
end