% Plot SI distributions and R profiles% Run many epidemics and get declarations across all scenarios
clearvars; clc;
close all; tic;

% Assumptions and notes
% - plots 3 R trajectories and SI distributions

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);
% New colour choices
cmap = linspecer(10, 'sequential');

% Save data and test figs
saveTrue = 1; thisDir = cd;

%% Choose distributions and R profiles

% Generation time distributions to plot
distNo = 1:3; nDist = length(distNo);
% Define all SI/generation time distributions - not bimodal
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
% Mean of SI distribution
distvals.omega = 14.2;

% Store probabilities
Pomega = cell(1, nDist); Pmean = zeros(1, nDist);
% For every chosen distribution get parameter
for i = 1:nDist
    % Set serial interval parameters
    distvals.type = distNo(i);
    
    % Define SI distribution (argument of 1 is dummy)
    [distvals.pm, ~, ~] = setupSIandR(distNo(i), 1);
    
    % Domain of generation time distribution
    tdist = 1:100;
    % Generation time distribution values
    serial = serialDistrTypes(max(tdist), distvals);
    % Single omega (mean) controls distribution
    Pomega{i} = serial(1/distvals.omega);
    
    % Mean of distribution
    Pmean(i) = sum(Pomega{i}.*tdist);
end

% R profile scenarios to plot
scenNo = [1 2 4]; nScen = length(scenNo);
% Define all dying epidemic scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};

% Store R changepoints and plotting values
Rch = cell(1, nScen); tch = Rch; t = Rch; R = Rch;
% Starting times (days)
tday0 = 1:300; nday0 = length(tday0);

% For every chosen R profile get trajectory
for i = 1:nScen
    % Define R parameters (argument of 1 is dummy)
    [~, Rch{i}, tch{i}] = setupSIandR(1, scenNo(i));
    
    % R trajectory
    [~, ~, R{i}, t{i}, ~] = epiSimDie(length(tday0), 1, tch{i}, Rch{i}, scenNo(i), distvals);
end


% Figure with 3 distributions and R profiles
if nDist == 3 && nScen == 3
    figure; 
    % Top row of all SI
    for i = 1:nDist
        subplot(2, 3, i);
        plot(tdist, Pomega{i}, 'LineWidth', 2);
        xlabel('$u$ (days)');
        ylabel('$w_u$');
        grid off; box off;
    end
    % Bottom row is all R
    for i = 1:nScen
        subplot(2, 3, nDist+i);
        plot(t{i}, R{i}, 'LineWidth', 2);
        xlabel('$s$ (days)');
        ylabel('$R_s$');
        grid off; box off;
    end
end
% Save composite figure
saveas(gcf, 'RandSI', 'fig');


