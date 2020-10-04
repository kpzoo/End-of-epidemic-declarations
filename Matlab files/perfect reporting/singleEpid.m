% Assess the elimination probability of a single epidemic
clearvars; clc;
close all; tic;

% Assumptions and notes
% - declarations are relative to tst; when last case was seen
% - computes estimated probability with NB uncertainty (z) or mean (z1)
% - computes true probability of elimination from known R (z0)
% - simulates dying epidemics and computes from first consistent 0
% - at each t compute zero sequence prob in future
% - computes p(I_{t+h} = 0 | I_{t+1}^{t+h-1} = 0, I_1^t)

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);
% New colour choices
cmap = linspecer(10);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving
saveFol = 'single data'; thisDir = cd;

% Confidence level for declaration
mu = 0.95;
disp(['Confidence = ' num2str(mu)]);

%% Simulate epidemic with decreasing Rt 

% Choose a scenario
scenNo = 1;
% Choose a generation time distribution
distNo = 4;
% Show whole epi-curve
wholeEp = 1;

% Define all dying epidemic scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};
scenChoice = scenNam{scenNo};
disp(['True R scenario: ' scenNam{scenNo}]);

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo};
disp(['True serial interval: ' distNam{distNo}]);

% Set serial interval parameters
distvals.type = distNo; 
distvals.omega = 14.2;

% Define SI distribution and R scenario
[distvals.pm, Rs, ts] = setupSIandR(distNo, scenNo);

% Starting times (days)
tday0 = 1:501; nday0 = length(tday0);
tcheck = max(ts);

% Simulate epidemic scenarios and truncate
Iwarn = 1; % ensure no warnings
while Iwarn 
    [Iday, Lam, Rtrue, tday, Iwarn] = epiSimDie(length(tday0), 1, ts, Rs, scenNo, distvals);
    % Ensure epidemic died but not prematurely
    Idayfir = Iday(1:tcheck); Idaysec = Iday(tcheck+1:end);
end
% Truncated observation period based on removed 0s
nday = length(tday);

% Window size (k+1 with current value)
k = 100; disp(['Window length: ' num2str(k+1)]);

% Set priors on R estimates, E[R] = ab
priors.a = 1; priors.b = 5;

%% Characterise elimination probability

% Start checking for epidemic end
tst = find(Iday, 1, 'last');
% If want to instead show whole epi-curve
if wholeEp
    % Now end time is absolute
    tst = 1;
end

% Times to interrogate for z
tindex = tst+1:length(Iday); 
lenind = length(tindex);

% Elimination probabilities and R
Rseq = cell(1, lenind); zseq = Rseq; z = zeros(1, lenind);
z0 = z; zseq0 = zseq; z1 = z; zseq1 = zseq0;

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered
    Itemp = Iday(1:tindex(i));
    Ltemp = Lam(1:tindex(i));
    
    % Estimates of elimination probability
    [z(i), zseq{i}, z1(i), zseq1{i}, Rseq{i}] = getProbElimin(k, Itemp, priors, distvals);
    % True elimination probabilities given R
    [z0(i), zseq0{i}] = getProbEliminTrue(Itemp, Rtrue(tindex(i)), distvals);
end
% Closeness of z curves
se = mean((z - z0).^2); se1 = mean((z1 - z0).^2);
disp(['SE [z z1]: [' [num2str(se) ' ' num2str(se1)] ']' ]);

% Find points of declaration
tdecs = [find(z0 > mu, 1, 'first'), find(z > mu, 1, 'first'), find(z1 > mu, 1, 'first')];
disp(['Declared at [z0 z z1]: ' num2str(tdecs)]);

% Check if epidemic actually ended
if any(isempty(tdecs))
    didEnd = 0;
    disp('Epidemic did not end');
else
    didEnd = 1;
end

% Domains of distributions
tdist = 1:100; Rdist = linspace(0, 30, 10000);
% Generation time distribution
serial = serialDistrTypes(max(tdist), distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);
% Priors on R
prior_R = gampdf(Rdist, priors.a, priors.b);

%% Visualisation and processing

% Plot from second time onwards
tplt = tday(2:end); Rplt = Rtrue(2:end); Iplt = Iday(2:end);
% Maximum time id
idmax = find(tindex == max(tdecs)+30);


% Generation time and R prior
figure;
subplot(2, 1, 1);
plot(tdist, Pomega, 'LineWidth', 2);
xlabel('$s$ (days)');
ylabel('$w_s$');
grid off; box off;
subplot(2, 1, 2);
plot(Rdist, prior_R, 'LineWidth', 2);
xlabel('$R_s$');
ylabel('$P(R_s)$');
grid off; box off;

% Elimination probabilities with time
if didEnd
    figure;
    plot(z, 'Color', 'r', 'LineWidth', 2);
    hold on;
    plot(z0, 'Color', 'k', 'LineWidth', 2);
    plot(1:length(z), mu*ones(size(z)), 'c--', 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$\Delta s$ (days)');
    ylabel('$z_s$');
    xlim([1 length(z)]);
    legend('$z$', '$z_0$', 'Location', 'Best');
    ylim([0 1.1]);
end

if wholeEp && didEnd
    % Examine how z changes when epidemic does not end
    figure;
    yyaxis right
    stairs(tindex, z0, '-', 'Color', grey2, 'LineWidth', 2);
    hold on;
    stairs(tindex, z, '-', 'Color', grey1, 'LineWidth', 2);
    plot(1:length(z), mu*ones(size(z)), 'c--', 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$s$ (days)');
    ylabel('$z_s$');
    xlim([tindex(1) max(tdecs)+30]);
    h = gca; h.YColor = h.XColor;
    h.YLim = [0 1.05];
    yyaxis left
    stairs(tindex(1:idmax), Iday(tindex(1:idmax)), 'Color', cmap(2,:), 'LineWidth', 2);
    grid off; box off;
    ylabel('$I_s$');
    h = gca; h.YColor = h.XColor;
    xlim([tindex(1) max(tdecs)+30]);
    %legend('$z_0$', '$z_1$', 'Location', 'Best');
    if saveTrue
        cd(saveFol);
        saveas(gcf, ['epEx_' scenChoice '_' num2str(k)], 'fig');
        cd(thisDir);
    end 
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save([scenChoice '_' num2str(nday) '_' num2str(k) '.mat']);
    cd(thisDir);
end