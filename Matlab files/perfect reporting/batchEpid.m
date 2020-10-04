% Run many epidemics and get declarations
clearvars; clc;
close all; tic;

% Assumptions and notes
% - R trajectory is fixed but various I trajectories
% - compare true tdec with estimates
% - fix best k, prior and examine different generation times

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving
saveFol = 'batch data'; thisDir = cd;

% Confidence level for declaration
mu = 0.95;

%% Setup epidemic parameters

% Choose scenario
scenNo = 4;
% Choose a generation time distribution
distNo = 2;
% Define all dying epidemic scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};
scenChoice = scenNam{scenNo};

% Window size and runs
k = 99; M = 1000;
% Set priors on R estimates, E[R] = ab
priors.a = 1;  priors.b = 5;

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo};
disp(['SI distribution: ' distNam{distNo}]);

% Set serial interval parameters
distvals.type = distNo; 
distvals.omega = 14.2;

% Define SI distribution and R scenario
[distvals.pm, Rs, ts] = setupSIandR(distNo, scenNo);

% Starting times (days)
tday0 = 1:601; nday0 = length(tday0);
tcheck = max(ts);

% Declaration quantiles of interest
qs = [0.05:0.05:0.95 0.99];
nQuant = length(qs);

%% Simulate M epidemics and find declaration time

% Store variables from each epidemic
Iday = cell(1, M); Lam = Iday; Rtrue = Iday; tday = Iday;
%z = Iday; z0 = Iday; z1 = Iday; R = Iday; tz = Iday; idtz = Iday;
nday = zeros(1, M); tdec0 = nday; win = nday; didEnd = nday;
tdecs = zeros(M, 2); tdecq = zeros(nQuant, M); tdec0q = tdecq; tdec1q = tdecq;

ii = 1;
while(ii <= M)
    % Simulate epidemic scenarios and truncate
    Iwarn = 1; % ensure no warnings
    while Iwarn
        [Iday{ii}, Lam{ii}, Rtrue{ii}, tday{ii}, Iwarn] = epiSimDie(length(tday0),....
            1, ts, Rs, scenNo, distvals);
    end
    % Truncated observation period based on removed 0s
    nday(ii) = length(tday{ii});
   
    % Main code for elimination probabilities and declarations
    [tdecs(ii, :), tdec0(ii), didEnd(ii), ~, ~, ~, win(ii), ~, ~, ~,...
        tdec0q(:, ii), tdecq(:, ii), tdec1q(:, ii)] = getSingleEpidFn(Iday{ii},...
        Lam{ii}, nday(ii), priors, distvals, Rtrue{ii}, k, qs, mu);
    
    % Ensure epidemic ends
    if didEnd(ii)        
        % Completed iteration
        disp(['Completed: ' num2str(ii) ' of ' num2str(M)]);
        ii = ii + 1;
    end
end

% Summarise declaration time quantiles across runs
tqrun = quantile(tdecq', [0.025, 0.5, 0.975]);
tq0run = quantile(tdec0q', [0.025, 0.5, 0.975]);
tq1run = quantile(tdec1q', [0.025, 0.5, 0.975]);

% Difference between estimators
tdecEstDiff = mean(tdecs(:,1) - tdecs(:, 2));
disp(['Diff between estimators: ' num2str(tdecEstDiff)]);

% Behaviour of estimate around true tdec
err = tdecs(:, 1) - tdec0';
Pearly = length(find(err < 0))/M; Pontime = length(find(err == 0))/M;
Plate = 1 - Pontime - Pearly; Pvals = [Pearly Pontime Plate];
disp(['Early, on-time and late fractions: ' num2str(Pvals)]);

%% Plotting and saving

% Confidence intervals across runs for plotting
e1 = tqrun(2,:) - tqrun(1, :); e1 = e1';
e2 = tqrun(3, :) - tqrun(2, :); e2 = e2';
e01 = tq0run(2,:) - tq0run(1, :); e01 = e01';
e02 = tq0run(3, :) - tq0run(2, :); e02 = e02';

% Quantiles of tdec true and estimated
figure;
plotCI(100*qs, tqrun(2,:)', e1, e2, 'c');
hold on;
plotCI(100*qs, tq0run(2,:)', e01, e02, 'g');
grid off; box off; hold off;
xlabel('$\mu \%$');
ylabel('$t_{\mu}$');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['tmu_' scenChoice '_' num2str(k) '_' num2str(M) '_' num2str(distvals.type)], 'fig');
    cd(thisDir);
end

% Histograms of estimated 95 declarations
figure;
ax(1) = subplot(2, 1, 1);
histogram(tdec0, 'Normalization', 'probability');
xlabel('true $t_{95}$');
ylabel('probability');
ax(2) = subplot(2, 1, 2);
histogram(tdecs(:, 1), 'Normalization', 'probability');
xlabel('estimated $t_{95}$');
ylabel('probability');
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['t95_' scenChoice '_' num2str(k) '_' num2str(M) '_' num2str(distvals.type)], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['TdecQ_' scenChoice '_' num2str(k) '_' num2str(M) '_' num2str(distvals.type) '.mat']);
    cd(thisDir);
end