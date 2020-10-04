% Run many epidemics and get declarations across all scenarios
clearvars; clc;
close all; tic;

% Assumptions and notes
% - does all distributions
% - removed quantiles and does all distributions at once
% - only uses disease distributions
% - R trajectory is fixed but various I trajectories
% - compare true tdec with estimates
% - fix best k, prior and examine different generation times

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
saveTrue = 1; testFig = 0;
% Folder for saving
saveFol = 'batch disease'; thisDir = cd;

%% Sim setup

% Choose a generation time distribution (all gamma)
distNos = 1:4;
% Define all SI/generation time distributions
distNam = {'Marburg', 'MERS', 'Measles', 'COVID-19'};

% Confidence level for declaration
mu = 0.95;
% Define all dying epidemic scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};

% Window size and runs
k = 100; M = 1000;
% Set priors on R estimates, E[R] = ab
priors.a = 1;  priors.b = 5;

% Starting times (days)
tday0 = 1:501; nday0 = length(tday0);

%% Loop all scenarios

% All possible scenarios at a fixed k and prior
scenExplore = [1 2 4]; nScen = length(scenExplore);
tsim = zeros(1, nScen);

for iii = 1:length(distNos)
    
    % Specific distribution considered
    distNo = distNos(iii);
    distChoice = distNam{distNo};
    disp(['SI distribution: ' distNam{distNo}]);
    
    for jj = 1:nScen
        
        % Current scenario
        scenChoice = scenNam{scenExplore(jj)};
        
        % Define gamma SI distribution and R scenario
        [distvals.pm, Rs, ts, distvals.omega] = setupSIandRGamma(distNo, scenExplore(jj));
        tcheck = max(ts);
        
        % Store variables from each epidemic
        Iday = cell(1, M); Lam = Iday; Rtrue = Iday; tday = Iday;
        nday = zeros(1, M); tdec0 = nday; didEnd = nday; tdecs = zeros(M, 2);
        
        ii = 1;
        while(ii <= M)
            % Simulate epidemic scenarios and truncate
            Iwarn = 1; % ensure no warnings
            while Iwarn
                % Use epiSimDie or 2 depending on how want to choose epidemics
                [Iday{ii}, Lam{ii}, Rtrue{ii}, tday{ii}, Iwarn] = epiSimDieDisease(length(tday0),....
                    ts, Rs, scenExplore(jj), distvals);
            end
            % Truncated observation period based on removed 0s
            nday(ii) = length(tday{ii});
            
            % Main code for elimination probabilities and declarations
            [tdecs(ii, :), tdec0(ii), didEnd(ii), ~, ~, ~, ~, ~, ~, ~] =...
                getSingleEpidFnDiseaseNoQ(Iday{ii}, Lam{ii}, nday(ii), priors, distvals, Rtrue{ii}, k, mu);
            
            % Ensure epidemic ends
            if didEnd(ii)
                % Completed iteration
                disp(['Completed: ' num2str(ii) ' of ' num2str(M)]);
                ii = ii + 1;
            end
        end
        
        % Difference between estimators
        tdecEstDiff = mean(tdecs(:,1) - tdecs(:, 2));
        disp(['Diff between estimators: ' num2str(tdecEstDiff)]);
        
        % Behaviour of estimate around true tdec
        err = tdecs(:, 1) - tdec0';
        Pearly = length(find(err < 0))/M; Pontime = length(find(err == 0))/M;
        Plate = 1 - Pontime - Pearly; Pvals = [Pearly Pontime Plate];
        disp(['Early, on-time and late fractions: ' num2str(Pvals)]);
        
        % Histograms of estimated 95 declarations
        figure;
        ax(1) = subplot(2, 1, 1);
        histogram(tdec0, 'Normalization', 'probability');
        xlabel('$t^*_{95}$');
        ylabel('probability');
        ax(2) = subplot(2, 1, 2);
        histogram(tdecs(:, 1), 'Normalization', 'probability');
        xlabel('$t_{95}$');
        ylabel('probability');
        linkaxes(ax, 'xy');
        if saveTrue
            cd(saveFol);
            saveas(gcf, ['t95_' scenChoice '_' num2str(k) '_' num2str(M) '_' distChoice], 'fig');
            cd(thisDir);
        end
        
        % Histograms of error
        figure; hold on;
        h = histogram(err, 'Normalization', 'probability');
        h.FaceAlpha = 0.15; h.EdgeColor = grey2;
        xlabel('$\delta t_{95}$', 'FontSize', 18);
        % Min and max for plotting limits
        tmin = min(err); tmax = max(err);
        tmin = tmin - 1; tmax = tmax + 1;
        xlim([tmin tmax]);
        box off; grid off;
        if saveTrue
            cd(saveFol);
            saveas(gcf, ['del_t95_' scenChoice '_' num2str(k) '_' num2str(M) '_' distChoice], 'fig');
            cd(thisDir);
        end
        
        if jj == nScen
            % Generation time distribution
            tdist = 1:100;
            Pomega = gammaDistrTypes(max(tdist), distvals);
            
            figure;
            plot(tdist, Pomega, 'LineWidth', 2);
            xlabel('$u$ (days)');
            ylabel('$w_u$');
            grid off; box off;
            if saveTrue
                cd(saveFol);
                saveas(gcf, ['SI_' distChoice '_' num2str(k) '_' num2str(M) '_' distChoice], 'fig');
                cd(thisDir);
            end
        end
        
        %close all;
        
        % Timing and data saving
        tsim(jj) = toc/60;
        disp(['Run time = ' num2str(tsim(jj))]);
        if saveTrue
            cd(saveFol);
            % All Rtrue same here
            Rtrue = Rtrue{1};
            save(['TdecQ_' scenChoice '_' num2str(k) '_' num2str(M) '_' distChoice '.mat']);
            cd(thisDir);
        end
        
        % Progress update
        clc; disp(['Completed: ' num2str(jj) ' of ' num2str(nScen)]);
        disp('******************************************************');
    end
end