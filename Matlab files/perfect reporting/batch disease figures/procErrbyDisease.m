% Replot batch of histograms for comparison
clearvars; clc;
close all; tic;

% Assumptions and notes
% - assumes data generated from allBatchDisease.m
% - process error in declaration times
% - groups by disease of interest

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
saveTrue = 0; thisDir = cd;
% Folder for loading
loadFol = './batch disease/k100data';

% Possible scenario names
scenNam = {'control', 'recovery' 'boom-bust'};
% Possible SI names
distNam = {'Marburg', 'MERS', 'Measles', 'COVID-19'};
% Size of figure in panels
nScen = length(scenNam); nDist = length(distNam);

% All declaration errors and true times
edec = cell(nScen, nDist); t0 = edec;
% Choices of k for checking and id for distributions
kval = zeros(nScen, nDist); idDist = kval;
% True R and serial interval
wdist = cell(1, nDist); Rtrue = cell(1, nScen);

% Process files into a combined histogram
cd ..; cd(loadFol);
for i = 1:nScen
    % All files of a single scenario
    files = dir(['TdecQ_' scenNam{i} '*']);
    nFiles = length(files);
    
    % Ensure correct file count
    if nFiles ~= nDist
        error('Incorrect number of data files');
    else
        % Load each file for given scenario
        for j = 1:nDist
            % Declaration time variables
            clearvars tdecs tdec0 tdec;
            load(files(j).name, 'tdecs', 'tdec0');
            
            % Which distribution is being considered
            whichDist = regexp(files(j).name, distNam);
            idDist(i, j) = find(~cellfun(@isempty, whichDist));
            
            % Only want main estimate (not exp approx)
            tdec = tdecs(:, 1); tdec = tdec';
            % Raw error in estimates and true declaration time
            edec{i, j} = tdec - tdec0; t0{i, j} = tdec0;
            
            % Serial interval only available from last nScen
            if i == nScen
                dat = load(files(j).name, 'Pomega');
                wdist{j} = dat.Pomega;
            end
            % True R available from last nDist
            if j == nDist
                dat = load(files(j).name, 'Rtrue');
                Rtrue{i} = dat.Rtrue;
            end
            % Choice of k
            dat = load(files(j).name, 'k');
            kval(i, j) = dat.k;
        end
    end
end
cd (thisDir);

% Statistics of runs
emean = cellfun(@mean, edec);
emed = cellfun(@(x) quantile(x, 0.5), edec);
elow = cellfun(@(x) quantile(x, 0.025), edec);
ehigh = cellfun(@(x) quantile(x, 1-0.025), edec);
eabsmean = cellfun(@(x) mean(abs(x)), edec);

% Statistics of true declaration times
t0min = cellfun(@min, t0);
t0max = cellfun(@max, t0);
t0range = cellfun(@range, t0);

% Check distribution order consistent
consistDist = zeros(1, nScen-1);
for i = 2:nScen
    if all(idDist(1, :) == idDist(i, :))
        consistDist(i-1) = 1;
    end
end
if sum(consistDist) == 2
    disp('Distribution choices consistent');
    distNam = distNam(idDist(1, :));
else
    error('Problem with distribution choices');
end


% Get length of runs
m = length(tdec); clearvars tdecs tdec0 tdec;
% Get k and number of repetitions
nReps = length(edec{1, 1}); k = unique(kval);

% % Plot error histograms for each disease
% for i = 1:nDist
%     % Name of disease considered
%     dname = distNam{i};
%     % Histogram construction
%     figure;
%     for j = 1:nScen
%         ax(j) = subplot(1, nScen, j);
%         h = histogram(edec{j, i}, 'Normalization', 'probability');
%         h.FaceAlpha = 0.15; h.EdgeColor = grey2;
%         hold on
%         %ksdensity(edec{j, i}, 'Kernel', 'normal')
%         xlabel('$\delta t_{95}$', 'FontSize', 18);
%         box off; grid off; hold off;
%         if j == round(nDist/2)
%             title(distNam{i});
%         end
%     end
%     linkaxes(ax, 'xy');
% end

% Colour vectors for plotting
cols = {'b', grey1, 'r', grey2};

% Plot true R values
figure;
subplot(2, 1, 1); hold on;
t = 1:length(Rtrue{1});
for i = 1:nScen
    plot(t, Rtrue{i}, 'Color', cols{i}, 'LineWidth', 2);
end
box off; grid off; hold off;
xlabel('Time, $s$ (days)', 'FontSize', 18);
ylabel('Reprod. number, $R_s$', 'FontSize', 18);
legend(scenNam{:}, 'Location', 'best', 'FontSize', 16, 'Box', 'off');
xlim([1 length(t)]);

% Plot serial interval distributions
subplot(2, 1, 2); hold on;
u = 1:length(wdist{1});
for i = 1:nDist
    stairs(u, wdist{i}, 'Color', cols{i}, 'LineWidth', 2);
end
box off; grid off; hold off;
xlabel('Elapsed time, $u$ (days)', 'FontSize', 18);
ylabel('Serial interval, $w_u$', 'FontSize', 18);
xlim([0 40]);
legend(distNam{:}, 'Location', 'best', 'FontSize', 16, 'Box', 'off');
if saveTrue
    saveas(gcf, ['SIandRdisease_' num2str(k) '_' num2str(m)], 'fig');
    close;
end

% Plot error histograms for each scenario
for i = 1:nScen
    % Name of disease considered
    sname = scenNam{i};
    % Histogram construction
    figure; hold on;
    for j = 1:nDist
        ax(j) = subplot(1, nDist, j);
        h = histogram(edec{i, j}, 'Normalization', 'probability');
        h.FaceAlpha = 0.15; h.EdgeColor = grey2;
        hold on
        %ksdensity(edec{j, i}, 'Kernel', 'normal')
        xlabel('$\delta t_{95}$', 'FontSize', 18);
        box off; grid off; hold off;
        ylabel(distNam{j}, 'FontSize', 18);
        if j == round(nDist/2)
            title(['Difference in $t_{95}$ for ' sname ' scenario'], 'FontSize', 18);
        end
    end
    linkaxes(ax, 'xy');
    if saveTrue
        saveas(gcf, [sname '_' num2str(k) '_' num2str(m)], 'fig');
        close;
    end
end

