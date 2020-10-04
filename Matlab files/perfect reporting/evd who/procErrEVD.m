% Replot batch of histograms for comparison
clearvars; clc;
close all; tic;

% Assumptions and notes
% - assumes data generated from allBatchDiseaseEVD.m
% - process error in declaration times and raw declaration times

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save data and test figs
saveTrue = 0; thisDir = cd;
% Folder for loading
loadFol = thisDir; saveTrue = 0;
% Choose disease
distChoice = 5; nDist = 1;

% WHO declaration time
twho = 42 + 5;

% Possible scenario names
scenNam = {'control', 'recovery' 'boom-bust'};
% Size of figure in panels
nScen = 1; 

% Possible SI names
distNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'Ebola virus'};
% Name of disease considered
dname = distNam{distChoice};

% All declaration errors and actual times times
edec = cell(nScen, nDist); t0 = edec; tx = edec;
% Choices of k for checking
kval = zeros(nScen, nDist);
% True R and serial interval
wdist = cell(1, nDist); Rtrue = cell(1, nScen);

% Process files into a combined histogram
cd(loadFol); 
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
            
            % Only want main estimate (not exp approx)
            tdec = tdecs(:, 1); tdec = tdec';
            % Raw error in estimates 
            edec{i, j} = tdec - tdec0; 
            % Absolute declaration times
            t0{i, j} = tdec0; tx{i, j} = tdec;
            
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

% Get length of runs
m = length(tdec); clearvars tdecs tdec0 tdec;
% Get k and number of repetitions
nReps = length(edec{1, 1}); k = unique(kval);

% Histogram construction for difference in declaration
figure;
for j = 1:nScen
    ax(j) = subplot(1, nScen, j);
    h = histogram(edec{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    hold on
    %ksdensity(edec{j, i}, 'Kernel', 'normal')
    xlabel('$\delta t_{95}$', 'FontSize', 18);
    ylabel(scenNam{j}, 'FontSize', 18);
    box off; grid off; hold off;
    if j == round(nScen/2)
        title(['Difference in $t_{95}$ for ' dname], 'FontSize', 18);
    end
end
if saveTrue
    saveas(gcf, ['delta_' dname '_' num2str(k) '_' num2str(m)], 'fig');
    close;
end

% Histogram construction for true declaration
figure;
for j = 1:nScen
    ax(j) = subplot(1, nScen, j);
    h = histogram(t0{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    hold on; h = gca;
    ksdensity(t0{j}, 'Kernel', 'normal')
    plot([twho twho], h.YLim, 'k--', 'LineWidth', 2);
    xlabel('$\Delta t_{95}^*$', 'FontSize', 18);
    ylabel(scenNam{j}, 'FontSize', 18);
    box off; grid off; hold off;
    if j == round(nScen/2)
        title(['True $t_{95}^*$ for ' dname], 'FontSize', 18);
    end
end
if saveTrue
    saveas(gcf, ['tz0_' dname '_' num2str(k) '_' num2str(m)], 'fig');
    close;
end

% Histogram construction for estimated declaration
figure;
for j = 1:nScen
    ax(j) = subplot(1, nScen, j);
    h = histogram(tx{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    hold on; h = gca;
    ksdensity(tx{j}, 'Kernel', 'normal')
    plot([twho twho], h.YLim, 'k--', 'LineWidth', 2);
    xlabel('$\Delta t_{95}$', 'FontSize', 18);
    ylabel(scenNam{j}, 'FontSize', 18);
    box off; grid off; hold off;
    if j == round(nScen/2)
        title(['Estimated $t_{95}$ for ' dname], 'FontSize', 18);
    end
end
if saveTrue
    saveas(gcf, ['tz_' dname '_' num2str(k) '_' num2str(m)], 'fig');
    close;
end

% Colour vectors for plotting
cols = {'b', grey1, 'r', grey2};

% Plot true R values and SI
figure; 
subplot(2, 1, 1); hold on;
t = 1:length(Rtrue{1});
for i = 1:nScen
    plot(t, Rtrue{i}, 'Color', cols{i}, 'LineWidth', 2);
end
box off; grid off; hold off;
xlabel('$s$ (days)', 'FontSize', 18);
ylabel('$R_s$', 'FontSize', 18);
legend(scenNam{:}, 'Location', 'best', 'FontSize', 16, 'Box', 'off');
xlim([1 length(t)]);

% Plot serial interval distributions
subplot(2, 1, 2); hold on;
u = 1:length(wdist{1});
for i = 1:nDist
    stairs(u, wdist{i}, 'Color', cols{i}, 'LineWidth', 2);
end
box off; grid off; hold off;
xlabel('$u$ (days)', 'FontSize', 18);
ylabel('$w_u$', 'FontSize', 18);
xlim([0 60]);
legend(dname, 'Location', 'best', 'FontSize', 16, 'Box', 'off');
if saveTrue
    saveas(gcf, ['SIandRdisease_' num2str(k) '_' num2str(m)], 'fig');
end


% Large figure that summarises a scenario for a disease
if nScen == 1 && nDist == 1    
    % Index of case
    j = 1;
    figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
    
    % True declaration times
    subplot(2, 2, 1);
    h = histogram(t0{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    hold on; h = gca;
    ksdensity(t0{j}, 'Kernel', 'normal')
    plot([twho twho], h.YLim, 'k--', 'LineWidth', 2);
    xlabel('True: $\Delta t_{95}^*$', 'FontSize', 18);
    box off; grid off; hold off;
    ylabel(distNam{distChoice}, 'FontSize', 18);
    
    % Estimated declaration times
    subplot(2, 2, 3);
    h = histogram(tx{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    hold on; h = gca;
    ksdensity(tx{j}, 'Kernel', 'normal')
    plot([twho twho], h.YLim, 'k--', 'LineWidth', 2);
    xlabel('Estimated: $\Delta t_{95}$', 'FontSize', 18);
    box off; grid off; hold off;
    ylabel(distNam{distChoice}, 'FontSize', 18);
    
    % Error in declaration times
    subplot(2, 2, 2);
    h = histogram(edec{j}, 'Normalization', 'probability');
    h.FaceAlpha = 0.15; h.EdgeColor = grey2;
    xlabel('Error: $\delta t_{95}$', 'FontSize', 18);
    box off; grid off; 
    ylabel(distNam{distChoice}, 'FontSize', 18);
    
    % Serial interval and R
    subplot(2, 2, 4);
    u = 1:length(wdist{j});
    stairs(u, wdist{j}, 'Color', 'b', 'LineWidth', 2);
    box off; grid off; 
    xlabel('Elapsed time, $u$ (days)', 'FontSize', 18);
    ylabel('Serial interval, $w_u$', 'FontSize', 18);
    xlim([0 60]);
    axes('Position',[.74 .28 .15 .15])
    t = 1:length(Rtrue{j});
    plot(t, Rtrue{j}, 'Color', 'r', 'LineWidth', 2);
    box off; grid off; 
    xlabel('Time, $s$ (days)', 'FontSize', 18);
    ylabel('Reprod. no., $R_s$', 'FontSize', 18);
    xlim([1 300]);
    
end


