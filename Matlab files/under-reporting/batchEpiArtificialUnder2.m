% Batch the sampling of an empirical epidemic with input data
clearvars; clc;
close all; tic;

% Assumptions and notes
% - always include rho = 1 for comparison
% - shows how binomial sampling at rho changes tdec
% - relative to tst; when last case was seen for each curve
% - computes estimated probability with NB uncertainty (z) or mean (z1)
% - uses input epidemics and computes from first consistent 0=

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Save data 
saveTrue = 1; thisDir = cd;
% Folder for saving and loading
loadFol = 'epi data/sars'; 
saveFol = 'epi results/sars';

% Confidence level for declaration and quantiles
mu = 0.95; aq = 0.025;
disp(['Confidence = ' num2str(mu)]);

%% Input data from EpiEstim or other package

% Load key data from other packages
cd(loadFol);

% Incidence curve (only local cases)
Iloc = csvread("Iday.csv", 1,1); Iloc = Iloc';
% Synchronise epidemic end to last 0 in Iloc
tst = find(Iloc, 1, 'last'); totLoc = sum(Iloc);

% Serial interval and total infectiousness
Lday = csvread("Lday.csv", 1,1); Lday = Lday';
Pomega = csvread("sidist.csv", 1,1); Pomega = Pomega';

% Days to consider and zeros from R
nday = length(Iloc); tday = 1:nday; nz = 100;
% Plot from second time onwards
tplt = tday(2:end); Iplt = Iloc(2:end);

cd(thisDir);

% Set priors on R estimates, E[R] = ab
priors.a = 1; priors.b = 5;

% Window size (k+1 with current value)
k = 100; win = k + 1;
disp(['Window length: ' num2str(win)]);

% Sampling fractions to consider
rho = linspace(0.1, 1, 10); lenE = length(rho);

% Artifically sample epi-curve and reconstruct Lam
nRep = 5000; IsampE = cell(1, lenE); IsampTemp = repmat(Iloc, [nRep 1]);
% Fractions of and total cases from false additions over replicates
fsampE = IsampE; fsampTemp = IsampTemp; ItotE = IsampE;


% Issue - decide where to start counting from e.g. first of 0s in Iloc?
% or first of 0s in any sequence? Chose latter here
for i = 1:lenE
    % Binomial samples with mean rho taken
    for j = 1:tst
        IsampTemp(:, j) = binornd(Iloc(j), rho(i), [nRep 1]);
    end
    % Merge with last value from original epi-curve
    %IsampTemp(:, tst) = Iloc(tst);
    % Assigned to cell - nRep replicate samples
    IsampE{i} = IsampTemp;
    % Fractions of cases and mean fraction
    fsampTemp = sum(IsampTemp, 2)/totLoc; 
    fsampTemp = fsampTemp'; fsampE{i} = fsampTemp;
end
% Mean fractions of imports
fracSampE = cellfun(@mean, fsampE);

% Find earliest time for last zero and align with that
tstmin = zeros(lenE, nRep);
for i = 1:lenE
    for j = 1:nRep
        tstmin(i, j) = find(IsampE{i}(j, :), 1, 'last');
    end
end
tstmin = min(min(tstmin));

%% End-of-epidemic with and without underreporting

% Elimination probabilities 
zsamp = cell(lenE, nRep); 
% Declaration times in relative time from last total case
tdecsamp = cell(1, lenE); 
% Check all epidemics ended and time references
didEndSamp = zeros(1, lenE); tstsamp = tdecsamp;

% For every boost in imports
for i = 1:lenE
    % Get sampled incidence for each replicate (local fixed)
    Isamp = IsampE{i};
    
    % Check if each replicate epidemic ended
    didEndSampTemp = zeros(1, nRep); 
    % Check on starting time reference and other declaration
    tdecCheckSamp = didEndSampTemp; tstsampj = didEndSampTemp;
    % Declaration times
    tdecsampTemp = didEndSampTemp; 
    
    % For every replicate at the boosted mean extraImp
    for j = 1:nRep
        % Compute elimination probability with underreporting
        [zsamp{i, j}, ~, ~, didEndSampTemp, ~, tdecsampTemp(j), tdecCheckSamp(j), tstsampj(j)]...
            = endIncid2(Isamp(j,:), Isamp(j,:), k, priors, Pomega, mu, tstmin);
    end
    
    % Store declaration times abd start times
    tdecsamp{i} = tdecsampTemp; tstsamp{i} = tstsampj; 
    
    % Check epidemics ended in all replicates
    didEndSamp(i) = all(didEndSampTemp); 
    % Check all declaration times match
    if ~all(tdecsampTemp == tdecCheckSamp)
        warning('Declaration times do not match');
    end
    
    disp(['Completed ' num2str(i) ' of ' num2str(lenE)]);
end

% Mean starting 0 time of all epidemics 
tstavg = cellfun(@mean, tstsamp);
tstsamp = cell2mat(tstsamp');

% Ensure all epidemics ended
if ~all(didEndSamp) 
    error('All epidemics did not end');
else
    % All checks cleared so variables removed
   clear didEndSampTemp tdecCheckSamp tstsampj;
end

% Convert cells and get quantiles
tdecsamp = cell2mat(tdecsamp'); 
tdecsampq = quantile(tdecsamp', [aq, 0.5, 1-aq]);

% Combine tst and tdec to get absolute time
tcomb = tstsamp + tdecsamp;

% Difference between declaration times (last is rho = 1)
tdiff = tdecsamp - tdecsamp(end, :);
tdiffq = quantile(tdiff', [aq, 0.5, 1-aq]);

% Difference in absolute time
tdcomb = tcomb - tcomb(end, :);
tdcombq = quantile(tdcomb', [aq, 0.5, 1-aq]);

%% Visualise local vs total

% Pick 4 parts of the range of rho
id = [1 3 6 lenE-1 lenE];

% Extract zsamp and get stats in ids
nid = length(id); zsampp = zsamp(id, :); 
% Combine all runs for each id
zsampid = cell(1, nid); 
for i = 1:nid
    % Collect raw values
    zsampid{i} = cell2mat(zsampp(i, :)');
    % Take quantiles
    zsampid{i} = quantile(zsampid{i}, [aq, 0.5, 1-aq]);
end


% Domains of distributions
tdist = 1:50; Rdist = linspace(0, 30, 10000);
% Generation time distribution
Pomegaplt = Pomega(tdist);
% Priors on R
prior_R = gampdf(Rdist, priors.a, priors.b);

% Generation time and R prior
figure;
subplot(2, 1, 1);
plot(tdist, Pomegaplt, 'LineWidth', 2);
xlabel('$s$ (days)');
ylabel('$w_s$');
grid off; box off;
subplot(2, 1, 2);
plot(Rdist, prior_R, 'LineWidth', 2);
xlabel('$R_s$');
ylabel('$P(R_s)$');
grid off; box off;

% Declaration times versus extra imports at ids 
figure; tzplt = tstmin+1:nday; %tzplt = 1:nz+1;
fracplt = round(fracSampE, 2, 'significant');
for i = 1:nid-1
    subplot(ceil((nid-1)/2), 2, i);
    zsampplt = zsampid{i}; zlocplt = zsampid{end};
    plot(tzplt, mu*ones(size(tzplt)), 'k--', 'LineWidth', 2);
    hold on;
    plotCIRaw(tzplt', zsampplt(2,:)', zsampplt(1,:)', zsampplt(3,:)', 'c');
    plotCIRaw(tzplt', zlocplt(2,:)', zlocplt(1,:)', zlocplt(3,:)', 'r');
    hold off; grid off; box off;
    xlim([tstmin+1 tstmin+100]); ylabel(['$z_s | \rho = $' num2str(fracplt(id(i)))]);
    if any(i == [nid-2 nid-1])
        xlabel('$s$ (days)');
    end
end
if saveTrue
    cd(saveFol);
    saveas(gcf, ['elimRepPlts_' num2str(k) '_' num2str(nRep) '_' num2str(lenE)], 'fig');
    cd(thisDir);
end

% Difference in declaration times in absolute terms (combined with tsts)
figure;
plotCI2(fracSampE', tdcombq(2,:)', tdcombq(2,:)'-tdcombq(1,:)', tdcombq(3,:)'-tdcombq(2,:)', 'b', 1);
grid off; box off;
ylabel('$\Delta t_{95}$ (days)');
xlabel('$\rho$ (sampled)');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['decAbsRep_' num2str(k) '_' num2str(nRep) '_' num2str(lenE)], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['sarsRep_' num2str(nday) '_' num2str(lenE) '_' num2str(nRep) '.mat']);
    cd(thisDir);
end





