% Assess the elimination probability of a single epidemic with input data
clearvars; clc;
close all; tic;

% Assumptions and notes
% - artifially increase the imports with Poiss(g)
% - shows how importations can change tdec
% - declarations are relative to tst; when last case was seen
% - computes estimated probability with NB uncertainty (z) or mean (z1)
% - uses input epidemics and computes from first consistent 0
% - if no 0 artifially include some at end

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Save data and test figs
saveTrue = 0; 
% Folder for saving and loading
loadFol = 'mers data'; thisDir = cd;
saveFol = 'mers results';

% Confidence level for declaration and quantiles
mu = 0.95; aq = 0.025;
disp(['Confidence = ' num2str(mu)]);

%% Input data from EpiEstim or other package

% Load key data from other packages
cd(loadFol);

% MERS incidence curve local and imported
Iloc = csvread("Iloc.csv", 1,1); Iloc = Iloc';
Iimp = csvread("Iimp.csv", 1,1); Iimp = Iimp';
% Original fraction of imported
totLoc = sum(Iloc); totImp = sum(Iimp);
fracImp = totImp/(totLoc + totImp);
disp(['Original import fraction: ' num2str(fracImp)]);

% MERS serial interval and total infectiousness
Lmers = csvread("Lmers.csv", 1,1); Lmers = Lmers';
genmers = csvread("genmers.csv", 1,1); genmers = genmers';

% Days to consider
if length(Iloc) == length(Iimp)
    nday = length(Iloc); tday = 1:nday;
else
    error('Input data inconsistent');
end
% Plot from second time onwards
tplt = tday(2:end); Iplt = Iloc(2:end); 

cd(thisDir);

% Set priors on R estimates, E[R] = ab
priors.a = 1; priors.b = 5;

% Window size (k+1 with current value)
k = 100; win = k + 1;
disp(['Window length: ' num2str(win)]);

% Non-zeros maintained (from R) and indices of incidence to be increased
nz = 100; idinc = 1:nday-nz; leninc = length(idinc);
% Check only zeros after idinc
if ~all(Iimp(idinc(end)+1:end) == 0) || ~all(Iloc(idinc(end)+1:end) == 0)
    error('Did not properly account for zeros');
end

% Additional imports (mean of Poisson) - local cases fixed
extraImp = linspace(0, 2, 10); lenE = length(extraImp);

% Artificially increase imports in non-zeros part
nRep = 1000; IimpE = cell(1, lenE); IimpTemp = zeros(nRep, nday);
% Fractions of and total cases from false additions over replicates
fimpE = IimpE; fimpTemp = IimpTemp; ItotE = IimpE;

for i = 1:lenE
    % Poisson samples with mean extraImp added
    randSamp = poissrnd(extraImp(i), [nRep, leninc]);
    IimpTemp(:, idinc) = randSamp +  Iimp(idinc);
    % Assigned to cell - nRep replicate samples
    IimpE{i} = IimpTemp;
    % Fractions of cases and mean fraction
    fimpTemp = sum(IimpTemp, 2)./(sum(IimpTemp, 2) + totLoc); 
    fimpTemp = fimpTemp'; fimpE{i} = fimpTemp;
    % Total number of cases
    ItotE{i} = IimpE{i} + Iloc;
end
% Mean fractions of imports
fracImpE = cellfun(@mean, fimpE);


%% End-of-epidemic ignoring vs accounting for local/imports

% Elimination probabilities locally (loc) and ignoring case type (tot)
ztot = cell(lenE, nRep); zloc = ztot;
% Declaration times in relative time from last total case
tdectot = cell(1, lenE); tdecloc = tdectot;
% Check all epidemics ended
didEndTot = zeros(1, lenE); didEndLoc = didEndTot;
% Check time references
tstloc = tdectot; tsttot = tdectot;

% For every boost in imports
for i = 1:lenE
    % Get total incidence for each replicate (local fixed)
    Itot = ItotE{i};
    
    % Check if each replicate epidemic ended
    didEndTotTemp = zeros(1, nRep); didEndLocTemp = didEndTotTemp;
    % Check on starting time reference and other declaration
    tdecCheckTot = didEndTotTemp; tdecCheckLoc = didEndTotTemp; 
    tstlocj = didEndTotTemp; tsttotj = didEndTotTemp;
    % Declaration times
    tdectotTemp = didEndTotTemp; tdeclocTemp = didEndTotTemp;
    
    % For every replicate at the boosted mean extraImp
    for j = 1:nRep
        % Compute elimination probability ignoring local/imports
        [ztot{i, j}, ~, ~, didEndTotTemp, ~, tdectotTemp(j), tdecCheckTot(j), tsttotj(j)]...
            = endIncid(Itot(j,:), Itot(j,:), k, priors, genmers, mu);
        % Compute elimination probability accounting for local/imports
        [zloc{i, j}, ~, ~, didEndLocTemp, ~, tdeclocTemp(j), tdecCheckLoc(j), tstlocj(j)]...
            = endIncid(Iloc, Itot(j,:), k, priors, genmers, mu);
    end
    
    % Store declaration times abd start times
    tdectot{i} = tdectotTemp; tdecloc{i} = tdeclocTemp;
    tsttot{i} = tsttotj; tstloc{i} = tstlocj;
    
    % Check epidemics ended in all replicates
    didEndTot(i) = all(didEndTotTemp); didEndLoc(i) = all(didEndLocTemp);
    % Check all declaration times match
    if ~all(tdeclocTemp == tdecCheckLoc) || ~all(tdectotTemp == tdecCheckTot)
        warning('Declaration times do not match');
    end
    
    disp(['Completed ' num2str(i) ' of ' num2str(lenE)]);
end

% Check all epidemics evaluated at same start time
tsttot = cell2mat(tsttot); tsttot = unique(tsttot);
tstloc = cell2mat(tstloc); tstloc = unique(tstloc);
if tsttot ~= tstloc
    error('Issue with starting times');
end

% Ensure all epidemics ended
if ~all(didEndTot) || ~all(didEndLoc)
    error('All epidemics did not end');
else
    % All checks cleared so variables removed
   clear didEndLocTemp didEndTotTemp tdecCheckLoc tdecCheckTot tstlocj tsttotj; 
end

% Convert cells and get quantiles
tdectot = cell2mat(tdectot'); tdecloc = cell2mat(tdecloc');
tdectotq = quantile(tdectot', [aq, 0.5, 1-aq]);
tdeclocq = quantile(tdecloc', [aq, 0.5, 1-aq]);

% Difference between declaration times
tdiff = tdectot - tdecloc;
tdiffq = quantile(tdiff', [aq, 0.5, 1-aq]);

%% Visualise local vs total

% Pick 4 parts of the range of extraImp
%id = round(quantile(1:lenE, [0 1/3 2/3 1]));
id = [1 3 6 lenE];

% Extract ztot and zloc and get stats in ids
nid = length(id); ztotp = ztot(id, :); zlocp = zloc(id, :);
% Combine all runs for each id
ztotid = cell(1, nid); zlocid = ztotid;
for i = 1:nid
    % Collect raw values
    ztotid{i} = cell2mat(ztotp(i, :)');
    zlocid{i} = cell2mat(zlocp(i, :)');
    % Take quantiles
    ztotid{i} = quantile(ztotid{i}, [aq, 0.5, 1-aq]);
    zlocid{i} = quantile(zlocid{i}, [aq, 0.5, 1-aq]);
end


% Domains of distributions
tdist = 1:100; Rdist = linspace(0, 30, 10000);
% Generation time distribution
Pomega = genmers(tdist);
% Priors on R
prior_R = gampdf(Rdist, priors.a, priors.b);

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

% Declaration times versus extra imports at ids
figure; tzplt = 1:nz; fracplt = round(fracImpE, 2, 'significant');
for i = 1:nid
    subplot(ceil(nid/2), 2, i);
    ztotplt = ztotid{i}; zlocplt = zlocid{i};
    plot(tzplt, mu*ones(size(tzplt)), 'k--', 'LineWidth', 2);
    hold on;
    plotCIRaw(tzplt', ztotplt(2,:)', ztotplt(1,:)', ztotplt(3,:)', 'c');
    plotCIRaw(tzplt', zlocplt(2,:)', zlocplt(1,:)', zlocplt(3,:)', 'r');
    hold off; grid off; box off;
    xlim([0 40]); ylabel(['$z_s | f(\epsilon) = $' num2str(fracplt(id(i)))]);
    if any(i == [nid-1 nid])
        xlabel('$\Delta s$ (days)');
    end
end
if saveTrue
    cd(saveFol);
    saveas(gcf, ['elimPlts_' num2str(k) '_' num2str(nRep) '_' num2str(lenE)], 'fig');
    cd(thisDir);
end

% Difference in declaration times vs noise
figure;
%plotCI2(extraImp', tdiffq(2,:)', tdiffq(2,:)'-tdiffq(1,:)', tdiffq(3,:)'-tdiffq(2,:)', 'b', 1);
plotCI2(fracImpE', tdiffq(2,:)', tdiffq(2,:)'-tdiffq(1,:)', tdiffq(3,:)'-tdiffq(2,:)', 'b', 1);
grid off; box off;
ylabel('Differeince in $t_{95}$, $\delta t_{95}$ (days)', 'FontSize', 18);
xlabel('Fraction of cases imported, $f_{\epsilon}$', 'FontSize', 18);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['decdiff_' num2str(k) '_' num2str(nRep) '_' num2str(lenE)], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['mers_' num2str(nday) '_' num2str(lenE) '.mat']);
    cd(thisDir);
end





