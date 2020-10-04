% Assess the elimination probability of a single epidemic with input data
clearvars; clc;
close all; tic;

% Assumptions and notes
% - shows how importations can change tdec, tst fixed to last Itot
% - declarations are relative to tst; when last case was seen
% - computes estimated probability with NB uncertainty (z) or mean (z1)
% - uses input epidemics and computes from first consistent 0
% - if no 0 artifially include some at end

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving and loading
loadFol = 'mers data'; thisDir = cd;
saveFol = 'mers results';

% Show whole epi-curve
wholeEp = 0;

% Confidence level for declaration
mu = 0.95;
disp(['Confidence = ' num2str(mu)]);

%% Input data from EpiEstim or other package

% Load key data from other packages
cd(loadFol);

% MERS incidence curve local and imported
Iloc = csvread("Iloc.csv", 1,1); Iloc = Iloc';
Iimp = csvread("Iimp.csv", 1,1); Iimp = Iimp';

% Artificially increase imports in non-zeros part
extraimp = 10; nz = 100; % nz known from R code
Iimp(1:end-nz) = Iimp(1:end-nz) + poissrnd(extraimp, [1 length(Iimp(1:end-nz))]);

% Total cases and fraction imported
Itot = Iloc + Iimp; fimp = sum(Iimp)/sum(Itot);
disp(['Fraction of cases imported = ' num2str(fimp)]);

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
tplt = tday(2:end); Iplt = Iloc(2:end); Iplttot = Itot(2:end);

cd(thisDir);

% Set priors on R estimates, E[R] = ab
priors.a = 1; priors.b = 5;

% Window size (k+1 with current value)
k = 100; win = k + 1;
disp(['Window length: ' num2str(win)]);

% Basic estimation with Iloc vs Iimp on long window
[pred, predInt, prob, R, RInt, ~] = getNegBinPriors(k, nday, Iloc, Lmers, priors);
[pred2, predInt2, prob2, R2, RInt2, ~] = getNegBinPriors(k, nday, Itot, Lmers, priors);

% Synchronise epidemic end based on total
tst = find(Itot, 1, 'last');
% If want to instead show whole epi-curve
if wholeEp
    % Now end time is absolute
    tst = 1;
end

% Times to interrogate for z
tindex = tst+1:length(Itot); 
lenind = length(tindex);

%% End-of-epidemic with no account for local/imports

% Elimination probabilities and R
Rtot = cell(1, lenind); zseqtot = Rtot; 
ztot = zeros(1, lenind); z1tot = ztot; 

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered is all cases
    Itemptot = Itot(1:tindex(i));
    
    % Estimates of elimination probability
    [ztot(i), zseqtot{i}, z1tot(i), ~, Rtot{i}] = getProbEliminSI(k, Itemptot, Itemptot, priors, genmers);
end

% Find points of declaration
tdecstot = [find(ztot > mu, 1, 'first'), find(z1tot > mu, 1, 'first')];
disp(['total declared at [z z1]: ' num2str(tdecstot)]);

% Check if epidemic actually ended
if any(isempty(tdecstot))
    didEnd = 0;
    disp('Epidemic did not end');
else
    didEnd = 1;
end

%% End-of-epidemic with account for local/imports

% Elimination probabilities and R
Rloc = cell(1, lenind); zseqloc = Rtot; 
zloc = zeros(1, lenind); z1loc = ztot; 

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered is all cases and just local
    Itemp = Iloc(1:tindex(i));
    Itemptot = Itot(1:tindex(i));
    
    % Estimates of elimination probability
    [zloc(i), zseqloc{i}, z1loc(i), ~, Rloc{i}] = getProbEliminSI(k, Itemp, Itemptot, priors, genmers);
end

% Find points of declaration
tdecsloc = [find(zloc > mu, 1, 'first'), find(z1loc > mu, 1, 'first')];
disp(['Local declared at [z z1]: ' num2str(tdecsloc)]);

% Check if epidemic actually ended
if any(isempty(tdecsloc))
    didEnd = 0;
    disp('Epidemic did not end');
else
    didEnd = 1;
end


%% Visualise local vs total

% Domains of distributions
tdist = 1:100; Rdist = linspace(0, 30, 10000);
% Generation time distribution
Pomega = genmers(tdist);
% Priors on R
prior_R = gampdf(Rdist, priors.a, priors.b);

% Local vs total incidence
figure;
stairs(tday, Itot, 'c', 'LineWidth', 2);
hold on;
stairs(tday, Iloc, 'Color', grey2, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)');
ylabel('$I_s$');
legend('$I_{total}$', '$I_{local}$', 'Location', 'Best');

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
    plot(zloc, 'Color', 'r', 'LineWidth', 2);
    hold on;
    plot(ztot, 'Color', 'k', 'LineWidth', 2);
    plot(1:length(zloc), mu*ones(size(zloc)), 'c--', 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$\Delta s$ (days)');
    ylabel('$z_s$');
    xlim([1 length(zloc)]);
    legend('$z_{local}$', '$z_{total}$', 'Location', 'Best');
    ylim([0 1.1]);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['mers_' num2str(nday) '_' num2str(100*extarimp) '.mat']);
    cd(thisDir);
end





