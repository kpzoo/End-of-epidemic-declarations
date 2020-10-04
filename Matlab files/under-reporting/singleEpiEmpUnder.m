% Elimination probability of a single epidemic with input data and underreporting
clearvars; clc;
close all; tic;

% Assumptions and notes
% - tst set to last time im fully sampled case
% - shows how constant underreporting can change tdec
% - declarations are relative to tst - last case time
% - computes estimated probability with NB uncertainty (z) or mean (z1)
% - uses input epidemics and computes from first consistent 0

% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving and loading
loadFol = 'epi data/sars'; thisDir = cd;
saveFol = 'epi results/sars';

% Show whole epi-curve
wholeEp = 0;

% Confidence level for declaration
mu = 0.95;
disp(['Confidence = ' num2str(mu)]);

%% Input data from EpiEstim or other package

% Load key data from other packages
cd(loadFol);

% Flu1918 incidence curve (only local cases)
Iloc = csvread("Iday.csv", 1,1); Iloc = Iloc';
% Serial interval and total infectiousness
Lday = csvread("Lday.csv", 1,1); Lday = Lday';
Pomega = csvread("sidist.csv", 1,1); Pomega = Pomega';

% Zeros and sample fraction
rho = 0.3; nz = 100; % nz known from R code
Isamp = zeros(size(Iloc));

% Synchronise epidemic end based on total
tst = find(Iloc, 1, 'last');
% If want to instead show whole epi-curve
if wholeEp
    % Now end time is absolute
    tst = 1;
end

% Only sample before last incidence of Iloc to align tst
%Isamp(1:tst-1) = binornd(Iloc(1:tst-1), rho); Isamp(tst) = Iloc(tst);
Isamp(1:tst) = binornd(Iloc(1:tst), rho); 
disp(['Underreporting = ' num2str(rho)]);

% Confirm fraction sampled
fimp = sum(Isamp)/sum(Iloc); 
disp(['Fraction of cases sampled = ' num2str(fimp)]);

% Days to consider
nday = length(Iloc); tday = 1:nday;
% Plot from second time onwards
tplt = tday(2:end); Iplt = Iloc(2:end); Ipltsamp = Isamp(2:end);

cd(thisDir);

% Set priors on R estimates, E[R] = ab
priors.a = 1; priors.b = 5;

% Window size (k+1 with current value)
k = 99; win = k + 1;
disp(['Window length: ' num2str(win)]);

% Undersampled total infectiousness
Lsamp = zeros(1, nday);
for i = 2:nday
    % Relevant part of SI: Pomega(1:i-1))
    Lsamp(i) = sum(Isamp(i-1:-1:1).*Pomega(1:i-1));
end

% Basic estimation with Iloc vs Iimp on long window
[pred, predInt, prob, R, RInt, ~] = getNegBinPriors(k, nday, Iloc, Lday, priors);
[pred2, predInt2, prob2, R2, RInt2, ~] = getNegBinPriors(k, nday, Isamp, Lsamp, priors);

% Times to interrogate for z
tindex = tst+1:length(Iloc); 
lenind = length(tindex);


%% End-of-epidemic without underreporting

% Elimination probabilities and R
Rloc = cell(1, lenind); zseqloc = Rloc; 
zloc = zeros(1, lenind); z1loc = zloc; 

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered is all cases and just local
    Itemp = Iloc(1:tindex(i));
    
    % Estimates of elimination probability
    [zloc(i), zseqloc{i}, z1loc(i), ~, Rloc{i}] = getProbEliminSI(k, Itemp, Itemp, priors, Pomega);
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

%% End-of-epidemic with underreporting

% Elimination probabilities and R
Rsamp = cell(1, lenind); zseqsamp = Rsamp; 
zsamp = zeros(1, lenind); z1samp = zsamp; 

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered is all cases
    Itempsamp = Isamp(1:tindex(i));
    
    % Estimates of elimination probability
    [zsamp(i), zseqsamp{i}, z1samp(i), ~, Rsamp{i}] = getProbEliminSI(k, Itempsamp, Itempsamp, priors, Pomega);
end

% Find points of declaration
tdecssamp = [find(zsamp > mu, 1, 'first'), find(z1samp > mu, 1, 'first')];
disp(['total declared at [z z1]: ' num2str(tdecssamp)]);

% Check if epidemic actually ended
if any(isempty(tdecssamp))
    didEnd = 0;
    disp('Epidemic did not end');
else
    didEnd = 1;
end


%% Visualise local vs total

% Domains of distributions
tdist = 1:100; Rdist = linspace(0, 30, 10000);
% Generation time distribution
Pomega = Pomega(tdist);
% Priors on R
prior_R = gampdf(Rdist, priors.a, priors.b);

% Local vs total incidence
figure;
stairs(tday, Iloc, 'c', 'LineWidth', 2);
hold on;
stairs(tday, Isamp, 'Color', grey2, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)');
ylabel('$I_s$');
legend('$I_{total}$', '$I_{sampled}$', 'Location', 'Best');

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
    plot(zsamp, 'Color', 'k', 'LineWidth', 2);
    plot(1:length(zloc), mu*ones(size(zloc)), 'c--', 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$\Delta s$ (days)');
    ylabel('$z_s$');
    xlim([1 length(zloc)]);
    legend('$z_{total}$', '$z_{sampled}$', 'Location', 'Best');
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





