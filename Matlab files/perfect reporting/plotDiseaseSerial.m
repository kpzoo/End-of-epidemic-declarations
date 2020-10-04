% Plot all serial disease distributions
type = 1:5; ntype = length(type);
shapepm = zeros(1, ntype);
scalepm = shapepm; sumgam = shapepm;

% Possible SI names
distNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'EVD'};

% Time to consider
t = 1:100; y = cell(1, ntype);

for i = type
    % Shape and mean hyerparameters of generation time
    switch(i)
        case 1
            % Marburg distribution from van Kerkhove 2015
            omega = 9; pm = (omega^2)/(5.4^2); % actually Marburg
        case 2
            % MERS distribution from Cauchemez 2016
           omega = 6.8; pm = (omega^2)/(4.1^2);
        case 3
            % Measles distribution from Cori 2013
            omega = 14.9; pm = (omega^2)/(3.9^2);
        case 4
            % COVID-19 distribution from Ferguson 2020
            omega = 6.5; pm = (1/0.65)^2;
        case 5
            % EVD distribution from van Kerkhove 2015
            omega = 15.3; pm = (omega^2)/(9.3^2);
    end
    % Gamma shape-scale
    shapepm(i) = pm; scalepm(i) = omega/shapepm(i);
    % Gamma distribution
    y{i} = gampdf(t, shapepm(i), scalepm(i));
    sumgam(i) = sum(y{i}); y{i} = y{i}/sumgam(i);
end

% Check mean and sd of serial distributions
[mu_g, sig_g] = gamstat(shapepm, scalepm);
sig_g = sqrt(sig_g);

figure; hold on;
for i = type
    stairs(t, y{i}, 'LineWidth', 2);
end
grid off; box off; hold off;
xlabel('$u$ (days)', 'FontSize', 18);
ylabel('$w_u$ (days)', 'FontSize', 18);
legend(distNam{:}, 'Location', 'best', 'FontSize', 16, 'Box', 'off');