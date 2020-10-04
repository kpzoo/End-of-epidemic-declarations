% Function to compute end-of-epidem with different incidence
function [z, z1, zseq, didEnd, R, tdecs, tdeccheck, tst] = endIncid(Iday, Ilam, k, priors, sidist, mu)

% Assumptions and notes
% - include tst time for checking the epidemic end
% - Ilam is total infected cases that informs Lam
% - enter either total or local infections only for Iday

% Start checking for epidemic end (total cases)
tst = find(Ilam > 0, 1, 'last');

% Times to interrogate for z
tindex = tst+1:length(Ilam); 
lenind = length(tindex);

% Elimination probabilities and R
R = cell(1, lenind); zseq = R; 
z = zeros(1, lenind); z1 = z; 

% For t > tst get elimination probability
for i = 1:lenind
    % True data considered
    Itemp = Iday(1:tindex(i));
    % Total cases input
    Itemptot = Ilam(1:tindex(i));
    
    % Estimates of elimination probability
    [z(i), zseq{i}, z1(i), ~, R{i}] = getProbEliminSI(k, Itemp, Itemptot, priors, sidist);
end

% Find points of declaration
tdecs = find(z > mu, 1, 'first'); tdeccheck = find(z1 > mu, 1, 'first');
tdecset = [tdecs tdeccheck];
%disp(['Declared at [z z1]: ' num2str(tdecset)]);

% Check if epidemic actually ended
if any(isempty(tdecset))
    didEnd = 0;
    disp('Epidemic did not end');
else
    didEnd = 1;
end