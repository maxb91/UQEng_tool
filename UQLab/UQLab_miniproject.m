clear; close all;

uqlab;

%% Define the model 
modelopts.mFile = 'TensionChordModel';
uqTCM = uq_createModel(modelopts);


%% Boundaries of uniform distributions 
% Steel stress at the crack
warning('To be changed')
sigsr_lim = [585, 595];

% Crack spacing
D = 16; % Reinforcing bar diameter in mm
rho = 0.02; % Reinforcement content As/Ac
sr0 = D/4*(1/rho-1); % Maximum crack spacing in mm
srmin = 0.5*sr0; % Minimum crack spacing in mm
srm_lim = [srmin, sr0];

% Nominal bond shear stress
alphaJCSS = 1;
lambdaJCSS = 0.96;
% Lognormal distribution for measured compressive strength of concrete
fc_mean = 40; % mean compressive strength of concrete, from an assumed material test series
fc_std = 5; % standard deviation of compressive strength of concrete, from an assumed material test series
fc_cov = fc_std/fc_mean;
fc_zeta = sqrt(log(1+fc_cov^2));
fc_lambda = log(fc_mean/sqrt(1+fc_cov^2));
% Lognormal distribution for Y1
Y1_mean = 1;
Y1_cov = 0.06;
Y1_std = Y1_mean*Y1_cov;
Y1_zeta = sqrt(log(1+Y1_cov^2));
Y1_lambda = log(Y1_mean/sqrt(1+Y1_cov^2));
% Lognormal distribution for Y2
Y2_mean = 1;
Y2_cov = 0.3;
Y2_std = Y1_mean*Y1_cov;
Y2_zeta = sqrt(log(1+Y2_cov^2));
Y2_lambda = log(Y2_mean/sqrt(1+Y2_cov^2));
% Transform to lognormal distribution of tensile stress
fct_zeta = 2/3*lambdaJCSS*fc_zeta+2/3*Y1_zeta+Y2_zeta;
fct_lambda = 2/3*lambdaJCSS*fc_lambda+2/3*Y1_lambda+Y2_lambda+log(0.3);
fct_mean = exp(fct_lambda+fct_zeta^2/2);
fct_std = exp(fct_lambda+fct_zeta^2/2)*sqrt(exp(fct_zeta^2)-1);
warning('to be discussed')
fct_lim = [fct_mean - fct_std,  fct_mean + fct_std];


%% Define the Input

% Steel stress at the crack
inputopts.Marginals(1).Name = 'sigsr';
inputopts.Marginals(1).Type = 'uniform';
inputopts.Marginals(1).Parameters = sigsr_lim;

% Crack spacing
inputopts.Marginals(2).Name = 'srm';
inputopts.Marginals(2).Type = 'uniform';
inputopts.Marginals(2).Parameters = srm_lim;

% Beam length in m
inputopts.Marginals(3).Name = 'taub1';
inputopts.Marginals(3).Type = 'uniform';
inputopts.Marginals(3).Parameters = fct_lim;

uqTCM_input = uq_createInput(inputopts);

%% Sample from the input

% p = 3; % polynomial degree
% M = 3; % 3 uncertain variables
% 
% % Compute cardinality:
% P = factorial(M+p)/(factorial(M)*factorial(p));
% 
% k = 2; % oversampling rate
% n = k*P; % resulting number of samples % evaluates to 40 -> low enough
% computational burden, so use n = 500 


n = 500; %with k=2 enough for p <= 9, relevant later on for convergence studies 
X_ED = uq_getSample(uqTCM_input, n, 'LHS');


% create a scatterplot of X2 against X1
% plot(X_ED(:,1), X_ED(:,2), 'o')
% 
% plot the sample: matrix of scatterplots
% plotmatrix(X_ED)

%% Evaluate the model on the sample X
Y_ED = uq_evalModel(uqTCM, X_ED);

%histogram(Y_ED)

%% PCE
metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.Method = 'OLS';

%use already created Experimental Design 
metaopts.ExpDesign.X = X_ED;
metaopts.ExpDesign.Y = Y_ED;

metaopts.Degree = 3;
metaopts.Input = uqTCM_input;

TCM_OLSPCE = uq_createModel(metaopts);
Y_PCE_ED = uq_evalModel(TCM_OLSPCE, X_ED);

% Check the fit
plot(Y_ED, Y_PCE_ED, 'o')
xlabel('model response')
ylabel('PCE response')

%% Validate model 
n=10E4;
X_val = uq_getSample(uqTCM_input, n, 'LHS');
Y_val = uq_evalModel(uqTCM, X_val);
Y_PCE_val = uq_evalModel(TCM_OLSPCE, X_val);

% Check the fit
plot(Y_val, Y_PCE_val, 'o')
xlabel('model response')
ylabel('PCE response')

%% Convergence study 

% set validation dataset
metaopts.ValidationSet.X = X_val;
metaopts.ValidationSet.Y = Y_val;

%% fixed n in ED, varying degree
degrees = 1:9;
hold on

for dd = 1:numel(degrees)
    metaopts.Degree = degrees(dd);
    TCM_OLSPCE(dd) = uq_createModel(metaopts); %model is always created with the same ED defined in metaopts
    LOO(dd) = TCM_OLSPCE(dd).Error.LOO;
    val_err(dd) = TCM_OLSPCE(dd).Error.Val;

    mu(dd) = TCM_OLSPCE(dd).PCE.Moments.Mean;
    var(dd) = TCM_OLSPCE(dd).PCE.Moments.Var;
    %plot(Y_val,Y_PCE_val,'.')
end

figure
subplot(2,2,1)
plot(degrees, mu, 'o-')
set(gca, 'yscale', 'log')
xlabel('Degree')
ylabel('mean')
subplot(2,2,2)
plot(degrees, var, 'o-')
set(gca, 'yscale', 'log')
xlabel('Degree')
ylabel('variance')
subplot(2,2,[3 4])
hold on 
plot(degrees, LOO, 'o-')
plot(degrees, val_err, '.-')
set(gca, 'yscale', 'log')
legend('LOO','val_err')  
xlabel('Degree')
ylabel('relatve error')


hold off 

%% fixed degree, varying n in ED
clear TCM_OLSPCE  LOO  val_err mu  var 
n = [50 100 200 500 1000];
metaopts.Degree = 3;

for dd = 1:numel(n)
    X_ED = uq_getSample(uqTCM_input, n(dd), 'LHS');
    Y_ED = uq_evalModel(uqTCM, X_ED);
    metaopts.ExpDesign.X = X_ED;
    metaopts.ExpDesign.Y = Y_ED;

    TCM_OLSPCE(dd) = uq_createModel(metaopts); %degree is kept constant 
    LOO(dd) = TCM_OLSPCE(dd).Error.LOO;
    val_err(dd) = TCM_OLSPCE(dd).Error.Val;  %note: validation set is kept the same

    mu(dd) = TCM_OLSPCE(dd).PCE.Moments.Mean;
    var(dd) = TCM_OLSPCE(dd).PCE.Moments.Var;
end

figure
subplot(2,2,1)
plot(n, mu, 'o-')
set(gca, 'yscale', 'log')
xlabel('N')
ylabel('mean')
subplot(2,2,2)
plot(n, var, 'o-')
set(gca, 'yscale', 'log')
xlabel('N')
ylabel('variance')
subplot(2,2,[3 4])
hold on 
plot(n, LOO, 'o-')
plot(n, val_err, '.-')
set(gca, 'yscale', 'log')
legend('LOO','val_err')  
xlabel('N')
ylabel('relatve error')


hold off 

