clear; close all;

uqlab;

%% Define the model 
modelopts.mFile = 'TensionChordModel';
uqTCM = uq_createModel(modelopts);


%% Boundaries of uniform distributions 
D = 16; % Reinforcing bar diameter in mm

% Steel stress at the crack
Qrange = 900; % kN Range of the load cell
acc = 0.001; % RO of the load cell
Qacc = Qrange*acc; % resolution of load cell
sigacc = 1e3*Qacc/(pi/4*D^2); % Assuming mu+sigacc = 95% quantile
stdev_sig = sigacc/norminv(0.95);
sigmeas = 590; % Measured value
sigsr_lim = [sigmeas-sigacc, sigmeas+sigacc];

% Crack spacing

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
b_tau = exp(norminv(0.95)*fct_zeta+fct_lambda);
a_tau = exp(norminv(0.05)*fct_zeta+fct_lambda);

tau_lim = [a_tau, b_tau];


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
inputopts.Marginals(3).Parameters = tau_lim;

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


n = 300; %with k=2 enough for p <= 9, relevant later on for convergence studies 
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

metaopts.Degree = 4 ;
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


%% Calculate "exact" Model Response
% n_large = 10E6;
% X_large = uq_getSample(uqTCM_input, n_large, 'LHS');
% Y_exact = uq_evalModel(uqTCM, X_large);
% mu_exact = mean(Y_exact);
mu_exact = 0.011556509939278;
% var_exact = var(Y_exact);
var_exact = 1.641797934414332e-05;

%% Convergence study 

% set validation dataset
metaopts.ValidationSet.X = X_val;
metaopts.ValidationSet.Y = Y_val;


%% fixed n in ED, varying degree

n = 300; %with k=2 enough for p <= 9, relevant later on for convergence studies 
X_ED = uq_getSample(uqTCM_input, n, 'LHS');
Y_ED = uq_evalModel(uqTCM, X_ED);
metaopts.ExpDesign.X = X_ED;
metaopts.ExpDesign.Y = Y_ED;


degrees = 1:8;
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
yline(mu_exact)
set(gca, 'yscale', 'log')
xlabel('Degree')
ylabel('mean')
subplot(2,2,2)
plot(degrees, var, 'o-')
yline(var_exact)
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
n = [35 53 70 88 105 123 140 300 500 700 1000];
metaopts.Degree = 4;

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
yline(mu_exact)
set(gca, 'yscale', 'log')
xlabel('N')
ylabel('mean')
subplot(2,2,2)
plot(n, var, 'o-')
yline(var_exact)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% RELEVANT CODE FOR REPORT %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% clear metaopts and start fresh
clear metaopts X_ED Y_ED
rng(1)
%%
%model settings
metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.Method = 'OLS';

%create ED
n = 70; %as determined in report 
X_ED = uq_getSample(uqTCM_input, n, 'LHS');
Y_ED = uq_evalModel(uqTCM, X_ED);
metaopts.ExpDesign.X = X_ED;
metaopts.ExpDesign.Y = Y_ED;

% set validation dataset
metaopts.ValidationSet.X = X_val;
metaopts.ValidationSet.Y = Y_val;

%use adaptive degrees
metaopts.Degree = 1:8;

TCM_OLSPCE = uq_createModel(metaopts);

uq_print(TCM_OLSPCE)

% -> UQLab also chooses a degree of 4, seems like we made the right choice 

%% See improvements if we let UQLab choose truncation (if any)
metaopts.Degree = 4; %fix Degree
metaopts.TruncOptions.qNorm = 0.5:0.05:1.0;

TCM_OLSPCE = uq_createModel(metaopts);
uq_print(TCM_OLSPCE)

% -> uses trunqation of q= 0.8


%% Use LARS, with same settings, see difference
metaopts.TruncOptions.qNorm = 0.8 % fix q-norm
metaopts.Method = 'LARS';

TCM_LARSPCE = uq_createModel(metaopts);
uq_print(TCM_LARSPCE)

% -> uses less coefficients than OLS (22 vs 23) 


%% Show coefficients for LARS vs OLS 

hold on 
plot(abs(TCM_OLSPCE.PCE.Coefficients), 'bo')
plot(abs(TCM_LARSPCE.PCE.Coefficients), 'ro')
set(gca,'yscale','log')
hold off

% -> all coefficients are rather similar/identical, which makes sense considering
% that the sparse matrix is close the the full matrix

%% Increase ED size, let UQLab have some fun, see improvements, if any 

clear metaopts 

metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.Method = 'LARS';

metaopts.Degree = 1:10;
metaopts.TruncOptions.qNorm = 0.5:0.05:1.0;

metaopts.Input = uqTCM_input;
metaopts.FullModel = uqTCM;
metaopts.ExpDesign.NSamples = 200; % stil a reasonable runtime

TCM_LARSPCE = uq_createModel(metaopts);
uq_print(TCM_LARSPCE)


% important observation -> Maximal degree as well as q-norm selected by
% UQLab is very dependent on the sample in the ED -> fixed by setting seed for rng 
% -problem? 






%% Convergence study: fixed N, increasing degree -- LARS vs OLS
% Additional experiment (not done in Live session)


degrees = 1:7;
N = 200;


% set validation dataset
metaopts.ValidationSet.X = X_val;
metaopts.ValidationSet.Y = Y_val;


for dd = 1:numel(degrees)
    metaopts.Degree = degrees(dd);
    
    metaopts.Method = 'OLS';
    myOLSPCE = uq_createModel(metaopts);
    mean(dd,1) = myOLSPCE.PCE.Moments.Mean;
    valError(dd,1) = myOLSPCE.Error.Val; 
    
    metaopts.Method = 'LARS';
    myLARSPCE = uq_createModel(metaopts);
    mean(dd,2) = myLARSPCE.PCE.Moments.Mean;
    valError(dd,2) = myLARSPCE.Error.Val; 
    
end



%% Visualize results
figure
plot(degrees, abs(mean - mu_exact)/abs(mu_exact), 'o-')
set(gca, 'yscale', 'log')
xlabel('Degree')
ylabel('Rel. error of mean')
legend('OLS', 'LARS')

figure
plot(degrees, valError, 'o-')
set(gca, 'yscale', 'log')
xlabel('Degree')
ylabel('validation error')
legend('OLS', 'LARS')


%% change input to lognormal & compare difference to all uniform

% Beam length in m
inputopts.Marginals(3).Name = 'taub1';
inputopts.Marginals(3).Type = 'lognormal';
inputopts.Marginals(3).Parameters = [fct_lambda, fct_zeta];

uqTCM_input = uq_createInput(inputopts);


n = 70; %with k=2 enough for p <= 9, relevant later on for convergence studies 
X_ED = uq_getSample(uqTCM_input, n, 'LHS');


% Evaluate the model on the sample X
Y_ED = uq_evalModel(uqTCM, X_ED);


% PCE
metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.Method = 'OLS';

%use already created Experimental Design 
metaopts.ExpDesign.X = X_ED;
metaopts.ExpDesign.Y = Y_ED;

metaopts.Degree = 4 ;
metaopts.Input = uqTCM_input;

TCM_OLSPCE = uq_createModel(metaopts);
Y_PCE_ED = uq_evalModel(TCM_OLSPCE, X_ED);

% Check the fit
plot(Y_ED, Y_PCE_ED, 'o')
xlabel('model response')
ylabel('PCE response')


% -> more comparsion necessary plots 

%metaopts.Degree = 1:8;
%metaopts.TruncOptions.qNorm = 0.5:0.1:1.0;
% Truncation options:
% metaopts.TruncOptions.qNorm = 0.5;
% metaopts.TruncOptions.MaxInteraction = 2;

