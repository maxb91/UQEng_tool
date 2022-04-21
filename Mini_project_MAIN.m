% UQ Eng mini project

clear
clc
close all
addpath('fcn')

%% 2.1 Choice and implementation of the model

% Fixed Parameters
D = 16; % Reinforcing bar diameter in mm
Es = 205000; % Elastic modulus of the steel in MPa
esu = 0.08; % Ultimate strain of the steel
fsu = 600; % Ultimate stress of the steel in MPa
fsy = 550; % Yield stress of the steel in MPa
rho = 0.02; % Reinforcement content As/Ac
% Material parameters for Ramberg-Osgood steel material law
ka = 0.002;
kb = 0.002;
alpha = log((esu-fsu/Es)/ka)/log(fsu/fsy);
kc =fsy/(kb^(1/alpha));

% Model - Alvarez 1998
% "Tension Chord Model" - calculate average strain in the reinforcing bar
M_avg_strain = @(sigsr,srm,taub1) sigsr/Es-taub1.*srm/(Es*D)+D./(2*taub1.*srm)*1./(alpha+1)./kc.^alpha.*...
    (sigsr.^(alpha+1)-(sigsr-2*taub1.*srm/D).^(alpha+1));

% Main uncertainties:
% sigma_sr (steel stress at the crack [MPa])
% s_rm (crack spacing [mm]) - distance between where cracks in the concrete
% occur
% tau_b1 (bond shear stress [MPa])


%% 2.2 Probabilistic input model

% Steel stress at the crack
warning('To be changed')
sigsr_lim = [585, 595];

% Crack spacing
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

sampling_limits = [sigsr_lim;
                   srm_lim;
                   fct_lim];

%%  Figure 1
Mini_project_FIG1

%% 2.3 Sampling the experimental design
% Latin hypercube sampling
warning('??')
p = 3; % polynomial degree
M = 3; % 3 uncertain variables

% Compute cardinality:
P = factorial(M+p)/(factorial(M)*factorial(p));

k = 2; % oversampling rate
n = k*P; % resulting number of samples

% Get sampled points in range [0,1]
u01 = lhsdesign(n,M);

% Transform -- Experimental design
x_sigsr = sigsr_lim(1)+u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
x_srm = srm_lim(1)+u01(:,2)*(srm_lim(2)-srm_lim(1));
x_taub1 = fct_lim(1)+u01(:,3)*(fct_lim(2)-fct_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
X = [x_sigsr, x_srm, x_taub1];

% Evaluation of the full model on the experimental design
y_ED = M_avg_strain(x_sigsr,x_srm,x_taub1);
y_ED_mean = mean(y_ED);
y_ED_std = std(y_ED);
y_ED_var = y_ED_std^2;

% for comparison: do the same thing with Monte Carlo:
u01MC = rand(n,M);
xMC_sigsr = sigsr_lim(1)+u01MC(:,1)*(sigsr_lim(2)-sigsr_lim(1));
xMC_srm = srm_lim(1)+u01MC(:,2)*(srm_lim(2)-srm_lim(1));
xMC_taub1 = fct_lim(1)+u01MC(:,3)*(fct_lim(2)-fct_lim(1)); 

%%  Figure 2
Mini_project_FIG2

%%  Figure 3
Mini_project_FIG3

%% 2.4 Polynomial Basis
%% 2.4.2 Isoprobabilistic transform
% Transform random variables back into range [-1,1]
% X_uniform corresponds to Xi in eq. (3) of the instructions
X_uniform = X;
for i=1:3
    X_uniform(:,i) = ((X(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
end

%% 2.4.3 Polynomial degree and basis truncation
% Get alphas for Legendre polynomials
[p_index, p_index_roots] = create_alphas(M, p);

% Set up Psi matrix
Psi = ones(n,P);
for i=1:n
    for j=1:P
        % For each element of the matrix, compute the multivariate basis by
        % a product of the univariate ones:
        Psi(i,j) = eval_legendre(X_uniform(i,1),p_index(j,1))*...
                   eval_legendre(X_uniform(i,2),p_index(j,2))*...
                   eval_legendre(X_uniform(i,3),p_index(j,3));
    end
end


%% 2.5 Compute polynomial coefficients


