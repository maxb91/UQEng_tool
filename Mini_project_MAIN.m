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
% PCE coefficients can be obtained from the following equation:
c=pinv(transpose(Psi)*Psi)*transpose(Psi)*y_ED;

%% 2.6 Numerical Experiments
%% 2.6.1 Validation: evaluation of the metamodel on a new set of points
% Value of PCE from ED
y_PCE=0;
for j=1:P
   y_PCE=y_PCE+c(j)*Psi(:,j);
end
% Number of new sample points
nv=10E4;
% Get sampled points in range [0,1]
u01_val = lhsdesign(nv,M);
%Create a large validation set XV
xv_sigsr = sigsr_lim(1)+u01_val(:,1)*(sigsr_lim(2)-sigsr_lim(1));
xv_srm = srm_lim(1)+u01_val(:,2)*(srm_lim(2)-srm_lim(1));
xv_taub1 = fct_lim(1)+u01_val(:,3)*(fct_lim(2)-fct_lim(1));
XV= [xv_sigsr, xv_srm, xv_taub1];
%% Evaluate the computational model M as well as the metamodel M_PCE on it
% Computational Model
y_ED_val=M_avg_strain(xv_sigsr,xv_srm,xv_taub1);
%PCE metalmodel
XV_uniform = XV;
for i=1:3
    XV_uniform(:,i) = ((XV(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
end
% Set up Psi matrix based on validation set
Psi_val = ones(nv,P);
for i=1:nv
    for j=1:P
        Psi_val(i,j) = eval_legendre(XV_uniform(i,1),p_index(j,1))*...
                   eval_legendre(XV_uniform(i,2),p_index(j,2))*...
                   eval_legendre(XV_uniform(i,3),p_index(j,3));
    end
end
% Value of PCE based on validation set 
y_PCE_val=0;
for j=1:P
   y_PCE_val=y_PCE_val+c(j)*Psi_val(:,j);
end
%% Compare the model responses with the responses of the metamodel for increasing polynomial degrees
figure
hold on
M_mean = zeros(1,3);
M_variance=zeros(1,3);
for p = 1:3 % polynomial degree
   
    P = factorial(M+p)/(factorial(M)*factorial(p));% Compute cardinality
    u01_val = lhsdesign(nv,M);% Get sampled points in range [0,1]
    
    %Create a large validation set XV
    xv_sigsr = sigsr_lim(1)+u01_val(:,1)*(sigsr_lim(2)-sigsr_lim(1));
    xv_srm = srm_lim(1)+u01_val(:,2)*(srm_lim(2)-srm_lim(1));
    xv_taub1 = fct_lim(1)+u01_val(:,3)*(fct_lim(2)-fct_lim(1));
    XV= [xv_sigsr, xv_srm, xv_taub1];

    % Computational Model Response
    y_ED_val=M_avg_strain(xv_sigsr,xv_srm,xv_taub1);
    %PCE metalmodel
    XV_uniform = XV;
    for i=1:3
        XV_uniform(:,i) = ((XV(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
    end
    % Get alphas for Legendre polynomials
    [p_index, p_index_roots] = create_alphas(M, p);
    % Set up Psi matrix based on validation set
    Psi_val = ones(nv,P);
    for i=1:nv
        for j=1:P
            Psi_val(i,j) = eval_legendre(XV_uniform(i,1),p_index(j,1))*...
                eval_legendre(XV_uniform(i,2),p_index(j,2))*...
                eval_legendre(XV_uniform(i,3),p_index(j,3));
        end
    end
    % Value of PCE based on validation set
    y_PCE_val=0;
    for j=1:P
        y_PCE_val=y_PCE_val+c(j)*Psi_val(:,j);
    end
    M_mean(p) = mean(y_ED_val);
    M_variance(p)=var(y_ED_val);
    plot(y_ED_val,y_PCE_val,'.')
end
legend('p=1','p=2','p=3')  
title('Comparison of models regarding increasing polynomial degrees')
ylabel('Computational Model Response')
xlabel('Metamodel Response')
hold off

%% Compute the relative leave-one-out (LOO) error using X_ED
A=Psi;
h=diag(A*pinv(transpose(A)*A)*transpose(A));
epsilon_LOO=0;
for i=1:n
    epsilon_LOO=epsilon_LOO+1/n*((y_ED(i)-y_PCE(i))/(1-h(i)))^2;
end 
%% Compute the relative mean-squared error (validation error) evaluated on the validation set XV
% Validation error
a=0;
b=0;
for i=1:nv
    a=a+(y_ED_val(i)-y_PCE_val(i)).^2;
    b=b+(y_ED_val(i)-mean(y_ED_val)).^2;
end
epsilon_val=a/b;
if epsilon_LOO >= epsilon_val
    fprintf('Leave-one-out error is greater than or equal to validation error.\n')
else
    fprintf('Leave-one-out error is smaller than validation error.\n')
end
%% 2.6.2 Validation: postprocessing of the coefficients
% check mean and variance of c
c_mean=mean(c);
c_variance=var(c);
% check mean and variance of M
%M_mean = mean(y_ED_val);
%M_variance=var(y_ED_val);
% PCE mean value
PCE_mean=c(1);
PCE_variance=sumsqr(c);
%% convergence plots for mean and variance regarding polynomial degrees
moment_mean=zeros(1,3);
moment_variance=zeros(1,3);
polynomials=[1,2,3];
% loop for various p
for p=1:3
    moment_mean(p)=M_mean(p)/PCE_mean-1;
    moment_variance(p)=M_variance(p)/PCE_variance-1;
end
figure
subplot(2,1,1)
plot(polynomials,moment_mean)
xlabel('polynomial degree')
title('Convergence plot for mean')
subplot(2,1,2)
plot(polynomials,moment_variance)
title('Convergence plot for variance')
xlabel('polynomial degree')
%% convergence plots for mean and variance regarding number of samples
i=1;
% loop for various sample values
for n_samples=[n,10,1e2,1e3,1e4]
    u01_n = lhsdesign(n_samples,M);
    xn_sigsr = sigsr_lim(1)+u01_n(:,1)*(sigsr_lim(2)-sigsr_lim(1));
    xn_srm = srm_lim(1)+u01_n(:,2)*(srm_lim(2)-srm_lim(1));
    xn_taub1 = fct_lim(1)+u01_n(:,3)*(fct_lim(2)-fct_lim(1));
    Xn= [xn_sigsr, xn_srm, xn_taub1];
    % Computational Model Response
    y_ED_n=M_avg_strain(xn_sigsr,xn_srm,xn_taub1);
    M_mean_n(i) = mean(y_ED_n);
    M_variance_n(i)=var(y_ED_n);
    i=i+1;
end
figure
subplot(2,1,1)
plot([n,10,1e2,1e3,1e4],M_mean_n)
xlabel('number of samples')
title('Convergence plot for mean')
subplot(2,1,2)
plot([n,10,1e2,1e3,1e4],M_variance_n)
title('Convergence plot for variance')
xlabel('number of samples')
