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
Qrange = 900; % kN Range of the load cell
acc = 0.001; % RO of the load cell
Qacc = Qrange*acc; % resolution of load cell
sigacc = 1e3*Qacc/(pi/4*D^2); % Assuming mu+sigacc = 95% quantile
stdev_sig = sigacc/norminv(0.95);
sigmeas = 590; % Measured value
sigsr_lim = [sigmeas-sigacc, sigmeas+sigacc];

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
b_tau = exp(norminv(0.95)*fct_zeta+fct_lambda);
a_tau = exp(norminv(0.05)*fct_zeta+fct_lambda);

tau_lim = [a_tau,  b_tau];

sampling_limits = [sigsr_lim;
                   srm_lim;
                   tau_lim];

%%  Figure 1
Mini_project_FIG1

%% 2.3 Sampling the experimental design
% Latin hypercube sampling
p = 4; % polynomial degree
M = 3; % 3 uncertain variables

% Compute cardinality:
P = factorial(M+p)/(factorial(M)*factorial(p));

k = 3; % oversampling rate
n = k*P; % resulting number of samples

% Get sampled points in range [0,1]
u01 = lhsdesign(n,M);

% Transform -- Experimental design
x_sigsr = sigsr_lim(1)+u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
x_srm = srm_lim(1)+u01(:,2)*(srm_lim(2)-srm_lim(1));
x_taub1 = tau_lim(1)+u01(:,3)*(tau_lim(2)-tau_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
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
xMC_taub1 = tau_lim(1)+u01MC(:,3)*(tau_lim(2)-tau_lim(1)); 
XMC = [xMC_sigsr, xMC_srm, xMC_taub1];

%%  Figure 2
Mini_project_FIG2
figure
plotmatrix(X)
grid on

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
xv_taub1 = tau_lim(1)+u01_val(:,3)*(tau_lim(2)-tau_lim(1));
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


% figure
% plot(y_ED_val,y_PCE_val,'.')
% 
% 
% figure
% plot(y_ED,y_PCE,'.')
%% Compare the model responses with the responses of the metamodel for increasing polynomial degrees
degrees = 1:5;
incdeg = arrayfun(@(x)struct('p',x),degrees);
for pp = degrees % polynomial degree
    tic;
    % FIRST CREATE PCE MODEL OF DEGREE p
    incdeg(pp).P = factorial(M+pp)/(factorial(M)*factorial(pp));% Compute cardinality
    incdeg(pp).n = k*incdeg(pp).P; % resulting number of samples

    % Get sampled points in range [0,1]
    incdeg(pp).u01 = lhsdesign(incdeg(pp).n,M);

    % Transform -- Experimental design
    incdeg(pp).x_sigsr = sigsr_lim(1)+incdeg(pp).u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
    incdeg(pp).x_srm = srm_lim(1)+incdeg(pp).u01(:,2)*(srm_lim(2)-srm_lim(1));
    incdeg(pp).x_taub1 = tau_lim(1)+incdeg(pp).u01(:,3)*(tau_lim(2)-tau_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
    incdeg(pp).X = [incdeg(pp).x_sigsr, incdeg(pp).x_srm, incdeg(pp).x_taub1];
    incdeg(pp).y_ED  =M_avg_strain(incdeg(pp).x_sigsr,incdeg(pp).x_srm,incdeg(pp).x_taub1);
    
    % Transform to [-1,1] for Legendre
    incdeg(pp).X_uniform = incdeg(pp).X;
    for i=1:3
        incdeg(pp).X_uniform(:,i) = ((incdeg(pp).X(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
    end

    % Get alphas for Legendre polynomials
    [incdeg(pp).p_index, incdeg(pp).p_index_roots] = create_alphas(M, pp);

    % Set up Psi matrix
    incdeg(pp).Psi = ones(incdeg(pp).n,incdeg(pp).P);
    for i=1:incdeg(pp).n
        for j=1:incdeg(pp).P
            % For each element of the matrix, compute the multivariate basis by
            % a product of the univariate ones:
            incdeg(pp).Psi(i,j) = eval_legendre(incdeg(pp).X_uniform(i,1),incdeg(pp).p_index(j,1))*...
                eval_legendre(incdeg(pp).X_uniform(i,2),incdeg(pp).p_index(j,2))*...
                eval_legendre(incdeg(pp).X_uniform(i,3),incdeg(pp).p_index(j,3));
        end
    end

    % PCE coefficients can be obtained from the following equation:
    incdeg(pp).c=pinv(transpose(incdeg(pp).Psi)*incdeg(pp).Psi)*transpose(incdeg(pp).Psi)*incdeg(pp).y_ED;
    
    incdeg(pp).t=toc;

    % Evaluate PCE on ED
    incdeg(pp).y_PCE=0;
    for j=1:incdeg(pp).P
        incdeg(pp).y_PCE=incdeg(pp).y_PCE+incdeg(pp).c(j)*incdeg(pp).Psi(:,j);
    end

    % SECOND: EVALUATE THE PCE MODEL ON VALIDTION SAMPLE   

    % Set up Psi matrix based on validation set (always the same validation
    % set)
    incdeg(pp).Psi_val = ones(nv,incdeg(pp).P);
    for i=1:nv
        for j=1:incdeg(pp).P
            incdeg(pp).Psi_val(i,j) = eval_legendre(XV_uniform(i,1),incdeg(pp).p_index(j,1))*...
                eval_legendre(XV_uniform(i,2),incdeg(pp).p_index(j,2))*...
                eval_legendre(XV_uniform(i,3),incdeg(pp).p_index(j,3));
        end
    end

    % Value of PCE based on validation set
    incdeg(pp).y_PCE_val=0;
    for j=1:incdeg(pp).P
        incdeg(pp).y_PCE_val=incdeg(pp).y_PCE_val+incdeg(pp).c(j)*incdeg(pp).Psi_val(:,j);
    end
    incdeg(pp).mu_PCE = mean(incdeg(pp).y_PCE_val);
    incdeg(pp).var_PCE=var(incdeg(pp).y_PCE_val);
    
    fprintf('PCE with degree %i done, mu = %4.4f, var = %4.4e, t = %4.2f \n',pp,incdeg(pp).mu_PCE,incdeg(pp).var_PCE,incdeg(pp).t)
end

for pp = 1:length(degrees)
    incdeg(pp).muhat_PCE = incdeg(pp).c(1);
    incdeg(pp).varhat_PCE =sum(incdeg(pp).c(2:end).*incdeg(pp).c(2:end));
    incdeg(pp).cvhat_PCE =sqrt(incdeg(pp).varhat_PCE)/incdeg(pp).muhat_PCE;
    incdeg(pp).h=diag(incdeg(pp).Psi*pinv(transpose(incdeg(pp).Psi)*incdeg(pp).Psi)*transpose(incdeg(pp).Psi));
    incdeg(pp).epsLOO = (1./incdeg(pp).n).*...
        sum(((incdeg(pp).y_ED-incdeg(pp).y_PCE)./(1-incdeg(pp).h)).^2)./...
        sum((incdeg(pp).y_ED-mean(incdeg(pp).y_ED)).^2);
    incdeg(pp).epsMSQ = sum((y_ED_val-incdeg(pp).y_PCE_val).^2)./...
        sum((y_ED_val-mean(y_ED_val)).^2);
    incdeg(pp).epsMSQ_ED = sum((incdeg(pp).y_ED-incdeg(pp).y_PCE).^2)./...
        sum((incdeg(pp).y_ED-mean(incdeg(pp).y_ED)).^2);
end

%% Figure 4 - Comparison degrees

Mini_project_FIG4
Mini_project_FIG4b
Mini_project_FIG4c % This one is the best, right?

%% Compute the relative leave-one-out (LOO) error using X_ED
A=Psi;
h=diag(A*pinv(transpose(A)*A)*transpose(A));
E_LOO=0;
d=0;
for i=1:n
    E_LOO=E_LOO+1/n*((y_ED(i)-y_PCE(i))/(1-h(i)))^2;
    d=d+(y_ED(i)-y_ED_mean)^2;
end 
epsilon_LOO=E_LOO/d;

epsLOO = (1./n).*...
        sum(((y_ED-y_PCE)./(1-h)).^2)./...
        sum((y_ED-mean(y_ED)).^2);

%% Compute the relative mean-squared error (validation error) evaluated on the validation set XV
% Validation error
a=0;
b=0;
mu_M = mean(y_ED_val);
for i=1:nv
    a=a+(y_ED_val(i)-y_PCE_val(i)).^2;
    b=b+(y_ED_val(i)-mu_M).^2;
end
epsilon_val=a/b;
if epsilon_LOO >= epsilon_val
    fprintf('Leave-one-out error is greater than or equal to validation error.\n')
else
    fprintf('Leave-one-out error is smaller than validation error.\n\n')
end
%% 2.6.2 Validation: postprocessing of the coefficients
% check mean and variance of M
% mu_M defined before
var_M=var(y_ED_val);
% PCE mean value
mu_PCE=c(1);
var_PCE=sum(c(2:end).*c(2:end));
fprintf('Comparison of the mean value\n')
fprintf('mu_M = %4.4f \n',mu_M)
fprintf('mu_PCE = %4.4f \n',mu_PCE)
fprintf('Comparison of the variance\n')
fprintf('var_M = %4.4e \n',var_M)
fprintf('var_PCE = %4.4e \n\n',var_PCE)
%% convergence plots for mean and variance regarding polynomial degrees

polynomials=degrees;


moment_mean=abs([incdeg.mu_PCE]/mu_M-1);
moment_variance=abs([incdeg.var_PCE]/var_M-1);


figure
hold on
box on
grid on
plot(polynomials,moment_mean,'k-o')
plot(polynomials,moment_variance,'-s','color',[0.7 0.7 0.7])
legend('Mean value','Variance')
xlabel('polynomial degree {\it p}')
ylabel('Relative error')
set(gca,'yscale','log')


%% convergence plots for mean and variance regarding number of samples
p = 4;
oversampling = 0.5:0.5:4;
incn = arrayfun(@(x)struct('oversampling',x),oversampling);
for kk = 1:length(oversampling) 
    tic;
    incn(kk).total_degree = p;

    % FIRST CREATE PCE MODEL OF DEGREE p
    incn(kk).P = factorial(M+p)/(factorial(M)*factorial(p));% Compute cardinality
    incn(kk).n = ceil(incn(kk).oversampling*incn(kk).P); % resulting number of samples

    % Get sampled points in range [0,1]
    incn(kk).u01 = lhsdesign(incn(kk).n,M);

    % Transform -- Experimental design
    incn(kk).x_sigsr = sigsr_lim(1)+incn(kk).u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
    incn(kk).x_srm = srm_lim(1)+incn(kk).u01(:,2)*(srm_lim(2)-srm_lim(1));
    incn(kk).x_taub1 = tau_lim(1)+incn(kk).u01(:,3)*(tau_lim(2)-tau_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
    incn(kk).X = [incn(kk).x_sigsr, incn(kk).x_srm, incn(kk).x_taub1];
    incn(kk).y_ED  =M_avg_strain(incn(kk).x_sigsr,incn(kk).x_srm,incn(kk).x_taub1);
    
    % Transform to [-1,1] for Legendre
    incn(kk).X_uniform = incn(kk).X;
    for i=1:3
        incn(kk).X_uniform(:,i) = ((incn(kk).X(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
    end

    % Get alphas for Legendre polynomials
    [incn(kk).p_index, incn(kk).p_index_roots] = create_alphas(M, p);

    % Set up Psi matrix
    incn(kk).Psi = ones(incn(kk).n,incn(kk).P);
    for i=1:incn(kk).n
        for j=1:incn(kk).P
            % For each element of the matrix, compute the multivariate basis by
            % a product of the univariate ones:
            incn(kk).Psi(i,j) = eval_legendre(incn(kk).X_uniform(i,1),incn(kk).p_index(j,1))*...
                eval_legendre(incn(kk).X_uniform(i,2),incn(kk).p_index(j,2))*...
                eval_legendre(incn(kk).X_uniform(i,3),incn(kk).p_index(j,3));
        end 
    end

    % PCE coefficients can be obtained from the following equation:
    incn(kk).c=pinv(transpose(incn(kk).Psi)*incn(kk).Psi)*transpose(incn(kk).Psi)*incn(kk).y_ED;

    incn(kk).t=toc;

    % Evaluate PCE on ED
    incn(kk).y_PCE=0;
    for j=1:incn(kk).P
        incn(kk).y_PCE=incn(kk).y_PCE+incn(kk).c(j)*incn(kk).Psi(:,j);
    end
    
    % SECOND: EVALUATE THE PCE MODEL ON VALIDTION SAMPLE   

    % Set up Psi matrix based on validation set (always the same validation
    % set)
    incn(kk).Psi_val = ones(nv,incn(kk).P);
    for i=1:nv
        for j=1:incn(kk).P
            incn(kk).Psi_val(i,j) = eval_legendre(XV_uniform(i,1),incn(kk).p_index(j,1))*...
                eval_legendre(XV_uniform(i,2),incn(kk).p_index(j,2))*...
                eval_legendre(XV_uniform(i,3),incn(kk).p_index(j,3));
        end
    end

    % Value of PCE based on validation set
    incn(kk).y_PCE_val=0;
    for j=1:incn(kk).P
        incn(kk).y_PCE_val=incn(kk).y_PCE_val+incn(kk).c(j)*incn(kk).Psi_val(:,j);
    end
    incn(kk).mu_PCE = mean(incn(kk).y_PCE_val);
    incn(kk).var_PCE=var(incn(kk).y_PCE_val);
    
    fprintf('PCE with experimental design of %i samples done, mu = %4.4f, var = %4.4e, t = %4.2f \n',incn(kk).n,incn(kk).mu_PCE,incn(kk).var_PCE,incn(kk).t)
end

for kk = 1:length(oversampling)
        % Evaluate PCE on ED
    incn(kk).muhat_PCE = incn(kk).c(1);
    incn(kk).varhat_PCE =sum(incn(kk).c(2:end).*incn(kk).c(2:end));
    incn(kk).cvhat_PCE =sqrt(incn(kk).varhat_PCE)/incn(kk).muhat_PCE;
    incn(kk).h=diag(incn(kk).Psi*pinv(transpose(incn(kk).Psi)*incn(kk).Psi)*transpose(incn(kk).Psi));
    incn(kk).epsLOO = (1./incn(kk).n).*...
        sum(((incn(kk).y_ED-incn(kk).y_PCE)./(1-incn(kk).h)).^2)./...
        sum((incn(kk).y_ED-mean(incn(kk).y_ED)).^2);
    incn(kk).epsMSQ = sum((y_ED_val-incn(kk).y_PCE_val).^2)./...
        sum((y_ED_val-mean(y_ED_val)).^2);
    incn(kk).epsMSQ_ED = sum((incn(kk).y_ED-incn(kk).y_PCE).^2)./...
        sum((incn(kk).y_ED-mean(incn(kk).y_ED)).^2);
end

rel_error_n_mu = abs([incn.mu_PCE]./mu_M-1);
rel_error_n_var = abs([incn.var_PCE]./var_M-1);

figure
hold on
box on
grid on
plot([incn.n],rel_error_n_mu,'k-o')
plot([incn.n],rel_error_n_var,'-s','color',[0.7 0.7 0.7])
legend('Mean value','Variance')
xlabel('Number of samples in the experimental design {\it n}')
ylabel('Relative error')
set(gca,'yscale','log')
Mini_project_FIG5

%% "Exact value" of mu, sigma
tic;
% Number of new sample points
ns=1E6;
% Get sampled points in range [0,1]
u01_s = lhsdesign(ns,M);
%Create a large validation set XV
xvs_sigsr = sigsr_lim(1)+u01_s(:,1)*(sigsr_lim(2)-sigsr_lim(1));
xvs_srm = srm_lim(1)+u01_s(:,2)*(srm_lim(2)-srm_lim(1));
xvs_taub1 = tau_lim(1)+u01_s(:,3)*(tau_lim(2)-tau_lim(1));
XVs= [xvs_sigsr, xvs_srm, xvs_taub1];
% Evaluate the computational model M as well as the metamodel M_PCE on it
% Computational Model
y_ED_s=M_avg_strain(xvs_sigsr,xvs_srm,xvs_taub1);
mu_exact = mean(y_ED_s);
var_exact = var(y_ED_s);
t_exact = toc;
%% Possible additional evaluation / postprocessing




Mini_project_ADD
