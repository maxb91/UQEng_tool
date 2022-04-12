% UQ Eng mini project

clear
clc
close all

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
M_avg_strain = @(sigsr,srm,taub1) sigsr/Es-taub1.*srm/(Es*D)+D./(2*taub1.*srm)*1./(alpha+1)./kc.^alpha.*...
    (sigsr.^(alpha+1)-(sigsr-2*taub1.*srm/D).^(alpha+1));



%% 2.2 Probabilistic input model

% Steel stress at the crack
warning('To be changed')
a_sigsr = 585;
b_sigsr = 595;

% Crack spacing
sr0 = D/4*(1/rho-1); % Maximum crack spacing in mm
srmin = 0.5*sr0; % Minimum crack spacing in mm
a_srm = srmin;
b_srm = sr0;

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
a_fct = fct_mean - fct_std;
b_fct = fct_mean + fct_std;


%% 2.2 - plots
F1.fig = figure('Units','Centimeters','Position',[2,2,21-5,6]);

F1.ax(1) = subplot(1,3,1);
grid on
box on
hold on
F1.p(1) = plot([a_sigsr:b_sigsr],M_avg_strain([a_sigsr:b_sigsr],mean([a_srm:b_srm]),mean([a_fct:b_fct])),'k');
xlabel('\sigma_{sr} [MPa]')
ylabel('\epsilon_{sm} [-]')

F1.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F1.p(2) = plot([a_srm:10:b_srm],M_avg_strain(mean([a_sigsr:b_sigsr]),[a_srm:10:b_srm],mean([a_fct:b_fct])),'k');
xlabel('s_{rm} [mm]')
ylabel('\epsilon_{sm} [-]')

F1.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F1.p(3) = plot([a_fct:0.1:b_fct],M_avg_strain(mean([a_sigsr:b_sigsr]),mean([a_srm:b_srm]),[a_fct:0.1:b_fct]),'k');
xlabel('\tau_{b1} [MPa]')
ylabel('\epsilon_{sm} [-]')

%% 2.3 Sampling the experimental design
% Latin hypercube sampling
warning('??')
p = 5;
M = 3;
P = factorial(M+p)/(factorial(M)*factorial(p));
k = 2; % oversampling rate
n = k*P;
u01 = lhsdesign(n,M);

% Transform -- Experimental design
x_sigsr = a_sigsr+u01(:,1)*(b_sigsr-a_sigsr);
x_srm = a_srm+u01(:,2)*(b_srm-a_srm);
x_taub1 = a_fct+u01(:,3)*(b_fct-a_fct);

% Evaluation of the full model on the experimental design
y_ED = M_avg_strain(x_sigsr,x_srm,x_taub1);
y_ED_mean = mean(y_ED);
y_ED_std = std(y_ED);
y_ED_var = y_ED_std^2;

