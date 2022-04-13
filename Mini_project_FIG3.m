% To be called within Mini_project_MAIN.m after section 2.3
%
% Experimental design from Monte Carlo sampling
% As a comparison with LHS

F3.fig = figure('Units','Centimeters','Position',[2,2,21-2,5.2],...
    'Name','Experimental design - Monte Carlo sampling');

F3.ax(1) = subplot(1,3,1);
grid on
box on
hold on
F3.p(1) = plot(xMC_sigsr,xMC_srm,'ko',...
    'Markersize',2);
xlabel('\sigma_{sr} [MPa]')
xlim([a_sigsr,b_sigsr])
ylabel('s_{rm} [mm]')
ylim([a_srm,b_srm]);

F3.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F3.p(2) = plot(xMC_srm,xMC_taub1,'ko',...
    'Markersize',2);
xlabel('s_{rm} [mm]')
xlim([a_srm,b_srm]);
ylabel('\tau_{b1} [MPa]')
ylim([a_fct,b_fct]);


F3.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F3.p(2) = plot(xMC_taub1,xMC_sigsr,'ko',...
    'Markersize',2);
xlabel('\tau_{b1} [MPa]')
xlim([a_fct,b_fct]);
ylabel('\sigma_{sr} [MPa]')
ylim([a_sigsr,b_sigsr])