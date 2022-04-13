% To be called within Mini_project_MAIN.m after section 2.3
%
% Experimental design from Latin hypercube sampling

F2.fig = figure('Units','Centimeters','Position',[2,2,21-2,5.2],...
    'Name','Experimental design - Latin hypercube sampling');

F2.ax(1) = subplot(1,3,1);
grid on
box on
hold on
F2.p(1) = plot(x_sigsr,x_srm,'ko',...
    'Markersize',2);
xlabel('\sigma_{sr} [MPa]')
xlim([a_sigsr,b_sigsr])
ylabel('s_{rm} [mm]')
ylim([a_srm,b_srm]);

F2.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F2.p(2) = plot(x_srm,x_taub1,'ko',...
    'Markersize',2);
xlabel('s_{rm} [mm]')
xlim([a_srm,b_srm]);
ylabel('\tau_{b1} [MPa]')
ylim([a_fct,b_fct]);


F2.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F2.p(2) = plot(x_taub1,x_sigsr,'ko',...
    'Markersize',2);
xlabel('\tau_{b1} [MPa]')
xlim([a_fct,b_fct]);
ylabel('\sigma_{sr} [MPa]')
ylim([a_sigsr,b_sigsr])