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
xlim([sigsr_lim(1),sigsr_lim(2)])
ylabel('s_{rm} [mm]')
ylim([srm_lim(1),srm_lim(2)]);

F2.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F2.p(2) = plot(x_srm,x_taub1,'ko',...
    'Markersize',2);
xlabel('s_{rm} [mm]')
xlim([srm_lim(1),srm_lim(2)]);
ylabel('\tau_{b1} [MPa]')
ylim([tau_lim(1),tau_lim(2)]);


F2.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F2.p(2) = plot(x_taub1,x_sigsr,'ko',...
    'Markersize',2);
xlabel('\tau_{b1} [MPa]')
xlim([tau_lim(1),tau_lim(2)]);
ylabel('\sigma_{sr} [MPa]')
ylim([sigsr_lim(1),sigsr_lim(2)])