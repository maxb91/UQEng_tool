% To be called within Mini_project_MAIN.m after section 2.2

F1.fig = figure('Units','Centimeters','Position',[2,2,21-5,6],...
    'Name','Probabilistic input model');

F1.ax(1) = subplot(1,3,1);
grid on
box on
hold on
F1.p(1) = plot([sigsr_lim(1):sigsr_lim(2)],M_avg_strain([sigsr_lim(1):sigsr_lim(2)],mean([srm_lim(1):srm_lim(2)]),mean([fct_lim(1):fct_lim(2)])),'k');
xlabel('\sigma_{sr} [MPa]')
ylabel('\epsilon_{sm} [-]')

F1.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F1.p(2) = plot([srm_lim(1):10:srm_lim(2)],M_avg_strain(mean([sigsr_lim(1):sigsr_lim(2)]),[srm_lim(1):10:srm_lim(2)],mean([fct_lim(1):fct_lim(2)])),'k');
xlabel('s_{rm} [mm]')
ylabel('\epsilon_{sm} [-]')

F1.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F1.p(3) = plot([fct_lim(1):0.1:fct_lim(2)],M_avg_strain(mean([sigsr_lim(1):sigsr_lim(2)]),mean([srm_lim(1):srm_lim(2)]),[fct_lim(1):0.1:fct_lim(2)]),'k');
xlabel('\tau_{b1} [MPa]')
ylabel('\epsilon_{sm} [-]')