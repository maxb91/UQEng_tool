% To be called within Mini_project_MAIN.m after section 2.2

F1.fig = figure('Units','Centimeters','Position',[2,2,21-5,6],...
    'Name','Probabilistic input model');

sigsr_mean = sigmeas;
srm_mean = 0.667*sr0;
taub_mean = fct_mean;

F1.ax(1) = subplot(1,3,1);
grid on
box on
hold on
F1.p(1) = plot([sigsr_lim(1):sigsr_lim(2) sigsr_lim(2)],M_avg_strain([sigsr_lim(1):sigsr_lim(2) sigsr_lim(2)],srm_mean,taub_mean),'k');
xlabel('\sigma_{sr} [MPa]')
ylabel('\epsilon_{sm} [-]')
xlim([sigsr_lim(1),sigsr_lim(2)])
ylim([0.005,0.02])

F1.ax(2) = subplot(1,3,2);
grid on
box on
hold on
F1.p(2) = plot([srm_lim(1):10:srm_lim(2) srm_lim(2)],M_avg_strain(sigsr_mean,[srm_lim(1):10:srm_lim(2) srm_lim(2)],taub_mean),'k');
xlabel('s_{rm} [mm]')
ylabel('\epsilon_{sm} [-]')
xlim([srm_lim(1),srm_lim(2)])
ylim([0.005,0.02])

F1.ax(3) = subplot(1,3,3);
grid on
box on
hold on
F1.p(3) = plot([tau_lim(1):0.1:tau_lim(2) tau_lim(2)],M_avg_strain(sigsr_mean,srm_mean,[tau_lim(1):0.1:tau_lim(2) tau_lim(2)]),'k');
xlabel('\tau_{b1} [MPa]')
ylabel('\epsilon_{sm} [-]')
xlim([tau_lim(1),tau_lim(2)])
ylim([0.005,0.02])