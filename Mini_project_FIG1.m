% To be called within Mini_project_MAIN.m after section 2.2

F1.fig = figure('Units','Centimeters','Position',[2,2,21-5,6],...
    'Name','Probabilistic input model');

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