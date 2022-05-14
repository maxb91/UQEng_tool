% To be called within Mini_project_MAIN.m after section 2.2

F5.fig = figure('Units','Centimeters','Position',[2,2,21-5,7],...
    'Name','Comparison increasing degree');


F5.ax(1)=subplot(1,2,1);
grid on
box on
hold on

plot(polynomials,moment_mean,'k-o')
plot(polynomials,moment_variance,'-s','color',[0.7 0.7 0.7])
legend('Mean value','Variance')
xlabel('polynomial degree {\it p}')
ylabel('Relative error')
set(gca,'yscale','log')


F5.ax(2)=subplot(1,2,2);
grid on
box on
hold on

plot([incn.n],rel_error_n_mu,'k-o')
plot([incn.n],rel_error_n_var,'-s','color',[0.7 0.7 0.7])
legend('Mean value','Variance')
xlabel('Number of samples in the experimental design {\it n}')
ylabel('Relative error')
set(gca,'yscale','log')


% set limits
ylim1 = get(F5.ax(1),'YLim');
ylim2 = get(F5.ax(2),'YLim');
ylimits = [min([ylim1(1),ylim2(1)]), max([ylim1(2),ylim2(2)])];
set(F5.ax(1),'YLim',ylimits);
set(F5.ax(2),'YLim',ylimits);
