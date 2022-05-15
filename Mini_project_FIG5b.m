% To be called within Mini_project_MAIN.m after section 2.6

F5b.fig = figure('Units','Centimeters','Position',[2,2,21-5,7],...
    'Name','Comparison increasing degree');


F5b.ax(1)=subplot(1,2,1);
grid on
box on
hold on
F5b.s1 = scatter(degrees_conv,abs(reshape([incdeg_conv.muhat_PCE],size(degrees_conv))/mu_M-1),...
    9,'o','MarkerEdgeColor',[0.3 0.3 0.3]);
F5b.s2 = scatter(degrees_conv,abs(reshape([incdeg_conv.varhat_PCE],size(degrees_conv))/var_M-1),...
    13,'s','MarkerEdgeColor',[0.7 0.7 0.7]);
F5b.p(1) = plot(degrees_conv(1,:),abs(mean(reshape([incdeg_conv.muhat_PCE],size(degrees_conv))/mu_M-1)),'k-o','markerfacecolor',[0.3 0.3 0.3]);
F5b.p(2) = plot(degrees_conv(1,:),abs(mean(reshape([incdeg_conv.varhat_PCE],size(degrees_conv))/var_M-1)),'k-s','markerfacecolor',[0.7 0.7 0.7],'markersize',7);
legend(F5b.p,{'Mean value','Variance'})
xlabel('polynomial degree {\it p}')
ylabel('\Delta')
set(gca,'yscale','log')

%%
F5b.ax(2)=subplot(1,2,2);
grid on
box on
hold on

F5b.s3 = scatter(oversampling_conv(:,2:end)*35,abs(reshape([incn_conv(:,2:end).muhat_PCE],size(oversampling_conv(:,2:end)))/mu_M-1),...
    9,'o','MarkerEdgeColor',[0.3 0.3 0.3]);
F5b.s4 = scatter(oversampling_conv(:,2:end)*35,abs(reshape([incn_conv(:,2:end).varhat_PCE],size(oversampling_conv(:,2:end)))/var_M-1),...
    13,'s','MarkerEdgeColor',[0.7 0.7 0.7]);
F5b.p2(1) = plot(oversampling_conv(1,2:end)*35,mean(reshape([incn_conv(:,2:end).muhat_PCE],size(oversampling_conv(:,2:end)))/mu_M-1),'k-o','markerfacecolor',[0.3 0.3 0.3]);
F5b.p2(2) = plot(oversampling_conv(1,2:end)*35,mean(reshape([incn_conv(:,2:end).varhat_PCE],size(oversampling_conv(:,2:end)))/var_M-1),'k-s','markerfacecolor',[0.7 0.7 0.7],'markersize',7);
legend(F5b.p2,{'Mean value','Variance'})
xlabel('Number of samples in the experimental design {\it n}')
ylabel('\Delta')
set(gca,'yscale','log')


% set limits
ylim1 = get(F5b.ax(1),'YLim');
ylim2 = get(F5b.ax(2),'YLim');
ylimits = [min([ylim1(1),ylim2(1)]), max([ylim1(2),ylim2(2)])];
set(F5b.ax(1),'YLim',ylimits);
set(F5b.ax(2),'YLim',ylimits);
