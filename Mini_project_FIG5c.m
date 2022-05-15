% To be called within Mini_project_MAIN.m after section 2.6

F5c.fig = figure('Units','Centimeters','Position',[2,2,21-5,14],...
    'Name','Comparison increasing degree');

xlimp = [1 5];
xlimn = [0 150];
ylimmu = [0.011 0.012];
ylimvar = [1 2.5]*10^(-5);


F5c.ax(1)=subplot(2,2,1);
grid on
box on
hold on
F5c.s1 = scatter(degrees_conv,reshape([incdeg_conv.muhat_PCE],size(degrees_conv)),...
    9,'o','MarkerEdgeColor',[0.3 0.3 0.3]);
F5c.p(1) = plot(degrees_conv(1,:),mean(reshape([incdeg_conv.muhat_PCE],size(degrees_conv))),'k-o','markerfacecolor',[0.3 0.3 0.3]);
F5c.mum = plot(xlimp, ones(size(xlimp))*mu_M,...
    'color',[0.3 0.3 0.3],'linewidth',1);
legend([F5c.s1(1),F5c.p(1),F5c.mum],{'20 runs','average','\mu_M'})
xlabel('polynomial degree {\it p}')
ylabel('\mu')
xlim(xlimp)
ylim(ylimmu)


F5c.ax(2)=subplot(2,2,3);
grid on
box on
hold on

F5c.s2 = scatter(degrees_conv,reshape([incdeg_conv.varhat_PCE],size(degrees_conv)),...
    13,'s','MarkerEdgeColor',[0.7 0.7 0.7]);
F5c.p(2) = plot(degrees_conv(1,:),mean(reshape([incdeg_conv.varhat_PCE],size(degrees_conv))),'k-s','markerfacecolor',[0.7 0.7 0.7],'markersize',7);
F5c.varm = plot(xlimp, ones(size(xlimp))*var_M,...
    'color',[0.7 0.7 0.7],'linewidth',1);
legend([F5c.s2(1),F5c.p(2),F5c.varm],{'20 runs','average','\sigma^2_M'})
xlabel('polynomial degree {\it p}')
ylabel('\sigma^2')
xlim(xlimp)
ylim(ylimvar)

%%
F5c.ax(3)=subplot(2,2,2);
grid on
box on
hold on

F5c.s3 = scatter(oversampling_conv(:,2:end)*35,reshape([incn_conv(:,2:end).muhat_PCE],size(oversampling_conv(:,2:end))),...
    9,'o','MarkerEdgeColor',[0.3 0.3 0.3]);
F5c.p2(1) = plot(oversampling_conv(1,2:end)*35,mean(reshape([incn_conv(:,2:end).muhat_PCE],size(oversampling_conv(:,2:end)))),'k-o','markerfacecolor',[0.3 0.3 0.3]);
F5c.mum = plot(xlimn, ones(size(xlimn))*mu_M,...
    'color',[0.3 0.3 0.3],'linewidth',1);
legend([F5c.s3(1),F5c.p2(2),F5c.mum],{'20 runs','average','\mu_M'})
xlabel('Number of samples in the ED {\it n}')
ylabel('\mu')
xlim(xlimn)
ylim(ylimmu)

F5c.ax(4)=subplot(2,2,4);
grid on
box on
hold on

F5c.s4 = scatter(oversampling_conv(:,2:end)*35,reshape([incn_conv(:,2:end).varhat_PCE],size(oversampling_conv(:,2:end))),...
    13,'s','MarkerEdgeColor',[0.7 0.7 0.7]);
F5c.p2(2) = plot(oversampling_conv(1,2:end)*35,mean(reshape([incn_conv(:,2:end).varhat_PCE],size(oversampling_conv(:,2:end)))),'k-s','markerfacecolor',[0.7 0.7 0.7],'markersize',7);
F5c.varm = plot(xlimn, ones(size(xlimn))*var_M,...
    'color',[0.7 0.7 0.7],'linewidth',1);
legend([F5c.s4(1),F5c.p2(2),F5c.varm],{'20 runs','average','\sigma^2_M'})
xlabel('Number of samples in the ED {\it n}')
ylabel('\sigma^2')
xlim(xlimn)
ylim(ylimvar)

% pos = get(F5c.fig,'Position');
% set(F5c.fig,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[pos(3), pos(4)])
% print('fig-conv1','-dpdf','-vector')

