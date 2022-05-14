% To be called within Mini_project_MAIN.m after section 2.2

F4c.fig = figure('Units','Centimeters','Position',[2,2,21-5,7],...
    'Name','Comparison increasing degree');


F4c.ax(1)=subplot(1,2,1);
grid on
box on
hold on

coltemp = bone(length(incdeg)+2);
cols = coltemp(end-1:-1:2,:);
leg = cell(1,length(incdeg));
plot([0 0.035],[0 0.035],'k')
for i = 1:length(incdeg)
    F4c.p(i)=plot(y_ED_val,incdeg(i).y_PCE_val,'.','color',cols(i,:));
    leg{i} = ['{\it p}',' = ',num2str(i)];
end


legend(F4c.p,leg,'Location','Southeast')  
ylabel('PCE Response {\it y_{PCE}}')
xlabel('Model Response {\it y_M}')
hold off

axis equal
xlim([0,0.035])


F4c.ax(2)=subplot(1,2,2);
grid on
box on
hold on

coltemp = bone(length(incdeg)+2);
cols = coltemp(end-1:-1:2,:);
leg2 = cell(1,length(incdeg));
plot([0 0.035],[0 0],'k')
for i = 1:length(incdeg)
    F4c.p2(i)=plot(y_ED_val,abs(incdeg(i).y_PCE_val./y_ED_val-1),'.','color',cols(i,:));
    leg2{i} = ['{\it p}',' = ',num2str(i)];
end
xlim([0,0.035])
%set(gca,'yscale','log')

legend(F4c.p2,leg2,'Location','north')  
ylabel('|{\it y_{PCE}} / {\it y_M} - 1|')
xlabel('Model Response {\it y_M}')
hold off


