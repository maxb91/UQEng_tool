% To be called within Mini_project_MAIN.m after section 2.2

F4b.fig = figure('Units','Centimeters','Position',[2,2,(21-5)/2,7],...
    'Name','Comparison increasing degree');


F4b.ax=axes;
grid on
box on
hold on

coltemp = bone(length(incdeg)+2);
cols = coltemp(end-1:-1:2,:);
leg = cell(1,length(incdeg));
plot([0 0.035],[0 0],'k')
for i = 1:length(incdeg)
    F4b.p(i)=plot(y_ED_val,abs(incdeg(i).y_PCE_val./y_ED_val-1),'.','color',cols(i,:));
    leg{i} = ['{\it p}',' = ',num2str(i)];
end
xlim([0,0.035])
%set(gca,'yscale','log')

legend(F4b.p,leg,'Location','north')  
ylabel('PCE Response')
xlabel('Model Response')
hold off


