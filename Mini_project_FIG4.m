% To be called within Mini_project_MAIN.m after section 2.2

F4.fig = figure('Units','Centimeters','Position',[2,2,(21-5)/2,7],...
    'Name','Comparison increasing degree');


F4.ax=axes;
grid on
box on
hold on

coltemp = bone(length(incdeg)+2);
cols = coltemp(end-1:-1:2,:);
leg = cell(1,length(incdeg));
plot([0 0.035],[0 0.035],'k')
for i = 1:length(incdeg)
    F4.p(i)=plot(y_ED_val,incdeg(i).y_PCE_val,'.','color',cols(i,:));
    leg{i} = ['{\it p}',' = ',num2str(i)];
end


legend(F4.p,leg,'Location','Southeast')  
ylabel('PCE Response')
xlabel('Model Response')
hold off

axis equal


