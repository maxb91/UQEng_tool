%% Coefficients
% Could be used together with polynomial degree variation plot to show from
% which degree on the model works well?
pol_degree = sum(incdeg(5).p_index,2);
markers = {'s','o','d','^','v'};

coltemp = bone(incdeg(5).p);
cols = coltemp(end:-1:1,:);

figure('Units','centimeters','Position',[1,1,21-5,7])
hold on
grid on
box on
for i=1:incdeg(5).p
    coeff_indices = find(pol_degree == incdeg(5).p+1-i);
    SC(incdeg(5).p+1-i)=scatter(coeff_indices,...
        abs(incdeg(5).c(coeff_indices)/incdeg(5).c(1)),...
        [],'MarkerEdgeColor',[0 0 0],'markerfacecolor',cols(i,:),...
        'marker',markers{i});
end
set(gca,'YScale','log')
legend(SC,{'First order','Second order','Third order','Fourth order', 'Fifth order'})

%% Run the whole thing multiple times for convergence plot
% TAKES VERY LONG TO RUN!!
% Mini_project_CONV
Mini_project_FIG5b
Mini_project_FIG5c


%% Computation times

% e.g. PCE compared to MC? 
% Dependent on number of MC samples. 
% Plot e.g. t on x-axis, accuracy on y-axis with different markers
% How to quantify accuracy?


%% In UQLab: OLS vs LARS? 
% with coefficient plot

%% In UQLab: Truncation?