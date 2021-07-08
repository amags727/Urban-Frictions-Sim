%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   m2_graphs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
if strcmp(getenv('USER'),'ntsivanidis')==0
    cd('D:/GitHub/slums-india/modelling/mumbai_modelsim/model_simulation/code');
else
    cd('/Users/ntsivanidis/Documents/GitHub/slums-india/modelling/mumbai_modelsim/model_simulation/code');
end
addpath(genpath('depends'));
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results
load('../output/Solve_Init.mat')
load('../output/Solve_PE.mat')
load('../output/Solve_GE.mat')

% Load parameters
load('../output/Parameters.mat');
struct2vars(Param)
mu_u_size = size(mu_u_set,2);
delta_u_size = size(delta_u_set,2);
dist_loc = t(:,loc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Calculate log changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars_GE = fieldnames(GE{1,1});
vars_PE = fieldnames(PE{1,1});
vars_init = fieldnames(Init{1,1});
var = {'Lr','rr','u', 'w_bar'};
var = matlab.lang.makeValidName(var);
dlogPE = cell(mu_u_size,delta_u_size);
dlogGE = cell(mu_u_size,delta_u_size);

for x=1:mu_u_size
    for y=1:delta_u_size
        for i = 1:length(vars_PE)
            dlogPE{x,y}.(var{i}) = log(mean(PE{x,y}.(vars_PE{i}),2))-log(mean(Init{x,y}.(vars_init{i}),2));
        end
    end
end

for x=1:mu_u_size
    for y=1:delta_u_size
        for i = 1:length(vars_GE)
            dlogGE{x,y}.(var{i}) = log(mean(GE{x,y}.(vars_GE{i}),2))-log(mean(Init{x,y}.(vars_init{i}),2));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: (i) compare shocks holds spillovers parameters contantant in x=2,
% y=2 and (ii) compare ammenities holds shock constant at s=2

%%
% Compare shocks
%%

x = 2;
y = 2; 

% Lr, rr, u
for i = 1:length(vars_PE)
    scatter(dist_loc,dlogPE{x,y}.(var{i})(:,:,1),'LineWidth',2)
    hold on 
    scatter(dist_loc,dlogPE{x,y}.(var{i})(:,:,2),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogPE{x,y}.(var{i})(:,:,3),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend('Shock 10%','Shock 20%','Shock 50%','Location','northeast')
    title(['\mu = ' num2str(mu_u_set(x)) ', \delta = ' num2str(delta_u_set(y))])
    grid on
    filename = ['../output/figures/PE/compare_shocks/' var{i} '_mu' num2str(mu_u_set(x)) '_delta' num2str(delta_u_set(y)) '.png'];
    print(filename, '-dpng');
    
    scatter(dist_loc,dlogGE{x,y}.(var{i})(:,:,1),'LineWidth',2);
    hold on 
    scatter(dist_loc,dlogGE{x,y}.(var{i})(:,:,2),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogGE{x,y}.(var{i})(:,:,3),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend('Shock 10%','Shock 20%','Shock 50%','Location','northeast')
    title(['\mu = ' num2str(mu_u_set(x)) ', \delta = ' num2str(delta_u_set(y))])
    grid on
    filename = ['../output/figures/GE/compare_shocks/' var{i} '_mu' num2str(mu_u_set(x)) '_delta' num2str(delta_u_set(y)) '.png'];
    print(filename, '-dpng');
end

% w_bar
scatter(dist_loc,dlogGE{x,y}.(var{4})(:,:,1),'LineWidth',2);
hold on 
scatter(dist_loc,dlogGE{x,y}.(var{4})(:,:,2),'LineWidth',2);
hold off
hold on
scatter(dist_loc,dlogGE{x,y}.(var{4})(:,:,3),'LineWidth',2);
hold off
xlabel('Distance to shock')
ylabel(['Log Change in ' var{4}])
legend('Shock 10%','Shock 20%','Shock 50%','Location','southeast')
title(['\mu = ' num2str(mu_u_set(x)) ', \delta = ' num2str(delta_u_set(y))])
grid on
filename = ['../output/figures/GE/compare_shocks/' var{4} '_mu' num2str(mu_u_set(x)) '_delta' num2str(delta_u_set(y)) '.png'];
print(filename, '-dpng');

%%
% Compare spillover
%%

% Different decay parameters
s=2;
x=2;

%Lr, rr, u
for i = 1:length(vars_PE)
    scatter(dist_loc,dlogPE{x,1}.(var{i})(:,:,s),'LineWidth',2)
    hold on 
    scatter(dist_loc,dlogPE{x,2}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogPE{x,3}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend(strcat('\delta = ', num2str(delta_u_set(1))),strcat('\delta = ', num2str(delta_u_set(2))), ...
        strcat('\delta = ', num2str(delta_u_set(3))), 'Location','northeast')
    title(['\mu = ' num2str(mu_u_set(x)) ', shock = ' num2str(u_shocks(s))])
    grid on
    filename = ['../output/figures/PE/compare_spillovers/' var{i} '_mu' num2str(mu_u_set(x)) '_shock' num2str(u_shocks(s)) '.png'];
    print(filename, '-dpng');
    
    scatter(dist_loc,dlogGE{x,1}.(var{i})(:,:,s),'LineWidth',2);
    hold on 
    scatter(dist_loc,dlogGE{x,2}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogGE{x,3}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend(strcat('\delta = ', num2str(delta_u_set(1))),strcat('\delta = ', num2str(delta_u_set(2))), ...
        strcat('\delta = ', num2str(delta_u_set(3))), 'Location','northeast')
    title(['\mu = ' num2str(mu_u_set(x)) ', shock = ' num2str(u_shocks(s))])
    grid on
    filename = ['../output/figures/GE/compare_spillovers/' var{i} '_mu' num2str(mu_u_set(x)) '_shock' num2str(u_shocks(s)) '.png'];
    print(filename, '-dpng');
end

% w_bar
scatter(dist_loc,dlogGE{x,1}.(var{4})(:,:,s),'LineWidth',2);
hold on 
scatter(dist_loc,dlogGE{x,2}.(var{4})(:,:,s),'LineWidth',2);
hold off
hold on
scatter(dist_loc,dlogGE{x,3}.(var{4})(:,:,s),'LineWidth',2);
hold off
xlabel('Distance to shock')
ylabel(['Log Change in ' var{4}])
legend(strcat('\delta = ', num2str(delta_u_set(1))),strcat('\delta = ', num2str(delta_u_set(2))), ...
    strcat('\delta = ', num2str(delta_u_set(3))), 'Location','southeast')
title(['\mu = ' num2str(mu_u_set(x)) ', shock = ' num2str(u_shocks(s))])
grid on
filename = ['../output/figures/GE/compare_spillovers/' var{4} '_mu' num2str(mu_u_set(x)) '_shock' num2str(u_shocks(s)) '.png'];
print(filename, '-dpng');

% Different spillovers 
s=2;
y=2;

%Lr, rr, u
for i = 1:length(vars_PE)
    scatter(dist_loc,dlogPE{1,y}.(var{i})(:,:,s),'LineWidth',2);
    hold on 
    scatter(dist_loc,dlogPE{2,y}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogPE{3,y}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend(strcat('\mu = ', num2str(mu_u_set(1))),strcat('\mu = ', num2str(mu_u_set(2))), ...
        strcat('\mu = ', num2str(mu_u_set(3))), 'Location','northeast')
    title(['\delta = ' num2str(delta_u_set(y)) ', shock = ' num2str(u_shocks(s))])
    grid on
    filename = ['../output/figures/PE/compare_spillovers/' var{i} '_delta' num2str(delta_u_set(y)) '_shock' num2str(u_shocks(s)) '.png'];
    print(filename, '-dpng');
    
    scatter(dist_loc,dlogGE{1,y}.(var{i})(:,:,s),'LineWidth',2)
    hold on 
    scatter(dist_loc,dlogGE{2,y}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    hold on
    scatter(dist_loc,dlogGE{3,y}.(var{i})(:,:,s),'LineWidth',2);
    hold off
    xlabel('Distance to shock')
    ylabel(['Log Change in ' var{i}])
    legend(strcat('\mu = ', num2str(mu_u_set(1))),strcat('\mu = ', num2str(mu_u_set(2))), ...
        strcat('\mu = ', num2str(mu_u_set(3))), 'Location','northeast')
    title(['\delta = ' num2str(delta_u_set(y)) ', shock = ' num2str(u_shocks(s))])
    grid on
    filename = ['../output/figures/GE/compare_spillovers/' var{i} '_delta' num2str(delta_u_set(y)) '_shock' num2str(u_shocks(s)) '.png'];
    print(filename, '-dpng');
end

% w_bar
scatter(dist_loc,dlogGE{1,y}.(var{4})(:,:,s),'LineWidth',2)
hold on 
scatter(dist_loc,dlogGE{2,y}.(var{4})(:,:,s),'LineWidth',2);
hold off
hold on
scatter(dist_loc,dlogGE{3,y}.(var{4})(:,:,s),'LineWidth',2);
hold off
xlabel('Distance to shock')
ylabel(['Log Change in ' var{4}])
legend(strcat('\mu = ', num2str(mu_u_set(1))),strcat('\mu = ', num2str(mu_u_set(2))), ...
    strcat('\mu = ', num2str(mu_u_set(3))), 'Location','southeast')
title(['\delta = ' num2str(delta_u_set(y)) ', shock = ' num2str(u_shocks(s))])
grid on
filename = ['../output/figures/GE/compare_spillovers/' var{4} '_delta' num2str(delta_u_set(y)) '_shock' num2str(u_shocks(s)) '.png'];
print(filename, '-dpng');
