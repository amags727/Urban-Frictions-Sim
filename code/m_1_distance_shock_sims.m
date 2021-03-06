%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   m1_SimulateModel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries

cd('/Volumes/GoogleDrive/My Drive/Github/Urban-Frictions-Sim/code');
addpath(genpath('depends'));
clear
rng(12345);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
alpha1 = 0.7; %firm cobb douglas exponent
beta1  = 0.7; %consumer cobb douglas exponent 
theta  = 3;  %frechet scale parameter
kappa  = 0.01; %commuting cost scalar
L_bar  = 10000; %city population 

 % right now set agglomeration forces to zero 
mu_a = 0; 
delta_a = 0;
mu_u_set = 0; %[0.2,0.3,0.4]
delta_u_set = 0; %[0.005,0.01,0.02]
distance_shock=[.9,.91,.92,.93,.94,.95,.96,.97,.98,.99];
xi_set=[0,.25,.5,.75, 1];
%%
% Monocentric city of NxN locations 
%%
n = 20; 
fun = (1:n); 

%create 400 ordered pairs 
    %create 1-20 20 times column vec
    coord_x = repmat(fun,1,n); 
    coord_x=coord_x'; %flip into column matrix 

    %create 1 20 times 2 20 times... column vec
    coord_y = repmat(fun,n,1); %1 20 times, 2 20 times...
    coord_y = coord_y(:); 

    %match the pairs
    A = [coord_y,coord_x]; 


% Defining distance between nodes 
D = pdist(A); 
Z = squareform(D);
t = 60*Z; %We set the distance between (i,j) and (i+1,j) equal to 60
d = exp(kappa*t); 

% Find central location (min distance to all other locations) 
tot_dist = sum(t,2);
[M,L] = min(tot_dist);
loc = L;

% Characteristics
I = n*n;
K = ones(I,1); % a columb vector of n^2 1s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Solve Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: procedure is: (i) solve model in GE to get w, A, (ii) feed in shock to amenities, 
%       (iii) solve for residential outcomes holding w,A fixed at initial values

%%
% Setup
%%

% Basic Settings
tol = 1e-7;
psi1 = 0.3;
N_sim = 1;
u_shocks = [1.1 1.2 1.5];
N_shocks = size(u_shocks,2);
mu_u_size = size(mu_u_set,2);
delta_u_size = size(delta_u_set,2);

% Matrix of shocks, one per simulation
u_bar_mat = ones(I,N_sim);%lognrnd(0,1,I,N_sim); 
A_bar_mat = ones(I,N_sim);%lognrnd(0,1,I,N_sim);
H_f_mat   = ones(I,N_sim);%lognrnd(0,1,I,N_sim);
H_r_mat   = ones(I,N_sim);%lognrnd(0,1,I,N_sim);

% Init values and arrays to store results
[resultsPE_L_r,resultsPE_rr,resultsPE_u] = deal(zeros(I,N_sim,N_shocks));
[resultsGE_L_r,resultsGE_rr,resultsGE_u,resultsGE_w_bar] = deal(zeros(I,N_sim,N_shocks));
[resultsInit_L_r,resultsInit_rr,resultsInit_u,resultsInit_w_bar] = deal(zeros(I,N_sim));

% Cells to store results for diff values of mu_u and delta_u
Init = cell(mu_u_size,delta_u_size);
PE   = cell(mu_u_size,delta_u_size);
GE   = cell(mu_u_size,delta_u_size);

%%
% Run simulations testing without frictions
        mu_u = mu_u_set(1,1);
        delta_u = delta_u_set(1,1);
        xi = 1;
        param    = v2struct(alpha1,beta1,theta,kappa,L_bar,I,mu_a,mu_u,delta_u,delta_a,xi);
        data     = v2struct(d,t,K);
        settings = v2struct(tol,psi1);
        
 
            %% (i) Solve for initial equilibrium 
            [w_init,rr_init,rf_init,A_init,u_init]   = deal(ones(I,1));
            U_bar_init = 1;
            open_city= 0;
            u_bar    = u_bar_mat(:,1);
            A_bar    = A_bar_mat(:,1);
            H_f      = H_f_mat(:,1);
            H_r      = H_r_mat(:,1);
            data     = struct('d',d,'t',t,'K',K,'u_bar',u_bar,'A_bar',A_bar,'H_f',H_f,'H_r',H_r);
            settings = struct('tol',tol,'psi1',psi1,'w_init',w_init,'rr_init',rr_init,'rf_init',rf_init,'A_init',A_init,'u_init',u_init,'U_bar_init',U_bar_init,'open_city',open_city);
            results_init_eqbm  = solveModel_GE(param,data,settings);
            save('../output/results_init.mat', 'results_init_eqbm');
            
    %% Test the effect of changing distances when there are no frictions  
    sample_size= 500; 
    no_friction_results= ones(sample_size,2 + length(distance_shock));
    %set the starting values based on the initial equilib. 
            open_city  = 0;
            U_bar_init = results_init_eqbm.U_bar;
            % Init values are from init eqbm
            w_init  = results_init_eqbm.w;
            rr_init = results_init_eqbm.rr;
            rf_init = results_init_eqbm.rf;
            A_init  = results_init_eqbm.A;
            u_init  = results_init_eqbm.u;
            settings = struct('tol',tol,'psi1',psi1,'w_init',w_init,'rr_init',rr_init,'rf_init',rf_init,'A_init',A_init,'u_init',u_init,'U_bar_init',U_bar_init,'open_city',open_city);
  
 
    sampled_nodes = transpose(randi([1, I],2,sample_size));
    for i =1:length(distance_shock)
        for j=1:sample_size
                test_t= t;
                test_t(sampled_nodes(j,1),sampled_nodes(j,2))=t(sampled_nodes(j,1),sampled_nodes(j,2))* distance_shock(i);
                test_t(sampled_nodes(j,2),sampled_nodes(j,1))=test_t(sampled_nodes(j,1),sampled_nodes(j,2));
                
               
                test_d=d;
                test_d(sampled_nodes(j,1),sampled_nodes(j,2))= exp(kappa*test_t(sampled_nodes(j,1),sampled_nodes(j,2)));
                test_d(sampled_nodes(j,2),sampled_nodes(j,1))= exp(kappa*test_t(sampled_nodes(j,1),sampled_nodes(j,2)));
                data     = struct('d',test_d,'t',test_t,'K',K,'u_bar',u_bar,'A_bar',A_bar,'H_f',H_f,'H_r',H_r);
                sim_eq    = solveModel_GE(param,data,settings);
                no_friction_results(j, 1)= log(t(sampled_nodes(j,1),sampled_nodes(j,2)));
                no_friction_results(j, 2)= log(results_init_eqbm.pi_ij(sampled_nodes(j,1),sampled_nodes(j,2)));
                if mod( j , sample_size/10 )==0
                    display("no frictions");
                    display(i);
                    display(j);
                    
               
                end 
                no_friction_results(j,i+2)= log(sim_eq.U_bar)- log(results_init_eqbm.U_bar);
        end
    end
    
save('../output/no_friction_results.mat', 'no_friction_results');
varnames= ["og_distance" "og_popularity" "d_shock_" + string(distance_shock)]
writematrix([varnames; no_friction_results],'../output/no_friction_results.csv');

%% Test the effect of changing distances when there are frictions  
    sample_size= 500; 
    friction_results= ones(sample_size,2 + length(xi_set));
    %set the starting values based on the initial equilib. 
            open_city  = 0;
            U_bar_init = results_init_eqbm.U_bar;
            % Init values are from init eqbm
            w_init  = results_init_eqbm.w;
            rr_init = results_init_eqbm.rr;
            rf_init = results_init_eqbm.rf;
            A_init  = results_init_eqbm.A;
            u_init  = results_init_eqbm.u;
            settings = struct('tol',tol,'psi1',psi1,'w_init',w_init,'rr_init',rr_init,'rf_init',rf_init,'A_init',A_init,'u_init',u_init,'U_bar_init',U_bar_init,'open_city',open_city);
  
 
    sampled_nodes = transpose(randi([1, I],2,sample_size));
    for i =1:length(xi_set)
        xi = xi_set(i);
        param  = v2struct(alpha1,beta1,theta,kappa,L_bar,I,mu_a,mu_u,delta_u,delta_a,xi);
        data     = struct('d',d,'t',t,'K',K,'u_bar',u_bar,'A_bar',A_bar,'H_f',H_f,'H_r',H_r);   
        results_init_eqbm  = solveModel_GE(param,data,settings);
        for j=1:sample_size
                test_t= t;
                test_t(sampled_nodes(j,1),sampled_nodes(j,2))=t(sampled_nodes(j,1),sampled_nodes(j,2))* distance_shock(1); %distance_shock is now set at 10 percent
                test_t(sampled_nodes(j,2),sampled_nodes(j,1))=test_t(sampled_nodes(j,1),sampled_nodes(j,2));
                test_d=d;
                test_d(sampled_nodes(j,1),sampled_nodes(j,2))= exp(kappa*test_t(sampled_nodes(j,1),sampled_nodes(j,2)));
                test_d(sampled_nodes(j,2),sampled_nodes(j,1))= exp(kappa*test_t(sampled_nodes(j,1),sampled_nodes(j,2)));
                data     = struct('d',test_d,'t',test_t,'K',K,'u_bar',u_bar,'A_bar',A_bar,'H_f',H_f,'H_r',H_r);
               
                %set the friction 
                sim_eq    = solveModel_GE(param,data,settings);
                friction_results(j, 1)= log(t(sampled_nodes(j,1),sampled_nodes(j,2)));
                friction_results(j, 2)= log(results_init_eqbm.pi_ij(sampled_nodes(j,1),sampled_nodes(j,2)));
                if mod( j , sample_size/10 )==0
                    display("frictions");
                    display(i);
                    display(j);
                end 
                friction_results(j,i+2)= log(sim_eq.U_bar)- log(results_init_eqbm.U_bar);
        end
    end
    
save('../output/friction_results.mat', 'friction_results');
varnames= ["og_distance" "og_popularity" "xi_level_" + string(xi_set)];
writematrix([varnames;friction_results],'../output/friction_results.csv');
