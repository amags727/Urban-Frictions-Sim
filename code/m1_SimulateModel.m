%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   m1_SimulateModel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
cd('/Volumes/GoogleDrive/My Drive/ grad school/Research/RA work/Gabriel /Nick stuff /old_model_simulation/code');
addpath(genpath('depends'));
clear
rng(12345);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
alpha1 = 0.7;
beta1  = 0.7;
theta  = 3;
kappa  = 0.01;
L_bar  = 10000;

mu_a = 0;
delta_a = 0;
mu_u_set = [0.2,0.3,0.4];
delta_u_set = [0.005,0.01,0.02];

%%
% Monocentric city of NxN locations 
%%
n = 20; 
fun = (1:n);
coord_x = repmat(fun,1,n);
coord_x=coord_x';
coord_y = repmat(fun,n,1);
coord_y = coord_y(:);
A = [coord_y,coord_x];

% Distance from i to j is the euclidean distance
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
K = ones(I,1);

% Understand how the decay fcuntion depends on parameters
explore_param = 0;
if explore_param==1
    mu_u = 0.3;
    delta_u=0.01;
    u_bar = ones(I,1);
    D_r_init = u_bar;
    u_shares = repmat(D_r_init',I,1).*exp(-delta_u*t);
    u_shares = u_shares./repmat(sum(u_shares,2),1,I);
    u_bar_init = u_bar;
    u_bar(loc,:)   = u_bar(loc,:).*1.25;
    u_hat =  (sum(u_shares.*repmat((u_bar./u_bar_init)',I,1),2)).^mu_u;
    scatter(t(:,loc),u_hat)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Solve Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: procedure is: (i) solve model in GE to get w, A, (ii) feed in shock to amenities, 
%       (iii) solve for residential outcomes holding w,A fixed at initial values

%%
% Setup
%%

% Basic Settings
tol = 1e-1;
psi1 = 0.3;
N_sim = 1;
u_shocks = [1.1 1.2 1.5];
N_shocks = size(u_shocks,2);
mu_u_size = size(mu_u_set,2);
delta_u_size = size(delta_u_set,2);

% Matrix of shcoks, one per simulation
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
% Run simulations
%%

% Loop over different value of delta_u and mu_u
for x=1:1
    for y=1:1
        mu_u = mu_u_set(1,x);
        delta_u = delta_u_set(1,y);
        param    = v2struct(alpha1,beta1,theta,kappa,L_bar,I,mu_a,mu_u,delta_u,delta_a);
        data     = v2struct(d,t,K);
        settings = v2struct(tol,psi1);
        
        parfor s = 1:N_sim
            %% (i) Solve for initial equilibrium 
            [w_init,rr_init,rf_init,A_init,u_init]   = deal(ones(I,1));
            U_bar_init = 1;
            open_city= 0;
            u_bar    = u_bar_mat(:,s);
            A_bar    = A_bar_mat(:,s);
            H_f      = H_f_mat(:,s);
            H_r      = H_r_mat(:,s);
            data     = struct('d',d,'t',t,'K',K,'u_bar',u_bar,'A_bar',A_bar,'H_f',H_f,'H_r',H_r);
            settings = struct('tol',tol,'psi1',psi1,'w_init',w_init,'rr_init',rr_init,'rf_init',rf_init,'A_init',A_init,'u_init',u_init,'U_bar_init',U_bar_init,'open_city',open_city);
            tic
            results_init_eqbm      = solveModel_GE(param,data,settings);
            toc
        end 
    end 
end 
 