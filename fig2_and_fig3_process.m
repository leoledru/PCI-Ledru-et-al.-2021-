%% FIG 2 PANEL A
clear all; clc;

% Number of runs and pre-allocation of variables
n_loop = 4;
out_simulation = cell(n_loop,1); 
out_symbiont = cell(n_loop,1);
out_symbiont_d = cell(n_loop,1);
out_time_speciation = zeros(n_loop,1);

% Host Germination and Competition
gamma = .2;
% gamma = 1;
host_germination = 1;
host_density_dependent = 0; % 0 = global, 1 = local
host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;

% Symbiont Colonization and Competition
symbiont_colonization = 1;
symbiont_density_dependent = 0; % 0 = global, 1 = local
coeff_compet = 0;
gamma_symbiont = .2;

% Mortality and Dispersal Cost
host_mortality = .06;
symbiont_mortality = .06;
coeff_disp = 0;
% coeff_disp = .45;
coeff_disp_h = 0;
% coeff_disp_h = .45;

% Initial conditions
load('data/cheater_eq.mat') 
initial_condis = total_new;

% Parallel loops on the number of iterations 
parfor i = 1:n_loop
    % Mutation Settings
    p_mut = 0.02; % probability of mutation
    muta_max = .5; % maximum mutation magnitude
    exp_mu = .25; % mu of exponential law
    muta_min = 0;
    muta_alpha = 1; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 1;
    muta_epsilon_symbiont = 1;
    
    % General Settings
    space_size = 100;
    tmax = 100000;
    % Neighbourhood 8 or 24
    neighbourhood = 8;
    % pre-allocation of host and symbiont matrices
    host = zeros(space_size^2,1);
    symbiont = zeros(space_size^2,1);
    host_d = zeros(space_size^2,1);
    symbiont_d = zeros(space_size^2,1);
    % application of the initial conditions
    total = initial_condis;
    h = sum(total>0,'all');
    m = sum(total==3,'all');
    initial_alpha = 0;
    h_alpha = ones(h,1).*initial_alpha;
    m_alpha = ones(m,1).*initial_alpha;
    host(total>0) = h_alpha;
    symbiont(total==3) = m_alpha;
    initial_epsilon = 1;
    h_epsilon = ones(h,1).*initial_epsilon;
    m_epsilon = ones(m,1).*initial_epsilon;
    host_d(total>0) = h_epsilon;
    symbiont_d(total==3) = m_epsilon;

    synchrone = ecoevo_model_time_to_speciation(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
             neighbourhood, host_mortality, symbiont_mortality,...
             host_density_dependent, host_germination, host_competition,...
             symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
             muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)

    out_time_speciation(i) = synchrone.time_speciation;
end

save(['time_to_speciation','.mat'],'out_time_speciation','-v7.3')
% save(['time_to_speciation_gamma02_dispcost45','.mat'],'out_time_speciation','-v7.3')
% save(['time_to_speciation_gamma1','mat'],'out_time_speciation','-v7.3')
% save(['time_to_speciation_gamma1_dispcost45','.mat'],'out_time_speciation','-v7.3')

%% FIG 2 PANEL B & FIG 3
clear all; clc;

n_loop = 1;
out_host = cell(n_loop,1);
out_host_d = cell(n_loop,1);
out_symbiont = cell(n_loop,1);
out_symbiont_d = cell(n_loop,1);
out_host_mean_trait = cell(n_loop,1);
out_symbiont_mean_trait = cell(n_loop,1);
out_time_speciation = zeros(n_loop,1);
out_symbiont_density = cell(n_loop,1);
out_host_density = cell(n_loop,1);
out_symbiont_mutualist_mean_trait = cell(n_loop,1);
out_symbiont_parasite_mean_trait = cell(n_loop,1);
out_parasite_density = cell(n_loop,1);
out_mutualiste_density = cell(n_loop,1);

% Host Germination and Competition
gamma = .2;
host_germination = 1;
host_density_dependent = 0; % 0 = global, 1 = local
host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;

% Symbiont Colonization and Competition
symbiont_colonization = 1;
symbiont_density_dependent = 0; % 0 = global, 1 = local
coeff_compet = 0;
gamma_symbiont = .2;

% Mortality and Dispersal Cost
host_mortality = .06;
symbiont_mortality = .06;
coeff_disp = 0;
coeff_disp_h = 0;

% Initial conditions
load('data/cheater_eq.mat') 
initial_condis = total_new;

for i = 1:n_loop
  % Mutation Settings
    p_mut = 0.02; % probability of mutation
    muta_max = .5; % maximum mutation magnitude
    exp_mu = .25; % mu of exponential law
    muta_min = 0;
    muta_alpha = 1; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 1;
    muta_epsilon_symbiont = 1;
    
    % General Settings
    space_size = 100;
    tmax = 10000;
    % Neighbourhood 8 or 24
    neighbourhood = 8;
    % pre-allocation of host and symbiont matrices
    host = zeros(space_size^2,1);
    symbiont = zeros(space_size^2,1);
    host_d = zeros(space_size^2,1);
    symbiont_d = zeros(space_size^2,1);
    % application of the initial conditions
    total = initial_condis;
    h = sum(total>0,'all');
    m = sum(total==3,'all');
    initial_alpha = 0;
    h_alpha = ones(h,1).*initial_alpha;
    m_alpha = ones(m,1).*initial_alpha;
    host(total>0) = h_alpha;
    symbiont(total==3) = m_alpha;
    initial_epsilon = 1;
    h_epsilon = ones(h,1).*initial_epsilon;
    m_epsilon = ones(m,1).*initial_epsilon;
    host_d(total>0) = h_epsilon;
    symbiont_d(total==3) = m_epsilon;
    
    synchrone = ecoevo_model(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                 neighbourhood, host_mortality, symbiont_mortality,...
                 host_density_dependent, host_germination, host_competition,...
                 symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)

    speciation_test = synchrone.speciation;
    out_symbiont_density{i} = synchrone.symbiont_density;
    out_parasite_density{i} = synchrone.parasite_density;
    out_mutualiste_density{i} = synchrone.mutualiste_density;
    out_host_density{i} = synchrone.host_density;
    out_host_mean_trait{i} = synchrone.host_mean_trait;
    out_symbiont_mean_trait{i} = synchrone.symbiont_mean_trait;
    out_symbiont_mutualist_mean_trait{i} = synchrone.symbiont_mutualist_mean_trait;
    out_symbiont_parasite_mean_trait{i} = synchrone.symbiont_parasite_mean_trait;
    out_time_speciation(i) = synchrone.time_speciation;
    out_host{i} = synchrone.host;
    out_host_d{i} = synchrone.host_d;
    out_symbiont{i} = synchrone.symbiont;
    out_symbiont_d{i} = synchrone.symbiont_d;
end

save(['simulations_gamma02_h','.mat'],'out_host_density','out_host_mean_trait',...
    'out_mutualiste_density','out_parasite_density','out_symbiont_mean_trait',...
    'out_symbiont_mutualist_mean_trait','out_symbiont_parasite_mean_trait',...
     'out_time_speciation','out_host','out_host_d','out_symbiont','out_symbiont_d','-v7.3')
