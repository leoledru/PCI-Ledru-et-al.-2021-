%% PANEL A
clear all; clc;

nb_linspace = 10;
nb_rep = 50;
host_density_dependent = 0; % 0 = global, 1 = local
symbiont_density_dependent = 0; % 0 = global, 1 = local
% Mortality and dispersal cost
host_mortality = .06;
symbiont_mortality = .06;
coeff_disp = .2;
coeff_disp_h = .2;
% Host Germination and Competition
gamma_loop = linspace(.1,2,nb_linspace);
host_germination = 1;
% Symbiont Colonization and Competition
symbiont_colonization = 1;
coeff_compet = 0;

params_name = {'host_germination','symbiont_colonization','coeff_compet'};
params_value = [host_germination,symbiont_colonization,coeff_compet];
for name_loop = 1:length(params_name)
    params{1,name_loop} = params_name{name_loop};
    params{2,name_loop} = params_value(name_loop);
end

% For save details
% out_simulation = cell(nb_linspace,nb_rep); 
% out_init = cell(nb_linspace,nb_rep); 
% out_host = cell(nb_linspace,nb_rep); 
% out_host_d = cell(nb_linspace,nb_rep); 
% out_symbiont = cell(nb_linspace,nb_rep);
% out_symbiont_d = cell(nb_linspace,nb_rep);
% out_persist = zeros(nb_linspace,nb_rep);
% out_time = zeros(nb_linspace,nb_rep);

load('data/cheater_eq.mat')
initial_condis = total_new;
for i = 1:nb_linspace
    space_size = 100;
    tmax = 10000;
    % Neighbourhood 8 or 24
    neighbourhood = 8;
    % Competition function
    gamma = gamma_loop(i);
    gamma_symbiont = gamma_loop(i);
    host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;
    % Mutation settings
    p_mut = 0.05; % probability of mutation
    muta_max = 0.5; % maximum mutation magnitude
    muta_min = 0;
    exp_mu = .25;
    muta_alpha = 1; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 1;
    muta_epsilon_symbiont = 1;
    mutualism_percent = zeros(1,nb_rep);
    parfor ii = 1:nb_rep
        host = zeros(space_size^2,1);
        symbiont = zeros(space_size^2,1);
        host_d = zeros(space_size^2,1);
        symbiont_d = zeros(space_size^2,1);
        total = initial_condis;
        h = sum(total>0,'all');
        s = sum(total==3,'all');
        % trait alpha
        alpha = 0;
        h_alpha = ones(h,1).*alpha;
        s_alpha = ones(s,1).*alpha;
        host(total>0) = h_alpha;
        symbiont(total==3) = s_alpha;
        % trait epsilon
        epsilon = 1;
        h_epsilon = ones(h,1).*epsilon;
        s_epsilon = ones(s,1).*epsilon;
        host_d(total>0) = h_epsilon;
        symbiont_d(total==3) = s_epsilon;
             
        synchrone = ecoevo_model(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                 neighbourhood, host_mortality, symbiont_mortality,...
                 host_density_dependent, host_germination, host_competition,...
                 symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)

%         For details
%         out_host{i,ii} = synchrone.host;
%         out_host_d{i,ii} = synchrone.host_d;
%         out_symbiont{i,ii} = synchrone.symbiont;
%         out_symbiont_d{i,ii} = synchrone.symbiont_d;
%         out_simulation{i,ii} = synchrone.simulation;
                              
        symbiont = synchrone.symbiont;
        simulation = synchrone.simulation;
        symbiont(simulation~=3) = NaN; % to not confound 0 because of empty cell with 0 because of a trait=0
        
        mutualism_percent(ii) = sum(symbiont>.475)/sum(~isnan(symbiont));
    end 
    Mutualism_percent(:,i) = mutualism_percent;
%     if mod(i,10)==0
%     save(['effect_of_competition','.mat'],'Mutualism_percent','host_density_dependent_loop','nb_linspace','nb_rep','params','-v7.3');
%     end
end
save(['effect_of_gamma','.mat'],'Mutualism_percent','gamma_loop','nb_linspace','nb_rep','params','-v7.3');

%% PANEL B
clear all; clc;

nb_linspace = 10;
nb_rep = 50;
host_density_dependent_loop = linspace(0,1,nb_linspace); % 0 = global, 1 = local
symbiont_density_dependent_loop = linspace(0,1,nb_linspace); % 0 = global, 1 = local

% Mortality and dispersal cost
host_mortality = .06;
symbiont_mortality = .06;
coeff_disp = .2;
coeff_disp_h = .2;
% Host Germination and Competition
gamma = .2;
host_germination = 1;
host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;
% Symbiont Colonization and Competition
symbiont_colonization = 1;
coeff_compet = 0;
gamma_symbiont = .2;

params_name = {'gamma','host_germination','symbiont_colonization','coeff_compet','gamma_symbiont'};
params_value = [gamma,host_germination,symbiont_colonization,coeff_compet,gamma_symbiont];
for name_loop = 1:length(params_name)
    params{1,name_loop} = params_name{name_loop};
    params{2,name_loop} = params_value(name_loop);
end

load('data/cheater_eq.mat')
initial_condis = total_new;
for i = 1:nb_linspace
    space_size = 100;
    tmax = 10000;
    % Neighbourhood 8 or 24
    neighbourhood = 8;
    % Competition scale
    host_density_dependent = host_density_dependent_loop(i);
    symbiont_density_dependent = symbiont_density_dependent_loop(i);
    % Mutation settings
    p_mut = 0.02; % probability of mutation
    muta_max = 0.5; % maximum mutation magnitude
    muta_min = 0;
    exp_mu = .25;
    muta_alpha = 1; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 1;
    muta_epsilon_symbiont = 1;
    
    mutualism_percent = zeros(1,nb_rep);
    parfor ii = 1:nb_rep
        host = zeros(space_size^2,1);
        symbiont = zeros(space_size^2,1);
        host_d = zeros(space_size^2,1);
        symbiont_d = zeros(space_size^2,1);
        total = initial_condis;
        h = sum(total>0,'all');
        s = sum(total==3,'all');
        % trait alpha
        alpha = 0;
        h_alpha = ones(h,1).*alpha;
        s_alpha = ones(s,1).*alpha;
        host(total>0) = h_alpha;
        symbiont(total==3) = s_alpha;
        % trait epsilon
        epsilon = 1;
        h_epsilon = ones(h,1).*epsilon;
        s_epsilon = ones(s,1).*epsilon;
        host_d(total>0) = h_epsilon;
        symbiont_d(total==3) = s_epsilon;
    
        synchrone = ecoevo_model(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                 neighbourhood, host_mortality, symbiont_mortality,...
                 host_density_dependent, host_germination, host_competition,...
                 symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)

        symbiont = synchrone.symbiont;
        simulation = synchrone.simulation;
        symbiont(simulation~=3) = NaN; % to not confound 0 because of empty cell with 0 because of a trait=0       
        mutualism_percent(ii) = sum(symbiont>0.475)/sum(simulation==3);
    end 
    Mutualism_percent(:,i) = mutualism_percent;
%     if mod(i,10)==0
%     save(['effect_of_competition','.mat'],'Mutualism_percent','host_density_dependent_loop','nb_linspace','nb_rep','params','-v7.3');
%     end
end
save(['effect_of_competition','.mat'],'Mutualism_percent','host_density_dependent_loop','nb_linspace','nb_rep','params','-v7.3');


%% PANEL C
% load a simulation with transition
load('mutualism_emerg_gamma02.mat')
Mutualist_emerg(simulation_emerg~=3) = NaN;       
Mutualism_percent = sum(Mutualist_emerg>0.475,1)./sum(simulation_emerg==3,1);
       
% continue the simulation but with a weak competition
% Host Germination and Competition
gamma = 2; % weak competition
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
tmax = 5000;
neighbourhood = 8;
% Initial conditions
initial_condis = simulation_emerg(:,end);
host = Host_emerg(:,end);
host_d = Host_d_emerg(:,end);
symbiont = Mutualist_emerg(:,end);
symbiont_d = Mutualist_d_emerg(:,end);
total = initial_condis;

synchrone = ecoevo_model(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
             neighbourhood, host_mortality, symbiont_mortality,...
             host_density_dependent, host_germination, host_competition,...
             symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
             muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)
% symbiont and simulation variables need to be "savec during dynamics"
% see : line 285 in ecoevo_model.m
out_symbiont = synchrone.symbiont;
out_simulation = synchrone.simulation;
out_symbiont(out_simulation~=3) = NaN;       
Mutualism_percent = [Mutualism_percent,...
    sum(out_symbiont>0.475,1)./sum(out_simulation==3,1)];

% Figure A3 is exactly the same procedure but saving the details of
% variables


