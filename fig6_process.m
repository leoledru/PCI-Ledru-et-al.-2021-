%% MEAN PERCENTAGE OF MUTUALISTIC SYMBIONTS
clear all; clc;

nb_linspace = 21;
nb_rep = 50;
% Mortality and dispersal cost for map
Host_mortality = flip(linspace(.005,.15,nb_linspace)); % flip to have an increase of mortality along vertical axis (see postprocess)
Host_mortality = repmat(Host_mortality,1,nb_linspace);
Symbiont_mortality = Host_mortality;
Coeff_disp_h = linspace(0,1,nb_linspace);
Coeff_disp_h = repelem(Coeff_disp_h,nb_linspace);
Coeff_disp = Coeff_disp_h;
n_loop = length(Symbiont_mortality);
% Host Germination and Competition
gamma = .2;
host_germination = 1;
host_density_dependent = 0; % 0=global, 1=local
host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;
% Symbiont Colonization and Competition
symbiont_colonization = 1;
symbiont_density_dependent = 0; % 0=global, 1=local
coeff_compet = 0;
gamma_symbiont = .2;

params_name = {'gamma','host_germination','host_density_dependent','symbiont_colonization','symbiont_density_dependent',...
    'coeff_compet','gamma_symbiont'};
params_value = [gamma,host_germination,host_density_dependent,symbiont_colonization,symbiont_density_dependent,...
    coeff_compet,gamma_symbiont];
for name_loop = 1:length(params_name)
    params{1,name_loop} = params_name{name_loop};
    params{2,name_loop} = params_value(name_loop);
end

out_mutualism_percent = zeros(n_loop,1);
out_persist = zeros(n_loop,1);
for i = 1:n_loop
    space_size = 100;
    tmax = 10000;
    neighbourhood = 8;
    % Mortality / dispersal cost
    host_mortality = Host_mortality(i);
    symbiont_mortality = Symbiont_mortality(i);
    coeff_disp_h = Coeff_disp_h(i);
    coeff_disp = Coeff_disp(i);
    % Mutation settings
    p_mut = 0.02; % probability of mutation
    muta_max = 0.5; % maximum mutation magnitude
    muta_min = 0;
    exp_mu = .25;
    muta_alpha = 1; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 1;
    muta_epsilon_symbiont = 1;
    
    mutualism_percent = [];
    persist = [];
    load('data/cheater_eq.mat')
    initial_condis = total_new;
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
                     muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu);
        
         mutualism_percent(ii) = synchrone.mutualiste_density(end);
         persist(ii) = synchrone.persist;
    end 
    out_mutualism_percent(i) = mean(mutualism_percent);
    out_persist(i) = mean(persist);
%     save(['map_morta_proba_emergence','.mat'],'out_mutualism_percent','-v7.3');
end
save(['map_morta_proba_emergence','.mat'],'out_mutualism_percent','out_persist','params','-v7.3');


%% VIABILITY OF THE PARASITIC SYSTEM
clear all; clc;

nb_linspace = 21;
% Mortality and dispersal cost for map
Host_mortality = flip(linspace(.005,.15,nb_linspace)); % flip to have an increase of mortality along vertical axis (see postprocess)
Host_mortality = repmat(Host_mortality,1,nb_linspace);
Symbiont_mortality = Host_mortality;
Coeff_disp_h = linspace(0,1,nb_linspace);
Coeff_disp_h = repelem(Coeff_disp_h,nb_linspace);
Coeff_disp = Coeff_disp_h;
n_loop = length(Symbiont_mortality);
% Host Germination and Competition
gamma = .2;
host_germination = 1;
host_density_dependent = 0; % 0=global, 1=local
host_competition = @(host_proportion) host_germination.*host_proportion.^gamma;
% Symbiont Colonization and Competition
symbiont_colonization = 1;
symbiont_density_dependent = 0; % 0=global, 1=local
coeff_compet = 1;
gamma_symbiont = .2;

out_simulation = cell(n_loop,1); 
out_host = cell(n_loop,1); 
out_host_d = cell(n_loop,1); 
out_symbiont = cell(n_loop,1);
out_symbiont_d = cell(n_loop,1);
out_persist = zeros(n_loop,1);

params_name = {'gamma','host_germination','host_density_dependent','symbiont_colonization','symbiont_density_dependent',...
    'coeff_compet','gamma_symbiont'};
params_value = [gamma,host_germination,host_density_dependent,symbiont_colonization,symbiont_density_dependent,...
    coeff_compet,gamma_symbiont];
for name_loop = 1:length(params_name)
    params{1,name_loop} = params_name{name_loop};
    params{2,name_loop} = params_value(name_loop);
end

parfor i = 1:n_loop
    space_size = 100;
    tmax = 5000;
    neighbourhood = 8;
    prop_init = 1;
    % Mortality / dispersal cost
    host_mortality = Host_mortality(i);
    symbiont_mortality = Symbiont_mortality(i);
    coeff_disp_h = Coeff_disp_h(i);
    coeff_disp = Coeff_disp(i);
    % Mutation settings
    p_mut = 0.02; % probability of mutation
    muta_max = 0.5; % maximum mutation magnitude
    muta_min = 0;
    exp_mu = .25;
    muta_alpha = 0; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 0;
    muta_epsilon_symbiont = 0;
    host = zeros(space_size^2,1);
    symbiont = zeros(space_size^2,1);
    host_d = zeros(space_size^2,1);
    symbiont_d = zeros(space_size^2,1);
    total = rand(space_size^2,1);
    total(total<prop_init) = 2;
    total(total<2) = 0;
    total(total==2) = 3;
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
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu);

    out_persist(i) = synchrone.persist;
    out_host{i} = synchrone.host;
    out_host_d{i} = synchrone.host_d;
    out_symbiont{i} = synchrone.symbiont;
    out_symbiont_d{i} = synchrone.symbiont_d;
    out_simulation{i} = synchrone.simulation;
    
end
save(['map_morta_dispcost_cheat','.mat'],'out_host','out_host_d','out_symbiont','out_persist',...
'out_symbiont_d','out_simulation','Symbiont_mortality','Coeff_disp_h','nb_linspace','params','-v7.3');


%% HOST DEPENDENCY
clear all; clc;

load('data/map_morta_dispcost_late_fbmin07.mat');
out_h = out_host;
out_h_d = out_host_d;
out_s = out_mutualist;
out_s_d = out_mutualist_d;
out_simu = out_simulation;
load('data/map_morta_dispcost_late_fbmin07_2.mat');
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_s = [out_s, out_mutualist];
out_s_d = [out_s_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];
load('data/map_morta_dispcost_late_fbmin07_3.mat');
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_s = [out_s, out_mutualist];
out_s_d = [out_s_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];
load('data/map_morta_dispcost_late_fbmin07_4.mat');
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_s = [out_s, out_mutualist];
out_s_d = [out_s_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];

out_host = out_h;
out_host_d = out_h_d;
out_symbiont = out_s;
out_symbiont_d = out_s_d;
out_simulation = out_simu;
nb_rep = size(out_simulation,2);
%%
% Settings for weak dependance
for row = 1:441
    for col = 1:nb_rep
    host = out_host{row,col};
    host_d = out_host_d{row,col};
    simu = out_simulation{row,col};
    host(simu==0) = NaN;
    host_d(simu==0) = NaN;
    Host_tot_ab(row,col) = sum(~isnan(host));
    Host_tot_d(row,col) = nanmean(host_d);
    Host_tot_alpha(row,col) = nanmean(host);
    end
end
Host_tot_ab = round(mean(Host_tot_ab,2));
Host_tot_d_min = min(Host_tot_d,[],2);
Host_tot_d = mean(Host_tot_d,2);
Host_tot_alpha = mean(Host_tot_alpha,2);
clearvars -except Host_tot_ab Host_tot_d_min Host_tot_d Host_tot_alpha

% Simulations
Host_mortality = flip(linspace(.005,.15,21)); % flip to have an increase of mortality along vertical axis (see postprocess)
Host_mortality = repmat(Host_mortality,1,21);
Symbiont_mortality = Host_mortality;
Coeff_disp_h = linspace(0,1,21);
Coeff_disp_h = repelem(Coeff_disp_h,21);
Coeff_disp = Coeff_disp_h;
n_loop = length(Symbiont_mortality);

% Host Germination and Competition
gamma = .2;
host_germination = 1;
host_density_dependent = 0; % 0=global, 1=local
coeff_disp_h = .05;
% Symbiont Colonization and Competition
symbiont_colonization = 1;
symbiont_density_dependent = 0; % 0=global, 1=local
coeff_compet = 1;
gamma_symbiont = .2;
coeff_disp = .05;

out_simulation = cell(n_loop,1); 
out_persist = cell(n_loop,1); 
out_host = cell(n_loop,1); 
out_host_d = cell(n_loop,1); 
params_name = {'gamma','host_germination','host_density_dependent','symbiont_colonization','symbiont_density_dependent',...
    'coeff_compet','gamma_symbiont'};
params_value = [gamma,host_germination,host_density_dependent,symbiont_colonization,symbiont_density_dependent,...
    coeff_compet,gamma_symbiont];
for name_loop = 1:length(params_name)
    params{1,name_loop} = params_name{name_loop};
    params{2,name_loop} = params_value(name_loop);
end

for i = 1:n_loop   
    space_size = 100;
    tmax = 5000;
    neighbourhood = 8;
    % Mortality / dispersal cost
    host_mortality = Host_mortality(i);
    symbiont_mortality = Symbiont_mortality(i);
    coeff_disp_h = Coeff_disp_h(i);
    coeff_disp = Coeff_disp_h(i);
    % Mutation settings
    p_mut = 0.02; % probability of mutation
    muta_max = 0.5; % maximum mutation magnitude
    muta_min = 0;
    exp_mu = .25;
    muta_alpha = 0; % 1 -> mutation on, 0-> mutation off
    muta_epsilon_host = 0;
    muta_epsilon_symbiont = 0;
    % setting for weak dependance % MEAN HOST
    host_abundance = Host_tot_ab(i);
    host_proportion = host_abundance/(space_size^2);
    host_epsilon = Host_tot_d(i);
    host_alpha = Host_tot_alpha(i);
    host_competition = host_germination.*host_proportion.^gamma;

    host = zeros(space_size^2,1);
    symbiont = zeros(space_size^2,1);
    host_d = zeros(space_size^2,1);
    symbiont_d = zeros(space_size^2,1);

    % Set only host, with the proper abundance : weak dependance
    total = zeros(space_size^2,1);
    idx = randperm(space_size^2);
    idx = idx(1:host_abundance);
    total(idx) = 1;

    h = sum(total>0,'all');
    s = sum(total==3,'all');
    % trait alpha
    alpha = host_alpha; % mean host alpha
    h_alpha = ones(h,1).*alpha;
    s_alpha = ones(s,1).*alpha;
    host(total>0) = h_alpha;
    symbiont(total==3) = s_alpha;
    % trait epsilon
    epsilon = host_epsilon; % mean host epsilon
    h_epsilon = ones(h,1).*epsilon;
    s_epsilon = ones(s,1).*epsilon;
    host_d(total>0) = h_epsilon;
    symbiont_d(total==3) = s_epsilon;

    % WITH FIXED COMPETITION : SEE LINES 153-157 IN ECOVO_MODEL.M
    % ALSO SEE LINES 328-349 ABOUT COMMENT/UNCOMMENT SECTION
    synchrone = ecoevo_model(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                 neighbourhood, host_mortality, symbiont_mortality,...
                 host_density_dependent, host_germination, host_competition,...
                 symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu);
    
    out_persist{i} = synchrone.persist;
    out_host{i} = synchrone.host;
    out_host_d{i} = synchrone.host_d;
    out_simulation{i} = synchrone.simulation;  
end
save(['map_morta_dispcost_dependance','.mat'],'out_host','out_host_d','out_simulation','out_persist','params','-v7.3');

