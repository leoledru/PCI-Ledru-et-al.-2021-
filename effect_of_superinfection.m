%% FIG 2 PANEL B & FIG 3
clear all; clc;

n_loop = 20;
superinfection_max_loop = linspace(.1,2,11);
n_rep = length(superinfection_max_loop);

out_host = cell(n_loop,n_rep);
out_host_d = cell(n_loop,n_rep);
out_symbiont = cell(n_loop,n_rep);
out_symbiont_d = cell(n_loop,n_rep);
out_host_mean_trait = cell(n_loop,n_rep);
out_symbiont_mean_trait = cell(n_loop,n_rep);
out_time_speciation = zeros(n_loop,n_rep);
out_symbiont_density = cell(n_loop,n_rep);
out_host_density = cell(n_loop,n_rep);
out_symbiont_mutualist_mean_trait = cell(n_loop,n_rep);
out_symbiont_parasite_mean_trait = cell(n_loop,n_rep);
out_parasite_density = cell(n_loop,n_rep);
out_mutualiste_density = cell(n_loop,n_rep);
out_simulation = cell(n_loop,n_rep);

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
coeff_disp = .45;
coeff_disp_h = .45;

% Initial conditions
load('~/Documents/MATLAB/Mycorhize/data/cheater_eq.mat') 
initial_condis = total_new;

for j = 1:length(superinfection_max_loop)
    superinfection_max = superinfection_max_loop(j);
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
        tmax = 5000;
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
        
        %     synchrone = ecoevo_model_beta_poisson(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
        %                  neighbourhood, host_mortality, symbiont_mortality,...
        %                  host_density_dependent, host_germination, host_competition,...
        %                  symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
        %                  muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)
        %     
        synchrone = ecoevo_model_beta_poisson_superinfection_b(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                     neighbourhood, host_mortality, symbiont_mortality,...
                     host_density_dependent, host_germination, host_competition,...
                     symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                     muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,...
                     coeff_disp_h,exp_mu,superinfection_max)
        
%         speciation_test = synchrone.speciation;
        %     out_symbiont_density{i} = synchrone.symbiont_density;
%         out_parasite_density{i,j} = synchrone.parasite_density;
%         out_mutualiste_density{i,j} = synchrone.mutualiste_density;
%         out_host_density{i,j} = synchrone.host_density;
%         out_host_mean_trait{i,j} = synchrone.host_mean_trait;
        %     out_symbiont_mean_trait{i} = synchrone.symbiont_mean_trait;
%         out_symbiont_mutualist_mean_trait{i,j} = synchrone.symbiont_mutualist_mean_trait;
%         out_symbiont_parasite_mean_trait{i,j} = synchrone.symbiont_parasite_mean_trait;
        out_time_speciation(i,j) = synchrone.time_speciation;
        out_host{i,j} = synchrone.host;
        out_host_d{i,j} = synchrone.host_d;
        out_symbiont{i,j} = synchrone.symbiont;
        out_symbiont_d{i,j} = synchrone.symbiont_d;
        out_simulation{i,j} = synchrone.simulation;
        
    end
    save(['test_superinfection_dispcost45_2','.mat'],'out_simulation',...
     'out_time_speciation','out_host','out_host_d','out_symbiont','out_symbiont_d','-v7.3')
end
save(['test_superinfection_dispcost45_2','.mat'],'out_simulation',...
 'out_time_speciation','out_host','out_host_d','out_symbiont','out_symbiont_d','-v7.3')

%% POSTPROCESS
load('test_superinfection_dispcost45.mat')
Host = out_host;
Host_d = out_host_d;
Symbiont = out_symbiont;
Symbiont_d = out_symbiont_d;
Simulation = out_simulation;
Time_speciation = out_time_speciation;

load('test_superinfection_dispcost45_2.mat')
Host = [Host;out_host];
Host_d = [Host_d;out_host_d];
Symbiont = [Symbiont;out_symbiont];
Symbiont_d = [Symbiont_d;out_symbiont_d];
Simulation = [Simulation;out_simulation];
Time_speciation = [Time_speciation;out_time_speciation];
%%
superinfection_max_loop = linspace(.1,2,11);
proba_transi = mean(~isnan(Time_speciation),1);

subplot(1,2,1)
plot(superinfection_max_loop,proba_transi,'-*k','LineWidth',3)
xlim([min(superinfection_max_loop),max(superinfection_max_loop)]);
xticks(superinfection_max_loop);
xlabel('intensity of superinfection','interpreter','latex','FontSize',20)
ylabel('probability of transition','interpreter','latex','FontSize',20)
%%
Mut_percent = [];
for i = 1:size(Symbiont,2)
    mut_percent = [];
    for j = 1:size(Symbiont,1)
        s = Symbiont{j,i};
        simu = Simulation{j,i};
        s(simu~=3) = NaN;
        mut_percent = [mut_percent,sum(s>=0.475)/sum(simu==3)];
    end
    Mut_percent = [Mut_percent, nanmean(mut_percent)];
end

subplot(1,2,2)
plot(superinfection_max_loop,Mut_percent.*100,'-*k','LineWidth',3)
xlim([min(superinfection_max_loop),max(superinfection_max_loop)]);
xticks(superinfection_max_loop);
xlabel('intensity of superinfection','interpreter','latex','FontSize',20)
ylabel('\% of mutualistic symbionts','interpreter','latex','FontSize',20)
%%
S_tot = [];
S_tot_d = [];
% h = [];
% h_d = [];
for i = 1:size(Symbiont,2)
    S = [];
    S_d = [];
    for j = 1:size(Symbiont,1)
        s = Symbiont{j,i};
        s_d = Symbiont_d{j,i};
        simu = Simulation{j,i};
        s(simu~=3) = NaN;
        s_d(simu~=3) = NaN;
        S = [S;s];
        S_d = [S_d;s_d];
    end
    S_tot = [S_tot,S];
    S_tot_d = [S_tot_d,S_d];
end
%%
subplot(2,5,1)
scatter(S_tot(:,1),S_tot_d(:,1),10,'filled')
xlabel('interaction trait','interpreter','latex','FontSize',15)
ylabel('dispersal trait','interpreter','latex','FontSize',15)
title(['superfection = ', num2str(superinfection_max_loop(1))])
subplot(2,5,2)
scatter(S_tot(:,2),S_tot_d(:,2),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(2))])
subplot(2,5,3)
scatter(S_tot(:,3),S_tot_d(:,3),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(3))])
subplot(2,5,4)
scatter(S_tot(:,4),S_tot_d(:,4),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(4))])
subplot(2,5,5)
scatter(S_tot(:,5),S_tot_d(:,5),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(5))])
subplot(2,5,6)
scatter(S_tot(:,6),S_tot_d(:,6),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(6))])
subplot(2,5,7)
scatter(S_tot(:,7),S_tot_d(:,7),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(7))])
subplot(2,5,8)
scatter(S_tot(:,8),S_tot_d(:,8),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(8))])
subplot(2,5,9)
scatter(S_tot(:,9),S_tot_d(:,9),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(9))])
subplot(2,5,10)
scatter(S_tot(:,10),S_tot_d(:,10),10,'filled')
title(['superfection = ', num2str(superinfection_max_loop(10))])
