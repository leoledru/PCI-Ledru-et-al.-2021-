%% FIGURE 2

%% panel a

load('data/time_to_speciation_gamma1.mat')
time = out_time_speciation;
proba_emergence = sum(1 - isnan(out_time_speciation))/length(out_time_speciation);
load('data/time_to_speciation_gamma1_dispcost45.mat')
time = [time,out_time_speciation];
proba_emergence = [proba_emergence,sum(1 - isnan(out_time_speciation))/length(out_time_speciation)];
load('data/time_to_speciation.mat')
time = [time,out_time_speciation];
proba_emergence = [proba_emergence,sum(1 - isnan(out_time_speciation))/length(out_time_speciation)];
load('data/time_to_speciation_gamma02_dispcost45.mat')
time = [time,out_time_speciation];
proba_emergence = [proba_emergence,sum(1 - isnan(out_time_speciation))/length(out_time_speciation)];
%% BOXPLOTS (OLD VERSION)
boxplot(time)
ylabel('transition time (log scale)','interpreter','latex','FontSize',20)
ax = gca;
ax.YAxis.Scale = "log";
proba_emergence % show the probablities
%% HISTROGRAMS (NEW VERSION)
subplot(1,2,1)
h = histogram(time(:,1),15)
% h = histogram(time(:,1),logspace(2,5));
h.FaceColor = 'k';
h.FaceAlpha = .4;

hold on
h2 = histogram(time(:,3),15)
% h2 = histogram(time(:,3),logspace(2,5))
h2.FaceAlpha = .2;
h.BinEdges = h2.BinEdges;
hold off

% ax = gca;
% ax.XAxis.Scale ="log";

yticks(linspace(0,150,11));
xlabel('time','interpreter','latex','Fontsize',20)
ylabel('number of transitions','interpreter','latex','Fontsize',20)
legend({['weak competition' newline '$\gamma_C = 1$'],['strong competition' newline '$\gamma_C = 0.2$']}...
    ,'interpreter','latex','Fontsize',15)
title('dispersal cost = 0','interpreter','latex','Fontsize',20)

subplot(1,2,2)
h = histogram(time(:,2),16);
% h = histogram(time(:,2),logspace(2,2));
h.FaceColor = 'k';
h.FaceAlpha = .4;
h.BinEdges = [100 213 326 439 552 665 778 891 1004 1117 1230 1343 1456 1569 1682 1795 1908];

hold on
h = histogram(time(:,4),16);
% h = histogram(time(:,4),logspace(2,2));
h.FaceAlpha = .2;
hold off

yticks(linspace(0,300,11));
xlabel('time','interpreter','latex','Fontsize',20)
ylabel('number of transitions','interpreter','latex','Fontsize',20)
legend({['weak competition' newline '$\gamma_C = 1$'],['strong competition' newline '$\gamma_C = 0.2$']}...
    ,'interpreter','latex','Fontsize',15)
title('dispersal cost = 0.45','interpreter','latex','Fontsize',20)

%% panel b
load('data/simulations_gamma02_h.mat')

out_time_speciation(out_time_speciation==0) = NaN; % no emergence
low_bound = 2000;
up_bound = 5000;
% conserve only simus with >2000 time steps before and after emergence
idx = find(out_time_speciation>low_bound & out_time_speciation<up_bound);

host_density = [];
symbiont_density = [];
host_mean_alpha = [];
symbiont_mean_alpha = [];
symbiont_parasite = [];
symbiont_mutualiste = [];
parasite_density = [];
mutualiste_density = [];
mutualiste_mean_trait = [];
parasite_mean_trait = [];

for j = 1:length(idx)
    i = idx(j);
    
    p_density = out_parasite_density{i};
    p_density = p_density(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    m_density = out_mutualiste_density{i};
    m_density = m_density(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    parasite_density = [parasite_density;p_density];
    mutualiste_density = [mutualiste_density;m_density];
    p_mean_trait = out_symbiont_parasite_mean_trait{i};
    p_mean_trait = p_mean_trait(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    parasite_mean_trait = [parasite_mean_trait;p_mean_trait];
    m_mean_trait = out_symbiont_mutualist_mean_trait{i};
    m_mean_trait = m_mean_trait(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    mutualiste_mean_trait = [mutualiste_mean_trait;m_mean_trait];
    
    s_density = out_symbiont_density{i};
    s_density = s_density(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    h_density = out_host_density{i};
    h_density = h_density(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    s_mean_trait = out_symbiont_mean_trait{i};
    s_mean_trait = s_mean_trait(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    h_mean_trait = out_host_mean_trait{i};
    h_mean_trait = h_mean_trait(out_time_speciation(i)-low_bound:out_time_speciation(i)+up_bound);
    s_parasite = out_symbiont_parasite_mean_trait{i};
    s_parasite = s_parasite(1:5000);
    s_mutualiste = out_symbiont_mutualist_mean_trait{i};
    s_mutualiste = s_mutualiste(1:5000);
    
    host_density = [host_density;h_density];
    symbiont_density = [symbiont_density;s_density];
    host_mean_alpha = [host_mean_alpha;h_mean_trait];
    symbiont_mean_alpha = [symbiont_mean_alpha;s_mean_trait];
    symbiont_parasite = [symbiont_parasite;s_parasite];
    symbiont_mutualiste = [symbiont_mutualiste;s_mutualiste];
end

close
N = 50;
mean_host_density = mean(host_density);
std_host_density = std(host_density);
mean_symbiont_density = mean(symbiont_density);
mean_parasite_density = mean(parasite_density);
std_parasite_density = std(parasite_density);
mean_mutualiste_density = mean(mutualiste_density);
std_mutualiste_density = std(mutualiste_density);
col = jet(100);

relative_mutualiste_density = mean_mutualiste_density./(mean_parasite_density+mean_mutualiste_density);
for i = 1:N:numel(mean_host_density(1,:))-N
%     plot(i,mean_host_density(i:i+N),'.','Color',[mean(host_mean_alpha(:,i))/max(mean(host_mean_alpha)),0,0])
%     plot(i,mean_host_density(i:i+N),'.','Color',[mean(host_mean_alpha(:,i)),.3,.3])
    subplot(1,2,1)
    col_idx = round(mean(host_mean_alpha(:,i))*100);
    plot(i,mean_host_density(i:i+N),'.','Color',col(col_idx,:),'MarkerSize',10)
    hold on
%     plot(i,mean_symbiont_density(i:i+N),'.','Color',[mean(symbiont_mean_alpha(:,i))/max(mean(symbiont_mean_alpha)),0,0])
%     plot(i,mean_symbiont_density(i:i+N),'.','Color',[mean(symbiont_mean_alpha(:,i)),.3,.3])
%     col_idx = round(mean(symbiont_mean_alpha(:,i))*100);
%     plot(i,mean_symbiont_density(i:i+N),'.','Color',col(col_idx,:),'MarkerSize',10)
    subplot(1,2,2)
    col_idx = round(mean(parasite_mean_trait(:,i))*100);
%     relative_parasite_density = mean_parasite_density(i:i+N)/(mean_parasite_density(i:i+N) + mean_mutualiste_density(i:i+N));
    plot(i,mean_parasite_density(i:i+N),'.','Color',col(col_idx,:),'MarkerSize',10) % density
%     plot(i,relative_parasite_density,'.','Color',col(col_idx,:),'MarkerSize',10) % relative density
    hold on
    col_idx = round(nanmean(mutualiste_mean_trait(:,i))*100);
%     relative_mutualiste_density = mean_mutualiste_density(i:i+N)/(mean_parasite_density(i:i+N) + mean_mutualiste_density(i:i+N));
    plot(i,mean_mutualiste_density(i:i+N),'.','Color',col(col_idx,:),'MarkerSize',10) % density
%     plot(i,relative_mutualiste_density,'.','Color',col(col_idx,:),'MarkerSize',10) % relative density

end
subplot(1,2,1)
mean_host_density = movmean(mean_host_density,100);
std_host_density = movmean(std_host_density,100);
x = 1:numel(mean_host_density(1,:));
x2 = [x, fliplr(x)];
inBetween = [mean_host_density-std_host_density, fliplr(mean_host_density+std_host_density)];
h = fill(x2, inBetween, 'k');
set(h,'facealpha',.05);
xlim([0,7000]);
xlabel('time','interpreter','latex','FontSize',22)
ylabel('host density','interpreter','latex','FontSize',22)
colormap(col)
c = colorbar('FontSize',12);
ylabel(c,'mean interaction trait','FontSize',22,'interpreter','latex')

subplot(1,2,2)
% yyaxis left
mean_parasite_density = movmean(mean_parasite_density,100);
std_parasite_density = movmean(std_parasite_density,100);
x = 1:numel(mean_host_density(1,:));
x2 = [x, fliplr(x)];
inBetween = [mean_parasite_density-std_parasite_density, fliplr(mean_parasite_density+std_parasite_density)];
h = fill(x2, inBetween, 'k');
set(h,'facealpha',.05);
mean_mutualiste_density = movmean(mean_mutualiste_density,100);
std_mutualiste_density = movmean(std_mutualiste_density,100);
x = 1:numel(mean_host_density(1,:));
x2 = [x, fliplr(x)];
inBetween = [mean_mutualiste_density-std_mutualiste_density, fliplr(mean_mutualiste_density+std_mutualiste_density)];
h = fill(x2, inBetween, 'k');
set(h,'facealpha',.05);
xlim([0,7000]);
xlabel('time','interpreter','latex','FontSize',22)
ylabel('symbiont density','interpreter','latex','FontSize',22)
colormap(col)
c = colorbar('FontSize',12);
ylabel(c,'mean interaction trait','FontSize',22,'interpreter','latex')
ylim([0 .35])
text(2800,0.065,'mutualistic symbionts','interpreter','latex','FontSize',15)
text(2800,0.33,'parasitic symbionts','interpreter','latex','FontSize',15)

yyaxis right
plot(movmean(relative_mutualiste_density,100),'--','LineWidth',3,'Color',[.5 0 .5])
plot([0,7000],[.1,.1],'k-','LineWidth',2) % transition threshold
text(2800,0.15,{'relative density','of mutualistic symbionts'},'interpreter','latex','FontSize',15)
ylim([0 .25])
ax = gca;
ax.YAxis(2).Color = [.5 .0 .5];
hold off

