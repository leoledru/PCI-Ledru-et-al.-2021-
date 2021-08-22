%% FIGURE 6

load('map_morta_dispcost_cheat_fbmin07.mat')
out_persist = reshape(out_persist,sqrt(length(out_persist)),sqrt(length(out_persist)));
out_persist = ~out_persist;

load('map_morta_dispcost_late_fbmin07.mat'); nb_rep = 12;
out_h = out_host;
out_h_d = out_host_d;
out_m = out_mutualist;
out_m_d = out_mutualist_d;
out_simu = out_simulation;
load('map_morta_dispcost_late_fbmin07_2.mat'); nb_rep = nb_rep + 12;
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_m = [out_m, out_mutualist];
out_m_d = [out_m_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];
load('map_morta_dispcost_late_fbmin07_3.mat'); nb_rep = nb_rep + 12;
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_m = [out_m, out_mutualist];
out_m_d = [out_m_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];
load('map_morta_dispcost_late_fbmin07_4.mat'); nb_rep = nb_rep + 12;
out_h = [out_h, out_host];
out_h_d = [out_h_d, out_host_d];
out_m = [out_m, out_mutualist];
out_m_d = [out_m_d, out_mutualist_d];
out_simu = [out_simu, out_simulation];

out_host = out_h;
out_host_d = out_h_d;
out_mutualist = out_m;
out_mutualist_d = out_m_d;
out_simulation = out_simu;

size = sqrt(length(out_simulation));
% speciation test
for i = 1:length(out_simulation)
    for j = 1:48
        symbiont = out_mutualist{i,j};
        community = out_simulation{i,j};
        symbiont(community~=3) = NaN;

        symbiont_collapse_test(i,j) = sum(community==3)==0; % 1 = symbiont collapse
        host_collapse_test(i,j) = sum(community>0)==0; % 1 = host collapse (hence symbiont collapse)

        speciation_test(i,j) = sum(symbiont>.475)/sum(~isnan(symbiont));
        speciation_bool(i,j) = speciation_test(i,j)>.1; % 1 if speciation / else 0
    end
end

map = jet(100);
map(1,:) = [.15,.15,.15];
symbiont_collapse = mean(symbiont_collapse_test,2);
symbiont_collapse = reshape(symbiont_collapse,size,size);
host_collapse = mean(host_collapse_test,2);
host_collapse = reshape(host_collapse,size,size);
symbiont_collapse(symbiont_collapse==1) = NaN;
host_collapse(host_collapse==1) = NaN;

speciation_test(speciation_bool==0) = NaN; 
speciation = nanmean(speciation_test,2);
speciation = reshape(speciation,size,size); % vertical axis = morta / horizontal axis = disp cost

h = imagesc(speciation,'AlphaData',~isnan(symbiont_collapse));
hold on
contour(out_persist,1,'--k','LineWidth',4)
c = colorbar('FontSize',15);
ylabel(c,'probability of mutualism emergence','FontSize',20,'interpreter','latex')
colormap(map)
set(gca,'Color','w')
grid on
xticks(1:21)
xticklabels(round(Coeff_disp_h(1:21:441),2))
yticks(1:21)
yticklabels(round(Mutualist_mortality(1:21),3))

% host dependency
load('map_morta_dispcost_cheat_fbmin07.mat')
clearvars -except out_persist symbiont_collapse
load('map_morta_dispcost_dependance_late_fbmin07_simple_cheater_epsilon05.mat')
nb_linspace = 21;
Host_mortality = flip(linspace(.005,.15,nb_linspace)); % flip to have an increase of mortality along vertical axis
Host_mortality = repmat(Host_mortality,1,nb_linspace);
Mutualist_mortality = Host_mortality;
Coeff_disp_h = linspace(0,1,nb_linspace);
Coeff_disp_h = repelem(Coeff_disp_h,nb_linspace);
persist = [];
for ii = 1:length(out_simulation)
    simu = out_simulation{ii};
persist = [persist;sum(simu)==0];
end
persist = reshape(persist,sqrt(length(persist)),sqrt(length(persist)));
out_persist = reshape(out_persist,sqrt(length(out_persist)),sqrt(length(out_persist)));
% out_persist = out_persist';
persist(out_persist==1) = NaN;

hold on
contour(persist,1,'--r','LineWidth',4)
hold off
