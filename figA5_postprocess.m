%% FIGURE A5

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


close
coeff_disp_value_list = [1,7,15];
cost_label = [0,0.30,0.75];
C_all_s = cell(length(coeff_disp_value_list),1);
C_all_h = cell(length(coeff_disp_value_list),1);
for jj = 1:length(coeff_disp_value_list)
    X_h = [];
    X_s = [];
    Symbiont = [];
    Symbiont_d = [];
    Host = [];
    Host_d = [];
    cost = coeff_disp_value_list(jj);
    for ii = 1:nb_rep
        col = ii;
        morta = 18;
        row = morta + 21*(cost-1);

        symbiont = out_mutualist{row,col};
        symbiont_d = out_mutualist_d{row,col};
        host = out_host{row,col};
        host_d = out_host_d{row,col};

        simu = out_simulation{row,col};

        symbiont(simu~=3) = [];
        symbiont_d(simu~=3) = [];
        host(simu==0) = [];
        host_d(simu==0) = [];

        X_h_ii = [host,host_d];
        X_s_ii = [symbiont,symbiont_d];
        X_h = [X_h;X_h_ii];
        X_s = [X_s;X_s_ii];
        Symbiont = [Symbiont;symbiont];
        Symbiont_d = [Symbiont_d;symbiont_d];
        Host = [Host;host];
        Host_d = [Host_d;host_d];

        mutualist = symbiont(symbiont>.475);
        mutualist_d = symbiont_d(symbiont>.475);
        parasite = symbiont(symbiont<=.475);
        parasite_d = symbiont_d(symbiont<=.475);

        subplot(1,3,jj)
        s = scatter(mutualist,mutualist_d,4,'or','filled');
        hold on
        xlim([0 1])
        ylim([0 1])
        s.MarkerFaceAlpha = .05;
        s = scatter(parasite,parasite_d,4,'ob','filled');
        s.MarkerFaceAlpha = .05;
        s = scatter(host,host_d,4,'og','filled');
        s.MarkerFaceAlpha = .05;
        xlabel('interaction trait $\alpha$','interpreter','latex','FontSize',25);
        ylabel('dispersal trait $\epsilon$','interpreter','latex','FontSize',25);
        title(['dispersal cost = ', num2str(cost_label(jj))],'interpreter','latex','FontSize',20)
    end
    [idx_h,C_h] = kmeans(X_h,1);
    if jj>1
        [idx_s,C_s] = kmeans(X_s,2);
        clust1 = [Symbiont(idx_s==1),Symbiont_d(idx_s==1)];
        clust2 = [Symbiont(idx_s==2),Symbiont_d(idx_s==2)];
        mean_s_1 = mean(clust1,1);
        mean_s_2 = mean(clust2,1);       
        std_s_1 = std(clust1,1);
        std_s_2 = std(clust2,1); 
        h = plotEllipses(mean_s_1,std_s_1);
        h.FaceColor = [0 0 1 .3];
        h.LineWidth = 3;
        h.LineStyle = '--';
        h.EdgeColor = [.3 .3 1];
        h = plotEllipses(mean_s_2,std_s_2);
        h.FaceColor = [1 0 0 .3];
        h.LineWidth = 3;
        h.LineStyle = '--';
        h.EdgeColor = [1 .3 .3];
    else
        [idx_s,C_s] = kmeans(X_s,1);
        CI_x = prctile(Symbiont, [(100-95)/2, 95+(100-95)/2]); 
        CI_y = prctile(Symbiont_d, [(100-95)/2, 95+(100-95)/2]); 
        std_s = [CI_x(2)-CI_x(1),CI_y(2)-CI_y(1)];
        h = plotEllipses(nanmean([Symbiont,Symbiont_d]),nanstd([Symbiont,Symbiont_d]));
        h.FaceColor = [0 0 1 .3];
        h.LineWidth = 3;
        h.EdgeColor = [.3 .3 1];
        h.LineStyle = '--';
    end
    C_all_s{jj} = C_s;
    C_all_h{jj} = C_h;
    CI_x = prctile(Host, [(100-95)/2, 95+(100-95)/2]); 
    CI_y = prctile(Host_d, [(100-95)/2, 95+(100-95)/2]); 
    std_h = [CI_x(2)-CI_x(1),CI_y(2)-CI_y(1)];
    h = plotEllipses(nanmean([Host,Host_d]),nanstd([Host,Host_d]));
    h.FaceColor = [0 1 0 .3];
    h.LineWidth = 3;
    h.EdgeColor = [.3 1 .3];
    h.LineStyle = '--';
    hold off
end
set(gcf, 'Position', get(0, 'Screensize'));

