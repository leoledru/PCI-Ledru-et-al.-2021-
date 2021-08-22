%% FIGURE A3
load('mutualism_emerg_gamma02_extinction_gamma2.mat')
Mutualist_tot(simulation_tot~=3) = NaN;
Mutualist_d_tot(simulation_tot~=3) = NaN;
Host_tot(simulation_tot==0) = NaN;
Host_d_tot(simulation_tot==0) = NaN;
host_mean = nanmean(Host_tot,1);
symbiont_mean = nanmean(Mutualist_tot,1);
mutualism_percent = sum(Mutualist_tot>.475,1)./sum(~isnan(Mutualist_tot),1);
host_density = sum(simulation_tot~=0,1);
symbiont_density = sum(simulation_tot==3,1);

subplot('Position',[0.06, .55, 0.4, 0.4]); % [left,bottom,width,height]
plot(host_mean,'-k','LineWidth',3)
hold on
plot(symbiont_mean,'-k','LineWidth',1)
plot((1:100:15001),mutualism_percent(1:100:end),'-ok','LineWidth',.1,'MarkerSize',3)
hold off
legend({'host','symbiont',['\% of mutualistic' 10 'symbiont']},'interpreter','latex','FontSize',12,'Location','northeast')
xlabel('time','interpreter','latex','FontSize',20)
ylabel('mean interaction trait','interpreter','latex','FontSize',20)
xlim([0 15000])
subplot('Position',[0.55, .55, 0.4, 0.4])
plot(host_density,'-k','LineWidth',3)
hold on
plot(symbiont_density,'-k','LineWidth',1)
hold off
legend('host','symbiont','interpreter','latex','FontSize',12,'Location','northeast')
xlabel('time','interpreter','latex','FontSize',20)
ylabel('density','interpreter','latex','FontSize',20)
xlim([0 15000])
subplot('Position',[0.06, .1, 0.28, 0.3]);
sc1 = scatter(Mutualist_tot(:,1000),Mutualist_d_tot(:,1000),'k','filled');
sc1.MarkerFaceAlpha = .15;
hold on
sc1 = scatter(Host_tot(:,1000),Host_d_tot(:,1000),'r','filled');
sc1.MarkerFaceAlpha = .15;
legend('symbiont','host','interpreter','latex','FontSize',15,'Location','southeast')
xlabel('dispersal trait','interpreter','latex','FontSize',20)
ylabel('interaction trait','interpreter','latex','FontSize',20)
subplot('Position',[0.38, .1, 0.28, 0.3]);
sc1 = scatter(Mutualist_tot(:,8000),Mutualist_d_tot(:,8000),'k','filled');
sc1.MarkerFaceAlpha = .15;
hold on
sc1 = scatter(Host_tot(:,8000),Host_d_tot(:,8000),'r','filled');
sc1.MarkerFaceAlpha = .15;
subplot('Position',[0.7, .1, 0.28, 0.3]);
sc1 = scatter(Mutualist_tot(:,13000),Mutualist_d_tot(:,13000),'k','filled');
sc1.MarkerFaceAlpha = .15;
hold on
sc1 = scatter(Host_tot(:,13000),Host_d_tot(:,13000),'r','filled');
sc1.MarkerFaceAlpha = .15;