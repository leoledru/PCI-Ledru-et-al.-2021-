%% FIGURE 4

load('assortment_measured_vs_rand.mat')
load('assortment_measured_vs_rand_interspe.mat')

close
% Mean and SD
% moving average to smooth out the stochasticity
window = 5;
assortment_host = movmean(assortment_h,window,2);
assortment_host(:,1) = 0;
assortment_symbiont = movmean(assortment_m,window,2);
assortment_symbiont(:,1) = 0;
assortment_inter = movmean(assortment,window,2);
assortment_inter(:,1) = 0;

p = stdshade(assortment_host,.4,'k');
p.LineWidth = 5;
p.LineStyle = '-';
hold on
ylabel('assortment index','interpreter','latex','FontSize',25)
xlabel('time','interpreter','latex','FontSize',25)
p = stdshade(assortment_symbiont,.2,'k');
p.LineWidth = 3;
p.LineStyle = '--';
p = stdshade(assortment_inter,.6,'k');
p.LineWidth = 2;
p.LineStyle = '-.';
legend('host intraspecific','symbiont intraspecific','interspecific',...
    'interpreter','latex','FontSize',15,'Location','northwest')
xlim([0 100])
plot([0,100],[0 0],'k');
xticks(0:10:10000);
xticklabels(0:1000:10000);

