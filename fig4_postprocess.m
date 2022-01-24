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

%% PANELS B & C

% load('data/simulations_gamma02_d.mat')

host = out_host{52};
symbiont = out_mutualist{52};
simulation = out_simulation{52};

before = 10;
after = 100;

simu = simulation(:,after);
h = host(:,after);
s = symbiont(:,after);

h(simu==0) = NaN; % do not confound trait=0 with empty cell
h(h<.475) = 0;
h(h>=.475) = 1;
s(simu~=3) = NaN;
s(s>=.475) = 3;
s(s<.475) = 1;
s(isnan(s)) = 0; % allow the addition of matrices h_reshape + s_reshape below

h_reshape = reshape(h,100,100);
s_reshape = reshape(s,100,100);
simu_reshape = reshape(simu,100,100);

to_show = h_reshape + s_reshape;
to_show(simu_reshape==1) = 0; % reset a value of 0 for mutualistic host ALONE
to_show = to_show + 1;
to_show(isnan(to_show)) = 0;

% 0 = empty cell ; 1 = host alone ; 2 = host non-mut with symbiont parasitic ;
% 3 = host mut with symbiont parasitic ; 4 = host non-mut with symbiont mutualistic ;
% 5 = host mutualistic with symbiont mutualistic
map = [0,0,0;
    0,1,0;
    0,0,1;
    .6,.5,.2;
    .5,0,.5;
    1,0,0];

imagesc(to_show)
imagesc(to_show(1:50,50:100))
colormap(map)
caxis([0 6]);


    







