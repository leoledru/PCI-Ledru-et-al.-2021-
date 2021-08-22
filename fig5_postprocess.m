%% FIGURE 5

%% panel b
load('effect_of_competition.mat')
mutualism_percent = Mutualism_percent;
load('effect_of_competition2.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_competition3.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_competition4.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_competition5.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];

% moving average to smooth out the stochasticity
mutualism_percent_mean = movmean(mutualism_percent,5,2);

p = stdshade(mutualism_percent_mean,.3,'k');
p.LineWidth = 3;
p.LineStyle = '-';
ylim([0 .40])
xlim([1 21])
xticklabels(.1:.1:1);
yticklabels(0:5:40);
ax = get(gca);
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xlabel('host competition','interpreter','latex','FontSize',26)
ylabel('\% of mutualist symbionts','interpreter','latex','Fontsize',26)
text(-.02,-.04,'global competition','interpreter','latex','FontSize',20)
text(17.5,-.04,'local competition','interpreter','latex','FontSize',20)

hold on
plot([0 21],[.1 0.1],'--r','LineWidth',2)
hold off

%% panel a
load('effect_of_gamma.mat');
mutualism_percent = Mutualism_percent;
load('effect_of_gamma2.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_gamma3.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_gamma4.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];
load('effect_of_gamma5.mat');
mutualism_percent = [mutualism_percent; Mutualism_percent];

mutualism_percent = mutualism_percent(:,3:end); % remove columns 1 and 2 -> NaN values because of non-viable conditions
% moving average to smooth out the stochasticity
mutualism_percent_mean = movmean(mutualism_percent,1,2);
p = stdshade(mutualism_percent_mean,.3,'k');
p.LineWidth = 3;
p.LineStyle = '-';
ylim([0 .40])
xlim([1 19])
xticks(1:2:21)
xticklabels(.2:.2:2);
yticklabels(0:5:40);
ax = get(gca);
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xlabel('host competition','interpreter','latex','FontSize',26)
ylabel('\% of mutualist symbionts','interpreter','latex','Fontsize',26)
text(0,-.04,'strong competition','interpreter','latex','FontSize',20)
text(16,-.04,'weak competition','interpreter','latex','FontSize',20)

hold on
plot([0 21],[.1 0.1],'--r','LineWidth',2)
hold off

%% panel c
close
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

close
windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smooth_line = filter(b,a,mutualism_percent(1:100:end));
smooth_line_2 = filter(b,a,smooth_line);

% raw
% plot((1:100:15001),mutualism_percent(1:100:end),'-ok','LineWidth',.1,'MarkerSize',3)
% smooth
plot((1:100:15001),smooth_line_2.*100,'k','LineWidth',4)

xlim([0 15000])
xlabel('time','interpreter','latex','FontSize',25)
ylabel('\% of mutualistic symbiont','interpreter','latex','FontSize',25)

