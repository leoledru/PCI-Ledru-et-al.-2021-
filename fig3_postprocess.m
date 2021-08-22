%% FIGURE 3

load('simulations_gamma02_g.mat')
idx = find(~cellfun('isempty',out_simulation));
idx_100 = idx(1:100); % 100 simus with emergence

host_alpha = [];
host_epsilon = [];
symbiont_alpha = [];
symbiont_epsilon = [];
for j = 1:length(idx)
    i = idx(j);
    simu = out_simulation{i};
    host = out_host{i};
    host_d = out_host_d{i};
    symbiont = out_mutualist{i};
    symbiont_d = out_mutualist_d{i};
    host(simu==0) = NaN;
    host_d(simu==0) = NaN;
    symbiont(simu~=3) = NaN;
    symbiont_d(simu~=3) = NaN;
  
    host_alpha = [host_alpha;host(:,end)];
    host_epsilon = [host_epsilon;host_d(:,end)];
    symbiont_alpha = [symbiont_alpha;symbiont(:,end)];
    symbiont_epsilon = [symbiont_epsilon;symbiont_d(:,end)];
end

close
symbiont_alpha(isnan(symbiont_alpha)) = [];
symbiont_epsilon(isnan(symbiont_epsilon)) = [];
host_alpha(isnan(host_alpha)) = [];
host_epsilon(isnan(host_epsilon)) = [];

sample = 200000;
% linear regression for host
x = host_alpha(1:sample);
x(isnan(x)) = [];
y = host_epsilon(1:sample);
y(isnan(y)) = [];
p = polyfit(x,y,1);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal

% plots
subplot(1,2,1)
p = scatter(symbiont_alpha(1:sample),symbiont_epsilon(1:sample),50,'ok','filled');
p.MarkerFaceAlpha = .01;
xlabel('symbiont interaction trait $\alpha_h$','interpreter','latex','FontSize',20)
ylabel('symbiont dispersal trait $\epsilon_h$','interpreter','latex','FontSize',20)
subplot(1,2,2)
p = scatter(host_alpha(1:sample),host_epsilon(1:sample),50,'ok','filled');
p.MarkerFaceAlpha = .1;
hold on
plot(x,yfit,'r')
xlabel('host interaction trait $\alpha_h$','interpreter','latex','FontSize',20)
ylabel('host dispersal trait $\epsilon_h$','interpreter','latex','FontSize',20)

