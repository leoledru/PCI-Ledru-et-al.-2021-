
clear all; clc;

%% POST PROCESS INTRASPECIFIC Moran INDEX
load('Aggregation_index.mat')


%% Figure
[nsimu,nt] = size(A_H); 

tt = 100*(1:nt);
tt_m = tt(4:end);
Color = [0.8,0.8,0.8;...
         0.3,0.3,0.3;...
         0.05,0.05,0.05];
     
A_H_med = quantile(A_H,[0.05,0.5,0.95]);
A_SM_med = quantile(A_SM(:,4:end),[0.05,0.5,0.95]);
A_SP_med = quantile(A_SP,[0.05,0.5,0.95]);

figure(1)
clf
hold on
plot(tt,A_H_med(2,:),'-','color','k','linewidth',3)
plot(tt_m,A_SM_med(2,:),'--','color','k','linewidth',2)
plot(tt,A_SP_med(2,:),'-.','color','k','linewidth',1.5)

%%% Symbiont mutualist
ttt = [tt_m, fliplr(tt_m)];
inBetween = [A_SM_med(1,:), fliplr(A_SM_med(3,:))];
f = fill(ttt, inBetween,Color(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt_m,A_SM_med(2,:),'--','color','k','linewidth',2)


%%% Symbiont parasit
ttt = [tt, fliplr(tt)];
inBetween = [A_SP_med(1,:), fliplr(A_SP_med(3,:))];
f = fill(ttt, inBetween,Color(2,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,A_SP_med(2,:),'-.','color','k','linewidth',1.5)


%%% Host
ttt = [tt, fliplr(tt)];
inBetween = [A_H_med(1,:), fliplr(A_H_med(3,:))];
f = fill(ttt, inBetween,Color(3,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,A_H_med(2,:),'-','color','k','linewidth',3)

line([2000,2000],[0 0.5],'linewidth',1.5,'color','r','linestyle','--')

fontsi = 16;
xlabel('time','fontsize',fontsi,'interpreter','latex');
ylabel('Aggregation index','fontsize',fontsi,'interpreter','latex');
legend('Host','Mutualist symbiont','Parasitic symbiont','location','southeast')
set(gca,'FontSize',fontsi,  'LineWidth'   , 1,...
    'FontName'   , 'Helvetica');
