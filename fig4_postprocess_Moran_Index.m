clear all; clc;

%% POST PROCESS INTRASPECIFIC Moran INDEX
load('Index_Moran.mat')

[nsimu,nt] = size(I_S(:,2:end)); 
tt = 100*(1:nt);

I_S_med = quantile(I_S(:,2:end),[0.05,0.5,0.95]);
I_H_med = quantile(I_H(:,2:end),[0.05,0.5,0.95]);

%% POST PROCESS INTERSPECIFIC Correlation coefficients
load('Corr_Host_Symbiont.mat')

[nsimu,nt] = size(C_HS(:,2:end)); 

tt = 100*(1:nt);

C_HS_med = quantile(C_HS(:,2:end),[0.05,0.5,0.95]);
P_HS_med = quantile(P_HS(:,2:end),[0.05,0.5,0.95]);


Color = [0.8,0.8,0.8;...
         0.3,0.3,0.3;...
         0.05,0.05,0.05];

figure(1)
clf
hold on
plot(tt,I_S_med(2,:),'--','color','k','linewidth',2)
plot(tt,I_H_med(2,:),'-','color','k','linewidth',3)
plot(tt,C_HS_med(2,:),'-.','color','k','linewidth',1.5)

%%% Symbiont
ttt = [tt, fliplr(tt)];
inBetween = [I_S_med(1,:), fliplr(I_S_med(3,:))];
f = fill(ttt, inBetween,Color(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,I_S_med(2,:),'--','color','k','linewidth',2)


%%% Host
ttt = [tt, fliplr(tt)];
inBetween = [I_H_med(1,:), fliplr(I_H_med(3,:))];
f = fill(ttt, inBetween,Color(2,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,I_H_med(2,:),'-','color','k','linewidth',3)



%%% Symbiont
ttt = [tt, fliplr(tt)];
inBetween = [C_HS_med(1,:), fliplr(C_HS_med(3,:))];
f = fill(ttt, inBetween,Color(3,:));
set(f,'EdgeColor','none','FaceAlpha', 0.6)
plot(tt,C_HS_med(2,:),'-.','color','k','linewidth',1.5)

line([tt(1),tt(end)],[0 0],'linewidth',1,'color','k')
line([2000,2000],[-0.2 1],'linewidth',1.5,'color','r','linestyle','--')

fontsi = 16;
xlabel('time','fontsize',fontsi,'interpreter','latex');
ylabel('Assortment index','fontsize',fontsi,'interpreter','latex');
legend('Symbiont Moran index','Host Moran index','Correlation coefficient host/symbiont','location','northeast')
set(gca,'FontSize',fontsi,  'LineWidth'   , 1,...
    'FontName'   , 'Helvetica');

% ,...
%     'TickDir'     , 'out', 'TickLength'  , [.02 .02] , ...
%     'XColor'      , [0 0 0],'YColor'      , [0 0 0], ...
%     'Box'         , 'off');


% %%% Host
% ttt = [tt, fliplr(tt)];
% inBetween = [P_HS_med(1,:), fliplr(P_HS_med(3,:))];
% f = fill(ttt, inBetween,Color(2,:));
% set(f,'EdgeColor','none','FaceAlpha', 0.5)
% plot(tt,P_HS_med(2,:),'-','color','k','linewidth',3)

