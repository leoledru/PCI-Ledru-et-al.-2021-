clear all; clc;

%% POST PROCESS INTRASPECIFIC Moran INDEX
load('Index_Moran.mat')

[nsimu,nt] = size(IM_S(:,2:end)); 
tt = 100*(1:nt);

IM_H =  -seuil_H + IM_H;
IM_S = -seuil_S + IM_S;

IG_H =  1- IG_H;
IG_S = 1- IG_S;

IL_H = IL_H-1;
IL_S = IL_S-1;

IM_S_med = quantile(IM_S(:,2:end),[0.05,0.5,0.95]);
IM_H_med = quantile(IM_H(:,2:end),[0.05,0.5,0.95]);
IG_S_med = quantile(IG_S(:,2:end),[0.05,0.5,0.95]);
IG_H_med = quantile(IG_H(:,2:end),[0.05,0.5,0.95]);
IL_S_med = quantile(IL_S(:,2:end),[0.05,0.5,0.95]);
IL_H_med = quantile(IL_H(:,2:end),[0.05,0.5,0.95]);

% %% WEIGHT MATRIX OF 8 NEIGHBORS
% iw = [];
% jw = [];
% ww = [];
% nW = 100;
% parfor ii=1:nW^2
%     [row,col] = ind2sub([nW nW],ii); % coordinates IJ of the cell ii
%         neighboor = [(row-1)*(row>1)+nW*(row==1),col; % coordinates IJ of the 8 neighboors (wrap space)
%             (row+1)*(row<nW)+1*(row==nW),col;
%             row,(col-1)*(col>1)+nW*(col==1);
%             row,(col+1)*(col<nW)+1*(col==nW);
%             (row-1)*(row>1)+nW*(row==1),(col-1)*(col>1)+nW*(col==1);
%             (row+1)*(row<nW)+1*(row==nW),(col+1)*(col<nW)+1*(col==nW);
%             (row-1)*(row>1)+nW*(row==1),(col+1)*(col<nW)+1*(col==nW);
%             (row+1)*(row<nW)+1*(row==nW),(col-1)*(col>1)+nW*(col==1)];
%     row = neighboor(:,1)';
%     col = neighboor(:,2)';
%     neighboor = sub2ind([nW nW],row,col);
%     voisins(ii,:) = neighboor;
%     ww = [ww ; ones(length(neighboor),1)];
%     iw = [iw ; ii*ones(length(neighboor),1)];
%     jw = [jw ; neighboor'];
% end
% W = sparse(iw,jw,ww);
% SIM_1 = sum(W.*W,'all');
% SIM_2 = sum(W,'all').^2;
% SIM_3 = sum(sum(W,2).^2);
% SIM = sqrt( nW^2.*SIM_1+3*SIM_2-nW*SIM_3./((nW^2-1)*SIM_2));
% 
% ZM_H = IM_H./SIM;
% ZM_S = IM_S./SIM;

%% POST PROCESS INTERSPECIFIC Correlation coefficients
load('Corr_Host_Symbiont.mat')

[nsimu,nt] = size(C_HS(:,2:end)); 

tt = 100*(1:nt);

C_HS_med = quantile(C_HS(:,2:end),[0.05,0.5,0.95]);
P_HS_med = quantile(P_HS(:,2:end),[0.05,0.5,0.95]);


Color = [0.8,0.8,0.8;...
         0.3,0.3,0.3;...
         0.05,0.05,0.05];

%% MORAN INDEX          
figure(1)
clf
hold on
plot(tt,IM_S_med(2,:),'--','color','k','linewidth',2)
plot(tt,IM_H_med(2,:),'-','color','k','linewidth',3)
plot(tt,C_HS_med(2,:),'-.','color','k','linewidth',1.5)

%%% Symbiont
ttt = [tt, fliplr(tt)];
inBetween = [IM_S_med(1,:), fliplr(IM_S_med(3,:))];
f = fill(ttt, inBetween,Color(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,IM_S_med(2,:),'--','color','k','linewidth',2)


%%% Host
ttt = [tt, fliplr(tt)];
inBetween = [IM_H_med(1,:), fliplr(IM_H_med(3,:))];
f = fill(ttt, inBetween,Color(2,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,IM_H_med(2,:),'-','color','k','linewidth',3)

%%% Symbiont/Host
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

%% MORAN INDEX          
figure(2)
clf
hold on
plot(tt,IG_S_med(2,:),'--','color','k','linewidth',2)
plot(tt,IG_H_med(2,:),'-','color','k','linewidth',3)
plot(tt,C_HS_med(2,:),'-.','color','k','linewidth',1.5)

%%% Symbiont
ttt = [tt, fliplr(tt)];
inBetween = [IG_S_med(1,:), fliplr(IG_S_med(3,:))];
f = fill(ttt, inBetween,Color(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,IG_S_med(2,:),'--','color','k','linewidth',2)


%%% Host
ttt = [tt, fliplr(tt)];
inBetween = [IG_H_med(1,:), fliplr(IG_H_med(3,:))];
f = fill(ttt, inBetween,Color(2,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(tt,IG_H_med(2,:),'-','color','k','linewidth',3)

%%% Symbiont/Host
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

