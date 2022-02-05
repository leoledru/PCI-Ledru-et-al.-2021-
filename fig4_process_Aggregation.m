clear all; clc;

%% INTRASPECIFIC INDEX
space_size = 100;
neighbourhood = 8;

load('data/simulations_gamma02_d_persist.mat')

voisins = zeros(space_size^2,neighbourhood);
iw = [];
jw = [];
ww = [];
parfor ii=1:space_size^2
    [row,col] = ind2sub([space_size space_size],ii); % coordinates IJ of the cell ii
    if neighbourhood==8
        neighboor = [(row-1)*(row>1)+space_size*(row==1),col; % coordinates IJ of the 8 neighboors (wrap space)
            (row+1)*(row<space_size)+1*(row==space_size),col;
            row,(col-1)*(col>1)+space_size*(col==1);
            row,(col+1)*(col<space_size)+1*(col==space_size);
            (row-1)*(row>1)+space_size*(row==1),(col-1)*(col>1)+space_size*(col==1);
            (row+1)*(row<space_size)+1*(row==space_size),(col+1)*(col<space_size)+1*(col==space_size);
            (row-1)*(row>1)+space_size*(row==1),(col+1)*(col<space_size)+1*(col==space_size);
            (row+1)*(row<space_size)+1*(row==space_size),(col-1)*(col>1)+space_size*(col==1)];
    end
    if neighbourhood==24
        neighboor = [(row-1)*(row>1)+space_size*(row==1),col; % coordinates IJ of the extended neighbourhood, 24 neighboors (wrap space)
            (row-2)*(row>2)+(space_size-1)*(row==1)+(space_size)*(row==2),col;
            (row+1)*(row<space_size)+1*(row==space_size),col;
            (row+2)*(row<space_size-1)+2*(row==space_size)+1*(row==space_size-1),col;
            row,(col-1)*(col>1)+space_size*(col==1);
            row,(col-2)*(col>2)+(space_size)*(col==2)+(space_size-1)*(col==1);
            row,(col+1)*(col<space_size)+1*(col==space_size);
            row,(col+2)*(col<space_size-1)+1*(col==space_size-1)+2*(col==space_size);
            (row-1)*(row>1)+space_size*(row==1),(col-1)*(col>1)+space_size*(col==1);
            (row+1)*(row<space_size)+1*(row==space_size),(col+1)*(col<space_size)+1*(col==space_size);
            (row-1)*(row>1)+space_size*(row==1),(col+1)*(col<space_size)+1*(col==space_size);
            (row+1)*(row<space_size)+1*(row==space_size),(col-1)*(col>1)+space_size*(col==1);
            (row-1)*(row>1)+space_size*(row==1),(col-2)*(col>2)+(space_size)*(col==2)+(space_size-1)*(col==1);
            (row-1)*(row>1)+space_size*(row==1),(col+2)*(col<space_size-1)+1*(col==space_size-1)+2*(col==space_size);
            (row+1)*(row<space_size)+1*(row==space_size),(col-2)*(col>2)+(space_size)*(col==2)+(space_size-1)*(col==1);
            (row+1)*(row<space_size)+1*(row==space_size),(col+2)*(col<space_size-1)+1*(col==space_size-1)+2*(col==space_size);
            (row-2)*(row>2)+(space_size-1)*(row==1)+(space_size)*(row==2),(col-1)*(col>1)+space_size*(col==1);
            (row-2)*(row>2)+(space_size-1)*(row==1)+(space_size)*(row==2),(col-2)*(col>2)+(space_size)*(col==2)+(space_size-1)*(col==1);
            (row-2)*(row>2)+(space_size-1)*(row==1)+(space_size)*(row==2),(col+1)*(col<space_size)+1*(col==space_size);
            (row-2)*(row>2)+(space_size-1)*(row==1)+(space_size)*(row==2),(col+2)*(col<space_size-1)+1*(col==space_size-1)+2*(col==space_size);                                 (row+2)*(row<space_size-1)+2*(row==space_size)+1*(row==space_size-1),(col-1)*(col>1)+space_size*(col==1);
            (row+2)*(row<space_size-1)+2*(row==space_size)+1*(row==space_size-1),(col-2)*(col>2)+(space_size)*(col==2)+(space_size-1)*(col==1);
            (row+2)*(row<space_size-1)+2*(row==space_size)+1*(row==space_size-1),(col+1)*(col<space_size)+1*(col==space_size);
            (row+2)*(row<space_size-1)+2*(row==space_size)+1*(row==space_size-1),(col+2)*(col<space_size-1)+1*(col==space_size-1)+2*(col==space_size)];
    end
    row = neighboor(:,1)';
    col = neighboor(:,2)';
    neighboor = sub2ind([space_size space_size],row,col);
    voisins(ii,:) = neighboor;
    ww = [ww ; ones(length(neighboor),1)];
    iw = [iw ; ii*ones(length(neighboor),1)];
    jw = [jw ; neighboor'];
end
W = sparse(iw,jw,ww);

%% Loop to compute the aggregation index per simulation
% Aggregation index in (0,1)
A_H = [];
A_SM= [];
A_SP = [];

seuil_mut = 0.475;
parfor jj = 1:length(out_host_persist)
    simulation = out_simulation_persist{jj};
    Host = out_host_persist{jj};
    Symbiont = out_mutualist_persist{jj};
    tmax = size(simulation,2);
    % Aggregation index for hosts and symbionts
    A_Ht = []; % host
    A_SMt = []; % symbiont mutualist
    A_SPt = []; % symbiont parasit
    
    for tt = 1:size(simulation,2)
        %%% Host
        simu = simulation(:,tt);
        W_host_simu = W((simu~=0)',:);
        W_host_simu = W_host_simu(:,(simu~=0));
        
        n_host = sum(simu~=0);
        a_host_theo = 4*n_host-ceil(6*sqrt(n_host));
        a_host_simu = sum(W_host_simu,'all')/2;
        a_host = a_host_simu./a_host_theo;
%         a_host = mean(sum(W_host_simu,2));
        
        
        %%% Symbionts
        symbiont = Symbiont(:,tt);
        %%% Mutualist
        I_mut = logical((simu==3).*(symbiont>seuil_mut));
        W_symbiont_mut_simu = W(I_mut',:);
        W_symbiont_mut_simu = W_symbiont_mut_simu(:,I_mut);
        n_s_mut = sum(I_mut);
        
        a_s_mut_theo = 4*n_s_mut-floor(6*sqrt(n_s_mut));
        a_s_mut_simu = sum(W_symbiont_mut_simu,'all')/2;
        a_s_mut = a_s_mut_simu./a_s_mut_theo;        
%         a_s_mut = mean(sum(W_symbiont_mut_simu,2));
        
        %%% Parasit
        I_par = logical((simu==3).*(symbiont<=seuil_mut));
        W_symbiont_par_simu = W(I_par',:);
        W_symbiont_par_simu = W_symbiont_par_simu(:,I_par);
        n_s_par = sum(I_par);

        a_s_par_theo = 4*n_s_par-floor(6*sqrt(n_s_par));
        a_s_par_simu = sum(W_symbiont_par_simu,'all')/2;
        a_s_par = a_s_par_simu./a_s_par_theo;
%         a_s_par = mean(sum(W_symbiont_par_simu,2));
             
        
        A_Ht = [A_Ht,a_host];
        A_SMt = [A_SMt,a_s_mut];
        A_SPt = [A_SPt,a_s_par];
        
    end
    A_H = [A_H;A_Ht];
    A_SM = [A_SM;A_SMt];
    A_SP = [A_SP;A_SPt];
    
end

save(['Aggregation_index','.mat'],'A_H','A_SM','A_SP','-v7.3')

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
plot(tt,A_H_med(2,:),'-')
plot(tt_m,A_SM_med(2,:),'--')
plot(tt,A_SP_med(2,:),'-.')

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

