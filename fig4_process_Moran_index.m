clear all; clc;

%% INTRASPECIFIC INDEX
space_size = 100;
neighbourhood = 8;
load('data/simulations_gamma02_d.mat'); out_persist = ~cellfun('isempty',out_mutualist);
%%
% Measurement only on simulations with transition
speciation_test = [];
parfor ii = 1:length(out_simulation)
    s = out_mutualist{ii};
    simu = out_simulation{ii};
    s(simu~=3) = NaN;
    speciation_test = [speciation_test ; sum(s>.475)/sum(~isnan(s))];
end
index = find(speciation_test>.15);
% out_host_persist = out_host(index);
% out_host_d_persist = out_host_d(index);
% out_host_density_persist = out_host_density(index);
% out_mutualist_persist = out_mutualist(index);
% out_mutualist_d_persist = out_mutualist_d(index);
% out_simulation_persist = out_simulation(index);
% out_symbiont_density_persist = out_symbiont_density(index);
%
% save('data/simulations_gamma02_d_persist.mat','out_host_d_persist','out_host_density_persist',...
%     'out_host_persist','out_mutualist_d_persist','out_mutualist_persist','out_simulation_persist','out_symbiont_density_persist')

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

% Loop over all persistent simulation
IL_H = [];  % Leo index
IL_S = [];
IM_H = []; % Moran index global
IM_S = [];
seuil_H = []; % Mean of Moran ondex under random assumption
seuil_S = [];
IG_H = [];    % Geary index local
IG_S = [];
parfor jj = 1:length(index)
    simulation = out_simulation{index(jj)};
    Host = out_host{index(jj)};
    Symbiont = out_mutualist{index(jj)};
    tmax = size(simulation,2);
    
    % Moran index for hosts and symbionts
    IM_Ht = []; % host
    IM_St = []; % symbiont
    IL_Ht = []; % host
    IL_St = []; % symbiont
    seuil_host_t = [];
    seuil_symbiont_t = [];
    IG_Ht = []; % host
    IG_St = []; % symbiont
    for tt = 1:size(simulation,2)
        if mod(tt,1)==0
            simu = simulation(:,tt);
            host = Host(:,tt);
            host_simu = host(simu~=0);
            host_simu_mean = mean(host_simu);
            W_host_simu = W((simu~=0)',:);
            W_host_simu = W_host_simu(:,(simu~=0));
            symbiont = Symbiont(:,tt);
            symbiont_simu = symbiont(simu==3);
            symbiont_simu_mean = mean(symbiont_simu);
            W_symbiont_simu = W((simu==3)',:);
            W_symbiont_simu = W_symbiont_simu(:,(simu==3));
            
            %%% MORAN INDEX
            I_host = (host_simu-host_simu_mean)'*W_host_simu*(host_simu-host_simu_mean)...
                ./var(host_simu)./sum(W_host_simu,'all');
            I_symbiont = (symbiont_simu-symbiont_simu_mean)'*W_symbiont_simu*(symbiont_simu-symbiont_simu_mean)...
                ./var(symbiont_simu)./sum(W_symbiont_simu,'all');
            IM_Ht = [IM_Ht,I_host];
            IM_St = [IM_St,I_symbiont];
            
            %%% GEARY INDEX
            IG_host = ((sum(W_host_simu)+sum(W_host_simu,2)')*(host_simu.*host_simu) ...
                      - 2*host_simu'*W_host_simu*host_simu)...
                      ./var(host_simu)./sum(W_host_simu,'all')/2;
            IG_symbiont = ((sum(W_symbiont_simu)+sum(W_symbiont_simu,2)')*(symbiont_simu.*symbiont_simu) ...
                      - 2*symbiont_simu'*W_symbiont_simu*symbiont_simu)...
                      ./var(symbiont_simu)./sum(W_symbiont_simu,'all')/2;
            IG_Ht = [IG_Ht,IG_host];
            IG_St = [IG_St,IG_symbiont];
            
            %%% LEO INDEX
            sW_host = full(sum(W_host_simu,2));        
            IL_host = sum((sW_host~=0).*(host_simu - W_host_simu*host_simu./(sW_host + (sW_host==0)) ).^2)...
                      ./var(host_simu)./(sum(sW_host~=0)-1);
            
            sW_symbiont = full(sum(W_symbiont_simu,2));
            IL_symbiont = sum((sW_symbiont~=0).*(symbiont_simu - W_symbiont_simu*symbiont_simu./(sW_symbiont + (sW_symbiont==0)) ).^2)...
                      ./var(symbiont_simu)./(sum(sW_symbiont~=0)-1);

            IL_Ht = [IL_Ht,IL_host];
            IL_St = [IL_St,IL_symbiont];
            
            %%%% POP SIZE
            sht = 1./(length(host_simu)-1);
            sst = 1./(length(symbiont_simu)-1);
            seuil_host_t = [seuil_host_t,sht];
            seuil_symbiont_t = [seuil_symbiont_t,sst];
            
        end
    end
    IM_H = [IM_H;IM_Ht];
    IM_S = [IM_S;IM_St];
    seuil_H = [seuil_H;seuil_host_t];
    seuil_S = [seuil_S;seuil_symbiont_t];
    IG_H = [IG_H;IG_Ht];
    IG_S = [IG_S;IG_St];
    IL_H = [IL_H;IL_Ht];
    IL_S = [IL_S;IL_St];
end
save(['Index_Moran','.mat'],'IM_H','IM_S','seuil_H','seuil_S','IG_H','IG_S','IL_H','IL_S','-v7.3')

%% INTERSPECIFIC INDEX of correlation
load('data/simulations_gamma02_d_persist.mat')

% Loop over all persistent simulation
C_HS = [];
P_HS = [];
parfor jj = 1:length(out_host_persist)
    simulation = out_simulation_persist{jj};
    Host = out_host_persist{jj};
    Symbiont = out_mutualist_persist{jj};
    tmax = size(simulation,2);
    
    % Correlation coefficient
    C_HSt = []; % host/symbiont
    P_HSt = [];
    for tt = 1:size(simulation,2)
        simu = simulation(:,tt);
        host = Host(:,tt);
        host_simu = host(simu==3);
        host_simu_mean = mean(host_simu);

        symbiont = Symbiont(:,tt);
        symbiont_simu = symbiont(simu==3);
        symbiont_simu_mean = mean(symbiont_simu);
        
        [c_hst,p_hst] = corrcoef(symbiont_simu,host_simu); 
        C_HSt = [C_HSt,c_hst(1,2)];
        P_HSt = [P_HSt,p_hst(1,2)];
    end
    C_HS = [C_HS;C_HSt];
    P_HS = [P_HS;P_HSt];
end
save(['Corr_Host_Symbiont','.mat'],'C_HS','P_HS','-v7.3')
