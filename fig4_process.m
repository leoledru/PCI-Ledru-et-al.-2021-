clear all; clc;

%% INTRASPECIFIC INDEX
space_size = 100;
neighbourhood = 8;
load('data/simulations_gamma02_d.mat'); out_persist = ~cellfun('isempty',out_mutualist);
%%
% Measurement only on simulations with transition
speciation_test = [];
for ii = 1:length(out_simulation)
    s = out_mutualist{ii};
    simu = out_simulation{ii};
    s(simu~=3) = NaN;
    speciation_test = [speciation_test ; sum(s>.475)/sum(~isnan(s))];
end
index = find(speciation_test>.15);

voisins = zeros(space_size^2,neighbourhood);
for ii=1:space_size^2
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
end

% Loop over all persistent simulation
for jj = 1:length(index)   
    simulation = out_simulation{index(jj)};
    Host = out_host{index(jj)};
    Symbiont = out_mutualist{index(jj)};
    tmax = size(simulation,2);
    
    % pre-allocation of the matrices that will contain the local trait means
    % and the difference between global and local means
    A_H = zeros(tmax*space_size,space_size); % host
    A_S = zeros(tmax*space_size,space_size); % symbiont
    R_H = zeros(tmax*space_size,space_size); % random host 
    R_S = zeros(tmax*space_size,space_size); % random symbiont
    A_H_diff = zeros(tmax*space_size,space_size);
    A_S_diff = zeros(tmax*space_size,space_size);

    % pre-allocation of the matrices that will contain the global average value
    % of trait 
    Globalmean_H = zeros(tmax,1); % host
    Globalmean_S = zeros(tmax,1); % symbiont
    Globalmean_R_H = zeros(tmax,1); % random host
    Globalmean_R_S = zeros(tmax,1); % random symbiont 

    % pre-allocation of vectors that will contain the pairing index values
    r_h = zeros(tmax,1);
    r_s = zeros(tmax,1);

    for tt = 1:size(simulation,2)
        if mod(tt,1)==0
            simu = simulation(:,tt);
            host = Host(:,tt);
            host(simu==0) = NaN;
            symbiont = Symbiont(:,tt);
            symbiont(simu~=3) = NaN;
%             Random matrices with alpha values of matrices obtained by 
%             reorganized simulations on random positions keeping the same empty cells
            random_matrix_s = zeros(space_size)*NaN;
            random_matrix_h = zeros(space_size)*NaN; 
            rnd_alpha_s = symbiont(simu==3);
            rnd_alpha_h = host(simu>0);    
            coord_alpha_s = find(simu==3);
            coord_alpha_h = find(simu>0);
            coord_alpha_s = coord_alpha_s(randperm(length(coord_alpha_s)));
            coord_alpha_h = coord_alpha_h(randperm(length(coord_alpha_h)));
            random_matrix_s(coord_alpha_s) = rnd_alpha_s;
            random_matrix_h(coord_alpha_h) = rnd_alpha_h;
            
            R_h = zeros(space_size);
            R_s = zeros(space_size);
            A_h = zeros(space_size);
            A_s = zeros(space_size);
           
%             Loop over all cells
            for ii=1:space_size^2       
                neighboor = voisins(ii,:);
                if simu(ii)==0 % empty cell
                    % assigning NaN value to not take in account these values
                    % in the average (thanks to nanmean())
                    A_h(ii) = NaN;
                    A_s(ii) = NaN;
                    R_h(ii) = NaN;
                    R_s(ii) = NaN;
                end
                if simu(ii)==1 % host alone  
                    % coordinates of neighboring cells also containing a host
                    coord = neighboor.*(simu(neighboor)>0);
                    coord = coord(coord>0);
                    % compute effective assortment
                    mean_alpha_h = nanmean(host(coord)); % mean of neighborhood alpha
                    cpr_h = 1 - abs(host(ii)-mean_alpha_h); % 1 - (alpha of focus cell - mean alpha of neighborhood)
                    % compute random assortment
                    mean_alpha_r_h = nanmean(random_matrix_h(coord)); % same thing for random matrix
                    cpr_r_h = 1 - abs(random_matrix_h(ii)-mean_alpha_r_h);
                    A_h(ii) = cpr_h;
                    A_s(ii) = NaN; % it is a cell without symbiont, so NaN value for A_s to prevent to have 0
                    R_h(ii) = cpr_r_h; 
                    R_s(ii) = NaN;
                end
                if simu(ii)==3 % host-symbiont pair
                    % For symbiont : coordinates of neighboring cells also containing a host-symbiont pair
                    coord = neighboor.*(simu(neighboor)==3)';
                    coord = coord(coord>0);
                    % compute effective assortment
                    mean_alpha_s = nanmean(symbiont(coord));
                    cpr_s = 1 - abs(symbiont(ii)-mean_alpha_s);
                    % compute random assortment
                    mean_alpha_r_s = nanmean(random_matrix_s(coord));
                    cpr_r_s = 1 - abs(random_matrix_s(ii)-mean_alpha_r_s);
                    % For host : coordinates of neighboring cells containing a host-symbiont pair, or a host alone
                    coord = neighboor.*(simu(neighboor)>0);
                    coord = coord(coord>0);
                    % effective assortment
                    mean_alpha_h = nanmean(host(coord));
                    cpr_h = 1 - abs(host(ii)-mean_alpha_h); 
                    % random assortment
                    mean_alpha_r_h = nanmean(random_matrix_h(coord));
                    cpr_r_h = 1 - abs(random_matrix_h(ii)-mean_alpha_r_h);

                    A_h(ii) = cpr_h;
                    A_s(ii) = cpr_s;
                    R_h(ii) = cpr_r_h;
                    R_s(ii) = cpr_r_s;
                end
           end
           % compute the global mean of comparisons over all the landscape (nanmean)
           globalmean_h = nanmean(A_h,'all');
           globalmean_s = nanmean(A_s,'all');
           globalmean_r_h = nanmean(R_h,'all');
           globalmean_r_s = nanmean(R_s,'all');

           % Difference between local means and global mean
           A_h_diff = abs(A_h-globalmean_h);
           A_s_diff = abs(A_s-globalmean_s);

           A_H((tt-1)*space_size+1:tt*space_size,:) = A_h;
           A_S((tt-1)*space_size+1:tt*space_size,:) = A_s;
           R_H((tt-1)*space_size+1:tt*space_size,:) = R_h;
           R_S((tt-1)*space_size+1:tt*space_size,:) = R_s;
           A_H_diff((tt-1)*space_size+1:tt*space_size,:) = A_h_diff;
           A_S_diff((tt-1)*space_size+1:tt*space_size,:) = A_s_diff;

           Globalmean_H(tt) = globalmean_h;
           Globalmean_S(tt) = globalmean_s;   
           Globalmean_R_H(tt) = globalmean_r_h;
           Globalmean_R_S(tt) = globalmean_r_s;

           % Compute the difference of global mean between real and random
           % >0 = positive pairing / 0 = random / <0 = negative pairing
           diff_h = globalmean_h-globalmean_r_h;
           diff_s = globalmean_s-globalmean_r_s;

           r_h(tt) = diff_h;
           r_s(tt) = diff_s;
        end  
    end
    % store the assortment index in function of time
    assortment_h(jj,:) = r_h;
    assortment_s(jj,:) = r_s;
end
% save(['assortment','.mat'],'assortment_h','assortment_m','-v7.3')


%% INTERSPECIFIC INDEX

clear all; clc;
load('data/simulations_gamma02_d.mat'); out_persist = ~cellfun('isempty',out_mutualist);
space_size = 100;
neighbourhood = 8;
tmax = 10000;
%%
speciation_test = []; 
for ii = 1:1000
    s = out_mutualist{ii};
    simu = out_simulation{ii};
    s(simu~=3) = NaN;
    speciation_test = [speciation_test ; sum(s>.475)/sum(~isnan(s))];
end
index = find(speciation_test>.15);

% Loop over all persistent simulation
for jj = 1:50 
    simulation = out_simulation{index(jj)};
    Host = out_host{index(jj)};
    Symbiont = out_mutualist{index(jj)};
    tmax = size(simulation,2);
    % pre-allocation of the matrices that will contain the local trait means
    % and the difference between global and local means
    A_H = zeros(tmax*space_size,space_size); % host
    A_S = zeros(tmax*space_size,space_size); % symbiont 
    R_H = zeros(tmax*space_size,space_size); % random host 
    R_S = zeros(tmax*space_size,space_size); % random symbiont 
    A_H_diff = zeros(tmax*space_size,space_size);
    A_S_diff = zeros(tmax*space_size,space_size);
    % pre-allocation of the matrices that will contain the global average value
    % of trait 
    Globalmean_H = zeros(tmax,1); % host
    Globalmean_S = zeros(tmax,1); % symbiont
    Globalmean_R_H = zeros(tmax,1); % random host
    Globalmean_R_S = zeros(tmax,1); % random symbiont 
    % pre-allocation of vectors that will contain the pairing index values
    r = zeros(tmax,1);

    for tt = 1:size(simulation,2)
        if mod(tt,1)==0
            simu = simulation(:,tt);
            host = Host(:,tt);
            host(simu==0) = NaN;
            symbiont = Symbiont(:,tt);
            symbiont(simu~=3) = NaN;
            random_matrix_s = zeros(space_size)*NaN;
            random_matrix_h = zeros(space_size)*NaN;
            rnd_alpha_s = symbiont(simu==3);
            rnd_alpha_h = host(simu>0);
            coord_alpha_s = find(simu==3);
            coord_alpha_h = find(simu>0);
            coord_alpha_s = coord_alpha_s(randperm(length(coord_alpha_s)));
            coord_alpha_h = coord_alpha_h(randperm(length(coord_alpha_h)));
            random_matrix_s(coord_alpha_s) = rnd_alpha_s;
            random_matrix_h(coord_alpha_h) = rnd_alpha_h;
            
            assort = zeros(space_size);
            assort_rand = zeros(space_size);

            for ii=1:space_size^2       
                if simu(ii)==0
                    % assigning NaN value to not take in account these values
                    % in the average (thanks to nanmean())
                    assort(ii) = NaN;
                    assort_rand(ii) = NaN;
                end
                if simu(ii)==3
                    diff = 1 - abs(symbiont(ii)-host(ii));
                    diff_rand = 1 - abs(random_matrix_s(ii)-random_matrix_h(ii));
                    assort(ii) = diff;
                    assort_rand(ii) = diff_rand;
                end
            end
           % compute the global mean of comparisons over all the landscape (nanmean)
           globalmean = nanmean(assort,'all');
           globalmean_rand = nanmean(assort_rand,'all');
           r(tt) = globalmean-globalmean_rand;
        end  
    end
    % store the assortment index
    assortment(jj,:) = r;
end
save(['assortment_interspe','.mat'],'assortment','-v7.3')

