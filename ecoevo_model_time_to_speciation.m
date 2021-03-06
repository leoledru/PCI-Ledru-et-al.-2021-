function [out] = ecoevo_model_time_to_speciation(tmax,total,host,symbiont,host_d,symbiont_d,space_size,...
                 neighbourhood, host_mortality, symbiont_mortality,...
                 host_density_dependent, host_germination, host_competition,...
                 symbiont_colonization, p_mut, muta_max, muta_min, muta_alpha, muta_epsilon_host,...
                 muta_epsilon_symbiont,symbiont_density_dependent,coeff_compet,gamma_symbiont,coeff_disp,coeff_disp_h,exp_mu)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fmut_min_s = 2.5;
    Fmut_min_h = .1;
    Fmut_max = 8;
    Fb_min_h = .7;
    Fb_min_s = .7;
    Fb_max = 1;
    gamma_1_s = 1;
    gamma_1_h = 4;
    gamma_2 = 1;
    gamma_h = 1;
    f_mutu_s = @(alpha_h) Fmut_min_s+(Fmut_max-Fmut_min_s).*alpha_h.^gamma_1_s;
    f_mutu_h = @(alpha_s) Fmut_min_h+(Fmut_max-Fmut_min_h).*alpha_s.^gamma_1_h;
    f_basal_s = @(alpha_s) Fb_max -(Fb_max-Fb_min_s).*alpha_s.^gamma_2;
    f_basal_h = @(alpha_h) Fb_max-(Fb_max-Fb_min_h).*alpha_h.^gamma_2; 
    fhseul_max = .5;
    fhseul_min = fhseul_max*Fb_min_h;
    fitness_hseul=@(alpha_h)fhseul_max-(fhseul_max-fhseul_min).*alpha_h.^gamma_h; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compuation of neighbourhood coordinates
    Neighboor = [];
    for ii = 1:space_size^2
      [row,col] = ind2sub([space_size space_size],ii); 
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
      Neighboor = [Neighboor;neighboor];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    total_old = total;
    total_new = total_old;
    test = 0;
    time = 0;
    
    while test==0 % run until speciation or time = tmax
        % Host alone mortality 
        h = sum(total_old==1,'all');
        proba = rand(h,1);
        total_new(total_old==1) = 1.*(proba>host_mortality) + 0.*(proba<host_mortality);

        % Host with symbiont mortality
        hm = sum(total_old==3,'all');
        proba = rand(hm,1);
        total_new(total_old==3) = 3.*(proba>host_mortality) + 0.*(proba<host_mortality);    

        % Symbiont mortality
        hm = sum(total_new==3,'all');
        proba = rand(hm,1);
        total_new(total_new==3) = 3.*(proba>symbiont_mortality) + 1.*(proba<symbiont_mortality);

        % update populations matrix after mortality
        symbiont(total_new==0) = 0;
        symbiont(total_new==1) = 0; 
        symbiont_d(total_new==0) = 0; 
        symbiont_d(total_new==1) = 0;
        host(total_new==0) = 0;
        host_d(total_new==0) = 0;

        % Reproduction, mutation, dispersal, colonization
        host_proportion_global = sum(total_new>0,'all')/(space_size^2); 
        symbiont_proportion_global = sum(total_new==3,'all')/(space_size^2);    
        % Host
        global_disp_h = zeros(space_size^2,1);
        local_disp_h = zeros(space_size^2,1);
        % fitness host alone
        number_h = sum(total_new==1,'all');
        fitness_h=fitness_hseul(host(total_new==1)).*ones(number_h,1);
        fitness_h = fitness_h.*(1 - host_d(total_new==1).*coeff_disp_h); % Apply dispersal cost on fitness
        fitness_host = fitness_h-mod(fitness_h,1)+binornd(1,mod(fitness_h,1)); % decimal part = parameter p of bernouilli trial
        global_disp_h(total_new==1) = binornd(fitness_host,host_d(total_new==1).*ones(size(fitness_host))); % number of offspring propagated by global dispersal
        local_disp_h(total_new==1) = fitness_host-global_disp_h(total_new==1); % number of offsprings propagated by local dispersal
        % fitness host with symbiont
        fitness_h= f_basal_h(host(total_new==3)).*f_mutu_h(symbiont(total_new==3));
        fitness_h = fitness_h.*(1 - host_d(total_new==3).*coeff_disp_h); % Apply dispersal cost on fitness
        fitness_host = fitness_h-mod(fitness_h,1)+binornd(1,mod(fitness_h,1));
        global_disp_h(total_new==3) = binornd(fitness_host,host_d(total_new==3).*ones(size(fitness_host))); % number of offsprings propagated by global dispersal
        local_disp_h(total_new==3) = fitness_host-global_disp_h(total_new==3); % number of offsprings propagated by local dispersal

        % Dispersal
        if sum(global_disp_h)+sum(local_disp_h)>0
            % arrival coordinates from global dispersal 
            disp_global = randi([1 space_size^2],sum(global_disp_h),1); % randomly choose cells for each offspring
            % arrival coordinates from local dispersal
            effective_neighboors = Neighboor.*(local_disp_h>0); % only coordinates of neighborhood of host which have local offspring
            effective_neighboors(~any(effective_neighboors,2),:) = [];
            disp_local = [];
            effective_local = local_disp_h(any(local_disp_h,2),:);
            for ii = 1:sum(local_disp_h>0)
                disp_local = [disp_local; datasample(effective_neighboors(ii,:),effective_local(ii))']; % select a random neihboors for each offspring 
            end
            disp = [disp_global;disp_local];
            germination = unique(disp); % only unique localisation which received offsprings
            proba = rand(length(germination),1);
            host_proportion_local = zeros(length(germination),1);
            for ii = 1:length(germination)
                jj = germination(ii);
                host_proportion_local(ii) = sum(total_new(Neighboor(jj,:))>0,'all')./neighbourhood;
            end
            host_proportion = host_density_dependent.*host_proportion_local+(1-host_density_dependent).*host_proportion_global;
            germination = germination.*(proba<(host_germination-host_competition(host_proportion))); % comment for fixed competition
    %         germination = germination.*(proba<(host_germination-host_competition)); % uncomment for fixed competition
            germination(germination==0) = [];
            germination = germination.*(total_new(germination)==0);
            germination(germination==0) = []; % only keeps the positions of successful germinations       
            total_new(germination) = 1;

            index = [repelem(host,global_disp_h);repelem(host,local_disp_h)]; % find the correspondance between host and localisation of germination
            index_d = [repelem(host_d,global_disp_h);repelem(host_d,local_disp_h)];
            for ii = 1:length(germination)
                if sum(disp==germination(ii))>1 % several offsprings at the same position : uniform lottery
                    participants_1 = index(disp==germination(ii))'; % find the host which disperse their offspring at the position germination(ii)
                    participants_2 = index_d(disp==germination(ii))'; % same thing but for trait epsilon
                    participants = [participants_1;participants_2];
                    lottery = participants(:,randperm(size(participants,2))); % random permutation of columns, first column win 
                    host(germination(ii)) = lottery(1,1);
                    host_d(germination(ii)) = lottery(2,1);
                else
                    host(germination(ii)) = index(disp==germination(ii)); % if only one offspring -> germination
                    host_d(germination(ii)) = index_d(disp==germination(ii));
                end
            end
           % Mutation
           if muta_alpha == 1
               % alpha
               proba = rand(size(germination));
               direction = randi([0 1],size(germination)); % random direction
               direction(direction==0) = -1;
               magnitude = (muta_max-muta_min).*exprnd(exp_mu,length(germination),1) + muta_min; % random mutation magnitude : EXP
               magnitude(magnitude>muta_max) = 0; % when a mutation exceed the muta_max it is cancelled           
               host(germination) = (proba<p_mut).*magnitude.*direction + host(germination); % set the new trait value for host-offspring (if mutation)
               % cancel mutation if out of trait domain
               host(germination(host(germination)>1)) = host(germination(host(germination)>1))-1.*magnitude(host(germination)>1);
               host(germination(host(germination)<0)) = host(germination(host(germination)<0))+1.*magnitude(host(germination)<0);  
           end
           if muta_epsilon_host == 1
               % epsilon
               proba = rand(size(germination));
               direction = randi([0 1],size(germination)); % random direction
               direction(direction==0) = -1;
               magnitude = (muta_max-muta_min).*exprnd(exp_mu,length(germination),1) + muta_min; % random mutation magnitude : EXP
               magnitude(magnitude>muta_max) = 0; % when a mutation exceed the muta_max it is cancelled    
               host_d(germination) = (proba<p_mut).*magnitude.*direction + host_d(germination); % set the new trait value for host-offspring (if mutation)
               % cancel mutation if out of trait domain
               host_d(germination(host_d(germination)>1)) = host_d(germination(host_d(germination)>1))-1.*magnitude(host_d(germination)>1);
               host_d(germination(host_d(germination)<0)) = host_d(germination(host_d(germination)<0))+1.*magnitude(host_d(germination)<0);   
           end
        end

        % Symbiont
        global_disp_s = zeros(space_size^2,1);
        local_disp_s = zeros(space_size^2,1);
        % fitness symbiont with host
        number_s = sum(total_new==3,'all');
        fitness_s = f_basal_s(symbiont(total_new==3)).*f_mutu_s(host(total_new==3));    
        fitness_symbiont = fitness_s-mod(fitness_s,1) + binornd(1,mod(fitness_s,1));
        % Competition on fitness for symbiont
        index = find(total_new==3);
        symbiont_proportion_local = zeros(number_s,1);
        for ii = 1:number_s
            jj = index(ii);
            symbiont_proportion_local(ii) = sum(total_new(Neighboor(jj,:))==3,'all')./neighbourhood;
        end
        symbiont_proportion = symbiont_density_dependent.*symbiont_proportion_local+(1-symbiont_density_dependent).*symbiont_proportion_global;

        fitness_symbiont = fitness_symbiont.*(1 - (symbiont_proportion.*coeff_compet).^gamma_symbiont); % Apply competition on fitness
        fitness_symbiont = fitness_symbiont.*(1 - symbiont_d(total_new==3).*coeff_disp); % Apply dispersal cost on fitness
        fitness_symbiont = fitness_symbiont-mod(fitness_symbiont,1)+binornd(1,mod(fitness_symbiont,1));

        global_disp_s(total_new==3) = binornd(fitness_symbiont,symbiont_d(total_new==3).*ones(size(fitness_symbiont))); % number of offsprings propagated by global dispersal
        local_disp_s(total_new==3) = fitness_symbiont-global_disp_s(total_new==3); % number of offsprings propagated by local dispersal
        % Dispersal
        if sum(global_disp_s)+sum(local_disp_s)>0
            % arrival coordinates from global dispersal 
            disp_global = randi([1 space_size^2],sum(global_disp_s),1); % randomly choose cells for each offspring
            % arrival coordinates from local dispersal
            effective_neighboors = Neighboor.*(local_disp_s>0); % only coordinates of neighborhood of host which have local offspring
            effective_neighboors(~any(effective_neighboors,2),:) = [];
            disp_local = [];
            effective_local = local_disp_s(any(local_disp_s,2),:);
            for ii = 1:sum(local_disp_s>0)
                disp_local = [disp_local; datasample(effective_neighboors(ii,:),effective_local(ii))']; % select a random neihboors for each offspring 
            end
            disp = [disp_global;disp_local];
            colonization = unique(disp); % only unique localisation which received offsprings
            proba = rand(length(colonization),1);
            colonization = colonization.*(proba<symbiont_colonization);
            colonization(colonization==0) = [];
            colonization = colonization.*(total_new(colonization)==1);
            colonization(colonization==0) = []; % only keeps the positions of successful germinations       
            total_new(colonization) = 3; 
            index = [repelem(symbiont,global_disp_s);repelem(symbiont,local_disp_s)]; % allow to find the correspondance between host and localisation of germination
            index_d = [repelem(symbiont_d,global_disp_s);repelem(symbiont_d,local_disp_s)];
            for ii = 1:length(colonization)
                if sum(disp==colonization(ii))>1 % several offsprings at the same position : uniform lottery
                    participants_1 = index(disp==colonization(ii))'; % find the host which disperse their offspring at the position germination(ii)
                    participants_2 = index_d(disp==colonization(ii))'; % same thing but for trait epsilon
                    participants = [participants_1;participants_2];
                    lottery = participants(:,randperm(size(participants,2))); % random permutation of columns, first column win 
                    symbiont(colonization(ii)) = lottery(1,1);
                    symbiont_d(colonization(ii)) = lottery(2,1);
                else
                    symbiont(colonization(ii)) = index(disp==colonization(ii)); % if only one offspring -> germination
                    symbiont_d(colonization(ii)) = index_d(disp==colonization(ii));
                end
            end
           % Mutation
           if muta_alpha==1
               % alpha
               proba = rand(size(colonization));
               direction = randi([0 1],size(colonization)); % random direction
               direction(direction==0) = -1;
               magnitude = (muta_max-muta_min).*exprnd(exp_mu,length(colonization),1) + muta_min; % random mutation magnitude : EXP
               magnitude(magnitude>muta_max) = 0; % when a mutation exceed the muta_max it is cancelled    
               symbiont(colonization) = (proba<p_mut).*magnitude.*direction + symbiont(colonization); % set the new trait value for host-offspring (if mutation)
               % cancel mutation if out of trait domain
               symbiont(colonization(symbiont(colonization)>1)) = symbiont(colonization(symbiont(colonization)>1))-1.*magnitude(symbiont(colonization)>1);
               symbiont(colonization(symbiont(colonization)<0)) = symbiont(colonization(symbiont(colonization)<0))+1.*magnitude(symbiont(colonization)<0);  
           end
           if muta_epsilon_symbiont == 1
               % epsilon
               proba = rand(size(colonization));
               direction = randi([0 1],size(colonization)); % random direction
               direction(direction==0) = -1;
               magnitude = (muta_max-muta_min).*exprnd(exp_mu,length(colonization),1) + muta_min; % random mutation magnitude : EXP
               magnitude(magnitude>muta_max) = 0; % when a mutation exceed the muta_max it is cancelled    
               symbiont_d(colonization) = (proba<p_mut).*magnitude.*direction + symbiont_d(colonization); % set the new trait value for host-offspring (if mutation)
               % cancel mutation if out of trait domain
               symbiont_d(colonization(symbiont_d(colonization)>1)) = symbiont_d(colonization(symbiont_d(colonization)>1))-1.*magnitude(symbiont_d(colonization)>1);
               symbiont_d(colonization(symbiont_d(colonization)<0)) = symbiont_d(colonization(symbiont_d(colonization)<0))+1.*magnitude(symbiont_d(colonization)<0);   
           end
        end
        % update populations matrices
        total_old = total_new;
        time = time + 1;

        if sum(total_new==3,'all')==0 || sum(total_new,'all')==0
            if sum(total_new==3,'all')==0 
                collapse = 1; % symbiont collapse
            end
            if sum(total_new,'all')==0
                collapse = 2; % full system collapse
            end

            break
        end

        % speciation test
        speciation_test = sum(symbiont>.475)/sum(total_new==3);
        if speciation_test>.1
            time_speciation = time;
            test = 1;
        end
        % tmax test
        if time >= tmax
            time_speciation = NaN;
            test = 1;
        end

    end
out.time_speciation = time_speciation;
end
