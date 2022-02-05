clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Approximation model of mutualistic/parasitic system             %
%                    with or without  mutation                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model parameters
f_hmin  = 0.1;     % minimal interaction host fecundity
f_hmax  =  8;      % maximal interaction host fecundity
f_hamax = 0.5;     % maximal host alone fecundity
f_smax  = f_hmax;  % maximal interaction symbiont fecundity
f_smin  = 2.5;     % minimal interaction symbiont fecundity
gammaf  = 4;       % selection strength on symbiont interaction trait
c       = 0.3;     % mutualism cost
mh      = 0.06;    % host mortality
ms      = 0.06;    % symbiont mortality
gammac  = 0.2;     % Competition stregth : gammac<1 = strong competition // gammac>=1 = weak competition
d       = 0;       % dispersal cost
e       = 0.02;    % mutation rate between symbionts

%% Dynamical system
%%% without symbiont
f_ha =@(alpha_h,e_h) (1-d*e_h).*(1-c*alpha_h).*f_hamax;

%%% With parasitic/mutualistic symbiont
f_h =@(alpha_h,alpha_s,e_h) (1-d*e_h).*(1-c*alpha_h).*(f_hmin+(f_hmax-f_hmin).*alpha_s.^gammaf);
f_s =@(alpha_h,alpha_s,e_s) (1-d*e_s).*(1-c*alpha_s).*(f_smin+(f_smax-f_smin).*alpha_h);

F1 = @(alpha_h,alpha_s,e_h,e_s,X) ...
    (1-mh).*X(1).*( 1 - f_s(alpha_h,alpha_s(1),e_s).*X(3).^2./(X(3)+X(4)) ...
    - f_s(alpha_h,alpha_s(2),e_s).*X(4).^2./(X(3)+X(4)) ...
    + f_ha(alpha_h,e_h).*( 1-(X(1)+X(2)).^gammac ) ) ...
    +(1-mh).*X(2).*( 1-(X(1)+X(2)).^gammac).*( f_h(alpha_h,alpha_s(1),e_h).*X(3)./(X(3)+X(4))...
    +f_h(alpha_h,alpha_s(2),e_h).*X(4)./(X(3)+X(4)));
F2 = @(alpha_h,alpha_s,e_h,e_s,X) (1-mh).*X(2) ...
    + (1-mh).*X(1).*( f_s(alpha_h,alpha_s(1),e_s).*X(3).^2./(X(3)+X(4)) ...
    + f_s(alpha_h,alpha_s(2),e_s).*X(4).^2./(X(3)+X(4)) );
F3 = @(alpha_h,alpha_s,e_h,e_s,X) (1-ms).*X(3).*(1+ (1-e).*f_s(alpha_h,alpha_s(1),e_s).*X(1).*X(3)./(X(3)+X(4)))...
                                  + (1-ms).*X(4).*e.*f_s(alpha_h,alpha_s(2),e_s).*X(1).*X(4)./(X(3)+X(4)) ;
F4 = @(alpha_h,alpha_s,e_h,e_s,X) (1-ms).*X(4).*(1+ (1-e).*f_s(alpha_h,alpha_s(2),e_s).*X(1).*X(4)./(X(3)+X(4)))...
                                  + (1-ms).*X(3).*e.*f_s(alpha_h,alpha_s(1),e_s).*X(1).*X(3)./(X(3)+X(4));                       
                              
F = @(alpha_h,alpha_s,e_h,e_s,X) [max(F1(alpha_h,alpha_s,e_h,e_s,X),0); ...
    F2(alpha_h,alpha_s,e_h,e_s,X); ...
    F3(alpha_h,alpha_s,e_h,e_s,X); ...
    F4(alpha_h,alpha_s,e_h,e_s,X)];

%% Traits of hosts and symbionts
alpha_h = 0;        % Host interaction trait
alpha_s = [0,0.55]; % Symbionts interaction traits 
e_h = 1;            % Host dispersal trait
e_s = 1;            % Symbionts dispersal trait

%% Initial proportion of Host and symbionts
X0ha = 0.1;   % Host alone density 
X0hs = 0.01;  % Host with a symbiont density
choice = 3;   % 1 = more parasitic / 2 = more mutualistic / 3 = coexistence p= p*
astar  = ((f_hamax-f_hmin)./(f_hmax-f_hmin)).^(1./gammaf);
psp_star = f_s(alpha_h,alpha_s(2),e_s)./(f_s(alpha_h,alpha_s(2),e_s)+f_s(alpha_h,alpha_s(1),e_s));

if (choice==1)       % more parasitic symbionts
    psp = 0.85;
elseif (choice == 2) % more mutualistic symbionts
    psp  = 0.42;
elseif(choice==3)    % coexistence p =p*
    psp = psp_star;
end
X0sp = psp*X0hs;            % Parasitic symbiont density
X0sm = (1-psp)*X0hs;        % Mutualistic symbiont density
X0 = [X0ha;X0hs;X0sp;X0sm]; % Initial densities

%% Computation of the dynamics
Tf = 500;                  % Final time
t = 0;                     % Initial time
X = X0;                    % Initial densities
Xnew = X;
while(t<Tf)
    Xold = Xnew;
    Xnew = F(alpha_h,alpha_s,e_h,e_s,Xold);
    t = t+1;
    X = [X,Xnew];
end

X_h = sum(X(1:2,:));    % Host density
X_sp = X(3,:);          % Parasitic symbiont density
X_sm = X(4,:);          % Mutualistic symbiont density
P_s= X_sm./(X_sp+X_sm); % Proportion of mutualistic symbiont


%% Equilibrium
%%% without symbiont
rho_ha =@(alpha_h) (1-mh/(1-mh)./f_ha(alpha_h,e_h)).^(1./gammac);
R_ha = rho_ha(0);

%%% With parasitic/mutualistic symbiont
p_ha =@(alpha_h,alpha_s) mh./(1-mh)./f_s(alpha_h,alpha_s,e_s);
rho_h =@(alpha_h,alpha_s,p) p.*(f_h(alpha_h,alpha_s,e_s).*(1-p.^gammac)-mh) ...
    + (1-mh).*p_ha(alpha_h,alpha_s).*(f_ha(alpha_h,e_h)-f_h(alpha_h,alpha_s,e_h)).*(1-p.^gammac);

%%% Parasitic system
R_hparasitic = fsolve(@(p) rho_h(0,0,p),1);
R_sparasitic = R_hparasitic - p_ha(0,0);

%%% Mutualistic system
R_hmutualistic = fsolve(@(p) rho_h(0,1,p),1);
R_smutualistic = R_hmutualistic - p_ha(0,1);

%% Transition threshold 
rho_ha_star = X(2,end);%0.15;
px = @(X) (1-e)*f_s(alpha_h,alpha_s(2),0)*rho_ha_star.*X.^2 ...
        + e*f_s(alpha_h,alpha_s(1),0)*rho_ha_star*(1-X).^2 - mh/(1-mh)*X;
pstar = fsolve(@(X) px(X),0.1);

%% Figures
Color = get(gca,'colororder');
Marker = ['o','*','d','^','v'];
tt = 0:Tf;

%%% Densities over time
figure(1)
clf
%%% Densities
plot(tt,X_h,'color',Color(5,:),'linewidth',2)
hold on
plot(tt,X_sp,'color',Color(1,:),'linewidth',2)
plot(tt,X_sm,'color',Color(2,:),'linewidth',2)
plot(tt,P_s,'--','color',Color(4,:),'linewidth',1.5)
%%% Equilbrium
if (choice==1)
    line([0,Tf],[R_hparasitic,R_hparasitic],'color',Color(5,:))
    line([0,Tf],[R_sparasitic,R_sparasitic],'color',Color(2,:))
else
    line([0,Tf],[R_hmutualistic,R_hmutualistic],'color',Color(5,:))
    line([0,Tf],[R_smutualistic,R_smutualistic],'color',Color(1,:))
end

xlabel('Time $t$','Interpreter','latex','FontSize',16)
ylabel('proportion of host and symbionts','Interpreter','latex','FontSize',16)
legend('Hosts','Parasitic symbionts','Mutualistic symbionts','Relative proportion of mutualistic symbiont','Interpreter','latex','FontSize',14)
axis([0,Tf,0,0.35])







