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
m       = 0.06;    % mortality probability
gammac  = 0.2;     % Competition stregth : gammac<1 = strong competition // gammac>=1 = weak competition
d       = 0;       % dispersal cost
e       = 0.008;    % mutation rate between symbionts

%% Dynamical system
%%% without symbiont
f_ha =@(alpha_h,e_h) (1-d*e_h).*(1-c*alpha_h).*f_hamax;

%%% With parasitic/mutualistic symbiont
f_h =@(alpha_h,alpha_s,e_h) (1-d*e_h).*(1-c*alpha_h).*(f_hmin+(f_hmax-f_hmin).*alpha_s.^gammaf);
f_s =@(alpha_h,alpha_s,e_s) (1-d*e_s).*(1-c*alpha_s).*(f_smin+(f_smax-f_smin).*alpha_h);

%%% Poisson distrib + space + choice of parasitic mutualistic
F1 = @(alpha_h,alpha_s,e_h,e_s,X) ...
    (1-m).*X(1).*exp(-(1-m)^2.*(f_s(alpha_h,alpha_s(1),e_s).*X(3) + f_s(alpha_h,alpha_s(2),e_s).*X(4))) ...
    +  (1-exp(-(1-m).*(f_h(alpha_h,alpha_s(1),e_h).*X(3) + f_h(alpha_h,alpha_s(2),e_h).*X(4) + f_ha(alpha_h,e_h).*X(1)) ) )...
    .*(1-(1-m)*(X(1)+X(2))).*( 1-((1-m)*(X(1)+X(2))).^gammac) + m.*(1-m).*X(2);

F2 = @(alpha_h,alpha_s,e_h,e_s,X) (1-m).*(1-m).*X(2) ...
    + (1-m).*X(1).*(1-exp(-(1-m)^2.*(f_s(alpha_h,alpha_s(1),e_s).*X(3) + f_s(alpha_h,alpha_s(2),e_s).*X(4))) ) ;

F3 = @(alpha_h,alpha_s,e_h,e_s,X) (1-m).*(1-m).*X(3) ...
    + (1-m).*X(1)...
    .*( (1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(1),e_s).*X(3)+e*f_s(alpha_h,alpha_s(2),e_s).*X(4))))...
    .*exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(2),e_s).*X(4)+e*f_s(alpha_h,alpha_s(1),e_s).*X(3))) ...
    +   (1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(1),e_s).*X(3)+e*f_s(alpha_h,alpha_s(2),e_s).*X(4))))...
    .*(1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(2),e_s).*X(4)+e*f_s(alpha_h,alpha_s(1),e_s).*X(3))))...
    .*((1-e)*f_s(alpha_h,alpha_s(1),e_s).*X(3)+e*f_s(alpha_h,alpha_s(2),e_s).*X(4))...
    ./(f_s(alpha_h,alpha_s(1),e_s).*X(3)+f_s(alpha_h,alpha_s(2),e_s).*X(4)));

F4 = @(alpha_h,alpha_s,e_h,e_s,X) (1-m).*(1-m).*X(4) ...
    + (1-m).*X(1)...
    .*( (1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(2),e_s).*X(4)+e*f_s(alpha_h,alpha_s(1),e_s).*X(3))))...
    .*exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(1),e_s).*X(3)+e*f_s(alpha_h,alpha_s(2),e_s).*X(4))) ...
    +   (1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(1),e_s).*X(3)+e*f_s(alpha_h,alpha_s(2),e_s).*X(4))))...
    .*(1-exp(-(1-m)^2.*((1-e)*f_s(alpha_h,alpha_s(2),e_s).*X(4)+e*f_s(alpha_h,alpha_s(1),e_s).*X(3))))...
    .*((1-e)*f_s(alpha_h,alpha_s(2),e_s).*X(4)+e*f_s(alpha_h,alpha_s(1),e_s).*X(3))...
    ./(f_s(alpha_h,alpha_s(1),e_s).*X(3)+f_s(alpha_h,alpha_s(2),e_s).*X(4)));
                              
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
X0ha = 0.2;   % Host alone density 
X0hs = 0.01;  % Host with a symbiont density
choice = 1;   % 1 = more parasitic / 2 = more mutualistic / 3 = coexistence p= p*

if (choice==1)       % more parasitic symbionts
    psp = 0.85;
    X0sp = psp*X0hs;            % Parasitic symbiont density
    X0sm = (1-psp)*X0hs;        % Mutualistic symbiont density
    X0 = [X0ha;X0hs;X0sp;X0sm]; % Initial densities
elseif (choice == 2) % more mutualistic symbionts
    psp  = 0.01;
    X0sp = psp*X0hs;            % Parasitic symbiont density
    X0sm = (1-psp)*X0hs;        % Mutualistic symbiont density
    X0 = [X0ha;X0hs;X0sp;X0sm]; % Initial densities
elseif(choice==3)    % coexistence p =p*
    X0 = [0.1151; 0.5219; 0.0255; 0.4964];
    psp = X0(3)./(X0(3)+X0(4)); %psp_star;
end


%% Computation of the dynamics
Tf = 1000;                  % Final time
t = 0;                     % Initial time
X = zeros(4,Tf+1);
X(:,1) = X0;                    % Initial densities
Xnew = X;
while(t<=Tf)
    Xold = Xnew;
    Xnew = F(alpha_h,alpha_s,e_h,e_s,Xold);
    t = t+1;
    X(:,t) = Xnew;
end

X_h = sum(X(1:2,:));    % Host density
X_sp = X(3,:);          % Parasitic symbiont density
X_sm = X(4,:);          % Mutualistic symbiont density
P_s= X_sm./(X_sp+X_sm); % Proportion of mutualistic symbiont


%% Equilibrium
%%% without symbiont Poisson + space competition
F0_ha =@(alpha_h,X)  -m.*X ...
    + (1-exp(-(1-m).* f_ha(alpha_h,e_h).*X ) )...
    .*(1-(1-m)*(X)).*( 1-((1-m)*(X)).^gammac);
R_ha = fsolve(@(X) F0_ha(alpha_h,X),0.5);
%%%%% Stability
Lambda_ha = (1-m)*(1-m)*(1+R_ha.*f_s(alpha_h,alpha_s,e_s));

A_h = 0:0.01:1;
A_s = 0:0.01:1;
[AH,AS] = meshgrid(A_h,A_s);
R_HA = zeros(length(A_h));
i = 1;
for a_h = A_h
    r_ha = fsolve(@(X) F0_ha(a_h,X),0.5);
    R_HA(:,i) = r_ha;
    i = i+1;
end
LAMBADA_ha = (1-m)*(1-m)*(1+R_HA.*f_s(AH,AS,e_s));



%%% With parasitic/mutualistic symbiont
%%% Parasitic
rho_ha_p = @(X) (1-(1-m)*(1-m))*X./((1-m)*(1-exp(-f_s(alpha_h,alpha_s(1),e_s)*(1-m)*X)));
F01 = @(X) (1-m).*X(1).*exp(-(1-m).*f_s(alpha_h,alpha_s(1),e_s).*X(2)) ...
         + (1-exp(-(1-m).*(f_h(alpha_h,alpha_s(1),e_s).*X(2) + f_ha(alpha_h,e_h).*X(1)) ) )...
        .*(1-(1-m)*(X(1)+X(2))).*( 1-((1-m)*(X(1)+X(2))).^gammac) + m.*(1-m).*X(2);
F0_hs = @(X) F01([rho_ha_p(X),X])-rho_ha_p(X);
Xstar = fsolve(@(X) F0_hs(X),0.5);

R_hs_parasitic_star = Xstar;
R_ha_parasitic_star = rho_ha_p(Xstar);
R_h_parasitic_star  = R_ha_parasitic_star + R_hs_parasitic_star;

%%% Mutualistic
rho_ha_m = @(X) (1-(1-m)*(1-m))*X./((1-m)*(1-exp(-f_s(alpha_h,alpha_s(2),e_s)*(1-m)*X)));
F01 = @(X) (1-m).*X(1).*exp(-(1-m).*f_s(alpha_h,alpha_s(2),e_s).*X(2)) ...
         + (1-exp(-(1-m).*(f_h(alpha_h,alpha_s(2),e_s).*X(2) + f_ha(alpha_h,e_h).*X(1)) ) )...
        .*(1-(1-m)*(X(1)+X(2))).*( 1-((1-m)*(X(1)+X(2))).^gammac) + m.*(1-m).*X(2);
F0_hs = @(X) F01([rho_ha_m(X),X])-rho_ha_m(X);
Xstar = fsolve(@(X) F0_hs(X),0.5);

R_hs_mutualistic_star = Xstar;
R_ha_mutualistic_star = rho_ha_m(Xstar);
R_h_mutualistic_star  = R_ha_mutualistic_star + R_hs_mutualistic_star;


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
    line([0,Tf],[R_h_parasitic_star,R_h_parasitic_star],'color',Color(5,:))
    line([0,Tf],[R_hs_parasitic_star,R_hs_parasitic_star],'color',Color(1,:))
else
    line([0,Tf],[R_h_mutualistic_star,R_h_mutualistic_star],'color',Color(5,:))
    line([0,Tf],[R_hs_mutualistic_star,R_hs_mutualistic_star],'color',Color(2,:))
end

xlabel('Time $t$','Interpreter','latex','FontSize',16)
ylabel('proportion of host and symbionts','Interpreter','latex','FontSize',16)
legend('Hosts','Parasitic symbionts','Mutualistic symbionts','Relative proportion of mutualistic symbiont','Interpreter','latex','FontSize',14)
xlim([0,Tf]) %,0,1])

%%% Stability of the parasitic system
figure(2)
clf
contourf(AH,AS,LAMBADA_ha)
xlabel('host trait $\alpha_h$','Interpreter','latex','FontSize',16)
ylabel('symbiont trait $\alpha_s$','Interpreter','latex','FontSize',16)






