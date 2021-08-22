%% FIGURE 1 
close
alpha_h = 0;
alpha_s = 0:0.01:1;

Fmut_min_s = 2.5;
Fmut_min_h = 0.1;
Fmut_max = 8;
Fb_min_h = .7;
Fb_min_s = .7;
Fb_max = 1;
gamma_1_h = 4;
gamma_1_s = 1;
gamma_2 = 1;
gamma_h = 1;
f_mutu_s = @(alpha_h) Fmut_min_s+(Fmut_max-Fmut_min_s).*alpha_h.^gamma_1_s;
f_mutu_h = @(alpha_s) Fmut_min_h+(Fmut_max-Fmut_min_h).*alpha_s.^gamma_1_h;
f_basal_s = @(alpha_s) Fb_max -(Fb_max-Fb_min_s).*alpha_s.^gamma_2;
f_basal_h = @(alpha_h) Fb_max-(Fb_max-Fb_min_h).*alpha_h.^gamma_2; 
fhseul_max = .5;
fhseul_min = fhseul_max*Fb_min_h;
fitness_hseul=@(alpha_h)fhseul_max-(fhseul_max-fhseul_min).*alpha_h.^gamma_h; 
f_seul = fitness_hseul(alpha_h); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_host_associated = f_basal_h(alpha_h).*f_mutu_h(alpha_s);

subplot(1,2,1)
plot(alpha_s,F_host_associated,'--k','LineWidth',3)
hold on
plot([0,1],[f_seul,f_seul],'-k','LineWidth',3)
plot([0.475,0.475],[0,8],'-.r','LineWidth',2)
text(.438,3,'parasitism','color','r','interpreter','latex','FontSize',20,'rotation',90);
text(.51,3,'mutualism','color','r','interpreter','latex','FontSize',20,'rotation',90);
legend({['host associated' 10 'with symbiont'],'host alone'},'interpreter','latex','FontSize',15,'Location','northwest')
xlabel('symbiont interaction trait $\alpha_s$','interpreter','latex','FontSize',25)
ylabel('host fecundity','interpreter','latex','FontSize',25)
hold off
%
subplot(1,2,2)
alpha_h = 0:0.01:1;
alpha_s = 0;
F_symbiont = f_basal_s(alpha_s).*f_mutu_s(alpha_h);
plot(alpha_h,F_symbiont,'--k','LineWidth',3)
xlabel('host interaction trait $\alpha_h$','interpreter','latex','FontSize',25)
ylabel('symbiont fecundity','interpreter','latex','FontSize',25)
