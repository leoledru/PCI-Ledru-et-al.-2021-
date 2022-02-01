%% FIGURE A4

load('data/mutualist_emerg_then_kill_center.mat')

%%
Mutualist_tot(simulation_tot~=3) = NaN;
simulation_tot(Mutualist_tot<0.475) = 2;
simulation_tot(Mutualist_tot>0.475) = 3;

cmap = [0,0,0;
        0,1,0;
        0,0,1;
        1,0,0];
    
subplot(2,3,1)
imagesc(reshape(simulation_tot(:,10000 - 1),100,100));
subplot(2,3,2)
imagesc(reshape(simulation_tot(:,10000),100,100));
subplot(2,3,3)
imagesc(reshape(simulation_tot(:,10000 + 5),100,100));
subplot(2,3,4)
imagesc(reshape(simulation_tot(:,10000 + 10),100,100));
subplot(2,3,5)
imagesc(reshape(simulation_tot(:,10000 + 100),100,100));
subplot(2,3,6)
imagesc(reshape(simulation_tot(:,10000 + 10000),100,100));
colormap(cmap)