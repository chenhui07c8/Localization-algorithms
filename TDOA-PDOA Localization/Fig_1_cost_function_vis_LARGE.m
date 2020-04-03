close all;
clear all;
clc;
%% 
target0 = [0 0];
radius = 2.5;
f = 20000:1000:23000;
v = 343;
sensors = [-1 -1 0;-1 0 0;-1 1 0; 0 -1 0; 0 1 0; 1 -1 0;1 0 0; 1 1 0]*radius;
N = length(sensors);
sensors = sensors(1:8,1:2);
N = length(sensors);
target_can = [];
area_range = -5:0.05:5;
for xa = area_range + target0(1)
    for ya = area_range + target0(2)
        target_can = [target_can; xa ya];
    end
end
errors = [];
cost1 = zeros(1,size(target_can,1));
cost2 = zeros(1,size(target_can,1));
   
range_std = 0.0;
phase_std = 0;

tic
for i = 1:size(target_can,1)
    if(mod(i,1000) == 0)
        disp(i);
    end
    target = target_can(i,:);
    distance = vecnorm(sensors - target0,2,2);
    comb = nchoosek(1:N,2);


    tdoa_all = distance(comb(:,1)) - distance(comb(:,2));   % tdoa of all the combinations
    tdoa_obs = [0; tdoa_all(1:(N-1))];      % take the first one as the reference receiver
    tdoa_all_noisy = tdoa_all + randn(size(tdoa_all))*range_std; 
    tdoa_obs_noisy = [0; tdoa_all_noisy(1:(N-1))];

    phase_all = wrapping(tdoa_all./(v./f)*2*pi);
    phase_all_noisy = phase_all + randn(size(phase_all))*phase_std;
    x = target;
    cost1(i) = sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);
%     ptemp = abs(wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy))*1/2/pi*v/f;
    ptemp = (wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2;
    cost2(i) = sum(sum(ptemp))/length(f);
end
toc

coe = 500;

%%
% cost 1, tdoa
figure;
imagesc(reshape(cost1,[length(area_range) length(area_range)]));
xticks([0 50 100 150 200]);
xticklabels({'-1','-0.5','0','0.5','1'});
yticks([0 50 100 150 200]);
yticklabels({'1','0.5','0','-0.5','-1'});
colorbar;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
set(gca,'FontSize',15)
set(gcf, 'Position',  [100, 100, 450, 350])

% cost 2
figure;
imagesc(reshape(cost2,[length(area_range) length(area_range)]));
xticks([0 50 100 150 200]);
xticklabels({'-1','-0.5','0','0.5','1'});
yticks([0 50 100 150 200]);
yticklabels({'1','0.5','0','-0.5','-1'});
% caxis([0 3]);
colorbar;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
set(gca,'FontSize',15)
set(gcf, 'Position',  [100, 100, 450, 350])

% cost 3
figure;
imagesc(reshape(cost1 + coe*cost2,[length(area_range) length(area_range)]));
xticks([0 50 100 150 200]);
xticklabels({'-1','-0.5','0','0.5','1'});
yticks([0 50 100 150 200]);
yticklabels({'1','0.5','0','-0.5','-1'});
% caxis([0 3]);
colorbar;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
set(gca,'FontSize',15)
set(gcf, 'Position',  [100, 100, 450, 350])

%%
% cost tdoa mesh
Z = reshape(cost1,[length(area_range) length(area_range)]);
[X,Y] = meshgrid(area_range);
figure;
mesh(X,Y,Z);
% set(gca, 'Zdir', 'reverse')
xlabel('x-axis [m]');
ylabel('y-axis [m]');
zlabel('c_t [m^2]');
view([55 20]);set(gca,'FontSize',15)


% cost pdoa mesh
Z = reshape(cost2,[length(area_range) length(area_range)]);
[X,Y] = meshgrid(area_range);
figure;
mesh(X,Y,Z);
% set(gca, 'Zdir', 'reverse')
xlabel('x-axis [m]');
ylabel('y-axis [m]');
zlabel('c_p [m^2]');
view([55 20]);set(gca,'FontSize',15)


% cost joint mesh
Z = reshape(cost1+coe*cost2,[length(area_range) length(area_range)]);
[X,Y] = meshgrid(area_range);
figure;
mesh(X,Y,Z);
% set(gca, 'Zdir', 'reverse')
xlabel('x-axis [m]');
ylabel('y-axis [m]');
zlabel('c_j [m^2]');
% view([50 45]);set(gca,'FontSize',15)
view([55 20]);set(gca,'FontSize',15)

%  print -painters -dpng -r150 fig-1-3.png

% % Z = cost2;
% regmax = imregionalmin(Z);
% disp('Local minima//////////////')
% sum(Z(regmax)<0.4*10^-3)
% sum(regmax(:))
% SS = sort(Z(regmax), 'descend');
% diff(-SS(end-1:end))*1000
