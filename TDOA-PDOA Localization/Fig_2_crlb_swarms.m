% test different swarm sizes
close all;
clear all;
clc;

%% Initialize system parameters
% load tdoa_std.mat;
% load sigma_phase;
rng(5); % For reproducibility
v = 343;
radius_all = 2.5;
bound = [-2.5 -2.5; 2.5 2.5];
target_all = [1 1];

N = 8;
% swarm_range = round(10.^(0.2:0.2:3.6));
swarm_range = round(10.^(0.2:0.2:3.0));
error_cell = cell(1,length(swarm_range));

errs = zeros(1,7);
% f = [20000 21000 22000 23000];
f = 21500;
radius = 2.5;
sensors = [-1 -1 0; -1 0 0;-1 1 0; 0 -1 0; 0 1 0; 1 -1 0; 1 0 0; 1 1 0]*radius;
sensors = sensors(:, 1:2);

% range_std = 0.05;
range_std = 0.002;
range_var = range_std^2;
phase_std_rad = 0.1265; %0.0161 20dB | 0.0125 10dB | 0.0396 0dB | 0.1265 -10dB | 0.4494 -20dB
phase_std = (phase_std_rad/2/pi*v/21500);

%% Main code
% tic
% scale_factor = 0.8;
% err1 = 0; % linear closed form
% err2 = 0; % SDP
% err3 = 0; % tdoa only pso, s = 10
% err4 = 0; % joint, M=4
% err5 = 0; % pdoa only, M = 4
% err6 = 0;
% err7 = 0;
% % parfor swarm_i = 1:length(swarm_range)
% parfor swarm_i = 1:length(swarm_range)
% 
%     swarmsize = swarm_range(swarm_i);
% %     margin = [range_std*scale_factor range_std*scale_factor];
%     margin = [0.02 0.02];
%     coe = (range_std / phase_std)^2;
%     for sim_i = 1:2000
%         if(mod(sim_i,10)==0)
%             disp(['sim_i:' num2str(sim_i)]);
%         end
%         for i = 1:size(target_all,1)
%             target = target_all(i,:);
%             distance = vecnorm(sensors - target,2,2);
%             distance_noisy = distance + randn(size(distance))*sqrt(range_var);
%             comb = combnk(1:N,2);
%             tdoa_all = distance(comb(:,1)) - distance(comb(:,2));
%             tdoa_obs = [0; tdoa_all(1:7)];
%             tdoa_all_noisy = distance_noisy(comb(:,1)) - distance_noisy(comb(:,2));
%             tdoa_obs_noisy = [0; tdoa_all_noisy(1:7)];
%             
%             distance_phase_noisy = distance + randn(size(distance))*phase_std;
%             tdoa_phase_noisy = distance_phase_noisy(comb(:,1)) - distance_phase_noisy(comb(:,2));
%             phase_all_noisy = wrapping(tdoa_phase_noisy./(v./f)*2*pi);
% 
% %             E1: tdoa-linear
%             [~, r_tdoa] = alg_tdoa_2d(sensors, -tdoa_obs_noisy', bound);
%             err1 = norm(r_tdoa - target);
%             
% %             E2: tdoa-sdp
%             r0 = alg_sdp_2d(tdoa_all_noisy, sensors', v);
%             err2 = norm(r0' - target);
% 
% %             E3: tdoa-pso
% %             clear x
%             fun = @(x)sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);  
%             lb = [-radius_all -radius_all];
%             ub = [radius_all radius_all];
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',20,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,~] = particleswarm(fun,nvars,lb,ub,options);
%             err3 = norm(x - target);
%             r_tdoa = x;
%             
%             % E4: phase + tdoa, multiple frequency, small swarmsize
% %             clear x
%             fun = @(x)coe*sum(sum((wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f) + ...
%                 sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);  
%             fix_margin = margin;
%             lb = max([-fix_margin + r_tdoa; bound(1,:)]);
%             ub = min([fix_margin + r_tdoa; bound(2,:)]);
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',swarmsize,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,~] = particleswarm(fun,nvars,lb,ub,options);
%             r_tp = x; 
%             err4 = norm(r_tp - target);
%             
%             % E5: pdoa only
% %             clear x
%             fun = @(x)sum(sum((wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f);  
%             fix_margin = margin;
%             lb = max([-fix_margin + r_tdoa; bound(1,:)]);
%             ub = min([fix_margin + r_tdoa; bound(2,:)]);
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',swarmsize,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,~] = particleswarm(fun,nvars,lb,ub,options);
%             r_tp = x; 
%             err5 = norm(r_tp - target);
% 
%             error_cell{swarm_i}(sim_i,1:7) = [err1 err2 err3 err4 err5 err6 err7];
%         end
%     end
%     
% end
% 
% toc
%% Last run results

% save crlb_fig1.mat
load crlb_fig1.mat

%% rmse [mm]
errors_rmse = zeros(length(swarm_range), length(errs));
p = zeros(length(swarm_range),length(errs));
for i = 1:size(p,1)
    temp = error_cell{i};
    p(i,:) = sqrt(mean(temp.^2))*1000;
%     p(i) = sqrt(mean(error_cell{i}.^2));
end
% figure;plot(p(1:20,:))

sigma = swarm_range;
figure; loglog(sigma, (p(:, 1)),'-c+', 'LineWidth',1.2);
hold on; loglog(sigma, (p(:,2)),'-g*', 'LineWidth',1.2);
hold on; loglog(sigma, p(:,3),'-ms', 'LineWidth',1.2);
hold on; loglog(sigma, (p(:,5)),'-r^', 'LineWidth',1.2);
hold on;loglog(sigma, (p(:,4)),'-bo', 'LineWidth',1.2);

xlabel('Swarm size')
ylabel('RMSE [mm]')
% legend({'TDOA Linear', 'TDOA SDP', 'TDOA Iter', 'TDOA-PDOA-1', 'TDOA-PDOA-2'},'Location','northeast')
grid on;
set(gca,'FontSize',12)
xlim([2 max(sigma)*1.5]);

% 0.1265/2/pi*v/21500*1000 = 0.32 mm

% print -dpng -r600 sim-2.png
%% CRLB
L = 8;
v = 343;
Evec = ones(L,1);
x = target_all(1);
y = target_all(2);

pairs0 = nchoosek([1:L],2);
pairs = pairs0(1:(L-1),:);  % non-redundant set
g = zeros(2, length(pairs));
for i = 1:length(pairs)
    a = pairs(i,1);
    b = pairs(i,2);
    g1 = ([x y]-sensors(a,:))/norm([x y] - sensors(a,:));
    g2 = ([x y]-sensors(b,:))/norm([x y] - sensors(b,:));
    g(:,i) = g2' - g1';
end
G = g;

varD = range_std^2; 
C = varD*(ones(length(G))+eye(length(G)));
GCG = G*C^-1*G';
F = GCG;
temp = inv(F);
crlb_t = trace(temp);

varD = phase_std^2;
C = varD*(ones(length(G))+eye(length(G)));
GCG = G*C^-1*G';
F = GCG;
temp = inv(F);
crlb_p = trace(temp);

coe = (range_std / phase_std)^2;
crlb_h = (coe^2*crlb_p + crlb_t)/(1+coe)^2;
% crlb_h = (crlb_p * crlb_t)/(crlb_p + crlb_t);

% end
hold on;semilogy(sigma, sqrt((crlb_t*ones(1,length(sigma))*1000^2)),'--+g','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt((crlb_p*ones(1,length(sigma))*1000^2)),'--xr','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt((crlb_h*ones(1,length(sigma))*1000^2)),'-..b','MarkerSize',10,'LineWidth',1.2);

legend({'TDOA Linear', 'TDOA SDP', 'TDOA Iter', 'TDOA-PDOA-1', 'TDOA-PDOA-2', 'CRLB-T', 'LB-P','LB-J'},'Location','northeast')

% figure;loglog(range_var, mse_all(1,:),'-o');
% hold on;loglog(range_var, mse_all(2,:),'-*');
% hold on; loglog(range_var, crlb_all);
% 
% legend('MSE-LS','MSE-SDP','CRLB')
% xlabel('Variance of Ranging Error [m]')
% ylabel('MSE of Location Estimations')

% print -dpng -r600 crlb-1.png

