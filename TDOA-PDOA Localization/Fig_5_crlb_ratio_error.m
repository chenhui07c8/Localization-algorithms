close all;
clear all;
clc;

%% Initialize system parameters
v = 343;
rng(1); % For reproducibility
radius_all = 2.5;
bound = [-2.5 -2.5; 2.5 2.5];
target_all = [1 1];

lambda = v/21500;
N = 8;
errs = zeros(1,7);

%sigma_phase(sigma_phase_ind); %0.0161 20dB | 0.0125 10dB | 0.0396 0dB | 0.1265 -10dB | 0.4494 -20dB
% f = [20000 21000 22000 23000];
f = 21500;
radius = 2.5;
sensors = [-1 -1 0; -1 0 0;-1 1 0; 0 -1 0; 0 1 0; 1 -1 0; 1 0 0; 1 1 0]*radius;
sensors = sensors(:, 1:2);

range_std = 0.005;
range_std_all = 10.^(-1:-0.25:-4);
phase_std_all = 10.^(-1:-0.25:-4)*0.5; 
error_cell = cell(1,length(phase_std_all));
%% Main code
% tic
% scale_factor = 0.8;
% parfor pdoa_i = 1:length(phase_std_all)
%     swarmsize = 500;
%     phase_std = phase_std_all(pdoa_i); 
%     range_std = range_std_all(pdoa_i);
%     % searching margin
%     margin0 = [range_std*scale_factor range_std*scale_factor];
%     margin0 = max([margin0; 0.01 0.01]);
%     coe = (range_std / phase_std)^2;
%     for sim_i = 1:500
% %         if(mod(sim_i,20)==0)
% %             clc;
% %             disp(['Range Error Index: ' num2str(swarm_i) ' | Sim:' num2str(sim_i)]);
% %         end
% %         for i = 1:2%
%         for i = 1:size(target_all,1)
%             target = target_all(i,:);
%             distance = vecnorm(sensors - target,2,2);
%             distance_noisy = distance + randn(size(distance))*range_std;
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
%             
%             % E1: tdoa-linear
%             [temp, r_tdoa] = alg_tdoa_2d(sensors, -tdoa_obs_noisy', bound);
%             err1 = norm(r_tdoa - target);
%             
%             % E2: tdoa-sdp
%             r0 = alg_sdp_2d(tdoa_all_noisy, sensors', v);
%             err2 = norm(r0' - target);
% 
%             % E3: tdoa-pso
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
%             fun = @(x)coe*sum(sum((wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f) + ...
%                 sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);  
%             fix_margin = margin0;
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
%             fun = @(x)sum(sum((wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f);  
%             fix_margin = margin0;
%             lb = max([-fix_margin + r_tdoa; bound(1,:)]);
%             ub = min([fix_margin + r_tdoa; bound(2,:)]);
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',swarmsize,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,~] = particleswarm(fun,nvars,lb,ub,options);
%             r_tp = x; 
%             err5 = norm(r_tp - target);
%            
%             error_cell{pdoa_i}(sim_i,1:5) = [err1 err2 err3 err4 err5];
%         end
%     end
% 
% end
% % errors
% 
% toc

%% last run result

load crlb_fig4.mat

% print -dpng -r600 sim-3.png

%% rmse [mm]
p = zeros(length(error_cell),5);
for i = 1:length(p)
    temp = error_cell{i};
    if(isempty(temp))
        p(i,:) = 0;
    else
        p(i,:) = sqrt(mean(temp.^2))*1000;
    end
%     p(i) = sqrt(mean(error_cell{i}.^2));
end
% figure;plot(p(1:20,:))
% sigma = 20:20:400;
ind = 1:(length(error_cell)-2);
% ind = 1:length(error_cell);
sigma = range_std_all*1000;
figure; loglog(sigma(ind), (p(ind, 1)),'-c+', 'LineWidth',1.2);
hold on; loglog(sigma(ind), (p(ind,2)),'-g*', 'LineWidth',1.2);
hold on; loglog(sigma(ind), (p(ind,3)),'-ms', 'LineWidth',1.2);
hold on; loglog(sigma(ind), (p(ind,5)),'-r^', 'LineWidth',1.2);
hold on; loglog(sigma(ind), (p(ind,4)),'-bd', 'LineWidth',1.2);
% hold on; semilogy(phase_std_all_mm, (p(:,7)),'-r^', 'LineWidth',1.2);
% hold on; semilogy(phase_std_all_mm, (p(:,6)),'-bo', 'LineWidth',1.2);

xlabel('TDOA standard deviation [mm]')
ylabel('RMSE [mm]')
% legend({'TDOA Linear', 'TDOA SDP', 'TDOA PSO', 'TDOA-PDOA-1 [s=100]', 'TDOA-PDOA-2 [s=100]'},'Location','southwest')
% legend({'TDOA Linear', 'TDOA SDP', 'TDOA-PDOA-1 [s=100]', 'TDOA-PDOA-2 [s=100]', 'TDOA-PDOA-1 [s=400]', 'TDOA-PDOA-2 [s=400]'},'Location','southeast')

grid on;
set(gca,'FontSize',12)
% axis([0 10 10^0 121])
axis([min(sigma(ind)) max(sigma(ind)) 10^-1.6 10^2.5])

% print -dpng -r600 sim-3.png

%%
L = 8;
v = 343;
Evec = ones(L,1);
x = target_all(1);
y = target_all(2);


pairs0 = nchoosek([1:L],2);
pairs = pairs0(1:(L-1),:);  % non-redundant set
g = zeros(2, length(pairs));
for i = 1:length(pairs)
%     g(1,i) = (x-sensors(i,1))/sqrt((x-sensors(i,1))^2 + (y-sensors(i,2))^2);
%     g(2,i) = (y-sensors(i,2))/sqrt((x-sensors(i,1))^2 + (y-sensors(i,2))^2);
    a = pairs(i,1);
    b = pairs(i,2);
    g1 = ([x y]-sensors(a,:))/norm([x y] - sensors(a,:));
    g2 = ([x y]-sensors(b,:))/norm([x y] - sensors(b,:));
    g(:,i) = g2' - g1';
end
G = g;

crlb_t = zeros(1,length(phase_std_all));
crlb_p = zeros(1,length(phase_std_all));
crlb_h = zeros(1,length(phase_std_all));

for i = 1:length(phase_std_all)
    % tdoa
    range_std = range_std_all(i);
    varD = range_std^2; 
    Ct = varD*(ones(length(G))+eye(length(G)));
    GCG = G*Ct^-1*G';
    F = GCG;
    temp = inv(F);
    crlb_t(i) = trace(temp);
    
    % pdoa
    phase_std = phase_std_all(i);
    varD = phase_std^2;
    Cp = varD*(ones(length(G))+eye(length(G)));
    GCG = G*Cp^-1*G';
    F = GCG;
    temp = inv(F);
    crlb_p(i) = trace(temp);

    % hybrid
    GCG = [G G]*[Ct zeros(size(Ct)); zeros(size(Cp)) Cp]^-1*[G G]';
    F = GCG;
    crlb_h(i) = trace(inv(F));

end

% end
% range_std
% sigma = phase_std_all;
hold on;semilogy(sigma, sqrt(crlb_t*1000^2),'--+g','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt(crlb_p*1000^2),'--xr','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt(crlb_h*1000^2),'-..b','MarkerSize',10,'LineWidth',1.2);

legend({'TDOA Linear', 'TDOA SDP', 'TDOA Iter', 'TDOA-PDOA-1', 'TDOA-PDOA-2', 'CRLB-T'...
   'LB-P','LB-J' },'Location','southeast')

% print -dpng -r600 crlb-4.png
