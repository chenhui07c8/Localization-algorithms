close all;
clear all;
clc;

%% Initialize system parameters
load tdoa_std.mat;
load sigma_phase;
sigma_phase_ind = 1; 

v = 343;
target_ind = 19;
target_all(19,:)

% range_errors = 0.5e-3:10e-3:60e-3;
% log10(10e-3)
rng(2); % For reproducibility
radius_all = 2.5;
bound = [-2.5 -2.5; 2.5 2.5];
target_all = [1 1];

N = 8;
swarmsize = 500;
tdoa_type = 1;
% sigma_range_ind = 1:5:length(sigma_range);
% sigma_range_ind = 0:5:100;
% sigma_range_ind(1) = 1;
errs1 = 0; % linear closed form
errs2 = 0; % SDP
errs3 = 0; % tdoa only, s = 10
errs4 = 0; % joint M=4, s = 100
errs5 = 0; % pdoa only M=4, s = 500

% % err7 = 0; pdoa only, M = 4, s = 500
% % err8 = 0; pdoa only, M = 4, s = 100
radius = 2.5;
sensors = [-1 -1 0; -1 0 0;-1 1 0; 0 -1 0; 0 1 0; 1 -1 0; 1 0 0; 1 1 0]*radius;
sensors = sensors(:, 1:2);

range_std_all = 10.^(-1:-0.5:-5);
phase_std_rad = 0.1265; %sigma_phase(sigma_phase_ind); %0.0161 20dB | 0.0125 10dB | 0.0396 0dB | 0.1265 -10dB | 0.4494 -20dB

p = zeros(5, length(range_std_all));
error_cell = cell(1,length(range_std_all));
phase_std = (phase_std_rad/2/pi*v/21500);

%% Main code
% f = 21500;
% scale_factor = 0.8;
% 
% tic
% parfor range_error_i = 1:length(range_std_all)
%     range_std = range_std_all(range_error_i);
%     % margin 1: based on TDOA
% %     margin = std_cell{sigma_range_ind(range_error_i)}(target_ind,1:2)*scale_factor;   
%     % margin 2: based on range_std
%     margin0 = [range_std*scale_factor range_std*scale_factor];
%     margin0 = max([margin0; 0.01 0.01]);
% %     margin
%     coe = (range_std / (phase_std_rad/2/pi*v/21500))^2;
%     % start...
%     simtimes = 1000;
%     errors = zeros(simtimes,5);
%     for sim_i = 1:simtimes
%         for i = 1:size(target_all,1)
%             if(mod(i,5)==0)
%                 disp(i);
%             end
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
%             % E1: tdoa-linear
%             [temp, r_tdoa] = alg_tdoa_2d(sensors, -tdoa_obs_noisy', bound);
%             errs1 = norm(r_tdoa - target);
%             
%             % E2: tdoa-sdp
%             r0 = alg_sdp_2d(tdoa_all_noisy, sensors', v);
%             errs2 = norm(r0' - target);
% 
%             % E3: tdoa-pso
%             fun = @(x)sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);  
%             lb = [-radius_all -radius_all];
%             ub = [radius_all radius_all];
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',20,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,fval] = particleswarm(fun,nvars,lb,ub,options);
%             errs3 = norm(x - target);
%             r_tdoa = x;
%             
%             % E4: pdoa only
%             fun = @(x)sum(sum((wrapping( (vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2))./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f);  
%             fix_margin = margin0;
%             lb = max([-fix_margin + r_tdoa; bound(1,:)]);
%             ub = min([fix_margin + r_tdoa; bound(2,:)]);
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',swarmsize,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,fval] = particleswarm(fun,nvars,lb,ub,options);
%             r_tp = x; 
%             errs4 = norm(r_tp - target);
%             
%             
%             % E5: phase + tdoa, multiple frequency, small swarmsize
%             fun = @(x)coe*sum(sum((wrapping(( vecnorm(sensors(comb(:,1),:)-x,2,2)-vecnorm(sensors(comb(:,2),:)-x,2,2) )./(v./f)*2*pi - phase_all_noisy)*1/2/pi*v./f).^2))/length(f) + ...
%                 sum((vecnorm(sensors(comb(:,1),:) - x, 2, 2) - vecnorm(sensors(comb(:,2),:)-x,2,2) - tdoa_all_noisy).^2);  
%             fix_margin = margin0;
%             lb = max([-fix_margin + r_tdoa; bound(1,:)]);
%             ub = min([fix_margin + r_tdoa; bound(2,:)]);
%             nvars = 2;
%             options = optimoptions('particleswarm','SwarmSize',swarmsize,'HybridFcn',@fmincon);
%             options.Display = 'off';
%             [x,fval] = particleswarm(fun,nvars,lb,ub,options);
%             r_tp = x; 
%             errs5 = norm(r_tp - target);
%             error_cell{range_error_i}(sim_i,1:5) = [errs1 errs2 errs3 errs4 errs5];
%         end
%     end
% end
% 
% toc
%% last run result

% save crlb_fig2.mat
load crlb_fig2.mat

%% rmse log[mm]
errors_rmse = zeros(length(range_std_all), 5);
p = zeros(length(error_cell),5);

for i = 1:length(p)
    temp = error_cell{i};
    p(i,:) = sqrt(mean(temp.^2))*1000;
end

% sigma = sigma_range*1000;
sigma = range_std_all*1000;
figure; loglog(sigma, p(:, 1),'-c+', 'LineWidth',1.2);
hold on; loglog(sigma, p(:,2),'-g*', 'LineWidth',1.2);
hold on; semilogy(sigma, p(:,3),'-ms', 'LineWidth',1.2);
hold on; loglog(sigma, p(:,4),'-r^', 'LineWidth',1.2);
hold on; loglog(sigma, p(:,5),'-bo', 'LineWidth',1.2);

legend({'TDOA Linear', 'TDOA SDP', 'TDOA PSO', 'TDOA-PDOA-1', 'TDOA-PDOA-2'},'Location','southeast')
xlabel('TDOA standard deviation [mm]')
ylabel('RMSE [mm]')
grid on;
set(gca,'FontSize',12)

% axis([0 100 10^-1.5 300])
% print -dpng -r600 sim-1.png


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

crlb_t = zeros(1,length(range_std_all));
crlb_p = zeros(1,length(range_std_all));
crlb_h = zeros(1,length(range_std_all));

for i = 1:length(range_std_all)
    % tdoa
    range_std = range_std_all(i);
    varD = range_std^2; 
    Ct = varD*(ones(length(G))+eye(length(G)));
    GCG = G*Ct^-1*G';
    F = GCG;
    temp = inv(F);
    crlb_t(i) = trace(temp);
    
    % pdoa
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
hold on;semilogy(sigma, sqrt(crlb_t*1000^2),'--+g','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt(crlb_p*1000^2),'--xr','MarkerSize',6,'LineWidth',1.2);
hold on;semilogy(sigma, sqrt(crlb_h*1000^2),'-..b','MarkerSize',10,'LineWidth',1.2);

legend({'TDOA Linear', 'TDOA SDP', 'TDOA Iter', 'TDOA-PDOA-1', 'TDOA-PDOA-2', 'CRLB-T'...
   'LB-P','LB-J' },'Location','northwest')

% print -dpng -r600 crlb-2.png
