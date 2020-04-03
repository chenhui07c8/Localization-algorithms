function [temp,output] = alg_tdoa_2d(sensors, tdoa_obs, bound)
%     tdoa_obs = [0; tdoa_obs];
    v = 343;
    t = tdoa_obs/v;

    x = sensors(:,1);
    y = sensors(:,2);
%     z = sensors(:,3);
    aa = [];
    bb = [];
    for m = 3:length(t)
        A = (-2*x(1)+2*x(m))/v/t(m) - (2*x(2)-2*x(1))/v/t(2);
        B = (-2*y(1)+2*y(m))/v/t(m) - (2*y(2)-2*y(1))/v/t(2);
%         C = 0;
        D = v*t(m) - v*t(2) + (x(1)^2+y(1)^2-x(m)^2-y(m)^2)/v/t(m) - (x(1)^2+y(1)^2-x(2)^2-y(2)^2)/v/t(2);
        aa = [aa; A B];
        bb = [bb; -D];
    end


    a2 = aa(:,1:2);
    b2 = bb;
    temp = (a2'*a2)^-1*a2'*b2;
    temp2 = sort([temp'; bound(:,1:2)]);
    x0 = temp2(2,1);
    y0 = temp2(2,2);
    
%     err_all = zeros(1,length(increment));
%     for i = 1:length(increment)
%         z0 = increment(i);
%         dis = vecnorm(sensors - [x0 y0 z0], 2, 2);
%         tdoa_est = dis - dis(1);
%         err_all(i) = norm(tdoa_est-t*v);
%     end
%     % figure;plot(err_all)
%     [~, q] = min(err_all);
%     z0 = increment(q);
    % err = norm([x0 y0 z0] - target)

    output = [x0 y0];
end