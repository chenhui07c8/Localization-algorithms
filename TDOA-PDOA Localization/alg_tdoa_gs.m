function [temp,output] = alg_tdoa_gs(sensors, tdoa_obs, bound)

    x = sensors(:,1);
    y = sensors(:,2);
    A = [];
    w = [];
    for m = 2:length(tdoa_obs)
        dm0 = tdoa_obs(m);
        wm0 = 1/2*(dm0^2 - x(m)^2 + x(1)^2-y(m)^2+y(1)^2);
        a = [x(1)-x(m) y(1)-y(m) dm0];
        A = [A; a];
        w = [w; wm0];
    end
    temp = (A'*A)^-1*A'*w;
   
	temp2 = sort([temp(1:2)'; bound(:,1:2)]);
    x0 = temp2(2,1);
    y0 = temp2(2,2);
    output = [x0 y0];

end