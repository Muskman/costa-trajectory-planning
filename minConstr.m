function [y, yeq] = minConstr(x,N,T,x_obs,r,opts)
    dt = opts.dt;
    y = -vecnorm(reshape(x-repmat(x_obs,N*T,1),2,N*T),2,1) + r;
    
    y_vel = zeros(1,N*(T-1)); yeq = zeros(1,4*N);
    e0 = zeros(2*N*T); % generate realization of noise
    [g, ~] = disturbance(x,e0,dt,opts);
    for j = 1:N
        y_vel((j-1)*(T-1)+1:j*(T-1)) = vecnorm(reshape(x(2*(j-1)*T+3:2*j*T)-x(2*(j-1)*T+1:2*j*T-2)-dt*g(2*(j-1)*T+1:2*j*T-2),2,T-1),2,1)-opts.v;
        yeq(4*j-3:4*j) = [x(2*(j-1)*T+1:2*(j-1)*T+2) - opts.x_start(2*j-1:2*j); x(2*j*T-1:2*j*T) - opts.x_goal(2*j-1:2*j)]';
    end

    % colission avoidance for two robots
    y_inter = -vecnorm(reshape(diff(reshape(x,2*T,N),1,2),2,T),2,1) + opts.r_a;
    y = [y y_vel y_inter];
end
