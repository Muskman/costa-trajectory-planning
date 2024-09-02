function y = minObj(x,N,K,T,std_noise,opts)
    dt = opts.dt;
    y = 0; 
    for k = 1:K
        e = std_noise*randn(2*N*T); % generate realization of noise
        [g, ~] = disturbance(x,e,dt,opts);
        for j = 1:N
            y = y + norm(x(2*(j-1)*T+3:2*j*T)-x(2*(j-1)*T+1:2*j*T-2)-dt*g(2*(j-1)*T+1:2*j*T-2))^2;
        end
    end
end