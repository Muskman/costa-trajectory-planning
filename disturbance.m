function [g, grad_f] = disturbance(x,e,dt,opts)
    % function to get ocean velocity / disturbance values at given set of waypoints x
    T = opts.T;
    N = opts.n_agents;

    p = opts.goal_weight;
    m = opts.m_scale;
    
    g = zeros(size(x));
    grad_f = zeros(size(x));

    for j = 1:N
        X = reshape(x(2*(j-1)*T+1:2*j*T),[2,T]); 
        E = reshape(e(2*(j-1)*T+1:2*j*T),[2,T]);
        L = zeros(T);

        if strcmp(opts.env,'test')
            g_temp = zeros(2,T);
            grad_temp = cell(T,1); 
        
            for i = 1:T
                a = X(1,i); b = X(2,i); 
                g_temp(:,i) = (opts.c_scale/m)*[(1-2*(a/m)^2); -2*(a/m)*(b/m)].*((exp(-(a^2+b^2)/m^2)*(1+E(:,i))));
                grad_temp{i} = (opts.c_scale/m^3)*[2*a*(2*(a/m)^2-3) 2*b*(2*(a/m)^2-1); 2*b*(2*(a/m)^2-1) 2*a*(2*(b/m)^2-1)].*(exp(-(a^2+b^2)/m^2)*(1+E(:,i)));
                L(i,i) = p*i/T;
            end
        else
            [g_temp, grad_temp] = find_ocean_vel(X,opts);

            for i = 1:T
                L(i,i) = p*i/T;
            end
        end

        g_temp = reshape(g_temp,[2*T,1]);
        g(2*(j-1)*T+1:2*j*T) = g_temp;

        grad_temp{T} = zeros(2);
        grad_temp = blkdiag(grad_temp{:});
        
        A = opts.A; B = opts.B; D = opts.D;
        
        xG = repmat(opts.x_goal(2*j-1:2*j),T,1);
        L = kron(L,eye(2));
        
        grad_f(2*(j-1)*T+1:2*j*T) = 2*L*(x(2*(j-1)*T+1:2*j*T)-xG) + A*x(2*(j-1)*T+1:2*j*T) + dt*B*g(2*(j-1)*T+1:2*j*T) + 2*(eye(2*T)-L)*dt*grad_temp'*(D*x(2*(j-1)*T+1:2*j*T)+dt*g(2*(j-1)*T+1:2*j*T));
    end

end