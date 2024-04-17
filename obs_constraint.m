function [C, d] = obs_constraint(x,x_obs,r,opts)
    % function to generate constraint over T poses
    
    T = opts.T;
    N = opts.n_agents;
    n_obs = opts.n_obs;
    
    C = cell(N,1); d = zeros(N*T*n_obs,1);
    % d = [];

    for j = 1:N
        for i = 1:n_obs
            x_obs = opts.x_obs(:,i);
            X = reshape(x(2*(j-1)*T+1:2*j*T),[2,T]);
        
            C_d = repmat(x_obs,1,T)-X;  % 2xT
            C_temp2 = num2cell(C_d,1);
            C{j} = [C{j}; blkdiag(C_temp2{:})'];    % [size(C{j}); Tx2T]
            
            d(n_obs*(j-1)*T+(i-1)*T+1:n_obs*(j-1)*T+i*T) = (x_obs'*x_obs - (vecnorm(X,2,1).^2)' - r(i)^2)/2;
        end
    end
    C = blkdiag(C{:});

    if N>1
        C_inter = zeros(T,length(x)); d_inter = zeros(T,1);
        for i = 1:T
            At = zeros(2,length(x)); At(:,2*i-1:2*i) = eye(2); At(:,2*i-1+T:2*i+T) = eye(2);
            Pt = At'*At;
            C_inter(i,:) = -2*x'*Pt; 
            d_inter(i) = -(x'*Pt*x + 4*opts.r_a^2);
        end
        try
            C = [C;C_inter]; d = [d; d_inter];
        catch
            % for debugging
            keyboard
        end
    end
end