function x = qcqp(b,C,d,R,g,s,x0,opts)
%%%    min  xTx -2bTx
%%%    s.t. Cx <= d
%%%         xTRtTRtx <= s   for all [T]
    T = opts.T;
    N = opts.n_agents;
    n_obs = opts.n_obs;
    f = -2*b;
    
    % min   0.5xTQx + fTx + c 
    % s.t.  0.5xTHix + kiTx + ddi <= 0 
    %       0.5xTJix + piTx + qi = 0
    Q = 2*eye(2*N*T);

    if N==1
        H = cell(N*(n_obs+1)*T,1);
        k = cell(N*(n_obs+1)*T,1);
        dd = cell(N*(n_obs+1)*T,1);
    else
        H = cell(N*(n_obs+1)*T+T,1);
        k = cell(N*(n_obs+1)*T+T,1);
        dd = cell(N*(n_obs+1)*T+T,1);
    end

    for j = 1:N
        for l = 1:n_obs
            for i = 1:T
                idx = (n_obs+1)*(j-1)*T+(l-1)*T+i;
                H{idx} = zeros(2*N*T,2*N*T);
                k{idx} = C(n_obs*(j-1)*T+(l-1)*T+i,:)';
                dd{idx} = - d(n_obs*(j-1)*T+(l-1)*T+i);
                if l == n_obs
                    H_temp = zeros(2*N*T,2*N*T);
                    H_temp(2*(j-1)*T+1:2*j*T,2*(j-1)*T+1:2*j*T) = 1*R(2*i-1:2*i,:)'*R(2*i-1:2*i,:);
                    H{idx+T} = H_temp;

                    k_temp = zeros(2*N*T,1);
                    k_temp(2*(j-1)*T+1:2*j*T,1) = -2*R(2*i-1:2*i,:)'*g(2*(j-1)*T+i:2*(j-1)*T+i+1)*opts.dt;
                    k{idx+T} = k_temp;

                    dd{idx+T} = norm(opts.dt*g(2*(j-1)*T+i:2*(j-1)*T+i+1))^2 - s;
                end
            end
        end
    end

    if N>1
        for i = 1:T
            H{N*(n_obs+1)*T+i} = zeros(2*N*T,2*N*T);
            k{N*(n_obs+1)*T+i} = C(N*n_obs*T+i,:)';
            dd{N*(n_obs+1)*T+i} = -d(N*n_obs*T+i,:);
        end
    end

    options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H),'Display','off'); % off/iter-detailed
    fun = @(x)quadobj(x,Q,f,0);
    nonlconstr = @(x)quadconstr(x,H,k,dd,opts);
    
    try
        [x,fval,eflag,output,lambda] = fmincon(fun,x0,...
            [],[],[],[],[],[],nonlconstr,options);
    catch
        keyboard
    end
end