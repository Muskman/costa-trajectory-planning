function [f, energy, constraint, xk] = minalgo(xk,std_noise,K,opts)

    T = opts.T; 
    N = opts.n_agents;
    x_obs = opts.x_obs;
    r_obs = opts.r_obs; r_a = opts.r_a;
        
    f = zeros(K,1);
    energy = zeros(K,1);
    constraint = zeros(K,1);
    
    if opts.full_opt == 0
        energy_fun = @(x) 0;% minObj(x,N,K,T,std_noise,opts);
    else
        energy_fun = @(x) minObj(x,N,K,T,std_noise,opts);
    end

    lb = -2.5*ones(2*N*T,1); ub = 2.5*ones(2*N*T,1);
    A = []; b = []; Aeq = []; beq = [];
    nonlcon = @(x) minConstr(x,N,T,x_obs,r_obs+r_a,opts);

    options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',4e3);%,'EnableFeasibilityMode',true,'SubproblemAlgorithm','cg');
    xk = fmincon(energy_fun,xk,A,b,Aeq,beq,lb,ub,nonlcon,options);
end





    