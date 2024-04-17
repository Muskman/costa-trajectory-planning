function [f, energy, constraint, xk] = costscaalgo(mu,kb,cb,w,xk,std_noise,K,opts)

dt = opts.dt;
T = opts.T; 
N = opts.n_agents;
x_obs = opts.x_obs;
r_obs = opts.r_obs; r_a = opts.r_a;
[C, d] = obs_constraint(xk,x_obs,r_obs+r_a,opts);

e = std_noise*randn(size(xk)); % generate realization of noise
e0 = zeros(size(xk));

[g, grad_f] = disturbance(xk,e,dt,opts); % disturbance prediction and gradient of objective function
[g0, grad_f0] = disturbance(xk,e0,dt,opts);

z = grad_f;

grad_f_old = zeros(size(xk));

f = zeros(K,1);
energy = zeros(K,1);
constraint = zeros(K,1);
L = opts.goal_weight*kron(diag([1/T:1/T:1-1/T]),eye(2));

for k = 1:K   
    for j = 1:N
        f(k) = f(k) + norm(sqrt(L)*(xk(2*(j-1)*T+1:2*j*T-2)-repmat(opts.x_goal(2*j-1:2*j),T-1,1)))+norm(sqrt(eye(2*T-2)-L)*(xk(2*(j-1)*T+3:2*j*T)-xk(2*(j-1)*T+1:2*j*T-2)-dt*g0(2*(j-1)*T+1:2*j*T-2)))^2;
        energy(k) = energy(k) + norm((xk(2*(j-1)*T+3:2*j*T)-xk(2*(j-1)*T+1:2*j*T-2)-dt*g0(2*(j-1)*T+1:2*j*T-2))/dt)^2;
    end

    try
        constraint(k) = sum(C*xk-d>0);
    catch
        % for debugging
        keyboard
    end

    w = w + mean((grad_f).^2);    
    eta = kb*(w^(-1/3));
    bt = cb.*(eta.^2);   

    z = grad_f + (1-bt).*(z - grad_f_old); 
    b = (2*mu*xk - z)/mu;

    P = opts.D'; P(1:2,1:2) = zeros(2);

    x = qcqp(b,C,d,P,g0,(opts.v*dt)^2,xk,opts);
    
    % next sample
    e = std_noise*randn(size(xk)); % generate realization of noise
    
    % CoSTORM SCA
    [gn, grad_f_old] = disturbance(xk,e,dt,opts);  % wrt next sample   
    xk = (1-eta)*xk + eta*x; % xk(1:2) = x_start; xk(end-1:end) = x_goal;
    [C, d] = obs_constraint(xk,x_obs,r_obs+r_a,opts);
    [g, grad_f] = disturbance(xk,e,dt,opts); % disturbance prediction and gradient of objective function
    [g0, grad_f0] = disturbance(xk,e0,dt,opts);
    
    % print optimization stats
    % fprintf('k: %3.0f, CoSTA: %2.6f, constraint: %d, eta: %1.6f, bt: %1.6f \n', k, f(k), constraint(k), norm(eta), norm(bt))
    
end
