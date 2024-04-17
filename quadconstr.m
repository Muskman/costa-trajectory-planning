function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d,opts)

T = opts.T;
N = opts.n_agents;

jj = length(H); % jj is the number of inequality constraints

y = zeros(1,jj);
for i = 1:jj
    try
        y(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
    catch
        keyboard
    end
end

yeq = zeros(1,4*N);
for j = 1:N
    yeq(4*j-3:4*j) = [x(2*(j-1)*T+1:2*(j-1)*T+2) - opts.x_start(2*j-1:2*j); x(2*j*T-1:2*j*T) - opts.x_goal(2*j-1:2*j)]';
end
    
if nargout > 2
    grady = zeros(length(x),jj);
    for i = 1:jj
        grady(:,i) = H{i}*x + k{i};
    end
end

gradyeq = zeros(length(x),4*N);
for j = 1:N
    gradyeq(2*(j-1)*T+1,4*j-3) = 1;
    gradyeq(2*(j-1)*T+2,4*j-2) = 1;
    gradyeq(2*j*T-1,4*j-1) = 1;
    gradyeq(2*j*T,4*j) = 1;
end

end