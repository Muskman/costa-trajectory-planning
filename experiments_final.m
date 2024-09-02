clc
close all
clear all

opts.env = 'test';

n_agents = 2;
n_obs = 1;
va_max = 1;         % va_max = 8.5; max_velocity of agent

m_scale = 1;        % scaling for environment size
c_scale = 0.8;      % scaling for current magnitude

nMC = 10;           % number of monte carlo runs 10
K  = 100;           % number of iterations

% generate environment with currents for plotting                
[X_lim,Y_lim] = meshgrid([-2.5:0.15:2.5].*m_scale);
Z = c_scale * X_lim .* exp(-(X_lim.^2 + Y_lim.^2)/m_scale^2)/m_scale;
[U,V] = gradient(Z,0.15,0.15);
mag = sqrt(U.^2+V.^2);

plotMap = true; plotEnergyCV = true; plotTime = true; 
test = 0;
if ~test
    N = [30 60];
    thetas = [-80 0 40 80];
    Tfs = [10 15 20];
    noise_levels = [0.05 0.1 0.2];
    saveFigs = true;
    saveMat = true;
else
    nMC = 1;
    N = [30];
    thetas = [80];
    Tfs = [15];
    noise_levels = 0;
    saveFigs = false;
    saveMat = false;
end

lo = length(thetas);
lp = length(N);
lq = length(Tfs);
lr = length(noise_levels);

res.SCA = zeros(lo,lr*lq*lp);
res.SCA_var = zeros(lo,lr*lq*lp);
res.CoSTORM = zeros(lo,lr*lq*lp);
res.CoSTORM_var = zeros(lo,lr*lq*lp);

res.SCA_cv = zeros(lo,lr*lq*lp);
res.SCA_cv_var = zeros(lo,lr*lq*lp);
res.CoSTORM_cv = zeros(lo,lr*lq*lp);
res.CoSTORM_cv_var = zeros(lo,lr*lq*lp);

res.SCA_t = zeros(lo,lr*lq*lp);
res.SCA_t_var = zeros(lo,lr*lq*lp);
res.CoSTORM_t = zeros(lo,lr*lq*lp);
res.CoSTORM_t_var = zeros(lo,lr*lq*lp);

cBar = lo*lp*lq*lr; cProg = 0;
oPB = textprogressbar(cBar , 'barlength', 20, ...
'updatestep', 1, ...
'startmsg', 'Overall Progress ',...
'endmsg', ' Done', ...
'showbar', true, ...
'showremtime', true, ...
'showactualnum', true, ...
'barsymbol', '+', ...
'emptybarsymbol', '-');

for o = 1:length(thetas)
    for p = 1:length(N)
        for q = 1:length(Tfs)
            for r = 1:length(noise_levels)
                oPB(cProg); 
                fprintf('\n');
                cProg = cProg+1;
                close all

                x_offset = 0; y_offset = 0;
                x_obs_offset = 0; y_obs_offset = 0;

                theta = thetas(o);
                R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

                T = N(p);        % number of waypoints
                Tf = Tfs(q);     % final time
                dt = Tf/T;       % spatial time difference between waypoints
                
                std_noise = noise_levels(r); % std of noise in disturbance perdictor
                
                % start, goal location and obstacle parameters
                x_start = zeros(2*n_agents,1); x_goal = zeros(2*n_agents,1);
                x_start(1:2) = (R*[-1; 1.7] + [x_offset; y_offset]).*m_scale; x_goal(1:2) = (R*[.9; -1.5] + [x_offset; y_offset]).*m_scale;
                if n_agents>1
                    thetaS = 0:90/(n_agents-1):90;
                    for i = 2:n_agents
                        Rs = [cosd(thetaS(i)) -sind(thetaS(i)); sind(thetaS(i)) cosd(thetaS(i))];
                        x_start(2*i-1:2*i) = Rs*x_start(1:2); x_goal(2*i-1:2*i) = Rs*x_goal(1:2); 
                    end
                end
                if n_obs == 1
                    x_obs = (R*[0.25; -0.5]+ [x_offset; y_offset] + [x_obs_offset; y_obs_offset]).*m_scale;
                    r_obs = 0.7*m_scale;
                elseif n_obs ==2 
                    x_obs = (R*[0.25 -0.25; -0.5 0.5]+ [x_offset; y_offset] + [x_obs_offset; y_obs_offset]).*m_scale;
                    r_obs = [0.6*m_scale; 0.2*m_scale];
                else
                    x_obs = (R*[0.25 -0.5; -0.25 0.5; 0.5 1]'+ [x_offset; y_offset] + [x_obs_offset; y_obs_offset]).*m_scale;
                    r_obs = [0.6*m_scale; 0.2*m_scale; 0.2*m_scale];
                end

                % tuning parameters  
                % CoSTORM SCA 
                mu = 0.001;   % strong convexity paramater
                kb = .365;      
                cb = 490;
                w = 700;

                % CSSCA algo
                mu1 = 0.001; 
                gam = 0.025;
                bet = 0.18;

                % straignt line initialization
                x_sl = [];
                for j = 1:n_agents
                    x_sl_1 = x_start(2*j-1):(x_goal(2*j-1)-x_start(2*j-1))/T:x_goal(2*j-1);
                    x_sl_2 = x_start(2*j):(x_goal(2*j)-x_start(2*j))/T:x_goal(2*j);
                    x_sl_temp = [x_sl_1; x_sl_2];
                    x_sl_temp(:,T/2) = []; 
                    x_sl_temp = reshape(x_sl_temp,2*T,1);
                    x_sl_temp(1:2) = x_start(2*j-1:2*j);
                    x_sl_temp(end-1:end) = x_goal(2*j-1:2*j);
                    x_sl = [x_sl; x_sl_temp];
                end

                % passing parameters to algorithm
                % environment parameters
                opts.c_scale = c_scale;
                opts.goal_weight = 0;
                opts.X_lim = X_lim; opts.Y_lim = Y_lim;
                opts.mag = mag;
                opts.x_start = x_start; opts.x_goal = x_goal;
                opts.x_obs = x_obs; opts.r_obs = r_obs;
                r_a = 0.1;
                opts.v = va_max-3*std_noise*c_scale;
                opts.r_a = r_a;     % agent radious
                opts.m_scale = m_scale;
                opts.dt = dt; opts.T = T; 
                opts.n_agents = n_agents; opts.n_obs = n_obs;

                % aux variables
                [A, B, D] = finitediff(T,opts);
                opts.A = A; opts.B = B; opts.D = D;

                % energy metrics
                energy_cost_mean = zeros(K,1); energy_cs_mean = zeros(K,1);
                energy_cost_var = zeros(K,1); energy_cs_var = zeros(K,1);
                constraint_cost_mean = zeros(K,1); constraint_cs_mean = zeros(K,1);
                constraint_cost_var = zeros(K,1); constraint_cs_var = zeros(K,1);
                t_cost_mean = zeros(K,1); t_cs_mean = zeros(K,1);
                t_cost_var = zeros(K,1); t_cs_var = zeros(K,1);

                % generate feasible initial guess
                opts.full_opt = 0;
                [~, ~, ~, x_init] = minalgo(x_sl,std_noise,K,opts);
                xk = x_init;
                
                % run algorithms
                iPB = textprogressbar(nMC, 'barlength', 20, ...
                    'updatestep', 1, ...
                    'startmsg', 'MC Progress      ',...
                    'endmsg', ' Done', ...
                    'showbar', true, ...
                    'showremtime', true, ...
                    'showactualnum', true, ...
                    'barsymbol', '+', ...
                    'emptybarsymbol', '-');
                for n = 1:nMC    
                    % tStart = tic; fprintf('Running CSSCA...')
                    [f_cs, energy_cs, constraint_cs, x_cs, t_cs] = csscaalgo(mu1,gam,bet,xk,std_noise,K,opts);
                    % t_cs = toc(tStart)-t_cost; fprintf('Done\n'); fprintf('\nRunning CoSTA...')
                    [f_cost, energy_cost, constraint_cost, x_cost, t_cost] = costscaalgo(mu,kb,cb,w,xk,std_noise,K,opts);
                    % t_cost = toc(tStart); fprintf('Done\n'); fprintf('Running IPM...')
                    % [f_min, energy_min, constraint_min, x_min] = minalgo(xk,std_noise,K,opts);
                    % t_min = toc(tStart)-t_cost-t_cs; fprintf('Done\n')

                    energy_cost_mean_old = energy_cost_mean; energy_cs_mean_old = energy_cs_mean; 
                    energy_cost_mean = ((n-1)*energy_cost_mean_old+energy_cost)/n;
                    energy_cs_mean = ((n-1)*energy_cs_mean_old+energy_cs)/n;

                    energy_cost_var = ((n-1)*energy_cost_var + (energy_cost-energy_cost_mean_old).*(energy_cost-energy_cost_mean))/n;
                    energy_cs_var = ((n-1)*energy_cs_var + (energy_cs-energy_cs_mean_old).*(energy_cs-energy_cs_mean))/n;

                    constraint_cost_mean_old = constraint_cost_mean; constraint_cs_mean_old = constraint_cs_mean; 
                    constraint_cost_mean = ((n-1)*constraint_cost_mean_old+constraint_cost)/n;
                    constraint_cs_mean = ((n-1)*constraint_cs_mean_old+constraint_cs)/n;

                    constraint_cost_var = ((n-1)*constraint_cost_var + (constraint_cost-constraint_cost_mean_old).*(constraint_cost-constraint_cost_mean))/n;
                    constraint_cs_var = ((n-1)*constraint_cs_var + (constraint_cs-constraint_cs_mean_old).*(constraint_cs-constraint_cs_mean))/n;
                    
                    t_cost_mean_old = t_cost_mean; t_cs_mean_old = t_cs_mean; 
                    t_cost_mean = ((n-1)*t_cost_mean_old+t_cost)/n;
                    t_cs_mean = ((n-1)*t_cs_mean_old+t_cs)/n;

                    t_cost_var = ((n-1)*t_cost_var + (t_cost-t_cost_mean_old).*(t_cost-t_cost_mean))/n;
                    t_cs_var = ((n-1)*t_cs_var + (t_cs-t_cs_mean_old).*(t_cs-t_cs_mean))/n;

                    iPB(n)
                end

                idx = lr*lq*(p-1)+lr*(q-1)+r;

                res.CoSTORM(o,idx) = energy_cost_mean(end); res.CoSTORM_var(o,idx) = energy_cost_var(end);
                res.SCA(o,idx) = energy_cs_mean(end); res.SCA_var(o,idx) = energy_cs_var(end);
                
                res.CoSTORM_cv(o,idx) = constraint_cost_mean(end); res.CoSTORM_cv_var(o,idx) = constraint_cost_var(end);
                res.SCA_cv(o,idx) = constraint_cs_mean(end); res.SCA_cv_var(o,idx) = constraint_cs_var(end);

                res.CoSTORM_t(o,idx) = t_cost_mean(end); res.CoSTORM_t_var(o,idx) = t_cost_var(end);
                res.SCA_t(o,idx) = t_cs_mean(end); res.SCA_t_var(o,idx) = t_cs_var(end);

                % opts.full_opt = 1;
                % t_startMin = tic; 
                % [~, ~, ~, x_minCon] = minalgo(x_init,std_noise,K,opts);
                % t_minCon = toc(t_startMin);

                if plotMap
                    figure; 
                    contourf(X_lim,Y_lim,mag,'LineColor','none','HandleVisibility','off');
                    
                    str = '#047495'; % #03719c (oceanblue) % #047495 (seablue)
                    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                    colmap = flipud(cmap('b',100,50,5));
                    caxis([0,max(max(mag))]); colormap (colmap); %parula
                    
                    htb = colorbar ;
                    htb.Label.String = 'Magnitude of ocean currents';
                    htb.Label.FontWeight = 'bold';
                    hold on;
                    quiver(X_lim,Y_lim,U,V,'k','LineWidth',0.8)
                    
                    % plotting trajectory
                    hold on; 
                    for i = 1:n_obs
                        circles(x_obs(1,i), x_obs(2,i), r_obs(i), 'linewidth',1.5 , 'facecolor','#BC9354','HandleVisibility','off')
                    end
                    for j = 1:n_agents
                        plot(x_init(2*(j-1)*T+1:2:2*j*T-1),x_init(2*(j-1)*T+2:2:2*j*T),'LineWidth',2,'Color','k','LineStyle','--')
                        plot(x_cost(2*(j-1)*T+1:2:2*j*T-1),x_cost(2*(j-1)*T+2:2:2*j*T),'LineWidth',2,'Color',[0, 1, 0])
                        plot(x_cs(2*(j-1)*T+1:2:2*j*T-1),x_cs(2*(j-1)*T+2:2:2*j*T),'LineWidth',2,'Color','m')
                        % plot(x_minCon(2*(j-1)*T+1:2:2*j*T-1),x_minCon(2*(j-1)*T+2:2:2*j*T),'LineWidth',2,'Color','r')
                        if j==1
                            shape_str = "square";
                        else
                            shape_str = "diamond";
                        end
                        plot(x_start(2*j-1),x_start(2*j),shape_str,'MarkerEdgeColor','k','MarkerFaceColor',"y",'MarkerSize',10)
                        plot(x_goal(2*j-1),x_goal(2*j),shape_str,'MarkerEdgeColor','k','MarkerFaceColor',"#A2142F",'MarkerSize',10)

                        x1 = x_cost((2*j-1)*T-1); y1 = x_cost((2*j-1)*T); u1 = diff(x_cost((2*j-1)*T-1:2:(2*j-1)*T+1)); v1 = diff(x_cost((2*j-1)*T:2:(2*j-1)*T+2));
                        x2 = x_cs((2*j-1)*T-1); y2 = x_cs((2*j-1)*T); u2 = diff(x_cs((2*j-1)*T-1:2:(2*j-1)*T+1)); v2 = diff(x_cs((2*j-1)*T:2:(2*j-1)*T+2));
                        q1 = quiver(x1,y1,u1,v1,'LineWidth',10,'Color',[0, 1, 0],'HandleVisibility','off'); 
                        q2 = quiver(x2,y2,u2,v2,'LineWidth',10,'Color','m','HandleVisibility','off');

                        axis equal
                        xlabel('Distance (m)','fontweight','bold'); ylabel('Distance (m)','fontweight','bold');
                        
                        pause(0.1)

                        q1.NodeChildren(2).Visible = 'off';
                        q2.NodeChildren(2).Visible = 'off';
                    end
                    if true %(thetas(o) == 0 || thetas(o) == 80)
                        legend({'Currents','Initial guess', 'CoSTA', 'CSSCA', 'Start','Goal', '', '', '',' ', ' ', ' ', ' '},'NumColumns',4,'fontweight','bold','fontsize',8,'Location','northwest'); legend boxoff
                        % legend({'Currents','Initial guess', 'CoSTA', 'CSSCA', 'IPM', 'Start','Goal', '', '','' , '',' ', ' ', ' ', ' '},'NumColumns',4,'fontweight','bold','fontsize',8,'Location','northwest'); legend boxoff
                    end
                    if saveFigs
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/map-T=',num2str(Tf),'noise=',num2str(std_noise),'.fig'),'fig')
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/map-T=',num2str(Tf),'noise=',num2str(std_noise),'.eps'),'epsc')
                    end
                end

                if plotEnergyCV
                    % plotting energy vs iters
                    figure;
                    % yyaxis left
                    plot(energy_cost_mean,'LineWidth',2,'Color','g','LineStyle', "-"); hold on
                    fill([1:K,fliplr(1:K)]', [energy_cost_mean+sqrt(energy_cost_var);flipud(energy_cost_mean-sqrt(energy_cost_var))],'b','FaceAlpha',0.3, 'LineStyle', "none" ,'HandleVisibility','off')
                    plot(energy_cs_mean,'LineWidth',2,'Color','m','LineStyle', "-")
                    fill([1:K,fliplr(1:K)]', [energy_cs_mean+sqrt(energy_cs_var);flipud(energy_cs_mean-sqrt(energy_cs_var))],'r','FaceAlpha',0.3, 'LineStyle', "none",'HandleVisibility','off')
                    xlabel('Iterations','fontweight','bold'); ylabel('Energy','fontweight','bold');
                    set(gca,'YColor','k');
                    % yyaxis right
                    % plot(constraint_cost_mean,'LineWidth',2,'Color','g','LineStyle', "--"); hold on
                    % fill([1:K,fliplr(1:K)]', [constraint_cost_mean+sqrt(constraint_cost_var);flipud(constraint_cost_mean-sqrt(constraint_cost_var))],'b','FaceAlpha',0.3, 'LineStyle', "none" ,'HandleVisibility','off')
                    % plot(constraint_cs_mean,'LineWidth',2,'Color','m','LineStyle', "--");
                    % fill([1:K,fliplr(1:K)]', [constraint_cs_mean+sqrt(constraint_cs_var);flipud(constraint_cs_mean-sqrt(constraint_cs_var))],'r','FaceAlpha',0.3, 'LineStyle', "none",'HandleVisibility','off')
                    % ylab = ylabel('Constraint Violation (CV)','fontweight','bold');
                    % ylab.Position(1) = 93;
                    % set(gca,'YColor','k');
                    % try
                    %     ylim([0 max([constraint_cost_mean; constraint_cs_mean])])
                    % catch
                    % end
                    if true %(thetas(o) == 0 || thetas(o) == 80)
                        % legend(["CoSTA Energy", "CSSCA Energy", "CoSTA CV", "CSSCA CV"],'Location', 'NorthEast','fontweight','bold','NumColumns',2); legend boxoff
                        legend(["CoSTA Energy", "CSSCA Energy"],'Location', 'NorthEast','fontweight','bold'); legend boxoff
                    end
                    if saveFigs
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/energy-cv-T=',num2str(Tf),'noise=',num2str(std_noise),'.fig'),'fig')
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/energy-cv-T=',num2str(Tf),'noise=',num2str(std_noise),'.eps'),'epsc')
                    end

                    % plotting evenergy vs time
                    figure;
                    % yyaxis left
                    plot(t_cost_mean,energy_cost_mean,'LineWidth',2,'Color','g','LineStyle', "-"); hold on
                    fill([t_cost_mean;flipud(t_cost_mean)], [energy_cost_mean+sqrt(energy_cost_var);flipud(energy_cost_mean-sqrt(energy_cost_var))],'b','FaceAlpha',0.3, 'LineStyle', "none" ,'HandleVisibility','off')
                    plot(t_cs_mean,energy_cs_mean,'LineWidth',2,'Color','m','LineStyle', "-")
                    fill([t_cs_mean;flipud(t_cs_mean)], [energy_cs_mean+sqrt(energy_cs_var);flipud(energy_cs_mean-sqrt(energy_cs_var))],'r','FaceAlpha',0.3, 'LineStyle', "none",'HandleVisibility','off')
                    xlabel('Time (s)','fontweight','bold'); ylabel('Energy','fontweight','bold');
                    set(gca,'YColor','k');
                    % yyaxis right
                    % plot(t_cost_mean, constraint_cost_mean,'LineWidth',2,'Color','g','LineStyle', "--"); hold on
                    % fill([t_cost_mean;flipud(t_cost_mean)], [constraint_cost_mean+sqrt(constraint_cost_var);flipud(constraint_cost_mean-sqrt(constraint_cost_var))],'b','FaceAlpha',0.3, 'LineStyle', "none" ,'HandleVisibility','off')
                    % plot(t_cs_mean,constraint_cs_mean,'LineWidth',2,'Color','m','LineStyle', "--");
                    % fill([t_cs_mean;flipud(t_cs_mean)], [constraint_cs_mean+sqrt(constraint_cs_var);flipud(constraint_cs_mean-sqrt(constraint_cs_var))],'r','FaceAlpha',0.3, 'LineStyle', "none",'HandleVisibility','off')
                    % ylab = ylabel('Constraint Violation (CV)','fontweight','bold');
                    % ylab.Position(1) = max([t_cost_mean; t_cs_mean])*0.95;
                    % set(gca,'YColor','k');
                    % try
                    %     ylim([0 max([constraint_cost_mean; constraint_cs_mean])])
                    % catch
                    % end
                    xlim([min([t_cost_mean; t_cs_mean]) max([t_cost_mean; t_cs_mean])])
                    if true %(thetas(o) == 0 || thetas(o) == 80)
                        % legend(["CoSTA Energy", "CSSCA Energy", "CoSTA CV", "CSSCA CV"],'Location', 'NorthEast','fontweight','bold','NumColumns',2); legend boxoff
                        legend(["CoSTA Energy", "CSSCA Energy"],'Location', 'NorthEast','fontweight','bold'); legend boxoff
                    end

                    if saveFigs
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/energy-cv-time-T=',num2str(Tf),'noise=',num2str(std_noise),'.fig'),'fig')
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/energy-cv-time-T=',num2str(Tf),'noise=',num2str(std_noise),'.eps'),'epsc')
                    end
                end

                if plotTime
                    % plot average time per iteration
                    figure
                    algo = categorical({'CSSCA','CoSTA'});
                    % algo = categorical({'CSSCA','CoSTA','IPM'});
                    % per iteration time
                    % data = [mean(diff(t_cost_mean)) mean(diff(t_cs_mean))];
                    % err = [std(diff(t_cost_mean)) std(diff(t_cs_mean))];
                    % total time required
                    data = [res.SCA_t/100 res.CoSTORM_t/100];
                    err = sqrt([res.SCA_t_var res.CoSTORM_t_var])/100;
                    % data = [res.SCA_t res.CoSTORM_t t_minCon];
                    % err = [res.SCA_t_var res.CoSTORM_t_var 0];
                    hBar = barh(algo,data); 
                    hold on; 
                    
                    x_e = hBar.YEndPoints;
                    errorbar(x_e,algo,err,'horizontal','k','linestyle','none','LineWidth',2)
                    ylabel('Algorithm', 'FontWeight','bold');
                    xlabel('Time per iteration (sec)')
                    title('Time', 'FontWeight','bold');
                    ax = gca;
                    ax.YAxis.FontWeight = 'bold';
                    ax.XAxis.FontWeight = 'bold';
                    ytickangle(90)

                    if saveFigs
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/avg-time-T=',num2str(Tf),'noise=',num2str(std_noise),'.fig'),'fig')
                        saveas(gcf,strcat('plots/N=',num2str(T),'/theta',num2str(theta),'/avg-time-T=',num2str(Tf),'noise=',num2str(std_noise),'.eps'),'epsc')
                    end
                end
            end
        end
    end
end


if saveMat
    save(strcat('plots/',date))
end

% fprintf('Time taken by Algorithms:\n')
% fprintf('CSSCA : %f\n',t_cs)
% fprintf('CoSTA : %f\n',t_cost)
% fprintf('IPM   : %f\n',t_min)
