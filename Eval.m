function d = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma_or_rho, c_or_alpha, alphay, betay,g, desiredtraj, ShowPlot, OneforDMPtwoforRMP, w1,w2,Nsaddles)
% %--------------------
goaltol = .005;
%ShowPlot = 1;
% % initialize or variable descriptions (basics that are the same for RMP/DMP)
%  Nt = 1000; % number of timesteps
%  y0 = .001 * [1,1]; %system initial position
%  y0dot= 0 * [1,1]; %system initial velocity
%   dt = .01; %timestep 
% tau = 1;   %time constant base system, also reused in canonical system        
% alphay = 4; %damping base system
% betay = alphay/4; %stiffness base system
% alphax = .5; %how fast cannonical system decays
% g = [1,0]; %goal position of base system
%  desiredtraj = [0:.05:.5, fliplr(.5*cos((-pi/2):-.1:-(pi))+1);
%                 0:.05:.5, fliplr(.5*sin((-pi/2):-.1:-(pi))+.5)]'; 

% %shared naming
% x0 = 1; %cannonical system initial state (scalar for DMP, vector for RMP)
% w  = weights that are learned  
    %for DMP w will be size 2 by Npsi (assuming planar trajectory)
% w  =   [0.8726   -4.8226   -5.2870   -5.6955   11.9854   -1.6538    0.2696   -0.6114    2.2657   -0.8382;
%        4.3342    1.8618    0.2612   -1.5936   -5.3375    5.2797   -3.1289    0.0051    0.4761   -1.3480];
% for RMP w will be size 
%
% 

% Nsaddles = 20;
if OneforDMPtwoforRMP == 1
% %for DMP
%   Npsi = 10; %number of basis functions (use this generate c and sigma)
 % c = logspace(log10(1), log10(.01), Npsi); % these are the positions of DMP in X, and they do not change. 
 % sigma = 2*logspace(log10(.3),log10(.002), Npsi);
 sigma = sigma_or_rho;
 c = c_or_alpha;
 
 x = x0;
 if ShowPlot
xall = nan*ones(Nt,1);
psiall =  nan*ones(Nt,length(sigma));
 end
 
 else
% %for RMP
% Nsaddles = 10;
alpha = 10;
beta = 1;
nu = 1.2;
rho = calculaterho_notacycle(alpha*ones(Nsaddles,1), beta*ones(Nsaddles,1), nu*ones(Nsaddles,1));
 
 
ep = 10^-9*ones(1,Nsaddles);
a0 = ep;
a0(1) = 1;
a0(2) = .01;
a = a0;
if ShowPlot
aall = nan*ones(Nt, Nsaddles);
end
% %%rho = []; %coupling matrix
% %%alpha = []; %growth factors
% % for now P and zd are zero since not measuring contact
% % and C is zero because there is no coupling and no measured z
% rho = sigma_or_rho;
% alpha = c_or_alpha;
% 
end

distance = 0;

ydot = y0dot;
y = y0;
fnum = 0;

maxtime = (Nt+1)*dt;
timetogoal = nan;

if ShowPlot
yall = nan*ones(Nt,2);
ydotall = nan*ones(Nt,2);
whichcolor = nan*ones(Nt,1);
alldistanceerror = nan*ones(Nt,1);
whichtraj = nan*ones(Nt,1);
allverror = nan*ones(Nt,1);
end

desiredtrajdir = diff(desiredtraj);
desiredtrajdir = .5*( [desiredtrajdir;desiredtrajdir(end,:)] + [desiredtrajdir(1,:);desiredtrajdir]);%centerdiffs
desiredtrajdir = desiredtrajdir ./sqrt(sum(desiredtrajdir.^2,2)); %normalize

trajpointdistancessq = Inf*ones(size(desiredtraj, 1), 1);

for i = 1:Nt
    % advance diff equation
    %DMP
    if OneforDMPtwoforRMP == 1        
%         w = [w(1,:); zeros(1,length(w))];
%         w = [zeros(1,length(w));w(2,:)];
    x = x + ((-alphax)*x)/tau*dt;
    psi = exp(-1./(2*sigma.^2).*(x-c).^2);
%     fnum = fnum + ((psi*w')*x);
    f = (psi*w') *x /sum(psi);
%     f = [0,f(2)];
%     f = [f(1),0];
%     %new weights
%         beta2 = ((alphay*(betay*(g-y)-ydot))*(dt^2)/tau) + ydot*dt + y;
%         beta1 = (psi*w'*(dt^2))*x/tau;
%         f = beta2; 
    end
    %RMP
    if OneforDMPtwoforRMP == 2
    dW = sqrt(dt)*randn(1,Nsaddles);
    da = a .* (alpha - a*rho) *dt + ep.*dW; %%%DID I forget tau??
    a = max(min(a + da, 1), .0005);
    f = a*w'/sum(a);
    end
    ydotdot = (alphay*(betay*(g-y)-ydot)+f)/tau;
    ydot = ydotdot*dt + ydot;
    y = ydot*dt+y;    
    
    
    if ShowPlot
        % find out which part of solution has max contribution
        if OneforDMPtwoforRMP == 1
        [~, whichcolor(i)] = max( [sum( abs(alphay*(betay*(g-y)-ydot)).^2)*.1;
                                   sum((abs((psi.*w) *x /sum(psi))').^2,2)]);
        else 
            [~, whichcolor(i)] = max( [sum( abs(alphay*(betay*(g-y)-ydot)).^2)*.1;
                                       sum((abs((a.*w)/sum(a))').^2,2)]);
                                  
        end
        
    end
    trajpointdistancessq = min ((y(1) - desiredtraj(:,1)).^2 + (y(2) - desiredtraj(:,2)).^2, ...
                                trajpointdistancessq);
       
    [dsquared, whichtrajpoint] = min( (y(1) - desiredtraj(:,1)).^2 + ...
        (y(2) - desiredtraj(:,2)).^2);
    distance = distance + ...
        sqrt( dsquared);
    desv = desiredtrajdir(whichtrajpoint, :);
    vnormalerror = (desv(1)*ydot(2) - desv(2)*ydot(1))/sqrt(ydot*ydot');
    
    % here if its at the goal, I record time
    if isnan(timetogoal) && sum((y-g).^2)<goaltol^2
        timetogoal = i*dt;
    end
    
    if ShowPlot
        alldistanceerror(i) = sqrt( dsquared);
        whichtraj(i) = whichtrajpoint;
        allverror(i) = vnormalerror;
        yall(i,:) = y;
        ydotall(i,:) = ydot;
        if OneforDMPtwoforRMP == 1
        xall(i) = x;
        psiall(i,:) = psi;
        end
        if OneforDMPtwoforRMP == 2
            if Nsaddles == 1
                aall = a;
            else
                aall(i,:) = a;
            end
        end
    end
end

% %evaluate trajectory
% distance = nan*ones(Nt,1);
% whichtraj = nan*ones(Nt,1);
% for i = 1:Nt
%     distance(i) = sqrt( min( (yall(i,1) - desiredtraj(:,1)).^2 + ... 
%                 (yall(i,2) - desiredtraj(:,2)).^2));
%     [~,whichtraj(i)] = min( (yall(i,1) - desiredtraj(:,1)).^2 + ... 
%                 (yall(i,2) - desiredtraj(:,2)).^2);
% end

if isnan(timetogoal) 
    timetogoal = maxtime-dt;
end

d = 1*distance+...
    10*timetogoal/Nt+...
    1*sum(sqrt(trajpointdistancessq));


if ShowPlot
%grayout colors that aren't used
colorstouse = setdiff(unique(whichcolor(1:i)-1), 0);
% colorstouse = setdiff(unique(whichcolor(1:(timetogoal/dt)))-1, 0);
colorstouse = setdiff(unique(whichcolor)-1, 0);

psicolors = ones(1, size(w,2))'*[0.7, 0.7, 0.7];
psicolors(colorstouse,:) = cool(length(colorstouse));


figure(3)
clf(3)
figure(3)
hold on 
title ('Basis functions times weights')
if OneforDMPtwoforRMP == 1
    for j = 1:size(psi,2)
        plot((1:Nt)*dt, psiall(:,j)*w(1,j), '-c', 'Color', psicolors(j,:))
        plot((1:Nt)*dt, psiall(:,j)*w(2,j), '-c.', 'Color', psicolors(j,:))
    end
else
    for j = 1:size(a,2)
        plot((1:Nt)*dt, aall(:,j)*w(1,j), '-c', 'Color', psicolors(j,:))
        plot((1:Nt)*dt, aall(:,j)*w(2,j), '-c.', 'Color', psicolors(j,:))
    end
end

figure(2)
 clf(2)
 figure(2)
hold on
plot(yall(:,1), yall(:,2), 'k-')
plot(desiredtraj(:,1), desiredtraj(:,2), '-g.');
plot(w2(1,:), w2(2,:), '-r.')
plot(w(1,:), w(2,:), '-ro', 'Color', [.6, .6, .6])

plot(yall(whichcolor==1,1), yall(whichcolor==1,2), 'k.')
for i = 1:size(w, 2)
    plot(yall(whichcolor==(1+i),1), yall(whichcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
    plot(w(1,i), w(2,i), 'k*', 'Color', psicolors(i,:))
end
axis equal

figure(5)
 clf(5)
 figure(5)
hold on
for i = 1:size(w, 2)
    plot(yall(whichcolor==(1+i),1), yall(whichcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
end

plot(desiredtraj(:,1), desiredtraj(:,2), '-k.');
axis equal

% figure(6)
% clf(6)
% figure(6)
% 
% hold on 
% for j = 1:size(psi,2)
%     plot((1:Nt)*dt, psiall(:,j), '-c', 'Color', psicolors(j,:))
%     plot((1:Nt)*dt, psiall(:,j), '-c.', 'Color', psicolors(j,:))
% end
% saddlecolors  = cool(Nsaddles);
% for j = 1:Nsaddles
%     plot((1:Nt)*dt, 1.1+aall(:,j), '-c', 'Color', saddlecolors(j,:))
%     plot((1:Nt)*dt, 1.1+aall(:,j), '-c.', 'Color', saddlecolors(j,:))
% end
end

