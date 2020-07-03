close all
clear all
clc
tic
rng(1)

% finalds = [];
OneforDMPtwoforRMP = 1

Usepresents = 1;
dt = .01;
Nt = 1000;
g = [1,0];
y0 = .001 * [1,1];
y0dot= 0 * [1,1];
x0 = 1;

% %ramp up bubble out
% desiredtraj = [0:.05:.5, .5*cos((pi/2):-.1:0)+.5;
%                0:.05:.5, .5*sin((pi/2):-.1:0)]'; 
           
% % %ramp up bubble in ------ 
%    desiredtraj = [0:.05:.5, fliplr(.5*cos((-pi/2):-.1:-(pi))+1);
%                   0:.05:.5, fliplr(.5*sin((-pi/2):-.1:-(pi))+.5)]'; 
%-------             
% %ramp up wrap around
   desiredtraj = [(0:.05:.5)*0, .5*cos(pi/2:.1:pi*2.6);
                   0:.05:.5, .5*sin(pi/2:.1:pi*2.6)]';
  g = [0,0];
% % ----- 
%% IGNORE THIS
%  bestw =[  3.2763   -7.9183    2.9071   -7.7999  -19.5263   -4.2977    7.2047   -0.1165  -1.7832    2.3749
%           6.9902   -0.5676   -2.2115    2.5635   11.9405   -6.1655   -0.0285    0.6614 -0.6912    0.6141];
%  bestw = [   2.5279   -6.3904    0.4768   -0.9810   -4.9774    3.7569 -0.7490    0.0487   -0.1125    0.3058;
%             7.3744   -2.2085   -0.0588   -0.1458    0.6147   -0.3531   -0.0641    0.0472   -0.0146   -0.0088];
% % After running all night! with all three factors
% % bestw =  [-2.9178   -4.4224    0.7990   -0.0852
%            2.1801    0.4769   -0.2220   -0.0464];
% two factors (no direction)      
%%
tau = 1;
alphay = 4; % if I decrease these numbers weird stuff happens % damping and stiffness
betay = alphay/4; %alphay/4; % if I decrease these numbers weird stuff happens % stiffnesses
alphax = .5; %how fast cannonical system decays
% Npsi = Nsaddles;
Npsi = 10;
c = logspace(log10(1), log10(.01), Npsi); % these are the positions of DMP in X, and they do not change. 
sigma = logspace(log10(.3),log10(.002), Npsi); %if I make this narrower (smaller sigma) it looks like it can get smaller turn radius, but best solution has a kink
%Cz = 0 no coupling

%-------------------
Ktrials = 500;
% w = 5-10*rand(2, Npsi); % these are the weights - they do change.
%Bestw for ramp up wrap around
%  w =[0.1682   -0.0909   -2.5292   -2.0515   -0.4892    1.4525    2.4341    1.8060   -0.0140   -2.3589
%     1.8293    2.7509    0.6173   -1.1173   -2.4039   -1.9603   -0.0670    1.5410    2.4332    1.5410];
% w =[    0.0990   -0.3065   -2.3560   -2.0671   -0.4662    1.4635    2.3638    1.8060   -0.0140   -1.7901;
%     2.0122    2.6248    0.6573   -1.1035   -2.4039   -1.9603   -0.2377    1.6481    2.3701    1.7008];

%What happens when w = desired traj (simple sampling)
% w = desiredtraj(round(linspace(1,length(desiredtraj),Npsi)) ,:)'*(alphay+betay); %when alphay = 4, betay = 1nice to multiply by 5
% w =   [  0.2060   -0.3065   -2.4163   -2.0671   -0.4646    1.4635    2.2877    1.8060   -0.1399   -2.0545;
%     1.7634    2.6421    0.6093   -1.1801   -2.4039   -1.8916   -0.2974    1.6481    2.4332    1.5173];
w = desiredtraj(round(linspace(1,length(desiredtraj),Npsi)) ,:)';
w1 = w;
w = desiredtraj(round(linspace(1,length(desiredtraj),Npsi)) ,:)'*(alphay+betay);
w2 = w;

%500 trials
% wbest =   [  0.1682   -0.0909   -2.5292   -2.0515   -0.4892    1.4525    2.4341    1.8060   -0.0140   -2.3589;
%     1.8293    2.7509    0.6173   -1.1173   -2.4039   -1.9603   -0.0670    1.5410    2.4332    1.5410];
%1000 trials
% wbest =  [   0.2060   -0.3065   -2.4163   -2.0671   -0.4646    1.4635    2.2877    1.8060   -0.1399   -2.0545;
%     1.7634    2.6421    0.6093   -1.1801   -2.4039   -1.8916   -0.2974    1.6481    2.4332    1.5173];

if Usepresents
%%%%%%% $$$ Here I am inserting linear algebra approach to weight finding
% where we know f = tau*ydotdot - (alphay*(betay*(g-y) - ydot);
desiredt =  linspace(0,10,length(desiredtraj))';% [1:length(desiredtraj)]'; % assumes constant velocity linspace(0,10,length(desiredtraj))
desiredydot = [0,0;diff(desiredtraj)./(diff(desiredt)*[1,1])]; % desired speeds
desiredydotdot = [0,0;diff(desiredydot)./(diff(desiredt)*[1,1])];
fdesired = tau * desiredydotdot - (alphay * (betay*(g-desiredtraj) - desiredydot));

if OneforDMPtwoforRMP == 1
    expected_x = x0*exp(-alphax/tau*(desiredt-desiredt(1)));
    expectedPsi = exp(ones(length(desiredt),1)*(-1./(2*sigma.^2)) .*(expected_x - c).^2); %n traj points by n basis 
    
    % proveXexpworking % run this code to make a plot
    betaw =  (expected_x .* expectedPsi) \ (sum(expectedPsi')'.* fdesired);
    errors = (expected_x .* expectedPsi)*betaw - (sum(expectedPsi')'.* fdesired);
else
    
    
Nsaddles = 10; %DUPLICATED IN EVAL
alpha = 10; %DUPLICATED IN EVAL
beta = 1;  %DUPLICATED IN EVAL
nu = 1.2; %DUPLICATED IN EVAL
rho = calculaterho_notacycle(alpha*ones(Nsaddles,1), beta*ones(Nsaddles,1), nu*ones(Nsaddles,1));
ep = 10^-9*ones(1,Nsaddles);  %DUPLICATED IN EVAL
a0 = ep;  %DUPLICATED IN EVAL
a0(1) = 1; %DUPLICATED IN EVAL
a0(2) = .01; %DUPLICATED IN EVAL
a = a0;
index = 1;
expected_SHCas = a;

for t = desiredt(1): dt: desiredt(end)  
    dW = sqrt(dt)*randn(1,Nsaddles);
    da = a .* (alpha - a*rho) *dt + ep.*dW; %wondering if I forgot tau
    a = max(min(a + da, 1), .0005);
    
    if t>=desiredt(index)
        index = index+1;
    else
        expected_SHCas(index, :) = a;
    end
end


    size(expected_SHCas)
    betaw =  ( expected_SHCas) \ (sum(expected_SHCas')'.* fdesired);
end


% meanerror = sum(sum(errors.^2))
w = betaw';
%%%%%%% $$$ end insertion
end


%TWO main functions here, one to evaluate weights
[d,dcomponents] = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1, OneforDMPtwoforRMP);

delta = 1;
allds = [];
dcomponentsplot = [];
for passes = 1:10
bestw = w;
bestd = d;
for i = 1:Ktrials
    if rand(1)>.5
        % try random stuff
       w = 1-.5*rand(2, Npsi);
       [d,dcomponents] = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP);
       if d<bestd
           bestw = w;
           bestd = d;
           delta = 1;
       end
    else
        %try different random stuff based on size delta 
       w = bestw + 3*delta*randn(2, Npsi);
       [d,dcomponents] = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP);
       if d<bestd
           bestw = w;
           bestd = d;
       end
       
    end

end
w = bestw;
d = bestd;
delta = delta*.99;
Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1,OneforDMPtwoforRMP);
allds = [allds,bestd];
dcomponentsplot = [dcomponentsplot;dcomponents];

% now for a gradient desecent on the best one. 
pause(1)
chosenew = Inf;
chosedel = Inf;
while chosenew+chosedel>=2
chosenew = 0;
chosedel = 0;
 for i = 1:size(w,1)
     for j = 1:size(w,2)
         wtrial = bestw;
         wnew = bestw;
         delt = delta*(1-2*rand(1));
         wtrial(i,j) = bestw(i,j) + delt;
         [dtrial,dtrialcomponents] =  Eval(wtrial, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP);
         wnew(i,j) = bestw(i,j) + max(min(bestd*delt/(bestd - dtrial), delta*5), -delta*5);
         [dnew,dnewcomponents] =  Eval(wnew, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP);
         if dnew < bestd
             bestw = wnew;
             bestd = dnew;
             dcomponents = dnewcomponents;
             chosenew = chosenew +1;
         end
         if dtrial < bestd
             bestw = wtrial;
             bestd = dtrial;
             dcomponents = dtrialcomponents;
             chosedel = chosedel +1;
         end
         
     end
 end
 chosenew
 chosedel
end

allds = [allds,bestd];
dcomponentsplot = [dcomponentsplot;dcomponents];
w = bestw;
d = bestd;
%ShowTraj(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj);
% if passes==1 || passes==5 || passes==10
%     pause
% end
end
Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1,OneforDMPtwoforRMP);

figure(4)
hold on
plot(allds,'k','Linewidth',2)
plot(dcomponentsplot(:,1),'b')
plot(dcomponentsplot(:,2),'m')
plot(dcomponentsplot(:,3),'r')
%  
% pause
% ShowTraj(w, Nt, y0dot, y0, x0, alphax, 2*tau, dt, sigma, c, alphay, betay,g, desiredtraj);
% w = weights
% example 
%     0.8726   -4.8226   -5.2870   -5.6955   11.9854   -1.6538    0.2696   -0.6114    2.2657   -0.8382
%     4.3342    1.8618    0.2612   -1.5936   -5.3375    5.2797   -3.1289    0.0051    0.4761   -1.3480
% Nt = number of timesteps
% y0 initial conditions
% x0 initial conditions
% alphax % keep same as training
% tau % changing this makes shape pretty different, sloppier if highter
%dt  %seems like I can change this without changing shape

% finalds = [finalds,allds(end)]
% figure(6)
% plot(finalds)
% end
% figure(7)
% hold on
% plot((1:Nt)*dt,ftrace(:,1),'Linewidth',2,'Color','r')
% plot((1:Nt)*dt,ftrace(:,2),'Linewidth',2,'Color','b')

toc

function [d,dcomponents,ftrace] = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma_or_rho, c_or_alpha, alphay, betay,g, desiredtraj, ShowPlot, OneforDMPtwoforRMP)
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
Nsaddles = 10;
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
% fnum = 0;
sumf = 0;
ftrace = [];

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

xall = [];
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
    da = a .* (alpha - a*rho) *dt + ep.*dW; %tau for scaling?
    a = max(min(a + da, 1), .0005); %enforcing boundaries on a
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
%     distance = distance + sqrt(dsquared);
distance = sqrt(dsquared);
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
            aall(i,:) = a;
        end
        fall(i,:) = f;
    end
    sumf = sumf + f;
    ftrace = [ftrace;f];
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
    1e5*timetogoal/Nt+...
    1*sum(sqrt(trajpointdistancessq));
dcomponents = [distance, timetogoal/Nt, sum(sqrt(trajpointdistancessq))];

if ShowPlot
%grayout colors that aren't used
colorstouse = setdiff(unique(whichcolor(1:i)-1), 0);
% colorstouse = setdiff(unique(whichcolor(1:(timetogoal/dt)))-1, 0);
colorstouse = setdiff(unique(whichcolor)-1, 0);

psicolors = ones(1, size(w,2))'*[0.7, 0.7, 0.7];
psicolors(colorstouse,:) = cool(length(colorstouse));

% majorpsi = exp(-1./(2*sum(sigma_or_rho.^2)).*(xall-sum(c_or_alpha)).^2);
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
    plot((1:Nt)*dt, fall(:,1), '-k')
    plot((1:Nt)*dt, fall(:,2), '-k.')
else
    for j = 1:size(a,2)
        plot((1:Nt)*dt, aall(:,j)*w(1,j), '-c', 'Color', psicolors(j,:))
        plot((1:Nt)*dt, aall(:,j)*w(2,j), '-c.', 'Color', psicolors(j,:))
    end
    plot((1:Nt)*dt, fall(:,1), '-k')
    plot((1:Nt)*dt, fall(:,2), '-k.')
end

figure(2)
 clf(2)
 figure(2)
hold on
plot(yall(:,1), yall(:,2), 'k-')
plot(desiredtraj(:,1), desiredtraj(:,2), '-g.');
% plot(w2(1,:), w2(2,:), '-r.')
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
end