%June 2020
%DMP and RMP structures using simple trajectory following (no point mass)

close all
clear all
clc
tic

rng(1) %seeding different random numbers to check validity of cost evaluation

%TRAJECTORY CHOICES
% %ramp up bubble out -----
% desiredtraj = [0:.05:.5, .5*cos((pi/2):-.1:0)+.5;
%                0:.05:.5, .5*sin((pi/2):-.1:0)]'; 
% g = [1,0]; %goal
% %ramp up bubble in ------ 
% desiredtraj = [0:.05:.5, fliplr(.5*cos((-pi/2):-.1:-(pi))+1);
%               0:.05:.5, fliplr(.5*sin((-pi/2):-.1:-(pi))+.5)]';
% g = [1,0];
%ramp up wrap around -----
desiredtraj = [(0:.05:.5)*0, .5*cos(pi/2:.1:pi*2.6);
               0:.05:.5, .5*sin(pi/2:.1:pi*2.6)]';
g = [0,0];

%system init
S = struct('desiredtraj',desiredtraj,'dt',0.01,'Nt',1000,'g',g);
%trajectory start, velocity at trajectory start, canonical state vector,
%size of timestep, number of timesteps, goal
y0 = [];
y0(1,:) = .001 * [1,1];
y0dot = [];
y0dot(1,:) = 0 * [1,1];
x0 = 1;
bfs = 10;

%dmp init
D = struct('tau',1,'alphay',4,'betay',1,'alphax',.5,... %how fast cannonical system decays
    'bfs',bfs,... %number of basis functions (psi activations)
    'c',logspace(log10(1), log10(.01), bfs),... %basis function centers 
    'sigma',logspace(log10(.3),log10(.002), bfs),... %basis function widths
    'w',zeros(1,length(desiredtraj)),...
    'psi',nan(S.Nt,bfs),'y',y0,'ydot',y0dot,'x',x0);

%rmp init
R = struct('tau',1,'alphay',4,'betay',1,'alpha',10,'beta',1,'nu',1.2,'bfs',bfs,...
    'rho',zeros(bfs),...
    'a',[],...
    'ep',[10^-9*ones(1,bfs)],...
    'w',zeros(1,length(desiredtraj)),'y',y0,'ydot',y0dot);

%Init
% D.w = S.desiredtraj(round(linspace(1,length(S.desiredtraj),D.bfs)) ,:)'*(D.alphay+D.betay);
% R.w = S.desiredtraj(round(linspace(1,length(S.desiredtraj),R.bfs)) ,:)'*(R.alphay+R.betay);

Ntrials = 100;
desiredt =  linspace(0,10,length(S.desiredtraj))';% assumes constant velocity
desiredydot = [0,0;diff(S.desiredtraj)./(diff(desiredt)*[1,1])]; % desired speed
desiredydotdot = [0,0;diff(desiredydot)./(diff(desiredt)*[1,1])];
fdesiredd = D.tau * desiredydotdot - (D.alphay * (D.betay*(S.g-S.desiredtraj) - desiredydot));

expected_x = x0*exp(-D.alphax/D.tau*(desiredt-desiredt(1)));
expectedPsi = exp(ones(length(desiredt),1)*(-1./(2*D.sigma.^2)) .*(expected_x - D.c).^2); %n traj points by n basis 
% proveXexpworking % run this code to make a plot
betawd =  (expected_x .* expectedPsi) \ (sum(expectedPsi')'.* fdesiredd);
errors = (expected_x .* expectedPsi)*betawd - (sum(expectedPsi')'.* fdesiredd);

D.w = betawd';

R.a(1,:) = [1,0.1,10^-9*ones(1,R.bfs-2)];
alpha_rho = R.alpha*ones(R.bfs,1);
beta_rho = R.beta*ones(R.bfs,1);
nu_rho = R.nu*ones(R.bfs,1);
R.rho = zeros(bfs);

for i = 1:bfs
    for j = 1:bfs
        if i==j
            R.rho(i,j) = alpha_rho(i)/beta_rho(i);
        else if 0 == j-1-i
                R.rho(i,j) = (alpha_rho(i) - alpha_rho(j)/nu_rho(j))./beta_rho(j);
            else
                R.rho(i,j) = (alpha_rho(i) + alpha_rho(j))/beta_rho(j);
            end
        end
    end
end

index = 1;
expected_SHCas = R.a;
for t = desiredt(1): S.dt: desiredt(end)    
    dW = sqrt(S.dt)*randn(1,R.bfs);
    da = R.a .* (R.alpha - R.a*R.rho) *S.dt + R.ep.*dW; %tau for scaling?
    R.a = max(min(R.a + da, 1), .0005); %enforcing boundaries on 'a'
    
    if t>=desiredt(index)
        index = index+1;
    else
        expected_SHCas(index, :) = R.a;
    end
end
R.a
pause
size(expected_SHCas)
fdesiredr = R.tau * desiredydotdot - (R.alphay * (R.betay*(S.g-S.desiredtraj) - desiredydot));
betawr =  (expected_SHCas) \ (sum(expected_SHCas')'.* fdesiredr);
% meanerror = sum(sum(errors.^2))
R.w = betawr';

'before optimization'
[dcost,rcost,D,R] = Eval(D,R,S);
ShowPlot(D,R,S)
pause

%DMP weight evaluation
delta = 1;
alldcosts = [];
allrcosts = [];
for passes = 1:10
    'enter optimization'
    passes
    bestdw = D.w;
    bestdcost = dcost;
    bestrw = R.w;
    bestrcost = rcost;
    for i = 1:Ntrials
        if rand(1)>.5
            % try random stuff
            D.w = 1-.5*rand(2, D.bfs);
            R.w = 1-.5*rand(2, R.bfs);
            [dcost,rcost,D,R] = Eval(D,R,S);
            ShowPlot(D,R,S)
            if dcost<bestdcost
               bestdw = D.w;
               bestdcost = dcost;
               delta = 1;
            end
            if rcost<bestrcost
               bestrw = R.w;
               bestrcost = rcost;
               delta = 1;
            end
        else
            %try different random stuff based on size delta 
           D.w = bestdw + 3*delta*randn(2, D.bfs);
           R.w = bestrw + 3*delta*randn(2, R.bfs);
           [dcost,rcost,D,R] = Eval(D,R,S);
           ShowPlot(D,R,S)
           if dcost<bestdcost
               bestdw = D.w;
               bestdcost = dcost;
           end
           if rcost<bestrcost
               bestrw = R.w;
               bestrcost = rcost;
           end
        end
    end
    D.w = bestdw;
    dcost = bestdcost;
    R.w = bestrw;
    rcost = bestrcost;
    delta = delta*.99;
%     [~,~,D,R] = Eval(D,R,S);
    
    alldcosts = [alldcosts,bestdcost];
    allrcosts = [allrcosts,bestrcost];
    
%     %Gradient Descent on the best weights
%     pause(1)
%     chosenew = Inf;
%     chosedel = Inf;
%     while chosenew+chosedel>=3
%     chosenew = 0;
%     chosedel = 0;
%      for i = 1:size(D.w,1)
%          for j = 1:size(D.w,2)
%              wtrial = bestw;
%              wnew = bestw;
%              delt = delta*(1-2*rand(1));
%              wtrial(i,j) = bestw(i,j) + delt;
%              [dcosttrial,~,D,~] =  Eval(D,R,S,wtrial);
%              ShowPlot(D,R,S)
%              wnew(i,j) = bestw(i,j) + max(min(bestcost*delt/(bestcost - dcosttrial), delta*5), -delta*5);
%              [dcostnew,~,D,~] =  Eval(D,R,S,wnew);
%              if dcostnew < bestcost
%                  bestw = wnew;
%                  bestcost = dcostnew;
%                  chosenew = chosenew +1;
%              end
%              if dcosttrial < bestcost
%                  bestw = wtrial;
%                  bestcost = dcosttrial;
%                  chosedel = chosedel +1;
%              end
%          end
%      end
%     end

%     %at the end of one pass, save progress
%     alldcosts = [alldcosts, bestcost];
%     D.w = bestw;
%     dcost = bestcost;

% %RMP weight evaluation
 
%     Gradient Descent on the best weights
%     pause(1)
%     chosenew = Inf;
%     chosedel = Inf;
%     while chosenew+chosedel>=3
%     chosenew = 0;
%     chosedel = 0;
%      for i = 1:size(R.w,1)
%          for j = 1:size(R.w,2)
%              wtrial = bestw;
%              wnew = bestw;
%              delt = delta*(1-2*rand(1));
%              wtrial(i,j) = bestw(i,j) + delt;
%              [~,rcosttrial,~,R] =  Eval(D,R,S,wtrial);
%              wnew(i,j) = bestw(i,j) + max(min(bestcost*delt/(bestcost - rcosttrial), delta*5), -delta*5);
%              [~,rcostnew,~,R] =  Eval(D,R,S,wnew);
%              if rcostnew < bestcost
%                  bestw = wnew;
%                  bestcost = rcostnew;
%                  chosenew = chosenew +1;
%              end
%              if rcosttrial < bestcost
%                  bestw = wtrial;
%                  bestcost = rcosttrial;
%                  chosedel = chosedel +1;
%              end
%          end
%      end
%     end
% 
%     at the end of one pass, save progress
%     allrcosts = [allrcosts, bestcost];
%     R.w = bestw;
%     rcost = bestcost;
ShowPlot(D,R,S)
end
ShowPlot(D,R,S)
Eval(D,R,S);
toc

%% Evaluate Cost
function [dcost,rcost,D,R] = Eval(D,R,S,wnew)
if nargin == 4
    wdmp = wnew;
    wrmp = wnew;
end
if nargin == 3
    wdmp = D.w;
    wrmp = R.w;
end

goaltol = 0.005;
distance = 0;
maxtime = (S.Nt+1)*S.dt;

%DMP
timetogoal = nan;
trajpointdistancessq = Inf*ones(size(S.desiredtraj, 1), 1);
for i = 1:S.Nt
    D.x(i+1,:) = D.x(i,:) + ((-D.alphax)*D.x(i,:))/D.tau*S.dt;
    D.psi(i,:) = exp(-1./(2*D.sigma.^2).*((D.x(i+1,:)-D.c).^2));
    f(i,:) = (D.psi(i,:)*wdmp') *D.x(i+1,:) /sum(D.psi(i,:));
%     f = forcing(i,:);
    
    ydotdot = (D.alphay*(D.betay*(S.g-D.y(i,:))-D.ydot(i,:))+f(i,:))/D.tau;
    D.ydot(i+1,:) = ydotdot*S.dt + D.ydot(i,:);
    D.y(i+1,:) = D.ydot(i+1,:)*S.dt+D.y(i,:);
    
    trajpointdistancessq = min ((D.y(i+1,1) - S.desiredtraj(:,1)).^2 + (D.y(i+1,2) - S.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);

    [dsquared, whichtrajpoint] = min( (D.y(i+1,1) - S.desiredtraj(:,1)).^2 + ...
        (D.y(i+1,2) - S.desiredtraj(:,2)).^2);
    distance = distance + sqrt(dsquared);
    
    %record the time if it's at the goal
    if isnan(timetogoal) && sum((D.y(i+1,:)-S.g).^2)<goaltol^2
        timetogoal = i*S.dt;
    end
end

if isnan(timetogoal) 
    timetogoal = maxtime-S.dt;
end
dcost = 1*distance + 10*timetogoal/S.Nt + 1*sum(sqrt(trajpointdistancessq));

%RMP
timetogoal = nan;
trajpointdistancessq = Inf*ones(size(S.desiredtraj, 1), 1);
for i = 1:S.Nt
    dW = sqrt(S.dt)*randn(1,R.bfs);
    da = R.a(i,:) .* (R.alpha - R.a(i,:)*R.rho) *S.dt + R.ep.*dW; %tau for scaling?
    R.a(i+1,:) = max(min(R.a(i,:) + da, 1), .0005); %enforcing boundaries on a
    f(i,:) = R.a(i+1,:)*wrmp'/sum(R.a(i+1,:));

    ydotdot = (R.alphay*(R.betay*(S.g-R.y(i,:))-R.ydot(i,:))+f(i,:))/R.tau;
    R.ydot(i+1,:) = ydotdot*S.dt + R.ydot(i,:);
    R.y(i+1,:) = R.ydot(i+1,:)*S.dt+R.y(i,:);
    
    trajpointdistancessq = min ((R.y(i+1,1) - S.desiredtraj(:,1)).^2 + (R.y(i+1,2) - S.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);

    [dsquared, whichtrajpoint] = min( (R.y(i+1,1) - S.desiredtraj(:,1)).^2 + ...
        (R.y(i+1,2) - S.desiredtraj(:,2)).^2);
    distance = distance + sqrt(dsquared);
    
    %record the time if it's at the goal
    if isnan(timetogoal) && sum((R.y(i+1,:)-S.g).^2)<goaltol^2
        timetogoal = i*S.dt;
    end
end

if isnan(timetogoal) 
    timetogoal = maxtime-S.dt;
end

rcost = 1*distance + 10*timetogoal/S.Nt + 1*sum(sqrt(trajpointdistancessq));
end

%% Plot
function ShowPlot(D,R,S)

xall = nan*ones(S.Nt,1);
psiall =  nan*ones(S.Nt,length(D.sigma));
aall = nan*ones(S.Nt, R.bfs);
yall = nan*ones(S.Nt,2);
ydotall = nan*ones(S.Nt,2);
whichcolor = nan*ones(S.Nt,1);
rmpcolor = whichcolor;
dmpcolor = whichcolor;
alldistanceerror = nan*ones(S.Nt,1);
whichtraj = nan*ones(S.Nt,1);
allverror = nan*ones(S.Nt,1);

% desiredtrajdir = diff(S.desiredtraj);
% desiredtrajdir = .5*( [desiredtrajdir;desiredtrajdir(end,:)] + [desiredtrajdir(1,:);desiredtrajdir]);%centerdiffs
% desiredtrajdir = desiredtrajdir ./sqrt(sum(desiredtrajdir.^2,2)); %normalize

for i = 1:S.Nt
    [~, dmpcolor(i)] = max( [sum( abs(D.alphay*(D.betay*(S.g-D.y(i,:))-D.ydot(i,:))).^2)*.1;
                               sum((abs((D.psi(i,:).*D.w) *D.x(i,:) /sum(D.psi(i,:)))').^2,2)]);
    [~, rmpcolor(i)] = max( [sum( abs(R.alphay*(R.betay*(S.g-R.y(i,:))-R.ydot(i,:))).^2)*.1;
                               sum((abs((R.a(i,:).*R.w) /sum(R.a(i,:)))').^2,2)]);
%     [dsquared, whichtrajpoint] = min( (D.y(1) - S.desiredtraj(:,1)).^2 + ...
%                                 (D.y(2) - S.desiredtraj(:,2)).^2);
%     [dsquared, whichtrajpoint] = min( (R.y(1) - S.desiredtraj(:,1)).^2 + ...
%                                 (R.y(2) - S.desiredtraj(:,2)).^2);
%     desv = desiredtrajdir(whichtrajpoint, :);
%     vnormalerror = (desv(1)*S.ydot(2) - desv(2)*S.ydot(1))/sqrt(S.ydot*S.ydot');
%     alldistanceerror(i) = sqrt(dsquared);
%     whichtraj(i) = whichtrajpoint;
%     allverror(i) = vnormalerror;
    yalldmp(i,:) = D.y(i,:);
    yallrmp(i,:) = R.y(i,:);
%     ydotall(i,:) = S.ydot;
%     xall(i) = S.x;
    psiall(i,:) = D.psi(i,:);
    aall(i,:) = R.a(i,:);
end
dmpcolorstouse = setdiff(unique(dmpcolor)-1, 0);
rmpcolorstouse = setdiff(unique(rmpcolor)-1, 0);
psicolors = ones(1, size(D.w,2))'*[0.7, 0.7, 0.7];
acolors = ones(1, size(R.w,2))'*[0.7, 0.7, 0.7];
psicolors(dmpcolorstouse,:) = cool(length(dmpcolorstouse));
acolors(rmpcolorstouse,:) = cool(length(rmpcolorstouse));


figure(1)
 clf(1)
 figure(1)
hold on
subplot(2,1,1)
title ('Basis functions times weights (dmp)')
for j = 1:size(D.psi,2)
    plot((1:S.Nt)*S.dt, psiall(:,j)*D.w(1,j), '-c', 'Color', psicolors(j,:))
    plot((1:S.Nt)*S.dt, psiall(:,j)*D.w(2,j), '-c.', 'Color', psicolors(j,:))
end
subplot(2,1,2)
title('Basis functions times weights (rmp)')
for j = 1:size(R.a,2)
    plot((1:S.Nt)*S.dt, aall(:,j)*R.w(1,j), '-c', 'Color', acolors(j,:))
    plot((1:S.Nt)*S.dt, aall(:,j)*R.w(2,j), '-c.', 'Color', acolors(j,:))
end

figure(2)
 clf(2)
 figure(2)
hold on
subplot(2,1,1)
plot(yalldmp(:,1), yalldmp(:,2), 'k-')
plot(S.desiredtraj(:,1), S.desiredtraj(:,2), '-g.');
plot(yalldmp(dmpcolor==1,1), yalldmp(dmpcolor==1,2), 'k.')
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
    plot(D.w(1,i), D.w(2,i), 'k*', 'Color', psicolors(i,:))
end
axis equal
subplot(2,1,2)
plot(yallrmp(:,1), yallrmp(:,2), 'k-')
plot(S.desiredtraj(:,1), S.desiredtraj(:,2), '-g.');
plot(yallrmp(rmpcolor==1,1), yallrmp(rmpcolor==1,2), 'k.')
for i = 1:size(R.w, 2)
    plot(yallrmp(rmpcolor==(1+i),1), yallrmp(rmpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
    plot(R.w(1,i), R.w(2,i), 'k*', 'Color', acolors(i,:))
end
axis equal

figure(3)
 clf(3)
 figure(3)
hold on
subplot(2,1,1)
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
end
plot(S.desiredtraj(:,1),S.desiredtraj(:,2), '-k.');
axis equal
subplot(2,1,2)
for i = 1:size(R.w, 2)
    plot(yallrmp(rmpcolor==(1+i),1), yallrmp(rmpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
end
plot(S.desiredtraj(:,1),S.desiredtraj(:,2), '-k.');
axis equal




end

