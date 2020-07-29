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
rbfs = 10;
sigma = logspace(log10(.3),log10(.002), bfs); %original
% sigma = 0.002*ones(1,bfs); %this gives bf #10 the same size as rmps
% sigma = logspace(log10(0.205),log10(0.002),bfs); %closest overall
% sigma = logspace(log10(0.3),log10(0.0005),bfs);

%dmp init
D = struct('tau',1,'alphay',4,'betay',1,'alphax',.5,... %how fast cannonical system decays
    'bfs',bfs,... %number of basis functions (psi activations)
    'c',logspace(log10(1), log10(.01), bfs),... %basis function centers 
    'sigma',sigma,... %basis function widths
    'w',zeros(1,length(desiredtraj)),...
    'f',zeros(S.Nt,2),...
    'psi',nan(S.Nt,bfs),'y',y0,'ydot',y0dot,'x',x0);

%rmp init
R = struct('tau',1,'alphay',4,'betay',1,'alpha',10,'beta',1,'nu',1.2,'bfs',rbfs,...
    'rho',zeros(rbfs),...
    'a0',[1,0.1,10^-9*ones(1,rbfs-2)],...
    'a',[],...
    'ep',10^-9*ones(1,rbfs),...
    'w',zeros(1,length(desiredtraj)),...
    'f',zeros(S.Nt,2),'y',y0,'ydot',y0dot);

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

R.a(1,:) = R.a0;
alpha_rho = R.alpha*ones(R.bfs,1);
beta_rho = R.beta*ones(R.bfs,1);
nu_rho = R.nu*ones(R.bfs,1);
% nu_rho(2) = 2.4;
R.rho = zeros(R.bfs);

for i = 1:R.bfs
    for j = 1:R.bfs
        if i==j
            R.rho(i,j) = alpha_rho(i)/beta_rho(i);
        else
            if 0 == j-1-i
                R.rho(i,j) = (alpha_rho(i) - alpha_rho(j)/nu_rho(j))./beta_rho(j);
            else
                R.rho(i,j) = (alpha_rho(i) + alpha_rho(j))/beta_rho(j);
            end
        end
    end
end

index = 1;
expected_SHCas = R.a;
a = R.a(1,:);

for t = desiredt(1): S.dt: desiredt(end)    
    dW = sqrt(S.dt)*randn(1,R.bfs);
    da = a .* (R.alpha - a*R.rho) *S.dt + R.ep.*dW; %tau for scaling?
    a = max(min(a + da, 1), .0005); %enforcing boundaries on 'a'
    
    if t>=desiredt(index)
        index = index+1;
    else
        expected_SHCas(index, :) = a;
    end
end

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
%     for i = 1:Ntrials
%         if rand(1)>.5
%             % try random stuff
%             D.w = 1-.5*rand(2, D.bfs);
%             R.w = 1-.5*rand(2, R.bfs);
%             [dcost,rcost,D,R] = Eval(D,R,S);
%             if dcost<bestdcost
%                bestdw = D.w;
%                bestdcost = dcost;
%                delta = 1;
%             end
%             if rcost<bestrcost
%                bestrw = R.w;
%                bestrcost = rcost;
%                delta = 1;
%             end
%         else
%             %try different random stuff based on size delta 
%            D.w = bestdw + 3*delta*randn(2, D.bfs);
%            R.w = bestrw + 3*delta*randn(2, R.bfs);
%            [dcost,rcost,D,R] = Eval(D,R,S);
%            if dcost<bestdcost
%                bestdw = D.w;
%                bestdcost = dcost;
%            end
%            if rcost<bestrcost
%                bestrw = R.w;
%                bestrcost = rcost;
%            end
%         end
%     end
%     D.w = bestdw;
%     dcost = bestdcost;
%     R.w = bestrw;
%     rcost = bestrcost;
%     delta = delta*.99;
%     [~,~,D,R] = Eval(D,R,S);
%     
%     alldcosts = [alldcosts,bestdcost];
%     allrcosts = [allrcosts,bestrcost];
    
%Gradient Descent on the best weights
pause(1)
chosenew = Inf;
chosedel = Inf;
while chosenew+chosedel>=2
    chosenew = 0;
    chosedel = 0;
     for i = 1:size(D.w,1)
         for j = 1:size(D.w,2)
             wtrial = bestdw;
             wnew = bestdw;
             delt = delta*(1-2*rand(1));
             wtrial(i,j) = bestdw(i,j) + delt;
             [dcosttrial,~,D,~] =  Eval(D,R,S,wtrial);
             wnew(i,j) = bestdw(i,j) + max(min(bestdcost*delt/(bestdcost - dcosttrial), delta*5), -delta*5);
             [dcostnew,~,D,~] =  Eval(D,R,S,wnew);
             if dcostnew < bestdcost
                 bestdw = wnew;
                 bestdcost = dcostnew;
                 chosenew = chosenew +1;
             end
             if dcosttrial < bestdcost
                 bestdw = wtrial;
                 bestdcost = dcosttrial;
                 chosedel = chosedel +1;
             end
         end
     end
end

%at the end of one pass, save progress
alldcosts = [alldcosts, bestdcost];
D.w = bestdw;
dcost = bestdcost;
 
% Gradient Descent on the best weights
pause(1)
chosenew = Inf;
chosedel = Inf;
while chosenew+chosedel>=2
    chosenew = 0;
    chosedel = 0;
     for i = 1:size(R.w,1)
         for j = 1:size(R.w,2)
             wtrial = bestrw;
             wnew = bestrw;
             delt = delta*(1-2*rand(1));
             wtrial(i,j) = bestrw(i,j) + delt;
             [~,rcosttrial,~,R] =  Eval(D,R,S,wtrial);
             wnew(i,j) = bestrw(i,j) + max(min(bestrcost*delt/(bestrcost - rcosttrial), delta*5), -delta*5);
             [~,rcostnew,~,R] =  Eval(D,R,S,wnew);
             if rcostnew < bestrcost
                 bestrw = wnew;
                 bestrcost = rcostnew;
                 chosenew = chosenew +1;
             end
             if rcosttrial < bestrcost
                 bestrw = wtrial;
                 bestrcost = rcosttrial;
                 chosedel = chosedel +1;
             end
         end
     end
end

% at the end of one pass, save progress
allrcosts = [allrcosts, bestrcost];
R.w = bestrw;
rcost = bestrcost;
% ShowPlot(D,R,S)
end
Eval(D,R,S);
ShowPlot(D,R,S)

figure(8)
hold on
title('All Costs')
plot(alldcosts,'-k','Linewidth',1.5)
plot(allrcosts,'--k','Linewidth',1.5)
legend('DMP','RMP')
xlabel('Iterations')
ylabel('Cost')

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
    D.f(i,:) = (D.psi(i,:)*wdmp') *D.x(i+1,:) /sum(D.psi(i,:));
%     f = forcing(i,:);
    
    ydotdot = (D.alphay*(D.betay*(S.g-D.y(i,:))-D.ydot(i,:))+D.f(i,:))/D.tau;
    D.ydot(i+1,:) = ydotdot*S.dt + D.ydot(i,:);
    D.y(i+1,:) = D.ydot(i+1,:)*S.dt+D.y(i,:);
    
    trajpointdistancessq = min ((D.y(i+1,1) - S.desiredtraj(:,1)).^2 + (D.y(i+1,2) - S.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);

    [dsquared,~] = min( (D.y(i+1,1) - S.desiredtraj(:,1)).^2 + ...
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
dcomponents = [1*distance, 10*timetogoal/S.Nt, 1*sum(sqrt(trajpointdistancessq))];
dcost = 1*distance + 10*timetogoal/S.Nt + 1*sum(sqrt(trajpointdistancessq));

%RMP
distance = 0;
timetogoal = nan;
trajpointdistancessq = Inf*ones(size(S.desiredtraj, 1), 1);
for i = 1:S.Nt
    dW = sqrt(S.dt)*randn(1,R.bfs);
    if i == 1
        da = R.a0 .* (R.alpha - R.a0*R.rho) *S.dt + R.ep.*dW;
        R.a(i,:) = max(min(R.a0 + da, 1), 0.0005);
    else
        da = R.a(i-1,:) .* (R.alpha - R.a(i-1,:)*R.rho) *S.dt + R.ep.*dW; %tau for scaling?
        R.a(i,:) = max(min(R.a(i-1,:) + da, 1), .0005); %enforcing boundaries on a
    end
    R.f(i,:) = (R.a(i,:)*wrmp')/sum(R.a(i,:));

    ydotdot = (R.alphay*(R.betay*(S.g-R.y(i,:))-R.ydot(i,:))+R.f(i,:))/R.tau;
    R.ydot(i+1,:) = ydotdot*S.dt + R.ydot(i,:);
    R.y(i+1,:) = R.ydot(i+1,:)*S.dt+R.y(i,:);
    
    trajpointdistancessq = min ((R.y(i+1,1) - S.desiredtraj(:,1)).^2 + (R.y(i+1,2) - S.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);

    [dsquared,~] = min( (R.y(i+1,1) - S.desiredtraj(:,1)).^2 + ...
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
rcomponents = [1*distance, 10*timetogoal/S.Nt, 1*sum(sqrt(trajpointdistancessq))];
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
psicolors = ones(1,size(D.w,2))'*[0.7, 0.7, 0.7];
acolors = ones(1, size(R.w,2))'*[0.7, 0.7, 0.7];
psicolors(dmpcolorstouse,:) = cool(length(dmpcolorstouse));
acolors(rmpcolorstouse,:) = cool(length(rmpcolorstouse));

figure(1)
 clf(1)
 figure(1)
subplot(2,1,1)
title ('Basis functions times weights (dmp)')
for j = 1:size(D.psi,2)
    hold on
    plot((1:S.Nt)*S.dt, psiall(:,j)*D.w(1,j), '-c', 'Color', psicolors(j,:))
    plot((1:S.Nt)*S.dt, psiall(:,j)*D.w(2,j), '-c.', 'Color', psicolors(j,:))
end
subplot(2,1,2)
title('Basis functions times weights (rmp)')
for j = 1:size(R.a,2)
    hold on
    plot((1:S.Nt)*S.dt, aall(:,j)*R.w(1,j), '-c', 'Color', acolors(j,:))
    plot((1:S.Nt)*S.dt, aall(:,j)*R.w(2,j), '-c.', 'Color', acolors(j,:))
end

figure(2)
 clf(2)
 figure(2)
subplot(2,1,1)
hold on
title('Weights plotted spatially (dmp)')
plot(D.w(1,:),D.w(2,:), '-ro', 'Color', [.6, .6, .6])
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
    plot(D.w(1,i), D.w(2,i), 'k*', 'Color', psicolors(i,:))
end
axis equal
subplot(2,1,2)
hold on
title('Weights plotted spatially (rmp)')
plot(R.w(1,:),R.w(2,:), '-ro', 'Color', [.6, .6, .6])
for i = 1:size(R.w, 2)
    plot(yallrmp(rmpcolor==(1+i),1), yallrmp(rmpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
    plot(R.w(1,i), R.w(2,i), 'k*', 'Color', acolors(i,:))
end
axis equal

figure(3)
 clf(3)
 figure(3)
subplot(2,1,1)
hold on
title('Desired traj (black) and produced traj (color) - dmp')
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
end
plot(S.desiredtraj(:,1),S.desiredtraj(:,2), '-k.');
axis equal
subplot(2,1,2)
hold on
title('Desired traj (black) and produced traj (color) - rmp')
for i = 1:size(R.w, 2)
    plot(yallrmp(rmpcolor==(1+i),1), yallrmp(rmpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
end
plot(S.desiredtraj(:,1),S.desiredtraj(:,2), '-k.');
axis equal

figure(4)
 clf(4)
 figure(4)
subplot(2,1,1)
hold on
title('DMP forcing function')
plot((1:S.Nt)*S.dt,D.f(:,1),'-r')
plot((1:S.Nt)*S.dt,D.f(:,2),'--r')
legend('x-component','y-component')
subplot(2,1,2)
hold on
title('RMP forcing function')
plot((1:S.Nt)*S.dt,R.f(:,1),'-b')
plot((1:S.Nt)*S.dt,R.f(:,2),'--b')
legend('x-component','y-component')

figure(5)
 clf(5)
 figure(5)
subplot(2,1,1)
hold on
title('Unweighted basis functions (dmp)')
for i = 1:size(psiall,2)
    plot((1:S.Nt)*S.dt, psiall(:,i),'-c', 'Color', psicolors(i,:))
end
subplot(2,1,2)
hold on
title('Unweighted basis functions (rmp)')
for i = 1:size(aall,2)
    plot((1:S.Nt)*S.dt, aall(:,i),'-c','Color',acolors(i,:))
end

figure(6)
 clf(6)
 figure(6)
subplot(2,1,1)
hold on
title('DMP Basis Functions')
plot((1:S.Nt)*S.dt, psiall(:,5),'-k','Linewidth', 1.5)
set(gca,'xtick',[])
subplot(2,1,2)
hold on
title('RMP Basis Functions')
plot((1:S.Nt)*S.dt,aall(:,5),'-k','Linewidth',1.5)
% legend('DMP','RMP')
set(gca,'xtick',[])
% dmparea5 = trapz(psiall(:,5))
% rmparea5 = trapz(aall(:,5))

figure(7)
 clf(7)
 figure(7)
hold on
title('Basis Functions')
plot((1:S.Nt)*S.dt, psiall(:,5),'-k','Linewidth', 1.5)
plot((1:S.Nt)*S.dt,aall(:,5),'--k','Linewidth',1.5)
legend('DMP','RMP')
set(gca,'xtick',[])

%Looking for area under curve
% dmparea = trapz(psiall,1)
% rmparea = trapz(aall,1)

end

