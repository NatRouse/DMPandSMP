%September 2020
%Dynamic movement primitive (DMP) and SHC-based movement primitive (SMP) 
%structures using simple trajectory following (no point mass).
%SMPs are based on a series of competitive Lotka-Volterra equations

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

% %ramp up wsap around -----
% desiredtraj = [(0:.05:.5)*0, .5*cos(pi/2:.1:pi*2.6);
%                0:.05:.5, .5*sin(pi/2:.1:pi*2.6)]';
% g = [0,0];

% %heart shape -----
% t = linspace(-pi,pi,78);
% desiredtraj = [16*sin(t).^3; 17 + 13*cos(t)- 5*cos(2*t) - 2*cos(3*t) - cos(4*t)]';
% desiredtraj = desiredtraj(2-end,:);
% g = desiredtraj(end,:);

% %ramp up bubble in adjusted ------ 
desiredtraj = [-1:.05:-.5, fliplr(.5*cos((-pi/2):-.1:-(pi))+1)-1;
              -1:.05:-.5, fliplr(.5*sin((-pi/2):-.1:-(pi))+0.5)-1]';
g = [2,0];


%SYSTEM PARAMETER INITIALIZATION
P = struct('desiredtraj',desiredtraj,'dt',0.001,'Nt',10000,'g',g);
%trajectory start, velocity at trajectory start, canonical state vector,
%size of timestep, number of timesteps, goal
y0 = [];
y0(1,:) = desiredtraj(1,:) + 0.001;
y0dot = [];
y0dot(1,:) = zeros(1,2);
x0 = 1;
bfs = 10;
sbfs = 10;
sigma = logspace(log10(.3),log10(.002), bfs); %original


%DMP INITIALIZATION
D = struct('tau',1,'alphay',4,'betay',1,'alphax',.5,... %how fast cannonical system decays
    'bfs',bfs,... %number of basis functions (psi activations)
    'c',logspace(log10(1), log10(.01), bfs),... %basis function centers 
    'sigma',sigma,... %basis function widths
    'w',zeros(2,bfs),...
    'f',zeros(P.Nt,2),...
    'psi',nan(P.Nt,bfs),'y',y0,'ydot',y0dot,'x',x0);

%SMP INITIALIZATION
S = struct('tau',1,'alphay',4,'betay',1,'alpha',10,'beta',1,'nu',1.2,'bfs',sbfs,...
    'rho',zeros(sbfs),...
    'a0',[1,0.1,10^-9*ones(1,sbfs-2)],...
    'a',[],...'ep',10^-9*ones(1,sbfs),...
    'ep',10^-9*ones(1,sbfs),...
    'w',zeros(2,bfs),...
    'f',zeros(P.Nt,2),'y',y0,'ydot',y0dot);


%INITIALIZE WEIGHTS
%Random
D.w = 5-10*rand(2,D.bfs);
S.w = 5-10*rand(2,S.bfs);
%Zero
% D.w = zeros(2,D.bfs);
% S.w = zeros(2,S.bfs);
%Scaled from trajectory
wd = P.desiredtraj(round(linspace(1,length(P.desiredtraj),D.bfs)) ,:)'...
    *(D.alphay+D.betay);
ws = P.desiredtraj(round(linspace(1,length(P.desiredtraj),S.bfs)) ,:)'...
    *(S.alphay+S.betay);
% D.w = wd;
% S.w = ws;


%BATCH LEARNING
Ntrials = 100;
desiredt = linspace(0,10,length(P.desiredtraj))';% assumes constant velocity
desiredydot = [0,0;diff(P.desiredtraj)./(diff(desiredt)*[1,1])]; % desired speed
desiredydotdot = [0,0;diff(desiredydot)./(diff(desiredt)*[1,1])];
fdesiredd = D.tau * desiredydotdot - (D.alphay * (D.betay*(P.g-P.desiredtraj) - desiredydot));

expected_x = x0*exp(-D.alphax/D.tau*desiredt);
expectedPsi = exp(ones(length(desiredt),1)*(-1./(2*D.sigma.^2)) .*(expected_x - D.c).^2);
betawd =  (expected_x .* expectedPsi) \ (sum(expectedPsi,2).* fdesiredd);
errors = (expected_x .* expectedPsi)*betawd - (sum(expectedPsi,2).* fdesiredd);

S.a(1,:) = S.a0;

%Parameters for rho matrix
alpha_rho = S.alpha*ones(S.bfs,1);
beta_rho = S.beta*ones(S.bfs,1);
nu_rho = S.nu*ones(S.bfs,1);
S.rho = zeros(S.bfs);

for i = 1:S.bfs
    for j = 1:S.bfs
        if i==j
            S.rho(i,j) = alpha_rho(i)/beta_rho(i);
        else
            if 0 == j-1-i % i == mod(j,S.bfs)+1
                S.rho(i,j) = (alpha_rho(i) - alpha_rho(j)/nu_rho(j))./beta_rho(j);
            else
                S.rho(i,j) = (alpha_rho(i) + alpha_rho(j))/beta_rho(j);
            end
        end
    end
end

index = 1;
expected_SHCas = S.a;
a = S.a(1,:);

%Run & plot for random weights
'random weights'
[dcostrand,rcostrand,D,R] = Eval(D,R,S);
ShowPlot(D,R,S,wd,ws)
alldcosts = dcostrand;
allrcosts = rcostrand;

for t = desiredt(1): P.dt: desiredt(end)    
    dW = sqrt(P.dt)*randn(1,S.bfs);
    da = a .* (S.alpha - a*S.rho) *P.dt + S.ep.*dW; %tau for scaling?
    a = max(min(a + da, 1), .0005); %enforcing boundaries on 'a'
    if t>=desiredt(index)
        index = index+1;
    else
        expected_SHCas(index, :) = a;
    end
end

fdesiredr = S.tau * desiredydotdot - (S.alphay *...
            (S.betay*(P.g-P.desiredtraj) - desiredydot));
betaws =  (expected_SHCas) \ fdesiredr;
meanerror = sum(sum(errors.^2))
S.w = betaws';
D.w = betawd';

%Run & plot after batch learning
'after batch learning'
[dcost,rcost,D,R] = Eval(D,R,S);
alldcosts = [alldcosts dcost];
allrcosts = [allrcosts rcost];
ShowPlot(D,R,S,wd,ws)
pause(5)

%Cost comparison (not necessary)
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
function [dcost,rcost,D,R] = Eval(D,S,P,wnew)
if nargin == 4
    wdmp = wnew;
    wsmp = wnew;
end
if nargin == 3
    wdmp = D.w;
    wsmp = S.w;
end
goaltol = 0.005;
maxtime = (P.Nt+1)*P.dt;

%DMP
distance = 0;
trajpointdistancessq = Inf*ones(size(P.desiredtraj, 1), 1);
ydotdotsum = 0;
for i = 1:P.Nt
    D.x(i+1,:) = D.x(i,:) + ((-D.alphax)*D.x(i,:))/D.tau*P.dt;
    D.psi(i,:) = exp(-1./(2*D.sigma.^2).*((D.x(i+1,:)-D.c).^2));
    D.f(i,:) = (D.psi(i,:)*wdmp') *D.x(i+1,:) /sum(D.psi(i,:));
    
    ydotdot = (D.alphay*(D.betay*(P.g-D.y(i,:))-D.ydot(i,:))+D.f(i,:))/D.tau;
    ydotdotsum = rssq(ydotdot) + ydotdotsum;
    D.ydot(i+1,:) = ydotdot*P.dt + D.ydot(i,:);
    D.y(i+1,:) = D.ydot(i+1,:)*P.dt+D.y(i,:);
    
    trajpointdistancessq = min ((D.y(i+1,1) - P.desiredtraj(:,1)).^2 +...
                            (D.y(i+1,2) - P.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);
    [dsquared,~] = min( (D.y(i+1,1) - P.desiredtraj(:,1)).^2 + ...
        (D.y(i+1,2) - P.desiredtraj(:,2)).^2);
    distance = distance + sqrt(dsquared);
end
dcost = 100*sum(sqrt((D.y(end,:)-P.g).^2)) ... %final goal error
    + sum(rssq(D.ydot,1))*P.dt ... %velocity
    + ydotdotsum*P.dt ... %acceleration
    + 1000*sum(sqrt(trajpointdistancessq)); %trajectory error

%SMP
distance = 0;
trajpointdistancessq = Inf*ones(size(P.desiredtraj, 1), 1);
ydotdotsum = 0;
for i = 1:P.Nt
    dW = sqrt(P.dt)*randn(1,S.bfs);
    if i == 1
        da = S.a0 .* (S.alpha - S.a0*S.rho) *P.dt + S.ep.*dW;
        S.a(i,:) = max(min(S.a0 + da, 1), 0.0005);
    else
        da = S.a(i-1,:) .* (S.alpha - S.a(i-1,:)*S.rho) *P.dt + S.ep.*dW;
        S.a(i,:) = max(min(S.a(i-1,:) + da, 1), .0005);
    end
    S.f(i,:) = (S.a(i,:)*wsmp');

    ydotdot = (S.alphay*(S.betay*(P.g-S.y(i,:))-S.ydot(i,:))+S.f(i,:))/S.tau;
    ydotdotsum = rssq(ydotdot) + ydotdotsum;
    S.ydot(i+1,:) = ydotdot*P.dt + S.ydot(i,:);
    S.y(i+1,:) = S.ydot(i+1,:)*P.dt+S.y(i,:);
    
    trajpointdistancessq = min ((S.y(i+1,1) - P.desiredtraj(:,1)).^2 +...
                            (S.y(i+1,2) - P.desiredtraj(:,2)).^2, ...
                            trajpointdistancessq);

    [dsquared,~] = min( (S.y(i+1,1) - P.desiredtraj(:,1)).^2 + ...
        (S.y(i+1,2) - P.desiredtraj(:,2)).^2);
    distance = distance + sqrt(dsquared);
end
rcost = 100*sum(sqrt((S.y(end,:)-P.g).^2)) ... %final goal error
    + sum(rssq(S.ydot,1))*P.dt ... %velocity
    + ydotdotsum*P.dt ... %acceleration
    + 1000*sum(sqrt(trajpointdistancessq)); %trajectory error
end

%% Plot
function ShowPlot(D,S,P,wd,ws)

xall = nan*ones(P.Nt,1);
psiall =  nan*ones(P.Nt,length(D.sigma));
aall = nan*ones(P.Nt, S.bfs);
yall = nan*ones(P.Nt,2);
ydotall = nan*ones(P.Nt,2);
whichcolor = nan*ones(P.Nt,1);
smpcolor = whichcolor;
dmpcolor = whichcolor;
alldistanceerror = nan*ones(P.Nt,1);
whichtraj = nan*ones(P.Nt,1);
allverror = nan*ones(P.Nt,1);

for i = 1:P.Nt
    [~, dmpcolor(i)] = max( [sum( abs(D.alphay*(D.betay*(P.g-D.y(i,:))-D.ydot(i,:))).^2)*.1;
                               sum((abs((D.psi(i,:).*D.w) *D.x(i,:) /sum(D.psi(i,:)))').^2,2)]);
    [~, smpcolor(i)] = max( [sum( abs(S.alphay*(S.betay*(P.g-S.y(i,:))-S.ydot(i,:))).^2)*.1;
                               sum((abs((S.a(i,:).*S.w) /sum(S.a(i,:)))').^2,2)]);
    yalldmp(i,:) = D.y(i,:);
    yallrmp(i,:) = S.y(i,:);
    psiall(i,:) = D.psi(i,:);
    aall(i,:) = S.a(i,:);
end

dmpcolorstouse = setdiff(unique(dmpcolor)-1, 0);
smpcolorstouse = setdiff(unique(smpcolor)-1, 0);
psicolors = ones(1,size(D.w,2))'*[0.7, 0.7, 0.7];
acolors = ones(1, size(S.w,2))'*[0.7, 0.7, 0.7];
psicolors(dmpcolorstouse,:) = cool(length(dmpcolorstouse));
acolors(smpcolorstouse,:) = cool(length(smpcolorstouse));

figure(1) %WEIGHTED KERNELS
 clf(1)
 figure(1)
subplot(1,2,1)
title ('DMP','Fontsize',20)
for j = 1:size(D.psi,2)
    hold on
    plot((1:P.Nt)*P.dt, psiall(:,j)*D.w(1,j), '-c', 'Color', psicolors(j,:))
    plot((1:P.Nt)*P.dt, psiall(:,j)*D.w(2,j), '-c.', 'Color', psicolors(j,:))
end
ylabel('canonical state, {\it x}','Fontsize',20)
xlabel('time (s)','Fontsize',20)
subplot(1,2,2)
title('SMP','Fontsize',20)
for j = 1:size(S.a,2)
    hold on
    plot((1:P.Nt)*P.dt, aall(:,j)*S.w(1,j), '-c', 'Color', acolors(j,:))
    plot((1:P.Nt)*P.dt, aall(:,j)*S.w(2,j), '-c.', 'Color', acolors(j,:))
end
plot((1:P.Nt)*P.dt,S.f(:,1),'-k')
plot((1:P.Nt)*P.dt,S.f(:,2),'--k')
ylabel('canonical state, {\it x}','Fontsize',20)
xlabel('time (s)','Fontsize',20)

figure(2) %SPATIAL WEIGHTS
 clf(2)
 figure(2)
subplot(1,2,1)
hold on
title('DMP')
plot(D.w(1,:),D.w(2,:), '-ro', 'Color', [.6, .6, .6])
plot(wd(1,:), wd(2,:), '-S.')
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
    plot(D.w(1,i), D.w(2,i), 'k*', 'Color', psicolors(i,:))
end
axis equal
% axes('Position',[.2 .5 .1 .1])
% box on
% for i = 1:size(D.w,2)
%     hold on
%     plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
% end
% set(gca,'YTick',[],'XTick',[])
subplot(1,2,2)
hold on
title('SMP')
plot(S.w(1,:),S.w(2,:), '-ro', 'Color', [.6, .6, .6])
plot(ws(1,:), ws(2,:), '-S.')
for i = 1:size(S.w, 2)
    plot(yallrmp(smpcolor==(1+i),1), yallrmp(smpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
    plot(S.w(1,i), S.w(2,i), 'k*', 'Color', acolors(i,:))
end
axis equal

figure(3) %TRAJECTORY
 clf(3)
 figure(3)
subplot(2,1,1)
hold on
title('DMP')
for i = 1:size(D.w, 2)
    plot(yalldmp(dmpcolor==(1+i),1), yalldmp(dmpcolor==(1+i),2), 'k.', 'Color', psicolors(i,:))
end
plot(P.desiredtraj(:,1),P.desiredtraj(:,2), '-k.');
axis image
subplot(2,1,2)
hold on
title('SMP')
for i = 1:size(S.w, 2)
    plot(yallrmp(smpcolor==(1+i),1), yallrmp(smpcolor==(1+i),2), 'k.', 'Color', acolors(i,:))
end
plot(P.desiredtraj(:,1),P.desiredtraj(:,2), '-k.');
axis image

figure(4) %FORCING FUNCTION
 clf(4)
 figure(4)
subplot(1,2,1)
hold on
title('DMP')
plot((1:P.Nt)*P.dt,D.f(:,1),'-k')
plot((1:P.Nt)*P.dt,D.f(:,2),'--k')
xlabel('time (s)')
ylabel('forcing function, {\it f(x)}')
% legend('x-component','y-component')
% axis image
subplot(1,2,2)
hold on
title('SMP')
plot((1:P.Nt)*P.dt,S.f(:,1),'-k')
plot((1:P.Nt)*P.dt,S.f(:,2),'--k')
xlabel('time (s)')
ylabel('forcing function, {\it f(x)}')
%legend('x-component','y-component')
% axis image

figure(5) %UNWEIGHTED KERNELS
 clf(5)
 figure(5)
subplot(2,1,1)
hold on
title('Unweighted basis functions (dmp)')
for i = 1:size(psiall,2)
    plot((1:P.Nt)*P.dt, psiall(:,i),'-c', 'Color', psicolors(i,:))
end
subplot(2,1,2)
hold on
title('Unweighted basis functions (smp)')
for i = 1:size(aall,2)
    plot((1:P.Nt)*P.dt, aall(:,i),'-c','Color',acolors(i,:))
end

figure(6) %SEPARATE KERNELS
 clf(6)
 figure(6)
subplot(2,1,1)
hold on
title('DMP Basis Functions')
plot((1:P.Nt)*P.dt, psiall(:,5),'-k','Linewidth', 1.5)
set(gca,'xtick',[])
axis square
subplot(2,1,2)
hold on
title('SMP Basis Functions')
plot((1:P.Nt)*P.dt,aall(:,5),'-k','Linewidth',1.5)
% legend('DMP','SMP')
set(gca,'xtick',[])
axis square

figure(7) %KERNEL COMPARISON
 clf(7)
 figure(7)
hold on
title('Basis Functions')
plot((1:P.Nt)*P.dt, psiall(:,5),'-k','Linewidth', 1.5)
plot((1:P.Nt)*P.dt,aall(:,5),'--k','Linewidth',1.5)
legend('DMP','RMP')
set(gca,'xtick',[])

end

