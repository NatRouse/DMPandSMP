close all
clear all
clc
tic
rng(1)

finalds = [];
for Nsaddles = 5:20
    OneforDMPtwoforRMP = 1; 
    
    Usepresents = 1; %for an 'if' statement later
    dt = .01;
    Nt = 1000;
    g = [1,0];
    y0 = .001 * [1,1];
    y0dot= 0 * [1,1];
    x0 = 1;
    
    desiredtraj = [(0:.05:.5)*0, .5*cos(pi/2:.1:pi*2.6);
        0:.05:.5, .5*sin(pi/2:.1:pi*2.6)]';
    g = [0,0];

    tau = 1; %time scaling factor, does nothing when tau=1
    alphay = 4;% damping and stiffness
    betay = alphay/4; % stiffnesses
         % if I decrease these numbers weird stuff happens
    alphax = .5; %how fast cannonical system decays
    Npsi = Nsaddles;
    % Npsi = 10;
    c = logspace(log10(1), log10(.01), Npsi); %center of the basis functions
        % these are the positions of DMP in X, and they do not change.
    sigma = logspace(log10(.3),log10(.002), Npsi); %"width" of basis function
        %if I make this narrower (smaller sigma) it looks like it can get smaller turn radius, but best solution has a kink
    %Cz = 0 no coupling
    
    Ktrials = 1000;
    
    %initializing trial function weights
    w = desiredtraj(round(linspace(1,length(desiredtraj),Npsi)) ,:)';
        %these are INDECES to the array "desiredtraj" defined earlier
    w1 = w;
    w = desiredtraj(round(linspace(1,length(desiredtraj),Npsi)) ,:)'*(alphay+betay);
    w2 = w;
    
    if Usepresents %Usepresents is a variable/boolean defined earlier
        %%%%%%% $$$ Here I am inserting linear algebra approach to weight finding
        % where we know f = tau*ydotdot - (alphay*(betay*(g-y) - ydot);
        desiredt =  linspace(0,10,length(desiredtraj))'; %desired length of motion, in time(sec?)
            % assumes constant velocity linspace(0,10,length(desiredtraj))
        desiredydot = [0,0;diff(desiredtraj)./(diff(desiredt)*[1,1])]; % desired speeds
        desiredydotdot = [0,0;diff(desiredydot)./(diff(desiredt)*[1,1])];
        fdesired = tau * desiredydotdot - (alphay * (betay*(g-desiredtraj) - desiredydot));
        
        if OneforDMPtwoforRMP == 1
            %'x' is physical position(?)
                %in schaal, 'x' is time scaled up or down by tau
            expected_x = x0*exp(-alphax/tau*(desiredt-desiredt(1)));
            expectedPsi = exp(ones(length(desiredt),1)*(-1./(2*sigma.^2)) .*(expected_x - c).^2); %n traj points by n basis
            
            % proveXexpworking % run this code to make a plot
            betaw =  (expected_x .* expectedPsi) \ (sum(expectedPsi')'.* fdesired);
            errors = (expected_x .* expectedPsi)*betaw - (sum(expectedPsi')'.* fdesired);
        else
            % Nsaddles = 10; %DUPLICATED IN EVAL
            %some variables that exist in "eval" function - consider
            %deleting
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
                %iterates by dt over the desired timeframe
                %canonical system??? don't touch this (yet)
                dW = sqrt(dt)*randn(1,Nsaddles);
                da = a .* (alpha - a*rho) *dt + ep.*dW; %wondering if I forgot tau
                a = max(min(a + da, 1), .0005);
                
                if t>=desiredt(index)
                    index = index+1;
                else
                    expected_SHCas(index, :) = a;
                    %SHC is the basis function? or something?
                end
            end
            
            
            size(expected_SHCas) %print the size of expectedSHCas to command line
            betaw =  ( expected_SHCas) \ (sum(expected_SHCas')'.* fdesired);
                %gain term in canonical system
        end
        
        
        % meanerror = sum(sum(errors.^2))
        w = betaw';
            %WEIGHTS OF BASIS FUNCTIONS
        %%%%%%% $$$ end insertion
    end
    
    
    %TWO main functions here, one to evaluate weights
    d = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1, OneforDMPtwoforRMP, w1,w2,Nsaddles);
        %d is the "cost"
    delta = 1;
    allds = [];
    for passes = 1:10
        bestw = w;
        bestd = d;
        for i = 1:Ktrials
            if rand(1)>.5
                % try random stuff
                w = 1-.5*rand(2, Npsi);
                d = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP, w1,w2,Nsaddles);
                if d<bestd
                    bestw = w;
                    bestd = d;
                    delta = 1;
                end
            else
                %try different random stuff based on size delta
                w = bestw + 3*delta*randn(2, Npsi);
                d = Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP, w1,w2,Nsaddles);
                if d<bestd
                    bestw = w;
                    bestd = d;
                end
                
            end
            
        end
        w = bestw;
        d = bestd;
        delta = delta*.99;
        Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1,OneforDMPtwoforRMP, w1,w2,Nsaddles);
        allds = [allds,bestd];
        
        
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
                    dtrial =  Eval(wtrial, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP,w1,w2,Nsaddles);
                    wnew(i,j) = bestw(i,j) + max(min(bestd*delt/(bestd - dtrial), delta*5), -delta*5);
                    dnew =  Eval(wnew, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 0,OneforDMPtwoforRMP,w1,w2,Nsaddles);
                    if dnew < bestd
                        bestw = wnew;
                        bestd = dnew;
                        chosenew = chosenew +1;
                    end
                    if dtrial < bestd
                        bestw = wtrial;
                        bestd = dtrial;
                        chosedel = chosedel +1;
                    end
                    
                end
            end
            %  chosenew
            %  chosedel
        end
        
        %at the end of one pass, save progress
        allds = [allds,bestd];
        w = bestw;
        d = bestd;
        %ShowTraj(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj);
        % if passes==1 || passes==5 || passes==10
        %     pause
        % end
    end
    Eval(w, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1,OneforDMPtwoforRMP,w1,w2,Nsaddles);
    
    % figure(4)
    % plot(allds)
    
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
    
    finalds = [finalds,allds(end)]
    figure(6)
    plot(finalds)
end

toc