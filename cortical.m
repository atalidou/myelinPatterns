%% This code computes the membrane potential along axons of cortical myelin
%% motifs. The percentage of myelin coverage has been set to 70%.
%
% The axons are constructed using the Poisson distribution. For the numerical
% scheme we used Finite Differences in space and backward Euler in time.
%%

clear;
close all;


%% ------------------------------------------------------------------------
% # of axons, time step, # of discrete points along each segment

numberOfAxons = 100;                        % number of axons
dt = 0.003;                                 % time step [us]
NX = 20;                                    % number of grid points along each segment

%% ------------------------------------------------------------------------
% axon setup

axonLength = 10000;                         % length of the myelinated axon (the reference axon) [um]
perMyel = 70;                               % percentage of myelin coverage (PMC) [OTHER CHOICES: 10, 30, 50, 70, 90]

%% -------------------------------------------------------------------------
% physiological parameters 

GL = 1*1e-8;                                % leak conductance density 
EL = -68.3;                                 % leak reversal potential [mV] 

Kth = 3.5;                                  % slope factor [mV] 
Ath = 520;                                  % (max(g_dep) = gL*Ath) 
Vth = -60.2;                                % spike threshold [mV]

Vsp = 10;                                   % spike-detecting threshold [mV]
Vre = EL + 20;                              % membrane potential for assuring repolarization [mV]

spT = 1.00;                                 % repolarization current time constant [ms]
spA = 90;                                   % (max(g_repolarize)= gL*spA) 

ee = exp(-dt/spT);                          % decay factor (alpha function) 
aa = spA*exp(1)/spT;                        % amplitude factor (alpha function)

a = 2;                                      % axon radius [um]
g = 0.6;                                    % g-ratio
r_int = 100;                                % axial resistivity [Mohm*um]
Rc = r_int/(pi*a^2);                        % cable resistance
Cc = 3.5*1e-9;                              % cable membrane capacitance [F/m]
eps = 1.4*1e-10;                            % permitivity 
Cmy = (2*pi*eps)/log(1/g);                  % capacitance of myelin sheath

%% ------------------------------------------------------------------------
% preallocate vectors to save conduction delays, velocities and conduction
% success/failures

tau = zeros(numberOfAxons, 1);              % conduction delays of all trials 
c = zeros(numberOfAxons, 1);                % conduction velocities of all trials
FailSuccess = NaN(numberOfAxons, 1);        % vector of conduction failures(0) and successes(1) 


%% ------------------------------------------------------------------------
% start the iteration for the number of trials

for trial = 1:numberOfAxons

    %% --------------------------------------------------------------------
    % call the function of the axon

    [numMyel, l0, lm, startAxon, endAxon, myel, node, interface, myelExtra, nodeExtra, interfaceAll] = axon_cortical10mm(axonLength, perMyel);

    %% --------------------------------------------------------------------
    % reference axon 

    numNode = numMyel;                      % number of ummyelinated segments
    m = numMyel + numNode + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAll(1)-l0)/NX:interfaceAll(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAll)
        xgrid(:, i) = interfaceAll(i-1):(interfaceAll(i)-interfaceAll(i-1))/NX:interfaceAll(i);
    end

    xgrid(:, m) = interfaceAll(end):(lm-interfaceAll(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)

    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1,m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end

    %% --------------------------------------------------------------------
    % FD matrix

    N = (NX+1)*m;                   % total number of unknowns

    A = zeros(N, N);                % full matrix

    % left boundary
    A(1,1) = 1 + (dt/(dx(1,1)^2))*kappa(1);
    A(1,2) = -(dt/(dx(1,1)^2))*kappa(1);

    % right boundary
    A(N,N)   = 1 + (dt/(dx(NX,m)^2))*kappa(m);
    A(N,N-1) = -(dt/(dx(NX,m)^2))*kappa(m);

    % mid compartmental points within UNmyelinated segments
    for i = 1:2:m
        for j = 2:NX
            r = (i-1)*(NX+1) + j;

            A(r,r-1) = -(dt/(dx(j,i)^2)) * kappa(i);
            A(r,r)   = 1 + 2*(dt/(dx(j,i)^2)) * kappa(i);
            A(r,r+1) = -(dt/(dx(j,i)^2)) * kappa(i);
        end
    end

    % mid compartmental points within Myelinated segments
    for i = 2:2:m
        for j = 2:NX
            r = (i-1)*(NX+1) + j;

            A(r,r-1) = -(dt/(dx(j,i)^2)) * kappa(i);
            A(r,r)   = 1 + 2*(dt/(dx(j,i)^2)) * kappa(i);
            A(r,r+1) = -(dt/(dx(j,i)^2)) * kappa(i);
        end
    end

    % compartmental point at the left interface of UNmyelinated segments
    for i = 3:2:m
        j = 1;
        r = (i-1)*(NX+1) + j;

        A(r,r-1) = -(dt/(dx(j,i)*dx(NX,i-1)))*kappa(i-1);  
        A(r,r)   = 1 + (dt/(dx(j,i)^2))*kappa(i) + (dt/(dx(j,i)*dx(NX,i-1)))*kappa(i-1);
        A(r,r+1) = -(dt/(dx(j,i)^2))*kappa(i);
    end

    % compartmental point at the left interface of Myelinated segments
    for i = 2:2:m
        j = 1;
        r = (i-1)*(NX+1) + j;

        A(r,r-1) = -(dt/(dx(j,i)*dx(NX,i-1)))*kappa(i-1);  
        A(r,r)   = 1 + (dt/(dx(j,i)^2))*kappa(i) + (dt/(dx(j,i)*dx(NX,i-1)))*kappa(i-1);
        A(r,r+1) = -(dt/(dx(j,i)^2))*kappa(i);
    end

    % compartmental point at the right interface of UNmyelinated segments
    for i = 1:2:m
        j = NX+1;
        r = (i-1)*(NX+1) + j;

        A(r,r-1) = -(dt/(dx(NX,i)^2))*kappa(i);
        A(r,r)   = 1 + (dt/(dx(NX,i)^2))*kappa(i) + (dt/(dx(NX,i)*dx(1,i+1)))*kappa(i+1);
        A(r,r+1) = -(dt/(dx(NX,i)*dx(1,i+1)))*kappa(i+1); 
    end

    % compartmental point at the right interface of Myelinated segments
    for i = 2:2:m-1
        j = NX+1;
        r = (i-1)*(NX+1) + j;

        A(r,r-1) = -(dt/(dx(NX,i)^2))*kappa(i);
        A(r,r)   = 1 + (dt/(dx(NX,i)^2))*kappa(i) + (dt/(dx(NX,i)*dx(1,i+1)))*kappa(i+1);
        A(r,r+1) = -(dt/(dx(NX,i)*dx(1,i+1)))*kappa(i+1); 
    end

    % check invertibility
    
    if rank(A) < N
        fprintf('axon # is %d | matrix is NOT invertible \n', trial)
        break
    else
        fprintf('axon # is %d | matrix is invertible \n', trial)
    end

    A = sparse(A);                  % sparse matrix

    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0 = @(x) ones(size(x));
    u = EL*u0(xgrid);    

    b = zeros(N,1);                 % vector to store u

    usoln = zeros(N,length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid,N,1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    % initialize currents     
    iL = zeros(N,1);                % leak current 
    iD = zeros(N,1);                % depolarizing current
    iR = zeros(N,1);                % repolarizing current  
    iI = zeros(N,1);                % injected current
    iTotal = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR = zeros(N,1);                % main
    yR = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady = ones(N,1);            % ready to spike 
    spStart = zeros(N,1);           % start spiking

    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b(1) = u(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b(r) = u(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b(r) = u(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b(r) = u(j,i);
        end

        b(end) = u(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI(3*(NX-4)) = 0;
        end

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL(r) = GL*(EL - b(r));
                iD(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b(r) - Vth)/Kth)));
                iR(r) = GL*xR(r)*(EL - b(r));

                iTotal(r) = (dt/Cc)*(iL(r) + iD(r) + iR(r) + iI(r));
            end
        end

        % calculate the numerical solution
        soln = A\(b + iTotal);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln(r) >= Vsp) 
                        spStart(r) = 1; % started spiking
                        spReady(r) = 0; % not ready for spiking
                    else 
                        spStart(r) = 0; % not started spiking
                        spReady(r) = 1; % ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln(r) <= Vre)
                        spStart(r) = 0; % not starting spiking 
                        spReady(r) = 1; % ready for next spike
                    else
                        spStart(r) = 0; % not starting spiking 
                        spReady(r) = 0; % not ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR = ee*yR + aa*spStart;
        xR = ee*xR + dt*yR;

        u = reshape(soln, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln(:,k) = soln;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        

    end % end of the while loop

    u = usoln;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex = zeros(1, length(interface)+2);    % vector of the indices where the AP crosses the threshold
    Tspike = zeros(1, length(interface)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel = zeros(1, numMyel);               % vector of the times along myelinated segments
    condSpeedMyel = zeros(1, numMyel);              % vector of the speeds along myelinated axons

    condTimeNode = zeros(1, numNode);               % vector of the times along exposed segments
    condSpeedNode = zeros(1, numNode);              % vector of the speeds along exposed axons



    TspikeIndexF = zeros(1, length(interface)+2);    
    TspikeF = zeros(1, length(interface)+2);     

    condTimeMyelF = zeros(1, numMyel);               
    condSpeedMyelF = zeros(1, numMyel);              

    condTimeNodeF = zeros(1, numNode);               
    condSpeedNodeF = zeros(1, numNode);              



    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u((length(interfaceAll)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess(trial) = -1;     % the "-1" indicates instability

    elseif isempty(find(u((length(interfaceAll)-2*length(myelExtra))*(NX+1),:)>=Vth & (u((length(interfaceAll)-2*length(myelExtra))*(NX+1),end))<Vth, 1)) == 0 
        
        fprintf('axon # is %d --> AP success \n', trial)

        FailSuccess(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interface) + 2
            TspikeIndex(i) = find(u((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike(i) = tspan(TspikeIndex(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyel
            condTimeMyel(ii) = Tspike(2*ii+1) - Tspike(2*ii);   % [ms]
            condSpeedMyel(ii) = myel(ii)/(condTimeMyel(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNode
            condTimeNode(ii) = Tspike(2*ii) - Tspike(2*ii-1);   % [ms]
            condSpeedNode(ii) = node(ii)/(condTimeNode(ii)*1000);  % [m/s]
        end

        tau(trial) = sum(condTimeMyel) + sum(condTimeNode);     % conduction delay of each trial [ms]
        c(trial) = axonLength/(tau(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay = %f ms', tau(trial))
        fprintf('\n conduction velocity = %f m/s \n', c(trial))

    else

        fprintf('axon # is %d --> AP failure \n', trial)

        FailSuccess(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interface) + 2

            if isempty(find(u((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF(i) = find(u((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF(i) = tspan(TspikeIndexF(i));

            else

                TspikeIndexF(i) = NaN;
                TspikeF(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyel
            condTimeMyelF(ii) = TspikeF(2*ii+1) - TspikeF(2*ii);   % [ms]
            condSpeedMyelF(ii) = myel(ii)/(condTimeMyelF(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNode
            condTimeNodeF(ii) = TspikeF(2*ii) - TspikeF(2*ii-1);   % [ms]
            condSpeedNodeF(ii) = node(ii)/(condTimeNodeF(ii)*1000);  % [m/s]
        end
     
    end

    nS = length(find(FailSuccess == 1));
    nF = length(find(FailSuccess == 0));

    if nS+nF == 100
        break;
    end

end % end of the for loop

tauSuccess = tau(tau ~= 0);         % save the conduction delays of only successful trials
cSuccess = c(c ~= 0);               % save the conduction velocities of only successful trials















%% ------------------------------------------------------------------------
%% define the function for the axon 

function [numMyel, l0, lm, startAxon, endAxon, myel, node, interface, myelExtra, nodeExtra, interfaceAll] = axon_cortical10mm(axonLength, perMyel)

% numMyel: number of myelinated segments [scalar]
% l0: startimg point of the entire axon [scalar]
% lm: ending point of the entire axon [scalar]
% startAxon: starting point of the reference axon [scalar]
% endAxon: ending point of the reference axon [scalar]
% myel: lengths of myelinated segments along the reference axon [vector]
% node: lengths of unmyelinated (exposed) segments along the reference axon
% [vector]
% interface: interfaces between myelinated/unmyelinated segments of the
% reference axon [vector]
% myelExtra: lengths of myelinated segments along the extra space [vector]
% nodeExtra: lengths of unmyelinated segments along the extra space [vector]
% interfaceAll: interfaces between myelinated/unmyelinated segments of the
% entire axon [vector]

%% ------------------------------------------------------------------------
% Define the total myelinated/unmyelinated length, and the number of 
% myelinated/unmyelinated segments of the reference axon

myelAxon = (perMyel/100)*axonLength;    % total covered (by myelin) length of the reference axon [um]
nodeAxon = axonLength - myelAxon;       % total uncovered length of the reference axon [um]


numAxon = 0;
alliterations = 1000;

numberOfMyel = zeros(1, alliterations);
nummerOfNode = zeros(1, alliterations);

for k = 1:alliterations

    numberOfMyel(k) = randi([1, 500],1);            % number of myelinated segments along the reference axon (sample a uniformly distributed pseudorandom integer)     
    nummerOfNode(k) = numberOfMyel(k);              % number of unmyelinated segments along the reference axon

    %% ------------------------------------------------------------------------
    % MYELINATED SEGMENTS
    
    z = rand(myelAxon,1);
    myel = (accumarray(ceil(z*numberOfMyel(k)),1,[numberOfMyel(k),1]))';
    
    %% ------------------------------------------------------------------------
    % UNMYELINATED SEGMENTS
    
    w = rand(nodeAxon,1);
    node = (accumarray(ceil(w*nummerOfNode(k)),1,[nummerOfNode(k),1]))';
    
    %%
    
    if (min(myel)>=1 && min(node)>=1 && mean(myel)>=30 && mean(myel)<=100) % if these conditions are satisfied save the axon
        
        numMyel = numberOfMyel(k);
        
        %% ------------------------------------------------------------------------
        % Create a vector with the first 7 myelin segments repeated 3 times, and a 
        % vector with the first 7 nodes repeated 3 times.
        
        myelExtra = repmat(myel(1:7),1,3);
        nodeExtra = repmat(node(1:7),1,3);
        
        %% --------------------------------------------------------------------
        % Define the parameters of the entire axon
        
        l0 = 0;                                         % starting point of the entire axon (including the extra space)
        startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
        endAxon = startAxon + axonLength;               % ending point of the reference axon
        lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 
        
        %% ------------------------------------------------------------------------
        % Create the reference axon
            
        axon = zeros(1, numberOfMyel(k) + nummerOfNode(k)); 
        axon(1:2:end) = node;                           % odd segments are unmyelinated
        axon(2:2:end) = myel;                           % even segments are myelinated
        
        %% ------------------------------------------------------------------------
        % interfaces along the reference axon: contact points between myelinated 
        % and unmyelinated segments
            
        interface = zeros(1, length(axon) - 1);
        
        interface(1) = startAxon + axon(1);             % first interface
        
        % mid and last interfaces
        for j = 2:length(interface)
            interface(j) = interface(j-1) + axon(j);
        end
            
        %% ------------------------------------------------------------------------
        % Create the extra space at the beginning and end of the reference axon
        
        axonExtra = zeros(1, length(myelExtra) + length(nodeExtra)); 
        
        axonExtra(1:2:end) = nodeExtra;               % odd segments are unmyelinated
        axonExtra(2:2:end) = myelExtra;               % even segments are myelinated
            
        %% ------------------------------------------------------------------------
        % interfaces of the EXTRA START axon
            
        interfaceStart = zeros(1, length(axonExtra) - 1);
        
        interfaceStart(1) = l0 + axonExtra(1);        % first interface of the EXTRA START axon
        
        % mid and last interfaces of the EXTRA START axon
        for j = 2:length(interfaceStart)
            interfaceStart(j) = interfaceStart(j-1) + axonExtra(j);
        end
            
        %% ------------------------------------------------------------------------
        % interfaces of the EXTRA END axon
            
        interfaceEnd = zeros(1, length(myelExtra)*2 - 1);
        
        interfaceEnd(1) = endAxon + axonExtra(1);      % first interface of the EXTRA END axon
        
        % mid and last interfaces of the EXTRA END axon
        for j = 2:length(interfaceEnd)
            interfaceEnd(j) = interfaceEnd(j-1) + axonExtra(j);
        end
            
        %% ------------------------------------------------------------------------
        % merge the axon pieces all the interfaces
        
        interfaceAll = [interfaceStart startAxon interface endAxon interfaceEnd];

        %%

        numAxon = numAxon + 1;

        if numAxon == 1
            break;
        end
            
    else
        continue
    end


end

end
