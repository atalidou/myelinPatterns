%% The code computes the membrane potential along healthy and demyelinated
%% axons of callosal myelin motifs.
%
% Initially, the percentage of myelin coverage is set to 90%. Demyelination
% is performed by removing a certain number of myelin sheaths. Once a
% myelin sheath is retracted from the axon (demyelinated segment) it
% becomes an excitable segment. The demyelinated segments together with the
% adjacent exposed segments form a longer exposed segment with a new spike
% threshold value VthD.
%
% A healthy axon undergoes demyelination six times (D1, ..., D6) where D6
% is the case of severe myelin damage.
%
% The outcome is the conduction delay and corresponding conduction
% velocity. If the action potential (AP) fails to conduct, no conduction
% delays/velocities are computed.



clear;
close all;


%%
%% ------------------------------------------------------------------------
% # of axons, time step, # of discrete points along each segment

numberOfAxons = 500;                          % number of axons
dt = 0.003;                                 % time step [us]
NX = 20;                                    % number of grid points along each segment

damage = 0.05;


%% ------------------------------------------------------------------------
% axon setup

axonLength = 10000;                         % length of the myelinated axon (the reference axon) [um]
perMyel = 90;                               % percentage of myelin coverage (PMC)

%% -------------------------------------------------------------------------
% physiological parameters 

GL = 1*1e-8;                                % leak conductance density 
EL = -68.3;                                 % leak reversal potential [mV] 

Kth = 3.5;                                  % slope factor [mV] 
Ath = 520;                                  % (max(g_dep) = gL*Ath) 
Vth = -60.2;                                % spike threshold [mV]

VthD = -55;                                 % spike threshold at demyelinated segments [mV]

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

%% at "healthy" axons
tau = zeros(numberOfAxons, 1);              % conduction delays of all trials 
c = zeros(numberOfAxons, 1);                % conduction velocities of all trials
FailSuccess = NaN(numberOfAxons, 1);      % vector of conduction failures(0) and successes(1) 

%% at "D1_demyelinated" axons
tau_D1 = zeros(numberOfAxons, 1);             
c_D1 = zeros(numberOfAxons, 1);               
FailSuccess_D1 = NaN(numberOfAxons, 1);

%% at "D2_demyelinated" axons
tau_D2 = zeros(numberOfAxons, 1);             
c_D2 = zeros(numberOfAxons, 1);               
FailSuccess_D2 = NaN(numberOfAxons, 1);

%% at "D3_demyelinated" axons
tau_D3 = zeros(numberOfAxons, 1);             
c_D3 = zeros(numberOfAxons, 1);               
FailSuccess_D3 = NaN(numberOfAxons, 1);

%% at "D4_demyelinated" axons
tau_D4 = zeros(numberOfAxons, 1);             
c_D4 = zeros(numberOfAxons, 1);               
FailSuccess_D4 = NaN(numberOfAxons, 1);

%% at "D5_demyelinated" axons
tau_D5 = zeros(numberOfAxons, 1);             
c_D5 = zeros(numberOfAxons, 1);               
FailSuccess_D5 = NaN(numberOfAxons, 1);

%% at "D6_demyelinated" axons
tau_D6 = zeros(numberOfAxons, 1);             
c_D6 = zeros(numberOfAxons, 1);               
FailSuccess_D6 = NaN(numberOfAxons, 1);



%% ------------------------------------------------------------------------
% start the iteration for the number of trials

for trial = 1:numberOfAxons

    %% --------------------------------------------------------------------
    % Start with the "healthy" axon. Call the function to create the axon.

    [numMyel, l0, lm, startAxon, endAxon, myel, node, interface, myelExtra, nodeExtra, interfaceAll] = axon_callosal10mm(axonLength, perMyel);

    %% --------------------------------------------------------------------
    % reference axon 
    
    numNode = numMyel                      % number of ummyelinated segments
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
    % FD matrix of "healthy" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

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
                        spStart(r) = 1; % now started spiking
                        spReady(r) = 0; % thus not ready for spiking
                    else 
                        spStart(r) = 0; % not started spiking yet
                        spReady(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln(r) <= Vre)
                        spStart(r) = 0; % not starting spiking, of course 
                        spReady(r) = 1; % and now ready for next spike
                    else
                        spStart(r) = 0; % not starting spiking, of course 
                        spReady(r) = 0; % not yet ready for next spike
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

    end

    pointer = NaN(1,length(interface) + 2);

    for i=1:length(interface) + 2

        if isempty(find(u((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end


    if nnz(~pointer) == 0

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



    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D1_demyelinated" axon

    [ix_affectedD1, ix_unaffectedD1, node_demyelD1, ix_demyelD1, numMyelD1, numNodeD1, myelD1, nodeD1, interfaceD1, interfaceAllD1] = axonDemyel_D1_callosal(axonLength, numMyel, numNode, myel, node, damage, myelExtra, nodeExtra);

    
    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD1 + numNodeD1 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D1 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD1 - 1;
    ix_unaffectedAxon_D1 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD1 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD1(1)-l0)/NX:interfaceAllD1(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD1)
        xgrid(:, i) = interfaceAllD1(i-1):(interfaceAllD1(i)-interfaceAllD1(i-1))/NX:interfaceAllD1(i);
    end

    xgrid(:, m) = interfaceAllD1(end):(lm-interfaceAllD1(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D1" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D1 = @(x) ones(size(x));
    u_D1 = EL*u0_D1(xgrid);    

    b_D1 = zeros(N, 1);                 % vector to store u

    usoln_D1 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D1 = zeros(N,1);                % leak current 
    iD_D1 = zeros(N,1);                % depolarizing current
    iR_D1 = zeros(N,1);                % repolarizing current  
    iI_D1 = zeros(N,1);                % injected current
    iTotal_D1 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D1 = zeros(N,1);                % main
    yR_D1 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D1 = ones(N,1);            % ready to spike 
    spStart_D1 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D1(1) = u_D1(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D1(r) = u_D1(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D1(r) = u_D1(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D1(r) = u_D1(j,i);
        end

        b_D1(end) = u_D1(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D1) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D1(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D1(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D1(r) = GL*(EL - b_D1(r));
                iD_D1(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D1(r) - Vth)/Kth)));
                iR_D1(r) = GL*xR_D1(r)*(EL - b_D1(r));

                iTotal_D1(r) = (dt/Cc)*(iL_D1(r) + iD_D1(r) + iR_D1(r) + iI_D1(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D1
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D1(r) = GL*(EL - b_D1(r));
                iD_D1(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D1(r) - VthD)/Kth)));
                iR_D1(r) = GL*xR_D1(r)*(EL - b_D1(r));

                iTotal_D1(r) = (dt/Cc)*(iL_D1(r) + iD_D1(r) + iR_D1(r) + iI_D1(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D1 = A\(b_D1 + iTotal_D1);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D1(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D1(r) >= Vsp) 
                        spStart_D1(r) = 1; % now started spiking
                        spReady_D1(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D1(r) = 0; % not started spiking yet
                        spReady_D1(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D1(r) <= Vre)
                        spStart_D1(r) = 0; % not starting spiking, of course 
                        spReady_D1(r) = 1; % and now ready for next spike
                    else
                        spStart_D1(r) = 0; % not starting spiking, of course 
                        spReady_D1(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D1 = ee*yR_D1 + aa*spStart_D1;
        xR_D1 = ee*xR_D1 + dt*yR_D1;

        u_D1 = reshape(soln_D1, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D1(:,k) = soln_D1;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D1 = usoln_D1;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D1 = zeros(1, length(interfaceD1)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D1 = zeros(1, length(interfaceD1)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D1 = zeros(1, numMyelD1);               % vector of the times along myelinated segments
    condSpeedMyel_D1 = zeros(1, numMyelD1);              % vector of the speeds along myelinated axons

    condTimeNode_D1 = zeros(1, numNodeD1);               % vector of the times along exposed segments
    condSpeedNode_D1 = zeros(1, numNodeD1);              % vector of the speeds along exposed axons



    TspikeIndexF_D1 = zeros(1, length(interfaceD1)+2);    
    TspikeF_D1 = zeros(1, length(interfaceD1)+2);     

    condTimeMyelF_D1 = zeros(1, numMyelD1);               
    condSpeedMyelF_D1 = zeros(1, numMyelD1);              

    condTimeNodeF_D1 = zeros(1, numNodeD1);               
    condSpeedNodeF_D1 = zeros(1, numNodeD1);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D1((length(interfaceAllD1)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D1(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD1) + 2);

    for i=1:length(interfaceD1) + 2

        if isempty(find(u_D1((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D1 demyelination \n', trial)

        FailSuccess_D1(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD1) + 2
            TspikeIndex_D1(i) = find(u_D1((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D1(i) = tspan(TspikeIndex_D1(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD1
            condTimeMyel_D1(ii) = Tspike_D1(2*ii+1) - Tspike_D1(2*ii);   % [ms]
            condSpeedMyel_D1(ii) = myelD1(ii)/(condTimeMyel_D1(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD1
            condTimeNode_D1(ii) = Tspike_D1(2*ii) - Tspike_D1(2*ii-1);   % [ms]
            condSpeedNode_D1(ii) = nodeD1(ii)/(condTimeNode_D1(ii)*1000);  % [m/s]
        end

        tau_D1(trial) = sum(condTimeMyel_D1) + sum(condTimeNode_D1);     % conduction delay of each trial [ms]
        c_D1(trial) = axonLength/(tau_D1(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D1 demyelinated axon = %f ms', tau_D1(trial))
        fprintf('\n conduction velocity of D1 demyelinated axon = %f m/s \n', c_D1(trial))

    else

        fprintf('axon # is %d --> AP failure after D1 demyelination \n', trial)

        FailSuccess_D1(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD1) + 2

            if isempty(find(u_D1((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D1(i) = find(u_D1((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D1(i) = tspan(TspikeIndexF_D1(i));

            else

                TspikeIndexF_D1(i) = NaN;
                TspikeF_D1(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD1
            condTimeMyelF_D1(ii) = TspikeF_D1(2*ii+1) - TspikeF_D1(2*ii);   % [ms]
            condSpeedMyelF_D1(ii) = myelD1(ii)/(condTimeMyelF_D1(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD1
            condTimeNodeF_D1(ii) = TspikeF_D1(2*ii) - TspikeF_D1(2*ii-1);   % [ms]
            condSpeedNodeF_D1(ii) = nodeD1(ii)/(condTimeNodeF_D1(ii)*1000);  % [m/s]
        end

    end
    
    
        
    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D2_demyelinated" axon

    [ix_affectedD2, ix_unaffectedD2, node_demyelD2, ix_demyelD2, numMyelD2, numNodeD2, myelD2, nodeD2, interfaceD2, interfaceAllD2] = axonDemyel_D2_callosal(axonLength, numMyelD1, numNodeD1, myelD1, nodeD1, damage, myelExtra, nodeExtra);
    
    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD2 + numNodeD2 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D2 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD2 - 1;
    ix_unaffectedAxon_D2 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD2 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD2(1)-l0)/NX:interfaceAllD2(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD2)
        xgrid(:, i) = interfaceAllD2(i-1):(interfaceAllD2(i)-interfaceAllD2(i-1))/NX:interfaceAllD2(i);
    end

    xgrid(:, m) = interfaceAllD2(end):(lm-interfaceAllD2(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D2" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D2 = @(x) ones(size(x));
    u_D2 = EL*u0_D2(xgrid);    

    b_D2 = zeros(N, 1);                 % vector to store u

    usoln_D2 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D2 = zeros(N,1);                % leak current 
    iD_D2 = zeros(N,1);                % depolarizing current
    iR_D2 = zeros(N,1);                % repolarizing current  
    iI_D2 = zeros(N,1);                % injected current
    iTotal_D2 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D2 = zeros(N,1);                % main
    yR_D2 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D2 = ones(N,1);            % ready to spike 
    spStart_D2 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D2(1) = u_D2(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D2(r) = u_D2(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D2(r) = u_D2(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D2(r) = u_D2(j,i);
        end

        b_D2(end) = u_D2(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D2) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D2(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D2(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D2(r) = GL*(EL - b_D2(r));
                iD_D2(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D2(r) - Vth)/Kth)));
                iR_D2(r) = GL*xR_D2(r)*(EL - b_D2(r));

                iTotal_D2(r) = (dt/Cc)*(iL_D2(r) + iD_D2(r) + iR_D2(r) + iI_D2(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D2
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D2(r) = GL*(EL - b_D2(r));
                iD_D2(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D2(r) - VthD)/Kth)));
                iR_D2(r) = GL*xR_D2(r)*(EL - b_D2(r));

                iTotal_D2(r) = (dt/Cc)*(iL_D2(r) + iD_D2(r) + iR_D2(r) + iI_D2(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D2 = A\(b_D2 + iTotal_D2);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D2(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D2(r) >= Vsp) 
                        spStart_D2(r) = 1; % now started spiking
                        spReady_D2(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D2(r) = 0; % not started spiking yet
                        spReady_D2(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D2(r) <= Vre)
                        spStart_D2(r) = 0; % not starting spiking, of course 
                        spReady_D2(r) = 1; % and now ready for next spike
                    else
                        spStart_D2(r) = 0; % not starting spiking, of course 
                        spReady_D2(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D2 = ee*yR_D2 + aa*spStart_D2;
        xR_D2 = ee*xR_D2 + dt*yR_D2;

        u_D2 = reshape(soln_D2, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D2(:,k) = soln_D2;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D2 = usoln_D2;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D2 = zeros(1, length(interfaceD2)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D2 = zeros(1, length(interfaceD2)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D2 = zeros(1, numMyelD2);               % vector of the times along myelinated segments
    condSpeedMyel_D2 = zeros(1, numMyelD2);              % vector of the speeds along myelinated axons

    condTimeNode_D2 = zeros(1, numNodeD2);               % vector of the times along exposed segments
    condSpeedNode_D2 = zeros(1, numNodeD2);              % vector of the speeds along exposed axons



    TspikeIndexF_D2 = zeros(1, length(interfaceD2)+2);    
    TspikeF_D2 = zeros(1, length(interfaceD2)+2);     

    condTimeMyelF_D2 = zeros(1, numMyelD2);               
    condSpeedMyelF_D2 = zeros(1, numMyelD2);              

    condTimeNodeF_D2 = zeros(1, numNodeD2);               
    condSpeedNodeF_D2 = zeros(1, numNodeD2);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D2((length(interfaceAllD2)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D2(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD2) + 2);

    for i=1:length(interfaceD2) + 2

        if isempty(find(u_D2((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D2 demyelination \n', trial)

        FailSuccess_D2(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD2) + 2
            TspikeIndex_D2(i) = find(u_D2((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D2(i) = tspan(TspikeIndex_D2(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD2
            condTimeMyel_D2(ii) = Tspike_D2(2*ii+1) - Tspike_D2(2*ii);   % [ms]
            condSpeedMyel_D2(ii) = myelD2(ii)/(condTimeMyel_D2(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD2
            condTimeNode_D2(ii) = Tspike_D2(2*ii) - Tspike_D2(2*ii-1);   % [ms]
            condSpeedNode_D2(ii) = nodeD2(ii)/(condTimeNode_D2(ii)*1000);  % [m/s]
        end

        tau_D2(trial) = sum(condTimeMyel_D2) + sum(condTimeNode_D2);     % conduction delay of each trial [ms]
        c_D2(trial) = axonLength/(tau_D2(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D2 demyelinated axon = %f ms', tau_D2(trial))
        fprintf('\n conduction velocity of D2 demyelinated axon = %f m/s \n', c_D2(trial))

    else

        fprintf('axon # is %d --> AP failure after D2 demyelination \n', trial)

        FailSuccess_D2(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD2) + 2

            if isempty(find(u_D2((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D2(i) = find(u_D2((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D2(i) = tspan(TspikeIndexF_D2(i));

            else

                TspikeIndexF_D2(i) = NaN;
                TspikeF_D2(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD2
            condTimeMyelF_D2(ii) = TspikeF_D2(2*ii+1) - TspikeF_D2(2*ii);   % [ms]
            condSpeedMyelF_D2(ii) = myelD2(ii)/(condTimeMyelF_D2(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD2
            condTimeNodeF_D2(ii) = TspikeF_D2(2*ii) - TspikeF_D2(2*ii-1);   % [ms]
            condSpeedNodeF_D2(ii) = nodeD2(ii)/(condTimeNodeF_D2(ii)*1000);  % [m/s]
        end

    end
    
   
     
    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D3_demyelinated" axon

    [ix_affectedD3, ix_unaffectedD3, node_demyelD3, ix_demyelD3, numMyelD3, numNodeD3, myelD3, nodeD3, interfaceD3, interfaceAllD3] = axonDemyel_D3_callosal(axonLength, numMyelD2, numNodeD2, myelD2, nodeD2, damage, myelExtra, nodeExtra);

    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD3 + numNodeD3 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D3 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD3 - 1;
    ix_unaffectedAxon_D3 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD3 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD3(1)-l0)/NX:interfaceAllD3(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD3)
        xgrid(:, i) = interfaceAllD3(i-1):(interfaceAllD3(i)-interfaceAllD3(i-1))/NX:interfaceAllD3(i);
    end

    xgrid(:, m) = interfaceAllD3(end):(lm-interfaceAllD3(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D3" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D3 = @(x) ones(size(x));
    u_D3 = EL*u0_D3(xgrid);    

    b_D3 = zeros(N, 1);                 % vector to store u

    usoln_D3 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D3 = zeros(N,1);                % leak current 
    iD_D3 = zeros(N,1);                % depolarizing current
    iR_D3 = zeros(N,1);                % repolarizing current  
    iI_D3 = zeros(N,1);                % injected current
    iTotal_D3 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D3 = zeros(N,1);                % main
    yR_D3 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D3 = ones(N,1);            % ready to spike 
    spStart_D3 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D3(1) = u_D3(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D3(r) = u_D3(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D3(r) = u_D3(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D3(r) = u_D3(j,i);
        end

        b_D3(end) = u_D3(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D3) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D3(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D3(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D3(r) = GL*(EL - b_D3(r));
                iD_D3(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D3(r) - Vth)/Kth)));
                iR_D3(r) = GL*xR_D3(r)*(EL - b_D3(r));

                iTotal_D3(r) = (dt/Cc)*(iL_D3(r) + iD_D3(r) + iR_D3(r) + iI_D3(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D3
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D3(r) = GL*(EL - b_D3(r));
                iD_D3(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D3(r) - VthD)/Kth)));
                iR_D3(r) = GL*xR_D3(r)*(EL - b_D3(r));

                iTotal_D3(r) = (dt/Cc)*(iL_D3(r) + iD_D3(r) + iR_D3(r) + iI_D3(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D3 = A\(b_D3 + iTotal_D3);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D3(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D3(r) >= Vsp) 
                        spStart_D3(r) = 1; % now started spiking
                        spReady_D3(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D3(r) = 0; % not started spiking yet
                        spReady_D3(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D3(r) <= Vre)
                        spStart_D3(r) = 0; % not starting spiking, of course 
                        spReady_D3(r) = 1; % and now ready for next spike
                    else
                        spStart_D3(r) = 0; % not starting spiking, of course 
                        spReady_D3(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D3 = ee*yR_D3 + aa*spStart_D3;
        xR_D3 = ee*xR_D3 + dt*yR_D3;

        u_D3 = reshape(soln_D3, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D3(:,k) = soln_D3;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D3 = usoln_D3;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D3 = zeros(1, length(interfaceD3)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D3 = zeros(1, length(interfaceD3)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D3 = zeros(1, numMyelD3);               % vector of the times along myelinated segments
    condSpeedMyel_D3 = zeros(1, numMyelD3);              % vector of the speeds along myelinated axons

    condTimeNode_D3 = zeros(1, numNodeD3);               % vector of the times along exposed segments
    condSpeedNode_D3 = zeros(1, numNodeD3);              % vector of the speeds along exposed axons



    TspikeIndexF_D3 = zeros(1, length(interfaceD3)+2);    
    TspikeF_D3 = zeros(1, length(interfaceD3)+2);     

    condTimeMyelF_D3 = zeros(1, numMyelD3);               
    condSpeedMyelF_D3 = zeros(1, numMyelD3);              

    condTimeNodeF_D3 = zeros(1, numNodeD3);               
    condSpeedNodeF_D3 = zeros(1, numNodeD3);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D3((length(interfaceAllD3)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D3(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD3) + 2);

    for i=1:length(interfaceD3) + 2

        if isempty(find(u_D3((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D3 demyelination \n', trial)

        FailSuccess_D3(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD3) + 2
            TspikeIndex_D3(i) = find(u_D3((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D3(i) = tspan(TspikeIndex_D3(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD3
            condTimeMyel_D3(ii) = Tspike_D3(2*ii+1) - Tspike_D3(2*ii);   % [ms]
            condSpeedMyel_D3(ii) = myelD3(ii)/(condTimeMyel_D3(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD3
            condTimeNode_D3(ii) = Tspike_D3(2*ii) - Tspike_D3(2*ii-1);   % [ms]
            condSpeedNode_D3(ii) = nodeD3(ii)/(condTimeNode_D3(ii)*1000);  % [m/s]
        end

        tau_D3(trial) = sum(condTimeMyel_D3) + sum(condTimeNode_D3);     % conduction delay of each trial [ms]
        c_D3(trial) = axonLength/(tau_D3(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D3 demyelinated axon = %f ms', tau_D3(trial))
        fprintf('\n conduction velocity of D3 demyelinated axon = %f m/s \n', c_D3(trial))


    else

        fprintf('axon # is %d --> AP failure after D3 demyelination \n', trial)

        FailSuccess_D3(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD3) + 2

            if isempty(find(u_D3((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D3(i) = find(u_D3((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D3(i) = tspan(TspikeIndexF_D3(i));

            else

                TspikeIndexF_D3(i) = NaN;
                TspikeF_D3(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD3
            condTimeMyelF_D3(ii) = TspikeF_D3(2*ii+1) - TspikeF_D3(2*ii);   % [ms]
            condSpeedMyelF_D3(ii) = myelD3(ii)/(condTimeMyelF_D3(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD3
            condTimeNodeF_D3(ii) = TspikeF_D3(2*ii) - TspikeF_D3(2*ii-1);   % [ms]
            condSpeedNodeF_D3(ii) = nodeD3(ii)/(condTimeNodeF_D3(ii)*1000);  % [m/s]
        end

    
    end
   
    
        
    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D4_demyelinated" axon

    [ix_affectedD4, ix_unaffectedD4, node_demyelD4, ix_demyelD4, numMyelD4, numNodeD4, myelD4, nodeD4, interfaceD4, interfaceAllD4] = axonDemyel_D4_callosal(axonLength, numMyelD3, numNodeD3, myelD3, nodeD3, damage, myelExtra, nodeExtra);

    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD4 + numNodeD4 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D4 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD4 - 1;
    ix_unaffectedAxon_D4 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD4 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD4(1)-l0)/NX:interfaceAllD4(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD4)
        xgrid(:, i) = interfaceAllD4(i-1):(interfaceAllD4(i)-interfaceAllD4(i-1))/NX:interfaceAllD4(i);
    end

    xgrid(:, m) = interfaceAllD4(end):(lm-interfaceAllD4(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D4" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D4 = @(x) ones(size(x));
    u_D4 = EL*u0_D4(xgrid);    

    b_D4 = zeros(N, 1);                 % vector to store u

    usoln_D4 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D4 = zeros(N,1);                % leak current 
    iD_D4 = zeros(N,1);                % depolarizing current
    iR_D4 = zeros(N,1);                % repolarizing current  
    iI_D4 = zeros(N,1);                % injected current
    iTotal_D4 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D4 = zeros(N,1);                % main
    yR_D4 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D4 = ones(N,1);            % ready to spike 
    spStart_D4 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D4(1) = u_D4(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D4(r) = u_D4(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D4(r) = u_D4(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D4(r) = u_D4(j,i);
        end

        b_D4(end) = u_D4(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D4) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D4(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D4(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D4(r) = GL*(EL - b_D4(r));
                iD_D4(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D4(r) - Vth)/Kth)));
                iR_D4(r) = GL*xR_D4(r)*(EL - b_D4(r));

                iTotal_D4(r) = (dt/Cc)*(iL_D4(r) + iD_D4(r) + iR_D4(r) + iI_D4(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D4
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D4(r) = GL*(EL - b_D4(r));
                iD_D4(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D4(r) - VthD)/Kth)));
                iR_D4(r) = GL*xR_D4(r)*(EL - b_D4(r));

                iTotal_D4(r) = (dt/Cc)*(iL_D4(r) + iD_D4(r) + iR_D4(r) + iI_D4(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D4 = A\(b_D4 + iTotal_D4);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D4(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D4(r) >= Vsp) 
                        spStart_D4(r) = 1; % now started spiking
                        spReady_D4(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D4(r) = 0; % not started spiking yet
                        spReady_D4(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D4(r) <= Vre)
                        spStart_D4(r) = 0; % not starting spiking, of course 
                        spReady_D4(r) = 1; % and now ready for next spike
                    else
                        spStart_D4(r) = 0; % not starting spiking, of course 
                        spReady_D4(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D4 = ee*yR_D4 + aa*spStart_D4;
        xR_D4 = ee*xR_D4 + dt*yR_D4;

        u_D4 = reshape(soln_D4, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D4(:,k) = soln_D4;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D4 = usoln_D4;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D4 = zeros(1, length(interfaceD4)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D4 = zeros(1, length(interfaceD4)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D4 = zeros(1, numMyelD4);               % vector of the times along myelinated segments
    condSpeedMyel_D4 = zeros(1, numMyelD4);              % vector of the speeds along myelinated axons

    condTimeNode_D4 = zeros(1, numNodeD4);               % vector of the times along exposed segments
    condSpeedNode_D4 = zeros(1, numNodeD4);              % vector of the speeds along exposed axons



    TspikeIndexF_D4 = zeros(1, length(interfaceD4)+2);    
    TspikeF_D4 = zeros(1, length(interfaceD4)+2);     

    condTimeMyelF_D4 = zeros(1, numMyelD4);               
    condSpeedMyelF_D4 = zeros(1, numMyelD4);              

    condTimeNodeF_D4 = zeros(1, numNodeD4);               
    condSpeedNodeF_D4 = zeros(1, numNodeD4);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D4((length(interfaceAllD4)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D4(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD4) + 2);

    for i=1:length(interfaceD4) + 2

        if isempty(find(u_D4((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D4 demyelination \n', trial)

        FailSuccess_D4(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD4) + 2
            TspikeIndex_D4(i) = find(u_D4((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D4(i) = tspan(TspikeIndex_D4(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD4
            condTimeMyel_D4(ii) = Tspike_D4(2*ii+1) - Tspike_D4(2*ii);   % [ms]
            condSpeedMyel_D4(ii) = myelD4(ii)/(condTimeMyel_D4(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD4
            condTimeNode_D4(ii) = Tspike_D4(2*ii) - Tspike_D4(2*ii-1);   % [ms]
            condSpeedNode_D4(ii) = nodeD4(ii)/(condTimeNode_D4(ii)*1000);  % [m/s]
        end

        tau_D4(trial) = sum(condTimeMyel_D4) + sum(condTimeNode_D4);     % conduction delay of each trial [ms]
        c_D4(trial) = axonLength/(tau_D4(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D4 demyelinated axon = %f ms', tau_D4(trial))
        fprintf('\n conduction velocity of D4 demyelinated axon = %f m/s \n', c_D4(trial))

    else

        fprintf('axon # is %d --> AP failure after D4 demyelination \n', trial)

        FailSuccess_D4(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD4) + 2

            if isempty(find(u_D4((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D4(i) = find(u_D4((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D4(i) = tspan(TspikeIndexF_D4(i));

            else

                TspikeIndexF_D4(i) = NaN;
                TspikeF_D4(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD4
            condTimeMyelF_D4(ii) = TspikeF_D4(2*ii+1) - TspikeF_D4(2*ii);   % [ms]
            condSpeedMyelF_D4(ii) = myelD4(ii)/(condTimeMyelF_D4(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD4
            condTimeNodeF_D4(ii) = TspikeF_D4(2*ii) - TspikeF_D4(2*ii-1);   % [ms]
            condSpeedNodeF_D4(ii) = nodeD4(ii)/(condTimeNodeF_D4(ii)*1000);  % [m/s]
        end

    end
    

        
    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D5_demyelinated" axon

    [ix_affectedD5, ix_unaffectedD5, node_demyelD5, ix_demyelD5, numMyelD5, numNodeD5, myelD5, nodeD5, interfaceD5, interfaceAllD5] = axonDemyel_D5_callosal(axonLength, numMyelD4, numNodeD4, myelD4, nodeD4, damage, myelExtra, nodeExtra);

    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD5 + numNodeD5 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D5 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD5 - 1;
    ix_unaffectedAxon_D5 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD5 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD5(1)-l0)/NX:interfaceAllD5(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD5)
        xgrid(:, i) = interfaceAllD5(i-1):(interfaceAllD5(i)-interfaceAllD5(i-1))/NX:interfaceAllD5(i);
    end

    xgrid(:, m) = interfaceAllD5(end):(lm-interfaceAllD5(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D5" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D5 = @(x) ones(size(x));
    u_D5 = EL*u0_D5(xgrid);    

    b_D5 = zeros(N, 1);                 % vector to store u

    usoln_D5 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D5 = zeros(N,1);                % leak current 
    iD_D5 = zeros(N,1);                % depolarizing current
    iR_D5 = zeros(N,1);                % repolarizing current  
    iI_D5 = zeros(N,1);                % injected current
    iTotal_D5 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D5 = zeros(N,1);                % main
    yR_D5 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D5 = ones(N,1);            % ready to spike 
    spStart_D5 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D5(1) = u_D5(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D5(r) = u_D5(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D5(r) = u_D5(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D5(r) = u_D5(j,i);
        end

        b_D5(end) = u_D5(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D5) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D5(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D5(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D5(r) = GL*(EL - b_D5(r));
                iD_D5(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D5(r) - Vth)/Kth)));
                iR_D5(r) = GL*xR_D5(r)*(EL - b_D5(r));

                iTotal_D5(r) = (dt/Cc)*(iL_D5(r) + iD_D5(r) + iR_D5(r) + iI_D5(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D5
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D5(r) = GL*(EL - b_D5(r));
                iD_D5(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D5(r) - VthD)/Kth)));
                iR_D5(r) = GL*xR_D5(r)*(EL - b_D5(r));

                iTotal_D5(r) = (dt/Cc)*(iL_D5(r) + iD_D5(r) + iR_D5(r) + iI_D5(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D5 = A\(b_D5 + iTotal_D5);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D5(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D5(r) >= Vsp) 
                        spStart_D5(r) = 1; % now started spiking
                        spReady_D5(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D5(r) = 0; % not started spiking yet
                        spReady_D5(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D5(r) <= Vre)
                        spStart_D5(r) = 0; % not starting spiking, of course 
                        spReady_D5(r) = 1; % and now ready for next spike
                    else
                        spStart_D5(r) = 0; % not starting spiking, of course 
                        spReady_D5(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D5 = ee*yR_D5 + aa*spStart_D5;
        xR_D5 = ee*xR_D5 + dt*yR_D5;

        u_D5 = reshape(soln_D5, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D5(:,k) = soln_D5;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D5 = usoln_D5;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D5 = zeros(1, length(interfaceD5)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D5 = zeros(1, length(interfaceD5)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D5 = zeros(1, numMyelD5);               % vector of the times along myelinated segments
    condSpeedMyel_D5 = zeros(1, numMyelD5);              % vector of the speeds along myelinated axons

    condTimeNode_D5 = zeros(1, numNodeD5);               % vector of the times along exposed segments
    condSpeedNode_D5 = zeros(1, numNodeD5);              % vector of the speeds along exposed axons



    TspikeIndexF_D5 = zeros(1, length(interfaceD5)+2);    
    TspikeF_D5 = zeros(1, length(interfaceD5)+2);     

    condTimeMyelF_D5 = zeros(1, numMyelD5);               
    condSpeedMyelF_D5 = zeros(1, numMyelD5);              

    condTimeNodeF_D5 = zeros(1, numNodeD5);               
    condSpeedNodeF_D5 = zeros(1, numNodeD5);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D5((length(interfaceAllD5)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D5(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD5) + 2);

    for i=1:length(interfaceD5) + 2

        if isempty(find(u_D5((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D5 demyelination \n', trial)

        FailSuccess_D5(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD5) + 2
            TspikeIndex_D5(i) = find(u_D5((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D5(i) = tspan(TspikeIndex_D5(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD5
            condTimeMyel_D5(ii) = Tspike_D5(2*ii+1) - Tspike_D5(2*ii);   % [ms]
            condSpeedMyel_D5(ii) = myelD5(ii)/(condTimeMyel_D5(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD5
            condTimeNode_D5(ii) = Tspike_D5(2*ii) - Tspike_D5(2*ii-1);   % [ms]
            condSpeedNode_D5(ii) = nodeD5(ii)/(condTimeNode_D5(ii)*1000);  % [m/s]
        end

        tau_D5(trial) = sum(condTimeMyel_D5) + sum(condTimeNode_D5);     % conduction delay of each trial [ms]
        c_D5(trial) = axonLength/(tau_D5(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D5 demyelinated axon = %f ms', tau_D5(trial))
        fprintf('\n conduction velocity of D5 demyelinated axon = %f m/s \n', c_D5(trial))

    else

        fprintf('axon # is %d --> AP failure after D5 demyelination \n', trial)

        FailSuccess_D5(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD5) + 2

            if isempty(find(u_D5((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D5(i) = find(u_D5((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D5(i) = tspan(TspikeIndexF_D5(i));

            else

                TspikeIndexF_D5(i) = NaN;
                TspikeF_D5(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD5
            condTimeMyelF_D5(ii) = TspikeF_D5(2*ii+1) - TspikeF_D5(2*ii);   % [ms]
            condSpeedMyelF_D5(ii) = myelD5(ii)/(condTimeMyelF_D5(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD5
            condTimeNodeF_D5(ii) = TspikeF_D5(2*ii) - TspikeF_D5(2*ii-1);   % [ms]
            condSpeedNodeF_D5(ii) = nodeD5(ii)/(condTimeNodeF_D5(ii)*1000);  % [m/s]
        end

    end
    

    
    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------

    %% --------------------------------------------------------------------



    %% --------------------------------------------------------------------
    % call the function for the "D6_demyelinated" axon

    [ix_affectedD6, ix_unaffectedD6, node_demyelD6, ix_demyelD6, numMyelD6, numNodeD6, myelD6, nodeD6, interfaceD6, interfaceAllD6] = axonDemyel_D6_callosal(axonLength, numMyelD5, numNodeD5, myelD5, nodeD5, damage, myelExtra, nodeExtra);

    %% --------------------------------------------------------------------
    % reference axon 
    
    m = numMyelD6 + numNodeD6 + 2*length(myelExtra) + 2*length(nodeExtra);  % total number of segments along the entire axon (unmyelinated + myelinated + extra)

    ix_affectedAxon_D6 = length(myelExtra) + length(nodeExtra) + 2*ix_affectedD6 - 1;
    ix_unaffectedAxon_D6 = length(myelExtra) + length(nodeExtra) + 2*ix_unaffectedD6 - 1;
    
    %% --------------------------------------------------------------------
    % geometric properties of finite difference method (FDM) 
    % (setup the spatial mesh and compartmental spacing)

    % spatial mesh
    xgrid = zeros(NX+1, m);                                             % rows: number of equally-spaced compartments, columns: total number of segments

    xgrid(:, 1) = l0:(interfaceAllD6(1)-l0)/NX:interfaceAllD6(1);           % segment #1 (left extra space)

    % segments #2, ... , #m-1 of the entire axon
    for i = 2:length(interfaceAllD6)
        xgrid(:, i) = interfaceAllD6(i-1):(interfaceAllD6(i)-interfaceAllD6(i-1))/NX:interfaceAllD6(i);
    end

    xgrid(:, m) = interfaceAllD6(end):(lm-interfaceAllD6(end))/NX:lm;       % segment #m (right extra space)

    dx = diff(xgrid);                                                   % grid spacing (compartmental spacing)
    
    
    %% --------------------------------------------------------------------
    % define the diffusivity constants
    kappa = zeros(1, m);

    for i=1:2:m
        kappa(i) = (a/(2*r_int))*(1/Cc);    % unmyelinated segments
    end

    for i=2:2:m
        kappa(i) = (1/Rc)*(1/Cmy);          % myelinated segments
    end
    
    
    %% --------------------------------------------------------------------
    % FD matrix of "D6" axon

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
    
    % if rank(A) < N
    %     fprintf('axon # is %d | matrix A is NOT invertible \n', trial)
    %     break
    % else
    %     fprintf('axon # is %d | matrix A is invertible \n', trial)
    % end

    A = sparse(A);                  % sparse matrix

    
    %% --------------------------------------------------------------------
    % time loop

    t = 0;                          % set initial time [ms]
    k = 1;                          % set a counter

    tend = 40;                      % final time [ms]
    tspan = 0:dt:tend;              % time discretization

    % define the initial condition
    u0_D6 = @(x) ones(size(x));
    u_D6 = EL*u0_D6(xgrid);    

    b_D6 = zeros(N, 1);                 % vector to store u

    usoln_D6 = zeros(N, length(tspan)); % numerical solution (action potential)

    xgrid = reshape(xgrid, N, 1);     % reshape the xgrid from a matrix to a column vector to perform the calculations

    
    % initialize currents     
    iL_D6 = zeros(N,1);                % leak current 
    iD_D6 = zeros(N,1);                % depolarizing current
    iR_D6 = zeros(N,1);                % repolarizing current  
    iI_D6 = zeros(N,1);                % injected current
    iTotal_D6 = zeros(N,1);            % all currents

    % repolarizing conductances (alpha function)
    xR_D6 = zeros(N,1);                % main
    yR_D6 = zeros(N,1);                % sub

    % define vectors for threshold checking
    spReady_D6 = ones(N,1);            % ready to spike 
    spStart_D6 = zeros(N,1);           % start spiking
    
    
    while t <= tend

        t = t + dt;                 % increase the time by one timestep

        b_D6(1) = u_D6(1,1);              % left boundary

        % mid compartmental points within each segment
        for i = 1:m
            for j = 2:NX
                r = (i-1)*(NX+1) + j;

                b_D6(r) = u_D6(j,i);
            end
        end

        % interface compartmental point (left)
        for i = 2:m
            j = 1;
            r = (i-1)*(NX+1) + j;

            b_D6(r) = u_D6(j,i);
        end

        % interface compartmental point (right)
        for i = 1:m-1
            j = NX+1;
            r = (i-1)*(NX+1) + j;

            b_D6(r) = u_D6(j,i);
        end

        b_D6(end) = u_D6(NX+1,m);           % right boundary


        %% define the ion currents

        % the injected current (iI_D6) is zero everywhere except at a particular 
        % time period

        if t>=10
            iI_D6(3*(NX-4)) = 1e-2;       % inject current at the 3rd segment [pA]
        end

        if t>=11
            iI_D6(3*(NX-4)) = 0;
        end

        
        % node segments along the entire axon
        for i=1:2:m %
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D6(r) = GL*(EL - b_D6(r));
                iD_D6(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D6(r) - Vth)/Kth)));
                iR_D6(r) = GL*xR_D6(r)*(EL - b_D6(r));

                iTotal_D6(r) = (dt/Cc)*(iL_D6(r) + iD_D6(r) + iR_D6(r) + iI_D6(r));
            end
        end
        
        % affected node segments along the entire axon
        for i=ix_affectedAxon_D6
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                iL_D6(r) = GL*(EL - b_D6(r));
                iD_D6(r) = GL*Kth*(Ath/(1 + Ath*exp(-(b_D6(r) - VthD)/Kth)));
                iR_D6(r) = GL*xR_D6(r)*(EL - b_D6(r));

                iTotal_D6(r) = (dt/Cc)*(iL_D6(r) + iD_D6(r) + iR_D6(r) + iI_D6(r));
            end
        end

        
        
        % calculate the numerical solution
        soln_D6 = A\(b_D6 + iTotal_D6);

        %% check for threshold crossing

        % UNmyelinated segments along the entire axon
        for i=1:2:m
            for j=1:NX+1
                r = (i-1)*(NX+1) + j;

                if (spReady_D6(r) == 1)    % if ready to spike, then check for spike occurrence
                    if (soln_D6(r) >= Vsp) 
                        spStart_D6(r) = 1; % now started spiking
                        spReady_D6(r) = 0; % thus not ready for spiking
                    else 
                        spStart_D6(r) = 0; % not started spiking yet
                        spReady_D6(r) = 1; % thus still ready for spiking
                    end
                else                    % if in repolarizing phase, then check whether voltage is back near rest
                    if (soln_D6(r) <= Vre)
                        spStart_D6(r) = 0; % not starting spiking, of course 
                        spReady_D6(r) = 1; % and now ready for next spike
                    else
                        spStart_D6(r) = 0; % not starting spiking, of course 
                        spReady_D6(r) = 0; % not yet ready for next spike
                    end
                end
            end
        end

        %% ----------------------------------------------------------------
        yR_D6 = ee*yR_D6 + aa*spStart_D6;
        xR_D6 = ee*xR_D6 + dt*yR_D6;

        u_D6 = reshape(soln_D6, NX+1, m); % reshape "soln" from a vector to a matrix

        usoln_D6(:,k) = soln_D6;          % store solution at requested values in tspan
        k = k + 1;                  % increase the counter

        % if the value of the counter k exceeds the length of tspan then
        % terminate the while loop

        if k == length(tspan)+1
            break;
        end

        
    end % end of the while loop
    
    u_D6 = usoln_D6;                      % action potential along the entire axon
    x = xgrid;                      % vector of the grid points in space

    %% --------------------------------------------------------------------
    % Compute the time where the AP crosses the threshold value at the 
    % interfaces and at the starting and ending points of the reference axon

    TspikeIndex_D6 = zeros(1, length(interfaceD6)+2);    % vector of the indices where the AP crosses the threshold
    Tspike_D6 = zeros(1, length(interfaceD6)+2);         % vector of the corresponding time at which the AP crosses the threshold

    condTimeMyel_D6 = zeros(1, numMyelD6);               % vector of the times along myelinated segments
    condSpeedMyel_D6 = zeros(1, numMyelD6);              % vector of the speeds along myelinated axons

    condTimeNode_D6 = zeros(1, numNodeD6);               % vector of the times along exposed segments
    condSpeedNode_D6 = zeros(1, numNodeD6);              % vector of the speeds along exposed axons



    TspikeIndexF_D6 = zeros(1, length(interfaceD6)+2);    
    TspikeF_D6 = zeros(1, length(interfaceD6)+2);     

    condTimeMyelF_D6 = zeros(1, numMyelD6);               
    condSpeedMyelF_D6 = zeros(1, numMyelD6);              

    condTimeNodeF_D6 = zeros(1, numNodeD6);               
    condSpeedNodeF_D6 = zeros(1, numNodeD6);              

    
    % Check for conduction success/failure: conduction success occurs if at 
    % the ending point of the reference axon the value of the potential is 
    % above the threshold. In this case, compute the index and the 
    % corresponding time at which u surpasses Vth. Otherwise, continue to 
    % the next iteration as conduction failure occurs.

    if isempty(find((u_D6((length(interfaceAllD6)-2*length(myelExtra))*(NX+1),end))>Vth, 1)) == 0

        fprintf('axon # is %d --> unstable \n', trial)

        FailSuccess_D6(trial) = -1;     % the "-1" indicates instability

    end

    
    pointer = NaN(1,length(interfaceD6) + 2);

    for i=1:length(interfaceD6) + 2

        if isempty(find(u_D6((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) ==0

            pointer(i) = 1;         % success
        
        else

            pointer(i) = 0;         % fail

        end

    end
    
    
    if nnz(~pointer) == 0

        fprintf('axon # is %d --> AP success after D6 demyelination \n', trial)

        FailSuccess_D6(trial) = 1;     % the "1" indicates success

        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD6) + 2
            TspikeIndex_D6(i) = find(u_D6((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
            Tspike_D6(i) = tspan(TspikeIndex_D6(i));
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD6
            condTimeMyel_D6(ii) = Tspike_D6(2*ii+1) - Tspike_D6(2*ii);   % [ms]
            condSpeedMyel_D6(ii) = myelD6(ii)/(condTimeMyel_D6(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD6
            condTimeNode_D6(ii) = Tspike_D6(2*ii) - Tspike_D6(2*ii-1);   % [ms]
            condSpeedNode_D6(ii) = nodeD6(ii)/(condTimeNode_D6(ii)*1000);  % [m/s]
        end

        tau_D6(trial) = sum(condTimeMyel_D6) + sum(condTimeNode_D6);     % conduction delay of each trial [ms]
        c_D6(trial) = axonLength/(tau_D6(trial)*1000);                % conduction velocity of each trial [m/s]

        fprintf('\n conduction delay of D6 demyelinated axon = %f ms', tau_D6(trial))
        fprintf('\n conduction velocity of D6 demyelinated axon = %f m/s \n', c_D6(trial))

    else

        fprintf('axon # is %d --> AP failure after D6 demyelination \n', trial)

        FailSuccess_D6(trial) = 0;     % the "0" indicates failure

        
        % index and corresponding time for which u surpasses Vth at the
        % interfaces and at the starting and ending points of the reference
        % axon

        for i=1:length(interfaceD6) + 2

            if isempty(find(u_D6((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1)) == 0
                
                TspikeIndexF_D6(i) = find(u_D6((length(myelExtra)*2+i-1)*(NX+1),:)>Vth, 1);
                TspikeF_D6(i) = tspan(TspikeIndexF_D6(i));

            else

                TspikeIndexF_D6(i) = NaN;
                TspikeF_D6(i) = NaN;

            end

        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along myelinated segments

        for ii=1:numMyelD6
            condTimeMyelF_D6(ii) = TspikeF_D6(2*ii+1) - TspikeF_D6(2*ii);   % [ms]
            condSpeedMyelF_D6(ii) = myelD6(ii)/(condTimeMyelF_D6(ii)*1000);  % [m/s]
        end

        %% ----------------------------------------------------------------
        % compute the conduction time and speed along exposed segments

        for ii=1:numNodeD6
            condTimeNodeF_D6(ii) = TspikeF_D6(2*ii) - TspikeF_D6(2*ii-1);   % [ms]
            condSpeedNodeF_D6(ii) = nodeD6(ii)/(condTimeNodeF_D6(ii)*1000);  % [m/s]
        end

    end
    
    
    nS = length(find(FailSuccess_D6 == 1));
    nF = length(find(FailSuccess_D6 == 0));

    if nS+nF == 100
        break;
    end
    
end % end of the for loop



tauSuccess = tau(tau ~= 0);                  % save the conduction delays of only successful trials
cSuccess = c(c ~= 0);                        % save the conduction velocities of only successful trials

tauSuccess_D1 = tau_D1(tau_D1 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D1 = c_D1(c_D1 ~= 0);               % save the conduction velocities of only successful trials

tauSuccess_D2 = tau_D2(tau_D2 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D2 = c_D2(c_D2 ~= 0);               % save the conduction velocities of only successful trials

tauSuccess_D3 = tau_D3(tau_D3 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D3 = c_D3(c_D3 ~= 0);               % save the conduction velocities of only successful trials

tauSuccess_D4 = tau_D4(tau_D4 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D4 = c_D4(c_D4 ~= 0);               % save the conduction velocities of only successful trials

tauSuccess_D5 = tau_D5(tau_D5 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D5 = c_D5(c_D5 ~= 0);               % save the conduction velocities of only successful trials

tauSuccess_D6 = tau_D6(tau_D6 ~= 0);         % save the conduction delays of only successful trials 
cSuccess_D6 = c_D6(c_D6 ~= 0);               % save the conduction velocities of only successful trials












%% ------------------------------------------------------------------------
%% define the function for the axon (perMyel = 90, case1)

function [numMyel, l0, lm, startAxon, endAxon, myel, node, interface, myelExtra, nodeExtra, interfaceAll] = axon_callosal10mm(axonLength, perMyel)

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

numMyel = randi([57, 280],1);            % number of myelinated segments along the reference axon (sample a uniformly distributed pseudorandom integer)     
numNode = numMyel;                      % number of unmyelinated segments along the reference axon
fprintf('numMyel # --> %d\n', numMyel)


%% ------------------------------------------------------------------------
% MYELINATED SEGMENTS
% sample the myelin segment lengths from the Gamma distribution

while true
    kappaMyel = randi([1,30],1); 
    thetaMyel = randi([1,15],1); 
    myel = round(gamrnd(kappaMyel, thetaMyel, [1 numMyel])); % myelin segment lengths
    
    if (sum(myel) == myelAxon && discretize(mean(myel), [20, 160]) == 1 && min(myel)>=1) 
        disp('myelin')
        break
    end 
end  

%% ------------------------------------------------------------------------
% UNMYELINATED SEGMENTS

while true
    kappaNode = randi([1,30],1);
    thetaNode = randi([1,15],1);
    node = round(gamrnd(kappaNode, thetaNode, [1 numNode])); % node segment lengths
   
    if (sum(node) == nodeAxon && min(node)>=1)
        
        disp('node')
        break
    end
end

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
    
axon = zeros(1, numMyel + numNode); 
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
    
end








%% ------------------------------------------------------------------------
% define the function for the D1 demyelination

function [ix_affectedD1, ix_unaffectedD1, node_demyelD1, ix_demyelD1, numMyelD1, numNodeD1, myelD1, nodeD1, interfaceD1, interfaceAllD1] = axonDemyel_D1_callosal(axonLength, numMyel, numNode, myel, node, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D1 (damage to the number of myelin segments of the "healthy" axon)
    
    %% pick the demyelinated segments
    numDemyelD1 = round(damage*numMyel);   % number of demyelinated segments
    ix_demyelD1 = sort(randperm(numMyel-1, numDemyelD1));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD1(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD1) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD1 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD1{i} = ix_demyelD1(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD1));     % index of affected nodes
    node_demyelD1 = zeros(1, length(seq_demyelD1));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD1)
        ix_nodeBeforeD1 = seq_demyelD1{i}(1);
        ix_nodeAfterD1 = seq_demyelD1{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD1:ix_nodeAfterD1;
        node_demyelD1(i) = sum(node(ix_nodeBeforeD1:ix_nodeAfterD1)) + sum(myel(seq_demyelD1{i}(1):seq_demyelD1{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD1 = setdiff(1:length(myel), ix_demyelD1);
    myelD1 = myel(ix_myelD1);
    numMyelD1 = length(myelD1);     % number of myelin sheaths at D1


    %% nodes (affected + unaffected)
    nodeD1 = zeros(1, numNode);

    for i=1:length(affected_nodes)        
        nodeD1(affected_nodes{i}(1)) = node_demyelD1(i);
        nodeD1(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD1 = nodeD1(~isnan(nodeD1));

    ix_affectedD1 = find(nodeD1); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(node), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD1(unaffected_nodes(i)) = node(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD1(affected_nodes{i}(1)) = node_demyelD1(i);
        nodeD1(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD1 = nodeD1(~isnan(nodeD1));

    numNodeD1 = length(nodeD1);     % number of nodes at D1
    ix_unaffectedD1 = setdiff(1:numNodeD1, ix_affectedD1);
    
   
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD1 = zeros(1, numMyelD1 + numNodeD1); 
    axonD1(1:2:end) = nodeD1;                           % odd segments are unmyelinated
    axonD1(2:2:end) = myelD1;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD1 = zeros(1, length(axonD1) - 1);

    interfaceD1(1) = startAxon + axonD1(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD1)
        interfaceD1(j) = interfaceD1(j-1) + axonD1(j);
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

    interfaceAllD1 = [interfaceStart startAxon interfaceD1 endAxon interfaceEnd];
    
end






%% ------------------------------------------------------------------------
% define the function for the D2 demyelination

function [ix_affectedD2, ix_unaffectedD2, node_demyelD2, ix_demyelD2, numMyelD2, numNodeD2, myelD2, nodeD2, interfaceD2, interfaceAllD2] = axonDemyel_D2_callosal(axonLength, numMyelD1, numNodeD1, myelD1, nodeD1, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D2 (damage to the number of myelin segments of the D1 axon)
    
    %% pick the demyelinated segments
    numDemyelD2 = round(damage*numMyelD1);   % number of demyelinated segments
    ix_demyelD2 = sort(randperm(numMyelD1 - 1, numDemyelD2));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD2(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD2) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD2 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD2{i} = ix_demyelD2(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD2));     % index of affected nodes
    node_demyelD2 = zeros(1, length(seq_demyelD2));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD2)
        ix_nodeBeforeD2 = seq_demyelD2{i}(1);
        ix_nodeAfterD2 = seq_demyelD2{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD2:ix_nodeAfterD2;
        node_demyelD2(i) = sum(nodeD1(ix_nodeBeforeD2:ix_nodeAfterD2)) + sum(myelD1(seq_demyelD2{i}(1):seq_demyelD2{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD2 = setdiff(1:length(myelD1), ix_demyelD2);
    myelD2 = myelD1(ix_myelD2);
    numMyelD2 = length(myelD2);     % number of myelin sheaths at D1


    %% nodes (affected + unaffected)
    nodeD2 = zeros(1, numNodeD1);

    for i=1:length(affected_nodes)        
        nodeD2(affected_nodes{i}(1)) = node_demyelD2(i);
        nodeD2(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD2 = nodeD2(~isnan(nodeD2));

    ix_affectedD2 = find(nodeD2); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(nodeD1), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD2(unaffected_nodes(i)) = nodeD1(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD2(affected_nodes{i}(1)) = node_demyelD2(i);
        nodeD2(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD2 = nodeD2(~isnan(nodeD2));

    numNodeD2 = length(nodeD2);     % number of nodes at D2
    ix_unaffectedD2 = setdiff(1:numNodeD2, ix_affectedD2);
  
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD2 = zeros(1, numMyelD2 + numNodeD2); 
    axonD2(1:2:end) = nodeD2;                           % odd segments are unmyelinated
    axonD2(2:2:end) = myelD2;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD2 = zeros(1, length(axonD2) - 1);

    interfaceD2(1) = startAxon + axonD2(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD2)
        interfaceD2(j) = interfaceD2(j-1) + axonD2(j);
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

    interfaceAllD2 = [interfaceStart startAxon interfaceD2 endAxon interfaceEnd];
    
end





%% ------------------------------------------------------------------------
% define the function for the D3 demyelination

function [ix_affectedD3, ix_unaffectedD3, node_demyelD3, ix_demyelD3, numMyelD3, numNodeD3, myelD3, nodeD3, interfaceD3, interfaceAllD3] = axonDemyel_D3_callosal(axonLength, numMyelD2, numNodeD2, myelD2, nodeD2, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D3 (damage to the number of myelin segments of the D2 axon)
    
    %% pick the demyelinated segments
    numDemyelD3 = round(damage*numMyelD2);   % number of demyelinated segments
    ix_demyelD3 = sort(randperm(numMyelD2 - 1, numDemyelD3));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD3(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD3) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD3 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD3{i} = ix_demyelD3(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD3));     % index of affected nodes
    node_demyelD3 = zeros(1, length(seq_demyelD3));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD3)
        ix_nodeBeforeD3 = seq_demyelD3{i}(1);
        ix_nodeAfterD3 = seq_demyelD3{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD3:ix_nodeAfterD3;
        node_demyelD3(i) = sum(nodeD2(ix_nodeBeforeD3:ix_nodeAfterD3)) + sum(myelD2(seq_demyelD3{i}(1):seq_demyelD3{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD3 = setdiff(1:length(myelD2), ix_demyelD3);
    myelD3 = myelD2(ix_myelD3);
    numMyelD3 = length(myelD3);     % number of myelin sheaths at D2


    %% nodes (affected + unaffected)
    nodeD3 = zeros(1, numNodeD2);

    for i=1:length(affected_nodes)        
        nodeD3(affected_nodes{i}(1)) = node_demyelD3(i);
        nodeD3(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD3 = nodeD3(~isnan(nodeD3));

    ix_affectedD3 = find(nodeD3); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(nodeD2), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD3(unaffected_nodes(i)) = nodeD2(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD3(affected_nodes{i}(1)) = node_demyelD3(i);
        nodeD3(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD3 = nodeD3(~isnan(nodeD3));

    numNodeD3 = length(nodeD3);     % number of nodes at D3
    ix_unaffectedD3 = setdiff(1:numNodeD3, ix_affectedD3);
  
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD3 = zeros(1, numMyelD3 + numNodeD3); 
    axonD3(1:2:end) = nodeD3;                           % odd segments are unmyelinated
    axonD3(2:2:end) = myelD3;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD3 = zeros(1, length(axonD3) - 1);

    interfaceD3(1) = startAxon + axonD3(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD3)
        interfaceD3(j) = interfaceD3(j-1) + axonD3(j);
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

    interfaceAllD3 = [interfaceStart startAxon interfaceD3 endAxon interfaceEnd];
    
end





%% ------------------------------------------------------------------------
% define the function for the D4 demyelination

function [ix_affectedD4, ix_unaffectedD4, node_demyelD4, ix_demyelD4, numMyelD4, numNodeD4, myelD4, nodeD4, interfaceD4, interfaceAllD4] = axonDemyel_D4_callosal(axonLength, numMyelD3, numNodeD3, myelD3, nodeD3, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D4 (damage to the number of myelin segments of the D3 axon)
    
    %% pick the demyelinated segments
    numDemyelD4 = round(damage*numMyelD3);   % number of demyelinated segments
    ix_demyelD4 = sort(randperm(numMyelD3 - 1, numDemyelD4));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD4(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD4) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD4 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD4{i} = ix_demyelD4(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD4));     % index of affected nodes
    node_demyelD4 = zeros(1, length(seq_demyelD4));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD4)
        ix_nodeBeforeD4 = seq_demyelD4{i}(1);
        ix_nodeAfterD4 = seq_demyelD4{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD4:ix_nodeAfterD4;
        node_demyelD4(i) = sum(nodeD3(ix_nodeBeforeD4:ix_nodeAfterD4)) + sum(myelD3(seq_demyelD4{i}(1):seq_demyelD4{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD4 = setdiff(1:length(myelD3), ix_demyelD4);
    myelD4 = myelD3(ix_myelD4);
    numMyelD4 = length(myelD4);     % number of myelin sheaths at D3


    %% nodes (affected + unaffected)
    nodeD4 = zeros(1, numNodeD3);

    for i=1:length(affected_nodes)        
        nodeD4(affected_nodes{i}(1)) = node_demyelD4(i);
        nodeD4(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD4 = nodeD4(~isnan(nodeD4));

    ix_affectedD4 = find(nodeD4); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(nodeD3), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD4(unaffected_nodes(i)) = nodeD3(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD4(affected_nodes{i}(1)) = node_demyelD4(i);
        nodeD4(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD4 = nodeD4(~isnan(nodeD4));

    numNodeD4 = length(nodeD4);     % number of nodes at D4
    ix_unaffectedD4 = setdiff(1:numNodeD4, ix_affectedD4);
  
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD4 = zeros(1, numMyelD4 + numNodeD4); 
    axonD4(1:2:end) = nodeD4;                           % odd segments are unmyelinated
    axonD4(2:2:end) = myelD4;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD4 = zeros(1, length(axonD4) - 1);

    interfaceD4(1) = startAxon + axonD4(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD4)
        interfaceD4(j) = interfaceD4(j-1) + axonD4(j);
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

    interfaceAllD4 = [interfaceStart startAxon interfaceD4 endAxon interfaceEnd];
    
end





%% ------------------------------------------------------------------------
% define the function for the D5 demyelination

function [ix_affectedD5, ix_unaffectedD5, node_demyelD5, ix_demyelD5, numMyelD5, numNodeD5, myelD5, nodeD5, interfaceD5, interfaceAllD5] = axonDemyel_D5_callosal(axonLength, numMyelD4, numNodeD4, myelD4, nodeD4, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D5 (damage to the number of myelin segments of the D4 axon)
    
    %% pick the demyelinated segments
    numDemyelD5 = round(damage*numMyelD4);   % number of demyelinated segments
    ix_demyelD5 = sort(randperm(numMyelD4 - 1, numDemyelD5));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD5(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD5) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD5 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD5{i} = ix_demyelD5(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD5));     % index of affected nodes
    node_demyelD5 = zeros(1, length(seq_demyelD5));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD5)
        ix_nodeBeforeD5 = seq_demyelD5{i}(1);
        ix_nodeAfterD5 = seq_demyelD5{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD5:ix_nodeAfterD5;
        node_demyelD5(i) = sum(nodeD4(ix_nodeBeforeD5:ix_nodeAfterD5)) + sum(myelD4(seq_demyelD5{i}(1):seq_demyelD5{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD5 = setdiff(1:length(myelD4), ix_demyelD5);
    myelD5 = myelD4(ix_myelD5);
    numMyelD5 = length(myelD5);     % number of myelin sheaths at D4


    %% nodes (affected + unaffected)
    nodeD5 = zeros(1, numNodeD4);

    for i=1:length(affected_nodes)        
        nodeD5(affected_nodes{i}(1)) = node_demyelD5(i);
        nodeD5(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD5 = nodeD5(~isnan(nodeD5));

    ix_affectedD5 = find(nodeD5); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(nodeD4), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD5(unaffected_nodes(i)) = nodeD4(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD5(affected_nodes{i}(1)) = node_demyelD5(i);
        nodeD5(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD5 = nodeD5(~isnan(nodeD5));

    numNodeD5 = length(nodeD5);     % number of nodes at D5
    ix_unaffectedD5 = setdiff(1:numNodeD5, ix_affectedD5);
  
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD5 = zeros(1, numMyelD5 + numNodeD5); 
    axonD5(1:2:end) = nodeD5;                           % odd segments are unmyelinated
    axonD5(2:2:end) = myelD5;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD5 = zeros(1, length(axonD5) - 1);

    interfaceD5(1) = startAxon + axonD5(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD5)
        interfaceD5(j) = interfaceD5(j-1) + axonD5(j);
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

    interfaceAllD5 = [interfaceStart startAxon interfaceD5 endAxon interfaceEnd];
    
end





%% ------------------------------------------------------------------------
% define the function for the D6 demyelination

function [ix_affectedD6, ix_unaffectedD6, node_demyelD6, ix_demyelD6, numMyelD6, numNodeD6, myelD6, nodeD6, interfaceD6, interfaceAllD6] = axonDemyel_D6_callosal(axonLength, numMyelD5, numNodeD5, myelD5, nodeD5, damage, myelExtra, nodeExtra)

    %% ---------------------------------------------------------------------
    % D6 (damage to the number of myelin segments of the D5 axon)
    
    %% pick the demyelinated segments
    numDemyelD6 = round(damage*numMyelD5);   % number of demyelinated segments
    ix_demyelD6 = sort(randperm(numMyelD5 - 1, numDemyelD6));   % from the myelin sheaths choose randomly (and without repetition) those that will be removed


    %% identify consecutive demyelinated segments
    ix_demyelD6(end + 1) = 1000000;   % adds new endpoint to very end so code picks up end of last group of consecutive values
    I = find(diff(ix_demyelD6) ~= 1);  % finds where sequences of consecutive numbers end
    [m, n] = size(I);   % finds dimensions of I, i.e. how many sequences of consecutive numbers you have
    start_point = 1;    % sets start index at the first value in your array

    seq_demyelD6 = cell(1,n);  % preallocate cell for the demyelinated segments
        
    for i = 1:n
        end_ix = I(i);   % set end index
        seq_demyelD6{i} = ix_demyelD6(start_point:end_ix);  % finds sequences of consecutive numbers and assigns to cell array
        start_point = end_ix + 1;   % update start index for the next consecutive sequence
    end

    
    %% find the lengths of the newly formed nodes
    affected_nodes = cell(1, length(seq_demyelD6));     % index of affected nodes
    node_demyelD6 = zeros(1, length(seq_demyelD6));     % lengths of the newly formed nodes 

    for i = 1:length(seq_demyelD6)
        ix_nodeBeforeD6 = seq_demyelD6{i}(1);
        ix_nodeAfterD6 = seq_demyelD6{i}(end) + 1;

        affected_nodes{i} = ix_nodeBeforeD6:ix_nodeAfterD6;
        node_demyelD6(i) = sum(nodeD5(ix_nodeBeforeD6:ix_nodeAfterD6)) + sum(myelD5(seq_demyelD6{i}(1):seq_demyelD6{i}(end)));
    end


    %% myelin sheaths (index + lengths) that are not affected by demyelination
    ix_myelD6 = setdiff(1:length(myelD5), ix_demyelD6);
    myelD6 = myelD5(ix_myelD6);
    numMyelD6 = length(myelD6);     % number of myelin sheaths at D5


    %% nodes (affected + unaffected)
    nodeD6 = zeros(1, numNodeD5);

    for i=1:length(affected_nodes)        
        nodeD6(affected_nodes{i}(1)) = node_demyelD6(i);
        nodeD6(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD6 = nodeD6(~isnan(nodeD6));

    ix_affectedD6 = find(nodeD6); % indices of affected segments along the new axon


    unaffected_nodes = setdiff(1:length(nodeD5), cell2mat(affected_nodes)); % indices of unaffected nodes

    for i=1:length(unaffected_nodes)
        nodeD6(unaffected_nodes(i)) = nodeD5(unaffected_nodes(i)); % lengths of unaffected nodes
    end

    for i=1:length(affected_nodes)
        nodeD6(affected_nodes{i}(1)) = node_demyelD6(i);
        nodeD6(affected_nodes{i}(2):affected_nodes{i}(end)) = NaN;
    end

    nodeD6 = nodeD6(~isnan(nodeD6));

    numNodeD6 = length(nodeD6);     % number of nodes at D6
    ix_unaffectedD6 = setdiff(1:numNodeD6, ix_affectedD6);
  
    
    %% --------------------------------------------------------------------
    % Define the parameters of the entire axon

    l0 = 0;                                         % starting point of the entire axon (including the extra space)
    startAxon = sum(myelExtra) + sum(nodeExtra);    % starting point of the reference axon (the reference axon is the myelinated axon where we compute the conduction time)
    endAxon = startAxon + axonLength;               % ending point of the reference axon
    lm = startAxon + endAxon;                       % ending point of the entire axon (including the extra space) 

    %% ------------------------------------------------------------------------
    % Create the reference axon

    axonD6 = zeros(1, numMyelD6 + numNodeD6); 
    axonD6(1:2:end) = nodeD6;                           % odd segments are unmyelinated
    axonD6(2:2:end) = myelD6;                           % even segments are myelinated
    
    %% ------------------------------------------------------------------------
    % interfaces along the reference axon: contact points between myelinated 
    % and unmyelinated segments

    interfaceD6 = zeros(1, length(axonD6) - 1);

    interfaceD6(1) = startAxon + axonD6(1);             % first interface

    % mid and last interfaces
    for j = 2:length(interfaceD6)
        interfaceD6(j) = interfaceD6(j-1) + axonD6(j);
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

    interfaceAllD6 = [interfaceStart startAxon interfaceD6 endAxon interfaceEnd];
    
end







