%
% Discrete reaction-diffusion on hexagonal lattice
%
% Szabolcs Zakany
% LANE, University of Geneva
% 2021

%% Reaction diffusion parameters

%(optimized for LL66 for n.n. stats):
%-------------------------
params.c1=-0.04;
params.c2=-0.056;
params.c3=0.382;
params.c4=-0.05;
params.c5=0;
params.c6=0.25;
params.c7=0.016;
params.c8=-0.0325; %!
params.c9=0.24;
params.cu=0.02;
params.cv=0.025;
params.cw=0.058; %!
params.Du=0.189; %!
params.Dv=0.189; %!
params.Dw=13.5;
params.P=0.0199; %!
params.epsilon=1;
params.S=20/3^(0.5);

dt = 0.4;

S=params.S;    % hexagon side
AreaHex = 3*3^(0.5)*S^2/2;

%% Steady state (analytical approximation)

% homogeneous unstable steady state
ffunc = @(x)(Fu(x(1),x(2),x(3),params).^2 + Fv(x(1),x(2),x(3),params).^2 + Fw(x(1),x(2),x(3),params).^2);
[hss_sol,~] = fminsearch(ffunc,[3 3 3]);
UVWanalytic_hss = hss_sol;

% main clusters
ffunc = @(x)(objfunc1stLevel(x,params));
[h1st_sol,~] = fminsearch(ffunc,[6 1 1 1 10 1]);
UVWanalyticB_1st = h1st_sol(1:3);
UVWanalyticG_1st = h1st_sol(4:6);

%  sub-clusters
for nG = 0:6
    ffunc = @(x)(objfunc2ndLevel(x,params,UVWanalyticB_1st,UVWanalyticG_1st,nG));
    [h2nd_sol,~] = fminsearch(ffunc,UVWanalyticB_1st);
    UVWanalyticB_2nd(nG+1,:) = h2nd_sol;
    [h2nd_sol,~] = fminsearch(ffunc,UVWanalyticG_1st);
    UVWanalyticG_2nd(nG+1,:) = h2nd_sol;
end



%% Periodic hexagonal lattice simulation


M = 110; %number of columns (~X coord)
N = 35; %number of rows (~Y coord)


% for juvenile-like initial conditions...
load('initialjuvenile_35x110.mat')
icJuvenile = 1; %1: use juvenile initial condition, 0: use small perturbations around hss

% L matrix (edge length matrix)
disp('Build matrices...')
L = zeros(N*M,N*M);
for n=0:N-1
    for m=0:M-1
        
        %for periodicity...
        nm1 = mod(n-1,N);
        np1 = mod(n+1,N);
        mm1 = mod((m-1),M);
        mp1 = mod((m+1),M);
        
        L( n*M+m+1, nm1*M+m+1) = S;
        L( n*M+m+1, np1*M+m+1) = S;
        L( n*M+m+1, n*M+mm1+1) = S;
        L( n*M+m+1, n*M+mp1+1) = S;
        if mod(m,2) == 0
            L( n*M+m+1, np1*M+mm1+1) = S;
            L( n*M+m+1, np1*M+mp1+1) = S;
        else
            L( n*M+m+1, nm1*M+mm1+1) = S;
            L( n*M+m+1, nm1*M+mp1+1) = S;
        end
    end
end

% M matrix
AreaHex = 3*3^(0.5)*S^2/2;
Mmatrix = (L-diag(sum(L,2)))./(AreaHex);

% time step matrices 
Mu = eye(N*M)+Mmatrix*dt*(params.Du*params.P/params.epsilon);
Mv = eye(N*M)+Mmatrix*dt*(params.Dv*params.P/params.epsilon);
Mw = eye(N*M)+Mmatrix*dt*(params.Dw*params.P/params.epsilon);


disp('Simulation...')

% build values of initial conditions
r=0.5; 
UVW_initB = UVWanalytic_hss+r*(UVWanalyticB_1st-UVWanalytic_hss);
UVW_initG = UVWanalytic_hss+r*(UVWanalyticG_1st-UVWanalytic_hss);
UicB = UVW_initB(1);
VicB = UVW_initB(2);
WicB = UVW_initB(3);
UicG = UVW_initG(1);
VicG = UVW_initG(2);
WicG = UVW_initG(3);


U = zeros(N*M,1);
V = U;
W = U;
for n=0:N-1
    for m=0:M-1
        
        if icJuvenile == 1
            % juvenile-like initial conditions...
            U(n*M+m+1) = initialJuvenile(n*M+m+1)*UicG + (1-initialJuvenile(n*M+m+1))*UicB;
            V(n*M+m+1) = initialJuvenile(n*M+m+1)*VicG + (1-initialJuvenile(n*M+m+1))*VicB;
            W(n*M+m+1) = initialJuvenile(n*M+m+1)*WicG + (1-initialJuvenile(n*M+m+1))*WicB;
        else
            %hss initial conditions...
            U(n*M+m+1) = Uhss+0.01*(2*rand()-1);
            V(n*M+m+1) = Vhss+0.01*(2*rand()-1);
            W(n*M+m+1) = Whss+0.01*(2*rand()-1);
        end
    end
end




% run the dRD simulation
%%%%%%%%%%%%%%%%%%%%%%%
tsteps = 5000;        % number of time steps before simulation ends
%%%%%%%%%%%%%%%%%%%%%%%
for step=1:tsteps
    disp(['timestep: '  num2str(step)])
    Unew = Mu*U + dt*Fu(U,V,W,params);
    Vnew = Mv*V + dt*Fv(U,V,W,params);
    Wnew = Mw*W + dt*Fw(U,V,W,params);
    U = Unew;
    V = Vnew;
    W = Wnew;
    output.Uhistory(:,step) = single(U);
    output.Vhistory(:,step) = single(V);
    output.Whistory(:,step) = single(W);
    output.thistory(step) = step*dt;
end

output.UVWanalytic_hss = UVWanalytic_hss;
output.UVWanalyticB_1st = UVWanalyticB_1st;
output.UVWanalyticG_1st = UVWanalyticG_1st;
output.UVWanalyticB_2nd = UVWanalyticB_2nd;
output.UVWanalyticG_2nd = UVWanalyticG_2nd;

save(['simulationResult_',datestr(now(),'yyyymmdd_HHMMSS'),'.mat'],'output','-v7.3');




%% plot result...

% primitive hexagon
xhex = cos( (0:5)*pi/3 );
yhex = sin( (0:5)*pi/3 );

sqrt3 = 3^(0.5);
%timepoint = 400;
timepoint = length(output.thistory); % final pattern...
U = output.Uhistory(:,timepoint);
V = output.Vhistory(:,timepoint);
W = output.Whistory(:,timepoint);


Xpatch = zeros(6,N*M);
Ypatch = Xpatch;
Cpatch = zeros(N*M,1,3);
for n=0:N-1
    for m=0:M-1
        
        c = n*M+m+1;
        xcenter = (m)*(3/2*params.S+sqrt3/2*params.epsilon);
        ycenter = -(n)*(sqrt3*params.S+params.epsilon)+mod(m,2)*(sqrt3/2*params.S+1/2*params.epsilon);
        
        Xpatch(:,c) = params.S*xhex+xcenter;
        Ypatch(:,c) = params.S*yhex+ycenter;
        Cpatch(c,1,:) = uvw2rgb(U(c),V(c),W(c));
        
    end
end
figure
subplot(2,1,1)
patch(Xpatch,Ypatch,Cpatch,'EdgeColor','none')
daspect([1 1 1])


subplot(2,1,2)
scatter(U,V,'.')
hold on
scatter(UVWanalyticB_2nd(:,1),UVWanalyticB_2nd(:,2),'rx')
scatter(UVWanalyticG_2nd(:,1),UVWanalyticG_2nd(:,2),'rx')




%% functions

function y = saturate(x)
    y = min(max(x,0),0.5);
end
function y = Fu(U,V,W,params)
    y = saturate(params.c1*V+params.c2*W+params.c3)-params.cu*U;
end
function y = Fv(U,V,W,params)
    y = saturate(params.c4*U+params.c5*W+params.c6)-params.cv*V;
end
function y = Fw(U,V,W,params)
    y = saturate(params.c7*U+params.c8*V+params.c9)-params.cw*W;
end

function rgbval = uvw2rgb(u,v,w)
    rgbval = [0 v/12 0];
end


% objective function for finding main cluster locations in uvw space
function y = objfunc1stLevel(x,params)

    % [x(1) x(2) x(3)] = [U V W] for black
    % [x(4) x(5) x(6)] = [U V W] for green

    S=params.S;    % hexagon side
    AreaHex = 3*3^(0.5)*S^2/2;
    alpha = params.P*S/(params.epsilon*AreaHex);
    
    vectB = [Fu(x(1),x(2),x(3),params)+alpha*params.Du*3*(x(4)-x(1))...
             Fv(x(1),x(2),x(3),params)+alpha*params.Dv*3*(x(5)-x(2))...
             Fw(x(1),x(2),x(3),params)+alpha*params.Dw*3*(x(6)-x(3))];
    vectG = [Fu(x(4),x(5),x(6),params)-alpha*params.Du*3*(x(4)-x(1))...
             Fv(x(4),x(5),x(6),params)-alpha*params.Dv*3*(x(5)-x(2))...
             Fw(x(4),x(5),x(6),params)-alpha*params.Dw*3*(x(6)-x(3))];     
         
    y = sum(vectB.*vectB)+sum(vectG.*vectG);
end

% objective function for finding main cluster locations in uvw space
function y = objfunc2ndLevel(x,params,UVW1stB,UVW1stG,nG)

    % [x(1) x(2) x(3)] = [U V W]

    S=params.S;    % hexagon side
    AreaHex = 3*3^(0.5)*S^2/2;
    alpha = params.P*S/(params.epsilon*AreaHex);
    
    vect = [Fu(x(1),x(2),x(3),params)+alpha*params.Du*(nG*UVW1stG(1)+(6-nG)*UVW1stB(1)-6*x(1))...
            Fv(x(1),x(2),x(3),params)+alpha*params.Dv*(nG*UVW1stG(2)+(6-nG)*UVW1stB(2)-6*x(2))...
            Fw(x(1),x(2),x(3),params)+alpha*params.Dw*(nG*UVW1stG(3)+(6-nG)*UVW1stB(3)-6*x(3))];
         
    y = sum(vect.*vect);
end




