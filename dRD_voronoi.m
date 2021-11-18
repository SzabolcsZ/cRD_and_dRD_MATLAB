%
% Discrete reaction-diffusion on arbitrary voronoi diagram
%
% Szabolcs Zakany
% LANE, University of Geneva
% 2021


%%

lizardDiagram = 0; 
%1: use a Voronoi diagram of an actual lizard skin-scale
%0: use a Voronoi diagram build from perturbed hexaognal lattice


%% Reaction diffusion parameters

%(optimized for LL66 for nearest neighbour stats):
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

dt = 2;

S=params.S;    % hexagon side
AreaHex = 3*3^(0.5)*S^2/2;


%% Prepare lattice and initial conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lizardDiagram==1  %for lizard voronoi 
    
    filestring = 'lizardVoronoi'
    load('LL38_data.mat')   
    xcenter = XYdata(:,1);
    ycenter = XYdata(:,2);
    
    kmax = size(XYdata,1);
    
    DT = delaunayTriangulation(XYdata);      %Delaunay triangulation
    [Vvect,C] = DT.voronoiDiagram;          %Voronoi Diagram (points become polygons):
    %(V: vertex coordinates)
    %(C: vertex list for each points~polygons)
    VA = vertexAttachments(DT);         %List of vertices for each points (minus the vertex at infinity... almost like C)
    EE = edges(DT);                     %List of edges in triangulation (~neighbourhood connectivity of polygons)
    
    %area list (area for each polygons, except infinite polygons for which the entry is 0)
    Avect = zeros(length(DT.Points),1);
    for k=1:length(DT.Points)
        if prod(C{k}~=1)==0
            Avect(k) = 1;
        else
            Avect(k) = area(polyshape(Vvect(C{k},1),Vvect(C{k},2)));
        end
    end
    
    %length matrix (lenght of common edge between two polygons, zero if
    %not touching, or when at least one of the polygons is infinite)
    L = zeros(length(DT.Points));
    for k=1:size(EE,1)
        triInds = VA{EE(k,1)}(find(sum(VA{EE(k,1)} == VA{EE(k,2)}')))+1;
        if (length(triInds) == 2) && (isInside(EE(k,1))>0) && (isInside(EE(k,2))>0)
            L(EE(k,1),EE(k,2)) = norm(Vvect(triInds(1),:)-Vvect(triInds(2),:));
        end
    end
    L = L + L';
    L = ((Avect~=0)*(Avect~=0)').*L; %for Neumann boundary condition
    Mmatrix = diag(1./Avect)*(L-diag(sum(L,2)));
    % Backward matrices
    disp('Invert matrices...')
    tic
    Mu = inv( eye(kmax)-Mmatrix*dt*(params.Du*params.P/params.epsilon) );
    toc
    Mv = inv( eye(kmax)-Mmatrix*dt*(params.Dv*params.P/params.epsilon) );
    toc
    Mw = inv( eye(kmax)-Mmatrix*dt*(params.Dw*params.P/params.epsilon) );
    toc
    
    
%     % Initial conditions: juvenile pattern
%     r=0.5; % ~juvenile colour difference intensity parameter
%     % homogeneous unstable steady state
%     ffunc = @(x)(Fu(x(1),x(2),x(3),params).^2 + Fv(x(1),x(2),x(3),params).^2 + Fw(x(1),x(2),x(3),params).^2);
%     [hss_sol,~] = fminsearch(ffunc,[3 3 3]);
%     UVWanalytic_hss = hss_sol;
%     % main clusters
%     ffunc = @(x)(objfunc1stLevel(x,params));
%     [h1st_sol,~] = fminsearch(ffunc,[6 1 1 1 10 1]);
%     UVWanalyticB_1st = h1st_sol(1:3);
%     UVWanalyticG_1st = h1st_sol(4:6);
%     UVW_initB = UVWanalytic_hss+r*(UVWanalyticB_1st-UVWanalytic_hss);
%     UVW_initG = UVWanalytic_hss+r*(UVWanalyticG_1st-UVWanalytic_hss);
%     %initial conditions
%     U = (juvstatesGreen==1)*UVW_initG(1)+(juvstatesGreen==0)*UVW_initB(1);
%     V = (juvstatesGreen==1)*UVW_initG(2)+(juvstatesGreen==0)*UVW_initB(2);
%     W = (juvstatesGreen==1)*UVW_initG(3)+(juvstatesGreen==0)*UVW_initB(3);
    

    % Initial conditions: homogeneous unstable steady state
    ffunc = @(x)(Fu(x(1),x(2),x(3),params).^2 + Fv(x(1),x(2),x(3),params).^2 + Fw(x(1),x(2),x(3),params).^2);
    [hss_sol,~] = fminsearch(ffunc,[3 3 3]);
    % initial conditions
    U = zeros(length(xcenter),1);
    V = U;
    W = U;
    for k=1:length(xcenter)
        U(k) = hss_sol(1)+0.01*(2*rand()-1);
        V(k) = hss_sol(2)+0.01*(2*rand()-1);
        W(k) = hss_sol(3)+0.01*(2*rand()-1);
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else        %for perturbed hex. voronoi
    
    % strength of noise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r0val = 0.2;
    filestring = 'r0p2'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=100;
    M=60;
    S = params.S;
    
    % build hex. lattice
    sqrt3 = 3^(0.5);
    xcenter = [];
    ycenter = [];
    for n=0:N-1
        for m=0:M-1
            c = n*M+m+1;
            xcenter(c) = (m)*(3/2*S);
            ycenter(c) = -(n)*(sqrt3*S)+mod(m,2)*(sqrt3/2*S);
            
        end
    end
    intermFactor = 0.975; % first border layer
    insideFactor = 0.925; % second border layer
    xMiddle = mean(xcenter);
    yMiddle = mean(ycenter);
    
    % label border layers
    for n=0:N-1
        for m=0:M-1      
            c = n*M+m+1;
            isInside(c) = 0;
            if abs(xcenter(c)-xMiddle)<intermFactor*abs(xMiddle) && abs(ycenter(c)-yMiddle)<intermFactor*abs(yMiddle)
                isInside(c) = 1;
                % these will be recomputed...
            end
            if abs(xcenter(c)-xMiddle)<insideFactor*abs(xMiddle) && abs(ycenter(c)-yMiddle)<insideFactor*abs(yMiddle)
                isInside(c) = 2;
                % these will be moved and recomputed...
            end
        end
    end
    % initial L matrix (before perturbing central region, with periodicity)
    L = zeros(N*M,N*M);
    for n=0:N-1
        for m=0:M-1
            %for periodicity...
            nm1 = mod(n-1,N);
            np1 = mod(n+1,N);
            mm1 = mod((m-1),M);
            mp1 = mod((m+1),M);
            %
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

    
    % perturb positions without changing topology (nearest neighbour number)
    r0 = r0val*S;
    %
    xcenterR = xcenter;
    ycenterR = ycenter;
    redo = 1:size(xcenter,2);
    while ~isempty(redo)
        xcenterR(redo) = xcenter(redo) + (isInside(redo)>1).*normrnd(0,r0,1,length(redo));
        ycenterR(redo) = ycenter(redo) + (isInside(redo)>1).*normrnd(0,r0,1,length(redo));
        [Vv,C] = voronoin([xcenterR' ycenterR']);
        TRI = delaunay(xcenterR', ycenterR');
        pairsDoubled = [[TRI(:,1);TRI(:,2);TRI(:,3)] [TRI(:,2);TRI(:,3);TRI(:,1)]];
        As = sparse(pairsDoubled(:,1),pairsDoubled(:,2),0*pairsDoubled(:,1)+1,N*M,N*M);
        nns = full(sum(As,2));
        redo = find( (nns ~= 6) & (isInside' > 1) ); % all inside points that have not the good num. of n.n.
        %length(redo)
    end
    
    % matrices for backward Euler integration
    % Lmatrix perturbed
    disp('Building matrices...')
    AreaHex = 3*3^(0.5)*S^2/2;
    Avect = AreaHex*ones(N*M,1);
    for n=0:N-1
        for m=0:M-1 
            c = n*M+m+1;
            if isInside(c) > 0              
                neighbourlabels = unique(reshape(TRI(find(sum(TRI==c,2)),:),[],1)); %with self
                neighbourlabels = neighbourlabels(neighbourlabels ~= c);              
                for np=0:N-1
                    for mp=0:M-1
                        cp = np*M+mp+1;
                        if sum(cp==neighbourlabels)==1
                            edgePoints = Vv(intersect(C{c},C{cp}),:);
                            L(c,cp) = norm(edgePoints(1,:)-edgePoints(2,:));
                            %A(c,cp) = 1;
                        else
                            L(c,cp) = 0;
                            %A(c,cp) = 0;
                        end
                    end
                end
                % Area vector
                Avect(c) = area(polyshape(Vv(C{c},:)));
            end
        end
    end
    disp( ['number of non-hexagons:' num2str(sum(sum(L>0)~=6))] )
    L=(L+L')/2;  %in case of numerical asymetries
    Mmatrix = diag(1./Avect)*(L-diag(sum(L,2)));
    % Backward matrices
    disp('Invert matrices...')
    tic
    Mu = inv( eye(N*M)-Mmatrix*dt*(params.Du*params.P/params.epsilon) );
    toc
    Mv = inv( eye(N*M)-Mmatrix*dt*(params.Dv*params.P/params.epsilon) );
    toc
    Mw = inv( eye(N*M)-Mmatrix*dt*(params.Dw*params.P/params.epsilon) );
    toc
    
    % Initial conditions: homogeneous unstable steady state
    ffunc = @(x)(Fu(x(1),x(2),x(3),params).^2 + Fv(x(1),x(2),x(3),params).^2 + Fw(x(1),x(2),x(3),params).^2);
    [hss_sol,~] = fminsearch(ffunc,[3 3 3]);
    % initial conditions
    U = zeros(N*M,1);
    V = U;
    W = U;
    for n=0:N-1
        for m=0:M-1
            U(n*M+m+1) = hss_sol(1)+0.01*(2*rand()-1);
            V(n*M+m+1) = hss_sol(2)+0.01*(2*rand()-1);
            W(n*M+m+1) = hss_sol(3)+0.01*(2*rand()-1);
        end
    end
    
end




%% Simulation

disp('Simulation...')
% run the simulation
tsteps = 3000;
for step=1:tsteps
    disp(step)
    Unew = Mu*(U+dt*Fu(U,V,W,params));
    Vnew = Mv*(V+dt*Fv(U,V,W,params));
    Wnew = Mw*(W+dt*Fw(U,V,W,params));
    diff = abs(Unew-U)+abs(Vnew-V)+abs(Wnew-W);
    sum(diff)      
    U = Unew;
    V = Vnew;
    W = Wnew;
    output.Uhistory(:,step) = single(U);
    output.Vhistory(:,step) = single(V);
    output.Whistory(:,step) = single(W);
    output.thistory(step) = step*dt;
end

output.L = L;
save([['outputVoronoi_' filestring '_'],datestr(now(),'yyyymmdd_HHMMSS'),'.mat'],'output','-v7.3');



%%

U = output.Uhistory(:,end)
V = output.Vhistory(:,end)
nG = (L>0)*(output.Vhistory(:,end)>4); % number of green neighbours
figure
scatter(U(isInside==2),V(isInside==2),'.') %only look at those inside the region 2 (central region)
hold on
scatter(U( (isInside'==2) & (nG==2) ),V( (isInside'==2) & (nG==2) ),'o') % for scales with nG=2 green neighbours



%%

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





