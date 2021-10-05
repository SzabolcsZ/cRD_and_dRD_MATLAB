%
% Continuous reaction-diffusion on 3D domain made of bell-shaped scales
%
% Szabolcs Zakany
% LANE, University of Geneva
% 2021

%%

gpuDevice(1)

% import 3D mesh of lattice of bell-shaped scales (generated with lizardskinLattice.m)
tic
mesh = getMeshStruct3D([pwd,'/lizardSkin.vtk']);
%mesh = getMeshStruct3D([pwd,'/lizardSkin100.msh']);
%mesh = getMeshStruct3D([pwd,'/lizardSkin400.msh']);
toc

save_all = 0; % 1: saves intermediate states, 0: overwrites...
filename = 'lizRD';
rng('shuffle');


%% Build matrices F00 (=A) ,F11 (=B)

tic

% Test building F00
disp('Computing F00...')
c=1;
clear iSP
clear jSP
clear vSP
for i=1:mesh.N_v

    %i=j
    iSP(c) = i;
    jSP(c) = i;
    vSP(c) = 1/10*sum(mesh.Vol(mesh.K{i}));
    c = c+1;

    %i~=j, A_ij=1;
    %(we use the fact that values with repeated (i,j,..) are added together when building sparse matrix)
    Tj = mesh.T(mesh.K{i},:);
    jSPnext = Tj(Tj~=i);
    iSPnext = repmat(i,length(jSPnext),1);
    kvals = repmat(mesh.K{i},1,4);
    vSPnext = 1/20*mesh.Vol(kvals(Tj~=i));
    cnext = length(jSPnext);
    
    jSP((c):(c+cnext-1))=jSPnext;
    iSP((c):(c+cnext-1))=iSPnext;
    vSP((c):(c+cnext-1))=vSPnext;
    c=c+cnext;
    
end
F00 = sparse(iSP,jSP,vSP);
F00 = 1/2*(F00+F00'); %maybe reduces numerical errors...
toc



%normal vectors for each tetraedra, ordered according to the opposite
%vertex labelled in T(k,:)
N(1,:,:) = 1/2*cross(mesh.V(mesh.T(:,3),:)-mesh.V(mesh.T(:,2),:),mesh.V(mesh.T(:,4),:)-mesh.V(mesh.T(:,2),:));
N(2,:,:) = 1/2*cross(mesh.V(mesh.T(:,1),:)-mesh.V(mesh.T(:,3),:),mesh.V(mesh.T(:,4),:)-mesh.V(mesh.T(:,3),:));
N(3,:,:) = 1/2*cross(mesh.V(mesh.T(:,1),:)-mesh.V(mesh.T(:,4),:),mesh.V(mesh.T(:,2),:)-mesh.V(mesh.T(:,4),:));
N(4,:,:) = 1/2*cross(mesh.V(mesh.T(:,3),:)-mesh.V(mesh.T(:,1),:),mesh.V(mesh.T(:,2),:)-mesh.V(mesh.T(:,1),:));

tic

% Test building F11
disp('Computing F11...')
c=1;
clear iSP
clear jSP
clear vSP
for i=1:mesh.N_v
    
    
    
    %i~=j, A_ij=1;
    %(we use the fact that values with repeated (i,j,..) are added together when building sparse matrix)
    Tj=mesh.T(mesh.K{i},:)';
    jSPnext = Tj(Tj~=i);
    iSPnext = repmat(i,length(jSPnext),1);
    knext = reshape(repmat(mesh.K{i},1,3)',[],1);
    labels = repmat(1:4,size(Tj,2),1)';
    labelsInext = reshape(repmat(labels(Tj==i),1,3)',[],1);
    labelsJnext = labels(Tj~=i);
    volnext = mesh.Vol(knext);
    
    Nnext = N(:,knext,:);
    N1next=Nnext(:,:,1);
    N2next=Nnext(:,:,2);
    N3next=Nnext(:,:,3);
    
    Nveci = [N1next(sub2ind(size(N1next),labelsInext',1:length(labelsInext)));
        N2next(sub2ind(size(N2next),labelsInext',1:length(labelsInext)));
        N3next(sub2ind(size(N3next),labelsInext',1:length(labelsInext)))];
    Nvecj = [N1next(sub2ind(size(N1next),labelsJnext',1:length(labelsJnext)));
        N2next(sub2ind(size(N2next),labelsJnext',1:length(labelsJnext)));
        N3next(sub2ind(size(N3next),labelsJnext',1:length(labelsJnext)))];
    
   
    Pijknext = sum(Nveci.*Nvecj)';
    vSPnext = 1/9 * Pijknext./volnext;    %here we should add average of D...
    cnext = length(jSPnext);
    
    jSP((c):(c+cnext-1))=jSPnext;
    iSP((c):(c+cnext-1))=iSPnext;
    vSP((c):(c+cnext-1))=vSPnext;
    c=c+cnext;

    %i=j
    iSP(c) = i;
    jSP(c) = i;
    vSP(c) =  1/9 * 1/3*sum(sum(Nveci.*Nveci)'./volnext);  %(1/3 for overcounting) here we should add average of D...
    c = c+1;
       
end
F11 = sparse(iSP,jSP,vSP);
F11 = 1/2*(F11+F11'); %maybe reduces numerical errors...

toc

%% selected locations on mesh to keep scale history

% this uses the ordering of the mesh points... (all points with same X and Y are consecutive...)
Itop = (mesh.V(2:end-1,3)>mesh.V(1:end-2,3)) & (mesh.V(2:end-1,3)>mesh.V(3:end,3));
Itop = [0; Itop; 1];
selectedVlabels = (mesh.V(:,3)>(max(mesh.V(:,3))-4)) & Itop;
selectedVxyz = [mesh.V(selectedVlabels,1) mesh.V(selectedVlabels,2) mesh.V(selectedVlabels,3)];


%% start simulation

disp('Setting initial conditions')
tic
% homogeneous steady-state values:
baseu1 = 3.51;
basev1 = 2.96;
basew1 = 3.45;
perturbation = 0.01;
u1pde = gpuArray( (baseu1).*ones(mesh.N_v, 1) + perturbation*randFunc2(mesh));
v1pde = gpuArray( (basev1).*ones(mesh.N_v, 1) + perturbation*randFunc2(mesh));
w1pde = gpuArray( (basew1).*ones(mesh.N_v, 1) + perturbation*randFunc2(mesh));
dw1pde = w1pde;
toc

dwsum = [1];
tarray = [0];
usHist = gather(u1pde(selectedVlabels));
vsHist = gather(v1pde(selectedVlabels));
wsHist = gather(w1pde(selectedVlabels));

Du1 = 0.189;
Dv1 = Du1;
Dw1 = 13.5;

nt = 10000;
dt = 5.0;


% Main simulation loop
disp('Building GPU matrices')
tic
MatU1 = F00 + dt*Du1*F11;
MatV1 = F00 + dt*Dv1*F11;
MatW1 = F00 + dt*Dw1*F11;
gpuMatU1 = gpuArray(MatU1);
gpuMatV1 = gpuArray(MatV1);
gpuMatW1 = gpuArray(MatW1);
%diagonal preconditioners:
gpuPREu = gpuArray(spdiags(diag(MatU1),0,size(MatU1,1),size(MatU1,2)));
gpuPREv = gpuArray(spdiags(diag(MatV1),0,size(MatV1,1),size(MatV1,2)));
gpuPREw = gpuArray(spdiags(diag(MatW1),0,size(MatW1,1),size(MatW1,2)));
toc
disp('Starting simulation')
tic
ii = 0;
for i = 1:nt
    
    % Render and write
    if mod(i, 50) == 1
        ii = ii+1;
        if save_all
            vtkfile = [filename '_' num2str(ii) '.vtk'];
        else
            vtkfile = [filename '.vtk'];            
        end
    
        disp(['Writing vtk for timestep ' num2str(i)] )
        vtkwrite(vtkfile,'UNSTRUCTURED_GRID',...
            mesh,'CELLS','TETRAHEDRON',mesh.T,...
            'SCALARS','u1',gather(u1pde),...
            'SCALARS','v1',gather(v1pde),...
            'SCALARS','w1',gather(w1pde),...
            'SCALARS','dw1',gather(dw1pde));        
        disp('Finished vtk')
        save('res.mat','u1pde','v1pde','w1pde','tarray','dwsum');
        save('selectedVhistory.mat','tarray','dwsum','usHist','vsHist','wsHist','selectedVxyz');
    end
    
    % Set reaction vectors
    FGHval = FGH(u1pde,v1pde,w1pde);
    % Solve RD equation (implicit)
    Mlu1 = F00*(u1pde+dt*FGHval(:,1));
    Mlv1 = F00*(v1pde+dt*FGHval(:,2));
    Mlw1 = F00*(w1pde+dt*FGHval(:,3));
    [u1pdeNEW, flagu, ~, iteru] = bicgstab(gpuMatU1, Mlu1,1e-5,3000,gpuPREu);
    [v1pdeNEW, flagv, ~, iterv] = bicgstab(gpuMatV1, Mlv1,1e-5,3000,gpuPREv);
    [w1pdeNEW, flagw, ~, iterw] = bicgstab(gpuMatW1, Mlw1,1e-5,3000,gpuPREw);
    disp(['Iter. ' num2str(i)])
    disp([flagu iteru; flagv iterv; flagw iterw])
    dw1pdeNEW = abs(w1pdeNEW-w1pde);
    if i>200 && log(gather(sum(dw1pdeNEW)))/log(gather(sum(dw1pde))) > 1.25
        [u1pdeNEW, flagu, ~, iteru] = bicgstab(gpuMatU1, Mlu1,1e-4,6000,gpuPREu);
        [v1pdeNEW, flagv, ~, iterv] = bicgstab(gpuMatV1, Mlv1,1e-4,6000,gpuPREv);
        [w1pdeNEW, flagw, ~, iterw] = bicgstab(gpuMatW1, Mlw1,1e-4,6000,gpuPREw);
        disp(['AGAIN Iter. ' num2str(i)])
        disp([flagu iteru; flagv iterv; flagw iterw])
        dw1pdeNEW = abs(w1pdeNEW-w1pde);
    end
    dw1pde = dw1pdeNEW;
    u1pde = u1pdeNEW;
    v1pde = v1pdeNEW;
    w1pde = w1pdeNEW;

    tarray = [tarray sum(tarray(end)+dt)];
    dwsum = [dwsum log(gather(sum(dw1pde))) ];
    usHist = [usHist gather(u1pde(selectedVlabels))];
    vsHist = [vsHist gather(v1pde(selectedVlabels))];
    wsHist = [wsHist gather(w1pde(selectedVlabels))];
    

end


save(['cRDsimulationResult_',datestr(now(),'yyyymmdd_HHMMSS'),'.mat'],...
    'tarray','dwsum','usHist','vsHist','wsHist','selectedVxyz');
toc


%% functions

% reaction function
function out = FGH(u,v,w)
    c1 = -0.04;
    c2 = -0.056;
    c3 = 0.382;
    c4 = -0.05;
    c5 = 0;
    c6 = 0.25;
    c7 = 0.016;
    c8 = -0.0325; %!
    c9 = 0.24;
    cu = 0.02;
    cv = 0.025;
    cw = 0.058; %!
    Fmax = 0.5;
    Gmax = 0.5;
    Hmax = 0.5;
    
    out = [max(0,min(Fmax,c1*v+c2*w+c3))-cu*u...
           max(0,min(Gmax,c4*u+c5*w+c6))-cv*v...
           max(0,min(Hmax,c7*u+c8*v+c9))-cw*w];
end


% function giving noise (1/f)
function out = randFunc2(mesh)
    
    Lx=max(mesh.V(:,1))-min(mesh.V(:,1));
    Ly=max(mesh.V(:,2))-min(mesh.V(:,2));
    out = zeros(size(mesh.V,1),1);
    disp('Building random function...')
    for n=1:200
        for t=1:5
            alpha = rand()*2*pi;
            out = out + (2*rand()-1)*1/sqrt(n)*sin( n*pi*(mesh.V(:,[1 2])*[cos(alpha) sin(alpha)]')/sqrt(Lx*Ly) );
        end
    end
    disp('...Finished.')
    out = 2*(out-min(out))./(max(out)-min(out))-1;
end



