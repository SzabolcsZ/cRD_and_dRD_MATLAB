%
% Create mesh for continuous reaction-diffusion 
% (3D domain made of bell-shaped scales)
%
% Szabolcs Zakany
% LANE, University of Geneva
% 2021
%%

% grid point numbers of mesh in X and Y directions:
Nx = 100;
Ny = 100;

%%
% Mesh of lizard skin

% scale lattice data
eps = 1;
Ni = 150; % (these should be enough for reasonable sized meshes)
Nj = 150;
S = 20/sqrt(3);
sigma = 6/6*S;
deltaZ = S/6; %height between bumps
Zheight = S;
[J,I] = meshgrid(-5:Nj,-5:Ni);
Xc = 3/2*S*J;
Yc = sqrt(3)*S*(I+1/2*mod(J-1,2) );

% mesh data
[Xgrid, Ygrid] = meshgrid(1:eps:Nx,1:eps:Ny);
ZgaussLin = gaussianLattice(Xgrid(:),Ygrid(:),Xc,Yc,sigma);
Zgauss = reshape(ZgaussLin,size(Xgrid,1),[]);
powerFactor = 1/3;
rFactor=3;
Zgrid = Zheight*rFactor*(Zgauss-min(Zgauss(:))) ./ ( rFactor*(Zgauss-min(Zgauss(:))) + 1 );
%Zgrid = (Zgauss-min(Zgauss(:)));

% figure
% surf(Xgrid,Ygrid,Zgrid)
% daspect([1 1 1])
%save('centers.mat','Xc','Yc')

disp('Building mesh')
Xall = [];
Yall = [];
Zall = [];
for k = 1:length(Zgrid(:))
    nn = max(1,ceil((deltaZ+Zgrid(k))/(eps)));
    Ztemp = (0:nn)/nn*(deltaZ+Zgrid(k));
    Xtemp = repmat(Xgrid(k),nn+1,1);
    Ytemp = repmat(Ygrid(k),nn+1,1);
    Zall= [Zall; Ztemp'];
    Xall= [Xall; Xtemp];
    Yall= [Yall; Ytemp];
end

disp('Building triangulation')
shp = alphaShape(Xall,Yall,Zall,1.5*eps);
tri = alphaTriangulation(shp);

figure
h = trimesh(tri,Xall,Yall,Zall,'FaceAlpha',0.3)
daspect([1 1 1])

mesh.V=[Xall Yall Zall];
vtkwrite('lizardSkin.vtk','UNSTRUCTURED_GRID',...
            mesh,'CELLS','TETRAHEDRON',tri,...
            'PRECISION',6);

% command to remesh using gmsh (https://gmsh.info/)        
command = 'gmsh lizardSkin.vtk -3 -o lizardSkin.msh';
system(command);


%%

function out = gaussianLattice(x,y,Xc,Yc,sigma)
    XcLin = Xc(:)';
    YcLin = Yc(:)';
    all = exp(-1/(2*sigma)*((x-XcLin).^2+(y-YcLin).^2) );
    out = sum(all,2);
end