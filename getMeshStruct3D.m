function mesh = getMeshStruct3D(file)

% mesh structure of 3D volume for finite element simulation code 

mesh.info = ['info: this information',newline,...
    'file: path of the .msh file (input)',newline,...
    'N_v: number of vertices',newline,...
    'V_id: vertex IDs (should be consecutive integers starting from 1)',newline,...
    'V: vertex coordinates',newline,...
    'N_t: number of volume elements (tetrahedra)',newline,...
    'T: labels of vertices in volume element (in right orientation)',newline,...
    'N_bf: nuber of border faces',newline,...
    'bF: labels of vertices in boundary face',newline,...
    ['K3: K3(a,:) set of rows [i j k t] with t the volume element label which is adjacent to the '...
    'face linking i and j and k in the positive direction'],newline,...
    'bK2: bK2(i,j) is the label of the border face adjacent to edge going from i to j (in right orientation), or 0',newline,...
    'A: adjacency matrix for vertices',newline,...
    'Ab: adjacency matrix for border faces',newline,...
    'J: labels of vertices sharing an edge with vertex, as cell',newline,...
    'K: labels of faces adjacent to vertex, as cell',newline,...
    'Vol: volumes of volume elements',newline,...
    'bNorm: unnormalized normals of border faces',newline,...
    'bnorm: normalized normals of border faces',newline,...
    'bAr: areas of border faces',newline,...
    ];


disp('  Reading mesh data...')


if strcmp(file(end-3:end),'.msh')
    % code to read .msh data, which can be produce by gmsh (https://gmsh.info)
    % no of nodes is mentioned in 5th row and first column
    N_v      = dlmread(file,'',[5-1 1-1 5-1 1-1]);          %number of vertices
    N_e      = dlmread(file,'',[7+N_v 0 7+N_v 0]);          %number of elements
    V_id     = dlmread(file,'',[5 0 4+N_v 0]);              %vertex ids
    V        = dlmread(file,'',[5 1 4+N_v 3]);               %vertex coordinates
    elements    = dlmread(file,'',[8+N_v 0 7+N_v+N_e 8]);   %all elements
    T = elements(elements(:,9)>0,6:9);                      %vertex labels for volumes (tetrahedra)
    N_t = size(T,1);                                        %number of faces
    bF = elements(elements(:,8)>0 & elements(:,9)==0,6:8);  %border faces (if structure available in .msh)
    N_bf = size(bF,1);                                      %number of border faces
end




if strcmp(file(end-3:end),'.vtk')
    % code to read .vtk data (e.g. produced by vtkwrite.m) for
    % UNSTRUCTURED_GRID format
    disp('  Assuming .vtk is UNSTRUCTURED_GRID and made of tetrahedron CELLS')
    fileID = fopen(file);
    fulltext = fileread(file);
    posPOINTS = strfind(fulltext,'POINTS');
    posCELLS = strfind(fulltext,'CELLS');
    posCELL_TYPES = strfind(fulltext,'CELL_TYPES');
    % extract points
    dpos = strfind(fulltext(posPOINTS+(0:100)),newline);
    text = fulltext((posPOINTS+dpos):(posCELLS-1));
    V = reshape( sscanf(text,'%f'),3,[] )';
    N_v = size(V,1);
    % extract triangles & tetrahedra
    % (assuming triangles are listed before tetrahedra)
    dpos3 = strfind(fulltext(posCELLS:(posCELL_TYPES-1)),[newline '3']);
    dpos4 = strfind(fulltext(posCELLS:(posCELL_TYPES-1)),[newline '4']);
    bF = [];
    T = [];
    if ~isempty(dpos3)
        text = fulltext( posCELLS+(dpos3(1):(dpos4(1)-2)) );
        temp = reshape( sscanf(text,'%f'),4,[] )';
        bF = temp(:,2:end)+1;
    end
    if ~isempty(dpos4)
        text = fulltext( (posCELLS+dpos4(1)):(posCELL_TYPES-1) );
        temp = reshape( sscanf(text,'%f'),5,[] )';
        T = temp(:,2:end)+1;
    end
    N_bf = size(bF,1);
    N_t = size(T,1);
    fclose(fileID);
    
end



% K3 list: K3 is a set of rows [i j k t]:
% t: the volume element label which is adjacent to the
% face linking i and j and k in the positive direction. 
disp(['  Computing adjacency data (' num2str(N_v) ' vertices)...'])
K3 = [T(:,1) T(:,2) T(:,3) (1:N_t)';
      T(:,2) T(:,1) T(:,4) (1:N_t)';
      T(:,3) T(:,4) T(:,1) (1:N_t)';
      T(:,4) T(:,3) T(:,2) (1:N_t)'];
% A adjacency matrix: symmetric matrix which is 1 when vertices share an edge
iSP = [T(:,1); T(:,1); T(:,1); T(:,2); T(:,2); T(:,3)];
jSP = [T(:,2); T(:,3); T(:,4); T(:,3); T(:,4); T(:,4)];
valSP = ones(6*N_t,1);
A = sparse(iSP,jSP,valSP,N_v,N_v);
A = (A + A')>0;




if N_bf>0
    % Ab adjacency matrix for border: symmetric matrix which is 1 when vertices share an edge
    iSP = [bF(:,1); bF(:,2); bF(:,3)];
    jSP = [bF(:,2); bF(:,3); bF(:,1)];
    valSP = [1:N_bf 1:N_bf 1:N_bf];
    %valSP = ones(3*N_bf,1);
    bK2 = sparse(iSP,jSP,valSP,N_v,N_v);
    bA = (bK2 + bK2')>0;

    % bNorm: non-normalized normal vectors for border faces
    % bAr: area values of border faces (triangles)
    % bnorm: normalized normals for border faces
    bNorm = cross(V(bF(:,2),:)-V(bF(:,1),:),V(bF(:,3),:)-V(bF(:,1),:) ) ;
    bAr =  1/2*vecnorm(bNorm,2,2);
    bnorm = bNorm ./ (2*bAr);
else
    bK2 = [];
    bA = [];
    bNorm = [];
    bAr = [];
    bnorm = [];
end



% K: list of faces labels adjacent to vertex.
% J: list of neighbour vertices to vertex
J = cell(N_v,1);
for i=1:N_v
    jset = find(A(:,i))';
    J{i} = sort( jset );
end
K = cell(N_v,1);
for k=1:N_t
    K{T(k,1)} = [K{T(k,1)};k];
    K{T(k,2)} = [K{T(k,2)};k];
    K{T(k,3)} = [K{T(k,3)};k];
    K{T(k,4)} = [K{T(k,4)};k];
end

% Vol: volume values of volume elements
Vol = 1/6 * sum( cross(V(T(:,2),:)-V(T(:,1),:),V(T(:,3),:)-V(T(:,1),:) ) .* ...
(V(T(:,4),:)-V(T(:,1),:)), 2 );


mesh.N_v = N_v;
%mesh.V_id = V_id;
mesh.V = V;
mesh.N_t = N_t;
mesh.T = T;
mesh.N_bf = N_bf;
mesh.bF = bF;
mesh.K3 = K3;
mesh.bK2 = bK2;
mesh.A = A;
mesh.bA = bA;
mesh.J = J;
mesh.K = K;
mesh.Vol = Vol;
mesh.bNorm = bNorm;
mesh.bnorm = bnorm;
mesh.bAr = bAr;




end