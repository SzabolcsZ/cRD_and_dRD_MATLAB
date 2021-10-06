function vtkwrite( filename,dataType,varargin )
% VTKWRITE Writes 3D Matlab array into VTK file format.
%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w) writes
%  a structured 3D vector data into VTK file, with name specified by the string
%  filename. (u,v,w) are the vector components at the points (x,y,z). x,y,z
%  should be 3-D matrices like those generated by meshgrid, where
%  point(ijk) is specified by x(i,j,k), y(i,j,k) and z(i,j,k).
%  The matrices x,y,z,u,v,w must all be the same size and contain
%  corrresponding position and vector component. The string title specifies
%  the name of the vector field to be saved. 
%
%  vtkwrite(filename,'structured_grid',x,y,z,'scalars',title,r) writes a 3D
%  scalar data into VTK file whose name is specified by the string
%  filename. r is the scalar value at the points (x,y,z). The matrices
%  x,y,z,r must all be the same size and contain the corresponding position
%  and scalar values. 
%
%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w,'scalars',
%  title2,r) writes a 3D structured grid that contains both vector and scalar values.
%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%
%  vtkwrite(filename, 'structured_points', title, m) saves matrix m (could
%  be 1D, 2D or 3D array) into vtk as structured points.
%
%  vtkwrite(filename, 'structured_points', title, m, 'spacing', sx, sy, sz)
%  allows user to specify spacing. (default: 1, 1, 1). This is the aspect
%  ratio of a single voxel. 
%
%  vtkwrite(filename, 'structured_points', title, m, 'origin', ox, oy, oz)
%  allows user to speicify origin of dataset. (default: 0, 0, 0).
%
%  vtkwrite(filename,'unstructured_grid',x,y,z,'vectors',title,u,v,w,'scalars',
%  title2,r) writes a 3D unstructured grid that contains both vector and scalar values.
%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%  
%  vtkwrite(filename, 'polydata', 'lines', x, y, z) exports a 3D line where
%  x,y,z are coordinates of the points that make the line. x, y, z are
%  vectors containing the coordinates of points of the line, where point(n)
%  is specified by x(n), y(n) and z(n).
%
%  vtkwrite(filename,'polydata','lines',x,y,z,'Precision',n) allows you to
%  specify precision of the exported number up to n digits after decimal
%  point. Default precision is 3 digits. 
%
%  vtkwrite(filename,'polydata','triangle',x,y,z,tri) exports a list of
%  triangles where x,y,z are the coordinates of the points and tri is an
%  m*3 matrix whose rows denote the points of the individual triangles.
%
%  vtkwrite(filename,'polydata','tetrahedron',x,y,z,tetra) exports a list
%  of tetrahedrons where x,y,z are the coordinates of the points
%  and tetra is an m*4 matrix whose rows denote the points of individual
%  tetrahedrons. 
%  
%  vtkwrite('execute','polydata','lines',x,y,z) will save data with default
%  filename ''matlab_export.vtk' and automatically loads data into
%  ParaView. 
%  
%  Version 2.3
%  Copyright, Chaoyuan Yeh, 2016
%  Codes are modified from William Thielicke and David Gingras's submission.    

if strcmpi(filename,'execute'), filename = 'matlab_export.vtk'; end
fid = fopen(filename, 'w'); 
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');


binaryflag = any(strcmpi(varargin, 'BINARY'));
if any(strcmpi(varargin, 'PRECISION'))
    precision = num2str(varargin{find(strcmpi(varargin, 'PRECISION'))+1});
else
    precision = '4';
end

switch upper(dataType)
    case 'STRUCTURED_POINTS'
        container = struct;
        if isstruct(varargin{1})
            container = varargin{1};
        else
            container.(varargin{1}) = varargin{2};
        end
        title = fieldnames(container);
        if any(strcmpi(varargin, 'spacing'))
            sx = varargin{find(strcmpi(varargin, 'spacing'))+1};
            sy = varargin{find(strcmpi(varargin, 'spacing'))+2};
            sz = varargin{find(strcmpi(varargin, 'spacing'))+3};
        else
            sx = 1;
            sy = 1;
            sz = 1;
        end
        if any(strcmpi(varargin, 'origin'))
            ox = varargin{find(strcmpi(varargin, 'origin'))+1};
            oy = varargin{find(strcmpi(varargin, 'origin'))+2};
            oz = varargin{find(strcmpi(varargin, 'origin'))+3};
        else
            ox = 0;
            oy = 0;
            oz = 0;
        end
        [nx, ny, nz] = size(container.(title{1}));
        setdataformat(fid, binaryflag);
        
        fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
        fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
        fprintf(fid, ['SPACING ', num2str(sx), ' ', num2str(sy), ' ', num2str(sz), '\n']);
        fprintf(fid, ['ORIGIN ', num2str(ox), ' ', num2str(oy), ' ', num2str(oz), '\n']);
        fprintf(fid, 'POINT_DATA %d\n', nx*ny*nz);
        for i=1:numel(title)
            data_type = class(container.(title{i}));
            switch class(container.(title{i}))
                case 'logical'
                    data_type = 'bit';
                case 'uint8'
                    data_type = 'unsigned_char';
                case 'uint16'
                    data_type = 'unsigned_short';
                case 'uint32'
                    data_type = 'unsigned_int';
                case 'uint64'
                    data_type = 'unsigned_long';
                case 'int8'
                    data_type = 'char';
                case 'int16'
                    data_type = 'short';
                case 'int32'
                    data_type = 'int';
                case 'int64'
                    data_type = 'long';                
                case 'single'
                    data_type = 'float';        
            end
            fprintf(fid, ['SCALARS ', title{i}, ' ',data_type,' 1\n']);
            fprintf(fid,'LOOKUP_TABLE default\n');
            if ~binaryflag
                spec = '%d ';
                if strcmp(data_type,'float')
                    spec = ['%0.', precision, 'f '];
                end
                fprintf(fid, spec, container.(title{i})(:)');

            else
                fwrite(fid, container.(title{i})(:)', class(container.(title{i})), 'b');
            end
        end
        
    case {'STRUCTURED_GRID','UNSTRUCTURED_GRID'}
        % 3. The format data proper is saved in (ASCII or Binary). Use
        % fprintf to write data in the case of ASCII and fwrite for binary.
%        if numel(varargin)<6, error('Not enough input arguments'); end
        setdataformat(fid, binaryflag);

        container = struct;
        if isstruct(varargin{1})
            container = varargin{1};
        else
            container.(varargin{1}) = varargin{2};
        end
        
%         fprintf(fid, 'BINARY\n');
        x = container.V(:,1);
        y = container.V(:,2);
        z = container.V(:,3);
        if sum(size(x)==size(y) & size(y)==size(z))~=length(size(x))
            error('Input dimensions do not match')
        end
        n_elements = numel(x);
        % 4. Type of Dataset ( can be STRUCTURED_POINTS, STRUCTURED_GRID,
        % UNSTRUCTURED_GRID, POLYDATA, RECTILINEAR_GRID or FIELD )
        % This part, dataset structure, begins with a line containing the
        % keyword 'DATASET' followed by a keyword describing the type of dataset.
        % Then the geomettry part describes geometry and topology of the dataset.
        if strcmpi(dataType,'STRUCTURED_GRID')
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
            fprintf(fid, 'DIMENSIONS %d %d %d\n', size(x,1), size(x,2), size(x,3));
        else
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
        end
        fprintf(fid, ['POINTS ' num2str(n_elements) ' float\n']);
        output = [x(:)'; y(:)'; z(:)'];
        
        if ~binaryflag
            spec = ['%0.', precision, 'f '];
            %fprintf(fid, spec, output);
            fprintf(fid, [spec spec spec ' \n'], output);
        else
            fwrite(fid, output, 'float', 'b');
        end
        
        cidx = find(strcmpi(varargin,'CELLS'));
        if cidx~=0
            switch upper(varargin{cidx+1})         
                case 'TRIANGLE'
                    ntri = length(varargin{cidx+2});
                    fprintf(fid,'\nCELLS %d %d\n',ntri,4*ntri);
                    if ~binaryflag
                        fprintf(fid,'3 %d %d %d\n',(varargin{cidx+2}-1)');
                    else
                        fwrite(fid,[ones(1,ntri)*3;(varargin{cidx+2}-1)'],'int','b');
                    end                    
                    fprintf(fid,'\nCELL_TYPES %d\n',ntri);
                    if ~binaryflag
                        fprintf(fid,'%d\n',5*ones(1,ntri));
                    else
                        fwrite(fid,ones(1,ntri)*5,'int','b');
                    end

                % added by S.ZAKANY
                case 'TETRAHEDRON'
                    ntri = length(varargin{cidx+2});
                    fprintf(fid,'\nCELLS %d %d\n',ntri,5*ntri);
                    if ~binaryflag
                        fprintf(fid,'4 %d %d %d %d\n',(varargin{cidx+2}-1)');
                    else
                        fwrite(fid,[ones(1,ntri)*4;(varargin{cidx+2}-1)'],'int','b');
                    end                    
                    fprintf(fid,'\nCELL_TYPES %d\n',ntri);
                % end added
                    if ~binaryflag
                        fprintf(fid,'%d\n',10*ones(1,ntri));
                    else
                        fwrite(fid,ones(1,ntri)*10,'int','b');
                    end
                    
            end
        end        
        
        % 5.This final part describe the dataset attributes and begins with the
        % keywords 'POINT_DATA' or 'CELL_DATA', followed by an integer number
        % specifying the number of points of cells. Other keyword/data combination
        % then define the actual dataset attribute values.
        %fprintf(fid, ['\nCELL_DATA ' num2str(ntri)]);
        fprintf(fid, ['\nPOINT_DATA ' num2str(n_elements)]);
        % Parse remaining argument.
        vidx = find(strcmpi(varargin,'VECTORS'));
        sidx = find(strcmpi(varargin,'SCALARS'));
        if vidx~=0
            for ii = 1:length(vidx)
                title = varargin{vidx(ii)+1};
                % Data enteries begin with a keyword specifying data type
                % and numeric format.
                fprintf(fid, ['\nVECTORS ', title,' float\n']);
                output = [varargin{ vidx(ii) + 2 }(:)';...
                          varargin{ vidx(ii) + 3 }(:)';...
                          varargin{ vidx(ii) + 4 }(:)'];

                if ~binaryflag
                    spec = ['%0.', precision, 'f '];
                    fprintf(fid, spec, output);
                else
                    fwrite(fid, output, 'float', 'b');
                end
%                 fwrite(fid, [reshape(varargin{vidx(ii)+2},1,n_elements);...
%                 reshape(varargin{vidx(ii)+3},1,n_elements);...
%                 reshape(varargin{vidx(ii)+4},1,n_elements)],'float','b');
            end
        end
        if sidx~=0
            for ii = 1:length(sidx)
                title = varargin{sidx(ii)+1};
                fprintf(fid, ['\nSCALARS ', title,' float 1\n']);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                if ~binaryflag
                    spec = ['%0.', precision, 'f '];
                    fprintf(fid, spec, varargin{ sidx(ii) + 2});
                else
                    fwrite(fid, varargin{ sidx(ii) + 2}, 'float', 'b');
                end
%                 fwrite(fid, reshape(varargin{sidx(ii)+2},1,n_elements),'float','b');
            end
        end
        
    case 'POLYDATA'

        fprintf(fid, 'ASCII\n');
        if numel(varargin)<4, error('Not enough input arguments'); end
        x = varargin{2}(:);
        y = varargin{3}(:);
        z = varargin{4}(:);
        if numel(varargin)<4, error('Not enough input arguments'); end
        if sum(size(x)==size(y) & size(y)==size(z))~= length(size(x))
            error('Input dimesions do not match')
        end
        n_elements = numel(x);
        fprintf(fid, 'DATASET POLYDATA\n');
        if mod(n_elements,3)==1
            x(n_elements+1:n_elements+2,1)=[0;0];
            y(n_elements+1:n_elements+2,1)=[0;0];
            z(n_elements+1:n_elements+2,1)=[0;0];
        elseif mod(n_elements,3)==2
            x(n_elements+1,1)=0;
            y(n_elements+1,1)=0;
            z(n_elements+1,1)=0;
        end
        nbpoint = numel(x);
        fprintf(fid, ['POINTS ' num2str(nbpoint) ' float\n']);
        
        spec = [repmat(['%0.', precision, 'f '], 1, 9), '\n'];
        
        output = [x(1:3:end-2), y(1:3:end-2), z(1:3:end-2),...
                  x(2:3:end-1), y(2:3:end-1), z(2:3:end-1),...
                  x(3:3:end), y(3:3:end), z(3:3:end)]';
              
        fprintf(fid, spec, output);
        
        switch upper(varargin{1})
            case 'LINES'
                if mod(n_elements,2)==0
                    nbLine = 2*n_elements-2;
                else
                    nbLine = 2*(n_elements-1);
                end
                conn1 = zeros(nbLine,1);
                conn2 = zeros(nbLine,1);
                conn2(1:nbLine/2) = 1:nbLine/2;
                conn1(1:nbLine/2) = conn2(1:nbLine/2)-1;
                conn1(nbLine/2+1:end) = 1:nbLine/2;
                conn2(nbLine/2+1:end) = conn1(nbLine/2+1:end)-1;
                fprintf(fid,'\nLINES %d %d\n',nbLine,3*nbLine);
                fprintf(fid,'2 %d %d\n',[conn1';conn2']);
            case 'TRIANGLE'
                ntri = length(varargin{5});
                fprintf(fid,'\nPOLYGONS %d %d\n',ntri,4*ntri);
                fprintf(fid,'3 %d %d %d\n',(varargin{5}-1)');
            case 'TETRAHEDRON'
                ntetra = length(varargin{5});
                fprintf(fid,'\nPOLYGONS %d %d\n',ntetra,5*ntetra);
                fprintf(fid,'4 %d %d %d %d\n',(varargin{5}-1)');
        end     
end
fclose(fid);
if strcmpi(filename,'matlab_export.vtk')
    switch computer
        case {'PCWIN','PCWIN64'}
            !paraview.exe --data='matlab_export.vtk' &
            % Exclamation point character is a shell escape, the rest of the
            % input line will be sent to operating system. It can not take
            % variables, though. The & at the end of line will return control to 
            % Matlab even when the outside process is still running. 
        case {'GLNXA64','MACI64'}
            !paraview --data='matlab_export.vtk' &
    end
end
end

function setdataformat(fid, binaryflag)

if ~binaryflag
    fprintf(fid, 'ASCII\n');
else
    fprintf(fid, 'BINARY\n');
end
end
    