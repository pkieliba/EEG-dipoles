% Visualizes the EEG equivalent-dipole as an elipsoid in the MNI average
% MRI head space.
%
%
% INPUTS
%
% centroid: a structure with computed centroid (ield posxyz) and directional 
% standard deviation (field loc_std) of the elipsoid.
%
% OPTONAL INPUTS
%
% color:   a char or RGB array with the color used for plotting the
%               elipsoid
% proj:    ['on' | 'off'] project the elipsoid onto the 2D planes 
%           (dafault: on)
% projcolor:  a char or RGB array with the color used for plotting the
%               2D projections
% mri:        mri template file to use
% zoom:       ammount of zoom for the plotting
% nvert:      number of vertices of the ellipsoid (default: 30)
%
% EXAMPLE
%
% ellipsoidPlot(centroidIC21, 'color', 'r', 'proj', 'on')
%
function  ellipsoidPlot(centroid, varargin)
        
    if nargin < 1
        help ellipsoidPlot;
        return
    end

    g = finputcheck( varargin, { 'color'    {'string','real'}  []  [1 0 0];
                                 'view'	'real'	[]	[0 0 1];                        
                                 'proj'	'string'	{'on', 'off'}	'on';
                                 'projcolor'    {'real', 'string'}  []  [];  
                                 'mri'	{'string', 'struct'} []	'';
                                 'zoom'	'integer'	[0 Inf]	1500;
                                 'nvert'	'integer'	[]  30 });

    if ischar(g)
        error(g); 
    end   

    g.coordformat = 'MNI';

    % look up an MRI file if necessary
    if isempty(g.mri)
        dipfitdefs;
        g.mri = template_models(1).mrifile;
    end
    
    % decode color if neccessary
    if ischar(g.color)
        g.color = strcol2real(g.color);
    end
    if isempty(g.projcolor)
        g.projcolor = g.color;
    elseif ischar(g.projcolor)
        g.projcolor = strcol2real(g.projcolor);
    end
        
    % read anatomical MRI using Fieldtrip and SPM2 functions
    if ischar(g.mri)
        try 
            g.mri = load('-mat', g.mri);
            g.mri = g.mri.mri;
        catch
            disp('Failed to read Matlab file. Attempt to read MRI file using function ft_read_mri');
            try
                g.mri = ft_read_mri(g.mri);
                g.mri.anatomy = round(gammacorrection( g.mri.anatomy, 0.8));
                g.mri.anatomy = uint8(round(g.mri.anatomy/max(reshape(g.mri.anatomy, prod(g.mri.dim),1))*255));
            catch
                error('Cannot load file using ft_read_mri');
            end
        end
    end
    

    dat.sph2spm = []; 
    dat.imgs = g.mri.anatomy;
    dat.transform = g.mri.transform;    
    
    % MRI coordinates for slices
    if ~isfield(g.mri, 'xgrid')
        g.mri.xgrid = 1:size(dat.imgs,1); 
        g.mri.ygrid = 1:size(dat.imgs,2);
        g.mri.zgrid = 1:size(dat.imgs,3);
    end

    dat.imgcoords = {g.mri.xgrid g.mri.ygrid g.mri.zgrid};            

    % plot head graph in 3D
    indx = ceil(dat.imgcoords{1}(end)/2);
    indy = ceil(dat.imgcoords{2}(end)/2);
    indz = ceil(dat.imgcoords{3}(end)/2);

    [planex, planey, planez] = plotimgs(dat, [indx indy indz], dat.transform);
    set(gca, 'color', 'k');
    axis equal;
    set(gca, 'cameraviewanglemode', 'manual'); 
    camzoom(1.2^2);
    view(g.view);
    axis off;

    if centroid.loc_std(1)<1
        centroid.loc_std(1) = 1;
    end
    if centroid.loc_std(2)<1
        centroid.loc_std(2) = 1;
    end
    if centroid.loc_std(3)<1
        centroid.loc_std(3) = 1;
    end
     
    % plot elipsoid
    [xx, yy, zz] = ellipsoid(centroid.posxyz(1), centroid.posxyz(2), centroid.posxyz(3), 3*centroid.loc_std(1), 3*centroid.loc_std(2), 3*centroid.loc_std(3), g.nvert);
    l = sqrt(xx.*xx+yy.*yy+zz.*zz);
    normals = reshape([xx./l yy./l zz./l],[g.nvert+1 g.nvert+1 3]);
    colorarray = repmat(reshape(g.color , 1,1,3), [size(zz,1) size(zz,2) 1]);
    hold on
    h = surf(xx, yy, zz, colorarray, 'tag', 'tmpmov', 'EdgeColor', g.color , 'VertexNormals', normals, ...
        'backfacelighting', 'reverselit', 'facelighting', 'gouraud', 'facecolor', 'interp', 'SpecularExponent', 5, ...
        'AmbientStrength', 0.5, 'DiffuseStrength', 0.7, 'SpecularColorReflectance', 0.7, 'EdgeLighting', 'gouraud');
    
    % plot projections

    if ~isempty(g.proj)
        colorarray  = repmat(reshape(g.projcolor, 1,1,3), [size(zz,1) size(zz,2) 1]);
        h(end+1) = surf((planex+1) * ones(size(xx)), yy, zz, colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none','FaceAlpha',0.8); 
        h(end+1) = surf(xx, (planey-1)  * ones(size(yy)), zz, colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none','FaceAlpha',0.8);

        h(end+1) = surf(xx, yy, (planez+1) * ones(size(zz)), colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none','FaceAlpha',0.8); 
        
    end
    
    camlight left;
    camlight right;
    rotate3d on;
    box off;
    axis equal;
    axis off;
    set(gca,'Clipping','off');
    
return;

end

% electrode space to MRI space
function [x,y,z] = transform(x, y, z, transmat)
    
    if isempty(transmat)
        return
    end
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            tmparray = transmat * [x(i,j) y(i,j) z(i,j) 1]';
            x(i,j) = tmparray(1);
            y(i,j) = tmparray(2);
            z(i,j) = tmparray(3);
        end
    end
end

% plot images (transmat is the uniform matrix MRI coords -> elec coords)
function [planex, planey, planez] = plotimgs(dat, mricoord, transmat)
   
    % loading images
    if ndims(dat.imgs) == 4 
        img1(:,:,3) = rot90(squeeze(dat.imgs(mricoord(1),:,:,3))); 
        img2(:,:,3) = rot90(squeeze(dat.imgs(:,mricoord(2),:,3))); 
        img3(:,:,3) = rot90(squeeze(dat.imgs(:,:,mricoord(3),3))); 
        img1(:,:,2) = rot90(squeeze(dat.imgs(mricoord(1),:,:,2))); 
        img2(:,:,2) = rot90(squeeze(dat.imgs(:,mricoord(2),:,2))); 
        img3(:,:,2) = rot90(squeeze(dat.imgs(:,:,mricoord(3),2))); 
        img1(:,:,1) = rot90(squeeze(dat.imgs(mricoord(1),:,:,1))); 
        img2(:,:,1) = rot90(squeeze(dat.imgs(:,mricoord(2),:,1))); 
        img3(:,:,1) = rot90(squeeze(dat.imgs(:,:,mricoord(3),1)));     
    else
        img1 = rot90(squeeze(dat.imgs(mricoord(1),:,:))); 
        img2 = rot90(squeeze(dat.imgs(:,mricoord(2),:))); 
        img3 = rot90(squeeze(dat.imgs(:,:,mricoord(3)))); 

        if ndims(img1) == 2
            img1(:,:,3) = img1; 
            img1(:,:,2) = img1(:,:,1); 
        end
        if ndims(img2) == 2
            img2(:,:,3) = img2; 
            img2(:,:,2) = img2(:,:,1); 
        end
        if ndims(img3) == 2
            img3(:,:,3) = img3; 
            img3(:,:,2) = img3(:,:,1); 
        end
    end
    
    % computing coordinates for planes
    wy1 = [min(dat.imgcoords{2}) max(dat.imgcoords{2}); min(dat.imgcoords{2}) max(dat.imgcoords{2})];
    wz1 = [min(dat.imgcoords{3}) min(dat.imgcoords{3}); max(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wx2 = [min(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) max(dat.imgcoords{1})];
    wz2 = [min(dat.imgcoords{3}) min(dat.imgcoords{3}); max(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wx3 = [min(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) max(dat.imgcoords{1})];
    wy3 = [min(dat.imgcoords{2}) min(dat.imgcoords{2}); max(dat.imgcoords{2}) max(dat.imgcoords{2})];

    wx1 =  [1 1; 1 1] * dat.imgcoords{1}(1);
    wy2 =  [1 1; 1 1] * dat.imgcoords{2}(end);
    wz3 =  [1 1; 1 1] * dat.imgcoords{3}(1);
    
    % transform MRI coordinates to electrode space
    [elecwx1, elecwy1, elecwz1] = transform(wx1, wy1, wz1, transmat);
    [elecwx2, elecwy2, elecwz2] = transform(wx2, wy2, wz2, transmat);
    [elecwx3, elecwy3, elecwz3] = transform(wx3, wy3, wz3, transmat);
    
    % ploting surfaces
    options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                'direct','tag','img', 'facelighting', 'none' };
    hold on;
    surface(elecwx1, elecwy1, elecwz1, img1(end:-1:1,:,:), options{:}); 
    surface(elecwx2, elecwy2, elecwz2, img2(end:-1:1,:,:), options{:});
    surface(elecwx3, elecwy3, elecwz3, img3(end:-1:1,:,:), options{:}); 

    rotate3d on 
    
    planex = elecwx1(1,1);
    planey = elecwy2(1,1);
    planez = elecwz3(1,1);
end

function color = strcol2real(color)
    switch color
     case 'r', color = [1 0 0];
     case 'g', color = [0 1 0];
     case 'b', color = [0 0 1];
     case 'c', color = [0 1 1];
     case 'm', color = [1 0 1];
     case 'y', color = [1 1 0];
     case 'k', color = [0 0 0];
     case 'w', color = [1 1 1];
     otherwise, error('Unknown color'); 
    end
end

function x = gammacorrection(x, gammaval)
    x = 255 * (double(x)/255).^ gammaval;
end
    
