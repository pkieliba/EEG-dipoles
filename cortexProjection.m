% Plots the cortex plot with chosen ICs (dipoles) depicted in different colors.
%
% INPUTS
%
% dipoles:  a cell array containing dipole structres, caluclating separately
%           for each IC, using EEGlab Dipfit plugin
% ics:      array containing the numbers of ICs to be plotted
% colors:   an array of colors (RGB values) used for projecting the ICs (optional)
%
% EXAMPLE
%
% cortexProjection(dipfit_all, [3 47], [0 0 1; 0 1 0]) 
%
function cortexProjection(dipoles, ics, colors)

    if nargin < 1
        help cortexProjection;
        return
    end

    if nargin<3
        % deafault color values
        colors = [0 0 1; 0 1 0; 1 0.843137 0; 1 0 0];
    end

    for n = 1:length(ics)

        i = ics(n);
        dipole = dipoles{i}.model;
        [dipObj, ~] = convert_dipole_structure_to_array(dipole);
        
        outX = isoutlier(dipObj.location(:,1));
        outY = isoutlier(dipObj.location(:,2));
        outZ = isoutlier(dipObj.location(:,3));
        outAll = sum([outX, outY, outZ],2);
        
        
        dipObj.coordinateFormat = 'mni';
        dipObj.headGrid = pr.headGrid();

        minValueNormalizedTo = 0.06;
        projectionParameter = pr.projectionParameter;
        projectionParameter.standardDeviationOfEstimatedDipoleLocation = 5;
        projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 3;  
        projectionParameter.normalizeInBrainDipoleDenisty = [];
        
        [~, dipoleDensity]= pr.meanProjection.getProjectionMatrix(dipObj, dipObj.headGrid, projectionParameter, dipObj.headGrid.insideBrainCube);

        fsf = load('MNImesh_dipfit.mat');
        cortexVertices = fsf.vertices;

        cortexPointDomainDenisty(:,n) = pr.project_domain_to_cortex(dipObj.headGrid.getPosition, cortexVertices, dipObj.headGrid.spacing, dipoleDensity);
        
    end

    pr.plot_cortex(cortexPointDomainDenisty, colors, 'minValueNormalizedTo', minValueNormalizedTo);

end