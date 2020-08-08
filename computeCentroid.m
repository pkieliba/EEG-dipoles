% Computes the mean and standard deviation of the dipoles' coordinates
% Takes as an argument a single dipfit dipole structure, computed from the 
% boottrapped distribution of the same IC.
%
function [centroid] = computeCentroid(dipoleBootstrap)

    centroid.posxyz = [0 0 0];
    centroid.momxyz = [0 0 0];
    
    [dipobj, ~] = convert_dipole_structure_to_array(dipoleBootstrap.model);
    
    outX = isoutlier(dipobj.location(:,1));
    outY = isoutlier(dipobj.location(:,2));
    outZ = isoutlier(dipobj.location(:,3));
    
    centroid.posxyz = [mean(dipobj.location(~outX,1)), mean(dipobj.location(~outY,2)), mean(dipobj.location(~outZ,3))];
    centroid.loc_std = [std(dipobj.location(~outX,1)), std(dipobj.location(~outY,2)), std(dipobj.location(~outZ,3))];
    
end