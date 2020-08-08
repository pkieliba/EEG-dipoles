% Computes the equivalent dipoles from a bootstrapped set of ICs
%
% INPUTS
%
% EEG:  EEG structure from the eeglab
% A_boot_percomp:   a 1xn cell array, where each cell contains (nxm) matrix of inverse ica weight 
%                   of a single bootstrapped IC
% W_boot_percomp:   a 1xn cell array, where each cell contains (mxn) matrix
%                   of ica weight of a single bootstrapped IC
% ics:              an array with the indices of the ICs, whose dipoles
%                   should be calculated
%
% m: number of bootstrapps run
% n: number of ICs
%
% OUTPUT
%
% dipolesBoot - a cell array, where each cell contains a dipole structure
%               (calculated by dipfit) of all bootstraps of a given IC
%
function [dipolesBoot] = calculateDipoles(dipolesBoot, EEG, A_boot_percomp, W_boot_percomp, ics)
    
    f = which('eegplugin_dipfit');f2 = f(1:end-19);
    chanfile = [f2 '/standard_BEM/elec/standard_1005.elc'];
    mrifile = [f2 '/standard_BEM/standard_mri.mat'];
    hdmfile = [f2 '/standard_BEM/standard_vol.mat'];
    coord_transform = EEG.etc.coord_transform;
    nbchan = EEG.nbchan;
    EEG = pop_dipfit_settings(EEG, 'hdmfile', hdmfile, 'coordformat', 'MNI', 'mrifile', mrifile, 'chanfile', chanfile, 'coord_transform', coord_transform ,'chansel', 3:nbchan);
    EEG = eeg_checkset(EEG);
    EEG = pop_reref( EEG, []);
    EEG = eeg_checkset(EEG);
    
    for ic = ics
        dipolesBoot{ic} = [];
        EEG = pop_dipfit_settings(EEG, 'hdmfile', hdmfile, 'coordformat', 'MNI', 'mrifile', mrifile, 'chanfile', chanfile, 'coord_transform', coord_transform ,'chansel', 3:nbchan);
        EEG.icawinv = A_boot_percomp{ic}(:,:);
        EEG.icaweights = W_boot_percomp{ic}(:,:);
        EEG.icasphere = eye(length(A_boot_percomp));
        EEG = pop_multifit(EEG, 1:size(EEG.icawinv,2), 'threshold', 100);
        dipolesBoot{ic} = EEG.dipfit;
    end

end

