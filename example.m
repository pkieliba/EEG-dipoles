
inputFile = 'data/sampleEEG.mat';
load(inputFile); % load data, containing EEG structure and bootstrapped ICA variables

%%%%%%%%%% dipfit settings %%%%%%%%%%
f = which('eegplugin_dipfit');f2 = f(1:end-19);
chanfile = [f2 '/standard_BEM/elec/standard_1005.elc'];
mrifile = [f2 '/standard_BEM/standard_mri.mat'];
hdmfile = [f2 '/standard_BEM/standard_vol.mat'];
coord_transform = EEG.etc.coord_transform;
nbchan = EEG.nbchan;
EEG = pop_dipfit_settings(EEG, 'hdmfile', hdmfile, 'coordformat', 'MNI', 'mrifile', mrifile, 'chanfile', chanfile, 'coord_transform', coord_transform ,'chansel', 3:nbchan);

%%%%%%%%%% plot the ICA components %%%%%%%%%%
EEG.icawinv = A;
EEG.icaweights = W;
EEG.icasphere = eye(nbchan);
pop_topoplot(EEG, 0, [1:size(EEG.icawinv,2)], 'ICA'); 

%%%%%%%%%% calculate equivalent dipoles of the chosen ICs based on the bootrapped ICA %%%%%%%%%%
ics = [3, 15];
dipolesBoot = cell(1,length(A_boot_percomp));
dipolesBoot = calculateDipoles(dipolesBoot, EEG, A_boot_percomp, W_boot_percomp, ics);

%%%%%%%%%% plot dippole density %%%%%%%%%%
ic = 20; 
if isempty(dipolesBoot{ic})
    dipolesBoot = calculateDipoles(dipolesBoot, EEG, A_boot_percomp, W_boot_percomp, [ic]);
end
[dens3d, mri] = dipoledensity(dipolesBoot{ic}.model, 'coordformat', 'mni');
mri3dplot(dens3d,mri)

%%%%%%%%%% plot dipole as an ellipsoid and project onto 2D MNI space %%%%%%%%%%
ic = 26; 
if isempty(dipolesBoot{ic})
    dipolesBoot = calculateDipoles(dipolesBoot, EEG, A_boot_percomp, W_boot_percomp, ic);
end
centroid = computeCentroid(dipolesBoot{ic});
figure()
ellipsoidPlot(centroid, 'color', 'b');

%%%%%%%%%% project chosen dipoles to cortical surface in 3D %%%%%%%%%%
ics = [3,15,20,26]; 
for ic=ics
    if isempty(dipolesBoot{ic})
        dipolesBoot = calculateDipoles(dipolesBoot, EEG, A_boot_percomp, W_boot_percomp, ic);
    end
end
cortexProjection(dipfit_total,ics);
