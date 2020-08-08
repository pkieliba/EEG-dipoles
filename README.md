# Visualisation of Dipolar Brain Sources

Many independent EEG components have scalp maps that nearly perfectly match the projection of a single equivalent brain dipole. The position of those dipoles can be easily computed using EEGLAB plugin [DIPFIT](https://github.com/sccn/dipfit). However a major obstacle to using this approach for the precise visualisation of the macroscopic brain dynamics lies in the fact that the dipole calculation depends heavily on the computed IC scalp maps whereas for several ICA algorithms, including Infomax, the IC scalp maps for even the same data set may differ slightly across runs. Moreover, as the calculated dipole is just a single point in a 3D space, it can be considered insufficient for visualisation of a brain source - as it is not likely that the observed brain activity originates from this single point only.

### This repository contains several Matlab scripts providing alternative ways of visualising dipolar brain sources in 3D.
To account for the differences between ICA runs and provide more accurate brain source position, these visualisation are based on the **bootstrap data decomposition procedure, called [RELICA](https://github.com/sccn/relica).** RELICA decomposes the data using Infomax with the ErpICASSO bootstrap approach and clusters the ICs from every bootstrap decomposition according to their similarity.

**ellipsoidPlot** calculates the position of the centroid of the dipoles’ cluster and its standard deviation along the x-, y- and z-direction and computes an equivalent ellipsoid. This ellipsoid is then plotted in the 3D space along with its 2D projections.

![ellipsoid plot of an equivalent dipole](https://github.com/pkieliba/EEG-dipoles/blob/master/images/ellipsoidDipole.png?raw=true)

**cortexProjection** uses [Measure Projection Analysis toolbox](https://bitbucket.org/bigdelys/measure-projection/src/default/) to model each of the dipoles by a Gaussian density function. It allows for simultaneous projection of several ICs onto the same 3D brain model. Each of the components can be plotted in different colour with the size of the patch reflecting on the spread of the given brain source

![cortex projection 4 equivalent dipoles](https://github.com/pkieliba/EEG-dipoles/blob/master/images/cortexProjection.png?raw=true)

### Required toolboxes/plugins

* [EEGlab](https://sccn.ucsd.edu/eeglab/index.php)
* [Dipfit](https://github.com/sccn/dipfit)
* [Measure Projection Analysis](https://bitbucket.org/bigdelys/measure-projection/src/default/)
