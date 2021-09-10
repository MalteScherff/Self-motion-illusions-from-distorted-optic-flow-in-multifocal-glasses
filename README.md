# Project PAL

This projects models a vision based self-motion estimation and explores the effects of how the distortion
of the visual input changes such estimations. The estimation process is based on the subspace algorithm by 
Heeger & Jepson ("Subspace Methods for Recovering Rigid Motion"), the distortion effects are coming from a
progressive addition lense (PAL) design (Carl Zeiss Vision GmbH, Aalen, Germany). This simulation is accompanied 
by a virtual reality based psychophysical experiment. Results are described in "Self-motion illusions from 
distorted optic flow in multifocal glasses" (not yet published).


#Simulation

The scripts and files here are the ones used to generate the data described in the aforementioned manuscript. 
Due to the deterministic nature of the model the exact values can be reconstructed, based on the predetermined 
spatial distribution of visible ground poins used. From that input (along with the simulation settings and 
heading space) the flowfields are calculated by first establishing the corresponding 3D environment that matches
the settings and the visible ground points and afterwards simulating the movement. The flowfields are then evaluated 
using an implementation of the subspace algorithm resulting in 50 heading estimates for each combination of
motion direction, distortion type used and effective field of view (FoV) size. Based on these estimation the 
PSE value for every distortion type and FoV size will be calculated, depicting the vertical heading component 
for which the the model most likely estimates an heading direction parallel to the ground. 

The simulation runs by defining the variable 'basefolder' to be the path to the folder 'Base' in this project, 
then executing the scripts 'Simulation_of_Movement.m' and 'Simulation_of_Heading_Estination.m'. 
'Simulation_Settings.m' can be used to alter the settings beforehand, although only in a limited way (Heading 
directions could be chosen freely but the distortion fits available only corresponds to the directions used). 
Results can be plotted using 'Plot_Barplot.m'.



# Sripts

1) Unpack_Project_PAL.m

Easy method to open the relevant scripts by categories.

2) Simulation_Settings.m

Input the settings for the simulation. Most relevant are  the relation between the simulated eye and the 
3D environment (angle between line of sight and ground plane, height of the eye, distance to the 
image plane), the different stimulus level (added vertical heading component) and the list of radi that 
includes the effective FoV sizes used later.

3) Simulation_of_Movement.m

A) Load the simulation settings, the predetermined spatial distribution of visible ground points in retinal 
coordinates and the fitted lense data. 

B) Calculate the 3D environment that fits the settings and visible points. 

C) Simulate the movement through that environment to get flowfields.

D) Include lense distortions for undistorted flowfields.

E) Save the results as 'Simulated_Movement.mat'.


4) Simulation_of_Heading_Estimation.m

A) Load the simulation settings, the calculated flowfields and the heading space. 

B) Calculate orthogonal complements for combinations of points in the flowfields and candidate directions
stored in the heading space.

C) Use orthogonal complements to evaluate flowfields, results in residual surfaces.

D) Combine residual surfaces with resect to the effective FoV sizes to determine candidate directions most 
fitting for the flowfields.

E) Calculate PSEs for the estimated heading directions.

F) Save the results as 'Heading Estimates.mat' and 'PSEs.mat'.


5) Calculate_Spatial_Coefficients.m

The matrices A(T) and B are calculated in reference to the aforementioned work "Subspace Methods for Recovering
Rigid Motion" by Heeger and Jepson. They are used for flowfield calculations as well as for the calculation of 
orthogonal complements.

6) Flowfield_Calculation.m

7) Calculate_Orthogonal_Complements.m

8) Evaluate_Flowfields.m

9) Process_Estimates.m

10) Stretch_Matrix.m

Functions to reshape the matrices containing the spatial coefficients. Needed at some point during the calculation
of orthogonal complements.

11) Plot_Barplot.m

# .mat files

1) Settings.mat

Contains the settings for the simulation.

2) Heading_Space.mat

Contains all the heading directions for which the eligibility for the flowfields is tested.

3) Pre_Grids.mat

Contains the spatial distribution of visible ground points in retinal coordinates.

4) Prepared_Lense_Data.mat

Contains the fitted functions describing the distortions when looking through certain points on a certai PAL.

5) Figure_Settings.mat

Settings for plotting.

6) Simulated_Movement.mat

Contains the flowfields generated after the first simulation step.

7) Heading_Estimates.mat

Contains the heading estimates.

8) PSEs.mat

Contains the points of subjective equalities.



# Sources

The estimation process is based on the work of Heeger and Jepson:

"Subspace Methods for Recovering Rigid Motion" (1992)

DOI: 10.1007/BF00128130


To calculate the PSE values Psignifit was used, developed by Schütt et al.

"Painfree and accurate bayesian estimation of psychometric functions for (potentially) overdispersed data" (2016)

DOI: 10.1016/j.visres.2016.02.002


The function to calculate theradius of circular receptive fields is drived from Albright and Desimone:

"Local precision of visuotopic organization in the middle temporal area (MT) of the macaque"

DOI: 10.1007/BF00235981

# Contact

malte.scherff@uni-muenster.de

