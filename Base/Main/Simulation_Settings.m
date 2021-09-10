%% Spatial relations between observer ('s eye) and the environment.
% Flowfields are calculated on an image plane perpendicular to the line of
% sight (Z-axis).
simulation_settings.distance_image_plane    = 1;

% Angle between line of sight and ground plane
simulation_settings.downwards_angle         = 20;

% Distance between "eye" and ground plane
simulation_settings.eye_height              = 1800;



%% Movement settings
% Movement is parallel to the ground plane + the values stored in
% stimulus_level
stimulus_level                              = -3.5:.5:3.5;

% Actual heading directions in retinal coordinates
simulation_settings.heading_directions      = [...
                                                repmat(-20,numel(stimulus_level),1),stimulus_level';...
                                                repmat(20,numel(stimulus_level),1),stimulus_level'];

% Flowfields are calculated for a movement in heading direction T with 
% |T| = observer_speed * distance_to_fixation_point                                               
simulation_settings.observer_speed          = 1;

% Used in approximating the derivatives of the functions describing the
% distortions.
simulation_settings.delta                   = 1e-04;


%% Data Analysis
% Heading estimates are calculated for flowfields of varying sizes, meaning
% the circular field of view gets restricted to a radius with values
% defined here and all flow vectors neglected that start outside that
% radius.
simulation_settings.distance_list           = 40:1:55;


%% Functions
% Some functions used throughout the simulation, including the conversion
% from retinal coordinates to cartesian coordinates and vice versa, as well
% as the function to calculate the radius of circular receptive fields depending on
% their eccentricity. 
simulation_settings.rc2cart                 = @(alpha,z) z*tand(alpha);
simulation_settings.cart2rc                 = @(l,z) atand(l/z);
simulation_settings.rfs                     = @(ecc) (1.04+.61.*ecc)/sqrt(pi);     




%% Save
cd(basefolder)
cd('Data')

save('Settings','simulation_settings')



