%% Grab data
cd(basefolder)
addpath(genpath(basefolder))

cd('Data')
load('Settings.mat')
load('Pre_Grids.mat')
cd('Lense_Data')
load('Prepared_Lense_Data.mat')

number_of_repetitions                                       = numel(iteration);
simulation_settings.number_of_repetitions                   = number_of_repetitions;


%% Repeat the simulation for every instance of the structure iteration
for iR = 1  : number_of_repetitions
  %% Step 1: Construct the environment
  % Define the environment based on 2D input in retinal coordinates   and
  % the simulation settings
  [ground_plane,visual_grid,testmatrix,info_struct]         = Recover_Spatial_Data(iteration(iR).pre_grid,...
                                                                simulation_settings,fit_cell);
  
  %% Step 2: Spatial coefficients
  x                                                         = ground_plane.xy(1,:);
  y                                                         = ground_plane.xy(2,:);
  sp_coeff                                                  = Calculate_Spatial_Coefficients(x,y,...
                                                                simulation_settings.distance_image_plane);
  
  %% Step 3: Simulate Movement
  motion                                                    = Calculate_Flowfields(ground_plane,...
                                                                simulation_settings,sp_coeff,fit_cell);
  
  simulated_movement(iR).ground_plane                       = ground_plane;
  simulated_movement(iR).visual_grid                        = visual_grid;
  simulated_movement(iR).testmatrix                         = testmatrix;
  simulated_movement(iR).info_struct                        = info_struct;
  
  simulated_movement(iR).sp_coeff                           = sp_coeff;
  simulated_movement(iR).motion                             = motion;
end

%% Save results
cd(basefolder)
cd('Data/Results')

save('Simulated_Movement','simulated_movement')




%% Functions
function [ground_plane,visual_grid,testmatrix,info_struct]  = Recover_Spatial_Data(pre_grid,simulation_settings,fit_cell)
  
  number_of_points                                      = size(pre_grid(1).xy_rc,2);
  number_of_fields                                      = numel(pre_grid);
  
  testmatrix                                            = zeros(number_of_fields, number_of_points);
  
  xy_rc                                                 = zeros(2,number_of_fields*number_of_points);
  max_point_distance_per_field                          = zeros(1,number_of_fields);
  min_point_distance_per_field                          = zeros(1,number_of_fields);
  RF_sizes                                              = zeros(1,number_of_fields);
  max_eccentricity_of_field                             = zeros(1,number_of_fields);
  center_point_distance                                 = zeros(1,number_of_fields);
  
  
  for iF = number_of_fields : -1 : 1
    indices                                             = (iF-1)*number_of_points+1:iF*number_of_points;
    visual_grid(iF).xy_rc                               = pre_grid(iF).xy_rc;
    visual_grid(iF).center_rc                           = pre_grid(iF).center_rc;
    visual_grid(iF).indices                             = indices;
    visual_grid(iF).RF_size                             = simulation_settings.rfs(norm(visual_grid(iF).center_rc));
    visual_grid(iF).max_eccentricity_of_field(1,iF)     = norm(visual_grid(iF).center_rc) + visual_grid(iF).RF_size;
    
    distances                                           = vecnorm(visual_grid(iF).xy_rc);
    
    xy_rc(:,indices)                                    = visual_grid(iF).xy_rc;
    testmatrix(iF,:)                                    = indices;
    
    
    %% Info Structure
    max_point_distance_per_field(1,iF)                  = max(distances);
    min_point_distance_per_field(1,iF)                  = min(distances);
    RF_sizes(1,iF)                                      = visual_grid(iF).RF_size;
    max_eccentricity_of_field(1,iF)                     = visual_grid(iF).max_eccentricity_of_field(1,iF);
    center_point_distance(1,iF)                         = norm(visual_grid(iF).center_rc);
  end
  
  ground_plane_tmp                                      = Recover_Ground_Plane(xy_rc,simulation_settings);
  
  ground_plane.xy                                       = ground_plane_tmp.xy;
  ground_plane.xy_rc                                    = xy_rc;
  ground_plane.XYZ                                      = ground_plane_tmp.XYZ;
  ground_plane.distance_fixation_point                  = ground_plane_tmp.distance_fp;
  
  info_struct.maximum_distance_of_point_per_field       = max_point_distance_per_field;
  info_struct.minimum_distance_of_point_per_field       = min_point_distance_per_field;
  info_struct.size_of_receptive_fields                  = RF_sizes;
  info_struct.maximal_eccentricity_covered_per_field    = max_eccentricity_of_field;
  info_struct.distances_of_receptive_field_centers      = center_point_distance;
  
  
  %%% Displacement
  for iL = 1 : size(fit_cell,2)
    p                                                   = ground_plane.xy;
    
    Fx_inv                                              = fit_cell{4,iL};
    Fy_inv                                              = fit_cell{5,iL};
    
    x                                                   = Fx_inv(p(1,:),p(2,:));
    y                                                   = Fy_inv(p(1,:),p(2,:));
    
    pd                                                  = [x;y];
    pd_rc(1:2,:)                                        = [simulation_settings.cart2rc(x,simulation_settings.distance_image_plane);...
                                                            simulation_settings.cart2rc(y,simulation_settings.distance_image_plane)]; 
        
    ground_plane.lense(iL).xy                           = pd;
    ground_plane.lense(iL).xy_rc                        = pd_rc;
    ground_plane_tmp                                    = Recover_Ground_Plane(ground_plane.lense(iL).xy_rc,...
                                                            simulation_settings);
    
    ground_plane.lense(iL).XYZ                          = ground_plane_tmp.XYZ;
    ground_plane.lense(iL).distance_fixation_point      = ground_plane_tmp.distance_fp;
  end
  
  
end

function [ground_plane]                                     = Recover_Ground_Plane(xy_rc,simulation_settings)
  %%% fov = [fov_horizontal_left,fov_horizontal_right; ...
  %%% ... fov_vertical_up, fov_vertical_down]
  %%% alpha is the angle between the line of sight and the plane
  downwards_angle                                       = simulation_settings.downwards_angle;
  eye_height                                            = simulation_settings.eye_height;
  distance_image_plane                                  = simulation_settings.distance_image_plane;
  
  % 1) Length: View line to FP
  a                                                     = eye_height/sind(downwards_angle);
  
  % 2) Z-Coordinate of Foot
  c                                                     = eye_height*sind(downwards_angle);
  
  % 3) y-coordinate of Foot
  d                                                     = cosd(downwards_angle)*eye_height*-1;
  
  % 4) Fixpoint
  FP_3D                                                 = [0; 0; a];
  
  % 5) Foot
  LP_3D                                                 = [0; d; c];
  
  % 6) Normal of the plane
  n_vec_3D                                              = -LP_3D/norm(LP_3D);
  
  % 7) Value of plane belonging
  w_3D                                                  = dot(n_vec_3D,LP_3D);
  
  
  for i = size(xy_rc,2): -1: 1
    x(i)                                                = simulation_settings.rc2cart(xy_rc(1,i),distance_image_plane);
    y(i)                                                = simulation_settings.rc2cart(xy_rc(2,i),distance_image_plane);
  end
  
  for i = numel(x) : -1: 1
    Z(i)                                                = w_3D/(n_vec_3D(1)*x(i)+n_vec_3D(2)*y(i)+n_vec_3D(3));
    X(i)                                                = x(i)*Z(i);
    Y(i)                                                = y(i)*Z(i);
  end
  
  ground_plane.XYZ                                      = [X;Y;Z];
  ground_plane.xy                                       = distance_image_plane*[ground_plane.XYZ(1,:)./ground_plane.XYZ(3,:);...
                                                            ground_plane.XYZ(2,:)./ground_plane.XYZ(3,:)];
  ground_plane.xy_rc                                    = xy_rc;
  ground_plane.distance_fp                              = FP_3D(3);
end

function [motion]                                           = Calculate_Flowfields(ground_plane,simulation_settings,sp_coeff,fit_cell)
  XYZ                                                   = ground_plane.XYZ;
  xy                                                    = ground_plane.xy;
  
  speed                                                 = simulation_settings.observer_speed;
  angle                                                 = simulation_settings.downwards_angle;
  f_r                                                   = simulation_settings.distance_image_plane;
  delta                                                 = simulation_settings.delta;
  
  heading_directions                                    = simulation_settings.heading_directions;
  
  number_of_heading_directions                          = size(heading_directions,1);
  number_of_points                                      = size(xy,2);
  
  motion.U                                              = zeros(3,number_of_heading_directions, number_of_points);
  motion.V                                              = zeros(3,number_of_heading_directions, number_of_points);
  
  for iDirection = 1 : size(heading_directions(:,1))
    direction                                           = heading_directions(iDirection,:);
    heading_v                                           = direction + [0, angle];
    
    T                                                   = [simulation_settings.rc2cart(heading_v(1),f_r),...
                                                            simulation_settings.rc2cart(heading_v(2),f_r), f_r];
    T                                                   = speed*ground_plane.distance_fixation_point*T./norm(T);
    Omega                                               = [T(2),-T(1),0]/ground_plane.distance_fixation_point;
    
    undist_motion                                       = Flowfield_Calculation(sp_coeff,XYZ(3,:),T,Omega);
    
    motion.U(1,iDirection,:)                            = undist_motion.U;
    motion.V(1,iDirection,:)                            = undist_motion.V;
    
    if direction(1) == -20
      key                                               = 'Right on lense';
    elseif direction(1) == 20
      key                                               = 'Left on lense';
    end
    
    
    for j = 1 : 2
      switch j
        case 1
          fit_cell_index                                = find(cellfun(@(x) isequal(x,key), fit_cell(1,:), 'UniformOutput', 1));
          row_index                                     = 2;
        case 2
          fit_cell_index                                = 3;
          row_index                                     = 3;
      end
      
      Fx                                                = fit_cell{2,fit_cell_index};
      Fy                                                = fit_cell{3,fit_cell_index};
      
      xy                                                = ground_plane.lense(fit_cell_index).xy;
      XYZ                                               = ground_plane.lense(fit_cell_index).XYZ;
      
      x                                                 = xy(1,:);
      y                                                 = xy(2,:);
      
      sp_coeff_redist                                   = Calculate_Spatial_Coefficients(x,y,...
                                                            simulation_settings.distance_image_plane);
      
      motion_dist                                       = Flowfield_Calculation(sp_coeff_redist,XYZ(3,:),...
                                                            T,Omega);
      
      U                                                 = motion_dist.U;
      V                                                 = motion_dist.V;
      xdpdelta                                          = Fx(x + delta,y);
      xdmdelta                                          = Fx(x - delta,y);
      distortion.fxdx                                   = 0.5*(xdpdelta - xdmdelta)./delta;
      
      xdpdelta                                          = Fx(x,y + delta);
      xdmdelta                                          = Fx(x,y - delta);
      distortion.fxdy                                   = 0.5*(xdpdelta - xdmdelta)./delta;
      
      ydpdelta                                          = Fy(x + delta,y);
      ydmdelta                                          = Fy(x - delta,y);
      distortion.fydx                                   = 0.5*(ydpdelta - ydmdelta)./delta;
      
      ydpdelta                                          = Fy(x,y + delta);
      ydmdelta                                          = Fy(x,y - delta);
      distortion.fydy                                   = 0.5*(ydpdelta - ydmdelta)./delta;
      
      Ud                                                = distortion.fxdx.*U + distortion.fxdy.*V;
      Vd                                                = distortion.fydx.*U + distortion.fydy.*V;
      
      motion.U(row_index,iDirection,:)                  = Ud;
      motion.V(row_index,iDirection,:)                  = Vd;
      
    end
  end
end









