%% Grab data
cd(basefolder)
cd('Data')
load('Heading_Space.mat')
load('Settings.mat')
cd('Results')
load('Simulated_Movement.mat')

number_of_repetitions                           = numel(simulated_movement);
simulation_settings.number_of_repetitions       = number_of_repetitions;
estimates                                       = zeros(3,numel(simulation_settings.distance_list),...
  size(simulation_settings.heading_directions,1),numel(iteration),2);



%% Initialise a waitbar
tic
h = waitbar(0,'Please wait...');


%% Repeat the simulation for every instance of the structure simulated movement
for iR = 1 : number_of_repetitions
  
  testmatrix                                    = simulated_movement(iR).testmatrix;
  info_struct                                   = simulated_movement(iR).info_struct;
  sp_coeff                                      = simulated_movement(iR).sp_coeff;
  motion                                        = simulated_movement(iR).motion;
  
  %% Step 1: Calculate orthogonal complements
  orthos                                        = Calculate_Orthogonal_Complements(sp_coeff, testmatrix,...
    heading_space);
  
  
  %% Step 2: Calculate residual maps
  residual_surfaces                             = Evaluate_Flowfields(motion,testmatrix,orthos,heading_space);
  
  %   fn = ['Worker_',num2str(iR)];
  %     save(fn,'residual_surfaces')
  
  %% Step 3: Heading estimation
  max_point_ecc                                 = info_struct.maximum_distance_of_point_per_field;
  
  for iM = 1 : numel(residual_surfaces)
    R_undistorted                               = squeeze(residual_surfaces(iM).R(1,:,:));
    R_peripheral                                = squeeze(residual_surfaces(iM).R(2,:,:));
    R_central                                   = squeeze(residual_surfaces(iM).R(3,:,:));
    
    for iD = 1 : numel(simulation_settings.distance_list)
      act_distance                              = simulation_settings.distance_list(iD);
      relevance_index                           = max_point_ecc <= act_distance;
      
      r_und                                     = sum(R_undistorted(relevance_index,:),1);
      r_per                                     = sum(R_peripheral(relevance_index,:),1);
      r_cen                                     = sum(R_central(relevance_index,:),1);
      
      [~,c_und]                                 = min(r_und,[],2);
      [~,c_per]                                 = min(r_per,[],2);
      [~,c_cen]                                 = min(r_cen,[],2);
      
      estimates(1,iD,iM,iR,1:2)                 = heading_space.candidates(c_und,:);
      estimates(2,iD,iM,iR,1:2)                 = heading_space.candidates(c_per,:);
      estimates(3,iD,iM,iR,1:2)                 = heading_space.candidates(c_cen,:);
    end
  end
  
  
  %% Update the waitbar
  
  waitbar(iR/number_of_repetitions)
  
  process_time                                  = toc;
  estimated_time                                = number_of_repetitions*process_time/iR-process_time;
  h.Children.Title.String                       = ['Please wait... Estimated Time ', num2str(round(estimated_time*10)/10),...
                                                    's, Run ', num2str(iR),' of ',num2str(number_of_repetitions),' completed'];
  
end
cd(basefolder)
cd('Data/Results')
save('Heading Estimates','estimates')

PSEs                                            = Process_Estimates(estimates,simulation_settings);
save('PSEs','PSEs')





