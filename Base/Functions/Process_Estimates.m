function [PSEs] = Process_Estimates(estimates, simulation_settings)
  stimulus_level                    = unique(simulation_settings.heading_directions(:,2)');
  downwards_angle                   = simulation_settings.downwards_angle;
  distance_list                     = simulation_settings.distance_list;
  number_of_distances               = numel(distance_list);
  number_of_stimulus_level          = numel(stimulus_level);
  number_of_repetitions             = simulation_settings.number_of_repetitions;
  
  
  estimates_peripheral(1:number_of_distances,1:number_of_stimulus_level, 1:number_of_repetitions)                           = ...
    squeeze(estimates(2,:,1:number_of_stimulus_level,:,2));
  estimates_peripheral(1:number_of_distances,1:number_of_stimulus_level, number_of_repetitions+1:2*number_of_repetitions)   = ...
    squeeze(estimates(2,:,number_of_stimulus_level+1:2*number_of_stimulus_level,:,2));
  
  estimates_central(1:number_of_distances,1:number_of_stimulus_level, 1:number_of_repetitions)                              = ...
    squeeze(estimates(3,:,1:number_of_stimulus_level,:,2));
  estimates_central(1:number_of_distances,1:number_of_stimulus_level, number_of_repetitions+1:2*number_of_repetitions)      = ...
    squeeze(estimates(3,:,number_of_stimulus_level+1:2*number_of_stimulus_level,:,2));
  
  
  percentage_upwards_answers_peripheral = sum(1*((estimates_peripheral-downwards_angle)>0) + .5*(estimates_peripheral-downwards_angle==0),3)/size(estimates_peripheral,3);
  percentage_upwards_answers_central = sum(1*((estimates_central-downwards_angle)>0) + .5*(estimates_central-downwards_angle==0),3)/size(estimates_peripheral,3);
  
  PSEs = zeros(numel(distance_list),2);
  for iD = 1 : numel(distance_list)
    
    output_evaluated_per.n_answered         = ones(number_of_stimulus_level,1)*number_of_repetitions*2;
    output_evaluated_per.angles             = stimulus_level';
    
    output_evaluated_cen                    = output_evaluated_per;
    
    output_evaluated_per.percent_up         = percentage_upwards_answers_peripheral(iD,:)';
    output_evaluated_per.n_up               = output_evaluated_per.percent_up.*output_evaluated_per.n_answered;
    
    output_evaluated_cen.percent_up         = percentage_upwards_answers_central(iD,:)';
    output_evaluated_cen.n_up               = output_evaluated_cen.percent_up.*output_evaluated_cen.n_answered;
    
    
    [~,PSEs_b]                              = Custom_fit_session(output_evaluated_per,output_evaluated_cen);
    
    PSEs(iD,:) = PSEs_b;
    
  end

  
end


function [output_fitted,PSEs] = Custom_fit_session(output_evaluated1,output_evaluated2)
  %% psignifit options
  options = struct;
  options.sigmoidName = 'norm';   % choose a cumulative Gauss as the sigmoid
  options.expType = 'equalAsymptote'; % free but euqal lower and upper asymptotes
  
  
  
  %% dataset 1
  angles = output_evaluated1.angles;
  n_up = output_evaluated1.n_up;
  n_answered = output_evaluated1.n_answered;
  percent_up = output_evaluated1.percent_up;
  
  index = isfinite(percent_up);
  x = angles(index);
  y1 = n_up(index);
  y2 = n_answered(index);
  res0 = psignifit([x,y1,y2],options);
  
  %% dataset 2
  angles = output_evaluated2.angles;
  n_up = output_evaluated2.n_up;
  n_answered = output_evaluated2.n_answered;
  percent_up = output_evaluated2.percent_up;
  
  index = isfinite(percent_up);
  x = angles(index);
  y1 = n_up(index);
  y2 = n_answered(index);
  res1 = psignifit([x,y1,y2],options);
  
  
  output_fitted = [res0 res1];
  PSEs = [res0.Fit(1) res1.Fit(1)];
  
end
