function [residual_values] = Evaluate_Flowfields(motion,testmatrix,orthogonal_complements,heading_space)
  
  number_of_selections = size(testmatrix,1);
  number_of_selected_points= size(testmatrix,2);
  number_of_candidates = heading_space.number_of_candidates;
  
  
  
  U = motion.U;
  V = motion.V;
  
  
  if numel(size(U)) == 2
    number_of_motions = size(U,1);
    number_motion_variants = 1;
  elseif numel(size(U)) == 3
    number_of_motions = size(U,2);
    number_motion_variants = size(U,1);
  end
  
  for iMV = 1 : number_motion_variants
    if number_motion_variants > 1
      UMV = squeeze(U(iMV,:,:));
      VMV = squeeze(V(iMV,:,:));
    else
      UMV = U;
      VMV = V;
    end
    
    for iS =  number_of_selections : -1 : 1
      selection_index = testmatrix(iS,:);
      
      U_sel = UMV(:,selection_index);
      V_sel = VMV(:,selection_index);
      
      theta = zeros(number_of_motions,number_of_selected_points*2);
      
      theta(:,1:2:end) = U_sel;
      theta(:,2:2:end) = V_sel;
      
      r =  (theta*orthogonal_complements(iS).OCCT).^2;
      
      pre_residual_values(iS,iMV).R = reshape(sum(reshape(r',number_of_selected_points-3,[]))',number_of_candidates,[])';
      
    end
  end
  
  for iMV = 1 : number_motion_variants
    for iS = 1 : number_of_selections
      for iM = number_of_motions : -1 : 1
        residual_values(iM).R(iMV,iS,:) =  pre_residual_values(iS,iMV).R(iM,:);
      end
    end
  end
end