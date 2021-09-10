function [orthogonal_complements]               = Calculate_Orthogonal_Complements(spatial_coefficients,testmatrix,heading_space)
  number_selections = size(testmatrix,1);
  for iS = number_selections :-1:1
    selected_points = testmatrix(iS,:);
    
    tmp1 = (selected_points-1)*2+1;
    tmp2 = tmp1+1;
    indices = zeros(1,numel(selected_points)*2);
    
    indices(:,1:2:end) = tmp1;
    indices(:,2:2:end) = tmp2;
    
    t_coeff = spatial_coefficients.translational(indices,:);
    r_coeff = spatial_coefficients.rotational(indices,:);
    
    n = size(t_coeff,1)/2;
    
    nRows = size(t_coeff,1);
    nCols = n-3;
    
    ocs = zeros(nRows,nCols*heading_space.number_of_candidates);
    
    for iT = heading_space.number_of_candidates :-1:1
      T = heading_space.candidates_c(iT,:);
      
      A = Stretch_Matrix(t_coeff*T');
      CT = [A,r_coeff];
      
      [P, ~] = qr(CT);
      [m2, nCT] = size(CT);
      
      % last m-nCT columns provide base of orthogonal complement of CT
      OCCT = P(:, nCT + 1 : m2);
      
      ocs(:,(iT-1)*nCols+1:iT*nCols) = OCCT;
    end
    orthogonal_complements(iS).OCCT = ocs;
  end
end