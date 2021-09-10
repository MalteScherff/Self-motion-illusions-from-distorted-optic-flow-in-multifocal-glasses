function [spatial_coefficients]         = Calculate_Spatial_Coefficients(x,y,f)
  % The translational and rotational spatial coefficients are used for the
  % flowfield calculation and, more importantely, the subspace algorithm
  % later. This is done to increase performance.
  k = numel(x);
  A = zeros(2*k, 3);
  B = zeros(2*k, 3);
  
  A(1:2:end,1:3) = [repmat(-f,k,1), zeros(k,1), x'];
  A(2:2:end,1:3) = [zeros(k,1), repmat(-f,k,1),  y'];
  
  B(1:2:end,1:3) = [x'.*y'/f, repmat(-f,k,1)-x'.^2/f, y'];
  B(2:2:end,1:3) = [repmat(f,k,1)+y'.^2/f, -x'.*y'/f,  -x'];
  
  spatial_coefficients.translational = A;
  spatial_coefficients.rotational = B;
end

