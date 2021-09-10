%% Flowfield
function [motion] = Flowfield_Calculation(spatial_coefficients,Z,T,Omega)
  k = numel(Z);
  p = 1./reshape([Z;Z],1,2*k)';
  
  A = p.*spatial_coefficients.translational*T';
  B = spatial_coefficients.rotational*Omega';
  
  UV = A+B;
  UV = reshape(UV,2,numel(UV)/2);
 
  motion.U = UV(1,:);
  motion.V = UV(2,:);
end