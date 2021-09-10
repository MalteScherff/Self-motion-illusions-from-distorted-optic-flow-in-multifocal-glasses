function [str_M]                                = Stretch_Matrix(M)
  % Turns a collection of (2xn) matrices A,B,C,D saved as one matrix
  % (A                  (A 0 0 0
  %  B          into     0 B 0 0
  %  C                   0 0 C 0
  %  D)                  0 0 0 D)
  
  
  num_dim = 2;
  
  [l,n] = size(M);
  k = l/num_dim;
  
  inds = [ones(1,n*k)+((1:n*k)-1)*l + reshape(repmat(((1:k)-1),n,1),1,[])*2;
    2*ones(1,n*k)+((1:n*k)-1)*l + reshape(repmat(((1:k)-1),n,1),1,[])*2];
  
  num_chunks = numel(inds)/prod([num_dim,n]);
  M2 = zeros(num_chunks*num_dim,n);
  
  for iChunk = 1 : num_chunks
    M2((iChunk-1)*num_dim+1:iChunk*num_dim,:) = inds(:,1:n);
    inds(:,1:n) = [];
  end
  inds = M2;
  str_M = zeros(l,n*k);
  str_M(inds) = M;
end