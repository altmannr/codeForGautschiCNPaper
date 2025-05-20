%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  [cosVal, H, V] = cosv_SaddleKrylov( tau, A, M, B, v, tol, m ) 
%
%  code for paper: 
%  R. Altmann, B. Dörich, C. Zimmer
%  Gautschi-type and implicit--explicit integrators for constrained wave-type systems
%
%  computes approximation of cos(tau*Om_ker)*v without computing Om_ker=sqrt(A_ker)
%  Krylov space is spanned by v, A_ker v, A_ker^2v,...
%
%  %% inputs:
%  tau = step size
%  A = matrix which defines A_ker
%  M = mass matrix  
%  B = constraint matrix
%  v = vector which is multiplied 
%  tol = NOT YET USED!
%  m = number of Arnoldi steps (dimension of Krylov subspace)
%
%  %% outputs:
%  cosVal = cos(tau*Om_ker)*v
%  H = Hessenberg matrix s.t. A*V ~ V*H
%  V = orthogonalized span of Krylov subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cosVal, H, V] = cosv_SaddleKrylov( tau, A, M, B, v, tol, max_iter ) 

[n,n] = size(A);
if nargin == 5
  tol = 1.0e-7;
  max_iter = min(n,15);
end
if nargin == 6
  max_iter = min(n,15);
end

% normalize first vector
vNorm = norm(v);
v = v/vNorm;

% build up V (for m=1)
V = v;              % matrix of orthonormal vectors 
H = zeros(max_iter,max_iter);     % Hessenberg matrix

R = [M, B'; B, sparse(size(B,1),size(B,1))];
% loop
for k = 2 : max_iter    
    % add new vector and normalize
    v_new = R\[A*V(:,end);zeros(size(B,1),1)];
    v_new = v_new(1:size(A,1));
    for j = 1 : k-1
        H(j,k-1) = V(:,j)'*v_new;
        v_new = v_new - H(j,k-1)*V(:,j);
    end
    H(k,k-1) = norm(v_new);
    V(:,end+1) = v_new/H(k,k-1);  
end

% up to now, V contains m orthogonal vectors; H is of size m*(m-1)
% we now add the last column
v_new = R\[A*V(:,end);zeros(size(B,1),1)];
v_new = v_new(1:size(A,1));

for j = 1 : max_iter
    H(j,max_iter) = V(:,j)'*v_new;
    v_new = v_new - H(j,max_iter)*V(:,j);
end

% use matrix-cos for Hessenbergmatrix
e1 = eye(max_iter,1); 
sqrtH = sqrtm(H);
cosVal = vNorm*V*cosm(tau*sqrtH)*e1;
