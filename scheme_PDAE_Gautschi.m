function u = scheme_PDAE_Gautschi(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,krylovDim,u_ref)
% Gautschi-type and implicit--explicit integrators for constrained wave-type systems
% R. Altmann, B. Dörich, C. Zimmer (2025)
% 
% Example: wave eqn with kinetic bc 
% 
% Numerical scheme in time:
% Gautschi for constrained PDAE formulation 
% 
% Output:
% u and u' at final time
%
% System equations:
% u'' - Laplace(u) + kappa*u = sin(t)  in Omega
% u'' - beta*LaplaceBeltrami(u) + partial_n*u = - u^3 + u  on Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TIME INTEGRATION
N = ceil(T/tau);

% time vector
t = (0:tau:T);
% errors at final time (u in L2- and energy norm)
errors = zeros(2,0);

%% matrices
Boundary = unique(Boundary_Edges);
%M_Om = M_Om; 
M_Ga = M_Ga(Boundary,Boundary);
A_Om = A_Om + kappa*M_Om;
A_Ga = beta*A_Ga(Boundary,Boundary);
% 
n = length(M_Om);
m = length(M_Ga);
M = [M_Om, sparse(n,m); sparse(m,n), M_Ga];
A = [A_Om, sparse(n,m); sparse(m,n), A_Ga];
%
B = sparse(m,n+m);
B(:,Boundary) = -M_Ga; 
B(:,n+1:end) = M_Ga;

%% starting values (as in HipK20, Sect 8.1)
x = Nodes(:,1);
y = Nodes(:,2);
u_old = exp(-20*((x-1).^2+y.^2));
p_old = u_old(Boundary);
w_old = 0.*u_old;                   % = dt_u_old
r_old = w_old(Boundary);

%% enlarge vectors to [u;p] to fit in structure of PDAE setting 
up_old = [u_old; p_old];
wr_old = [w_old; r_old];

%% compute u_1
up_curr = up_old + tau*wr_old;
[f_Om,f_Ga] = func_rhs_waveKin(0,Nodes);    
b_0 = [M_Om*f_Om; M_Ga*f_Ga(Boundary)];
nl_0 = [sparse(n,1); -M_Ga*(p_old.^3-p_old)];
% compute A_ker*(u_0-Bg_0)
%tmp = [M, B'; B, sparse(m,m)] \ [A*up_old;sparse(m,1)];
%ddotup = M \ (b_0 + nl_0 - tmp(1:n+m));
%up_curr = up_old + tau*wr_old + 0.5*tau^2*ddotup;  % !!! noch einschalten
tmp = [M, B'; B, sparse(m,m)] \ [b_0 + nl_0 - A*up_old; sparse(m,1)];
ddotup = tmp(1:n+m);
up_curr = up_old + tau*wr_old + 0.5*tau^2*ddotup;

%% iteration loop 
for n_timestep = 2:N    
    t_curr = t(n_timestep);
    p_curr = up_curr(n+1:end);

    %%%%%
    % rhs f_\Om and f_\Ga (at curr time step)
    [f_Om_curr,f_Ga_curr] = func_rhs_waveKin(t_curr,Nodes);    
    b_curr = [M_Om*f_Om_curr; M_Ga*f_Ga_curr(Boundary)];  
   
    %%%%%
    % Gautschi step 
    nl_curr = [sparse(n,1); -M_Ga*(p_curr.^3-p_curr)];
    tmp = [A, B'; B, sparse(m,m)] \ [b_curr + nl_curr; sparse(m,1)];
    bn = tmp(1:n+m);

    %%%%% 
    % Krylov  
    cosProduct = cosv_SaddleKrylov(tau, A, M, B, up_curr-bn, 0, krylovDim);
    up_new = - up_old + 2*cosProduct + 2*bn;  

    %%%%%
    % update past values of the solution
    up_old = up_curr;
    up_curr = up_new;
end

u = up_new(1:n); 

% compare to exact/reference solution (if specified)
if nargin > 13
    err_u = u - u_ref;
    
    % L2-L2 error at t=T
    errors(1) = sqrt(err_u.'*(M_Om*err_u));  
    % energy error at t=T
    errors(2) = sqrt(err_u.'*(A_Om*err_u));
    
    % save error 
    save(strcat(['errors_PDAE_Gautschi/error_PDAE_Gautschi','_n',num2str(meshno),'_tau',num2str(tau),'.txt']), 'errors', '-ASCII');
end

end
