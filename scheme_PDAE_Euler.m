function u = scheme_PDAE_Euler(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,u_ref)
% Gautschi-type and implicit--explicit integrators for constrained wave-type systems
% R. Altmann, B. Dörich, C. Zimmer (2025)
% 
% Example: wave eqn with kinetic bc 
% 
% Numerical scheme in time:
% IMEX Euler for constrained PDAE formulation (explicit in nonlinearity) 
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

% matrix for Euler iteration
iterMatrix = [M, -tau*M, sparse(n+m,m); ...
              tau*A, M, tau*B'; ...
              B, sparse(m,n+m), sparse(m,m)];

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

%% iteration loop 
for n_timestep = 1:N    
    t_new = t(n_timestep+1);
    
    %%%%%
    % rhs at new time step
    [f_Om_new,f_Ga_new] = func_rhs_waveKin(t_new,Nodes);    
    b_new = [M_Om*f_Om_new; M_Ga*f_Ga_new(Boundary)];
    p_old = up_old(n+1:end);
    nl_old = [sparse(n,1); -M_Ga*(p_old.^3-p_old)];

    %%%%%
    % Euler (nonlinearity explicit)
    z = iterMatrix \ [M*up_old; M*wr_old + tau*b_new + tau*nl_old; sparse(m,1)];

    %%%%%
    % update past values of the solution
    up_old = z(1:(n+m));
    wr_old = z(n+m+1:2*(n+m));
end

u = z(1:n); 

% compare to exact/reference solution (if specified)
if nargin > 12
    err_u = u - u_ref;
    
    % L2-L2 error at t=T
    errors(1) = sqrt(err_u.'*(M_Om*err_u));  
    % energy error at t=T
    errors(2) = sqrt(err_u.'*(A_Om*err_u));
    
    % save error 
    save(strcat(['errors_PDAE_Euler/error_PDAE_Euler','_n',num2str(meshno),'_tau',num2str(tau),'.txt']), 'errors', '-ASCII');
end

end
