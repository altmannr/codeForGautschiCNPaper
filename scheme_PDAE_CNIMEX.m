function u = scheme_PDAE_CNIMEX(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,u_ref)
% Gautschi-type and implicit--explicit integrators for constrained wave-type systems
% R. Altmann, B. Dörich, C. Zimmer (2025)
% 
% Example: wave eqn with kinetic bc (with kappa=0, beta=1)
% 
% Numerical scheme in time:
% IMEX Crank--Nicolson for constrained PDAE formulation (explicit in nonlinearity) 
% no damping --> 2-step formulaiton
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
M_Om = M_Om; 
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

%% enlarge vectors to [u;p] to fit in structure of PDAE setting 
up_old = [u_old; p_old];
w_old = 0.*up_old;

%% first step - equation (4.3) in the paper
[f_Om_0,f_Ga_0] = func_rhs_waveKin(0,Nodes);    % time-dep part of rhs without nonlinearity

b_0 = [M_Om*f_Om_0; M_Ga*f_Ga_0(Boundary)];
nl_0 = [sparse(n,1); M_Ga*(p_old.^3-p_old)];
z = [M+0.25*tau^2*A, B'; B, sparse(m,m)] \ ...
    [M*up_old + tau*M*w_old - 0.25*tau^2*A*up_old + 0.5*tau^2*(b_0-nl_0); sparse(m,1)];
up_curr = z(1:n+m);

%% iteration loop - starting from step 2
for n_timestep = 2:N    
    t_curr = t(n_timestep); 
    
    %%%%%
    % rhs f_\Om and f_\Ga (at old and new time step)
    [f_Om_curr,f_Ga_curr] = func_rhs_waveKin(t_curr,Nodes);
    b_curr = [M_Om*f_Om_curr; M_Ga*f_Ga_curr(Boundary)];  
    p_curr = up_curr(n+1:end);
    nl_curr = [sparse(n,1); M_Ga*(p_curr.^3-p_curr)];
    
    %%%%%
    % CN (as presented in paper) - two-step formulation 
    z = [M+0.25*tau^2*A, B'; B, sparse(m,m)] \ ...
        [2*M*up_curr - M*up_old - 0.5*tau^2*A*up_curr - 0.25*tau^2*A*up_old + tau^2*(b_curr-nl_curr); sparse(m,1)];

    %%%%%
    % update past values of the solution
    up_old = up_curr;
    up_curr = z(1:n+m);
end

u = up_curr(1:n); 

% compare to exact/reference solution (if specified)
if nargin > 12
    err_u = u - u_ref;
    
    % L2-L2 error at t=T
    errors(1) = sqrt(err_u.'*(M_Om*err_u));  
    % energy error at t=T
    errors(2) = sqrt(err_u.'*(A_Om*err_u));
    
    % save error 
    save(strcat(['errors_PDAE_CNIMEX/error_PDAE_CNIMEX','_n',num2str(meshno),'_tau',num2str(tau),'.txt']), 'errors', '-ASCII');
end

end
