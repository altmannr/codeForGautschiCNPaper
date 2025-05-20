function run_kineticWave
% Run-script serving as a unified interface for various convergence tests
% for parabolic problems with dynamic boundary conditions.
% 
% Uses a formulation based on partial differential algebraic equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;

%% time integrator, meshes, time steps, and coefficients, etc.
% meshes
meshno_vect=[4];      % range: 0,...,9  
% time steps
tau_vect=2.^(-1:-1:-9); % zB (-4:-1:-10)

% time interval
T = 1; 

% coefficients appearing in the problem 
kappa = 1; % weight of u in Omega
beta = 1;  % weight of Laplace-Beltrami

%% a loop over the meshes
for ih=1:length(meshno_vect)
    meshno = meshno_vect(ih)
    
    % loading current mesh from the directory 'mesh'
    Nodes=load(['meshes/disk_nodes',num2str(meshno),'.txt']);
    Elements=load(['meshes/disk_elements',num2str(meshno),'.txt']);
    
    %% loading boundary edges or preprocessing and saving them
    if exist(['meshes/disk_boundary',num2str(meshno),'.txt'], 'file')==2
        Boundary_Edges=load(['meshes/disk_boundary',num2str(meshno),'.txt']);
    else
        [Boundary_Edges]=preprocess_mesh(Nodes,Elements);
        save(['meshes/disk_boundary',num2str(meshno),'.txt'], 'Boundary_Edges', '-ASCII');
    end
    
    %% loading existing matrices or assembleing them
    if exist(['meshes/matrices_disk',num2str(meshno),'.mat'], 'file')==2
        matrices=load(['meshes/matrices_disk',num2str(meshno),'.mat']);
        A_Om=matrices.A_Om;
        M_Om=matrices.M_Om;
        A_Ga=matrices.A_Ga;
        M_Ga=matrices.M_Ga;
    else
        % BULK mass and stiffness matrix assembly
        [A_Om,M_Om]=assembly_bulk(Nodes,Elements);

        % BOUNDARY mass and stiffness matrix assembly 
        [A_Ga,M_Ga]=assembly_surface(Nodes,Boundary_Edges);
        
        s1.A_Om=A_Om;
        s1.M_Om=M_Om;
        s1.A_Ga=A_Ga;
        s1.M_Ga=M_Ga;
        save(['meshes/matrices_disk',num2str(meshno),'.mat'],'-struct','s1');
    end
        
    %% use finest time mesh to produce reference solution 
    tau_ref = 2^(-11);  
    u_ref = scheme_PDAE_Gautschi(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau_ref,beta,kappa,Elements,10);  

    %% a loop over step sizes
    for jtau = 1:length(tau_vect)
        tau = tau_vect(jtau)        
        % 
        scheme_PDAE_Euler(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,u_ref); 
        scheme_PDAE_CNIMEX(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,u_ref); 
        scheme_PDAE_Gautschi(Nodes,Boundary_Edges,A_Om,M_Om,A_Ga,M_Ga,T,meshno,tau,beta,kappa,Elements,5,u_ref); 
    end
end

end
