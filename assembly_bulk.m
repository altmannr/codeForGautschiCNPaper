function [S,M]=assembly_bulk(Nodes,Elements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness and mass matrix assembly in the bulk 
% using first order finite element functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% inicialization
dof=length(Nodes);
M = spalloc(dof,dof,5*dof);
S = spalloc(dof,dof,5*dof);

%% matrices on the reference element
Mtx_m = [1/12    1/24    1/24
        1/24    1/12    1/24
        1/24    1/24    1/12];

% \pa_x v^2
Mtxx = [1/2   -1/2     0
       -1/2    1/2     0
         0      0      0];

% \pa_y v^2
Mtxy = [0     0       0
        0    1/2    -1/2
        0   -1/2     1/2 ];


% \pa_x v \pa_y v
Mtxxy = [0    -1/2    1/2 
         0    1/2     -1/2
         0     0       0  ];

no_Elements=size(Elements,1);

% M_Cell = cell(no_Elements,1);
% A_Cell = cell(no_Elements,1);

%% assembly using reference element
for i=1:no_Elements
    AA=[Nodes(Elements(i,1),1) Nodes(Elements(i,1),2)];
    BB=[Nodes(Elements(i,2),1) Nodes(Elements(i,2),2)];
    CC=[Nodes(Elements(i,3),1) Nodes(Elements(i,3),2)];

    M_A=[AA(1)-BB(1) CC(1)-BB(1) ; AA(2)-BB(2) CC(2)-BB(2)];
    
    C=inv(M_A);
    determ=det(M_A);
    
    S_loc=(Mtxx*(C(1,1)^2+C(1,2)^2) + ...
          (Mtxxy+Mtxxy')*(C(1,1)*C(2,1)+C(1,2)*C(2,2)) +...
           Mtxy*(C(2,1)^2+C(2,2)^2))*abs(determ);
       
    M_loc=Mtx_m*abs(determ);
    
    %% mass matrix
    M(Elements(i,:),Elements(i,:))=M(Elements(i,:),Elements(i,:))+M_loc;
    
    %% stiffness matrix
    S(Elements(i,:),Elements(i,:))=S(Elements(i,:),Elements(i,:))+S_loc;
end