clearvars
close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Practice 2.5 Poisson Equation (2D Linear Triangles)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Poisson example
%
% We consider the Poisson equation on a triangle domain 
% with f(x,y) constant over the domain

% Coefficients:
% a11 = a22 = 1, a12 = a21 = a00 = 0, f = 1

% The domain is meshed by linear triangles. We recall that,
% when the coefficients of the model equafion are constant
% over the elment, there exist "closed" formulas that give
% the components of the local stiffness matrix. The same 
% works for the "source" term. 

% See the documents with the Theory's presentation at the
% Numerical Factory

% Geometry: nodes and elements
nodes=[
    0,0;
    0.5,0;
    0.5,0.5;
    1,0;
    1,0.5;
    1,1;
    ];

elem=[1,2,3;
      5,3,2;
      2,4,5;
      3,5,6;
      ];

numNod=size(nodes,1); %find out the number of nodes
numElem=size(elem,1); %find out the nimber of elements

% Plot the meshed domain
numbering=1; %if you want to see number of nodes and elements
plotElementsOld(nodes,elem,numbering);

% Coefficients of the modal equation 
a11=1;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=1;
coeff=[a11,a12,a21,a22,a00,f];

K = zeros(numNod);
Q = zeros(numNod,1);
F = zeros(numNod,1);

for e = 1:numElem
    %
    % Assemble the elements
    %
    [Ke,Fe] = linearTriangElement(coeff, nodes, elem, e);
    rows = [elem(e,1), elem(e,2), elem(e,3)];
    cols = rows;
    K(rows,cols) = K(rows,cols) + Ke;
    if (coeff(6) ~= 0)
        F(rows) = F(rows) + Fe;
    end
end
%
% Boundary Conditions
%
fixedNodes = [4,5,6]; %fixed nodes (global numbering)
freeNodes = setdiff(1:numNod,fixedNodes);

% Natural B.C.
Q(freeNodes) = 0; %In *this* case redundant, since Q
                  %Q was initialised to zero
% Essential B.C.
u = zeros(numNod,1);
u(fixedNodes) = 0;

% Reduced system 
Qm = F(freeNodes) + Q(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);

% Solution of the reduced system
um = Km\Qm;
u(freeNodes) = um;

% Post-process
Q = K*u - F;

% Table with solutions
sols = [(1:numNod)',nodes, u,Q];

%format short e
%format compact
fprintf("%5s%9s%14s%14s%14s\n","node","X", "Y", "U", "Q")
fprintf("%5d%14.6e%14.6e%14.6e%14.6e\n",sols');

% Color map for the solution 
titol='Poisson solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale)
