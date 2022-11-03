clearvars
close all

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
numNod=size(nodes,1);
numElem=size(elem,1);

numbering=1; %if you want to see number values
plotElements(nodes,elem,numbering);

%Coefficients 
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
    [Ke,Fe] = linearTriangElement(coeff, nodes, elem, e);
    rows = [elem(e,1), elem(e,2), elem(e,3)];
    cols = rows;
    K(rows,cols) = K(rows,cols) + Ke;
    if (coeff(6) ~= 0)
        F(rows) = F(rows) + Fe;
    end
end

% Boundary Conditions
fixedNodes = [4,5,6];
freeNodes = setdiff(1:numNod,fixedNodes);

%Essential B.C.
Q(freeNodes) = 0; %Redudant
%Natural B.C
u = zeros(numNod,1);
u(fixedNodes) = 0;

%Reduced system 
Qm = F(freeNodes) + Q(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);

um = Km\Qm;
u(freeNodes) = um;

Q = K*u - F;

format short e
[u, Q]

titol='Poisson solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale)
