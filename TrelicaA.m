
%---------------------------
%  control input data
%---------------------------

clear all; clf
nel=9;           % number of elements
nnel=2;          % number of nodes per element
ndof=2;          % number of dofs per node
nnode=6;         % total number of nodes in system
sdof=nnode*ndof; % total system dofs

%---------------------------
%  nodal coordinates
%---------------------------

gcoord(1,1)=0.0;  gcoord(1,2)=0.0;   % x, y-coordinate of node 1
gcoord(2,1)=20.0; gcoord(2,2)=0.0;   % x, y-coordinate of node 2
gcoord(3,1)=40.0;  gcoord(3,2)=0.0;  % x, y-coordinate of node 3
gcoord(4,1)=0.0;  gcoord(4,2)=15.0;   % x, y-coordinate of node 4
gcoord(5,1)=20.0; gcoord(5,2)=15.0;   % x, y-coordinate of node 5
gcoord(6,1)=40.0;  gcoord(6,2)=15.0;  % x, y-coordinate of node 6

%------------------------------------------
%  material and geometric properties
%------------------------------------------

elprop(1,1)=200000000000;  % elastic modulus of 1st element
elprop(1,2)=2500*10^-6;       % cross-section of 1st element
elprop(2,1)=200000000000;  % elastic modulus of 2nd element
elprop(2,2)=2500*10^-6;       % cross-section of 2nd element
elprop(3,1)=200000000000;  % elastic modulus of 3rd element
elprop(3,2)=2500*10^-6;       % cross-section of 3rd element
elprop(4,1)=200000000000;  % elastic modulus of 4th element
elprop(4,2)=2500*10^-6;       % cross-section of 4th element
elprop(5,1)=200000000000;  % elastic modulus of 5th element
elprop(5,2)=2500*10^-6;       % cross-section of 5th element
elprop(6,1)=200000000000;  % elastic modulus of 6th element
elprop(6,2)=2500*10^-6;       % cross-section of 6th element
elprop(7,1)=200000000000;  % elastic modulus of 7th element
elprop(7,2)=2500*10^-6;       % cross-section of 7th element
elprop(8,1)=200000000000;  % elastic modulus of 8th element
elprop(8,2)=2500*10^-6;       % cross-section of 8th element
elprop(9,1)=200000000000;  % elastic modulus of 9th element
elprop(9,2)=2500*10^-6;       % cross-section of 9th element

%-----------------------------
%  nodal connectivity
%-----------------------------

nodes(1,1)=1;  nodes(1,2)=2;   % nodes associated with element 1
nodes(2,1)=2;  nodes(2,2)=3;   % nodes associated with element 2
nodes(3,1)=1;  nodes(3,2)=4;   % nodes associated with element 3
nodes(4,1)=1;  nodes(4,2)=5;   % nodes associated with element 4
nodes(5,1)=2;  nodes(5,2)=5;   % nodes associated with element 5
nodes(6,1)=2;  nodes(6,2)=6;   % nodes associated with element 6
nodes(7,1)=3;  nodes(7,2)=6;   % nodes associated with element 7
nodes(8,1)=4;  nodes(8,2)=5;   % nodes associated with element 8
nodes(9,1)=5;  nodes(9,2)=6;   % nodes associated with element 9
%-----------------------------
%  applied constraints
%-----------------------------

bcdof = [2; 6];
bcval = [0; 0];


%----------------------------
%  initialization to zero
%----------------------------

ff=zeros(sdof,1);              % system force vector
kk=zeros(sdof,sdof);           % system stiffness matrix
index=zeros(nnel*ndof,1);      % index vector
elforce=zeros(nnel*ndof,1);    % element force vector
eldisp=zeros(nnel*ndof,1);     % element nodal displacement vector
k=zeros(nnel*ndof,nnel*ndof);  % element stiffness matrix
stress=zeros(nel,1);           % stress vector for every element
%-----------------------------
%  applied nodal force
%-----------------------------

ff(8)=-10000;     % 2nd node has 1000 lb in downward direction
ff(10)=-10000;     % 2nd node has 1000 lb in downward direction
ff(12)=-10000;     % 2nd node has 1000 lb in downward direction
%--------------------------
%  loop for elements
%--------------------------

for iel=1:nel    % loop for the total number of elements

nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);  % coordinate of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);  % coordinate of 2nd node

leng=sqrt((x2-x1)^2+(y2-y1)^2);  % element length

if (x2-x1)==0;
if y2>y1;
   beta=2*atan(1);       % angle between local and global axes
else
   beta=-2*atan(1);
end
else
beta=atan((y2-y1)/(x2-x1));
end

el=elprop(iel,1);               % extract elastic modulus
area=elprop(iel,2);             % extract cross-sectional area

index=feeldof(nd,nnel,ndof);  % extract system dofs for the element

k=fetruss2(el,leng,area,0,beta,1); % compute element matrix

kk=feasmbl1(kk,k,index);           % assemble into system matrix

end

%---------------------------------------------------
%  apply constraints and solve the matrix
%---------------------------------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);  % apply the boundary conditions

disp=kk\ff;   % solve the matrix equation to find nodal displacements

%--------------------------------------------------
%  post computation for stress calculation
%--------------------------------------------------

for iel=1:nel         % loop for the total number of elements

nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);  % coordinate of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);  % coordinate of 2nd node

leng=sqrt((x2-x1)^2+(y2-y1)^2);  % element length

if (x2-x1)==0;
if y2>y1;
   beta=2*atan(1);       % angle between local and global axes
else
   beta=-2*atan(1);
end
else
beta=atan((y2-y1)/(x2-x1));
end

el=elprop(iel,1);               % extract elastic modulus
area=elprop(iel,2);             % extract cross-sectional area

index=feeldof(nd,nnel,ndof);  % extract system dofs for the element

k=fetruss2(el,leng,area,0,beta,1); % compute element matrix

for i=1:(nnel*ndof)           % extract displacements associated with
eldisp(i)=disp(index(i));     % (iel)-th element
end

elforce=k*eldisp;             % element force vector
stress(iel)=sqrt(elforce(1)^2+elforce(2)^2)/area; % stress calculation

if ((x2-x1)*elforce(3)) < 0;
stress(iel)=-stress(iel);
end

if ((x2 - x1) == 0)
  if ((y2-y1)*elforce(3)) < 0;
    stress(iel)=-stress(iel);
  end
end


end

%----------------------------
% print fem solutions
%----------------------------

num=1:1:sdof;
displ=[num' disp]          % print displacements

numm=1:1:nel;
stresses=[numm' stress]    % print stresses

%--------------------------------------------------------------------

kk_A = kk / (200000000000*2500*10^-6)

disp_B = disp * (200000000000*2500*10^-6) / 2000

forcas_nas_barras = stress * 2500*10^-6

figure;
hold on;
title('Forças Axiais nas Barras (compressão x tração)');
xlabel('x [m]');
ylabel('y [m]');
axis equal;

% Define o mapa de cores de compressão (vermelho) a tração (azul)
cmap = jet(256);   % azul (valores altos) → vermelho (valores baixos)

% Centraliza a barra de cores no zero (compressão negativa, tração positiva)
max_abs_force = max(abs(forcas_nas_barras));
caxis([-max_abs_force, max_abs_force]);  % centraliza a colorbar no zero
colormap(cmap);
cb = colorbar;
ylabel(cb, 'Força axial [N]');

% Loop para plotar cada barra
for iel = 1:nel
    n1 = nodes(iel,1);
    n2 = nodes(iel,2);

    x = [gcoord(n1,1), gcoord(n2,1)];
    y = [gcoord(n1,2), gcoord(n2,2)];

    % Cor baseada no valor real da força (negativa = compressão)
    cor_idx = round( (forcas_nas_barras(iel) + max_abs_force) / (2*max_abs_force) * 255 ) + 1;
    cor_idx = max(min(cor_idx, 256), 1);  % garante que está dentro dos limites
    cor = cmap(cor_idx, :);

    % Largura da linha proporcional à força (valor absoluto)
    linewidth = 1 + 4 * abs(forcas_nas_barras(iel)) / max_abs_force;

    plot(x, y, 'Color', cor, 'LineWidth', linewidth);

    % Escrever o valor da força no meio da barra
    xm = mean(x);
    ym = mean(y);
    text(xm, ym, sprintf('%.1f N', forcas_nas_barras(iel)), ...
        'FontSize', 8, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'w', 'Margin', 1);
end

hold off;

