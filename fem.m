%% 2-D GRID STRUCTURE - FEM ANALYSIS
clear;clc;close all

% Material, area and load
E=210e9; A=5e-3; P=500;
scale=500;                         % DOUBLE deformation visualization

% Nodes [x y]
N=[0 0;2 0;4 0;6 0;8 0;10 0;12 0;
   0 1;2 1;4 1;6 1;8 1;10 1;12 1;
   0 2;2 2;0 3;2 3;0 4;2 4;0 5;2 5;
   0 6;2 6;0 7;2 7;0 8;2 8;0 9;2 9];

% Elements
e=[1 2;2 3;3 4;4 5;5 6;6 7;
   8 9;9 10;10 11;11 12;12 13;13 14;
   1 8;2 9;3 10;4 11;5 12;6 13;7 14;
   1 9;8 2;2 10;9 3;3 11;10 4;4 12;11 5;
   5 13;12 6;6 14;13 7;
   8 15;15 17;17 19;19 21;21 23;23 25;25 27;27 29;
   9 16;16 18;18 20;20 22;22 24;24 26;26 28;28 30;
   8 9;15 16;17 18;19 20;21 22;23 24;25 26;27 28;29 30;
   8 16;15 9;15 18;17 16;17 20;19 18;19 22;21 20;
   21 24;23 22;23 26;25 24;25 28;27 26;27 30;29 28];

% Global stiffness matrix
nd=2*size(N,1);
K=zeros(nd);

for i=1:size(e,1)
    n1=e(i,1); n2=e(i,2);
    dx=N(n2,1)-N(n1,1);
    dy=N(n2,2)-N(n1,2);
    L=hypot(dx,dy);
    c=dx/L;
    s=dy/L;

    ke=E*A/L*[c^2 c*s -c^2 -c*s;
               c*s s^2 -c*s -s^2;
              -c^2 -c*s c^2 c*s;
              -c*s -s^2 c*s s^2];

    d=[2*n1-1 2*n1 2*n2-1 2*n2];
    K(d,d)=K(d,d)+ke;
end

% Load
F=zeros(nd,1);
F(14)=-P;                         % 500 N downward at Node 7

% Supports
% Node 30 fixed in X and Y
% Node 29 restrained in Y
fix=[58 59 60];
free=setdiff(1:nd,fix);

% FEM solution
U=zeros(nd,1);
U(free)=K(free,free)\F(free);

% Displacements
Ux=U(1:2:end);
Uy=U(2:2:end);
D=hypot(Ux,Uy);

[dm,nm]=max(D);

% Deformed coordinates
Nd=N+scale*[Ux Uy];

%% Plot
figure('Color','w','Position',[100 50 1100 650])
hold on
grid on
box on

% Original and deformed members
for i=1:size(e,1)

    n1=e(i,1);
    n2=e(i,2);

    % Original - Black
    plot(N([n1 n2],1),N([n1 n2],2), ...
        'k-','LineWidth',1.1)

    % Deformed - Red
    plot(Nd([n1 n2],1),Nd([n1 n2],2), ...
        'r-','LineWidth',1.8)
end

% Original nodes
plot(N(:,1),N(:,2),'ko', ...
    'MarkerFaceColor','k','MarkerSize',5)

% Deformed nodes
plot(Nd(:,1),Nd(:,2),'ro', ...
    'MarkerFaceColor','r','MarkerSize',5)

% Node numbers
for i=1:size(N,1)

    text(N(i,1)+0.08,N(i,2)+0.07, ...
        num2str(i),'Color','k','FontSize',8)

    text(Nd(i,1)+0.08,Nd(i,2)-0.10, ...
        num2str(i),'Color','r','FontSize',8)
end

% Fixed support
plot(N(30,1),N(30,2),'ks', ...
    'MarkerFaceColor','k','MarkerSize',10)

text(2.25,9.15,'FIXED', ...
    'FontWeight','bold')

% Load arrow
quiver(12,0.75,0,-0.5,0, ...
    'b','LineWidth',2,'MaxHeadSize',0.7)

text(11.65,0.9,'500 N', ...
    'FontWeight','bold')

% Axis labels
xlabel('X (m)','FontSize',12)
ylabel('Y (m)','FontSize',12)

% Title
title(sprintf( ...
    'Original (black) vs Deformed (red) - Scale = %g',scale), ...
    'FontSize',14,'FontWeight','bold')

% Legend
legend('Original','Deformed (scaled)', ...
    'Location','northeast')

axis equal
xlim([-1 13.5])
ylim([-1 10])

% Maximum displacement
text(0.3,9.55, ...
    sprintf('Max displacement = %.3e m at node %d',dm,nm), ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'FontSize',10)

% Command window result
fprintf('\n============================================\n')
fprintf('       2-D GRID FEM ANALYSIS\n')
fprintf('============================================\n')
fprintf('Maximum displacement = %.6e m\n',dm)
fprintf('Maximum displacement at Node = %d\n',nm)
fprintf('Visualization scale = %g\n',scale)
fprintf('============================================\n')