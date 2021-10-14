%% ---------------------------PA signal-------------------------
%% Montecarlo results:
clear; clc;

%Loading of mat file containing the power density (W/cm^3):
load('MAT_files/Fluence_and_heat_rate.mat');

%Loading of mat file with temperature data:
%load("PA+PTT/temperature2D.mat");

%Definir potencia del laser: (actual = 2 W)
global QL
QL(3,:)=0.05*QL(3,:);
%Change to make it look the temperature matrix
pixel_size = 0.05; %0.5 mm
x = 0:pixel_size:2;
y = -3:pixel_size:3;
[~,a] = size(x);
[~,b] = size(y);
Q_L = zeros(3,a*b);
counter = 1;
for i=1:a
    actual_x = x(1,i);
    for j=1:b
        actual_y = y(1,j);
         [~,nodo]= min((QL(1,:) - actual_x).^2 + (QL(2,:) - actual_y).^2);
         Q_L(:,counter)= [actual_x,actual_y,QL(3,nodo)]; 
         counter = counter+1;
    end
end

%Matrix for PA 0:
PA_0 = zeros(3,a*b);
PA_0(1:2,:)=Q_L(1:2,:);
PA_0(3,:)=1*Q_L(3,:);

%Matrix for PA:
[~,t_size] = size(tlist);
PA_t = zeros(3,a*b,t_size);
PA_t(:,:,1) = PA_0;
for j=1:t_size
    PA_t(1:2,:,j) = PA_0(1:2,:);
end

%Function to find the closest node to an specific coordinate:
getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);

%Matrix to find nodes from T matrix with flujo coordinates:
col = a*b;
nodes = zeros(1,col);
for i=1:col
    x = PA_0(1,i);
    y = PA_0(2,i);
    [~,node]= getClosestNode(mesh.Nodes,x,y);
    nodes(1,i)=node;
end

% For to complete PA matrix in time: PA = PA_0*T/T_0 con T_0=37C
for k=1:col
    PA_t(3,k,2:t_size)=PA_0(3,k)*T(nodes(k),2:t_size)/37;
end

%Adding of noise to PA_t signal:
PA_tn = PA_t;
for k=1:col
    PA_tn(3,k,1:t_size)=awgn(PA_t(3,k,1:t_size),30,'measured');
end

%% Adding of movement to PA_t signal:
% First, define PA window: This is equivalent to PA_0 in the window of
% interest:
x = 0:pixel_size:1;
y = -1:pixel_size:1;
[~,a] = size(x);
[~,b] = size(y);
PA_window = zeros(3,a*b);
counter = 1;
for i=1:a
    actual_x = x(1,i);
    for j=1:b
        actual_y = y(1,j);
         [~,nodo]= min((PA_0(1,:) - actual_x).^2 + (PA_0(2,:) - actual_y).^2);
         PA_window(:,counter)= [actual_x,actual_y,PA_0(3,nodo)]; 
         counter = counter+1;
    end
end

%Define PA signal with movement:
%First create a bigger movement window. Window will be moved a maximum of 2
%mm, therefore, I set it to be a square with x=[0,1.6]
PA_tm = zeros(3,a*b,t_size);
PA_tm(:,:,1) = PA_window;
