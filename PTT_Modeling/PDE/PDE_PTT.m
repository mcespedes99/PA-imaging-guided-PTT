% Ejemplo matlab esfera dentro de cubo: (0.5,0,0), (0,0,0), (0,0.7,0)
%% ---------------------------Montecarlo results for light distribution inside the tissue-------------------------
clear; clc;

%Loading of mat file containing the power density (W/cm^3):
load('MAT_files/Fluence_and_heat_rate.mat');
%Define laser power: (actual = 1 W)
global QL
QL(3,:)=0.1*QL(3,:);

%% -------------------------------PDE Solution------------------------------
% 1) Creation of the PDE's that we need to solve the problem:
%Main PDE
thermalmodel = createpde(1);

%% 2) Creation of Geometry:
load('MAT_files/Geo_PDE.mat')
geometryFromEdges(thermalmodel,dl);

%% 3) Plot of Cells y Faces of the Geometry
figure('Position',[10,10,800,400]);
subplot(1,2,1)
pdegplot(thermalmodel,'FaceLabel','on')
title('Geometry with Cell Labels')
subplot(1,2,2)
pdegplot(thermalmodel,'EdgeLabel','on')
title('Geometry with Edge Labels')

%% 4) Definition of the mesh:
mesh=generateMesh(thermalmodel,'Hmax',0.05,'Hmin', 0.04);
%figure
%pdemesh(mesh)

%% 5) Themal properties needed for the simulation to run:
%Laser on (L=1) or Laser off (L=0)
global L
L = 1;
% Skin:
%Constante "d": Cp*rho
d_skin = 0.0012*3800;
%Conductividad térmica: constante "c"
k_skin = 0.0053;
specifyCoefficients(thermalmodel,'Face',4,'m',0,'d',d_skin,"c",k_skin,"a",@a_skin,"f",@f_skin);
specifyCoefficients(thermalmodel,'Face',2,'m',0,'d',d_skin,"c",k_skin,"a",@a_skin,"f",@f_skin);

% Adipose:
%Constante "d": Cp*rho
d_sub = 0.0023*850;
%Conductividad térmica: constante "c"
k_sub = 0.0016;
specifyCoefficients(thermalmodel,'Face',5,'m',0,'d',d_sub,"c",k_sub,"a",@a_sub,"f",@f_sub);
specifyCoefficients(thermalmodel,'Face',6,'m',0,'d',d_sub,"c",k_sub,"a",@a_sub,"f",@f_sub);

% Muscle:
%Constante "d": Cp*rho
d_muscle =  0.00127*3500;
%Conductividad térmica: constante "c"
k_muscle = 0.0053;
specifyCoefficients(thermalmodel,'Face',1,'m',0,'d',d_muscle,"c",k_muscle,"a",@a_muscle,"f",@f_muscle);

% Damage tissue:
%Constante "d": Cp*rho
d_tumor_cte = 0.001*3500;
%Conductividad térmica: constante "c"
k_tumor_cte = 0.00642;
specifyCoefficients(thermalmodel,'Face',3,'m',0,'d',d_tumor_cte,"c",k_tumor_cte,"a",@a_tumor,"f",@f_tumor);
specifyCoefficients(thermalmodel,'Face',7,'m',0,'d',d_tumor_cte,"c",k_tumor_cte,"a",@a_tumor,"f",@f_tumor);
specifyCoefficients(thermalmodel,'Face',8,'m',0,'d',d_tumor_cte,"c",k_tumor_cte,"a",@a_tumor,"f",@f_tumor);

%% 6) Boundary conditions Laser:
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',8,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',9,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',10,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',1,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',5,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',6,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Edge',7,'u',37);

%Function to find the closest node to an specific coordinate:
getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);

% iii) Initial conditions:
%Initial temperature:
R = 37;
setInitialConditions(thermalmodel,R);

%% 7) PDE solution:
%thermalmodel.SolverOptions.MinStep = 0.5;
%model.SolverOptions.AbsoluteTolerance = 1.0e-4;
tlist = [0:0.25:300];
fprintf("solución");
R = solvepde(thermalmodel,tlist);
T = R.NodalSolution;

figure;
pdeplot(thermalmodel,'XYData',T(:,end),'Contour','on') %Índice 2 = 10 segundos
title(['Temperature at Time 300 s']);
colorbar
xlim([0 2])
ylim([-3 3])

%% ----------------------Tissue coefficients functions--------------------
%% Skin:
% a = -Cp*W(T)
function skin_a = a_skin(location,state)
skin_a = zeros(size(location.x));    
ids_piel_menor_44 = state.u<=44;
skin_a(ids_piel_menor_44) = (3500*0.45*(1+9.2*exp(-(state.u(ids_piel_menor_44)-44).^2 / 10)))/1000000;
ids_piel_mayor_44 = state.u>44;
skin_a(ids_piel_mayor_44) = (3500*0.45*10.2)/1000000;
end

% Coeficiente "f": QL+QM-Cp*Tb*W(T)
function f_s = f_skin(location,state)
global QL;
f_s = nan(size(location.x));
[~,c_s] = size(location.x);
x_array_s = [location.x];
y_array_s = [location.y];
t_array_s = [state.u];
%display(state.time);
%cond = isnan(state.time);
%if(cond == 1)
%    f_coef_s = nan(size(location.x));
%else
for i=1:c_s
    x = x_array_s(i);
    y = y_array_s(i);
    t = t_array_s(i);
    
    [~,nodo]= min((QL(1,:) - x).^2 + (QL(2,:) - y).^2); 
    P_laser = QL(3,nodo);
    %Primero, se coloca la convección con el ambiente en zona expuesta de la piel:
    if(x == 0)
        if(t <= 44)
           f_s(i) = P_laser + 0.001091 + (3500*(37)*0.45*(1+9.2*exp(-(t-44)^2)))/1000000 - 0.0005*(25-t);
        else
           f_s(i) = P_laser + 0.001091 + (3500*(37)*0.45*10.2)/1000000 - 0.0005*(25-t);
        end
        
    %Para la piel no expuesta:
    else
        if(t <= 44)
           f_s(i) = P_laser + 0.001091 + (3500*(37)*0.45*(1+9.2*exp(-(t-44)^2)))/1000000;
        else
           f_s(i) = P_laser + 0.001091 + (3500*(37)*0.45*10.2)/1000000;
        end
    end
end
end

%% Subcutaneous tissue:
% a = -Cp*W(T)
function sub_a = a_sub(location,state)
sub_a = zeros(size(location.x));    
ids_adipose_menor_45 = state.u<=45;
sub_a(ids_adipose_menor_45) = (3500*(0.36+0.36*exp(-(state.u(ids_adipose_menor_45)-45).^2 / 12)))/1000000;
ids_adipose_mayor_45 = state.u>45;
sub_a(ids_adipose_mayor_45) = (3500*0.72)/1000000;
end

% Coeficiente "f": QL+QM-Cp*Tb*W(T)
function f_sub = f_sub(location,state)
global QL;
f_sub = nan(size(location.x));
[~,c_s] = size(location.x);
x_array_s = [location.x];
y_array_s = [location.y];
t_array_s = [state.u];
%display(state.time);
%cond = isnan(state.time);
%if(cond == 1)
%    f_coef_s = nan(size(location.x));
%else
for i=1:c_s
    x = x_array_s(i);
    y = y_array_s(i);
    t = t_array_s(i);
    
    [~,nodo]= min((QL(1,:) - x).^2 + (QL(2,:) - y).^2); 
    P_laser = QL(3,nodo);
    if(t <= 45)
       f_sub(i) = P_laser + 0.001091 + (3500*(37)*(0.36+0.36*exp(-(t-45)^2/12)))/1000000;
    else
       f_sub(i) = P_laser + 0.001091 + (3500*(37)*0.72)/1000000;
    end
end
end

%% Muscle:
% a = -Cp*W(T)
function a_m = a_muscle(location,state)
a_m = zeros(size(location.x));    
ids_muscle_menor_45 = state.u<=45;
a_m(ids_muscle_menor_45) = (3500*(0.45+3.55*exp(-(state.u(ids_muscle_menor_45)-45).^2 / 12)))/1000000;
ids_muscle_mayor_45 = state.u>45;
a_m(ids_muscle_mayor_45) = (3500*4)/1000000;
end

% Coeficiente "f": QL+QM-Cp*Tb*W(T)
function f_m = f_muscle(location,state)
global QL;
f_m = nan(size(location.x));
[~,c_s] = size(location.x);
x_array_s = [location.x];
y_array_s = [location.y];
t_array_s = [state.u];
%display(state.time);
%cond = isnan(state.time);
%if(cond == 1)
%    f_coef_s = nan(size(location.x));
%else
for i=1:c_s
    x = x_array_s(i);
    y = y_array_s(i);
    t = t_array_s(i);
    
    [~,nodo]= min((QL(1,:) - x).^2 + (QL(2,:) - y).^2); 
    P_laser = QL(3,nodo);
    if(t <= 45)
       f_m(i) = P_laser + 0.001091 + (3500*(37)*(0.45+3.55*exp(-(t-45)^2 / 12)))/1000000;
    else
       f_m(i) = P_laser + 0.001091 + (3500*(37)*4)/1000000;
    end
end
end

%% Tumor:
% a = -Cp*W(T)
function a_t = a_tumor(location,state)
a_t = zeros(size(location.x));    
ids_t_menor_37 = state.u<37;
a_t(ids_t_menor_37) = (3500*(0.833))/1000000;
ids_t_entre_37_42 = state.u>=37 & state.u<=42;
a_t(ids_t_entre_37_42) = (3500*(0.833*(state.u(ids_t_entre_37_42)-37).^4.8 / 5438))/1000000;
ids_t_mayor_42 = state.u>42;
a_t(ids_t_mayor_42) = (3500*0.416)/1000000;
end

% Coeficiente "f": QL+QM-Cp*Tb*W(T)
function f_t = f_tumor(location,state)
global QL;
f_t = nan(size(location.x));
[~,c_s] = size(location.x);
x_array_s = [location.x];
y_array_s = [location.y];
t_array_s = [state.u];
%display(state.time);
%cond = isnan(state.time);
%if(cond == 1)
%    f_coef_s = nan(size(location.x));
%else
for i=1:c_s
    x = x_array_s(i);
    y = y_array_s(i);
    t = t_array_s(i);
    
    [~,nodo]= min((QL(1,:) - x).^2 + (QL(2,:) - y).^2); 
    P_laser = QL(3,nodo);
    %Primero, se coloca la convección con el ambiente en zona expuesta del tumor:
    if(x == 0)
        % Convección Q_c = -h*(Tamb-T) = -5 W/(K*m²)*(25-T) = 
        % modelo_ecuaciones.pdf= Thermal dosage investigation for optimal
        % temperature...-Yatao Ren
        if(t<37)
            f_t(i)= P_laser+0.001091+(3500*37*(0.833))/1000000 - 0.0005*(25-t);
        elseif(t>=37) && (t<=42)
            f_t(i)= P_laser+0.001091+(3500*37*(0.833*(t-37)^4.8 / 5438))/1000000 - 0.0005*(25-t);
        else
            f_t(i)= P_laser+0.001091+(3500*(37)*0.416)/1000000 - 0.0005*(25-t);
        end
        
    %Para el tumor no expuesto:
    else
        if(t<37)
            f_t(i)= P_laser+0.001091+(3500*37*(0.833))/1000000;
        elseif(t>=37) && (t<=42)
            f_t(i)= P_laser+0.001091+(3500*37*(0.833*(t-37)^4.8 / 5438))/1000000;
        else
            f_t(i)= P_laser+0.001091+(3500*(37)*0.416)/1000000;
        end
    end
end
end
