%% ---------------------------Acople con Montecarlo-------------------------
%Limpieza de consola y del workspace:
clear; clc;

%Carga de matriz de interés:
load('MAT_files/Montecarlo_fluence.mat');
%En cada fila de H están las coordenadas de cada tetrahedro. ValoMC usa vectores para definir el tetraedro entonces cada fila de 
%H contiene las filas del array "r" en las que están contenidas los vectores que forman el tetrahedro. En r, lo que están son las
%coordenadas de estos vectores que forman los tetraedros. https://inverselight.github.io/ValoMC/structures3d.html
H = vmcmesh.H;
r = vmcmesh.r;

%Definición del tamaño de matriz de interés: El número de filas de H es
%equivalente al número de tetrahedros o subdivisiones de la geometria.
[m,n] = size(H);

%Variable fluence para guardar el fluence y sus coordenadas. Tiene 4 filas
%porque, para cada elemento de la geometría se van a guardar sus
%coordenadas (x,y) y el fluence asociado a ese elemento.
global fluence;
fluence = zeros(3,m);

%Bucle para recorrer todas las filas y columnas de H:
%El for está pensado para que vaya fila por fila de H (tetraedro por tetraedro),
%sacando las coordenadas de cada uno de los vectores que forman el
%tetraedro. Luego se hace un promedio en cada uno de los ejes x,y,z para
%lograr calcular una coordenada promedio del centro del tetraedro. Luego
%estas coordenadas las pase a esfericas para calcular la tasa de calor
%generada por el laser, dada por la irradiancia (fluence) por mu_a.
%Finalmente, las guarde en el array flujo que creé antes.
for i=1:m
    coordenadas = zeros(3,2);
    contador = 1;
    for j=1:n
       fila_r = H(i,j);
       coordenadas(contador,:)= r(fila_r,:);
       contador = contador + 1;
    end
    xc = mean(coordenadas(:,1))/10;
    yc = mean(coordenadas(:,2))/10;
    fluence(:,i)=[xc,yc,solution.element_fluence(i)*100];
end
fluence(1,:) =  round(fluence(1,:),4);
fluence(2,:) =  round(fluence(2,:),4);
fluence(3,:) =  round(fluence(3,:),4);

%Esto es para plotearlo bonito. Si quiere entenderlo, me escribe y le
%ayudo, si no logra sacar.
[x,y]=meshgrid(0:0.025:2, -3:0.025:3);
[a,b]=size(x);
fluenceXY = zeros(size(x));
for i=1:b
    actual_x = x(1,i);
    for j=1:a
        actual_y = y(j,1);
        [~,nodo]= min((fluence(1,:) - actual_x).^2 + (fluence(2,:) - actual_y).^2);
        fluenceXY(j,i)= fluence(3,nodo);
    end
end

%%Calculation of QL:
global QL;
QL = zeros(3,m);
for i=1:m
    coordenadas = zeros(3,2);
    contador = 1;
    for j=1:n
       fila_r = H(i,j);
       coordenadas(contador,:)= r(fila_r,:);
       contador = contador + 1;
    end
    xc = mean(coordenadas(:,1))/10;
    yc = mean(coordenadas(:,2))/10;
    [~,radio] = cart2pol(xc,yc);
    if(radio <= 0.505)
        if(ismember(i,chosen_nodes))
           node = chosen_nodes(1,:)==i;
           fluence_Ne = solution.element_fluence(i)*(chosen_nodes(2,node))*1000;
        else
           fluence_Ne = solution.element_fluence(i)*0.021433576*1000; 
        end
    elseif(xc<=0.13)
        fluence_Ne = solution.element_fluence(i)*0.0207638*1000;
    elseif(xc <= 0.305) && (xc > 0.13)
         fluence_Ne = solution.element_fluence(i)*0.008352283*1000;
    elseif(xc > 0.305)
         fluence_Ne = solution.element_fluence(i)*0.021433576*1000;  
    end
    QL(:,i)=[xc,yc,fluence_Ne];
end
QL(1,:) =  round(QL(1,:),4);
QL(2,:) =  round(QL(2,:),4);
QL(3,:) =  round(QL(3,:),4);

QLXY = zeros(size(x));
for i=1:b
    actual_x = x(1,i);
    for j=1:a
        actual_y = y(j,1);
        [~,nodo]= min((QL(1,:) - actual_x).^2 + (QL(2,:) - actual_y).^2);
        QLXY(j,i)= QL(3,nodo);
    end
end
