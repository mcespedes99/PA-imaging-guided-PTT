%% Working with NetGen: netgentest.m
%
% This example demonstrates how to import a mesh from Netgen. The python
% source code for Netgen that generates the mesh can found in the
% examples/netgen_square_with_two_circles.py. The Python source code can be
% viewed <netgen_square_with_two_circles.py here>

%% Import the NetGen file
% Netgen meshes can be imported using 'importNetGenMesh'. In addition to
% the mesh structure, it returns the regions and boundaries in the vol file
% as cell arrays. If the second argument is set to 'false', a new boundary
% will be generated and the one in the file will not be used. It is
% recommended since the original boundary oftain contains boundary elements
% that are between normal elements. This is currently not supported in
% ValoMC.

clear all;

if(exist('square_with_circle.vol', 'file') ~= 2)
    error('Could not find the mesh data file.');
end

[vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh('square_with_circle.vol', false);


%% Find indices

healthy = cell2mat(regions(find(strcmp(region_names,'healthy'))));
tumor = cell2mat(regions(find(strcmp(region_names,'tumor'))));
indices_for_lightsource = cell2mat(boundaries(find(strcmp(boundary_names,'lightsource'))));


%% Set optical parameters and light sources using the indices

%% Set optical coefficients
%For muscle:
vmcmedium.absorption_coefficient(healthy) = 0.0214336;   
vmcmedium.scattering_coefficient(healthy) = 8.2729053;  
vmcmedium.scattering_anisotropy(healthy) = 0.933957;                            
vmcmedium.refractive_index(healthy) = 1.37;  

% Use the indices
skin = findElements(vmcmesh,'rectangle',[0.77 0],[1.54],[60]);
vmcmedium.absorption_coefficient(skin) = 0.0207638;   
vmcmedium.scattering_coefficient(skin) = 5.1004766;  
vmcmedium.scattering_anisotropy(skin) = 0.715;                            
vmcmedium.refractive_index(skin) = 1.3773113; 

adipose = findElements(vmcmesh,'rectangle',[2.3 0],[2.2],[60]);
vmcmedium.absorption_coefficient(adipose) = 0.0083523;   
vmcmedium.scattering_coefficient(adipose) = 3.7088950;  
vmcmedium.scattering_anisotropy(adipose) = 0.715;                            
vmcmedium.refractive_index(adipose) = 1.44; 


vmcmedium.absorption_coefficient(tumor) = 0.0214336;   
vmcmedium.scattering_coefficient(tumor) = 16.755845;  
vmcmedium.scattering_anisotropy(tumor) = 0.933957;                            
vmcmedium.refractive_index(tumor) = 1.37; 
%For tumor, 80% has nanoparticles, in different concentrations.
r = rand;
[num,~]=size(tumor);
withnp = int16(num*8);
range_mu_a= (0.089544-0.0350556);
range_mu_s= (16.85308-16.791084);
%Matrix with modified nodes:
chosen_nodes = zeros(3,num);
for i=1:num
   r = rand;
   node = randsample(tumor,1);
   if(ismember(node,chosen_nodes))
      while(ismember(node,chosen_nodes)==1)
          node = randsample(tumor,1);
      end
   end
   chosen_nodes(1,i) = node;
   mu_a_nodo = r*range_mu_a + 0.0350556;
   mu_s_nodo = r*range_mu_s + 16.791084;
   vmcmedium.absorption_coefficient(node)=mu_a_nodo;
   vmcmedium.scattering_coefficient(node)=mu_s_nodo;
   chosen_nodes(2,i) = mu_a_nodo;
   chosen_nodes(3,i) = mu_s_nodo;
end

vmcboundary.lightsource(indices_for_lightsource) = {'direct'};   


%% Run the Monte Carlo simulation
options.photon_count=1e8;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

%% Plot the solution

figure('rend','painters','pos',[10 10 1200 400])

%h = subplot(1,2,1);
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', vmcmedium.absorption_coefficient(:), 'FaceColor', 'flat', 'EdgeColor','none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       
hold off
title('Absorption coefficient [1/mm]');

%h=subplot(1,2,2);
figure('rend','painters','pos',[10 10 1200 400]);
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       
title('Fluence [W/mm^2]');
hold off

