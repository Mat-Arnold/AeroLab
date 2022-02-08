%% Script For Aerolab 1

clear, clc, close all

%% Raw Data

% Geometry

shape.plate.dia_m = 75.62e-3;
shape.plunger.dia_m = 75.62e-3;
shape.mushroom.dia_m = 75.62e-3;
shape.sphere.dia_m = 75.62e-3;
shape.stream.dia_m = 57.46e-3;

% results
frequencies_hz = [16 26 36 46];

fixture.force_lbf = [-0.007 -0.023 -0.046 -0.078];
fixture.velocity_mps = [5.52 9.20 13.60 16.96];

shape.sphere.force_lbf = [-0.016	-0.047	-0.094	-0.155];
shape.sphere.velocity_mps = [5.53	9.24	13.11	17.01];

shape.plate.force_lbf = [-0.027	-0.075	-0.149	-0.253];
shape.plate.velocity_mps = [5.58	9.20	13.00	16.91];

shape.mushroom.force_lbf = [-0.014	-0.042	-0.083	-0.142];
shape.mushroom.velocity_mps = [5.59	9.18	13.07	16.96];

shape.plunger.force_lbf = [-0.028	-0.080	-0.159	-0.267];
shape.plunger.velocity_mps = [5.60	9.25	13.08	16.98];

shape.stream.force_lbf = [-0.008	-0.024	-0.050	-0.084];
shape.stream.velocity_mps = [5.53	9.25	13.05	16.95];

% standard atmosphere
% from fundamentals of aerodynamics
pressure_pa = 1.01325e5;
density_kgpm3 = 1.2250;
temperaure_K = 288.16;
viscosity_kgpms = 1.7894e-5;



%% Calculations
shapeList = fieldnames(shape);
n = length(shapeList);

for k = 1:n
    name = shapeList{k};
    % adjust for fixure, and make positive, and convert to newtons
    shape.(name).adjForce_N = abs(shape.(name).force_lbf - fixture.force_lbf)*4.44822;
    
    % calculate cross-sectional area (all objects are circular)
    shape.(name).area_m2 = (shape.(name).dia_m/2)^2 * pi;
    
    % calculate the drag coefficient (lab manual)
    shape.(name).cd = shape.(name).adjForce_N./((1/2) * density_kgpm3 * (shape.(name).velocity_mps.^2) * shape.(name).area_m2);
    
    % calculate average cd
    shape.(name).avgCd = mean(shape.(name).cd);
    
    % calculate the reynolds number (fundamentals of aerodynamics)
    % using diameter as our characteristic length 
    shape.(name).re = (density_kgpm3 * shape.(name).velocity_mps * shape.(name).dia_m)./viscosity_kgpms;
end

%% Results

fprintf('\t\tSphere \t\tPlate \t\tMushroom \tPlunger \tStreamline\n')
fprintf('Area (mm^2): \t%0.2f \t\t%0.2f \t\t%0.2f \t\t%0.2f \t\t%0.2f \n', shape.sphere.area_m2, shape.plate.area_m2, shape.mushroom.area_m2, shape.plunger.area_m2, shape.stream.area_m2);
fprintf('Average Cd: \t%0.2f \t\t%0.2f \t\t%0.2f \t\t%0.2f \t\t%0.2f \n', shape.sphere.avgCd, shape.plate.avgCd, shape.mushroom.avgCd, shape.plunger.avgCd, shape.stream.avgCd);

figure('Name','Cd vs Re');
for k = 1:n
   name = shapeList{k};
   plot(shape.(name).re, shape.(name).cd, '*-', 'LineWidth',2)
   hold on
end
grid on
ax = gca;
ax.FontSize = 25;
xlabel('Reynolds Number');
ylabel('Drag Coefficient');
title('Drag Coefficient vs Reynolds Number')
legend({'Flat Plate','Plunger','Mushroom','Sphere','Streamlined'},'Location','Best')

