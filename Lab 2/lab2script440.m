%% Script For Aerolab 2

clear, clc, close all

%% Raw Data
table = readtable('MAE 440 Spring 2022 Experiment 2 -02.xlsx');
angle_deg = table.Var1(2:end);
run1Hi_inH2O = table.Var2(2:end);
run1Lo_inH2O = table.Var3(2:end);
run2Hi_inH2O = table.Var6(2:end);
run2Lo_inH2O = table.Var7(2:end);

diameter_in = 2.63;
diameter_m = diameter_in*0.0254;
bias_inH2O = -0.181;

pRun1_inH2O = ((run1Hi_inH2O + run1Lo_inH2O)./2) - bias_inH2O;
pRun2_inH2O = ((run2Hi_inH2O + run2Lo_inH2O)./2) - bias_inH2O;

%% Calcs
% convert angle to radians
angle_rad = deg2rad(angle_deg);

% calculate Cp
cp1 = pRun1_inH2O./pRun1_inH2O(1);
cp2 = pRun2_inH2O./pRun2_inH2O(1);

% graph cp vs angle
figure('Name','Run 1')
plot(angle_deg, cp1)
hold on, grid on
xlabel('Angle [deg]')
ylabel('Pressure Coefficient')
title('Cp vs Angle for Run 1')
xlim([0 360])

figure('Name','Run 2')
plot(angle_deg, cp2)
hold on, grid on
xlabel('Angle [deg]')
ylabel('Pressure Coefficient')
title('Cp vs Angle for Run 2')
xlim([0 360])

% calculate the drag coefficients
cd1 = (1/2)*sum(cp1.*cos(angle_rad).*deg2rad(10));
cd2 = (1/2)*sum(cp2.*cos(angle_rad).*deg2rad(10));

% calculate the air velocity using bernouli's 
airDensity_kgpm3 = 1.225;
airViscosity_pas = 1.81e-5;

% covert inH20 to pascals using rho*g*h, convering inH20 to mH20
dp1_pa = 997*9.81*(pRun1_inH2O(1)*0.0254);
dp2_pa = 997*9.81*(pRun2_inH2O(1)*0.0254);

v1_mps = sqrt((2*dp1_pa)/airDensity_kgpm3);
v2_mps = sqrt((2*dp2_pa)/airDensity_kgpm3);

% calculate reynold's nubmer
re1 = (airDensity_kgpm3*v1_mps*diameter_m)/airViscosity_pas;
re2 = (airDensity_kgpm3*v2_mps*diameter_m)/airViscosity_pas;

%% Print Results
fprintf('\t Run 1 \t\t Run 2\n')
fprintf('Cd: \t %0.3f \t\t %0.3f\n',cd1,cd2);
fprintf('Re: \t %0.3f \t %0.3f\n',re1, re2)











