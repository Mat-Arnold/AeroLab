%% Lab 3 Script
clear, clc, close all
%% Inputs
% calibration parameters
a0 = 0.161918555250023;
a1 = 3.01930069200996;


airDensity_kgpm3 = 1.225;
airViscosity_pas = 1.789e-5;

% data
stepSize_m = ((0.5 * 0.0254)/4000) * (100);
pitotHalf_m = (0.01*0.0254)/2;

t = readtable('lab3rawdata.csv');

[~,n] = size(t);

% all different lengths so we slap em in a cell array
data_volts{1} = t.run10cm_volt(~isnan(t.run10cm_volt));
data_volts{2} = t.run20cm_volt(~isnan(t.run20cm_volt));
data_volts{3} = t.run30cm_volt(~isnan(t.run30cm_volt));
data_volts{4} = t.run40cm_volt(~isnan(t.run40cm_volt));
data_volts{5} = t.run50cm_volt(~isnan(t.run50cm_volt));

%% Calculations

% convert to pressures and velocities
for run = 1:n
    data_inH2o{run} = a0 + a1.*data_volts{run};
    data_pa{run} = (data_inH2o{run}*0.0254)*997.77*9.81;
    data_mps{run} = sqrt((2*data_pa{run})./airDensity_kgpm3);
end


% find free stream velocity for each run by averaging the final 10 readings
uInf_mps = zeros(1,5);
for run = 1:n
    uInf_mps(run) = mean(data_mps{run}(end-9:end));
end

% find the delta at which U = 0.99 U Infinity
delta_m = zeros(1,5);
for run = 1:n
    idx = find((data_mps{run} >= 0.99*uInf_mps(run)),1,'first');
    delta_m(run) = idx * stepSize_m + pitotHalf_m;
end

% find d* and theta, and H-Shape
deltaStar_m = zeros(1,5);
theta_m = zeros(1,5);
hShape= zeros(1,5);
for run = 1:n
    deltaStar_m(run) = sum((1-data_mps{run}./uInf_mps(run))*stepSize_m);
    theta_m(run)= sum((data_mps{run}./uInf_mps(run)).*(1-data_mps{run}./uInf_mps(run)).*stepSize_m);
    hShape(run) = deltaStar_m(run)/theta_m(run);
end

% perform Cf Calcs
for run = 1:n
    % shear stress method
    tau0_pa(run) = airViscosity_pas * (data_mps{run}(1)/pitotHalf_m);
    cfTau(run) = tau0_pa(run)/(0.5*airDensity_kgpm3*(uInf_mps(run)^2));

    % von karman method
    p = polyfit(0.1:0.1:0.5, theta_m,1);
    dtheta_dx = p(1);
    cfVonKarman = 2 * dtheta_dx;

    % reynold's method
    re_x(run) = (airDensity_kgpm3*uInf_mps(run)*(run/10))/airViscosity_pas;
    if run <= 1
        deltaRe_m(run) = (5*(run/10))/(sqrt(re_x(run)));
        cfRe(run) = (0.664)/sqrt(re_x(run));
    else
        deltaRe_m(run) = (0.382*(run/10))/(re_x(run)^(1/5));
        cfRe(run) = (0.592)/(re_x(run)^(1/5));
    end
end
% find delta percent difference
dDiff = abs(deltaRe_m - delta_m)./deltaRe_m * 100;


% grab needed averages
cfTau_avg = mean(cfTau);
cfRe_avg = mean(cfRe);

%% Plots
% Velocity vs Height Plots
height_m = (1:1:150) * stepSize_m + pitotHalf_m;
figure('Name','Velocities')
for run = 1:n
    plot(data_mps{run},height_m(1:length(data_mps{run})),'LineWidth',2)
    hold on, grid on;
end
% xlim([0,10.5])
xlabel('Velocity [m/s]');
ylabel('Height [m]');
title('Velocity vs Height');
legend({'10 cm','20 cm','30 cm','40 cm','50 cm'},'Location','Best')
ax = gca;
ax.FontSize = 25;

%% Make Tables
% raw data + conversions to velocities
for k = 1:n
    name = ['run_',num2str(k*10), '_cm.csv'];
    tempTable = table;
    tempTable.volts = data_volts{k};
    tempTable.inH20 = data_inH2o{k};
    tempTable.pa = data_pa{k};
    tempTable.mps = data_mps{k};
    writetable(tempTable,fullfile('Tables',name))
end

% output table
lab3output = table;
lab3output.delta_m = delta_m';
lab3output.deltaStar_m = deltaStar_m';
lab3output.theta_m = theta_m';
lab3output.hShape = hShape';
lab3output.tau0_pa = tau0_pa';
lab3output.cfTau = cfTau';
lab3output.cfRe  = cfRe';
lab3output.deltaRe_m = deltaRe_m';
lab3output.deltaDiff = dDiff';
writetable(lab3output,fullfile('Tables','lab3output.csv'))













