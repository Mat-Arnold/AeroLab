%% lab5script440
clear, clc, close all

%% import data
rawData = readtable("rawdata.xlsx");
xd = rawData.x_d(~isnan(rawData.x_d));
delY_m = rawData.deltay_in(~isnan(rawData.deltay_in))*0.0254;
n = 1;
vel_mps{n} = rawData.xd_0_1_mps(~isnan(rawData.xd_0_1_mps)); n = n+1;
vel_mps{n} = rawData.xd_1_mps(~isnan(rawData.xd_1_mps)); n = n+1;
vel_mps{n} = rawData.xd_2_5_mps(~isnan(rawData.xd_2_5_mps)); n = n+1;
vel_mps{n} = rawData.xd_5_mps(~isnan(rawData.xd_5_mps)); n = n+1;
clear n

%% begin calcs
n = length(xd);

% preallocation block
uinf_mps = zeros(n,1);
um_mps = zeros(n,1);
umIdx = zeros(n,1);
centerIdx = zeros(n,1);

for k = 1:n
    uinf_mps(k) = mean([vel_mps{k}(1:5);vel_mps{k}(end-4:end)]);
    [um_mps(k), umIdx(k)] = max(vel_mps{k});
    
    % calculate the flow center
    upperIdx(k) = find(diff(vel_mps{k}) > 100*delY_m(k),1,'first');
    lowerIdx(k) = find(abs(diff(vel_mps{k})) > 100*delY_m(k),1,'last');
    centerIdx(k) = round(mean([upperIdx(k),lowerIdx(k)]));

    % calculate r
    r_m(k) = (lowerIdx(k)-upperIdx(k))*delY_m(k);

end