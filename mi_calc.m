clear
clc
close all

load('source_pressure_4C1_30V_M1_9_100MHz_single.mat')

figure()
for i = 1:100
    imagesc(elPos, latPos, squeeze(pressure(i,:,:)));
    colorbar()
    title(t(i))
    caxis([0 10000])


    pause(0.05)
end

% how to calculate focal depth of the surface pressure
