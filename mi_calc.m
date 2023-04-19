clear
clc
close all

load('source_pressure_4C1_30V_M1_9_100MHz_single.mat')
% 
% figure()
% for i = 1:100
%     imagesc(elPos, latPos, squeeze(pressure(i,:,:)));
%     colorbar()
%     title(t(i))
%     caxis([0 10000])
% 
%     pause(0.05)
% end

%% water
% how to calculate focal depth from the surface pressure, analytically?
loadfile_water = load("pressure_focal_dist_alpha_0.00217_f0_3.5MHz_.mat");

loadfile = loadfile_water
p = loadfile.apaz;
prop_dist = loadfile.prop_dist;
fc = loadfile.f0;
alpha0 = 0.3;
actual_alpha = loadfile.alpha0;
t_max = 7e-6;

% construct the taxis bc you didn't do it before you dumbass
dT = ((1/fc)/3)/5;
taxis = 0:dT:t_max+dT;

% how to calc MI from the ending surface pressure frame
x_midind = round(length(p(:,1,1))/2);
y_midind = round(length(p(1,:,1))/2);

p_waveform = squeeze(p(x_midind, y_midind, :));

% derate pressure using linear derating scheme
% in a way that you can then impose the attenuation of
% tissue to the same system.
p_derate = derate(p_waveform, alpha0, fc, prop_dist);
pnp_derate = abs(min(p_derate));

% ok don't derate it
p_not_derated = p_waveform;

c0 = 1480;
rho = 1000;

pd = estimate_pd(p_derate, taxis);
isppa = compute_isppa(p_derate,taxis, c0, rho);

% divide by freq to get MI
mi = (pnp_derate/(1e6)) / sqrt(fc/(1e6));

figure()
visual_fft(p_derate, dT, fc)
hold on

%% tissue, alpha = 0.3
loadfile_tissue = load("pressure_focal_dist_alpha_0.3_f0_3.5MHz_.mat");

loadfile = loadfile_tissue;
p = loadfile.apaz;
prop_dist = loadfile.prop_dist;
fc = loadfile.f0;
alpha0 = loadfile.alpha0;

if ~exist('loadfile.taxis', 'var')
    dT = ((1/fc)/3)/5;
    taxis = 0:dT:(7e-6)+dT;
end

% how to calc MI from the ending surface pressure frame
x_midind = round(length(p(:,1,1))/2);
y_midind = round(length(p(1,:,1))/2);

p_waveform = squeeze(p(x_midind, y_midind, :));

% dont need to derate, we are assuming that this calc was done with
% alpha = 0.3
p_derate = p_waveform;
pnp_derate = abs(min(p_derate));

c0 = 1540;
rho = 1000;

%pd = estimate_pd(p_derate, taxis);
%isppa = compute_isppa(p_derate,taxis, c0, rho)

% divide by freq to get MI
mi = (pnp_derate/(1e6)) / sqrt(fc/(1e6));
y = 9;

visual_fft(p_derate, dT, fc)
legend('\alpha = 2.17e-3 but derated', '\alpha = 0.3')

%% for alpha = 0..?
loadfile = load('pressure_focal_dist_alpha_0_f0_3MHz_.mat');

p = loadfile.apaz;
prop_dist = loadfile.prop_dist;
fc = loadfile.f0;
alpha0 = loadfile.alpha0;
dT = mean(diff(loadfile.taxis));

p_derate = p_waveform;
visual_fft(p_derate, dT, fc)


%% TO DO:

% simulate alpha values 0.3, 2.17e-3, and ..0?
% at frequencies f = 2.5, 3.5, and 4 MHz

% time delay undo and re-calculate - impose your own time delays for focal
% depth?



%% p waveform to look at and observe higher harmonics

function visual_fft(p, dT, fc)
    % for now we assume that
    
    p_fft = abs(fftshift(fft(p)));
    p_fft = p_fft/max(p_fft);
    n_fft = length(p_fft);
    fs = 1/dT;
    p_freq = fs*1e-6/2*linspace(0,1,n_fft/2+1);
    
    plot(p_freq,p_fft(round(n_fft/2):end))
    xlabel("Pressure Frequency (MHz)");
    ylabel("Normalized FT Magnitude")
    title("Pulse Waveform Harmonics, Derated")
    subtitle_str = ['f_c = ' num2str(fc*1e-6) 'MHz'];
    subtitle(subtitle_str)

end

function visual_pulse(p,taxis)
   
    p = p*1e-6;
    plot(taxis,p)
    ylim([0,])
    
end



