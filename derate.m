function [p_derated] = derate(p, alpha_0, fc, prop_dist)
% Compute the spatial-peak pulse-average intensity
% 
% INPUTS:
%   alpha_0 = attenuation to apply (dB/cm/MHz)
%   fc = center frequency (Hz)
%   prop_dist = propagation distance at pulse (m)
% 
% OUTPUTS:
%   p_derated - P (Pa)
%
% p_derated = derate(p_waveform, alpha, fc, prop_dist);

np_to_dB = 0.115;

alpha_np = alpha_0 * np_to_dB; % need alpha in Np
prop_dist = prop_dist * 1e2;   % need it in cm
fc = fc * 1e-6; % need fc in MHz

% using a linear derating scheme
p_derated = exp(-alpha_np*prop_dist*fc) .* p;

end
