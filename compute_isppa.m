function [isppa,tpulse] = compute_isppa(p_derated, t, c, rho)
% Compute the spatial-peak pulse-average intensity
% 
% INPUTS:
%     p_derated - Derated pressure trace (Pa)
%     t - Time trace (seconds)
%     c - Speed of sound (m/s)
%     rho - Density (kg/m^3)
% 
% OUTPUTS:
%     isppa - ISPPA (W/cm^2)
%
% isppa = COMPUTE_ISPPA(p_derated, t);

if ~exist('c','var'), c = 1480; end
if ~exist('rho','var'), rho = 1000; end

tpulse = estimate_pd(p_derated, t);
isppa  = trapz(t, p_derated.^2)./(rho*c)/(100^2)/tpulse; 