function t_pulse = estimate_pd(p, t)
% Estimate the pulse duration
% 
% INPUTS:
%     p - Pressure trace
%     t - Time trace (seconds)
% 
% OUTPUTS:
%     t_pulse - Pulse duration (seconds)
% 
% t_pulse = ESTIMATE_PD(p, t);

p_int = cumsum(p.^2);
p_int = p_int./max(p_int);
t1 = t(find(p_int >= 0.1,1));
t2 = t(find(p_int >= 0.9,1));

t_pulse = 1.25*(t2-t1);

end
