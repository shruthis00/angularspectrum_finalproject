function [isppa] = calc_isppa(p, alpha_0, alpha_true, prop_dist, t_axis)
    
    % using a linear derating scheme
    p_derated = exp(alpha_0*prop_dist) * exp(-1*alpha_true*prop_dist) * p;
    I_onaxis = (p_derated* 1e-6).^2;

    % isolate the actual pulse
    dT = mean(diff(t_axis));
    pd = 1e-6;       % cheating here, find a better way to put this in
    
    isppa = sum(I_onaxis)/(t_axis(end) - t_axis(1));

end
