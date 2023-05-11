function [delta_m, delta_k] = decay_fn(t, a_k, e_k, b_k, b_m, e_m)
%Generates a sawtooth function and cos curve of the form:

    if nargin <2
        a_k = 0.0004; 
        e_k = 0.05;
        b_k = 0.02; 
        b_m = 0.15;
        e_m = 0.25;
    end

    delta_m = [e_m * sawtooth(b_m*t-pi)];
    delta_k = [exp(-a_k*t).*((1+e_k*cos(b_k*t))/(1+e_k)) - 1];
end