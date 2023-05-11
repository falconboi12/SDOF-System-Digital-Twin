function [a_e,b_e] = error_incorp(a,b,mu,sigma)
    a_e = a + randn(1,length(a))*sigma + mu;
    b_e = b + randn(1,length(b))*sigma + mu;
end
