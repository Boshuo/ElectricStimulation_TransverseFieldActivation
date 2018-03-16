function J = JError_Lapicque(x,PW,I) % relative error for Lapicque S-D durve
% x(1) = t_ch  x(2) = I_rh

I_fit=x(2)./(1 - 2.^( -PW / x(1)));

J=sqrt(nanmean((abs(I_fit-I)./abs(I_fit)).^2));

if any(x < 0)
    J = NaN;
end

end