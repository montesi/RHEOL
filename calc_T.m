function T=calc_T(z,Ts,Ti,G,del,thid);
% T=calc_T(z,Ts,Ti,G,del,thid);
% Compute temperature from geothermal parameters
%   thid: 1 for linear; 2 for erf; 3 for exponential;
%      G: surface geotherm
%    del: coefficient for erf
%     Ts: surface temperature
%     Ti: Temperature at infinity (for erf)
switch thid
case 2
    T=Ts+(Ti-Ts)*erf(z/del);
otherwise
    T=Ts+z*G;
end