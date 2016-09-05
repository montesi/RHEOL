function r=MixCstStressGeom(C,rh,s,T,P,g)
r=1;
for i=1:numel(rh)
    r=r.*((RateSimple(rh(i),s,T,P,g(i))).^(C(i)));
end