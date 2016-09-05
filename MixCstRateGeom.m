function s=MixCstRateGeom(C,rh,r,T,P,g)
s=1;
for i=1:numel(rh)
    s=s.*((RheolSimple(rh(i),r,T,P,g(i))).^(C(i)));
end