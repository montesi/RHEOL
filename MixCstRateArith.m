function s=MixCstRateArith(C,rh,r,T,P,g)
s=0;
for i=1:numel(rh)
    s=s+C(i)*RheolSimple(rh(i),r,T,P,g(i));
end