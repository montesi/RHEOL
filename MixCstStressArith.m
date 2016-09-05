function r=MixCstStressArith(C,rh,s,T,P,g)
r=0;
for i=1:numel(rh)
    r=r+C(i)*RateSimple(rh(i),s,T,P,g(i));
end