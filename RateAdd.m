function r=RateAdd(rh,s,T,P,g)
r=0;
for i=1:numel(rh)
    r=r+RateSimple(rh(i),s,T,P,g);
end