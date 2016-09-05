function s=StressAdd(rh,r,T,P,g)
s=0;
for i=1:numel(rh)
    s=s+RheolSimple(rh(i),r,T,P,g);
end