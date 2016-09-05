function Stotal=IntegrateStrength(model,e,thid,Temperature,rock,did);

nlayer=numel(model);
calc_strength;

nstep=100;
Stotal=0;
in=0;
for il=1:nlayer;
    for im=1:model(il).nrock;
        for is=1:model(il).rock(im).nstr;
            in=in+1;
            if model(il).rock(im).str(is).law<0; %brittle law
                z=[model(il).rock(im).str(is).ztop,model(il).rock(im).str(is).zbot];
            else
                z=linspace(model(il).rock(im).str(is).ztop,model(il).rock(im).str(is).zbot,nstep);
            end
            stress=model(il).rock(im).str(is).s(z,e);
            Stotal=Stotal+(sum(stress)-(stress(1)+stress(end))/2)*(z(end)-z(1))/(numel(z)-1);
        end
    end
end
disp(sprintf('Total strength: %g',Stotal));