%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rheol.m
% 
% Laurent Montesi, 01/06/2015
% Strength profile  with evolving grain size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load planet data
load planet
% load rock data
load ./rock.mat
% Some constant
Celsius=273.15;
%%

% Define model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choose_planet
choose_th
choose_def
choose_model

%%
for il=1:nlayer
    rhog=model(il).rhog;
    nrock=model(il).nrock;
    model(il).thid=thid;
    model(il).Temperature=Temperature;
    %model(il).zbot=model(il).ztop+model(il).thick;
    Pressure=@(z)model(il).Ptop+(z-model(il).ztop)*model(il).rhog;
    model(il).Pressure=Pressure;
    %    model(il).str=[]; %initialize stratification;
   
    for im=1:nrock
        ztop=model(il).ztop;
        zbot=model(il).zbot;
        nrheol=model(il).rock(im).nrheol;
        model(il).rock(im).str=[];        
        %% initialize laws 
        rh=[];Euse=[]
        for ir=1:nrheol
            irh(ir)=model(il).rock(im).irheol(ir);
            if irh(ir)<0 %brittle law
                [S(ir),f(ir)]=byerlee(did,-irh(ir));
                rh(ir).s=@(z,e,g)S(ir)+f(ir)*Pressure(z);
                rh(ir).Drate=@(z,s,g)0; %dislocation rate;
            else % ductile law
                rh(ir).s=@(z,e,g)RheolSimple(rock(model(il).irock(im)).rheol(irh(ir)),e,...
                    Temperature(z),Pressure(z),g);
                Euse(ir)=or(or(or(...
                    rock(model(il).irock(im)).rheol(irh(ir)).type==1,...
                    rock(model(il).irock(im)).rheol(irh(ir)).type==2),...
                    rock(model(il).irock(im)).rheol(irh(ir)).type==5),...
                    rock(model(il).irock(im)).rheol(irh(ir)).type==6);
                rh(ir).Drate=@(z,s,g)RheolSimple(rock(model(il).irock(im)).rheol(irh(ir)),s,...
                    Temperature(z),Pressure(z),g);*Euse(ir); %dislocation rate;
            end
        end
        npz=numel(rock(model(il).irock(im)).piezo);
        ipz=npz+1;
        while (ipz>npz)|(ipz<0)
            disp(sprint('%4d: $s',0,'Final grain size'));
            for i=1:npz
                disp(sprintf('%4d: %s',i,rock(model(il).irock(im)).piezo(i).ref));
            end
            ipz=input('Enter constraint on final grain size: ');
        end
        if ipz==0
            gf=10e-6;
            disp(sprintf('Default final grain size: %g microns',gf*1e6));
            ans=input(sprintf(...
                'Enter desired final grain size or leave blank for default (in m):',il));
            if ~isempty(ans);
                gf=ans;%model(il).gs(im)=ans;
            end
            gevolve=@(z,s,g)sum(rh(ir).Drate(z,s,g)).*(gf-g(z))/2;
        else
            gevolve=@(z,s,g)sum(rh(ir).Drate(z,s,g)).*(rock(model(il).irock(im)).piezo(ipz).geq(s)-g(z))/2;
        end
    end
end

%%
nz=100;
z=linspace(ztop,zbot,nz);
g=z*0+model(il).rock(im).gs;
nrh=numel(rh);
s=NaN(nrh,nz);
for ir=1:numel(rh);
    s(ir,:)=rh(ir).s(z,e,g);
end
%%
figure(2); clf;
plot(min(s)/1e6,z/1e3); set(gca,'ydir','reverse');
xlim([0,1000]);
xlabel('Stress (MPa)');
ylabel('Depth (km)');




        % Define mechanical transitions; return
        %        model(il).str(im,is).law=irt;
        %        model(il).str(im,is).ztop=zn;
        %        model(il).str(im,is).zbot=zbot;
%         [model(il).rock(im).str,model(il).rock(im).nstr]=...
%             Define_Layers(rh,model(il).rock(im).irheol,...
%             ztop,zbot,e,Pressure,Temperature);
%%


%%
% calc_strength
% ifig=1; plot_strength