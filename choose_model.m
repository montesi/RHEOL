%Choose_model
il=0;
Pbot=P0;ztop=0;zbot=0;
newlayer='yes';
while or(strncmpi(newlayer,'y',1),newlayer==1);
    il=il+1;
    disp(sprintf('For layer %d: ',il));
    for ir=1:nrock
        disp(sprintf('%4d: %s',ir,rock(ir).name));
    end
    model(il).irock=1;
    ans=input('Desired rock type (can be a vector): ');
    if ~isempty(ans)
        model(il).irock=ans;
    end
    model(il).nrock=length(model(il).irock);
    model(il).Ci=ones(1,model(il).nrock)/model(il).nrock; %Concentration initialize
    if model(il).nrock>1;
        ans=input('Enter the volume proportion of each rock type');
        if ~isempty(ans)
            model(il).Ci=ans/sum(ans);
        end
    end
    if ~isempty(model(il).irock);
        model(il).thick=input('Layer thickness (in km): ')*1000;
        model(il).ztop=zbot;
        zbot=zbot+model(il).thick;
        model(il).zbot=zbot;
        
        %ibrit=[];
        %disp('  1: Byerlee, low P branch');
        %disp('  2: Byerlee, high P branch');
        %disp('  3: Tensile strength');
        %ans=input('Desired brittle rheology (can be a vector): ');
        %if ~isempty(ans)
        %    ibrit=ans;
        %end
        ibrit=[1,2];
        if model(il).irock==10;
            ans=input('Serpentine detected: Is it lizardite (low F, default = no): ');
            if ~isempty(ans);
                if or(strncmpi(ans,'yes',1),ans==1);
                    ibrit=[3,4];
                end
            end
        end
        
        model(il).pf='p';
        disp('Default pore fluid pressure: hydrostatic');
        ans=input('Enter pore fluid pressure: lambda or keep hydrostatic');
        if ~isempty(ans)
            model(il).pf=ans;
        end
                
%         disp('Do you want to saturate at 300 MPa?');
        ans=input('Do you want to saturate at 300 MPa? (default is no)');
        if ~isempty(ans)
            if or(strncmpi(ans,'yes',1),ans==1);
                ibrit=[ibrit,5];
            end
        end
        
        
        
        model(il).Ptop=Pbot;
        rhoav=0;
        for im=1:model(il).nrock;
            rhoav=rhoav+rock(model(il).irock(im)).density*model(il).Ci(im);
        end        
        if model(il).pf=='p'
            model(il).rhog=(rhoav-1e3)*...
                planet(ip).global.gravity;
        else
            model(il).rhog=rhoav*(1-model(il).pf)*...
                planet(ip).global.gravity;
        end
        Pbot=model(il).Ptop+model(il).thick*model(il).rhog;    
        model(il).Pbot=Pbot;


        iduct=[];
        model(il).irheol=[];
        for im=1:model(il).nrock
            nrheol=size(rock(model(il).irock(im)).rheol,2);
            for ir=1:nrheol
                disp(sprintf('%4d: %s: %s',...
                    ir,...
                    rock(model(il).irock(im)).rheol(ir).name,...
                    rock(model(il).irock(im)).rheol(ir).ref));
            end
            ans=input('Desired ductile rheology (can be a vector; enter 0 for brittle only): ');
            if ~isempty(ans)
                if ans==0;
                    nduct=0;
                    iduct=[];
                else
                    nduct=length(ans);
                    iduct=ans;
                end
                
                %model(il).gs(im)=0;
                model(il).rock(im).gs=0;
                model(il).rock(im).gdep=zeros([nrock,nduct+2]);
                for ir=1:nduct
                    if (rock(model(il).irock(im)).rheol(iduct(ir)).m~=0);
                        model(il).rock(im).gdep(ir)=1;%(model(il).gs(im)==0);
                        model(il).rock(im).gs=10e-3; %model(il).gs(im)=100e-6;
%                     if (rock(model(il).irock(im)).rheol(iduct(ir)).m~=0)&...
%                             (model(il).rock(im).gs==0);%(model(il).gs(im)==0);
%                         model(il).rock(im).gs=10e-3; %model(il).gs(im)=100e-6;
                    end
                end
                if model(il).rock(im).gs~=0;
                    %disp(sprintf('Detecting grain-size-dependent laws; default grain size: %g microns',model(il).rock(im).gs*1e6));
                    disp(sprintf('Detecting grain-size-dependent laws\n default behavior: grain size fixed at %g microns',...
                        model(il).rock(im).gs*1e6));
                    npiez=numel(rock(model(il).irock(im)).piezo);
                    if npiez==0;
                        ans=input(sprintf(...
                            'Enter new grain size for layer %d (in m):',il));
                        if ~isempty(ans);
                            model(il).rock(im).gs=ans;%model(il).gs(im)=ans;
                        end
                    else
                        ans=input(sprintf(...
                            'Enter 0 for choosing a piezometer or new grain size for layer %d (in m):',il));
                        if ~isempty(ans);
                            if ans~=0;
                                model(il).rock(im).gs=ans;%model(il).gs(im)=ans;
                            else
                                for ipiez=1:npiez;
                                    disp(sprintf('%2g: %s',ipiez,rock(model(il).irock(im)).piezo(ipiez).ref))
                                end
                                ans=input('Enter piezometer choice:');
                                if ~isempty(ans);
                                    if (ans<=npiez)&(ans>0);
                                        % imposes limits to 1 micron and 1 meter
                                        model(il).rock(im).gs=@(s)min(max(rock(model(il).irock(im)).piezo(ans).geq(s),1e-6),1);%model(il).gs(im)=ans;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
            model(il).rock(im).nrheol=nduct+length(ibrit);
            %model(il).irheol(im,[1:model(il).nrheol(im)])=[iduct(im,[1:nduct(im)]),-ibrit];
               
            model(il).rock(im).irheol=[iduct,-ibrit];
            model(il).rock(im).Wk=1;
        end
                      
        model(il).rock(im).Ty=NaN;
        if ~isempty(find(ibrit==0));
            disp(sprintf('Default tensile strength: Ty=%g MPa',model(il).rock(im).Ty));
            ans=input('Enter desired tensile strength (in MPa): ');
            if ~isempty(ans)
                model(il).rock(im).Ty=ans*1e6;
            end
        end
    end
    newlayer=input('Do you want to add a new layer?');
end
nlayer=size(model,2);
% model.nlayer=nlayer;