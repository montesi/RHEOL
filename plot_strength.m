%plot_strength
set(0,'DefaultAxesLineWidth',1,...
    'DefaultAxesFontSize',12,...
    'DefaultAxesColor','none');
collist=[0,0,0;1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1;...
    0.75,0.25,0.25;0.25,0.75,0.25;0.25,0.25,0.75];
ncolor=10;
icolor=0;
nz=100;

figure (ifig)
clf
hold on

axst=gca;% plotting strength envelope
hold on
axth=axes; % plotting temperature profile
hold on
axgs=axes; % plotting grain size profile
hold on
% axab=axes;
% hold on

labels=[];

nstep=100;
Stotal=0;
%plot strength and a-b
icolor=0;in=0;
for il=1:nlayer;
    for im=1:model(il).nrock
        for is=1:model(il).rock(im).nstr;
            in=in+1;
            if model(il).rock(im).str(is).law<0; %brittle law
                z=[model(il).rock(im).str(is).ztop,model(il).rock(im).str(is).zbot];
                labels(in).text='Brittle';
            else
                z=linspace(model(il).rock(im).str(is).ztop,model(il).rock(im).str(is).zbot,nstep);
                labels(in).text=rock(model(il).irock(im)).rheol(model(il).rock(im).str(is).law).name;
            end
            if isa(model(il).rock(im).gs,'function_handle'); 
                %can't use vector for z
                stress=z*0;
                grain=stress;
                for iz=1:numel(z)
                    stress(iz)=model(il).rock(im).str(is).s(z(iz),e);
                    grain(iz)=model(il).rock(im).gs(stress(iz));
                end
            else %use vector for z
                stress=model(il).rock(im).str(is).s(z,e);
                grain=model(il).rock(im).gs;
            end
            Stotal=Stotal+(sum(stress)-(stress(1)+stress(end))/2)*(z(end)-z(1))/(numel(z)-1);
            
            
            icolor=icolor+1;
            plot(stress/1e6,z/1000,'color',collist(icolor,:),'Parent',axst);
            plot(grain*1e6,z/1000,'color',collist(icolor,:),'Parent',axgs);
            
            %plot thermal profile
            z=linspace(model(il).rock(im).str(is).ztop,model(il).rock(im).str(is).zbot,100);
            T=model(il).Temperature(z)-Celsius;
            %         ztop=model(1).str(1).ztop;
            %         zbot=model(end).scomb(end).zbot;
            %         z=ztop+(zbot-ztop)*[0:nz-1]/(nz-1);
            %         T=calc_T(z,Ts,Ti,G,del,thid)-Celsius;
            plot(T,z/1000,'color',collist(icolor,:),'linewidth',2,'parent',axth);
        end
        

        %     z0=model(il).ztop;
        %     Ptop=model(il).Ptop;
        %     rhog=model(il).rhog;
        %     for icomb=1:model(il).ncomb;
        %         ztop=model(il).scomb(icomb).ztop;
        %         zbot=model(il).scomb(icomb).zbot;
        %         z=ztop+(zbot-ztop)*[0:nz-1]'/(nz-1);
        %         stress=z*0;
        %         S=[];f=[];B=[];n=[];QRn=[];
        %         for im=1:model(il).nrock;
        %             irh=model(il).irheol(im,model(il).scomb(icomb).laws(im));
        %             if irh<0 %brittle law
        %                 [S,f]=byerlee(did,-irh);
        %                 stress=stress+...
        %                     model(il).Ci(im)*...
        %                     calc_stress(z,1,irh,...
        %                         S,f,z0,Ptop,rhog,...
        %                         B,n,m,p,QRn,e,Ts,Ti,G,del,thid);
        %             else % ductile law
        %                 n=rock(model(il).irock(im)).rheol(irh).n;
        %                 m=rock(model(il).irock(im)).rheol(irh).m;
        %                 p=rock(model(il).irock(im)).rheol(irh).p;
        %                 B=(rock(model(il).irock(im)).rheol(irh).A*...
        %                     model(il).gs(im)^(-m)).^(-1/n);
        % %                 2.^((1+(3/n))/2).*...
        % %                     (rock(model(il).irock(im)).rheol(irh).A*...
        % %                     model(il).gs(im)^(-m)).^(-1/n);
        %                 QRn=rock(model(il).irock(im)).rheol(irh).QR/n;
        %                 for iz=1:nz;
        %                     stress(iz)=stress(iz)+...
        %                     model(il).Ci(im)*...
        %                     calc_stress(z(iz),1,irh,...
        %                         S,f,z0,Ptop,rhog,...
        %                         B,n,m,p,QRn,e,Ts,Ti,G,del,thid);
        %                 end
        %             end
        %         end
        %         icolor=icolor+1;
        %         plot(stress/1e6,z/1000,'color',collist(icolor,:),'Parent',axst);
    end
end


% 
% %plot a-b
% plot([0,0],[ztop,zbot]/1000,'k');
% Tbot=calc_T(zbot,Ts,Ti,G,del,thid)-Celsius;
% zlast=ztop;
% Tlast=calc_T(zlast,Ts,Ti,G,del,thid)-Celsius;
% zt=spline(T,z,150);
% if zt>zlast
%     if zt>zbot
%         plot(-0.01-[Tlast-150,Tbot-150]*1.5e-4,[zlast,zbot]/1000,'parent',axab);
%     else
%         plot(-0.01-[Tlast-150,0]*1.5e-4,[zlast,zt]/1000,'parent',axab);
%         zlast=zt;
%         Tlast=calc_T(zlast,Ts,Ti,G,del,thid)-Celsius;
%     end
% end
% zt=spline(T,z,400);
% if (zt>zlast)
%     if zt>zbot
%         plot(-0.01*[1,1],[zlast,zbot]/1000,'parent',axab);
%     else
%         plot(-0.01*[1,1],[zlast,zt]/1000,'parent',axab);
%         zlast=zt;
%         plot(-0.01+[0,Tbot-400]*4e-4,[zlast,zbot]/1000,'parent',axab);
%     end
% end

axstpos=get(axst,'Position');
axthpos=axstpos;
axthpos(3)=axstpos(3)/2;
axstpos(1)=axthpos(1)+axthpos(3);
axstpos(3)=axthpos(3);
axgspos=axstpos;
% axabpos=axthpos;

set (axst,'Position',axstpos,...
    'Ydir','reverse',...
    'Xdir','normal',...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'Xcolor','k',...
    'Ycolor','k',...
    'TickDir','in',...
    'YTickLabel','');
set (axth,'Position',axthpos,...
    'Ydir','reverse',...
    'Xdir','reverse',...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xcolor','k',...
    'Ycolor','k',...
    'TickDir','in');
set (axgs,'Position',axthpos,...
    'Ydir','reverse',...%    'Xdir','reverse',...
    'xscale','log',...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none',...
    'Xcolor','k',...
    'Ycolor','k',...
    'TickDir','in');
% set (axab,'Position',axabpos,...
%     'Ydir','reverse',...
%     'Xdir','reverse',...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'Xcolor','k',...
%     'Ycolor','k',...
%     'TickDir','in',...
%     'YTickLabel','');

set(get(axst,'xlabel'),'string','Stress (MPa)');
set(get(axst,'ylabel'),'string','Depth (km)');
set(get(axth,'xlabel'),'string','Temperature (\circ C)');
set(get(axgs,'xlabel'),'string','Grain size (\mum)');
% set(get(axab,'xlabel'),'string','a-b');


legend(labels.text,'Location','NorthWest')

orient landscape