load rock; nrock=numel(rock);
%%
%Choose_rock
for ir=1:nrock
    disp(sprintf('%4d: %s',ir,rock(ir).name));
end
irs=input('Rock Type: ');

disp('Rheologies:');
nrheol=size(rock(irs).rheol,2);
for ir=1:nrheol
    disp(sprintf('%4d: %s: %s',...
        ir,...
        rock(irs).rheol(ir).name,...
        rock(irs).rheol(ir).ref));
end
irh=input('Desired rheologies ');
for ir=1:numel(irh)
    rh(ir)=rock(irs).rheol(irh(ir));
end

disp('Piezometer: ')
npiez=size(rock(irs).piezo,2);
for ip=1:npiez
    disp(sprintf('%4d: %s',...
        ip,...
        rock(irs).piezo(ip).ref));
end

ipz=input('Desired piezometer: ');
pz=rock(irs).piezo(ipz).geq;
%%
gm=10.^(-6);gM=10.^(-2); %grain size (m)
sm=fzero(@(s)pz(s)-gm,[1e-20,1e20]);
sM=fzero(@(s)pz(s)-gM,[1e-20,1e20]);
seq=linspace(sm,sM,10);geq=seq*0;
for is=1:numel(seq)
    geq(is)=pz(seq(is));
end
%% Localization plot
P=20000*10*3000; %pressure (Pa)
Tmin=400; Tmax=1200; nT=100;
Tall=linspace(Tmin, Tmax, nT);
gi=0.01; %initial grain size
gf=10e-6; %final grain size (when imposed)
si=NaN(size(Tall));rf=si;ge=si;re=si;
ri=1e-15;
ssat=1e20;
for iT=1:numel(Tall);
    Tk=Tall(iT)+273;
    if  RateAdd(rh,ssat,Tk,P,gi)>=ri;
        si(iT)=fzero(@(s)RateAdd(rh,s,Tk,P,gi)-ri,[0,ssat]);
        rf(iT)=RateAdd(rh,si(iT),Tk,P,gf);
        ge(iT)=pz(si(iT));
        re(iT)=RateAdd(rh,si(iT),Tk,P,ge(iT));
    end
end
%%
figure(5); clf; 
%orient tall;
% subplot 311; hold on
hold on;
plot([Tmin,Tmax],ri*[1,1],'k')
plot(Tall,[rf;re]);

legend(sprintf('Initial grain size %4.1e micron',gi*1e6),...
    sprintf('Final grain size %4.1e micron',gf*1e6),...
    sprintf('%s piezometer',rock(irs).piezo(ipz).ref))
title(sprintf(repmat('%s +',[1,numel(rh)]),rh.name)...
    ,'fontsize',18);
set(gca,'fontsize',12,'box','on','yscale','log',...
    'xlim',[400,1200],'ylim',ri*[1,1e9]);
xlabel('Temperature (C)','fontsize',18)
ylabel('Strain rate (s^{-1})','fontsize',18);

figure(6)
subplot 211; hold on
plot([Tmin,Tmax],gi*[1,1]*1e6,'k');
% plot([Tmin,Tmax],gf*[1,1]*1e6,'b');
plot(Tall,[[Tall*0+gf];ge]*1e6)
% legend(sprintf('Initial grain size %4.1e micron',gi*1e6),...
%     sprintf('Final grain size %4.1e micron',gf*1e6),...
%     sprintf('%s piezometer',rock(irs).piezo(ipz).ref))
title(sprintf(repmat('%s +',[1,numel(rh)]),rh.name)...
    ,'fontsize',18);
set(gca,'fontsize',12,'box','on','yscale','log',...
    'xlim',[Tmin,Tmax],'ylim',[1e-0;1e5]);
xlabel('Temperature (C)','fontsize',18)
ylabel('Grain size (microns)','fontsize',18);

subplot 212; hold on
plot(Tall,si/1e6,'k');
% plot([Tmin,Tmax],gf*[1,1]*1e6,'b');
% plot(Tall,[[Tall*0+gf];ge]*1e6)
% legend(sprintf('Initial grain size %4.1e micron',gi*1e6),...
%     sprintf('Final grain size %4.1e micron',gf*1e6),...
%     sprintf('%s piezometer',rock(irs).piezo(ipz).ref))
title(sprintf(repmat('%s +',[1,numel(rh)]),rh.name)...
    ,'fontsize',18);
set(gca,'fontsize',12,'box','on','yscale','log',...
    'xlim',[Tmin,Tmax],'ylim',[1e-1;1e4]);
xlabel('Temperature (C)','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);


%%
P=20000*10*3000; %pressure (Pa)
T=800; %Temperature (°C)
Tk=T+273; %Temperature (K)
rall=10.^([-15:-12]); %strain rate (1/s)
gall=10.^([-6:0.1:-2]); %grain size (m)
ss=NaN(numel(rall),numel(gall));
for ir=1:numel(rall);
    r=rall(ir); %current strain rate
    for ig=1:numel(gall)
        gn=gall(ig);
        ss(ir,ig)=fzero(@(s)RateAdd(rh,s,Tk,P,gn)-r,[0,1e10]);
    end
end
%
figure(1);
clf;
plot(gall*1e6,ss/1e6,'linewidth',2); set(gca,'yscale','log');
hold on;
plot(geq*1e6,seq/1e6,'k','linewidth',1)
legend(num2str(rall'));
title(sprintf([repmat('%s +',[1,numel(rh)]),'\n%gK'],rh.name,T),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale','log','ylim',[1e-1,1e4]);
xlabel('Grain size (micron)','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

%%
P=20000*10*3000; %pressure (Pa)
Tall=[400:100:1000]; %Temperature (°C)
r=10.^(-15); %strain rate (1/s)
gall=10.^([-6:0.1:-2]); %grain size (m)
ss=NaN(numel(Tall),numel(gall));
for iT=1:numel(Tall);
    Tk=Tall(iT)+273; %Temperature (K)
    for ig=1:numel(gall)
        gn=gall(ig);
        ss(iT,ig)=fzero(@(s)RateAdd(rh,s,Tk,P,gn)-r,[1e-20,1e20]);
    end
end
%
figure(2);
clf;
plot(gall*1e6,ss/1e6,'linewidth',2); set(gca,'yscale','log');
hold on;
plot(geq*1e6,seq/1e6,'k','linewidth',1);
legend(num2str(Tall'));
title(sprintf([repmat('%s +',[1,numel(rh)]),'\n%g s^{-1}'],rh.name,r),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale','log','ylim',[1e-1,1e4]);
xlabel('Grain size (micron)','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

%%
P=20000*10*3000; %pressure (Pa)
T=800; %Temperature (°C)
Tk=T+273;
rall=10.^([-16:0.25:-8]); %strain rate (1/s)
gall=[1,2,5,10,20,50,100,200,500,1000]/1e6;%10.^([-6:-2]); %grain size (m)
ss=NaN(numel(gall),numel(rall));
for ig=1:numel(gall);
    gn=gall(ig);
    for ir=1:numel(rall)
        ss(ig,ir)=fzero(@(s)RateAdd(rh,s,Tk,P,gn)-rall(ir),[0,1e20]);
    end
end
%
figure(7);
clf;
plot(rall,ss/1e6,'linewidth',2); set(gca,'yscale','log');
hold on;
% plot(geq*1e6,seq/1e6,'k','linewidth',1);
legend(num2str(gall'*1e6),'location','northwest');
title(sprintf([repmat('%s +',[1,numel(rh)]),'\n%g s^{-1}'],rh.name,r),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale',...
    'log','ylim',[1e-1,1e4]);
xlabel('Strain rate (s^{-1})','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

%%
P=20000*10*3000; %pressure (Pa)
T=800; %Temperature (°C)
Tk=T+273; %Temperature (K)
sall=10.^[5:10]; %stress (Pa)
gall=10.^([-6:0.1:-2]); %grain size (m)
rs=NaN(numel(sall),numel(gall));
for is=1:numel(sall);
    rs(is,:)=RateAdd(rh,sall(is),Tk,P,gall);
end
%
figure(3);
clf;
plot(gall*1e6,rs,'linewidth',2); set(gca,'yscale','log');
hold on;
plot(geq*1e6,RateAdd(rh,seq,Tk,P,geq),'k','linewidth',1)
legend(num2str((sall/1e6)'));
title(sprintf([repmat('%s +',[1,numel(rh)]),'\n%g K'],rh.name,T),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale','log','ylim',[1e-17,1e-10]);
xlabel('Grain size (micron)','fontsize',18)
ylabel('Strain rate (s^{-1})','fontsize',18);

%%
P=20000*10*3000; %pressure (Pa)
Tall=[400:100:1000]; %Temperature (°C)
ri=10.^(-15); %initial strain rate (1/s)
gi=1e-2; %initial grain size (m)
gall=10.^(linspace(-6,log10(gi),30)); %grain size (m)
rs=NaN(numel(Tall),numel(gall));
si=NaN(size(Tall)); gf=si; rf=si;
for iT=1:numel(Tall);
    Tk=Tall(iT)+273; %Temperature (K)
    if  RateAdd(rh,ssat,Tk,P,gi)>=ri;
        si(iT)=fzero(@(s)RateAdd(rh,s,Tk,P,gi)-ri,[1e-20,1e20]); %stress level
        rs(iT,:)=RateAdd(rh,si(iT),Tk,P,gall);
        gf(iT)=pz(si(iT));
        rf(iT)=RateAdd(rh,si(iT),Tk,P,gf(iT));
    end
end
%%
figure(4);
clf; orient tall
subplot 211
plot(gall*1e6,rs,'linewidth',2); set(gca,'yscale','log');
hold on;
%plot(geq*1e6,RateAdd(rh,seq,Tk,P,geq),'k','linewidth',1)
plot(gf*1e6,rf,'k-o','linewidth',1)
legend(num2str(Tall'),'location','NorthEast');
% title(sprintf([repmat('%s +',[1,numel(rh)])],rh.name),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale','log',...
    'xlim',10.^[0,4],'ylim',[1e-16,1e-8]);
xlabel('Grain size (micron)','fontsize',18)
ylabel('Strain rate (s^{-1})','fontsize',18);
%
subplot 212
plot([repmat(gi,1,numel(si));gf]*1e6,repmat(si/1e6,2,1),...
    'linewidth',2); set(gca,'yscale','log');
hold on;
%plot(geq*1e6,RateAdd(rh,seq,Tk,P,geq),'k','linewidth',1)
plot(geq*1e6,seq/1e6,'k','linewidth',1)
plot(gf*1e6,si/1e6,'ko','linewidth',1)
% legend(num2str(Tall'));
% title(sprintf([repmat('%s +',[1,numel(rh)])],rh.name),'fontsize',18);

set(gca,'fontsize',12,'xscale','log','yscale','log',...
    'xlim',10.^[0,4],'ylim',[1e-1,1e4],'xlim',[1,1e4]);
xlabel('Grain size (micron)','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

%%%
%print(1,sprintf('%s_%s_%gto%g.pdf',rhd.name,rhg.name,gi,ge),'-dpdf')
    

