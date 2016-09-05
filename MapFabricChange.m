load rock; nrock=numel(rock);
%%
%Choose_rock
for ir=1:nrock
    disp(sprintf('%4d: %s',ir,rock(ir).name));
end
irs=input('Desired strong mineral phase: ');
irw=input('Desired weak mineral phase: ');
C=input('Weak phase abundance');

disp('Rheologies for strong phase');
nrheol=size(rock(irs).rheol,2);
for ir=1:nrheol
    disp(sprintf('%4d: %s: %s',...
        ir,...
        rock(irs).rheol(ir).name,...
        rock(irs).rheol(ir).ref));
end
ihs=input('Desired ductile rheology for the strong phase: ');
rhs=rock(irs).rheol(ihs);
%
disp('Rheologies for weak phase');
nrheol=size(rock(irw).rheol,2);
for ir=1:nrheol
    disp(sprintf('%4d: %s: %s',...
        ir,...
        rock(irw).rheol(ir).name,...
        rock(irw).rheol(ir).ref));
end
ihw=input('Desired ductile rheology for the weak phase: ');
rhw=rock(irw).rheol(ihw);

%% Analytical Potential
P=20000*10*3000;
r=1e-15;
Tmin=400; Tmax=1200;
T=linspace(Tmin, Tmax, 100);
Tk=T+273;

ss=RheolSimple(rhs,r,Tk,P);
sw=RheolSimple(rhw,r,Tk,P);
ns=rhs.n;
nw=rhw.n;
R=ss./sw;

L=(1-C)*((1-C)+C./R).^(ns)+C*((1-C)*R+C).^nw;
Ls=C*((1-C)*R).^nw;
Lw=(1-C)*(C./R).^ns;
%%
figure(1);
clf;
plot(T,[L;Ls;Lw]); set(gca,'yscale','log'); legend('L','Ls','Lw');

title(sprintf('%g%s%s,\n%g%s%s',(1-C)*100,'% ',rhs.ref,C*100,'% ',rhw.ref),'fontsize',18);

% if max(L)>1e12; 
    axis ([Tmin,Tmax,1,1e50]);
% else
%     axis ([Tmin,Tmax,1,1e6]);
% end
set(gca,'fontsize',12);
xlabel('Temperature (°C)','fontsize',18)
ylabel('Localization potential','fontsize',18);


%% Form rheology structure
clear rh;
rh(1)=rhs;
rh(2)=rhw;
Cm=[1-C,C];
%% Localization potential
P=20000*10*3000;
r=1e-15;
Tall=[400:100:1200]
g=ones(size(Cm)); % Need to replace with input request if GSS is used
%
s=NaN(size(Tall)); rf=s;
for iT=1:numel(Tall)
    Tk=Tall(iT)+273;
    s(iT)=MixCstRateArith(Cm,rh,r,Tk,P,g); 
    rf(iT)=MixCstStressArith(Cm,rh,s(iT),Tk,P,g); 
end

figure(1); hold on; 
plot(Tall, rf./r, 'ko')

%% Setup two-level mix
%foliation parameter F
slev2G=@(F,C,rh,r,T,P,g)(1-F).*MixCstRateArith(C,rh,r,T,P,g)+...
        F.*fzero(@(s)MixCstStressArith(C,rh,s,T,P,g)-r,[0,1e10]);
rlev2G=@(F,C,rh,s,T,P,g)fzero(@(r)(slev2G(F,C,rh,r,T,P,g)-s),[1e-20,1e20]);


%% Rheology demo
T=1000; %temperature (°C)
Tk=T+273;
rall=10.^([-17:-10]); %strain rate
gn=1; %grain size (m)
g=[gn,gn];

figure (3);clf; 
subplot 211; hold on;
plot(rall,RheolSimple(rhs,rall,Tk,P,gn)/1e6,'r')
plot(rall,RheolSimple(rhw,rall,Tk,P,gn)/1e6,'color',[0.5,0,0])
sCRA=MixCstRateArith(Cm,rh,rall,Tk,P,g)
sCRG=MixCstRateGeom(Cm,rh,rall,Tk,P,g)
;
rCRA=NaN(size(sCRA));rCRG=rCRA; %Constant Rates
for ir=1:numel(sCRA);
    rCRA(ir)=fzero(@(r)(MixCstRateArith(Cm,rh,r,Tk,P,g)/sCRA(ir))-1,[0,1e-10]);
    rCRG(ir)=fzero(@(r)(MixCstRateGeom(Cm,rh,r,Tk,P,g)/sCRG(ir))-1,[0,1e-10]);
end

sCSA=NaN(size(rall)); sCSG=sCSA;%rcheck=sall;
for ir=1:numel(rall)
    sCSA(ir)=fzero(@(s)MixCstStressArith(Cm,rh,s,Tk,P,g)-rall(ir),[0,1e10]);
    sCSG(ir)=fzero(@(s)MixCstStressGeom(Cm,rh,s,Tk,P,g)-rall(ir),[0,1e10]);
% rcheck=MixCstStressArith(C,rh,sall,T,P,g)
end
% rcheck=MixCstStressArith(Cm,rh,sall,Tk,P,g)
%
plot(rall,sCRA/1e6,'g')
plot(rall,sCRG/1e6,'color',[0,0.5,0])
%
plot(rall,sCSA/1e6,'b')
plot(rall,sCSG/1e6,'color',[0,0,0.5])


% plot(rall,sall/1e6,'m')
plot(rCRA,sCRA/1e6,'go')
plot(rCRG,sCRG/1e6,'color',[0,0.5,0],'marker','o')
% plot(rcheck,sall/1e6,'mo')
% legend('Strong','Weak','Cst Rate','cst stress','check cst rate','check cst stress');
% title(sprintf([repmat('%s +',[1,numel(rh)])]),'fontsize',18);

set(gca,'fontsize',12,'yscale','log','xscale','log',...
    'xlim',[1e-16,1e-10],'ylim',[0.1,1e3],'box','on');
xlabel('Strain Rate (s^{-1})','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

subplot 212;
hold on;
Fall=[0:.25:1];
plot(rall,RheolSimple(rhs,rall,Tk,P,gn)/1e6,'k')
plot(rall,RheolSimple(rhw,rall,Tk,P,gn)/1e6,'k')

sF=NaN(numel(Fall),numel(rall));
ns=8; scheck=NaN(numel(Fall),ns); rF=sF; 
for iF=1:numel(Fall);
    for ir=1:numel(rall)
        sF(iF,ir)=slev2G(Fall(iF),Cm,rh,rall(ir),Tk,P,g);
    end
    scheck(iF,:)=linspace(min(sF(iF,:)),max(sF(iF,:)),ns);
    for is=1:ns;
        rF(iF,is)=rlev2G(Fall(iF),Cm,rh,scheck(iF,is),Tk,P,g);
    end
end
plot(rall,sF/1e6,'x-');
plot(rF',scheck'/1e6,'o');

% sM=max(sCRA); sm=min(sall);
% sa=linspace(sm,sM,10);
% rF=NaN(numel(Fall),numel(sa));
% for is=1:numel(sa);
%     for iF=1:numel(Fall);
%         rF(iF,is)=rlev2G(Fall(iF),Cm,rh,sa(is),Tk,P,g);
%     end
% end

% plot(geq*1e6,seq/1e6,'k','linewidth',1);
legend('Strong','Weak','various F');
% title(sprintf([repmat('%s +',[1,numel(rh)])]),'fontsize',18);

set(gca,'fontsize',12,'yscale','log','xscale','log',...
    'xlim',[1e-16,1e-10],'ylim',[0.1,1e3],'box','on');
xlabel('Strain Rate (s^{-1})','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);

%% Localization
r=1e-15;
Tall=[500:100:1000];
Fall=linspace(0,1,35);
si=NaN(size(Tall));
rf=NaN(numel(Tall),numel(Fall));sc=rf;
for iT=1:numel(Tall);
    Tk=Tall(iT)+273;
    si(iT)=slev2G(Fall(1),Cm,rh,r,Tk,P,g);
    for iF=1:numel(Fall);
        rf(iF,iT)=rlev2G(Fall(iF),Cm,rh,si(iT),Tk,P,g);
        sc(iF,iT)=slev2G(Fall(iF),Cm,rh,rf(iF,iT),Tk,P,g);
    end
end
%
figure(4); clf; orient tall; hold on;
subplot(211)
plot(Fall,rf);

legend(num2str(Tall'),'location','northwest')
set(gca,'fontsize',12,'yscale','log','xlim',[0,1],'ylim',10.^[-16,-5],'box','on');
xlabel('Foliated Fraction','fontsize',18);
ylabel('Strain Rate (s^{-1})','fontsize',18)
title(sprintf('%g%s%s,\n%g%s%s',(1-C)*100,'% ',rhs.ref,C*100,'% ',rhw.ref),'fontsize',18);


subplot(212)
plot(sc/1e6,rf);
% legend(num2str(Tall'))
set(gca,'fontsize',12,'xscale','log','yscale','log','xlim',10.^[0.1,3],'ylim',10.^[-16,-10],'box','on');
xlabel('Foliated Fraction','fontsize',18);
ylabel('Strain Rate (s^{-1})','fontsize',18)





% print(1,sprintf('%s%g%s.pdf',rhs.name,C,rhw.name),'-dpdf')
    

