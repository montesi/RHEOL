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
ihd=input('Desired Dislocation creep rheology: ');
rhd=rock(irs).rheol(ihd);
ihg=input('Desired Diffusion creep rheology: ');
rhg=rock(irs).rheol(ihg);

disp('Piezometer: ')
npiez=size(rock(irs).piezo,2);
for ip=1:npiez
    disp(sprintf('%4d: %s',...
        ip,...
        rock(irs).piezo(ip).ref));
end
ipz=input('Desired Dislocation creep rheology: ');
pz=rock(irs).piezo(ipz).geq;

%%
P=20000*10*3000;
r=1e-15;
Tmin=400;Tmax=1200;
T=linspace(Tmin,Tmax,30);
gi=0.001;ge=10e-6;
Tk=T+273;
ss=[];L=[];geq=[];
for iT=1:numel(Tk)
    ss(iT)=fzero(@(s)RateSimple(rhd,s,Tk(iT),P,gi)+...
        RateSimple(rhg,s,Tk(iT),P,gi)-r,[1,1e20]);
    L(iT)=(RateSimple(rhd,ss(iT),Tk(iT),P,ge)+...
        RateSimple(rhg,ss(iT),Tk(iT),P,ge))/r;
    
    %geq(iT)=0.015*(ss(iT)/1e6)^(-1.33); %van der wal piezometer
    geq(iT)=pz(ss(iT));
    Leq(iT)=(RateSimple(rhd,ss(iT),Tk(iT),P,geq(iT))+...
        RateSimple(rhg,ss(iT),Tk(iT),P,geq(iT)))/r;
    
    gbd(iT)=fzero(@(g)RheolSimple(rhg,r,Tk(iT),P,g)-ss(iT),geq(iT));
    Lbd(iT)=(RateSimple(rhd,ss(iT),Tk(iT),P,gbd(iT))+...
        RateSimple(rhg,ss(iT),Tk(iT),P,gbd(iT)))/r;
end
%%
figure(1);
clf;
plot(T,L,'linewidth',2); set(gca,'yscale','log');
hold on;
plot(T,Leq,'r','linewidth',1); set(gca,'yscale','log');
plot(T,Lbd,'g','linewidth',1); set(gca,'yscale','log');
legend('Imposed','piezometer','boundary')
%legend(num2str(Qall'));
% 
title(sprintf('%s%s,\n%gto%gm',rhd.ref,rhg.ref,gi,ge),'fontsize',18);

% if max(L)>1e12; 
%     axis ([Tmin,Tmax,1,1e50]);
% else
    axis ([Tmin,Tmax,1,1e6]);
% end
set(gca,'fontsize',12);
xlabel('Temperature (°C)','fontsize',18)
ylabel('Localization potential','fontsize',18);
%%
figure(3);
clf;
plot(geq*1e6,ss/1e6,'r'); hold on
plot(gbd*1e6,ss/1e6,'g'); hold on
axis ([1,1e6,1e-3,1e3]);
set(gca,'fontsize',12,'xscale','log','yscale','log');
xlabel('Grain size (micron)','fontsize',18)
ylabel('Stress (MPa)','fontsize',18);
legend ('Equilibrium','Boundary')
%% 
figure(2);
clf;hold on;
plot(T,T*0+gi*1e6,'b','linewidth',2);
plot(T,T*0+ge*1e6,'b');
plot(T,geq*1e6,'r'); hold on
plot(T,gbd*1e6,'g'); hold on
    axis ([Tmin,Tmax,1,1e6]);
set(gca,'fontsize',12,'yscale','log');
xlabel('Temperature (°C)','fontsize',18)
ylabel('Grain size (micron)','fontsize',18);
legend ('initial','final','piezometer','boundary')

%%
print(1,sprintf('%s_%s_%gto%g.pdf',rhd.name,rhg.name,gi,ge),'-dpdf')
    

