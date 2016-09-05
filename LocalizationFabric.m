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

%%
P=20000*10*3000;
r=1e-15;
Tmin=400;Tmax=1200;
T=linspace(Tmin,Tmax,100);
Tk=T+273;
%
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

%%
% print(1,sprintf('%s%g%s.pdf',rhs.name,C,rhw.name),'-dpdf')
    

