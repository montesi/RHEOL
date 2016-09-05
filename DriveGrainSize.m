modelSafe=model;

RunName='Run';
ans=input(sprintf('Enter run name (default: %s)',RunName));
if ~isempty(ans);
    RunName=ans;
end
eall=10.^[-18:0.5:-6];
% eall=10.^[-18:3:-9];
Sall=zeros(size(eall)+[2,0]);

ifig=1;
gname(ifig).name='Piezometer';
sequence_integrate

ifig=2; 
for inow=1:numel(model)
model(inow).rock.gs=1e-5; %repeat with 10 micron grain size
end
gname(ifig).name='10 microns';
sequence_integrate

ifig=3; 
for inow=1:numel(model)
model(inow).rock.gs=1e-2; %repeat with 10 mm grain size
end
gname(ifig).name='10 mm';
sequence_integrate;
Eall=Sall.*repmat(eall,[3,1]); %energy
%%

ifig=ifig+1;
figure(ifig); clf;
subplot(211)
hold on;
loglog(Sall/1e6/1e3,eall,'linewidth',2);
et=10.^[-18:-15]; %initial strain rate
ii=size(Sall,1); %code for initial condition (here: last calculated)
ep=10.^[linspace(-18,-9,20)]; %strain rate for display
%initialize final strain rates
ef=zeros(size(Sall,1),size(et,2)); %For constant energy
es=ef; %for constant integrated stress

for it=1:numel(et);
    st=interp1(eall,Sall(ii,:),et(it)); %initial integrated stress
    Et=st.*et(it); %initial energy
    for ia=1:size(Eall,1); %final sets
        ef(ia,it)=interp1(Eall(ia,:),eall,Et); %final strain rate / energy
        es(ia,it)=interp1(Sall(ia,:),eall,st); %final strain rate / stress
    end
    loglog(Et./ep/1e6/1e3,ep,'k','linewidth',1);
end

legend(gname.name,'Location','NorthWest')
set(get(gca,'xlabel'),'string','Integrated Stress (MPa*km)','fontSize',18);
set(get(gca,'ylabel'),'string','Strain Rate (s^{-1})','fontSize',18);
set(gca,...
    'fontSize',12,...
    'xscale','log',...
    'yscale','log',...
    'xlim',[1e-0,1e3]*100,...
    'ylim',10.^[-18,-6],...
    'box','on');

subplot(212); hold on;
loglog(et,ef./repmat(et,[3,1]),'linewidth',2);
loglog(et,es./repmat(et,[3,1]),'linewidth',1);
set(gca,...
    'fontSize',12,...
    'xscale','log',...
    'yscale','log',...'xlim',[1e-0,1e3]*100,...
    'xlim',10.^[-18,-15],...
    'box','on');

legend(gname.name,'Location','NorthWest')
set(get(gca,'xlabel'),'string','Strain Rate (s^{-1})','fontSize',18);
set(get(gca,'ylabel'),'string','Strain rate enhancement','fontSize',18);

orient tall;
print(ifig,'-dpdf',sprintf('Summ%s',RunName));

model=modelSafe;


