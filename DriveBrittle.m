modelSafe=model;

RunName='Run';
ans=input(sprintf('Enter run name (default: %s)',RunName));
if ~isempty(ans);
    RunName=ans;
end
eall=10.^[-20:0.5:-9];
% eall=10.^[-18:3:-6];
Wkall=[1,0.5,0.25,0.1]; nWk=numel(Wkall);
Sall=zeros(nWk,size(eall,2));
Eall=Sall;

%% 
for ifig=1:nWk;
    for irock=1:1; %limited to 1st layer
        model(irock).rock.Wk=Wkall(ifig);
    end
    %gname(ifig).name=sprintf('%4.2f strength',model.rock.Wk)
    gname(ifig).name=sprintf('%4.2f strength',Wkall(ifig));
    sequence_integrate;
end

%%

ifig=ifig+1;
figure(ifig); clf;
subplot(211)
hold on;
loglog(Sall/1e6/1e3,eall,'linewidth',2);
et=10.^[-18:-15]; %initial strain rate
ii=1;%size(Sall,1); %code for initial condition (here: first calculated)
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
    'ylim',10.^[-20,-10],...
    'box','on');

subplot(212); hold on;
loglog(et,ef'./repmat(et',[1,nWk]),'linewidth',2);
loglog(et,es'./repmat(et',[1,nWk]),'linewidth',1);
set(gca,...
    'fontSize',12,...
    'xscale','log',...
    'yscale','log',...'xlim',[1e-0,1e3]*100,...
    'xlim',10.^[-20,-15],...
    'box','on');

legend(gname.name,'Location','NorthWest')
set(get(gca,'xlabel'),'string','Strain Rate (s^{-1})','fontSize',18);
set(get(gca,'ylabel'),'string','Strain rate enhancement','fontSize',18);
%
orient tall;
print(ifig,'-dpdf',sprintf('Summ%s',RunName));

model=modelSafe;


