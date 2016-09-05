para.A=[1.1e5, 3.98e7, 6.5e3, 6.31e4, 5.7e11]; %preexponential factor MPa^{-n}/s or 1/s
para.Q=[530,370,400,445, 535]; %activation energy kJ/mol
para.n=[3.5,1,3.5,2.9,2]; %stress exponent
para.m=[0,3,2,0.7,0]; %grain size exponent
para.P=[0,0,0,0,8500]; %Peierls stress
R=8.32; %Gas constant, J/K

for i=1:numel(para.A);
    if para.P(i)~=0; % exponential creep
        rh(i).r=@(s,T)para.A(i)...
            .*exp(-(para.Q(i)*1000./(R*T)).*(1-s./para.P(i)).^(para.n(i)));
        rh(i).type='Exponential creep';
    elseif para.m(i)==0; %grain size insensitive creep;
        rh(i).r=@(s,T)para.A(i).*exp(-para.Q(i)*1000./(R*T))...
            .*(s.^para.n(i));
        rh(i).type='Dislocation creep';
    else %grain size sensitive creep
        rh(i).r=@(s,d,T)para.A(i).*exp(-para.Q(i)*1000./(R*T))...
            .*(s.^para.n(i)).*(d.^(-1*para.m(i)));
        if para.n(i)==1;
            rh(i).type='Diffusion creep';
        else
            rh(i).type='disGBS';
        end
    end
end
%%
rhChoice=[1,2,3];
r=@(s,d,T) 0;
%construct combined rheology
for i=1:numel(rhChoice)
    ir=rhChoice(i);
    if strcmpi(rh(ir).type,'Exponential creep')|strcmpi(rh(ir).type,'Dislocation creep');
        %         disp('s,T');
        %         plot(rh(ir).r(s,T),s)
        r=@(s,d,T)r(s,d,T)+rh(ir).r(s,T);
    else
        %         disp('D,G');
        %         plot(rh(ir).r(s,d,T),s)
        r=@(s,d,T)r(s,d,T)+rh(ir).r(s,d,T);
    end
end
% implicit definitions
Tr=@(s,d,rn)bisection(@(T)(r(s,d,T)./rn)-1,[1e3])
dr=@(s,T,rn)bisection(@(d)(r(s,d,T)./rn)-1,100)
sr=@(d,T,rn)bisection(@(s)(r(s,d,T)./rn)-1,10)

%% evaluate for various conditions
s=10.^(linspace(0,4,30));
dall=10.^[-1:5]; %microns
T=900+273.15;
figure(1);clf; hold on;
for id=1:numel(dall);
    dn=dall(id);
    rn=nan(numel(rhChoice)+1,numel(s));
    for i=1:numel(rhChoice)
        ir=rhChoice(i);
        if strcmpi(rh(ir).type,'Exponential creep')|strcmpi(rh(ir).type,'Dislocation creep');
            rn(i,:)=rh(ir).r(s,T);
        else
            rn(i,:)=rh(ir).r(s,dn,T);
        end
        plot(rn,s);
    end
    plot(r(s,dn,T),s,'k','linewidth',2); %combined rheology
end
legend(rh(rhChoice).type,'Location','southeast')
set(gca,'xscale','log','yscale','log','box','on','fontSize',12)
set(gca,'xlim',10.^[-18,-10]);
xlabel('Strain rate (s^{-1})','fontSize',18)
ylabel('Stress (MPa)','fontSize',18)

%%
figure(2); clf; hold on;
Tall=400:100:1000;
dall=10.^linspace(1,4,30);
rn=1e-15;
sn=nan(numel(dall),numel(Tall));
for iT=1:numel(Tall);
    Tn=Tall(iT);
    for id=1:numel(dall)
        sn(id,iT)=sr(dall(id),Tn,rn);
    end
end
plot(dall,sn');

    
legend(num2str(Tall'))
set(gca,'xscale','log','yscale','log','box','on','fontSize',12)
set(gca,'xlim',10.^[1,4],'ylim',10.^[-3,3]);
xlabel('Grain size (\mu m)','fontSize',18)
ylabel('Stress (MPa)','fontSize',18)

    