% Choose thermal properties
thid=-1;
while (thid<0)|(thid>2)
    disp('Available profile types');
    disp('  1: linear profile');
    disp('  2: erf profile');
    thid=input('Enter ID for the type of thermal profile');
end
if isempty(thid); thid=1; end
ans=input(sprintf('Default surface temperature is %gC; Enter chosen Ts: ',Ts-Celsius));
if ~isempty(ans)
    Ts=ans+Celsius;
end

ans=input(sprintf('Default surface T gradient is %g K/km; Enter chosen G: ',G*1000));
if ~isempty(ans)
    G=ans/1000;
end
%if thid==2;
    ans=input(sprintf('Default adiabatic temperature is %gC; Enter chosen Ti: ',Ti-Celsius))+Celsius;
    if ~isempty(ans)
        Ti=ans+Celsius;
    end
%end
del=2*(Ti-Ts)/(G*sqrt(pi));

% setup temperature
switch thid 
case 2
    Temperature=@(z)Ts+(Ti-Ts)*erf(z/del);
otherwise
    Temperature=@(z)min(Ts+z*G,Ti);
end

