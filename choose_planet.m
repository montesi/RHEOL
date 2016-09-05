% Choose planet
ip=nplanets+1;
while (ip>nplanets)|(ip<1)
    for i=1:nplanets
        disp(sprintf('%4d: %s',i,planet(i).name))
    end
    ip=input('Enter Planet ID: ');
end
Ts=planet(ip).env.Ts;
Ti=planet(ip).env.Ti;
G=planet(ip).env.G;
P0=planet(ip).env.P0;
g=planet(ip).global.gravity;