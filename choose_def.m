%Parameter for type of deforamtion 
disp(sprintf('%4d: %s',1,'Compression'));
disp(sprintf('%4d: %s',2,'Extension'));
disp(sprintf('%4d: %s',3,'Strike-slip'));
did=input('Enter style of deformation: ');
e=1e-15;
ans=input(sprintf('Strain rate (default is %g/s): ',e));
if ~isempty (ans)
    e=ans;
end