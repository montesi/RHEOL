for ie=1:numel(eall);
    e=eall(ie);
    calc_strength;
%     ifig=ie;
    plot_strength_integrate;
    Sall(ifig,ie)=Stotal;
    Eall(ifig,ie)=Stotal.*eall(ie);
    print(ifig,'-dpdf',sprintf('%s%s E%4.1f.pdf',RunName,gname(ifig).name,log10(e)));
end