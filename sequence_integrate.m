for ie=1:numel(eall);
    e=eall(ie);
    calc_strength;
%     ifig=ie;
    plot_strength_integrate;
    Sall(ifig,ie)=Stotal;
    Eall(ifig,ie)=Stotal.*eall(ie);
    orient portrait
    set(gcf,'PaperPosition',[1.25,3,6,5]);
    print(ifig,sprintf('%s%s E%4.1f.pdf',RunName,gname(ifig).name,log10(e)),'-dpdf');
end