function [S,f]=byerlee(did,branch);
%[S,f]=byerlee(did,branch);
%gives coefficient for SII=S+f Pov
%    branch= 1: Byerlee low pressure 
%    branch= 2: Byerlee high pressure
%    branch= 3: tensile strength
%      did = 1: Compression
%      did = 2: Extension
%      did = 4: Strike-slip

switch branch
case 1 %low pressure branch
    C=0;
    mu=0.85;
    %convert for invariants
    Ci=C/sqrt(1+mu^2);
    fi=mu/sqrt(1+mu^2);
case 2 %high pressure branch
    C=50e6;
    mu=0.6;
    %convert for invariants
    Ci=C/sqrt(1+mu^2);
    fi=mu/sqrt(1+mu^2);
case 3 %Ice friction
    C=8.3e6;
    mu=0.2;
    %convert for invariants
    Ci=C/sqrt(1+mu^2);
    fi=mu/sqrt(1+mu^2);
case 4 %Tension
    Ci=1;
    fi=1;
end

switch did
case 1% compression
    S=Ci/(1-fi);
    f=fi/(1-fi);
case 2% extension
    S=Ci/(1+fi);
    f=fi/(1+fi);
case 3%strike-slip
    S=Ci;
    f=fi;
end