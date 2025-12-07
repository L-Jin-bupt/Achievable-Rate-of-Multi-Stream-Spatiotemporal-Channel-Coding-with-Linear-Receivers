function [CequIC,RequIC,CequNW,RequNW]=equallayer(D)
global t;
CequIC=0;
RequIC=0;
CequNW=0;
RequNW=0;
if(mod(t, D)==0)
    layer=t/D*ones(D,1);
    %% IC receiver
    [CequIC,RequIC]=ICreceiver(layer);
    %% NW receover
    [CequNW,RequNW]=NWreceiver(layer);
end
