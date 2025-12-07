function [CranIC,RranIC,CranNW,RranNW]=randomlayer(D)
global t;  
if(D==1)
    layer=[t];
    %% IC receiver
    [CranIC,RranIC]=ICreceiver(layer);
    %% NW receiver
    [CranNW,RranNW]=NWreceiver(layer);
end
if(D>1)
    for mont=1:100
        layer=zeros(D,1);   
        k=t-D+1;
        for i=1:D-1
            layer(i)=randi([1,k],1,1);
            k=t-sum(layer)-(D-i)+1;
        end
        layer(i+1)=t-sum(layer);
        %% IC receiver
        [CranICt(mont),RranICt(mont)]=ICreceiver(layer);
        %% NW receiver
        [CranNWt(mont),RranNWt(mont)]=NWreceiver(layer);
    end
    CranIC=mean(CranICt);
    RranIC=mean(RranICt);
    CranNW=mean(CranNWt);
    RranNW=mean(RranNWt);
end
end