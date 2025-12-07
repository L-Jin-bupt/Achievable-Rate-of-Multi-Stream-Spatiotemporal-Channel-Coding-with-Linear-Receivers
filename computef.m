function [Rft,Cft,Rt,Ct]=computef(D)
global t H sigma I_r p n epsi;  
for mont=1:10
    layer=zeros(D,1);   
    k=t-D+1;
    for i=1:D-1
        layer(i)=randi([1,k],1,1);
        k=t-sum(layer)-(D-i)+1;
    end
    layer(i+1)=t-sum(layer);
    
    for d=1:length(layer)
        I_l=eye(layer(d),layer(d));
        H_d=H(:,sum(layer(1:d))-layer(d)+1:sum(layer(1:d)));
        tt=H;
        H(:,sum(layer(1:d))-layer(d)+1:sum(layer(1:d)))=[];
        H_dbar=H;
        H=tt;
        K_d=sigma*I_r+p^2*H_dbar*H_dbar';
        Hda=K_d^(-1/2)*H_d;
        C(d)=real(log2(det(I_l+p^2*Hda'*Hda)));
        R(d)=C(d)-real(sqrt((layer(d)-trace((I_l+p^2*Hda'*Hda)^(-2)))/n)*log2(exp(1))*qfuncinv(epsi));
      %% orthogonal
        if (d==1)
            s=1;
        else
            s=sum(layer(1:d-1))+1;
        end 
        Cf(d)=0;
        Rlossf(d)=0;
        for k=s:sum(layer(1:d))
            Cf(d)=Cf(d)+abs(log2(1+p^2*H(:,k)'*K_d^(-1)*H(:,k)));
            Rlossf(d)=Rlossf(d)+abs(1-(1+p^2*H(:,k)'*K_d^(-1)*H(:,k))^(-2));
        end
        Rf(d)=Cf(d)-sqrt(Rlossf(d))*log2(exp(1))*qfuncinv(epsi)/sqrt(n);
    end
    Rft(mont)=sum(Rf);
    Cft(mont)=sum(Cf);
    Ct(mont)=sum(C);
    Rt(mont)=sum(R);
end
end