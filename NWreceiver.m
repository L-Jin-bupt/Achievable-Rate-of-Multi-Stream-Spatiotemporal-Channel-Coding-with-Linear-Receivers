function [Ct,Rt]=NWreceiver(layer)
global  H sigma I_r p n epsi;
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
end
Ct=sum(C);
Rt=sum(R);
end
