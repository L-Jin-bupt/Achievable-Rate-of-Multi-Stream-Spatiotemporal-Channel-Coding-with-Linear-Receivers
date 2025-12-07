function [Ct,Rt]=ICreceiver(layer)
global t H p n epsi r sigma;
for d=1:length(layer)
    I_l=eye(layer(d),layer(d));
    H_d=H(:,sum(layer(1:d))-layer(d)+1:sum(layer(1:d)));
    tt=H;
    H(:,sum(layer(1:d))-layer(d)+1:sum(layer(1:d)))=[];
    H_dbar=H;
    H=tt;
    [U,S,V] = svd(H_dbar);
    U0=U(:,((t-layer(d)))+1:r);
    M=U0'*H_d*H_d'*U0;
    [V, Lamb] = eig(M);      
    [lambda_sorted, idx] = sort(diag(Lamb), 'descend');  
    V_ld = V(:, idx(1:layer(d)));  
    A=U0*V_ld;
    Hdabar=A'*H_dbar;
    Hda=A'*H_d;
    C(d)=real(log2(det(I_l+p^2*Hda'*Hda/sigma)));
    R(d)=C(d)-real(sqrt((layer(d)-trace((I_l+p^2*Hda'*Hda/sigma)^(-2)))/n)*log2(exp(1))*qfuncinv(epsi));
end
Ct=sum(C);
Rt=sum(R);
end
