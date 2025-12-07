function [CMMSE,RMMSE,CZF,RZF,CU,RU]=TCC(H)
global t p n epsi sigma I_r;
    %% MMSE
for q=1:t
    ht=H;
    H(:,q)=[];
    Hbar=H;
    H=ht;
    temp=sigma*I_r+p^2*Hbar*Hbar';
    SINR(q)=p^2*H(:,q)'*temp^(-1)*H(:,q);
    V(q)=1-(1+SINR(q))^(-2);
    C_T(q)=real(log2(1+SINR(q)));
    R_T(q)=C_T(q)-real(sqrt(V(q)/n)*log2(exp(1))*qfuncinv(epsi));
end
CMMSE=sum(C_T);
RMMSE=sum(R_T);
    %% ZF
for q=1:t
    temp=pinv(H'*H)*sigma;
    SINR_ZF(q)=p^2/(temp(q,q));
    V(q)=1-(1+SINR_ZF(q))^(-2);
    CZF_T(q)=real(log2(1+SINR_ZF(q)));
    RZF_T(q)=CZF_T(q)-real(sqrt(V(q)/n)*log2(exp(1))*qfuncinv(epsi));
end
CZF=sum(CZF_T);
RZF=sum(RZF_T);
%% TCC-upper bound
lambda=p^2*eig(H'*H)/sigma;
for i=1:t
    Vu(i)=1-(1+lambda(i))^(-2);
    C_Tu(i)=real(log2(1+lambda(i)));
    R_Tu(i)=C_Tu(i)-real(sqrt(Vu(i)/n)*log2(exp(1))*qfuncinv(epsi));
end
CU=sum(C_Tu);
RU=sum(R_Tu);     
end
