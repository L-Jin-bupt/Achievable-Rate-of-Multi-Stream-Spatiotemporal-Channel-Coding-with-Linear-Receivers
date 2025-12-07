clc;
clear all;
close all;
global t H delta sigma I_r p n epsi r;
t=12;
r=20;
m=12;%B/15k
n=96;%blocklength:m*0.534 ms/Ts
epsi=10^-7;
Ts=66.7*10^-6;
delta=1;
I_r=eye(r,r);
D=[5];
PdB=10:-2:-10;
for i=1:length(PdB)
    P(i)=10^(PdB(i)/10);
end
for monte=1:1000
    H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
    for x=1:length(P)
        p=sqrt(P(x));
        sigma=1;
        %% random layer
        [CranIC(x,monte,1),RranIC(x,monte,1),CranNW(x,monte,1),RranNW(x,monte,1)]=randomlayer(D(1));
        %% equal layer
%         [CequIC(x,monte,1),RequIC(x,monte,1),CequNW(x,monte,1),RequNW(x,monte,1)]=equallayer(D(1));
        %% TCC
       [CMMSE(x,monte),RMMSE(x,monte),CZF(x,monte),RZF(x,monte),CU(x,monte),RU(x,monte)]=TCC(H);
       sigma=1+p^2*0.2;
        %% random layer
        [imCranIC(x,monte,1),imRranIC(x,monte,1),imCranNW(x,monte,1),imRranNW(x,monte,1)]=randomlayer(D(1));
        %% TCC
       [imCMMSE(x,monte),imRMMSE(x,monte),imCZF(x,monte),imRZF(x,monte),imCU(x,monte),imRU(x,monte)]=TCC(H);
    end
end
for i=1:length(D)
    for x=1:length(P)
        RranICm(x,i)=mean(RranIC(x,:,i))*m/Ts;
        RranNWm(x,i)=mean(RranNW(x,:,i))*m/Ts;
        RMMSEm(x,i)=mean(RMMSE(x,:))*m/Ts;
        RZFm(x,i)=mean(RZF(x,:))*m/Ts;
        imRranICm(x,i)=mean(imRranIC(x,:,i))*m/Ts;
        imRranNWm(x,i)=mean(imRranNW(x,:,i))*m/Ts;
        imRMMSEm(x,i)=mean(imRMMSE(x,:))*m/Ts;
        imRZFm(x,i)=mean(imRZF(x,:))*m/Ts;
    end
end



figure;
hold on;

plot(PdB, RranNWm(:,1), '-','Color',[205,85,85]/255,'LineWidth',1.5);
plot(PdB, imRranNWm(:,1), 'p','Color',[205,85,85]/255,'LineWidth',1.5);
% plot(PdB, RranNWm(:,2), ':p','Color',[205,85,85]/255,'LineWidth',1.5);
plot(PdB, RranICm(:,1), '-','Color',[70,130,180]/255,'LineWidth',1.5);
plot(PdB, imRranICm(:,1), '+','Color',[70,130,180]/255,'LineWidth',1.5);
% plot(PdB, RranICm(:,2), ':+','Color',[70,130,180]/255,'LineWidth',1.5);
plot(PdB, RMMSEm(:,1), ':','Color',[205,85,85]/255,'LineWidth',1.5);
plot(PdB, imRMMSEm(:,1), 'o','Color',[205,85,85]/255,'LineWidth',1.5);
plot(PdB, RZFm(:,1), ':','Color',[70,130,180]/255,'LineWidth',1.5);
plot(PdB, imRZFm(:,1), 'd','Color',[70,130,180]/255,'LineWidth',1.5);

xlabel('Transmit SNR, {\it p}^2/\sigma^2 (dB)','FontName','Times New Roman');
ylabel('Achievable rate / bps','FontName','Times New Roman');
legend('STCC-NW perfect CSI','STCC-NW \sigma_i^2=0.2','STCC-IC perfect CSI','STCC-IC \sigma_i^2=0.2',...
       'TCC-MMSE perfect CSI','TCC-MMSE \sigma_i^2=0.2','TCC-ZF perfect CSI','TCC-ZF \sigma_i^2=0.2','FontName','Times New Roman','Location','best');


grid on; box on;


