clc;
clear all;
close all;
global t H delta sigma I_r p n epsi r;
B=180*10^3;
m=12;%B/15k
n=96;%blocklength:m*0.534 ms/Ts
fc=3.5*10^9;
dr=200;brlos=10^((-32.4-30*log10(dr)-20*log10(3.5))/10);
P=10*10^-3;
noi=10^-3*10^(-17.4)*B;
p=sqrt(P*brlos/noi);
epsi=10^-7;
delta=1;
sigma=1;
rx=20:20:200;
t=4;
D=3;
for x=1:length(rx)
    r=rx(x);
    I_r=eye(r,r);
    Ccount(x)=0;
    Rcount(x)=0;
    for monte=1:1000
        H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
        [RfNW(x,:,monte),CfNW(x,:,monte),RNW(x,:,monte),CNW(x,:,monte)]=computef(D);
        %% TCC
       [CMMSE(x,monte),RMMSE(x,monte),CZF(x,monte),RZF(x,monte),CU(x,monte),RU(x,monte)]=TCC(H);
       Ctemp=CNW(x,:,monte)-CMMSE(x,monte);
       Ccount(x)=Ccount(x)+sum(Ctemp>0);
       Rtemp=RNW(x,:,monte)-RMMSE(x,monte);
       Rcount(x)=Rcount(x)+sum(Rtemp>0);
       disp(t);disp(x);disp(monte);
    end
    CP(x)=Ccount(x)/(size(CNW,2)*size(CNW,3));
    RP(x)=Rcount(x)/(size(CNW,2)*size(CNW,3));
    P(x)=CP(x)*RP(x);
    Ptemp=RNW(x,:,:)./RfNW(x,:,:);
    PH(x)=mean(Ptemp(:));
        Ctemp=CNW(x,:,:)./CfNW(x,:,:);
    CH(x)=mean(Ctemp(:));
end
Ps1=P;PHs1=PH;

t=12;
D=11;
for x=1:length(rx)
    r=rx(x);
    I_r=eye(r,r);
    Ccount(x)=0;
    Rcount(x)=0;
    for monte=1:1000
        H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
        [RfNW(x,:,monte),CfNW(x,:,monte),RNW(x,:,monte),CNW(x,:,monte)]=computef(D);
        %% TCC
       [CMMSE(x,monte),RMMSE(x,monte),CZF(x,monte),RZF(x,monte),CU(x,monte),RU(x,monte)]=TCC(H);
       Ctemp=CNW(x,:,monte)-CMMSE(x,monte);
       Ccount(x)=Ccount(x)+sum(Ctemp>0);
       Rtemp=RNW(x,:,monte)-RMMSE(x,monte);
       Rcount(x)=Rcount(x)+sum(Rtemp>0);
       disp(t);disp(x);disp(monte);
    end
    CP(x)=Ccount(x)/(size(CNW,2)*size(CNW,3));
    RP(x)=Rcount(x)/(size(CNW,2)*size(CNW,3));
    P(x)=CP(x)*RP(x);
    Ptemp=RNW(x,:,:)./RfNW(x,:,:);
    PH(x)=mean(Ptemp(:));
end
Ps=P;PHs=PH;



t=32;
D=31;
for x=1:length(rx)
    r=rx(x);
    I_r=eye(r,r);
    Ccount(x)=0;
    Rcount(x)=0;
    for monte=1:1000
        H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
        [RfNW(x,:,monte),CfNW(x,:,monte),RNW(x,:,monte),CNW(x,:,monte)]=computef(D);
        %% TCC
       [CMMSE(x,monte),RMMSE(x,monte),CZF(x,monte),RZF(x,monte),CU(x,monte),RU(x,monte)]=TCC(H);
       Ctemp=CNW(x,:,monte)-CMMSE(x,monte);
       Ccount(x)=Ccount(x)+sum(Ctemp>0);
       Rtemp=RNW(x,:,monte)-RMMSE(x,monte);
       Rcount(x)=Rcount(x)+sum(Rtemp>0);
       disp(t);disp(x);disp(monte);
    end
    CP(x)=Ccount(x)/(size(CNW,2)*size(CNW,3));
    RP(x)=Rcount(x)/(size(CNW,2)*size(CNW,3));
    P(x)=CP(x)*RP(x);
    Ptemp=RNW(x,:,:)./RfNW(x,:,:);
    PH(x)=mean(Ptemp(:));
        Ctemp=CNW(x,:,:)./CfNW(x,:,:);
    CH(x)=mean(Ctemp(:));
end
Ps2=P;PHs2=PH;



rx=rx.';
figure;
hold on
set(0,'DefaultAxesFontName','Times New Roman')
plot(rx,Ps1,'o','Color',[70,130,180]/255,'LineWidth',1.5);
plot(rx,Ps,'k+','LineWidth',1.5);
plot(rx,Ps2,'-m','LineWidth',1.5);
plot(rx,PHs1,'-*','Color',[205,85,85]/255,'LineWidth',1.5);
plot(rx,PHs,'-p','Color',[205,149,12]/255,'LineWidth',1.5);
plot(rx,PHs2,'-cx','LineWidth',1.5);
 xlim([20 200])
xlabel('Number of receive antennas, \it{r}');
ylabel('Proportion')
legend('{\it P_CP_R} , {\it t}=4, {\it D}=3','{\it P_CP_R} , {\it t}=12, {\it D}=11','{\it P_CP_R} , {\it t}=32, {\it D}=31',...
    '{\it P_H} ,{\it t}=4, {\it D}=3','\it{P_H} , {\it t}=12, {\it D}=11', '{\it P_H} , {\it t}=32, {\it D}=31');
legend('location','best');
grid on;
box on;
