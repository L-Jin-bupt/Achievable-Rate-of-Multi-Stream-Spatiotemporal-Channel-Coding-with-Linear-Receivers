clc;
clear all;
close all;
global t H delta sigma I_r p n epsi r;
t=12;
r=20;
B=180*10^3;
m=12;%B/15k
n=96;%blocklength:m*0.534 ms/Ts
Ts=66.7*10^-6;
fc=3.5*10^9;
dr=200;brlos=10^((-32.4-30*log10(dr)-20*log10(3.5))/10);
P=0.1*10^-3;
noi=10^-3*10^(-17.4)*B;
p=sqrt(P*brlos/noi);
epsi=10^-7;
delta=1;
sigma=1;
epsi=10^-7;
I_r=eye(r,r);
D=1:1:t;
for monte=1:1000
    H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
    for x=D
    %% random layer
    [CranIC(x,monte),RranIC(x,monte),CranNW(x,monte),RranNW(x,monte)]=randomlayer(x);
    %% equal layer
   [CequIC(x,monte),RequIC(x,monte),CequNW(x,monte),RequNW(x,monte)]=equallayer(x);
    end
    %% TCC
   [CMMSE(monte),RMMSE(monte),CZF(monte),RZF(monte),CU(monte),RU(monte)]=TCC(H);
end
for x=D
    CranICm(x)=mean(CranIC(x,:))*m/Ts;
    RranICm(x)=mean(RranIC(x,:))*m/Ts;
    CranNWm(x)=mean(CranNW(x,:))*m/Ts;
    RranNWm(x)=mean(RranNW(x,:))*m/Ts;
    CequICm(x)=mean(CequIC(x,:))*m/Ts;
    RequICm(x)=mean(RequIC(x,:))*m/Ts;
    CequNWm(x)=mean(CequNW(x,:))*m/Ts;
    RequNWm(x)=mean(RequNW(x,:))*m/Ts;
end
CMMSEm=mean(CMMSE)*ones(t,1)*m/Ts;RMMSEm=mean(RMMSE)*ones(t,1)*m/Ts;
CZFm=mean(CZF)*ones(t,1)*m/Ts;RZFm=mean(RZF)*ones(t,1)*m/Ts;
CUm=mean(CU)*ones(t,1)*m/Ts;RUm=mean(RU)*ones(t,1)*m/Ts;

valid_idx = CequICm ~= 0;
x_valid = D(valid_idx);
CequICm_valid = CequICm(valid_idx);

valid_idx = RequICm ~= 0;
RequICm_valid = RequICm(valid_idx);

valid_idx = CequNWm ~= 0;
CequNWm_valid = CequNWm(valid_idx);

valid_idx = RequNWm ~= 0;
RequNWm_valid = RequNWm(valid_idx);

figure;
hold on
plot(D,CranNWm,'-p','Color',[205,85,85]/255,'LineWidth',1.5);
plot(x_valid,CequNWm_valid,':p','Color',[205,85,85]/255,'LineWidth',1.5);
plot(D,CranICm,'-+','Color',[70,130,180]/255,'LineWidth',1.5);
plot(x_valid,CequICm_valid,':+','Color',[70,130,180]/255,'LineWidth',1.5);
plot(D,CUm,'-*','Color',[205,149,12]/255,'LineWidth',1.5);
plot(D,CMMSEm,'-o','Color',[205,85,85]/255,'LineWidth',1.5);
plot(D,CZFm,'-d','Color',[70,130,180]/255,'LineWidth',1.5);
xlim([1 D(end)])
xlabel('Number of streams, \it{D}','FontName','Times New Roman');
ylabel('Capacity / bps','FontName','Times New Roman')
legend('STCC-NW-R','STCC-NW-E','STCC-IC-R','STCC-IC-E','TCC-U','TCC-MMSE','TCC-ZF','FontName','Times New Roman');
legend('location','best');
grid on;
box on;
set(gca,'DefaultAxesFontName','Times New Roman')

figure;
hold on
plot(D,RranNWm,'-p','Color',[205,85,85]/255,'LineWidth',1.5);
plot(x_valid,RequNWm_valid,':p','Color',[205,85,85]/255,'LineWidth',1.5);
plot(D,RranICm,'-+','Color',[70,130,180]/255,'LineWidth',1.5);
plot(x_valid,RequICm_valid,':+','Color',[70,130,180]/255,'LineWidth',1.5);
plot(D,RUm,'-*','Color',[205,149,12]/255,'LineWidth',1.5);
plot(D,RMMSEm,'-o','Color',[205,85,85]/255,'LineWidth',1.5);
plot(D,RZFm,'-d','Color',[70,130,180]/255,'LineWidth',1.5);
xlim([1 D(end)])
xlabel('Number of streams, \it{D}','FontName','Times New Roman');
ylabel('Achievable rate / bps','FontName','Times New Roman')
legend('STCC-NW-R','STCC-NW-E','STCC-IC-R','STCC-IC-E','TCC-U','TCC-MMSE','TCC-ZF','FontName','Times New Roman');
legend('location','best');
grid on;
box on;
set(gca,'DefaultAxesFontName','Times New Roman')