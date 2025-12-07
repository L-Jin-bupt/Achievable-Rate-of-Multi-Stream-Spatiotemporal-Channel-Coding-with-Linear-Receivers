clc;
clear all;
close all;
global t H delta sigma I_r p n epsi r;
t=12;
r=20;
m=12;%B/15k
B=180*10^3;
Ts=66.7*10^-6;
fc=3.5*10^9;
dr=200;brlos=10^((-32.4-30*log10(dr)-20*log10(3.5))/10);
P=0.1*10^-3;
noi=10^-3*10^(-17.4)*B;
p=sqrt(P*brlos/noi);
delta=1;
sigma=1;
I_r=eye(r,r);
%% D=11
D=11;
temp=1;
for io=-7:1:-1
    error(temp)=10^io;
    temp=temp+1;
end
infmount=3.5*10^3;
for monte=1:1000
     H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
     ind=1;
    for epsi=error
        indiIC=0; indiNW=0; indiZF=0; indiMMSE=0;
        for n=130:-1:20
            %% random layer
            [CranIC,RranIC,CranNW,RranNW]=randomlayer(D);
            infIC=RranIC*n;
            if infIC-infmount<0 && indiIC==0
                nIC(ind,monte)=n;
                indiIC=1;
            end
            infNW=RranNW*n;
            if infNW-infmount<0 && indiNW==0
                nNW(ind,monte)=n;
                indiNW=1;
            end
      
            %% TCC
             [CMMSE,RMMSE,CZF,RZF,CU,RU]=TCC(H);
              infZF=RZF*n;
            if infZF-infmount<0 && indiZF==0
                nZF(ind,monte)=n;
                indiZF=1;
            end
            infMMSE=RMMSE*n;
            if infMMSE-infmount<0 && indiMMSE==0
                nMMSE(ind,monte)=n;
                indiMMSE=1;
            end
            if indiIC==1 && indiNW==1
                break;
            end
            if indiIC==1 && indiNW==1 && indiZF==1 && indiMMSE==1
                break;
            end
        end
        ind=ind+1;
    end
end
nICm=mean(nIC,2);
nNWm=mean(nNW,2);
nZFm=mean(nZF,2);
nMMSEm=mean(nMMSE,2);
TIC=nICm*Ts/12/10^-3;
TNW=nNWm*Ts/12/10^-3;
TZF=nZFm*Ts/12/10^-3;
TMMSE=nMMSEm*Ts/12/10^-3;

clear nIC nNW;
%% D=5
D=5;
for monte=1:1000
     H=delta/sqrt(2)*(randn(r,t)+randn(r,t)*j);
     ind=1;
    for epsi=error
        indiIC=0; indiNW=0; indiZF=0; indiMMSE=0;
        for n=120:-1:20
            %% random layer
            [CranIC,RranIC,CranNW,RranNW]=randomlayer(D);
            infIC=RranIC*n;
            if infIC-infmount<0 && indiIC==0
                nIC(ind,monte)=n;
                indiIC=1;
            end
            infNW=RranNW*n;
            if infNW-infmount<0 && indiNW==0
                nNW(ind,monte)=n;
                indiNW=1;
            end
            if indiIC==1 && indiNW==1
                break;
            end
        end
        ind=ind+1;
    end
end
nICm=mean(nIC,2);
nNWm=mean(nNW,2);
TIC2=nICm*Ts/12/10^-3;
TNW2=nNWm*Ts/12/10^-3;

error1 = logspace(-7, -1, 7);

figure; hold on;
semilogx(error1, TNW, '-p', 'Color', [205,85,85]/255, 'LineWidth', 1.5);
semilogx(error1, TNW2, ':p', 'Color', [205,85,85]/255, 'LineWidth', 1.5);
semilogx(error1, TIC, '-+', 'Color', [70,130,180]/255, 'LineWidth', 1.5);
semilogx(error1, TIC2, ':+', 'Color', [70,130,180]/255, 'LineWidth', 1.5);
semilogx(error1, TMMSE, '-o', 'Color', [205,85,85]/255, 'LineWidth', 1.5);
semilogx(error1, TZF, '-d', 'Color', [70,130,180]/255, 'LineWidth', 1.5);


set(gca, 'XScale', 'log');           
xlim([1e-7 1e-1]);                   
set(gca, 'XTick', 10.^(-7:-1));      
set(gca, 'XMinorTick','on');         
set(gca, 'XTickLabel', {'10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
set(gca, 'TickLabelInterpreter', 'tex', 'FontSize', 12);
grid on; box on;
ax = gca;
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridAlpha = 0.3; ax.MinorGridAlpha = 0.2;


xlabel('Error probability, $\epsilon$', 'Interpreter','latex', 'FontName','Times New Roman');
ylabel('Latency / ms', 'Interpreter','latex', 'FontName','Times New Roman');


lgd = legend('STCC-NW {\it D}=11','STCC-NW {\it D}=5','STCC-IC {\it D}=11','STCC-IC {\it D}=5',...
             'TCC-MMSE','TCC-ZF','Location','best');
set(lgd,'FontName','Times New Roman','Interpreter','latex');



