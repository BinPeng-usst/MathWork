clear all;clc;close all;
%%Preparation of the ambinet test data.S is the original signal,TS is the truncated&detrended signal,NS is SVD-processed signal,FS is EKF filtered signal
SamFreq=200;Rsl=4e-3;% Ambient test parameters
ImDataFile1='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-1_20141023172144190.mat'
ImDataFile4='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-4_20141023172144190.mat'
ImDataFile7='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-7_20141023172144190.mat'
importdata(ImDataFile1);
S1=ans.Data;
importdata(ImDataFile4);
S4=ans.Data;
importdata(ImDataFile7);
S7=ans.Data;
%delete small voltage turbulents
j=1;
for i=1:size(S1)
    if abs(S1(i))>Rsl
        TS1(j,1)=S1(i);
        j=j+1;
    end
end
j=1;
for i=1:size(S4)
    if abs(S4(i))>Rsl
        TS4(j,1)=S4(i);
        j=j+1;
    end
end
j=1;
for i=1:size(S7)
    if abs(S7(i))>Rsl
        TS7(j,1)=S7(i);
        j=j+1;
    end
end
Num=min([size(TS1,1);size(TS4,1);size(TS7,1)]);

%detrend
DtrdOdr=2;
t=[0:1/SamFreq:(Num-1)/SamFreq]';
TS1=dtrend(TS1(1:Num)-polyval(polyfit(t,TS1(1:Num),DtrdOdr),t));
TS4=dtrend(TS4(1:Num)-polyval(polyfit(t,TS4(1:Num),DtrdOdr),t));
TS7=dtrend(TS7(1:Num)-polyval(polyfit(t,TS7(1:Num),DtrdOdr),t));
% SVD decompositon
NS=[TS1,TS4,TS7];
[U,e,V] = svd(NS,0);
NS=U(:,1)*e(1,1)*V(1,:);
NS=(NS(:,1)+NS(:,2)+NS(:,3))/3;

%% Preparation of structural matricies
ImMassFile='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-1_20141023172144190.mat'
ImDampFile='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-4_20141023172144190.mat'
ImStiffnessFile='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-7_20141023172144190.mat'
ImOutputMFile='C:\Users\Bin Peng\OneDrive\�������������Σ�\����\����ԭʼ����\W2\W2��1-1��\AI1-7_20141023172144190.mat'
importdata(ImMassFile);
M=ans.Data; 
importdata(ImDampFile);
K=ans.Data;
importdata(ImStiffnessFile);
C=ans.Data; 
importdata(ImOutputMFile);
D=ans.Data; 

%% EKF construction
% Omiga=14.55*(2*pi);Amp=mean(abs(TS7));Phase=0;% priors
% StateTranF=@(T) (T+sign(T)*sign((1-(T/Amp)^2))*(Amp*Omiga*sqrt(sign((1-(T/Amp)^2))*(1-(T/Amp)^2)))*(1/SamFreq));
DeltaT=1/SamFreq;

Phi=[0,1;(-1)*(M^-1)*K,(-1)*(M^-1)*C];
Psi=[1;(-1)*M];

A=expm(Phi*DeltaT);
B=(Phi^-1)*(1-(expm(Phi*DeltaT))^-1)*Psi;

StateTranF=@(X,U) A*X+B*U
MeasurementF=@(Y) D*Y;
obj = extendedKalmanFilter(StateTranF,MeasurementF,[NS(1)],'StateCovariance',5);
obj.ProcessNoise=0.618;
obj.MeasurementNoise=1; 

for k = 1:size(NS)
  [CorrectedState,CorrectedStateCovariance] = correct(obj,NS(k)); 
  [PredictedState,PredictedStateCovariance] = predict(obj);
  FS(k,1)=PredictedState;
end
% FS=polyval(polyfit(t,FS,DtrdOdr),t);

%% PSD of the filtered signal
[Psd0,f0]=pwelch(TS7,[],[],1024,SamFreq);%Hamming����Ĭ�ϴ����ȡ��ص����Ⱥ�DFT����
[Psd1,f1]=pwelch(NS,[],[],1024,SamFreq);%Hamming����Ĭ�ϴ����ȡ��ص����Ⱥ�DFT����
[Psd2,f2]=pwelch(FS,[],[],1024,SamFreq);%Hamming����Ĭ�ϴ����ȡ��ص����Ⱥ�DFT����
Psd0(1:20)=0;
Psd1(1:20)=0;
Psd2(1:20)=0;

%% Output
figure(1)
plot(t,TS7); hold on;
plot(t,NS); hold on;
plot(t,FS); 
legend('TS','NS','FS');
xlabel('\fontname{����}ʱ��\fontname{Times new Roman}(s)','FontSize',10);
ylabel('\fontname{����}���ٶ�\fontname{Times new Roman}(cm/s^{2})','FontSize',10);

% ax=gca;outerpos=ax.OuterPosition;ti=ax.TightInset; 
% left=outerpos(1)+ti(1);bottom=outerpos(2)+ti(2);ax_width=outerpos(3)-ti(1)-ti(3);ax_height=outerpos(4)-ti(2)-ti(4);
% ax.Position=[left bottom ax_width ax_height];
set(gca,'FontName','Times new Roman','FontSize',11);
set(gcf,'Units','centimeters','Position',[0 0 16 16],'Resize','off');
RFile='C:\Users\Bin Peng\Desktop\T';
print('-f1',RFile,'-painters','-dmeta','-r600');

figure(2);
[pks0,loc0]=max(Psd0);BsFreq0=f0(loc0);
[pks1,loc1]=max(Psd1);BsFreq1=f1(loc1);
[pks2,loc2]=max(Psd2);BsFreq2=f2(loc2);
plot(f0,Psd0);hold on;
plot(f1,Psd1);hold on;
plot(f2,Psd2); 
legend(['PSD of TS',' ',num2str(BsFreq0)],['PSD of NS',' ',num2str(BsFreq1)],['PSD of FS',' ',num2str(BsFreq2)]);
