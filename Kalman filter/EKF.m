clear all;clc;close all;
%%Preparation of the ambinet test data.S is the original signal,TS is the detrended&truncated signal,NS is SVD-processed signal,FS is EKF filtered signal
SamFreq=200;Rsl=1e-3;% Ambient test parameters
ImDataFile1='C:\Users\pengbin\OneDrive - usst.edu.cn\桌面\Publication\振动与冲击\data1.mat';
DRepos=importdata(ImDataFile1);
IO=[DRepos.input,DRepos.output1,DRepos.output2,DRepos.output3,DRepos.output4,DRepos.output5,DRepos.output6,DRepos.output7];

%detrend
DtrdOdr=2;
t=[0:1/SamFreq:(size(IO,1)-1)/SamFreq]';
for i=1:size(IO,2)
    IO(:,i)=dtrend(IO(:,i)-polyval(polyfit(t,IO(:,i),DtrdOdr),t));
end

%delete small voltage turbulents
for i=1:size(IO,2)
    for j=1:size(IO(:,i))
      if abs(IO(j,i))>Rsl
        TIO(j,i)=IO(j,i);
        j=j+1;
       else
         TIO(j,i)=0;
         j=j+1;
        end
    end
end

% SVD decompositon the output records
for j=2:8;
  NS(:,j-1)=TIO(:,j);
end
[U,e,V] = svd(NS,0);
NS=U(:,1)*e(1,1)*V(1,:);

%% Preparation of structural matricies
M=DRepos.M; 
C=DRepos.C;
K=DRepos.K; 
D=DRepos.juzhen; 

%% EKF construction
% Omiga=14.55*(2*pi);Amp=mean(abs(TS7));Phase=0;% priors
% StateTranF=@(T) (T+sign(T)*sign((1-(T/Amp)^2))*(Amp*Omiga*sqrt(sign((1-(T/Amp)^2))*(1-(T/Amp)^2)))*(1/SamFreq));
DeltaT=1/SamFreq;

Phi=[zeros(size(M)),eye(size(M));(-1)*(M^-1)*K,(-1)*(M^-1)*C];
Psi=[zeros(size(M));(-1)*M];

A=expm(Phi*DeltaT);
B=(Phi^-1)*(1-(expm(Phi*DeltaT))^-1)*Psi;

StateFcn=@(X,U)(A*X+B*U);
MeasurementFcn=@(Y)(D*Y);
obj = extendedKalmanFilter(StateFcn,MeasurementFcn,NS(1,1)*ones(size(A,2),1),'StateCovariance',5);
obj.ProcessNoise=0.618;
obj.MeasurementNoise=1; 

for k = 1:size(NS)
  [CorrectedState,CorrectedStateCovariance] = correct(obj,NS(k,1:7)'); 
  [PredictedState,PredictedStateCovariance] = predict(obj,TIO(k,1)*ones(size(B,2),1));
  FS(1,1:7)=(D*CorrectedState)';
end
% FS=polyval(polyfit(t,FS,DtrdOdr),t);

%% PSD of the filtered signal
[Psd0,f0]=pwelch(TS7,[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
[Psd1,f1]=pwelch(NS,[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
[Psd2,f2]=pwelch(FS,[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
Psd0(1:20)=0;
Psd1(1:20)=0;
Psd2(1:20)=0;

%% Output
figure(1)
plot(t,TS7); hold on;
plot(t,NS); hold on;
plot(t,FS); 
legend('TS','NS','FS');
xlabel('\fontname{宋体}时间\fontname{Times new Roman}(s)','FontSize',10);
ylabel('\fontname{宋体}加速度\fontname{Times new Roman}(cm/s^{2})','FontSize',10);

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
