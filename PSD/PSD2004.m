clear all;clc;close all;
% %% Data input
% % Sampling Parameters
% % ImportSourceFile1='C:\Users\Administrator\Desktop\zhong\1－通道 3.dat';
% % ImportSourceFile2='C:\Users\Administrator\Desktop\zhong\1－通道 3.dat';
% % ImportSourceFile3='C:\Users\Administrator\Desktop\zhong\1－通道 3.dat';
% ImportResponseFile1='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% ImportResponseFile2='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% ImportResponseFile3='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% 
% %OutPutSourceAccTh='C:\Users\Administrator\Desktop\zhong\第3通道InputAcceleration';
% OutPutResponseAccTh='C:\Users\2013-FWQ\Desktop\第1通道Acceleration';
% OutPutResponsePSD='C:\Users\2013-FWQ\Desktop\第1通道PSD';
% %OutPutResponseVelTh='C:\Users\Bin Peng\Desktop\201901李余阳动测\191221\第3通道Velocity';
% %OutPutResponseDisTh='C:\Users\Bin Peng\Desktop\201901李余阳动测\191221\第3通道Displacement';
% % OutPutXCorr='C:\Users\Bin Peng\OneDrive\Temp\第4通道XCorrelation'
% %OutPutCoheret='C:\Users\Bin Peng\OneDrive\Temp\第4通道Coherent';
% 
% SigLen=12803893-12;
% NumSeg=1;
% SamFreq=50;
% %施卫星老师设备
% % SenseOfAcc=(1/49)*980 ;%Gal/v
% % SenseOfInstru=23.8; %Gal/v
% % GainOfInstru=1;
% %东华设备
% SenseOfAcc=32.8 ;%Gal/v
% SenseOfInstru=34.2; %Gal/v
% GainOfInstru=1;
% R=1/((SenseOfInstru/SenseOfAcc)*GainOfInstru);
% th=[1/SamFreq:1/SamFreq:SigLen*(1/SamFreq)]';
% DtrOdr=1;
% 
% %Input data
% % importdata(ImportSourceFile1,' ',1);
% % A1=ans.data*R;
% % importdata(ImportSourceFile2,' ',1);
% % A2=ans.data*R;
% % importdata(ImportSourceFile3,' ',1);
% % A3=ans.data*R;
% % Asum=[A1;A2;A3];
% % for i=1:NumSeg;
% %    for j=1:SigLen;
% %     A(j,i)=Asum((i-1)*SigLen+j);
% %    end
% % end
% % for k=1:NumSeg;
% %      A(:,k)=detrend(A(:,k));去线性趋势项
% %      A(:,k)=A(:,k)-polyval(polyfit(th,A(:,k),DtrOdr),th);%去一般趋势项
% % end
% % In_th=mean(A,2);
% 
% importdata(ImportResponseFile1,'\t',12);
% B1=ans.data*R;
% importdata(ImportResponseFile2,'\t',12);
% B2=ans.data*R;
% importdata(ImportResponseFile3,'\t',12);
% B3=ans.data*R;
% Bsum=[B1;B2;B3];
% for i=1:NumSeg
%    for j=1:SigLen
%      B(j,i)=Bsum((i-1)*SigLen+j);
%    end
% end
% for k=1:NumSeg
%      %B(:,k)=detrend(B(:,k));%仅去线性趋势项
%      B(:,k)=B(:,k)-polyval(polyfit(th,B(:,k),DtrOdr),th);%去一般趋势项
% end
% Acc_th=mean(B,2);
% load('C:\Users\Bin Peng\Desktop\MAT20200420-1254（非安静）\MAT20200420-1254原始数据\202004201254.mat');
load('C:\Users\Bin Peng\Desktop\1-4号场地选取数据\1号场地晚上（8个测点）.mat');
SamFreq=500;
mkr=['o','o','o','+','+','+','*','*','*','x','x','x','s','s','s','d','d','d','^','^','^','v','v','v'];
% for i=1:3:22
for i=3:3:24
% OutPutResponsePSD=['C:\Users\Bin Peng\Desktop\',num2str(i)];
OutPutResponsePSD=['C:\Users\Bin Peng\Desktop\','方向3'];
%% show the In_ch
% figure(1);
% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% plot(th,In_th);
% title('Input accelartion');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Acceleration (cm/s^2)');
% ylim([1.1*min(In_th),1.1*max(In_th)]);
% print('-f1',OutPutSourceAccTh,'-painters','-dmeta','-r600');
%% show the Acc_ch
% figure(2);
% plot(th,Acc_th)
% title('Accelartion');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Acceleration (cm/s^2)');
% ylim([1.1*min(Acc_th),1.1*max(Acc_th)]);
% print('-f2',OutPutResponseAccTh,'-painters','-dmeta','-r600');
%% PSD maker of Acc_th
% figure(3);
% Acc_f=fft(Acc_th,SigLen);
% Psd=Acc_f.*conj(Acc_f)/SigLen;
% f=SamFreq/SigLen*(1:(0.5*SigLen))';
% plot(f,Psd(1:(0.5*SigLen)));
% [Psd,f]=periodogram(Acc_th,[],[],SamFreq);
% plot(f,Psd);
% [Psd,f]=pwelch(Acc_th,[],[],[],SamFreq);%默认窗、分段数和重叠长度
% plot(f,Psd);
% title('Power spectral density');
% xlabel('Frequency (Hz)');
% xlim([0,100]);
% ylim([0,0.05]);
% %(《结构动力学》p.297上的单位为|Acc|^{2}/2Hz,与此处不同）The units of the PSD estimate are in squared magnitude units of the time series data per unit frequency.
% %For example, if the input data is in volts, the PSD estimate is in units of squared volts per unit frequency. 
% %For a time series in volts, if you assume a resistance of 1 Ω and specify the sample rate in hertz, the PSD estimate is in watts per hertz.
% ylabel('|Acc|^{2}/Hz');
% ylim([0,1.1*max(Psd)]);
% [~,S] = max(Psd);
% Smax=S*(SamFreq/SigLen);
% sSmax=num2str(Smax,2);
% text(1,0.9*max(Psd),['\rightarrow',sSmax]);
% print('-f3',OutPutResponsePSD,'-painters','-dmeta','-r600');
%% Displacement time history maker of Acc_th
% figure(5);
% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% th=th';
% D_th=cumtrapz(th,V_th)*(1/SamFreq);
% % D_th(1)=0;
% % for j=2:SigLen
% %     D_th(j)=D_th(j-1)+(V_th(j)+V_th(j-1))*(1/SamFreq)*0.5;
% % end   
% plot(th,D_th);
% title('Displacement');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Displacement (cm)');
% ylim([1.1*min(D_th),1.1*max(D_th)]);
%% Cross correlation  
% figure(6);
% [acor,lag] = xcorr(In_th,Acc_th);
% acor=acor/max(acor);
% [~,tI] = max(abs(acor));
% tlagDiff = lag(tI);
% timeDiff = tlagDiff/SamFreq;
% plot(lag/SamFreq,acor);
% title('Cross correlation function');
% xlim([-0.5*((length(lag)-1)/SamFreq),0.5*((length(lag)-1)/SamFreq)]);
% xlabel('Time lag (s)');
% ylim([1.1*min(acor),1.1*max(acor)]);
% ylabel('|Acc|^2');
%% Cross coherent  
% figure(7);
% fs=[0:1/SamFreq:SamFreq/2];
% [Cxy,F] = mscohere(In_th,Acc_th,[],[],fs,SamFreq);
% [Cxy,F] = mscohere(In_th,Acc_th,[],[],[],SamFreq);%默认窗、分段数和重叠长度
% [~,fI] = max(abs(Cxy));
% flagDiff = F(fI);
% frequencyDiff = flagDiff/SamFreq;
% plot(F,Cxy);
% title('Cross coherent function');
% xlim([0,(length(F)-1)/SamFreq]);
% xlabel('frequency (Hz)');
% ylim([0,1]);
% print('-f7',OutPutCoheret,'-painters','-dmeta','-r600');
%% divide the In_ch to n sections and show the time-history
figure(1);
dD=reshape(Datas(:,i),300000,[]);
bp=0:6000:300000;
for k=1:size(dD,2)
dD(:,k)=dD(:,k)-mean(dD(:,k));
dD(:,k)=detrend(dD(:,k),'linear',bp);
end

SigLen=size(dD,1);
th=[0:1/SamFreq:(SigLen-1)/SamFreq];

for j=1:(size(dD,2)-1)
DD_th(:,j)=cumtrapz(th,dD(:,j))*(1/SamFreq);
end
plot(th,DD_th(:,j));
hold on;
end
title('Displacement');
xlabel('time (s)');
xlim([0,600]);
ylabel('Dis(mm)');
legend('1','2','3','4','5','6','7','8','location','northeastoutside')
print('-f1', OutPutResponsePSD,'-painters','-djpeg','-r600');

%% divide the In_ch to n sections and make the root-mean-square (RMS) level or mean PSD
% figure(1);
% dD=reshape(Datas(:,i),30000,[]);
% bp=0:6000:30000;
% for k=1:size(dD,2)
% dD(:,k)=dD(:,k)-mean(dD(:,k));
% dD(:,k)=detrend(dD(:,k),'linear',bp);
% end
% 
% SigLen=size(dD,1);
% th=0:1:(SigLen-1);
% 
% for j=1:(size(dD,2))
% DD_th(:,j)=cumtrapz(th',dD(:,j))*(1/SamFreq);
% [Psd(:,j),f]=pwelch(DD_th(:,j),[],[],[],SamFreq);%默认窗、分段数和重叠长度
% % DD(j)=rms(DD_th(:,j));
% end
% 
% MPsd=mean(Psd,2);
% % scatter(0:1:119,DD,'marker',mkr(i));
% semilogx(f,MPsd);
% hold on;
% end
% title('PSD ');
% xlabel('Frequence(Hz)');
% xlim([1,200]);
% % ylabel('RMS(mm)');
% legend('1','2','3','4','5','6','7','8','location','northeastoutside')
% print('-f1', OutPutResponsePSD,'-painters','-djpeg','-r600');
%% divide the In_ch to n sections and make the histogram of the RMS of the displacement
% figure(1);
% dD=reshape(Datas(:,i),300000,[]);
% bp=0:6000:300000;
% for k=1:size(dD,2)
% dD(:,k)=dD(:,k)-mean(dD(:,k));
% dD(:,k)=detrend(dD(:,k),'linear',bp);
% end
% 
% SigLen=size(dD,1);
% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% 
% for j=1:(size(dD,2))
% DD_th(:,j)=cumtrapz(th,dD(:,j))*(1/SamFreq);
% R_S(j)=rms(DD_th(:,j));
% end
% [htct,edges] = histcounts(R_S,100,'Normalization','probability');
% Bw=edges(2)-edges(1);
% E=0.5*Bw:0.5*Bw:100*0.5*Bw;
% plot(E,htct,'marker',mkr(i));
% hold on;
% end
% title('histogram of the RMS of the displacement');
% % xlim([0,E(end)]);
% % ylim([0,1.1*max(Psd)]);
% xlabel('Displacement (mm)');
% % ylabel('|Dis|^{2}/Hz');
% %(《结构动力学》p.297上的单位为|Acc|^{2}/2Hz,与此处不同）The units of the PSD estimate are in squared magnitude units of the time series data per unit frequency.
% %For example, if the input data is in volts, the PSD estimate is in units of squared volts per unit frequency. 
% %For a time series in volts, if you assume a resistance of 1 Ω and specify the sample rate in hertz, the PSD estimate is in watts per hertz.
% % [~,S] = max(MP);
% % Smax=S*(SamFreq/SigLen);
% % sSmax=num2str(Smax,2);
% % text(1,0.9*max(MP),['\rightarrow',sSmax]);
% legend('1','2','3','4','5','6','7','8','location','northeastoutside')
% print('-f1', OutPutResponsePSD,'-painters','-djpeg','-r600');
% % end
%% divide the In_ch to n sections and list the histogram of the peak-to-peak value
% figure(1);
% dD=reshape(Datas(:,i),300000,[]);
% bp=0:6000:300000;
% for k=1:size(dD,2)
% dD(:,k)=dD(:,k)-mean(dD(:,k));
% dD(:,k)=detrend(dD(:,k),'linear',bp);
% end
% 
% SigLen=size(dD,1);
% th=[0:1:(SigLen-1)];
% 
% for j=1:size(dD,2)
% DD_th(:,j)=cumtrapz(th,dD(:,j))*(1/SamFreq);
% Peak(j)=abs(max(DD_th(:,j))-min(DD_th(:,j)));
% end
% 
% [htct,edges] = histcounts(Peak,100,'Normalization','probability');
% Bw=edges(2)-edges(1);
% E=0.5*Bw:0.5*Bw:100*0.5*Bw;
% plot(E,htct,'marker',mkr(i));
% hold on;
% % ylabel('level');
% % ylim([0,1]);
% end
% title('Histogram of the peak-to-peak value of the displacement');
% xlabel('Peak-to-peak value of the displacement(mm)');
% xlim([0,E(end)]);
% legend('1','2','3','4','5','6','7','8','location','northeastoutside')
% print('-f1', OutPutResponsePSD,'-painters','-djpeg','-r600');
% % end
