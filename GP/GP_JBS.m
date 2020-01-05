%% data prepare
clear all;

DataFile='C:\Users\Bin Peng\OneDrive - usst.edu.cn\桌面\Publication\建筑结构学报（赵文昊）\Draft';   
CollectI1=xlsread(DataFile,'collected results (1)','B2:C40');
CollectI2=xlsread(DataFile,'collected results (2)','K2:L14');
TestI1=xlsread(DataFile,'Test results','B2:C7');
CollectO1=xlsread(DataFile,'collected results (1)','D2:D40');    
CollectO2=xlsread(DataFile,'collected results (2)','N2:N14');
TestO1=xlsread(DataFile,'Test results','G2:G7');

% TestI=TestI1([1 2 3 4 5],[1 2]);
% TestO=TestO1([1 2 3 4 5]);

%Training data
X=[CollectI1;CollectI2([1:8],[1 2])];
R=[CollectO1;CollectO2([1:8])];
% Sort=[X R];
% Sort=sortrows(Sort);
% X=Sort(:,[1 2]);
% R=Sort(:,[3]);

%Testing data
% xs=[0.01:0.01:1.0]'; 
% xs(:,2)=[0.02:0.02:2.0];
xs=[CollectI2([9:13],[1 2]);TestI1]; 
rs=[CollectO2([9:13]);TestO1];

fa=xlsread(DataFile,'collected results (2)','H10:H14');
fb=xlsread(DataFile,'Test results','L11:L16');
f=[fa*0.45;fb];
cs=[CollectI2([9:13],[1 2]);TestI1];
cr=[CollectO2([9:13]);TestO1];
fv=0.42*cs(:,1);
for i=1:size(f)
    if 1.2*cs(i,2)>0.8*f(i)
        cs(i,2)=(0.8*f(i))/1.2;
    end
end
code12=(fv+0.6*(0.26-0.082*(1.2)*cs(:,2)./f).*cs(:,2)*(1.2))/(0.42);
code135=(fv+0.64*(0.23-0.065*(1.35)*cs(:,2)./f).*cs(:,2)*(1.35))/(0.42); 

%% calculation
% switch n
%     case 1
      hfcn1=@(Me) 0.8723*Me(:,1)+0.4146*Me(:,2);
%       gprMdl1 = fitrgp(X,R,'BasisFunction',hfcn1,'Beta',[1],'KernelFunction','ardsquaredexponential');
        gprMdl1 = fitrgp(X,R,'BasisFunction',hfcn1,'Beta',[1],'KernelFunction','ardsquaredexponential',...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
         struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres1=resubPredict(gprMdl1);
      [ypred1,ysd1,yint1] = predict(gprMdl1,xs,'Alpha',0.01);
%       Weight1=sum((yres1-R).^2)/size(yres1,1);
%       Delta1=ypred1-rs;
%       MSE1=(sum(Delta1.^2)/size(Delta1,1))/var(rs);
%       Max1=max(abs(Delta1));
%       Min1=min(abs(Delta1));

%     case 2
       hfcn2=@(Me) 0.2809*Me(:,1).*sqrt(1+20.93*Me(:,2)./Me(:,1));
       gprMdl2 = fitrgp(X,R,'BasisFunction',hfcn2,'Beta',[1],'KernelFunction','ardsquaredexponential',...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
         struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
%        gprMdl2 = fitrgp(X,R,'BasisFunction',hfcn2,'Beta',[1],'KernelFunction','ardsquaredexponential','Leaveout','on');
%        yres2=kfoldPredict(gprMdl2);

        yres2=resubPredict(gprMdl2);       
        [ypred2,ysd2,yint2] = predict(gprMdl2,xs,'Alpha',0.01);
%       Weight2=sum((yres2-R).^2)/size(yres2,1);
%       Delta2=ypred2-rs;
%       MSE2=(sum(Delta2.^2)/size(Delta2,1))/var(rs);
%       Max2=max(abs(Delta2));
%       Min2=min(abs(Delta2));

%     case 3
      hfcn3=@(Me) 0.1818+0.3875*Me(:,2);
      gprMdl3 = fitrgp(X,R,'BasisFunction',hfcn3,'Beta',[1],'KernelFunction','ardsquaredexponential',...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
         struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres3=resubPredict(gprMdl3);
      [ypred3,ysd3,yint3] = predict(gprMdl3,xs,'Alpha',0.01);
%       Weight3=sum((yres3-R).^2)/size(yres3,1);
%       Delta3=ypred3-rs;
%       MSE3=(sum(Delta3.^2)/size(Delta3,1))/var(rs);
%       Max3=max(abs(Delta3));
%       Min3=min(abs(Delta3));

%     otherwise
      gprMdl4 = fitrgp(X,R,'BasisFunction','none','KernelFunction','ardsquaredexponential',...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
         struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres4=resubPredict(gprMdl4);
      [ypred4,ysd4,yint4] = predict(gprMdl4,xs,'Alpha',0.05);
%       Weight4=sum((yres4-R).^2)/size(yres4,1);
%       Delta4=ypred4-rs;
%       MSE4=(sum(Delta4.^2)/size(Delta4,1))/var(rs);
%       Max4=max(abs(Delta4));
%       Min4=min(abs(Delta4));

%     Total probability       
      S=ysd1./abs(ypred1)+ysd2./abs(ypred2)+ysd3./abs(ypred3)+ysd4./abs(ypred4);
      W1=(ysd1./abs(ypred1))./S;
      W2=(ysd2./abs(ypred2))./S;
      W3=(ysd3./abs(ypred3))./S;
      W4=(ysd4./abs(ypred4))./S;
      T=ypred1.*W1+ypred2.*W2+ypred3.*W3+ypred4.*W4;
%       DeltaT=T-rs;
%       MSET=(sum(DeltaT.^2)/size(DeltaT,1))/var(rs);
%       MaxT=max(abs(DeltaT));
%       MinT=min(abs(DeltaT));
% end
%% plot 2-D
  %% 曲线图
     %% 回归结果(形式1）
% figure(1) 
% subplot(2,2,1);
% plot(yres1,'-^r','linewidth',2.5,'markersize',3);
% hold on;
% scatter([1:1:size(R,1)]',R,'b');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 size(R,1)],'xtick',[0:5:size(R,1)],'ylim',[0,1.2],'ytick',[0:0.2:1.2]);
% xlabel('墙体编号','FontName','{宋体}','FontSize',16);
% ylabel('\fontname{宋体}名义抗剪强度\fontname{Times new Roman}{\it\bf\tau} (N/mm^{2})');
% legend('\fontname{宋体}本文模型推断','\fontname{宋体}文献实验结果');
% legend('boxoff');
% box off;
% 
% subplot(2,2,2);
% plot(yres2,'-^r','linewidth',2.5,'markersize',3);
% hold on;
% scatter([1:1:size(R,1)]',R,'b');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 size(R,1)],'xtick',[0:5:size(R,1)],'ylim',[0,1.2],'ytick',[0:0.2:1.2]);
% xlabel('墙体编号','FontName','{宋体}','FontSize',16);
% ylabel('\fontname{宋体}名义抗剪强度\fontname{Times new Roman}{\it\bf\tau} (N/mm^{2})');
% legend('\fontname{宋体}本文模型推断','\fontname{宋体}文献实验结果');
% legend('boxoff');
% box off;
% 
% subplot(2,2,3);
% plot(yres3,'-^r','linewidth',2.5,'markersize',3);
% hold on;
% scatter([1:1:size(R,1)]',R,'b');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 size(R,1)],'xtick',[0:5:size(R,1)],'ylim',[0,1.2],'ytick',[0:0.2:1.2]);
% xlabel('墙体编号','FontName','{宋体}','FontSize',16);
% ylabel('\fontname{宋体}名义抗剪强度\fontname{Times new Roman}{\it\bf\tau} (N/mm^{2})');
% legend('\fontname{宋体}本文模型推断','\fontname{宋体}文献实验结果');
% legend('boxoff');
% box off;
% 
% subplot(2,2,4);
% plot(yres4,'-^r','linewidth',2.5,'markersize',3);
% hold on;
% scatter([1:1:size(R,1)]',R,'b');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 size(R,1)],'xtick',[0:5:size(R,1)],'ylim',[0,1.2],'ytick',[0:0.2:1.2]);
% xlabel('墙体编号','FontName','{宋体}','FontSize',16);
% ylabel('\fontname{宋体}名义抗剪强度\fontname{Times new Roman}{\it\bf\tau} (N/mm^{2})');
% legend('\fontname{宋体}本文模型推断','\fontname{宋体}文献实验结果');
% legend('boxoff');
% box off;
% box off;
% set(gcf,'units','normalized','position',[0.1 0.1 0.9 0.9]);
     %% 回归结果(形式2）
% figure(1) 
% subplot(1,4,1);
% scatter(R,yres1,10,'d','MarkerEdgeColor','r');
% hold on;
% plot([0:0.2:1.2]',[0:0.2:1.2]','r');
% set(gca,'FontName','Times new Roman','FontSize',11,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0,1.2],'ytick',[0:0.3:1.2],'DataAspectRatio',[1 1 1]);
% xlabel('\fontname{宋体}回归结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% box on;
% % ylabel('\fontname{宋体}实验结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% 
% subplot(1,4,2);
% scatter(R,yres2,10,'d','MarkerEdgeColor','r');
% hold on;
% plot([0:0.1:1.2]',[0:0.1:1.2]','r');
% set(gca,'FontName','Times new Roman','FontSize',11,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0,1.2],'ytick',[0:0.3:1.2],'DataAspectRatio',[1 1 1]);
% xlabel('\fontname{宋体}回归结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% box on;
% % ylabel('\fontname{宋体}实验结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% 
% subplot(1,4,3);
% scatter(R,yres3,10,'d','MarkerEdgeColor','r');
% hold on;
% plot([0:0.1:1.2]',[0:0.1:1.2]','r');
% set(gca,'FontName','Times new Roman','FontSize',11,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0,1.2],'ytick',[0:0.3:1.2],'DataAspectRatio',[1 1 1]);
% xlabel('\fontname{宋体}回归结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% box on;
% % ylabel('\fontname{宋体}实验结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% 
% subplot(1,4,4);
% scatter(R,yres4,10,'d','MarkerEdgeColor','r');
% hold on;
% plot([0:0.1:1.2]',[0:0.1:1.2]','r');
% set(gca,'FontName','Times new Roman','FontSize',11,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0,1.2],'ytick',[0:0.3:1.2],'DataAspectRatio',[1 1 1]);
% xlabel('\fontname{宋体}回归结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% box on;
% % ylabel('\fontname{宋体}实验结果\fontname{Times new Roman}(N/mm^{2})','FontSize',10);
% set(gcf,'units','centimeters','position',[1 1 23 12.5/2]);    
     %% 预测结果
figure (1)
subplot(2,2,1)
    yM1= 0.8723*xs(:,1)+0.4146*xs(:,2);
    plot([1:1:size(yM1)]',yM1,'k-.*','markersize',6); 
    hold on;
    plot([1:1:size(ypred1)]',ypred1,'r-^','linewidth',1.5,'markersize',5); 
    hold on;
    scatter([1:size(rs,1)-6]',rs(1:size(rs,1)-6),60,'MarkerEdgeColor','b');
    hold on;
    scatter([size(rs,1)-5:size(rs,1)]',rs(size(rs,1)-5:size(rs,1)),60,'filled','MarkerFaceColor','b');
    hold on;
    f = [ypred1+ysd1; flipdim(ypred1-ysd1,1)];
    fill([[1:1:size(ypred1)]'; flipdim([1:1:size(ypred1)]',1)], f, 'g','FaceAlpha',0.3,'EdgeColor','g','EdgeAlpha',0.5);
    set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 (size(ypred2,1))],'xtick',[0:1:(size(ypred2,1))],'ylim',[0,1.8],'ytick',[0:0.3:1.8]);
    box off;

subplot(2,2,2)
    yM2=0.2809*xs(:,1).*sqrt(1+20.93*xs(:,2)./xs(:,1));
    plot([1:1:size(yM2)]',yM2,'k-.*','markersize',6); 
    hold on;
    plot([1:1:size(ypred2)]',ypred2,'r-^','linewidth',1.5,'markersize',5); 
    hold on;
    scatter([1:size(rs,1)-6]',rs(1:size(rs,1)-6),60,'MarkerEdgeColor','b');
    hold on;
    scatter([size(rs,1)-5:size(rs,1)]',rs(size(rs,1)-5:size(rs,1)),60,'filled','MarkerFaceColor','b');
    hold on;
    f = [ypred2+ysd2; flipdim(ypred2-ysd2,1)];
    fill([[1:1:size(ypred2)]'; flipdim([1:1:size(ypred2)]',1)], f, 'g','FaceAlpha',0.3,'EdgeColor','g','EdgeAlpha',0.5);
    set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 (size(ypred2,1))],'xtick',[0:1:(size(ypred2,1))],'ylim',[0,1.8],'ytick',[0:0.3:1.8]);
    box off;
 
subplot(2,2,3);
    yM3=0.1818+0.3875*xs(:,2);
    plot([1:1:size(yM3)]',yM3,'k-.*','markersize',6); 
    hold on;
    plot([1:1:size(ypred3)]',ypred3,'r-^','linewidth',1.5,'markersize',5); 
    hold on;
    scatter([1:size(rs,1)-6]',rs(1:size(rs,1)-6),60,'MarkerEdgeColor','b');
    hold on;
    scatter([size(rs,1)-5:size(rs,1)]',rs(size(rs,1)-5:size(rs,1)),60,'filled','MarkerFaceColor','b');
    hold on;
    f = [ypred3+ysd3; flipdim(ypred3-ysd3,1)];
    fill([[1:1:size(ypred3)]'; flipdim([1:1:size(ypred3)]',1)], f, 'g','FaceAlpha',0.3,'EdgeColor','g','EdgeAlpha',0.5);
    set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 (size(ypred3,1))],'xtick',[0:1:(size(ypred3,1))],'ylim',[0,1.8],'ytick',[0:0.3:1.8]);
    box off;

subplot(2,2,4);
    plot([1:1:size(ypred4)]',ypred4,'r-^','linewidth',1.5,'markersize',5); 
    hold on;
    scatter([1:size(rs,1)-6]',rs(1:size(rs,1)-6),60,'MarkerEdgeColor','b');
    hold on;
    scatter([size(rs,1)-5:size(rs,1)]',rs(size(rs,1)-5:size(rs,1)),60,'filled','MarkerFaceColor','b');
    hold on;
    f = [ypred4+ysd4; flipdim(ypred4-ysd4,1)];
    fill([[1:1:size(ypred4)]'; flipdim([1:1:size(ypred4)]',1)], f, 'g','FaceAlpha',0.3,'EdgeColor','g','EdgeAlpha',0.5);
    set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 (size(ypred4,1))],'xtick',[0:1:(size(ypred4,1))],'ylim',[0,1.8],'ytick',[0:0.3:1.8]);
    box off;
% xlabel('墙体编号','FontName','{宋体}','FontSize',10);
% ylabel('\fontname{宋体}名义抗剪强度\fontname{Times new Roman}{\it\bf\tau} (N/mm^{2})');
% L1=legend('\fontname{宋体}本文模型推断','\fontname{宋体}文献实验结果','\fontname{宋体}本文实验结果',,'Location','Best');
% set(L1,'FontSize',15,'Box','off');
box off;
set(gcf,'units','normalized','position',[0.1 0.1 0.9 0.9]);
  %% 玫瑰图
% figure(1);
% DeltaT=T-rs;
% MSET=(sum(DeltaT.^2)/size(DeltaT,1))/var(rs);
% MaxT=max(abs(DeltaT));
% MinT=min(abs(DeltaT));
% s1=size(T,1);
% s2=size(CollectO2([9:13]),1);
% s3=size(TestO1,1);
% delta=2*pi/s1;
% polarplot([0:delta:(2*pi-delta)],T,'r-^','linewidth',1.5,'markersize',5);
% hold on;
% polarplot([0:delta:(2*pi-delta)],code12,'k--');
% hold on;
% polarplot([0:delta:(2*pi-delta)],code135,'k--');
% hold on;
% polarscatter([0:delta:(s2-1)*delta]',CollectO2([9:13]),45,'MarkerEdgeColor','b');
% hold on;
% polarscatter([s2*delta:delta:(2*pi-delta)]',TestO1,45,'filled','MarkerFaceColor','b');
% hold on;
% TTick=[0:rad2deg(delta):rad2deg(2*pi-delta)];
% TLable={'\fontname{宋体}墙体1';'2';'3';'4';'\fontname{宋体}墙体5';'6';'7';'8';'\fontname{宋体}墙体9';'10';'11';'12';'\fontname{宋体}墙体13';'14';'15';'16'};
% % L1=legend('','','','Location','Best');
% % set(L1,'FontSize',15,'Box','off');
% set(gca,'FontName','Times new Roman','FontSize',9,'ThetaTick',TTick,'ThetaTickLabel',TLable,'Rlim',[0,1.5],'Rtick',[0:0.3:1.5],'RAxisLocation',45,'RTickLabel',{});
% set(gcf,'units','centimeters','position',[1 1 8.5 8.5]); 
%% plot 3-D
  %% 管状图;
%    figure(1)
%  TubeLike(xs(:,1),xs(:,2),ymu,ys2);
% % % set(gca,'xlim',[0 0.6],'xtick',[0:0.2:0.6],'ylim',[0 1.4],'ytick',[0:0.2:1.4],'zlim',[0 1.2],'ztick',[0:0.2:1.2],'DataAspectRatio',[1 1 1]);
% % % xlabel('Cohesion {\it\bff}_{v0}(N/mm^{2})','Rotation',23,'Position',[12.58,15.88,-11.66]);
% % % ylabel('Nominal compressive stress {\it\bfs}_{v0}(N/mm^{2})','Rotation',-35,'Position',[11.95,16.46,-11.52]);
% % % zlabel('Nominal shear strength {\it\bf\tau} (N/mm^{2})');
% % % shading interp;
% % % light;
% % % lighting gouraud;
% % % material metal ;
% % % axis vis3d;
  %% 曲面图
% figure(1);
% [Xcor Ycor]=meshgrid(xs(:,1),xs(:,2));
% tri=delaunay(Xcor,Ycor);
% Zcor=griddata(xs(:,1),xs(:,2),T,Xcor,Ycor,'v4');
% trisurf(tri,Xcor,Ycor,Zcor,'FaceAlpha',0.5,'EdgeAlpha',0.3);
% hold on;
% 
% % ZCode12=griddata(xs(:,1),xs(:,2),code12,Xcor,Ycor,'v4');
% % trisurf(tri,Xcor,Ycor,ZCode12,'FaceAlpha',0.1,'EdgeAlpha',0.05);
% % hold on;
% % ZCode135=griddata(cs(:,1),cs(:,2),code135,Xcor,Ycor,'v4');
% % trisurf(tri,Xcor,Ycor,ZCode135,'FaceAlpha',0.5,'EdgeAlpha',0.3);
% % hold on;
% % scatter3(cs([1:size(cr,1)-6],[1]),cs([1:size(cr,1)-6],[2]),cr([1:size(cr,1)-6]),60,'MarkerEdgeColor','b');
% % hold on;
% % scatter3(cs([size(cr,1)-5:size(cr,1)],[1]),cs([size(cr,1)-5:size(cr,1)],[2]),cr([size(cr,1)-5:size(cr,1)]),60,'filled','MarkerFaceColor','b');
% % hold on
% colormap(copper);
% set(gca,'FontName','Times new Roman','FontSize',9,'xlim',[0 0.4],'xtick',[0:0.1:0.4],'ylim',[0,1.0],'ytick',[0:0.2:1.0],'zlim',[0,0.6],'ztick',[0:0.2:0.6]);
% set(gcf,'units','centimeters','position',[1 1 8.5 8.5]); 
% view(58,12);
%% export data
RFile='C:\Users\Bin Peng\OneDrive\MathWork\Results\fig5';
print('-f1',RFile,'-painters','-dmeta','-r600');
% xlswrite('testdata.xls',X,1,'C:D');