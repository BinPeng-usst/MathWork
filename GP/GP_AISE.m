%%%%%%%%%%%data prepare
clear all;
% DataDir='C:\Users\Administrator\OneDrive\Advances in Structural Engineering\Initial submission\Draft.xls';                                      
DataFile='C:\Users\Bin Peng\OneDrive\Advances in Structural Engineering\Initial submission\Draft.xls';   
X1=xlsread(DataFile,'collected results (1)','A3:B45');
X2=xlsread(DataFile,'collected results (2)','J2:K36');
R1=xlsread(DataFile,'collected results (1)','C3:C45');
R2=xlsread(DataFile,'collected results (2)','M2:M36');
TestI1=xlsread(DataFile,'Test results','B2:C7');
TestO1=xlsread(DataFile,'Test results','G2:G7');

TestI=TestI1([1 2 3 4 5],[1 2]);
TestO=TestO1([1 2 3 4 5]);

X=[X1;X2;TestI];
R=[R1;R2;TestO];

% X50=X(1:50,1:2);
% R50=R(1:50);

% Frange=[0.0001:0.012:1.2001];
% Crange=[0.0001:0.02:2.0001];
% for j=0:size(Frange,2)-1
%     for k=1:100
%     xs(j*100+k,1)=Frange(1,k); 
%     end
% end
% xs(:,2)=repmat(Crange',100,1); 

xs=TestI1([6],[1 2]); 

%%%prior models
%  n = input('Enter a number: ');
% n=4;
% switch n
%     case 1
      hfcn1=@(Me) Me(:,1)+0.19*Me(:,2);
      gprMdl1 = fitrgp(X,R,'BasisFunction',hfcn1,'Beta',[1],'KernelFunction','ardsquaredexponential',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres1=resubPredict(gprMdl1);
      [ypred1,ysd1,yint1] = predict(gprMdl1,xs,'Alpha',0.1);
      L1 = resubLoss(gprMdl1)/var(R);
      Sig1=sqrt(gprMdl1.Sigma^2+(gprMdl1.KernelInformation.KernelParameters(end))^2);
          
%       gprMdl1_p50 = fitrgp(X50,R50,'BasisFunction',hfcn1,'Beta',[1],'KernelFunction','ardsquaredexponential');
%       yres1_p50=resubPredict(gprMdl1_p50);
%       [ypred1_p50,ysd1_p50,yint1_p50] = predict(gprMdl1_p50,xs,'Alpha',0.1);
%       L1_p50 = resubLoss(gprMdl1_p50)/var(R50);
%       
%       L1pv =sum(((X(:,1)+0.19*X(:,2))-R).^2)/(size(X,1)*var(R));
%       
%       [yp1,ysd1,yint1] = predict(gprMdl1,X,'Alpha',0.1);
%       Rs1=1-(sum((R-yp1).^2))/var(R);
%       
%       [yp1_p50,ysd1_p50,yint1_p50] = predict(gprMdl1,X50,'Alpha',0.1);
%       Rs1_p50=1-(sum((R50-yp1_p50).^2))/var(R50);
%       
%       Rs1pv=1-sum(((X(:,1)+0.19*X(:,2))-R).^2)/var(R);
      
     %     case 2
       hfcn2=@(Me) 1.0*Me(:,1).*sqrt(1+1.0*Me(:,2)./Me(:,1));
       gprMdl2 = fitrgp(X,R,'BasisFunction',hfcn2,'Beta',[1],'KernelFunction','ardsquaredexponential',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
       yres2=resubPredict(gprMdl2);
       [ypred2,ysd2,yint2] = predict(gprMdl2,xs,'Alpha',0.1);
       L2 = resubLoss(gprMdl2)/var(R);
       Sig2=sqrt(gprMdl2.Sigma^2+(gprMdl2.KernelInformation.KernelParameters(end))^2);
%        gprMdl2_p50 = fitrgp(X50,R50,'BasisFunction',hfcn2,'Beta',[1],'KernelFunction','ardsquaredexponential');
%        yres2_p50=resubPredict(gprMdl2_p50);
%        [ypred2_p50,ysd2_p50,yint2_p50] = predict(gprMdl2_p50,xs,'Alpha',0.1);
%        L2_p50 = resubLoss(gprMdl2_p50)/var(R50);
%        
%        L2pv =sum((R-1.0*X(:,1).*sqrt(1+1.0*X(:,2)./X(:,1))).^2)/(size(X,1)*var(R));
%        
%       [yp2,ysd2,yint2] = predict(gprMdl2,X,'Alpha',0.1);
%       Rs2=1-(sum((R-yp2).^2))/var(R);
%       
%       [yp2_p50,ysd2_p50,yint2_p50] = predict(gprMdl2,X50,'Alpha',0.1);
%       Rs2_p50=1-(sum((R50-yp2_p50).^2))/var(R50);
%       
%       Rs2pv=1-sum((R-1.0*X(:,1).*sqrt(1+1.0*X(:,2)./X(:,1))).^2)/var(R);
      
%     case 3
      hfcn3=@(Me) 0.25+0.45*Me(:,2);
      gprMdl3 = fitrgp(X,R,'BasisFunction',hfcn3,'Beta',[1],'KernelFunction','ardsquaredexponential',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres3=resubPredict(gprMdl3);
      [ypred3,ysd3,yint3] = predict(gprMdl3,xs,'Alpha',0.1);
      L3 = resubLoss(gprMdl3)/var(R);
      Sig3=sqrt(gprMdl3.Sigma^2+(gprMdl3.KernelInformation.KernelParameters(end))^2);
%       gprMdl3_p50 = fitrgp(X50,R50,'BasisFunction',hfcn3,'Beta',[1],'KernelFunction','ardsquaredexponential');
%       yres3_p50=resubPredict(gprMdl3_p50);
%       [ypred3_p50,ysd3_p50,yint3_p50] = predict(gprMdl3_p50,xs,'Alpha',0.1);
%       L3_p50 = resubLoss(gprMdl3_p50)/var(R50);   
%       
%       L3pv =sum((R-(0.25+0.45*X(:,2))).^2)/(size(X,1)*var(R));
%       
%       [yp3,ysd3,yint3] = predict(gprMdl3,X,'Alpha',0.1);
%       Rs3=1-(sum((R-yp3).^2))/var(R);
%       
%       [yp3_p50,ysd3_p50,yint3_p50] = predict(gprMdl3,X50,'Alpha',0.1);
%       Rs3_p50=1-(sum((R50-yp3_p50).^2))/var(R50);
%       
%       Rs3pv=1-sum((R-(0.25+0.45*X(:,2))).^2)/var(R);
      
%     otherwise
      gprMdl4 = fitrgp(X,R,'BasisFunction','none','KernelFunction','ardsquaredexponential',...
    'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'Verbose',0);
      yres4=resubPredict(gprMdl4);
      [ypred4,ysd4,yint4] = predict(gprMdl4,xs,'Alpha',0.1);
      L4 = resubLoss(gprMdl4)/var(R);
      Sig4=sqrt(gprMdl4.Sigma^2+(gprMdl4.KernelInformation.KernelParameters(end))^2);
%       gprMdl4_p50 = fitrgp(X50,R50,'BasisFunction','none','KernelFunction','ardsquaredexponential');
%       yres4_p50=resubPredict(gprMdl4_p50);
%       [ypred4_p50,ysd4_p50,yint4_p50] = predict(gprMdl4_p50,xs,'Alpha',0.1);
%       L4_p50 = resubLoss(gprMdl4_p50)/var(R50);
%       
%       [yp4,ysd4,yint4] = predict(gprMdl4,X,'Alpha',0.1);
%       Rs4=1-(sum((R-yp4).^2))/var(R);
%       
%       [yp4_p50,ysd4_p50,yint4_p50] = predict(gprMdl4,X50,'Alpha',0.1);
%       Rs4_p50=1-(sum((R50-yp4_p50).^2))/var(R50);
     
%    end
         


%%%%%%%%%%%%%plot 2-D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1);
% subplot(1,2,1);
% % f = [ypred+2*ysd; flipdim(ypred-2*ysd,1)];
% f=[yint(:,1); flipdim(yint(:,2),1)];
% fill([xs(:,2); flipdim(xs(:,2),1)], f,'g','FaceAlpha',0.2,'EdgeColor','g','EdgeAlpha',0.2,'linestyle', '- .');
% hold on; 
% % f = [ymu+2*sqrt(ys2); flipdim(ymu-2*sqrt(ys2),1)];
% % patch([xs(:,2); flipdim(xs(:,2),1)], f, 'b','FaceAlpha',0.2);
% % hold on; 
% f = [ypred+1*ysd; flipdim(ypred-1*ysd,1)];
% fill([xs(:,2); flipdim(xs(:,2),1)], f, 'r','FaceAlpha',0.5,'EdgeColor','r','EdgeAlpha',0.5);
% hold on; 
% plot(xs(:,2), ypred,'k--','linewidth',2.5); 
% scatter(X(:,2), R, 56, 'ko');
% scatter(TestI(:,2), TestO, 56,'r*');
% % plot(TestI(:,2), TestSO, 'r+','MarkerSize',8);
% set(gca,'XDir','reverse','xlim',[0 2.0],'xtick',[0:0.2:2.0],'ylim',[0,1.2],'ytick',[0:0.2:1.2],'YColor','none','DataAspectRatio',[1 1 1],'position',[0.1 0.1 0.35 0.5]);
% xlabel('Nominal compressive stress {\it\bfs}_{v0}(N/mm^{2})');
% % ylabel('Nominal shear strength {\it\bf\tau} (N/mm^{2})','Rotation',0,'Position',[0 1.02 0],'Color','k');
% box off;
% 
% % figure(2);
% subplot(1,2,2);
% % f = [ypred+2*ysd; flipdim(ypred-2*ysd,1)];
% f=[yint(:,1); flipdim(yint(:,2),1)];
% fill([xs(:,1); flipdim(xs(:,1),1)], f,'g','FaceAlpha',0.2,'EdgeColor','g','EdgeAlpha',0.2,'linestyle', '- .');
% hold on; 
% % f = [ymu+2*sqrt(ys2); flipdim(ymu-2*sqrt(ys2),1)];
% % fill([xs(:,1); flipdim(xs(:,1),1)], f, 'b','FaceAlpha',0.2);
% % hold on; 
% f = [ypred+1*ysd; flipdim(ypred-1*ysd,1)];
% fill([xs(:,1); flipdim(xs(:,1),1)], f, 'r','FaceAlpha',0.5,'EdgeColor','r','EdgeAlpha',0.5);
% % errorbar(xs(:,1),ymu,0.5*ys2);
% hold on;
% plot(xs(:,1), ypred,'k--','linewidth',2.5); 
% scatter(X(:,1), R, 56, 'ko');
% scatter(TestI(:,1), TestO, 56,'r*');
% % plot(TestI(:,1), TestSO, 'r+','MarkerSize',8);
% set(gca,'xlim',[0 1.2],'xtick',[0:0.2:1.2],'ylim',[0,1.2],'ytick',[0:0.2:1.2],'DataAspectRatio',[1 1 1],'position',[0.365 0.1 0.35 0.5]);
% xlabel('Cohesion {\it\bff}_{v0}(N/mm^{2})');
% % legend('Interval of ¡À3{\sigma}','Interval of ¡À1{\sigma}','Mean','prior mean','Collected result for training','Experimental result','Location','westoutside');
% % legend('boxoff');
% box off;
% set(gcf,'units','normalized','position',[0.3 0.1 0.7 0.5]);


% %%%%%%%%%%%%%plot 3-D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% TubeLike(xs(:,1),xs(:,2),ypred,ysd.*ysd);//¹Ü×´Í¼
% set(gca,'xlim',[0 0.6],'xtick',[0:0.2:0.6],'ylim',[0 1.4],'ytick',[0:0.2:1.4],'zlim',[0 1.2],'ztick',[0:0.2:1.2],'DataAspectRatio',[1 1 1]);
% xlabel('Cohesion {\it\bff}_{v0}(N/mm^{2})','Rotation',45);
% ylabel('Nominal compressive stress {\it\bfs}_{v0}(N/mm^{2})','Rotation',-45);
% zlabel('Nominal shear strength {\it\bf\tau} (N/mm^{2})');
% shading interp;
% light;
% lighting gouraud;
% material metal ;
% axis vis3d;
% subplot(1,3,3);
% plot3(xs(:,1), xs(:,2),ymu,'-.','linewidth',2);
% hold on; 
% plot3(X(:,1), X(:,2), R, 'ko');
% hold on;
% plot3(TestI(:,1), TestI(:,2),TestO,'k*');
% hold on;
% plot3(TestI(:,1), TestI(:,2),TestSO,'R+');
% grid on;
% box off;
% figure(1);
% subplot(2,2,1);%ÇúÃæÍ¼
% [Xcora Ycora]=meshgrid(Frange,Crange);
% ypv1a=Xcora+0.19*Ycora;
% surf(Xcora,Ycora,ypv1a,'FaceAlpha',0.7,'EdgeAlpha',0);
% hold on;
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% ypv1=Xcor+0.19*Ycor;
% % L1pv =sum(((X(:,1).*sqrt(1+X(:,2)./X(:,1))-R).^2))/(size(X,1)*var(R));
% m=0;n=0;
% for i=1:size(R)
%     if ypv1(i,i)>R(i)
%     n=n+1;
%     Rp1l(n)=R(i);
%     Xp1l(n)=Xcor(1,i);
%     Yp1l(n)=Ycor(i,1);        
%     elseif ypv1(i,i)<=R(i)
%     m=m+1;
%     Rp1u(m)=R(i);
%     Xp1u(m)=Xcor(1,i);
%     Yp1u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(Xp1l,Yp1l,Rp1l,60,'b','o');
% hold on;
% scatter3(Xp1u,Yp1u,Rp1u,60,'r','*');
% hold on;
% view(51,45);
% % legend('','','','','location','southoutside');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% subplot(2,2,2);%ÇúÃæÍ¼
% [Xcora Ycora]=meshgrid(Frange,Crange);
% ypv2a=Xcora.*sqrt(1+Ycora./Xcora);
% surf(Xcora,Ycora,ypv2a,'FaceAlpha',0.7,'EdgeAlpha',0.0);
% 
% hold on;
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% ypv2=Xcor.*sqrt(1+Ycor./Xcor);
% % L2pv =sum(((X(:,1).*sqrt(1+X(:,2)./X(:,1))-R).^2))/(size(X,1)*var(R));
% m=0;n=0;
% for i=1:size(R);
%     if ypv2(i,i)>R(i)
%     n=n+1;
%     Rp2l(n)=R(i);
%     Xp2l(n)=Xcor(1,i);
%     Yp2l(n)=Ycor(i,1);        
%     elseif ypv2(i,i)<=R(i)
%     m=m+1;
%     Rp2u(m)=R(i);
%     Xp2u(m)=Xcor(1,i);
%     Yp2u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(Xp2l,Yp2l,Rp2l,60,'b','o');
% hold on;
% scatter3(Xp2u,Yp2u,Rp2u,60,'r','*');
% hold on;
% view(51,45);
% % legend('','','','','location','southoutside');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% subplot(2,2,3);%ÇúÃæÍ¼
% [Xcora Ycora]=meshgrid(Frange,Crange);
% ypv3a=0.25+0.45*Ycora;
% surf(Xcora,Ycora,ypv3a,'FaceAlpha',0.7,'EdgeAlpha',0.0);
% 
% hold on;
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% ypv3=0.25+0.45*Ycor;
% % L3pv =sum(((X(:,1).*sqrt(1+X(:,2)./X(:,1))-R).^2))/(size(X,1)*var(R));
% m=0;n=0;
% for i=1:size(R)
%     if ypv3(i,i)>R(i)
%     n=n+1;
%     Rp3l(n)=R(i);
%     Xp3l(n)=Xcor(1,i);
%     Yp3l(n)=Ycor(i,1);        
%     elseif ypv3(i,i)<=R(i)
%     m=m+1;
%     Rp3u(m)=R(i);
%     Xp3u(m)=Xcor(1,i);
%     Yp3u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(Xp3l,Yp3l,Rp3l,60,'b','o');
% hold on;
% scatter3(Xp3u,Yp3u,Rp3u,60,'r','*');
% hold on;
% view(51,45);
% % legend('','','','','location','southoutside');
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% set(gcf,'units','normalized','position',[0.1 0.1 0.9 0.9]);
% 
% figure(2);
% subplot(2,2,1);%ÇúÃæÍ¼
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% Zcor=griddata(xs(:,1),xs(:,2),ypred1,Xcor,Ycor,'v4');
% m=0;n=0;
% for i=1:size(R)
%     if Zcor(i,i)>R(i)
%     n=n+1;
%     R1l(n)=R(i);
%     X1l(n)=Xcor(1,i);
%     Y1l(n)=Ycor(i,1);        
%     elseif Zcor(i,i)<=R(i)
%     m=m+1;
%     R1u(m)=R(i);
%     X1u(m)=Xcor(1,i);
%     Y1u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(X1l,Y1l,R1l,60,'b','o');
% hold on;
% scatter3(X1u,Y1u,R1u,60,'r','*');
% hold on;
% [XcorA YcorA]=meshgrid(Frange,Crange);
% ZcorA=griddata(xs(:,1),xs(:,2),ypred1,XcorA,YcorA,'v4');
% % tri=delaunay(Xcor,Ycor);
% % trisurf(tri,Xcor,Ycor,Zcor,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% surf(XcorA,YcorA,ZcorA,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% hold on;
% view(56,10);
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% subplot(2,2,2);%ÇúÃæÍ¼
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% Zcor=griddata(xs(:,1),xs(:,2),ypred2,Xcor,Ycor,'v4');
% m=0;n=0;
% for i=1:size(R);
%     if Zcor(i,i)>R(i)
%     n=n+1;
%     R2l(n)=R(i);
%     X2l(n)=Xcor(1,i);
%     Y2l(n)=Ycor(i,1);        
%     elseif Zcor(i,i)<=R(i)
%     m=m+1;
%     R2u(m)=R(i);
%     X2u(m)=Xcor(1,i);
%     Y2u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(X2l,Y2l,R2l,60,'b','o');
% hold on;
% scatter3(X2u,Y2u,R2u,60,'r','*');
% hold on;
% [XcorA YcorA]=meshgrid(Frange,Crange);
% ZcorA=griddata(xs(:,1),xs(:,2),ypred2,XcorA,YcorA,'v4');
% % tri=delaunay(Xcor,Ycor);
% % trisurf(tri,Xcor,Ycor,Zcor,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% surf(XcorA,YcorA,ZcorA,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% hold on;
% view(56,10);
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% 
% subplot(2,2,3);%ÇúÃæÍ¼
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% Zcor=griddata(xs(:,1),xs(:,2),ypred3,Xcor,Ycor,'v4');
% m=0;n=0;
% for i=1:size(R);
%     if Zcor(i,i)>R(i)
%     n=n+1;
%     R3l(n)=R(i);
%     X3l(n)=Xcor(1,i);
%     Y3l(n)=Ycor(i,1);        
%     elseif Zcor(i,i)<=R(i)
%     m=m+1;
%     R3u(m)=R(i);
%     X3u(m)=Xcor(1,i);
%     Y3u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(X3l,Y3l,R3l,60,'b','o');
% hold on;
% scatter3(X3u,Y3u,R3u,60,'r','*');
% hold on;
% [XcorA YcorA]=meshgrid(Frange,Crange);
% ZcorA=griddata(xs(:,1),xs(:,2),ypred3,XcorA,YcorA,'v4');
% % tri=delaunay(Xcor,Ycor);
% % trisurf(tri,Xcor,Ycor,Zcor,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% surf(XcorA,YcorA,ZcorA,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% hold on;
% view(56,10);
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% 
% subplot(2,2,4);%ÇúÃæÍ¼
% [Xcor Ycor]=meshgrid(X(:,1),X(:,2));
% Zcor=griddata(xs(:,1),xs(:,2),ypred4,Xcor,Ycor,'v4');
% m=0;n=0;
% for i=1:size(R);
%     if Zcor(i,i)>R(i)
%     n=n+1;
%     R4l(n)=R(i);
%     X4l(n)=Xcor(1,i);
%     Y4l(n)=Ycor(i,1);        
%     elseif Zcor(i,i)<=R(i)
%     m=m+1;
%     R4u(m)=R(i);
%     X4u(m)=Xcor(1,i);
%     Y4u(m)=Ycor(i,1); 
%     end 
% end
% scatter3(X4l,Y4l,R4l,60,'b','o');
% hold on;
% scatter3(X4u,Y4u,R4u,60,'r','*');
% hold on;
% [XcorA YcorA]=meshgrid(Frange,Crange);
% ZcorA=griddata(xs(:,1),xs(:,2),ypred4,XcorA,YcorA,'v4');
% % tri=delaunay(Xcor,Ycor);
% % trisurf(tri,Xcor,Ycor,Zcor,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% surf(XcorA,YcorA,ZcorA,'FaceColor','interp','FaceAlpha',0.7,'EdgeAlpha',0.0);
% hold on;
% view(56,10);
% set(gca,'FontName','Times new Roman','FontSize',16,'xlim',[0 1.2],'xtick',[0:0.3:1.2],'ylim',[0 2.1],'ytick',[0:0.3:2.1],'zlim',[0 1.2],'ztick',[0:0.3:1.2]);
% 
% set(gcf,'units','normalized','position',[0.1 0.1 0.9 0.9]);

%%%%%%%%%%%%%export data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FigFileDir='C:\Users\Bin Peng\OneDrive\MathWork\Figures\';
% print('-f1',[FigFileDir,'f1'],'-painters','-dmeta','-r600');
% print('-f2',[FigFileDir,'f2'],'-painters','-dmeta','-r600');
% xlswrite('testdata.xls',ymu,1,'A');
% xlswrite('testdata.xls',ys2,1,'B');
% xlswrite('testdata.xls',X,1,'C:D');
% xlswrite('testdata.xls',R,1,'E');
% xlswrite('testdata.xls',xs,1,'F:G');
% xlswrite('testdata.xls',ys,1,'H');
% xlswrite('testdata.xls',TestI,1,'I:J');
% xlswrite('testdata.xls',TestO,1,'K');
% xlswrite('testdata.xls',TestSO,1,'L');
% xlswrite('testdata.xls',fmu,1,'M');
% xlswrite('testdata.xls',fs2,1,'N');
% xlswrite('testdata.xls',lp,1,'O');





