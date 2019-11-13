clear all;clc;close all;

N=1e3;
Emeanf=-1.16;
Edf=sqrt(0.02);
% meanf=exp(Emeanf+0.5*Edf*Edf);
% df=exp(2*Emeanf+Edf*Edf)*(exp(Edf*Edf)-1);
meanc=0.5;
dc=0.06;
% A=477600;
% Vk=18340;
% DeltaT=0.26;%实验结果变异系数
TestI=[ 0.26 	0.43 
 0.34 	0.54 
 
 0.34 	0.55 
 0.36 	0.54 
 0.36 	0.52 
 0.25 	0.42 
];

TestO=[0.74 
 0.86 
 0.66 
 0.67 
 0.72 
 0.34 
];

% f2s=rand(1,N);
% f2=logninv(f2s,Emeanf,Edf);
fv=lognrnd(Emeanf,Edf,[N,1]);
% cs=rand(1,N);
% c=norminv(cs,meanc,dc);
c=normrnd(meanc,dc,[N,1]);
% Vs=rand(1,N);
% V=gevinv(Vs,1/2.35,0.385*Vk/2.35,0.385*Vk);
% V=gevrnd(1/2.35,0.385*Vk/2.35,0.385*Vk,[N,1]);
% Rtemp=0.1*A*sqrt(fv+3.6*c.*sqrt(fv));
% R=Rtemp-mean(Rtemp)*DeltaT+mean(Rtemp)*DeltaT*2*rand(N,1);

n=1;
switch n
    
    case 1
      a=1.001;b=0.6949;
      R=a*fv+b*c;
      
    
    case 2
      a=0.6007;b=7.09;
      R=a*fv.*sqrt(1+b*(c./fv));
     
    case 3
      a=-0.1856;b=1.701;
      R=a+b*c;
end

% H=(R<V);
% Pf=sum(H)/N;
% Beta=norminv(1-Pf);

%三维极限状态曲面图
% figure(1);
subplot(1,3,1);
[Xcor Ycor]=meshgrid(fv,c);
tri=delaunay(Xcor,Ycor);
Zcor=griddata(fv,c,R,Xcor,Ycor,'v4');
% trimesh(tri,Xcor,Ycor,Zcor);
trisurf(tri,Xcor,Ycor,Zcor,'FaceAlpha',0.5,'EdgeAlpha',0.1);
hold on;
% ZcorUp=griddata(xs(:,1),xs(:,2),ymu+sqrt(ys2),Xcor,Ycor,'v4');
% trisurf(tri,Xcor,Ycor,ZcorUp,'FaceAlpha',0.5,'EdgeAlpha',0.2);
% hold on;
% ZcorDown=griddata(xs(:,1),xs(:,2),ymu-sqrt(ys2),Xcor,Ycor,'v4');
% trisurf(tri,Xcor,Ycor,ZcorDown,'FaceAlpha',0.5,'EdgeAlpha',0.2);
% hold on;
% set(gca,'FontName','Times new Roman','xlim',[0 1.0],'xtick',[0:0.2:1.0],'ylim',[0 1.4],'ytick',[0:0.2:1.4],'zlim',[0 1.4],'ztick',[0:0.2:1.4],'DataAspectRatio',[1 1 1]);
% xlabel('Cohesion {\itf}_{v0}(N/mm^{2})','Rotation',23,'Position',[23.53,29.87,-21.85]);
% ylabel('Nominal compressive stress {\itc}(N/mm^{2})','Rotation',-35,'Position',[22.74,30.56,-21.74]);
% zlabel('Nominal shear strength {\it\tau} (N/mm^{2})');
% hold on;
% Zcor=RegularizeData3D(xs(:,1),xs(:,2),ymu,xs(:,1),xs(:,2));
% mesh(xs(:,1),xs(:,2),Zcor);
plot3(fv,c, R, 'k+');
hold on;
plot3(TestI(:,1), TestI(:,2),TestO,'r*');
% hold on;
% plot3(xs(:,1), xs(:,2),ymu,'ro');
% hold on;
% view(-37.5,10);
% for i=1:36
%    camorbit(10,0);
%    drawnow;
% end

subplot(1,3,2);
xn0=0:0.02:1.0;
yn0=zeros(1,51);
zn0=griddata(fv,c,R,xn0,yn0,'v4');
% xn1=0:0.02:1.0;
% yn1=0.7*ones(1,51);
% zn1=griddata(fv,c,R,xn1,yn1,'v4');
% xn2=0:0.02:1.0;
% yn2=1.4*ones(1,51);
% zn2=griddata(fv,c,R,xn2,yn2,'v4');
plot(xn0,zn0);
% hold on;
% plot(xn1,zn1);
% hold on;
% plot(xn2,zn2);
% hold on;
% set(gca,'FontName','Times new Roman','xlim',[0 1.0],'xtick',[0:0.2:1.0],'ylim',[0 1.4],'ytick',[0:0.2:1.4]);
% xlabel('Cohesion {\itf}_{v0}(N/mm^{2})');
% ylabel('Nominal shear strength {\it\tau} (N/mm^{2})');

subplot(1,3,3);
xn0=zeros(1,51);
yn0=0:0.02:1.0;
zn0=griddata(fv,c,R,xn0,yn0,'v4');
% xn1=0:0.02:1.0;
% yn1=0.7*ones(1,51);
% zn1=griddata(fv,c,R,xn1,yn1,'v4');
% xn2=0:0.02:1.0;
% yn2=1.4*ones(1,51);
% zn2=griddata(fv,c,R,xn2,yn2,'v4');
plot(yn0,zn0);


%%%%%%%%%%%%export data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  print('-f1','H:\OneDrive\MathWork\Figures\f3(a)','-painters','-dmeta','-r600');