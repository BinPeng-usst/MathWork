%% Brick, normal prior
pdf = @(x)(0.01^(-3))*exp((-1/(2*0.01))*((1.7-x)^2+(1.9-x)^2+(1.7-x)^2+(1.8-x)^2+(1.9-x)^2+(1.5-x)^2))*(((0.01)^(-0.5))*exp((-1/(2*0.01))*(x-2.1)^2));
%proppdf = @(x,y) normpdf(y,x,0.05);
proppdf = @(x,y) y+normpdf(x,0,(2.7-1.1)/2);
proprnd = @(x)  unifrnd(1.1,2.7); %the range is 2 times of the one set by the  minimun and maximum measured values
nsamples = 5500;
[x,accept]= mhsample(1.1,nsamples,'pdf',pdf,'proprnd',proprnd,'proppdf',proppdf,'thin',21,'burnin',2000);
h=histfit(x,100);
X=get(h(2),'Xdata');
Y=get(h(2),'Ydata');

