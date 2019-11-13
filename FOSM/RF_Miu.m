function g_out= RF_Miu(xi,PdofS) %非正态分布的当量正态化
switch PdofS.DistributionName
    case 'GeneralizedExtremeValue'
      if xi<=0
          g_out=0;
          return;
      end
      SigmaN=normpdf(norminv(cdf(PdofS,xi),0,1))/pdf(PdofS,xi);
      g_out=xi-SigmaN*norminv(cdf(PdofS,xi),0,1);
    case 'ExtremeValue'
       SigmaN=normpdf(norminv(1-cdf(PdofS,-1*xi),0,1))/pdf(PdofS,-1*xi);
       g_out=xi-SigmaN*norminv(1-cdf(PdofS,-1*xi),0,1);
    otherwise
       SigmaN=normpdf(norminv(cdf(PdofS,xi),0,1))/pdf(PdofS,xi);
       g_out=xi-SigmaN*norminv(cdf(PdofS,xi),0,1);
end
end