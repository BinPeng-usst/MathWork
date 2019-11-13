function g_out= RF_Sigma(xi,PdofS) %非正态分布的当量正态化
switch PdofS.DistributionName
    case 'GeneralizedExtremeValue'
      if xi<=0
          g_out=1;
          return;
      end
      g_out=normpdf(norminv(cdf(PdofS,xi),0,1))/pdf(PdofS,xi);
     case 'ExtremeValue'
      g_out=normpdf(norminv(1-cdf(PdofS,-1*xi),0,1))/pdf(PdofS,-1*xi);
     otherwise
         g_out=normpdf(norminv(cdf(PdofS,xi),0,1))/pdf(PdofS,xi);
end
end

