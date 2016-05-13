function h=CompareLogN(meanx,covx,x)


v=(covx.*meanx).^2;
[mux,sigmax] = logninvstat(meanx,v);

% X = logninv(.99,mux,sigmax) %for future reference, this is the 99th
% percentile of this distribution

[rHist.N,rHist.c]=hist(x,35);

xval=linspace(0.9*min(rHist.c),1.1*max(rHist.c),100);
yval=lognpdf(xval,mux,sigmax);

h=plotyy(rHist.c,rHist.N,xval,yval);

return
