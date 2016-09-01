function MakeFigs (D,Truth,Prior,C,E,Err,AllTruth,DAll,AllObs)

figure(1)
% if C.Estimateq,
%     n=3;
% else
%     n=2;
% end
n=3;
subplot(n,1,1); 
plot(C.thetaA0'); grid on;
title('Baseflow cross-sectional area, m^2')
subplot(n,1,2)
plot(C.thetana'); grid on;
title('Roughness coefficient parameter n0')
subplot(n,1,3)
plot(C.thetax1'); grid on;
title('Roughness coefficient parameter x1')
    
figure(2)
h1=errorbar(D.xkm./1000,E.A0hat,E.stdA0Post,'LineWidth',2,'Color','r'); hold on;
hp=errorbar(D.xkm./1000,Prior.meanA0.*ones(D.nR,1),Prior.stdA0.*ones(D.nR,1),'LineWidth',2,'Color','g');
h2=plot(D.xkm./1000,Truth.A0,'LineWidth',2); hold off;
set(gca,'FontSize',14)
xlabel('Flow distance, km')
ylabel('Cross-sectional area, m^2')
legend([h1; hp; h2;],'Estimate','Prior','True')

figure(3)
meanA0=Prior.meanA0;
covA0=Prior.stdA0./meanA0;
vA0=(covA0.*meanA0).^2;
[muA0,sigmaA0] = logninvstat(meanA0,vA0);

for i=1:D.nR,    
    subplot(1,D.nR,i)
    x=0:max(C.thetaA0(i,:));
    y=lognpdf(x,muA0(i),sigmaA0(i));
    plot(x,y/max(y)*C.N/20,'r--','LineWidth',2); hold on;
    hist(C.thetaA0(i,C.Nburn+1:end),50); 
    set(gca,'FontSize',14)
    plot(Truth.A0(i)*ones(2,1),get(gca,'YLim'),'g-','LineWidth',2);    
    plot(E.A0hat(i)*ones(2,1),get(gca,'YLim'),'k--','LineWidth',2); hold off;    
    title(['Reach ' num2str(i)])
    xlabel('A_0, m^2')
    ylabel('Frequency')
end

figure(4)
if isnan(Truth.n),
    Truth.n=nan(size(Truth.A0));
end
    
for i=1:D.nR,
    subplot(1,D.nR,i)
    x=linspace(0,max(C.thetana(i,:)),100);
    y=normpdf(x,Prior.meanna(i),Prior.stdna(i));
    plot(x,y/max(y)*C.N/20,'r--','LineWidth',2); hold on;
    
    hist(C.thetana(i,C.Nburn+1:end),50); 
    plot(Truth.n(i)*ones(2,1),get(gca,'YLim'),'g-','LineWidth',2); 
    set(gca,'FontSize',14)
    plot(E.nahat(i)*ones(2,1),get(gca,'YLim'),'k--','LineWidth',2); hold off;    
    title(['Reach ' num2str(i)])
    xlabel('n0, [-]')
    ylabel('Frequency')    
end

figure(5)
for i=1:D.nR,
    subplot(1,D.nR,i)
    x=linspace(0,max(C.thetax1(i,:)),100);
    y=normpdf(x,Prior.meanx1(i),Prior.stdx1(i));
    plot(x,y/max(y)*C.N/20,'r--','LineWidth',2); hold on;
    
    hist(C.thetax1(i,C.Nburn+1:end),50);     
    set(gca,'FontSize',14)
    plot(E.x1hat(i)*ones(2,1),get(gca,'YLim'),'k--','LineWidth',2); hold off;    
    title(['Reach ' num2str(i)])
    xlabel('x1, [-]')
    ylabel('Frequency')    
end

% Qbar=squeeze(mean(mean(C.thetaAllQ(:,:,C.Nburn+1:end) )));
Qbar=squeeze(mean(mean(C.thetaAllQ(:,:,:) )));

figure(6)
plot(Qbar); grid on;
set(gca,'FontSize',14)
xlabel('Iteration')
ylabel('Discharge, m^3/s')

hold on;
plot([1 C.N],mean(mean(Truth.Q))*ones(2,1));
hold off;


figure(7)
h=plot(D.t,Truth.Q,D.t,E.QhatPostf','LineWidth',2); 
set(h(1:D.nR),'Color','b');  set(h(D.nR+1:end),'Color','r');
set(gca,'FontSize',14)
xlabel('Time, days')
ylabel('Discharge, m^3/s')
legend(h([1 end]),'True','MetroMan')

figure(8)
plot(1:D.nR,Err.QRelErrPrior,'.-',1:D.nR,Err.QRelErrPost,'.-');
xlabel('Reach'); ylabel('Relative error');
legend('Prior','Posterior');

figure(9)
plot(Qbar,C.LogLike,'o')
set(gca,'FontSize',14)
xlabel('Average discharge, m^3/s')
ylabel('Log of likelihood')


figure(10)
plot(DAll.t,mean(AllTruth.Q,1),DAll.t,mean(E.AllQ,1),DAll.t,mean(E.QhatAllPrior,1),'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Discharge, m^3/s')
legend('True','Estimate','Prior')
if DAll.t(1) > datenum(1900,0,0,0,0,0),
    datetick('x','mmm-dd-hh:MM','keepticks')
end

r=D.nR;
figure(11)
subplot(121)
loglog(E.AllQ(r,:)',AllObs.w(r,:)','+');
xlabel('Estimated Discharge, m^3/s')
ylabel('Width, m')
title(['AHG for Reach #' num2str(r)])

subplot(122)
plot(AllObs.h(r,:)',AllObs.w(r,:)','+');
xlabel('Height, m')
ylabel('Width, m')
title(['Stage-area for Reach #' num2str(r)])

iPos=AllObs.S>1E-5; %should make this a variable to be passed in... and should get the "true" slopes added to input file
nTrue=nan(size(AllTruth.Q));
ATrue=AllTruth.A0'*ones(1,DAll.nt)+AllTruth.dA;
nTrue(iPos)=1./AllTruth.Q(iPos).*ATrue(iPos).^(5/3).*AllTruth.W(iPos).^(-2/3).*(AllObs.S(iPos)).^.5;

figure(15)
r=1:DAll.nR;
subplot(121),
plot(E.AllQ(r,:)',E.nhatAll(r,:)','o')
set(gca,'FontSize',14)
xlabel('Estimated Discharge, m^3/s')
ylabel('Estimated n, [-]');
subplot(122)
plot(AllTruth.Q(r,:)',nTrue(r,:)','o')
set(gca,'FontSize',14)
xlabel('True Discharge, m^3/s')
ylabel('"True" n, [-]');
subplot(122)


return