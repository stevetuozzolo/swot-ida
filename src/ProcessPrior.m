<<<<<<< HEAD
function [Prior,jmp,AllObs]=ProcessPrior(Prior,Obs,D,AllObs,jmp,DAll)

%% 1 determine minimum A0 values: note these are different for all and estimation window
% 1.1 calculate minimum A0 values -- this is for "All"
for i=1:D.nR,
=======
function [Prior,jmp,AllObs]=ProcessPrior(Prior,AllObs,DAll,Obs,D,ShowFigs,BjerklienOpt)


%% 1 handle input prior information
% note that A0min is refined for inclusion in the "jmp" variable at the
% bottom
allA0min=nan(DAll.nR,1);
for i=1:DAll.nR,
>>>>>>> varnQbarPrior
    if min(AllObs.dA(i,:))>=0
        allA0min(i,1)=0;
    else
        allA0min(i,1)=-min(AllObs.dA(i,:));
    end
end

<<<<<<< HEAD
% 1.2 calculate minimum values for A0 for the estimation window
for i=1:D.nR,
    if min(Obs.dA(i,:))>=0
        estA0min(i,1)=0;
    else
        estA0min(i,1)=-min(Obs.dA(i,:));
    end
end

% 1.3 shift that the "all" A0 into the estimate window
i1=find(DAll.t==D.t(1));
AllObs.A0Shift=AllObs.dA(:,i1); %this is area at first estimation time > than time 0

% 1.4 for future reference save the more restrictive
jmp.A0min=max(allA0min+AllObs.A0Shift,estA0min);

%% 2 simulated distribution of A0, from the AllObs perspective
N=1E4;

v=(Prior.covQbar*Prior.meanQbar)^2;
[mu,sigma] = logninvstat(Prior.meanQbar,v);
Qbar=lognrnd(mu,sigma,1,N);

dAbar=mean(AllObs.dA,2);
Sbar=mean(AllObs.S,2);
wbar=mean(AllObs.w,2);

for i=1:N,
    A0bar(:,i)= ( Qbar(i).*Prior.meann.*wbar.^(2/3).*Sbar.^(-.5) ).^(3/5)-dAbar;
end

j=all(A0bar>allA0min*ones(1,N));

%switch everything to estimation window perspective
A0bar=A0bar+AllObs.A0Shift*ones(1,N);

Prior.meanA0=mean(A0bar(:,j),2);
Prior.stdA0=std(A0bar(:,j),[],2);

%% 3 move A0 back to discharge to get prior for baseflow
Qavg=nan(N,D.nt);

for i=1:N,
    Q(:,:,i) = 1./(Prior.meann*ones(1,D.nt)) .* ...
        (A0bar(:,i)*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S) ;
    Qavg(i,:)=squeeze(mean(Q(:,:,i),1));
    Qbase(i)=min(Qavg(i,:));
    QstAvg(i)=mean(Qavg(i,:));
end

Prior.meanQbase=mean(Qbase(j));
Prior.stdQbase=std(Qbase(j));

return

% the beginning of an attempt to do part 2 better... abandoned for now...
% see notes on October 14, 2015
% out!
% r=1;
% i=1;
% A0max=10*max(Obs.dA(1,:));
% Ntry=100;
% A0try=linspace(allA0min(1),A0max,Ntry);
% k=1;
% QbarTry(k)=mean(1./Prior.meann(r).* (A0try(k)*ones(1,D.nt)+Obs.dA(r,:)).^(5/3) .* ...
%     Obs.w(r,:).^(-2/3).* Obs.S(r,:).^(1/2) );
=======
%% 2 Bjerklie calcs

Prior.Wa=mean(AllObs.w,2);
Prior.Ha=mean(AllObs.h,2);

for r=1:DAll.nR,
    DB.Sa(r)=mean(AllObs.S(r,:));

    Prior.Wa(r)=mean(AllObs.w(r,:));
    Prior.Ha(r)=mean(AllObs.h(r,:));
    DB.chiW(r)=std(AllObs.w(r,:))/Prior.Wa(r);
    DB.chiH(r)=std(AllObs.h(r,:))/Prior.Ha(r);

    if BjerklienOpt < 3,
        c1=0.85;
        meanx1(r)=2.257+1.308*log10(DB.chiH(r))+0.99*log10(DB.chiW(r))+0.435*log10(DB.Sa(r));
        meanna(r)=0.22*Sa^0.18; %this is "na" in Bjerklie's notation
    elseif BjerklienOpt == 3,
        c1=nan;
        meanx1(r)=-0.09; %these values computed across 10 rivers, all reaches
        meanna(r)=0.04;
    end
end

covna=.05;
covx1=.25;

%% 3 initial probability calculations
%n calcs
v=(covna.*meanna).^2; %ok, could just do Prior.stdn^2... 
[mun,sigman] = logninvstat(meanna,v);

%x1 calcs: note the pdf is actually for -x1 
% covx1=0.25; %move this to be a prior input
% covx1=1; %move this to be a prior input
v=(covx1.*meanx1).^2;
[mux1,sigmax1] = logninvstat(meanx1,v);

%Q calcs
v=(Prior.covQbar*Prior.meanQbar)^2;
[muQbar,sigmaQbar] = logninvstat(Prior.meanQbar,v);


%% chain setup

N=1E4; %chain length

A0u=allA0min; %arbitrary -- factor should not matter
A0u(A0u==0)=1; 
nau=meanna; %ok as long as all r
x1u=meanx1;

z1=randn(DAll.nR,N); %A0
z2=randn(DAll.nR,N); %n
z3=randn(DAll.nR,N); %x1
u1=rand(DAll.nR,N);
u2=rand(DAll.nR,N);
u3=rand(DAll.nR,N);
na1(DAll.nR)=0;
na2(DAll.nR)=0;
na3(DAll.nR)=0;

thetaA0=nan(DAll.nR,N);
thetaA0(:,1)=A0u;
thetana=nan(DAll.nR,N);
thetana(:,1)=nau;
thetax1=nan(DAll.nR,N);
thetax1(:,1)=x1u;
thetaQ=nan(DAll.nR,N);
f=nan(DAll.nR,N);

%% 3 chain calculations
tic
disp('Processing prior for each reach...')
for j=1:DAll.nR,
    
    if j==11,
        stop=1;
    end
            
    A0u=3*thetaA0(j,1); %the initial value is the minimum 
    nau=thetana(j,1);
    x1u=thetax1(j,1);
    
    %preliminary
    pjmp.stdA0=A0u; % these are just trial-and-error
    pjmp.stdna=nau;
    pjmp.stdx1=x1u;
    
    pjmp.target=0.5; %since each reach is hanlded individually, goal is 50%
        
    pu1=1;
    pu2=lognpdf(nau,mun(j),sigman(j));     
    pu3=lognpdf(-x1u,mux1(j),sigmax1(j));     
    
%     nhatu=c1.*( AllObs.w(j,:).*AllObs.h(j,:)./Prior.Wa(j)./Prior.Ha(j) ).^x1u .* nau; %updated if either nau or x1u change
    nhatu = calcnhat(AllObs.w(j,:),AllObs.h(j,:),Prior.Wa(j),Prior.Ha(j),c1,x1u,nau,BjerklienOpt);
    
    Qu=mean( 1./nhatu.*(A0u+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
    
    fu=lognpdf(Qu,muQbar,sigmaQbar);

    for i=1:N,    
        
        % let the jumps adapt to target efficiencies
        if i<N*0.2 && i~=1 && mod(i,100)==0,
            pjmp.stdA0=mean(pjmp.record.stdA0(j,1:i-1))/pjmp.target*(na1(j)/i);
            pjmp.stdna=mean(pjmp.record.stdna(j,1:i-1))/pjmp.target*(na2(j)/i);
            pjmp.stdx1=mean(pjmp.record.stdx1(j,1:i-1))/pjmp.target*(na3(j)/i);                        
        end                
        
        pjmp.record.stdA0(j,i)=pjmp.stdA0;
        pjmp.record.stdna(j,i)=pjmp.stdna;
        pjmp.record.stdx1(j,i)=pjmp.stdx1;

        %A0
        A0v=A0u+z1(j,i).*pjmp.stdA0;

        if A0v<allA0min(j),
            pv1=0; fv=0;
        else
            pv1=1;
            Qv = mean( 1./nhatu.*(A0v+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
            fv=lognpdf(Qv,muQbar,sigmaQbar);
        end

        MetRatio=fv/fu*pv1/pu1;

        if MetRatio>u1(j,i),
            na1(j)=na1(j)+1;
            A0u=A0v; Qu=Qv;
            fu=fv; pu1=pv1;
        end

        %na
        nav=nau+z2(j,i).*pjmp.stdna;   
        if nav<=0,
            pv2=0;
        else
            pv2=lognpdf(nav,mun(j),sigman(j));
        end
        
%         nhatv=c1.*( AllObs.w(j,:).*AllObs.h(j,:)./Prior.Wa(j)./Prior.Ha(j) ).^x1u .* nav; %updated if either nau or x1u change
        nhatv = calcnhat(AllObs.w(j,:),AllObs.h(j,:),Prior.Wa(j),Prior.Ha(j),c1,x1u,nav,BjerklienOpt);

        Qv = mean( 1./nhatv.*(A0u+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
        fv=lognpdf(Qv,muQbar,sigmaQbar);

        MetRatio=fv/fu*pv2/pu2;

        if MetRatio>u2(j,i),
            na2(j)=na2(j)+1;
            nau=nav; 
            Qu=Qv;
            fu=fv; pu2=pv2;
            
        end    

        %x1
        x1v=x1u+z3(j,i).*pjmp.stdx1;   
        if x1v>=0,
            pv3=0;
        else
            pv3=lognpdf(-x1v,mux1(j),sigmax1(j));
        end
        
%         nhatv=c1.*( AllObs.w(j,:).*AllObs.h(j,:)./Prior.Wa(j)./Prior.Ha(j) ).^x1v .* nau; %updated if either nau or x1u change
        nhatv = calcnhat(AllObs.w(j,:),AllObs.h(j,:),Prior.Wa(j),Prior.Ha(j),c1,x1v,nau,BjerklienOpt);

        Qv = mean( 1./nhatv.*(A0u+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
        fv=lognpdf(Qv,muQbar,sigmaQbar);

        MetRatio=fv/fu*pv3/pu3;

        if MetRatio>u3(j,i),
            na3(j)=na3(j)+1;
            x1u=x1v; 
            Qu=Qv;
            fu=fv; pu3=pv3;
            
        end    
        
        
        %store
        thetaA0(j,i)=A0u;    
        thetana(j,i)=nau;    
        thetax1(j,i)=x1u;    
        thetaQ(j,i)=Qu;
        f(j,i)=fu;
    end
    
    stop=1;
end

disp(['... Done. Elapsed time=' num2str(toc) 'sec.'])

if min(na1)/N < 0.15 || min(na2)/N < 0.15 || min(na3)/N < 0.15 || max(na1)/N > 0.8 || max(na2)/N > 0.8 || max(na3)/N > 0.8
    disp('Calculation of prior A0 & n failed; adjust jump standard deviations')
    disp('Execution killed...')
    disp(['Acceptance for A0: ' num2str(na1/N*100) '%'])
    disp(['Acceptance for na: ' num2str(na2/N*100) '%'])
    disp(['Acceptance for x1: ' num2str(na3/N*100) '%'])

    clear Prior
    return
end

iUse=N/5+1:N;

%% 4 posterior Q estimation
Prior.meanA0=mean(thetaA0(:,iUse),2);
Prior.stdA0=std(thetaA0(:,iUse),[],2);
Prior.meanna=mean(thetana(:,iUse),2);  %should check these parameters actually fit the posterior...
Prior.stdna=std(thetana(:,iUse),[],2);
Prior.meanx1=mean(thetax1(:,iUse),2);  %should check these parameters actually fit the posterior...
Prior.stdx1=std(thetax1(:,iUse),[],2);

for r=1:DAll.nR,
%     nhat=c1.*( AllObs.w(r,:).*AllObs.h(r,:)./Prior.Wa(r)./Prior.Ha(r) ).^Prior.meanx1(r) .* Prior.meanna(r); %updated if either nau or x1u change
    nhat = calcnhat(AllObs.w(j,:),AllObs.h(j,:),Prior.Wa(j),Prior.Ha(j),c1,Prior.meanx1(r),Prior.meanna(r),BjerklienOpt);
    QPrior(r,:)=1./nhat.*(Prior.meanA0(r)+AllObs.dA(r,:)).^(5/3).*AllObs.w(r,:).^(-2/3).*sqrt(AllObs.S(r,:)) ;
end

if ShowFigs,

    %check validity of the n parameterization ...
    r=3; %reach to check out.

    Nuse=length(iUse);
    
    figure(11)
    h=CompareLogN(Prior.meanna(r),Prior.stdna(r)./Prior.meanna(r),thetana(r,iUse));
    
    set(gca,'FontSize',14)
    ylabel(h(1),'Histogram of the posterior')
    ylabel(h(2),'Probability')
    xlabel('Roughness coefficient parameter, na, [-]')
    title(['Reach #' num2str(r) ' for na'])

    %check validity of the A0 parameterization ...
    r=2; %reach to check out.
    
    figure(12)
    h=CompareLogN(Prior.meanA0(r),Prior.stdA0(r)./Prior.meanA0(r),thetaA0(r,iUse));
    set(gca,'FontSize',14)
    ylabel(h(1),'Histogram of the posterior')
    ylabel(h(2),'Probability')
    xlabel('A_0, m^2')
    title(['Reach #' num2str(r) ' for A0'])
    
    figure(13)
    h=CompareLogN(Prior.meanQbar,Prior.covQbar,thetaQ(r,iUse)); 
    ylabel(h(1),'Histogram of the posterior')
    ylabel(h(2),'Probability')
    set(gca,'FontSize',14)
    ylabel('Histogram of the posterior')
    xlabel('Discharge, m^3/s')
    
    figure(14)
    hist(thetax1(r,iUse),35)
    
    h=CompareLogN(-Prior.meanx1(r),Prior.stdx1(r)./Prior.meanx1(r),-thetax1(r,iUse)); 
    ylabel(h(1),'Histogram of the posterior')
    ylabel(h(2),'Probability')
    set(gca,'FontSize',14)
    ylabel('Histogram of the posterior')
    xlabel('x1')

    
end

%% 5 

% 5.1 calculate minimum values for A0 for the estimation window
for i=1:D.nR,
    if min(Obs.dA(i,:))>=0
        estA0min(i,1)=0;
    else
        estA0min(i,1)=-min(Obs.dA(i,:));
    end
end

% 5.3 shift that the "all" A0 into the estimate window
i1=find(DAll.t==D.t(1));
AllObs.A0Shift=AllObs.dA(:,i1); %this is area at first estimation time > than time 0

% 5.4 for future reference save the more restrictive
jmp.A0min=max(allA0min+AllObs.A0Shift,estA0min); %this is in the "estimation" window
jmp.nmin=.001; %temporary!

return
>>>>>>> varnQbarPrior
