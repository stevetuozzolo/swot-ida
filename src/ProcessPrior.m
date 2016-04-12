function [Prior,jmp,AllObs]=ProcessPrior(Prior,AllObs,jmp,DAll,Obs,D)

N=1E4;

%% 1 handle input prior information
% note that A0min is refined for inclusion in the "jmp" variable at the
% bottom
allA0min=nan(DAll.nR,1);
for i=1:DAll.nR,
    if min(AllObs.dA(i,:))>=0
        allA0min(i,1)=0;
    else
        allA0min(i,1)=-min(AllObs.dA(i,:));
    end
end

meann=Prior.meann;
covn=Prior.stdn/meann;

%% 2 initial probability calculations
%n calcs
v=(covn*meann)^2;
[mun,sigman] = logninvstat(meann,v);

%Q calcs
v=(Prior.covQbar*Prior.meanQbar)^2;
[muQbar,sigmaQbar] = logninvstat(Prior.meanQbar,v);

A0u=allA0min; %arbitrary -- factor should not matter
nu=meann; %ok as long as all r

z1=randn(DAll.nR,N);
z2=randn(DAll.nR,N);
u1=rand(DAll.nR,N);
u2=rand(DAll.nR,N);
na1(DAll.nR)=0;
na2(DAll.nR)=0;

thetaA0=nan(DAll.nR,N);
thetaA0(:,1)=A0u;
thetan=nan(DAll.nR,N);
thetan(:,1)=nu;
thetaQ=nan(DAll.nR,N);
f=nan(DAll.nR,N);

%% 3 chain calculations
tic
disp('Processing prior for each reach...')
for j=1:DAll.nR,
            
    A0u=3*thetaA0(j,1); %the initial value is the minimum 
    nu=thetan(j,1);
    
    %preliminary
    pjmp.stdA0=A0u; % these are just trial-and-error
    pjmp.stdn=0.25*nu;
    
    pjmp.target=0.5; %since each reach is hanlded individually, goal is 50%
        
    pu1=1;
    pu2=lognpdf(nu,mun,sigman);       
    
    Qu=mean( 1./nu.*(A0u+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
    
    fu=lognpdf(Qu,muQbar,sigmaQbar);

    for i=1:N,    
        
        % let the jumps adapt to target efficiencies
        if i==N*0.2,
            pjmp.stdA0=pjmp.stdA0/pjmp.target*(na1(j)/N/0.2);
            pjmp.stdn=pjmp.stdn/pjmp.target*(na2(j)/N/0.2);
        end                

        %A0
        A0v=A0u+z1(j,i).*pjmp.stdA0;

        if A0v<allA0min(j),
            pv1=0; fv=0;
        else
            pv1=1;
            Qv = mean( 1./nu.*(A0v+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
            fv=lognpdf(Qv,muQbar,sigmaQbar);
        end

        MetRatio=fv/fu*pv1/pu1;

        if MetRatio>u1(j,i),
            na1(j)=na1(j)+1;
            A0u=A0v; Qu=Qv;
            fu=fv; pu1=pv1;
        end

        %n
        nv=nu+z2(j,i).*pjmp.stdn;   
        if nv<=0,
            pv2=0;
        else
            pv2=lognpdf(nv,mun,sigman);
        end

        Qv = mean( 1./nv.*(A0u+AllObs.dA(j,:)).^(5/3).*AllObs.w(j,:).^(-2/3).*sqrt(AllObs.S(j,:)) );
        fv=lognpdf(Qv,muQbar,sigmaQbar);

        MetRatio=fv/fu*pv2/pu2;

        if MetRatio>u2(j,i),
            na2(j)=na2(j)+1;
            nu=nv; Qu=Qv;
            fu=fv; pu2=pv2;
        end    

        %store
        thetan(j,i)=nu;    
        thetaA0(j,i)=A0u;    
        thetaQ(j,i)=Qu;
        f(j,i)=fu;
    end
end

disp(['... Done. Elapsed time=' num2str(toc) 'sec.'])

if min(na1)/N < 0.2 || min(na2)/N < 0.2 || max(na1)/N > 0.8 || max(na2)/N > 0.8
    disp('Calculation of prior A0 & n failed; adjust jump standard deviations')
    disp('Execution killed...')
    disp(['Acceptance for A0: ' num2str(na1/N*100) '%'])
    disp(['Acceptance for n: ' num2str(na2/N*100) '%'])

    clear Prior
end

iUse=N/5*4+1:N;

%% 4 posterior Q estimation
Prior.meanA0=mean(thetaA0(:,iUse),2);
Prior.stdA0=std(thetaA0(:,iUse),[],2);
Prior.meann=mean(thetan(:,iUse),2);  %should check these parameters actually fit the posterior...
Prior.stdn=std(thetan(:,iUse),[],2);

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