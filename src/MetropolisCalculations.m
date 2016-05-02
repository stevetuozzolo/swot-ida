function [C] = MetropolisCalculations(Prior,D,Obs,jmp,C,R,DAll,AllObs,BjerklienOpt)

[Delta,DeltaA,B,C,thetauA0,thetauna,thetaux1,thetauq,R]=InitializeMetropolis (D,C,Prior,R);

%TEMPORARY ...
jmp.stdx1=0.03.*mean(thetaux1);

% log-normal probabilty calculations
meanA0=Prior.meanA0;
covA0=Prior.stdA0./meanA0;
vA0=(covA0.*meanA0).^2;
[muA0,sigmaA0] = logninvstat(meanA0,vA0);

meanna=Prior.meanna;
covna=Prior.stdna./meanna;
vna=(covna.*Prior.meanna).^2;
[muna,sigmana] = logninvstat(Prior.meanna,vna);

meanx1=Prior.meanx1;
covx1=Prior.stdx1./meanx1;
vx1=(covx1.*Prior.meanx1).^2;
[mux1,sigmax1] = logninvstat(-Prior.meanx1,vx1);

%1) initial probability calculations
pu1=lognpdf(thetauA0,muA0,sigmaA0);
pu2=lognpdf(thetauna,muna,sigmana);
pu3=lognpdf(-thetaux1,mux1,sigmax1);
% if C.Estimateq,
%     pu4=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
% end

[fu,dQdx,dAdt]=CalcLklhd(Obs,thetauA0,thetauna,thetaux1,D,Prior,Delta,DeltaA,B,thetauq,BjerklienOpt);

% if Prior.meanq==-1,
%     Prior.meanq=dQdx+dAdt;
%     C.thetaq(:,1)=Prior.meanq;
%     thetauq=C.thetaq(:,1);
%     pu4=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
% end

%2) The loop
tic

C.n_a1=0; %A0
C.n_a2=0; %na
C.n_a3=0; %x1

for i=1:C.N,
    if mod(i,C.Nburn/2)==0, 
        disp(['Iteration #' num2str(i) '/' num2str(C.N) '.']); 
    end    
        
    %A0
    thetavA0=thetauA0+jmp.stdA0.*R.z1(:,i);   
    thetavA0(thetavA0<jmp.A0min)=jmp.A0min(thetavA0<jmp.A0min);
%     pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));    
    pv1=lognpdf(thetavA0,muA0,sigmaA0);
    fv=CalcLklhd(Obs,thetavA0,thetauna,thetaux1,D,Prior,Delta,DeltaA,B,thetauq,BjerklienOpt);    

    MetRatio=exp(fv-fu)*exp(sum(log(pv1))-sum(log(pu1)));
    if MetRatio>R.u1(i),
        C.n_a1=C.n_a1+1; %increment
        thetauA0=thetavA0; fu=fv; pu1=pv1; %update u->v     
    end    
    C.thetaA0(:,i)=thetauA0;
    
    %na
    thetavna=thetauna+jmp.stdn.*R.z2(:,i);
    thetavna(thetavna<jmp.nmin)=jmp.nmin;
    pv2=lognpdf(thetavna,muna,sigmana);
    fv=CalcLklhd(Obs,thetauA0,thetavna,thetaux1,D,Prior,Delta,DeltaA,B,thetauq,BjerklienOpt);    

    MetRatio=exp(fv-fu)*exp(sum(log(pv2))-sum(log(pu2)));
    if MetRatio>R.u2(i),
        C.n_a2=C.n_a2+1; %increment
        thetauna=thetavna; fu=fv; pu2=pv2; %update u->v     
    end    
    C.thetana(:,i)=thetauna;    

    %x1   
    thetavx1=thetaux1+jmp.stdx1.*R.z3(:,i);
%     thetavna(thetavna<jmp.nmin)=jmp.nmin; %use no limits for now
    pv3=lognpdf(-thetavx1,mux1,sigmax1);
    fv=CalcLklhd(Obs,thetauA0,thetauna,thetavx1,D,Prior,Delta,DeltaA,B,thetauq,BjerklienOpt);    

    MetRatio=exp(fv-fu)*exp(sum(log(pv3))-sum(log(pu3)));
    if MetRatio>R.u3(i),
        C.n_a3=C.n_a3+1; %increment
        thetaux1=thetavx1; fu=fv; pu3=pv3; %update u->v     
    end    
    C.thetax1(:,i)=thetaux1; 
    
%     if C.Estimateq,
%         thetavq=thetauq+jmp.stdq.*R.z3(:,i);
%         thetavq(thetavq<jmp.qmin)=jmp.qmin;
%         pv3=exp(-0.5.*(thetavq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetavq-Prior.meanq));
%         fv=CalcLklhd(Obs,thetauA0,thetaun,D,Prior,Delta,DeltaA,B,thetavq);    
% 
%         MetRatio=exp(fv-fu)*pv3/pu4;
%         if MetRatio>R.u3(i),
%             C.n_a3=C.n_a3+1; %increment
%             thetauq=thetavq; fu=fv; pu4=pv3; %update u->v     
%         end    
%         C.thetaq(:,i)=thetauq;  
%     end
    
    C.Like(i)=exp(fu);
    C.LogLike(i)=fu;
end

disp(['Elapsed time: ' num2str(toc) ' seconds.'])

disp(['A0: Acceptance rate =' num2str(C.n_a1/C.N*100) ' pct.'])
disp(['na: Acceptance rate =' num2str(C.n_a2/C.N*100) ' pct.'])
disp(['x1: Acceptance rate =' num2str(C.n_a3/C.N*100) ' pct.'])
% if C.Estimateq,
%     disp(['q: Acceptance rate =' num2str(C.n_a3/C.N*100) ' pct.'])
% end

return