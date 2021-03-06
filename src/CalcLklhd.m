function [f,dQdxv,dAdtv,Cf,Qv]=CalcLklhd(Obs,AllObs,A0,na,x1,D,Prior,Delta,DeltaA,B,qhatv,BjerklienOpt)

%All vectors ordered "space-first"
% theta(1)=theta(r1,t1)
% theta(2)=theta(r1,t2)
% ... 
% theta(nt)=theta(r1,nt)
% theta(nt+1)=theta(r2,t1)

%prep...
N=D.nR*(D.nt-1); %total number of "equations" / constraints
M=D.nR*D.nt;

A0v=reshape((A0*ones(1,D.nt))',D.nR*D.nt,1);

c1=0.85;
for r=1:D.nR,
    nhat(r,:) = calcnhat(Obs.w(r,:),Obs.h(r,:),AllObs.hmin(r),Prior.Wa(r),Prior.Ha(r),c1,x1(r),na(r),BjerklienOpt);
end

nv=reshape(nhat',D.nR*D.nt,1);
% nv=
Qv=1./nv.*(A0v+Obs.dAv).^(5/3).*Obs.wv.^(-2/3).*sqrt(Obs.Sv);

if any(Obs.hv)<0 || any(A0v)<0 || any(Obs.Sv)<0,
    f=0;
    return
end

%1) Calculate dQdx, dQdt, and q for channel mass balance 
dQdxv=Delta*Qv;
dAdtv=(DeltaA*Obs.hv) ./ D.dt .* (B*Obs.wv);

%2) Calculate covariance matrix of theta

%2.1) Calculate covariance matrix of dQdx

%2.1.1) Handle Jacobian of dQdx w.r.t. slope, dA, width
TSv=Obs.Sv.^-1;
TdAv=1./(A0v+Obs.dAv);
Tw=Obs.wv.^-1;

JS=.5.*Delta.*(ones(N,1)*Qv').*(ones(N,1)*TSv');
JdA=5/3.*Delta.*(ones(N,1)*Qv').*(ones(N,1)*TdAv');
Jw=-2/3.*Delta.*(ones(N,1)*Qv').*(ones(N,1)*Tw');


%2.1.3) Covariance matrix calculation
J=[JS JdA Jw];
CdQ=J*Obs.CSdAw*J';

% %2.1.4) Calculate error covariance due to Manning's error
% CdQm=Delta*(eye(M).*Prior.eQm^2)*Delta'; 

%2.3) Final covariance matrix calculation
Cf=Obs.CA+CdQ+Prior.Cqf; %+CdQm;

%3) Calculate likelihood
Theta=dQdxv+dAdtv; %-qhatv; %omit q estimation for now

if isnan(rcond(Cf))
    f=0;
else
    f=-0.5.*Theta'/Cf*Theta;
end



return