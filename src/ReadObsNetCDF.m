function [Domain,Obs] = ReadObsNetCDF(fname)

gdrch=ncread(fname,'/River_Info/gdrch');
rchbnd=ncread(fname,'/River_Info/rch_bnd');
Domain.nR=length(ncread(fname,'/River_Info/gdrch'));
if strcmp(fname(end-6:end-3),'sity')
    groupname='/Reach_Timeseries_Simulator/';
else
    groupname='/Reach_Timeseries/';
end

for i=1:length(gdrch)
    Domain.xkm(i)=(rchbnd(gdrch(i)+1)+rchbnd(gdrch(i)))/2/1e3;
    Domain.L(i)=(rchbnd(gdrch(i)+1)-rchbnd(gdrch(i)))/1e3;    
end

allH=ncread(fname,[groupname 'H']);
allS=ncread(fname,[groupname 'S'])./1E5;
allW=ncread(fname,[groupname 'W']);

Domain.nt=size(allH,2);
Domain.t=ncread(fname,'River_Info/orbit')'/86400;
Domain.dt=reshape( (diff(Domain.t)'.*86000*ones(1,Domain.nR)),Domain.nR*(Domain.nt-1),1);

Obs.h=allH(gdrch,:);
Obs.h0=Obs.h(:,1);
Obs.S=allS(gdrch,:);
Obs.w=allW(gdrch,:);

Obs.sigS=1e-5;
Obs.sigh=10;
Obs.sigw=1;

Obs.hv=reshape(Obs.h',Domain.nR*Domain.nt,1);
Obs.wv=reshape(Obs.w',Domain.nR*Domain.nt,1);
Obs.Sv=reshape(Obs.S',Domain.nR*Domain.nt,1);

return