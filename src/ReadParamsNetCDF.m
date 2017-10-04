function [Chain,Prior,R,Exp] = ReadParamsNetCDF(fname)

fid=fopen(fname,'r');
%% a bunch of former values in params.txt will be hard-wired constants
Chain.N=1e4; Chain.Nburn=2e3; R.Seed=9; Prior.covQbar=0.5; Exp.nOpt=6;
Prior.meanQbar=ncread(fname,'/River_Info/QWBM');

orbit=round(ncread(fname,'/River_Info/orbit')/86400);
orbitlen=21;
timesteps=size(ncread(fname,'/Reach_Timeseries_Simulator/H'),2);
tUse1=orbit(1);
times_mat=[orbit orbit+[1:floor(timesteps/orbitlen)]*orbitlen];
veclen=size(times_mat,1)*size(times_mat,2);
Exp.tUse=reshape(times_mat,[veclen 1]);
Exp.tStep=1;
% fgetl(fid); Chain.N=fscanf(fid,'%f',1); fscanf(fid,'\n');
% fgetl(fid); Chain.Nburn=fscanf(fid,'%f',1); fscanf(fid,'\n');
% fgetl(fid); R.Seed=fscanf(fid,'%f',1); fscanf(fid,'\n');
%fgetl(fid); tUse1=fscanf(fid,'%f',1); tUseEnd=fscanf(fid,'%f',1); fscanf(fid,'\n');
%fgetl(fid); Prior.meanQbar=ncread(fname,'/River_Info/QWBM');
%fgetl(fid); Exp.nOpt=fscanf(fid,'%f',1); fscanf(fid,'\n');

% if fgetl(fid)~=-1 %if extra line in param file for time step size
%     Exp.tStep=fscanf(fid,'%f',1); fscanf(fid,'\n'); %time step as fraction of day
% else
%     Exp.tStep=1; %default to one day
% end

%Exp.tUse=(round(tUse1/Exp.tStep):1:round(tUseEnd/Exp.tStep))*Exp.tStep;
Exp.Est_nt=length(Exp.tUse);

fclose(fid);

return