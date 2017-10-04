function [Chain,Prior,R,Exp] = ReadParams(fname)

fid=fopen(fname,'r');
fgetl(fid); Chain.N=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Chain.Nburn=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); R.Seed=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); tUse1=fscanf(fid,'%f',1); tUseEnd=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.meanQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.covQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.nOpt=fscanf(fid,'%f',1); fscanf(fid,'\n');

Exp.tUse=tUse1:tUseEnd;
Exp.Est_nt=length(Exp.tUse);

fclose(fid);

return