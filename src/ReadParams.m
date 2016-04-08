function [Chain,Prior,jmp,R,Exp] = ReadParams(fname,D)

fid=fopen(fname,'r');
fgetl(fid); Chain.N=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Chain.Nburn=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); R.Seed=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Chain.Estimateq=fscanf(fid,'%f',1); fscanf(fid,'\n');

fgetl(fid); Exp.Est_nt=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.tUse=fscanf(fid,'%f',Exp.Est_nt)'; fscanf(fid,'\n');

fgetl(fid); Prior.meanA0=fscanf(fid,'%f',D.nR); fscanf(fid,'\n');
fgetl(fid); Prior.stdA0=fscanf(fid,'%f',D.nR); fscanf(fid,'\n');
fgetl(fid); Prior.meann=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.stdn=fscanf(fid,'%f',1); fscanf(fid,'\n');

fgetl(fid); Prior.meanq=fscanf(fid,'%f',D.nR*(D.nt-1)); fscanf(fid,'\n');
fgetl(fid); Prior.stdq=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); jmp.stdA0=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); jmp.rhoA0=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); jmp.stdn=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); jmp.stdq=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); jmp.qmin=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.meanQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.covQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fclose(fid);

return