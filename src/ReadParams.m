function [Chain,Prior,R,Exp] = ReadParams(fname)

fid=fopen(fname,'r');
fgetl(fid); Chain.N=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Chain.Nburn=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); R.Seed=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.Est_nt=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.tUse=fscanf(fid,'%f',Exp.Est_nt); tUseEnd=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.meanQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.covQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.nOpt=fscanf(fid,'%f',1); fscanf(fid,'\n');

if isempty(Exp.nOpt)
    Exp.nOpt=3;
end

if fgetl(fid)~=-1 %if extra line in param file for time step size
    Exp.tStep=fscanf(fid,'%f',1); fscanf(fid,'\n'); %time step as fraction of day
    Exp.tUse=(round(tUse1/Exp.tStep):1:round(tUseEnd/Exp.tStep))*Exp.tStep;
else
    Exp.tStep=1; %default to one day
end

%Exp.Est_nt=length(Exp.tUse); %move this up to be in params.txt

fclose(fid);

return
