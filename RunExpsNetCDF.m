%clear all

%uselib('SWOTQAlg')
addpath src

ShowFigs=true;
RunBjerklie=false;
RunMetroMan=true;

fid=fopen('RunFile.txt');
while ~feof(fid),
   RunName=fgetl(fid);
   if strcmp(RunName(1:2),'//')
       disp(['Skipping ' RunName(3:end)])
       continue
   end
   if RunMetroMan,
       RunExpNetCDF(RunName,ShowFigs);
   end
   if RunBjerklie,
       RunBjerklieExp(RunName,ShowFigs)
   end
end
