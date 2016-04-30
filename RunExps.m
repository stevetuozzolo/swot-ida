clear all

uselib('SWOTQAlg')

ShowFigs=true;
RunBjerklie=false;
ReRunPrior=true; %need to re-run every time jump params change!
RunMetroMan=true;
Smin=1E-5;

fid=fopen('RunFile.txt');
while ~feof(fid),
   RunName=fgetl(fid);
   if strcmp(RunName(1:2),'//')
       disp(['Skipping ' RunName(3:end)])
       continue
   end
   if RunMetroMan,
       RunExp(RunName,ShowFigs,ReRunPrior);
   end
   if RunBjerklie,
       RunBjerklieExp(RunName,ShowFigs)
   end
end
