clear all

uselib('SWOTQAlg')

ShowFigs=true;
RunBjerklie=false;
ReRunPrior=false; %need to re-run every time jump params change!
RunMetroMan=true;
Smin=1E-5; %should move this to parameter file

%TEMP -- need to move to parameter file
BjerklienOpt=2; %1 = use standard form; 2 = use height-only form

fid=fopen('RunFile.txt');
while ~feof(fid),
   RunName=fgetl(fid);
   if strcmp(RunName(1:2),'//')
       disp(['Skipping ' RunName(3:end)])
       continue
   end
   if RunMetroMan,
%        RunExp(RunName,ShowFigs,ReRunPrior);
        RunExpVarn(RunName,ShowFigs,ReRunPrior,Smin,BjerklienOpt); %variable roughness!
   end
   if RunBjerklie,
       RunBjerklieExp(RunName,ShowFigs);
   end
end
