clear all

addpath ./src

ShowFigs=true;
RunBjerklie=false;
ReRunPrior=true; %need to re-run every time jump params change! just use for debugging etc.
RunMetroMan=true;
Smin=1E-5; %should move this to parameter file

%TEMP -- need to move to parameter file
BjerklienOpt=3; %1 = use standard form; 2 = use height-only form; 3 = use rating-curve form

fid=fopen('RunFile.txt');
while ~feof(fid),
   RunName=fgetl(fid);
   if strcmp(RunName(1:2),'//')
       disp(['Skipping ' RunName(3:end)])
       continue
   end
   if RunMetroMan,
        RunExpVarn(RunName,ShowFigs,ReRunPrior,Smin,BjerklienOpt); %variable roughness!
   end
   if RunBjerklie,
       RunBjerklieExp(RunName,ShowFigs);
   end
end
