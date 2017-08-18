% this script runs a single experiment, taking Run Directory as input

function RunExpVarn(RunDir,ShowFigs,ReRunPrior,Smin)

% RunDir='.';

disp(['Running ' RunDir])

ObsFile=[RunDir '/SWOTobs.txt'];
[DAll,Obs] = ReadObs(ObsFile);

ParamFile=[RunDir '/params.txt'];
[Chain,Prior,R,Exp] = ReadParams(ParamFile);

TruthFile=[RunDir '/truth.txt'];
AllTruth=ReadTruth (TruthFile,DAll);

[D,Obs,AllObs,DAll,Truth]=SelObs(DAll,Obs,Exp,AllTruth);

[Obs] = CalcdA(D,Obs);
[AllObs] = CalcdA(DAll,AllObs);

%limit slopes to Smin... should just move this above the SelObs...
Obs.S(Obs.S<Smin)=Smin;
AllObs.S(AllObs.S<Smin)=Smin;

AllObs.hmin=min(AllObs.h,[],2);

if ReRunPrior,
    [Prior,jmp,AllObs]=ProcessPrior(Prior,AllObs,DAll,Obs,D,ShowFigs,Exp.BjerklienOpt); 
    save([ RunDir '/Prior.mat'],'Prior','jmp','AllObs');
else    
    load([ RunDir '/Prior.mat'],'Prior','jmp','AllObs');
end

[Obs,Prior] = GetCovMats(D,Obs,Prior);

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R,DAll,AllObs,Exp.BjerklienOpt);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior,DAll,AllObs,Exp.BjerklienOpt);

[Estimate] = FilterEstimate(Estimate,Chain,D,Obs);

Err=CalcErrorStats(AllTruth,Estimate,DAll);

Err=DispRMSEStats(Err,Truth,Prior,Estimate);

if ShowFigs,
    MakeFigs(D,Truth,Prior,Chain,Estimate,Err,AllTruth,DAll,AllObs);
end

% WriteSummary (R,Err,Estimate,RunDir);

% CheckBals(Truth,Obs,D,Prior,Chain,R,Estimate)

save([RunDir '/RunData.mat'])

return