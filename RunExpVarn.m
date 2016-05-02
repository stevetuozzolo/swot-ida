% this script runs a single experiment, taking Run Directory as input

function RunExpVarn(RunDir,ShowFigs,ReRunPrior,Smin,BjerklienOpt)

% RunDir='.';

disp(['Running ' RunDir])

ObsFile=[RunDir '/SWOTobs.txt'];
[DAll,Obs] = ReadObs(ObsFile);

ParamFile=[RunDir '/params.txt'];
[Chain,Prior,jmp,R,Exp] = ReadParams(ParamFile,DAll);

TruthFile=[RunDir '/truth.txt'];
AllTruth=ReadTruth (TruthFile,DAll);

[D,Obs,AllObs,DAll,Truth]=SelObs(DAll,Obs,Exp,AllTruth);

[Obs] = CalcdA(D,Obs);
[AllObs] = CalcdA(DAll,AllObs);

if ReRunPrior,
    [Prior,jmp,AllObs]=ProcessPrior(Prior,AllObs,jmp,DAll,Obs,D,ShowFigs,BjerklienOpt); 
    save([ RunDir '/Prior.mat'],'Prior','jmp','AllObs');
else    
    load([ RunDir '/Prior.mat'],'Prior','jmp','AllObs');
end

[Obs,Prior] = GetCovMats(D,Obs,Prior);

%limit slopes to Smin...
Obs.S(Obs.S<Smin)=Smin;
AllObs.S(AllObs.S<Smin)=Smin;

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R,DAll,AllObs,BjerklienOpt);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior,DAll,AllObs,BjerklienOpt);

[Estimate] = FilterEstimate(Estimate,Chain,D,Obs);

Err=CalcErrorStats(AllTruth,Estimate,DAll);

Err=DispRMSEStats(Err,Truth,Prior,Estimate);

if ShowFigs,
    MakeFigs(D,Truth,Prior,Chain,Estimate,Err,AllTruth,DAll);
end

% WriteSummary (R,Err,Estimate,RunDir);

% CheckBals(Truth,Obs,D,Prior,Chain,R,Estimate)

save([RunDir '/RunData.mat'])

return