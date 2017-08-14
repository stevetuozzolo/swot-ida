%% OlentangyTester.m
% rewrite params.txt file in Olentangy.m to determine RMSEs + error stats
% for different initial conditions. Save output files of successful
% algorithm retrievals
Qrng=[1:20];
for i=1:length(Qrng)
    ReWriteParams(Qrng(i))
    RunExpsVarn;
    load Olentangy\RunData
    if isstruct(Prior)
        save(['Olentangy\Opt4_VarQ\Q_' num2str(Qrng(i))])
    end
end

%%
for i=1:length(Qrng)
    load(['Olentangy\Q_' num2str(Qrng(i))]);
    if isstruct(Prior)
        rel_err(i,1)=mean(Err.QRelErrPrior);
        rel_err(i,2)=mean(Err.QRelErrPost);
    end
end
