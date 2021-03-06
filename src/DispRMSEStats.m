function [Err]=DispRMSEStats (Err,Truth,Prior,E)

%Quick Error stats
Err.RelErrA0=mean(E.A0hat'-Truth.A0)./mean(Truth.A0);
% Err.RelErrN=mean(E.nhat'-Truth.n')./mean(Truth.n);
% RMSq=sqrt(mean( (E.qhat-Truth.q).^2 ));
% RMSqPrior=sqrt(mean( (Prior.meanq-Truth.q).^2 ));

disp(['Relative Error in A0:' num2str(Err.RelErrA0)])
disp(['Relative Uncertainty in A0: ' num2str(mean(E.stdA0Post'./E.A0hat'))])

% disp(['Relative Error in n:' num2str(Err.RelErrN)])
% disp(['Relative Uncertainty in n: ' num2str(mean(E.stdnPost'./E.nhat'))])

% disp(['RMS for q posterior:' num2str(RMSq)])
% disp(['Relative Uncertainty in q: ' num2str(mean(E.stdqpost./E.qhat))])

RMSQPost=sqrt(mean( (E.QhatPostf'-Truth.Q').^2 ));
RMSQPrior=sqrt(mean( (E.QhatPrior'-Truth.Q').^2 ));
ratio=((RMSQPrior-RMSQPost)./RMSQPrior)*100;

nR=length(E.A0hat);

fprintf('RMS for Q prior: ');
for i=1:nR,
    fprintf('%.1f,', RMSQPrior(i))
end
fprintf('\n');

fprintf('RMS for Q posterior: ')
for i=1:nR,
    fprintf('%.1f, ', RMSQPost(i))
end

fprintf('\n')

Err.QRelErrPrior=RMSQPrior./mean(Truth.Q');
Err.QRelErrPost=RMSQPost./mean(Truth.Q');

fprintf('Average RMS for Q posterior: %.3f\n', mean(Err.QRelErrPost))

fprintf('Average relative Q uncertainty: %.3f\n', ...
    mean(mean(E.QstdPost./E.QhatPost)) )

fprintf('For entire timeseries, rRMSE= %.2f, and relative bias= %.2f\n',[Err.Stats.rRMSE Err.Stats.meanRelRes]);

fprintf('For entire timeseries, RMSE/Qbart= %.2f, and bias/Qbart= %.2f\n',[Err.Stats.RMSE/Err.Stats.Qbart Err.Stats.bias/Err.Stats.Qbart]);

fprintf('\n')

return