function [choiceModelStats, PFxdata, PFydata, freqBias, freqSensitivity] = logistic_multiFreq(performance, binaryResponseSide, stimFreq, categoryBoundary);

choicePredictors = []; sideChosen = [];
for n = 1:length(performance)
    if ~isnan(performance(n))
    choicePredictors(end+1) = stimFreq(n) - categoryBoundary;
    sideChosen(end+1) = binaryResponseSide(n);
    end
end

[choiceCoef, ~, choiceModelStats] = glmfit(choicePredictors',sideChosen','binomial','link','logit');

xdata = min(stimFreq)-categoryBoundary:0.5:max(stimFreq)-categoryBoundary;
PFydata = 1./(1+exp(-(choiceCoef(1) + choiceCoef(2)*xdata)));
PFxdata = xdata+12;

mu = (log(1) - choiceCoef(1))/choiceCoef(2);
% myXval = 0;
% s = (mu- myXval)/(choiceCoef(1) + choiceCoef(2)*myXval);
% sigma = sqrt((s^2*pi^2)/3);



freqBias = categoryBoundary + (log(1) - choiceCoef(1))/choiceCoef(2);
freqSensitivity = abs(choiceCoef(2)/4); %According to the divide by 4 rule when evaluating slope at the steepest point.
%Solve for the case of p = 0.5 , only use positive values 


end