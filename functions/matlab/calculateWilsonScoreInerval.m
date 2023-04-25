function [lowerBound, upperBound] = calculateWilsonScoreInerval(data);
% [lowerBound, upperBound] = calculateWilsonScoreInerval(data);
%
% Calculate the Wilson score confidence interval for binomial
% distributions. This confidence intervall as asymmetric, bounded by 0 and
% 1 and "pulls" towards 0.5. Currently, this function only handles an alpha
% level of 0.05.
% Fromula was drawn from:
% https://www.ucl.ac.uk/english-usage/staff/sean/resources/binomialpoisson.pdf
%
% INPUT: -data: A vector of 0s and 1s containing the observed outcomes.
%               Nans can be tolerated.
%
% OUTPUTS: -lowerBound/upperBound: The lower and upper bound of the
%                                  confidence interval.
%
% LO, 5/4/2021
%
%--------------------------------------------------------------------------
p = nanmean(data); %take the average fo all the specified values
n = sum(~isnan(data)); %number of valid observations
z = 1.95996; %the critical value for an alpha level of 0.05

lowerBound = (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/(1 + z^2/n);
upperBound = (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/(1 + z^2/n);

end
