function dPrime = dPrime2AFC(correctSide, responseSide);
%dPrime = dPrime2AFC(correctSide, responseSide);
%
% Calculate the d' score for a 2AFC task.
% 
% INPUTS: -correctSide: The assigned correect side
%         -responseSide: The side chosen by the animal (NaN tolerated)
%
% OUTPUT: -dPrime: The dPrime score, currently without 1/2^(1/2)
%                  correction.
%
% LO, 5/5/2021
%--------------------------------------------------------------------------

%Input check
correctSide = correctSide(~isnan(responseSide)); %Make sure to only consider valid trials
responseSide = responseSide(~isnan(responseSide));

%Determine the side to look at
sideCodes = unique(responseSide); %To make it neutral to the input code
theSide = sideCodes(1); %First take the lower one (usually left)


hitRate = sum(correctSide == theSide & responseSide == theSide)/sum(correctSide == theSide);
%Is defined as the proportion the animal responded correctly with left side
%when left side was true (so, out of all left side is true)
falseAlarmRate = sum(correctSide ~= theSide & responseSide == theSide)/sum(correctSide ~= theSide);
%Is defined as the proportion the animal responded incorrectly with left
%side when the right side was true (so, out of all right side is true)

dPrime = norminv(hitRate) - norminv(falseAlarmRate); %Norminv is the inverse normal distribution fucntion
%dPrime = (1/sqrt(2))*(norminv(hitRate) - norminv(falseAlarmRate)); %This
%is when comapring a 2AFC to a detection task.
% 
% % Handle degenerate cases
% if dPrime == Inf || dPrime == -Inf %Check for misclassifications
%     if ~(hitRate == 1 && falseAlarmrate == 0) %If it is not the ideal observer
%       if hitrate == 1 %The 
%           dPrime =  norminv(1-falseAlaramRate);
%        
%     if (hitRate-falseAlarmRate) < 0.5 %This is only the case when animals are so biased that they only choose one side!
%         dPrime = 0; %Force 0 discrimination
%     end
% end
