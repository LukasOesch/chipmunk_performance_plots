function [performanceMetrics, figureHandle,outputLoc] = performanceOverviewDecisionTask(sessionFile,outputLoc, sessionsBack);
%
%
%
%
%--------------------------------------------------------------------------
%Find the directory where the running function is located and add the
%functions folder
function_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(function_dir,'functions')); %Sorry for the abusive naming...
%
performanceMetrics = [];
if ~exist('sessionsBack') || isempty(sessionsBack)
sessionsBack = 5; %how many sessions back should be analysed
end
trialTimeBinSize = 4; %minutes

%Convenience
if ~exist('sessionFile') || isempty(sessionFile) %let the user pick a file if nothing has been specified
   [fileName, fileFolder] = uigetfile('*.mat','Select the session');
   sessionFile = fullfile(fileFolder, fileName);
   [fileFolder, fileName, fileExt] = fileparts(sessionFile);
   [~, parentDir] = fileparts(fileFolder); %Get the name of the parent directory
else
   [fileFolder, fileName, fileExt] = fileparts(sessionFile); %Get the location of the specified file
   [~, parentDir] = fileparts(fileFolder); %Get the name of the parent directory
    %This relies on all the data being stored in one folder together
    %instead of within distinct folders!
end

%List the files inside the containing directory
listFiles = dir(fullfile(fileFolder,['*' fileExt(2:end)])); %Get all the files in the folder
%Unfortunately dir(*.ext) does not return an exact match in the style of
%strcmp but rather of the sort of strfind. That is why .obsmat files are
%mixed into the list when looking for .mat...

%Determine whether the chipmunk data live in the churchland lab directory
%structure or whether they are all grouped together in one folder per mouse
if strcmp(parentDir, 'chipmunk') & (length(listFiles)==1)
    parts = strsplit(fileFolder, filesep);
    sessionDate = parts{end-1}; %Retrieve the session data
    
    %Get the directory containing the session directories
    if ispc %Sensitive to operating system
        sessionsFolder = [];
        for k = 1:length(parts) - 2
            sessionsFolder = fullfile(sessionsFolder, parts{k});
        end
    else
        sessionsFolder = filesep;
        for k = 2:length(parts) - 2
            sessionsFolder = fullfile(filesep, sessionsFolder, parts{k});
        end
    end
    
    dirList = struct2cell(dir(sessionsFolder)); %Where rows are: name, folder, date, bytes, isdir and datenum
    directories = dirList(1, cell2mat(dirList(5,:)) == 1); %Find the entries that are directories and retrieve folder names for these only
    %Sort the directories in ascending order
    directories = flip(sort(directories));
    
    %Locate the selected session in the session files
    for k=1:length(directories)
        if strcmp(directories{k}, sessionDate)
            sessionIdx = k;
        end
    end
    
    includeFiles = [];
    temp = dir(fullfile(sessionsFolder, directories{sessionIdx}, 'chipmunk', ['*' fileExt(2:end)]));
    includeFiles{1} = fullfile(temp.folder, temp.name);
    for k=1:sessionsBack
        temp = dir(fullfile(sessionsFolder, directories{sessionIdx + k}, 'chipmunk',['*' fileExt(2:end)]));
        includeFiles{k+1} = fullfile(temp.folder, temp.name);
    end
    includeFiles = flip(includeFiles); %Reverse the order so that later sessions appear on higher indices
    
    
    %When the input is inside one folder
else
    
    %Ugly and bulky, but works..
    delIdx = [];
    for k=1:length(listFiles)
        [~,~, ext] = fileparts(listFiles(k).name);
        if ~strcmp(ext, fileExt)
            delIdx(end+1) = k;
        end
    end
    listFiles(delIdx) = [];
    
    %%%%--------Sort the files according to the date in their name-----------
    sessionFileDateTime = []; %Get the date and time of the files in their name to sort  them
    for k = 1:size(listFiles,1)
        pos = strfind(listFiles(k).name,'202'); %Look for the 20ies.
        if ~isempty(pos)
            if length(pos) > 1
                dateStart = [];
                for n=1:length(pos)
                    if  length(listFiles(k).name)-pos(n) >= 7 %Make sure to pick up the entire date if it is one!
                        if ~isempty(str2num(listFiles(k).name(pos(n):pos(n)+7))) %Check whether convertible to number
                            dateStart = n;
                        end
                    end
                end
                if length(dateStart) == 1
                    pos = pos(dateStart);
                else
                    error('Found multiple 7-digit number inside you file name, date could not be extracted.')
                end
            end
            sessionFileDateTime(k) = datenum(datetime(listFiles(k).name(pos:pos+14),'InputFormat','yyyyMMdd_HHmmss'));
        else
            sessionFileDateTime(k) = 0;
        end
    end
    [~,fidx] = sort(sessionFileDateTime,'descend');
    listFiles = listFiles(fidx);
    
    %%%---- find your file and load the outones obtained x sessions before
    myFileIdx = 0;
    for k = 1:size(listFiles,1)
        if strcmp(listFiles(k).name, [fileName fileExt])
            myFileIdx = k;
        end
    end
    
    % for k=1:size(listFiles,1)
    %     if ~isempty(strfind(listFiles(k).name, fileName))
    %         sessionIdx = k;
    %     end
    % end
    
    includeFiles = []; %Check which files need to be loaded according to how many past sessions should be plotted
    for s = 1:sessionsBack+1
        includeFiles{s} = fullfile(fileFolder,listFiles(myFileIdx+(sessionsBack+1-s)).name);
    end
end

if ~exist('outputLoc') || isempty(outputLoc) %let the user pick a file if nothing has been specified
   outputLoc = uigetdir(fileFolder,'Select the output folder');
end

%--------------------------------------------------------------------------
%% Start with the within session analysis
%data to be displayed within session potentially
experimentName = cell(1,sessionsBack+1); %Type of experiment defined by a state matrix
dates = cell(1,sessionsBack+1); %The dates when the sessions were run
motionTime = cell(1,sessionsBack+1); %The time the animal spends moving from the center to the response port
absoluteWaitTime = cell(1,sessionsBack+1); %The effective time the animal waits inside the center port
targetWaitTime = cell(1,sessionsBack+1); %The required time to wait in the center
performance = cell(1,sessionsBack+1); %Here, a vector of correct = 1, error = 0 and early withdrawal and no choice = NaN
earlyWithdrawals = cell(1,sessionsBack+1); %Vector of early withdrawals = 1, rest = 0
gapTime = cell(1,sessionsBack+1); %The time between the offset of the stimulus and the animal reporting its decision
choiceModelStats = cell(1, sessionsBack+1); %The stats for the GLM fit for choice history and stimulus
modelFitWarnings = cell(1, sessionsBack+1); %Warnings for unconverged model fits or perfect separation between cases
stimFreq = cell(1, sessionsBack+1);
trialDuration = cell(1, sessionsBack+1); %The duration between initiation of the trial and outcome (Reward, Punishment, Early Withdrawal, noChoice)

%data summarized across sessions
numTrials = []; %
noChoice = [];
extendedStim = []; %Extended stimulus presentation time
totalLeftRewards = [];
totalRightRewards = [];
totalReward = [];
bodyWeight = [];
leftwardBias = []; %A side preference metric based on a given side being falsely normalized to the sum of erroneous choices on both sides, further divided by the anmal's performance
choiceCoef = NaN(3,sessionsBack+1); %The choice history model coefficients
performanceEasiest = [];
%--------------------------------------------------------------------------
% Start looping
for s = 1:sessionsBack+1
load(includeFiles{s});

newChipmunk = false;
if isfield(SessionData.TrialSettings(1),'smaAssembler'); %Check whether this is the newer implementatio of chipmunk from June 2021
    newChipmunk = true;
end

%%%%%%%%%%TEMPORARY%%%%%%%%%%%%%%
if newChipmunk
    if isfield(SessionData,'TaskPase')
        SessionData.TaskPhase = SessionData.TaskPase;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine the phase of the learning to display the data accordingly but
%also to process the data adequately
if isfield(SessionData,'TaskPhase')
    experimentName{s} = SessionData.TaskPhase;
elseif SessionData.TrialSettings(1).initCenterStimUntilCorrect == 1;
    experimentName{s} = 'initCenterStimUntilCorrect';
elseif SessionData.TrialSettings(1).holdCenterStimUntilCorrect == 1;
    experimentName{s} = 'holdCenterStimUntilCorrect';
else
    experimentName{s} = 'Habituation';
end

dates{s} = string(datetime(datestr(SessionData.Info.SessionDate),'Format','MM/dd'));
% try
% dates{s} = string(datetime(includeFiles{s}(end-18:end-11),'InputFormat','yyyyMMdd','Format','MM/dd'));
% catch
%     dates{s} = string(datetime(includeFiles{s}(7:14),'InputFormat','yyyyMMdd','Format','MM/dd'));
% end
%--------------------------------------------------------------------------
%General values of interest

if isfield(SessionData.TrialSettings,'BodyWeight')
    bodyWeight(s) = SessionData.TrialSettings(end).BodyWeight;
elseif isfield(SessionData.TrialSettings,'demonWeight')
    if ~isempty(SessionData.TrialSettings(end).demonWeight)
        if ~isempty(str2num(num2str(SessionData.TrialSettings(end).demonWeight))) %This can happen when the user inserts a string
            bodyWeight(s) = SessionData.TrialSettings(end).demonWeight;
        else
            bodyWeight(s) = NaN;   
        end
    else
        bodyWeight(s) = NaN;
    end
else
    bodyWeight(s) = NaN;
end

if newChipmunk
    numTrials(s) = sum(SessionData.OutcomeRecord > -2); %Initiated
else
    numTrials(s) = SessionData.nTrials;
end
noChoice(s) = sum(SessionData.DidNotChoose);
extendedStim(s) = mean(SessionData.ExtraStimDuration);
if newChipmunk
fractionLeftAssigned(s) = mean(SessionData.CorrectSide == 0); %Note to change if sides change to 0 and 1
else
fractionLeftAssigned(s) = mean(SessionData.CorrectSide == 1); %Note to change if sides change to 0 and 1
end

earlyWithdrawals{s} = SessionData.EarlyWithdrawal;
%--------------------------------------------------------------------------
%Reward parameters
leftRewards = []; rightRewards = [];
for k = 1:SessionData.nTrials
    if SessionData.Rewarded(k) == 1
        if SessionData.CorrectSide(k) == 1
            leftRewards(end+1) = SessionData.TrialSettings(k).leftRewardVolume;
        else
            rightRewards(end+1) = SessionData.TrialSettings(k).rightRewardVolume;
        end
    end
end
totalLeftRewards = sum(leftRewards);
totalRightRewards = sum(rightRewards);
totalReward(s) = totalLeftRewards + totalRightRewards;

%--------------------------------------------------------------------------
%Get the correct choices and the side chosen for different experimental
%condition
correctChoice = []; responseSide = [];
if newChipmunk
    correctChoice = nan(1,SessionData.nTrials);
    correctChoice(SessionData.OutcomeRecord==1) = 1;
    correctChoice(SessionData.OutcomeRecord==0) = 0;
     responseSide = SessionData.ResponseSide;
elseif ~isempty(strfind(SessionData.TaskPhase, 'holdCenterStimUntilCorrect')) || ~isempty(strfind(SessionData.TaskPhase,'Fixation'))
    for n=1:SessionData.nTrials
        P1visits = []; P3visits = [];
        if SessionData.EarlyWithdrawal(n) == 0 & SessionData.DidNotChoose(n) == 0
            if isfield(SessionData.RawEvents.Trial{1,n}.Events, 'Port1In')
                P1visits = SessionData.RawEvents.Trial{1,n}.Events.Port1In(find(SessionData.RawEvents.Trial{1,n}.Events.Port1In > SessionData.RawEvents.Trial{1,n}.States.WaitCenter(1)));
            end
            if isfield(SessionData.RawEvents.Trial{1,n}.Events, 'Port3In')
                P3visits = SessionData.RawEvents.Trial{1,n}.Events.Port3In(find(SessionData.RawEvents.Trial{1,n}.Events.Port3In > SessionData.RawEvents.Trial{1,n}.States.WaitCenter(1)));
            end
            if isempty(P3visits)
                responseSide(n) = 1;
                correctChoice(n) = 1;
            elseif isempty(P1visits)
                responseSide(n) = 2;
                correctChoice(n) = 1;
            else
                responseSide(n) = (P1visits(1) > P3visits(1))+1; % the case when the animal changes its mind after trying
                correctChoice(n) = responseSide(n) == SessionData.CorrectSide(n);
            end
        else
            responseSide(n) = NaN;
            correctChoice(n) = NaN;
        end
    end
else
    responseSide = SessionData.ResponseSide;
    correctChoice = SessionData.Rewarded;
    correctChoice(SessionData.EarlyWithdrawal == 1 | SessionData.DidNotChoose == 1) = NaN;
end

%--------------------------------------------------------------------------
% Get within session performance metrics
performance{s} = correctChoice;

if newChipmunk
wrongChoiceInd = find(SessionData.OutcomeRecord == 0);
else
wrongChoiceInd = find(correctChoice == 0 & SessionData.EarlyWithdrawal == 0 & SessionData.DidNotChoose ==0);
end
if newChipmunk
    wrongOnLeftTrial = sum(SessionData.CorrectSide(wrongChoiceInd) == 0)/sum(SessionData.CorrectSide == 0);
wrongOnRightTrial = sum(SessionData.CorrectSide(wrongChoiceInd) == 1)/sum(SessionData.CorrectSide == 1);
else
wrongOnLeftTrial = sum(SessionData.CorrectSide(wrongChoiceInd) == 1)/sum(SessionData.CorrectSide == 1);
wrongOnRightTrial = sum(SessionData.CorrectSide(wrongChoiceInd) == 2)/sum(SessionData.CorrectSide == 2);
end
%leftwardBias(s) = (wrongOnRightTrial/(wrongOnLeftTrial+wrongOnRightTrial)-0.5)*2;
leftwardBias(s) = (wrongOnRightTrial/(wrongOnLeftTrial+wrongOnRightTrial)-0.5)/nanmean(performance{s});

%--------------------------------------------------------------------------
% Get the aboulute wait time and the set wait times to compute the
% difference, the motion speed and the reaction time...
% ....also ge the stim frequency and the trial duration as we are already
% looping through all the trials
if newChipmunk
    for n = 1:SessionData.nTrials
        if ~isnan(SessionData.RawEvents.Trial{n}.States.PlayStimulus(1)) %General case where the mouse waits until real or virtual observer pokes
            targetWaitTime{s}(end+1) = SessionData.RawEvents.Trial{n}.States.PlayStimulus(1)- SessionData.RawEvents.Trial{n}.States.DemonInitFixation(1) + SessionData.SetWaitTime(n);
        else
            if isfield(SessionData.TrialDelays, 'virtualObsInitDelay') %Early withdrawals before play stimulus when learning to wait for a virtual observer
                targetWaitTime{s}(end+1) = SessionData.SetWaitTime(n) + SessionData.PreStimDelay(n) + SessionData.TrialDelays(n).virtualObsInitDelay;
            elseif isfield(SessionData, 'ObsOutcomeRecord') %Early withdrawals happpening before a real observer pokes, assume obsInitiationWindow as required
                targetWaitTime{s}(end+1) = SessionData.SetWaitTime(n) + SessionData.PreStimDelay(n) + SessionData.TrialSettings(n).obsInitiationWindow;
            else %For early withdrawals in older implementations
                targetWaitTime{s}(end+1) = SessionData.SetWaitTime(n) + SessionData.PreStimDelay(n);
            end
        end
        absoluteWaitTime{s}(end+1) = SessionData.ActualWaitTime(n);
        %stimFreq{s}(end+1) = SessionData.StimulusRate(n,SessionData.Modality(n));
        if ~isnan(SessionData.StimulusRate(n,1)) & isnan(SessionData.StimulusRate(n,2)) %Visual only
            stimFreq{s}(end+1,:) = [SessionData.StimulusRate(n,1) nan nan];
        elseif isnan(SessionData.StimulusRate(n,1)) & ~isnan(SessionData.StimulusRate(n,2)) %Auditory only
            stimFreq{s}(end+1,:) = [nan SessionData.StimulusRate(n,2) nan];
        elseif ~isnan(SessionData.StimulusRate(n,1)) & ~isnan(SessionData.StimulusRate(n,2)) %Audiovisual 
            stimFreq{s}(end+1,:) = [nan nan SessionData.StimulusRate(n,1)];
        end
        if SessionData.ValidTrials(n) == 1
            
            %%%%%%%%%%TEMPORARY%%%%%%%%%%%%%%
            if ~isfield(SessionData.RawEvents.Trial{n}.States, 'DemonWaitForWithdrawalFromCenter')
                idx = min(find(SessionData.RawEvents.Trial{n}.Events.Port2Out > SessionData.RawEvents.Trial{n}.States.DemonWaitForResponse(1)));
                motionTime{s}(end+1) = SessionData.RawEvents.Trial{n}.States.DemonWaitForResponse(2)- SessionData.RawEvents.Trial{n}.Events.Port2Out(idx);
                gapTime{s}(end+1) = (SessionData.RawEvents.Trial{n}.States.DemonWaitForResponse(2) - SessionData.RawEvents.Trial{n}.States.PlayStimulus(2)) - (SessionData.ExtraStimDuration(n) + SessionData.TotalStimDuration(n));
                if gapTime{s}(end) < 0
                    gapTime{s}(end) = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                motionTime{s}(end+1) = SessionData.RawEvents.Trial{n}.States.DemonWaitForResponse(2)- SessionData.RawEvents.Trial{n}.States.DemonWaitForWithdrawalFromCenter(2);
                gapTime{s}(end+1) = (SessionData.RawEvents.Trial{n}.States.DemonWaitForResponse(2) - SessionData.RawEvents.Trial{n}.States.PlayStimulus(2)) - (SessionData.ExtraStimDuration(n) + SessionData.TotalStimDuration(n));
                if gapTime{s}(end) < 0
                    gapTime{s}(end) = 0;
                end
            end
        else
            motionTime{s}(end+1) = NaN;
            gapTime{s}(end+1) = NaN;
        end
    end
    
else
    
    for n = 1:SessionData.nTrials
        targetWaitTime{s}(end+1) = SessionData.SetWaitTime(n) + SessionData.PreStimDelay(n);
        if SessionData.EarlyWithdrawal(n) == 0
            absoluteWaitTime{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.WaitForWithdrawalFromCenter(2) - SessionData.RawEvents.Trial{n}.States.PlayWaitStartCue(1);
            if SessionData.DidNotChoose(n) == 0
                if ~isempty(strfind(experimentName{s}, 'Task'))
                    motionTime{s}(end+1) = SessionData.RawEvents.Trial{n}.States.WaitForResponse(2)- SessionData.RawEvents.Trial{n}.States.WaitForResponse(1);
                    gapTime{s}(end+1) = (SessionData.RawEvents.Trial{n}.States.WaitForResponse(2) - SessionData.RawEvents.Trial{n}.States.PlayStimulus(2)) - (SessionData.ExtraStimDuration(n) + SessionData.TotalStimDuration(n)/1000);
                    if gapTime{s}(end) < 0
                        gapTime{s}(end) = 0;
                    end
                    
                elseif ~isempty(strfind(experimentName{s}, 'holdCenterStimUntilCorrect')) %check the movement time for the first response
                    timeOfResponse = [];
                    if responseSide(n) == 1 %left response
                        timeOfResponse = min(SessionData.RawEvents.Trial{1,n}.Events.Port1In(SessionData.RawEvents.Trial{1,n}.Events.Port1In > SessionData.RawEvents.Trial{1,n}.States.WaitCenter(1)));
                    else  %right response
                        timeOfResponse = min(SessionData.RawEvents.Trial{1,n}.Events.Port3In(SessionData.RawEvents.Trial{1,n}.Events.Port3In > SessionData.RawEvents.Trial{1,n}.States.WaitCenter(1)));
                    end
                    motionTime{s}(end+1) =  timeOfResponse - SessionData.RawEvents.Trial{n}.States.WaitForResponse(1);
                    gapTime{s}(end+1) = (timeOfResponse - SessionData.RawEvents.Trial{n}.States.PlayStimulus(2)) - (SessionData.ExtraStimDuration(n) + SessionData.TotalStimDuration(n)/1000);
                    if gapTime{s}(end) < 0
                        gapTime{s}(end) = 0;
                    end
                end
            else
                motionTime{s}(end+1) = NaN;
                gapTime{s}(end+1) = NaN;
            end
        else
            absoluteWaitTime{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.EarlyWithdrawal(1) - SessionData.RawEvents.Trial{n}.States.PlayWaitStartCue(1);
            motionTime{s}(end+1) = NaN;
            gapTime{s}(end+1) = NaN;
        end
        stimFreq{s}(end+1) = length(SessionData.stimEventList{n})-1; %Access the list of the stimulus times
    end
end

%--------------------------------------------------------------------------
%Loop again for trial time on different outcomes
%For differentiating between outcomes
if newChipmunk
for n = 1:SessionData.nTrials
    if SessionData.OutcomeRecord(n) > -2 && SessionData.OutcomeRecord(n) < 2 %All the trials was an initiation except for no choice trials
        try %For the revised chipmunk version
        trialDuration{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.DemonWaitForResponse(2) - SessionData.RawEvents.Trial{1,n}.States.DemonTrialStart(1);
        catch
        trialDuration{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.DemonWaitForResponse(2) - SessionData.RawEvents.Trial{1,n}.States.PlayStimulus(1);   
        end
    end
end
    
else       
    
for n = 1:SessionData.nTrials
    if ~isnan(SessionData.RawEvents.Trial{1,n}.States.Reward(1))
        trialDuration{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.Reward(1) - SessionData.RawEvents.Trial{1,n}.States.GoToCenter(2);
    elseif ~isnan(SessionData.RawEvents.Trial{1,n}.States.WrongChoice(1))
        trialDuration{s}(end+1) = SessionData.RawEvents.Trial{1,n}.States.WrongChoice(1) - SessionData.RawEvents.Trial{1,n}.States.GoToCenter(2);
    else
        trialDuration{s}(end+1) = NaN;
    end
end
end
    
%---------------------------------------------------------------------------
%Construct predictor matrix for animal choice and oucome vector and fit a
%logistic regression model

%----------version for constant versus shifting bias--------------------
choicePredictors = NaN(SessionData.nTrials-1,2);
%The model takes the form:
%choice ~ constantSideBias + a1*Stimulus + a2*ChoosingSamePortAgain
for k = 2:SessionData.nTrials %here only
    if ~isnan(correctChoice(k)) & ~isnan(correctChoice(k-1))
        %First column is the stimulus. Here, we look at the case for
        %the easiest condition with only two stimuli, thus we can set
        %the left indicating to 0 and right indicating one to 1
        choicePredictors(k-1,1) = SessionData.CorrectSide(k)-1;
        if ~isempty(strfind(experimentName{s},'holdCenterStimUntilCorrect')) || (isfield(SessionData,'reviseChoiceFlag') && SessionData.reviseChoiceFlag)
            if correctChoice(k-1) == 0 & SessionData.Rewarded(k-1) == 1 %If the mouse revises then it will be at the other port to start a new trial
                if responseSide(k-1) == 2
                    choicePredictors(k-1,2) = 0; %If animal responded incorrectly to the right but harvested on the left count as response on the left
                else
                    choicePredictors(k-1,2) = 1; %Vice versa if incorrect on left
                end
            elseif correctChoice(k-1) == 0 & SessionData.Rewarded(k-1) == 1
                choicePredictors(k-1,2) = responseSide(k-1)-1;
            end
        else
            choicePredictors(k-1,2) = responseSide(k-1)-1;
        end
    else
        choicePredictors(k-1,1:2) = NaN;
    end
end
%-------

sidesChosen = (responseSide(2:end) - 1)';

if newChipmunk
    sidesChosen = sidesChosen+1;
    choicePredictors = choicePredictors+1;
end

lastwarn(''); %clear the last warning
[choiceCoef(:,s), dev, choiceModelStats{s}] = glmfit(choicePredictors,sidesChosen,'binomial','link','logit');
modelFitWarnings{s} = lastwarn;
if ~isempty(lastwarn)
    choiceCoef(:,s) = NaN(3,1);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%--------------------------------------------------------------------------
%Trial number per time histogram
trialEnds = SessionData.TrialEndTimestamp/60; % in minutes
binEdges = trialTimeBinSize/2:trialTimeBinSize:trialEnds(end);
trialDist = histcounts(trialEnds,length(binEdges));

%--------------------------------------------------------------------------
%Inter-trial intervals, defined as the time between end of last trial and
%initiation of the next trial by poking in the center
interTrialInterval = [];
if newChipmunk
    lastStartTime = SessionData.TrialStartTimestamp(1); %Set to first trial
    for k=1:SessionData.nTrials
        if SessionData.OutcomeRecord(k) > -2 %Consider initiated trials
            interTrialInterval(end+1) = (SessionData.TrialStartTimestamp(k) + SessionData.RawEvents.Trial{1,k}.States.DemonInitFixation(1)) - lastStartTime;
            lastStartTime = SessionData.TrialStartTimestamp(k) + SessionData.RawEvents.Trial{1,k}.States.DemonInitFixation(1);
        end
    end
    
else
for k=1:SessionData.nTrials    
interTrialInterval(k) = SessionData.RawEvents.Trial{1,k}.States.GoToCenter(2);
end
end

itiDistrObj = fitdist(interTrialInterval','kernel'); 
itiXdata = (0:1:300)';
itiYdata = pdf(itiDistrObj,itiXdata); 

%%%----Get initiation delay for the demonstrator
initiationDelay = NaN(SessionData.nTrials,1);
for k=1:SessionData.nTrials
    if isfield(SessionData, 'PacedFlag')
        if SessionData.PacedFlag
            initiationDelay(k) = SessionData.RawEvents.Trial{1,k}.States.DemonTrialStart(2) - SessionData.RawEvents.Trial{1,k}.States.DemonTrialStart(1);
        else
            try
                initiationDelay(k) = SessionData.RawEvents.Trial{1,k}.States.GoToCenter(2) - SessionData.RawEvents.Trial{1,k}.States.GoToCenter(1);
            catch
                initiationDelay(k) = SessionData.RawEvents.Trial{1,k}.States.DemonWaitForCenterFixation(2) - SessionData.RawEvents.Trial{1,k}.States.DemonWaitForCenterFixation(1);
            end
        end
    else
        initiationDelay(k) = SessionData.RawEvents.Trial{1,k}.States.GoToCenter(2) - SessionData.RawEvents.Trial{1,k}.States.GoToCenter(1);
    end
end

%Generate histogram of initiation delay
initDelayBinEdges = 0:0.2:10;
initiationDelayCounts = histcounts(initiationDelay,initDelayBinEdges);
initiationDelayBinCenters = initDelayBinEdges(1:end-1) + 0.1;

%--------------------------------------------------------------------------
% Compute within session graphs
%performance across session
 window = 24; trialVect = []; withinPerformance = []; withinEarlyWithdrawals = [];
 for k=window/2:1:SessionData.nTrials-(window/2)
     withinPerformance(:,end+1) = correctChoice(k-window/2+1:k+window/2)';
     withinEarlyWithdrawals(end+1) = nanmean(earlyWithdrawals{end}(k-window/2+1:k+window/2));
     trialVect(end+1) = k;
 end
%--------------------------------------------------------------------------
% Generate the psychometric plot
% trialStimFreq = [];
% for n = 1:SessionData.nTrials; % extract the list of unique stimulus frequencies
%     trialStimFreq(end+1) = length(SessionData.stimEventList{n})-1; %Disregard modality for now
% end
% uniqueFreqVals = unique(trialStimFreq); %returns values in ascending order

%uniqueFreqVals = unique(stimFreq{end}); %returns values in ascending order
uniqueFreqVals = unique(stimFreq{end}(:)); %Get the unique frequencies for all the modalities combined
uniqueFreqVals(isnan(uniqueFreqVals)) = []; % Remove NaN in the data
rightSideChoices = cell(3,length(uniqueFreqVals)); %Find the choices sides for all the frequencies presented
for k = 1:SessionData.nTrials
    if newChipmunk
        %modality_idx = find(~isnan(stimFreq{end}(k,:)));
%         modality_idx = tmp(1);
        [currentFreq,inCols,~]= find(uniqueFreqVals == stimFreq{end}(k,:)); %Get the current frequency and localize the presented modalities
       % rightSideChoices{find(uniqueFreqVals == stimFreq{end}(k,modality_idx))}(end+1) = responseSide(k);
         rightSideChoices{inCols, currentFreq(1)}(end+1) = responseSide(k); %Assign to the proper modality
    else
        rightSideChoices{find(uniqueFreqVals == stimFreq{end}(k))}(end+1) = responseSide(k)-1;
    end
end

avRightSideChoice = nan(size(rightSideChoices)); lowerBound = nan(size(rightSideChoices)); upperBound = nan(size(rightSideChoices));
for k=1:length(uniqueFreqVals)
    for n=1:3 %The sensory modalities
        if ~isempty(rightSideChoices{n,k})
            avRightSideChoice(n,k) = nanmean(rightSideChoices{n,k});
            [lowerBound(n,k), upperBound(n,k)] = calculateWilsonScoreInerval(rightSideChoices{n,k});
        end
    end
end

%Get the average performance on the easiest condtions for all the sessions
easiestStimuli = [min(uniqueFreqVals) max(uniqueFreqVals)]; %Find the easiest conditions assuming that there will always be presentation of these
for k=1:sessionsBack+1
    tmp = nanmean([performance{k}(find(stimFreq{k}(:,1) == easiestStimuli(1))),performance{k}(find(stimFreq{k}(:,1) == easiestStimuli(2)))]);
    tmp = [tmp; nanmean([performance{k}(find(stimFreq{k}(:,2) == easiestStimuli(1))),performance{k}(find(stimFreq{k}(:,2) == easiestStimuli(2)))])];
    tmp = [tmp; nanmean([performance{k}(find(stimFreq{k}(:,3) == easiestStimuli(1))),performance{k}(find(stimFreq{k}(:,3) == easiestStimuli(2)))])];
    
    performanceEasiest(k) = nanmean(tmp);
    
%     performanceEasiest(k) = nanmean([performance{k}(find(stimFreq{k} == easiestStimuli(1))),...
%     performance{k}(find(stimFreq{k} == easiestStimuli(2)))]);
end

%Calculate d' value for easiest
correctSideEasiest = SessionData.CorrectSide(sum(stimFreq{k} == easiestStimuli(1),2) | sum(stimFreq{k} == easiestStimuli(2),2));
responseSideEasiest = responseSide(sum(stimFreq{k} == easiestStimuli(1),2) | sum(stimFreq{k} == easiestStimuli(2),2)); %Get the side picked by the animal only on either the easy high or low rate

dPrime = dPrime2AFC(correctSideEasiest, responseSideEasiest); %Compute dPrime, no 1/2^(1/2) correction here!

%Fit psychometric curve if more than two frequencies are present
fitPsychometric = length(uniqueFreqVals) > 2;
if fitPsychometric
    % [PMparameters, PFxdata, PFydata] = fitPsychometricFunction(uniqueFreqVals, rightSideChoices, SessionData.TrialSettings(1).highRateSide);
    perceptualModelStats = cell(1,3); PFxdata = []; PFydata = []; freqBias = []; freqSensitivity = [];
    if newChipmunk
        for n=1:3
            if ~isempty(rightSideChoices{n,1})
                tmpResp = []; tmpStim = [];
                for q=1:size(rightSideChoices,2)
                    tmpResp = [tmpResp rightSideChoices{n,q}]; %Rebuild a continuous representation, arrgh!
                    tmpStim = [tmpStim ones(1,size(rightSideChoices{n,q},2))*uniqueFreqVals(q)]; %Match the stim strength
                end
                 [perceptualModelStats{n}, PFxdata(n,:), PFydata(n,:), freqBias(n), freqSensitivity(n)] = logistic_multiFreq(tmpResp, tmpResp, tmpStim, 12);%The first tmpResp is provided instead of performance. Just to check for nan
            end
        end
      else
        [perceptualModelStats, PFxdata, PFydata, freqBias, freqSensitivity] = logistic_multiFreq(performance{end}, responseSide-1, stimFreq{end}, 12);
    end
else
    freqBias = [];
    freqSensitivity = [];
end
%------
%Wait time difference and reaction time current session
waitTimeDiff = absoluteWaitTime{end} - targetWaitTime{end};
goCueReactionTimes = absoluteWaitTime{end}(earlyWithdrawals{s}==0) - targetWaitTime{s}(earlyWithdrawals{end}==0);

watiTimeDistrObj = fitdist(absoluteWaitTime{end}','kernel'); 
waitTimeXdata = (0:0.01:3)';
waitTimeYdata = pdf(watiTimeDistrObj,waitTimeXdata); 

%--------
%Motion time distribution
motionTimeDistrObj = fitdist(motionTime{end}','kernel'); 
motionTimeXdata = (0:0.01:3)';
motionTimeYdata = pdf(motionTimeDistrObj,motionTimeXdata); 

%-----
%Over multiple sessions
earlyWithdrawalSessions = []; performanceSessions = []; movementSessions = []; gapTimeSessions = []; 
for s = 1:sessionsBack+1
    if newChipmunk
        earlyWithdrawalSessions(s) = sum(earlyWithdrawals{s})/numTrials(s);
    else
    earlyWithdrawalSessions(s) = nanmean(earlyWithdrawals{s});
    end
    performanceSessions(s) = nanmean(performance{s});
    movementSessions(s) = nanmedian(motionTime{s});
    gapTimeSessions(s) = nanmedian(gapTime{s});
end

%-----------------------------------
%%%%%%
%% Plotting business
%Before plotting
if strfind(experimentName{s},'Task')
    performancePlotsColor = [0.9 0.4 0.1];
elseif ~isempty(strfind(experimentName{s},'holdCenterStimUntilCorrect')) || ~isempty(strfind(experimentName{s},'Fixation'))
    performancePlotsColor = [0.9 0.7 0.1];
end
psychometricPlotsColor = [0.02 0.26 0.54]; %plot visual stimuli on dark blue
psychometricPlotsColor(2,:) = [0.02 0.56 0.4]; %plot auditory stimuli in dark green
psychometricPlotsColor(3,:) = [0.3 0.3 0.3]; %plot multi-sensory stimuli in dark gray

%--plotting
fi = figure('Position',[49, 42, 0.9*1280, 0.9*720],'NumberTitle','off','MenuBar','none');
%
%Start the panel structure
figPanels = struct();
   figPanels.parameters = uipanel('Parent', fi, 'Units', 'Normal',...
            'Position', [0.05, 0.75, 0.3, 0.25],'Title','Parameters','FontWeight','bold','BackgroundColor','w');
        figPanels.acrossTrials = uipanel('Parent', fi, 'Units', 'Normal',...
            'Position', [0.05, 0, 0.3, 0.75],'Title','Across trials','FontWeight','bold','BackgroundColor','w');
         figPanels.withinSession = uipanel('Parent', fi, 'Units', 'Normal',...
            'Position', [0.35, 0, 0.3, 1],'Title','Within session','FontWeight','bold','BackgroundColor','w');
         figPanels.acrossSessions = uipanel('Parent', fi, 'Units', 'Normal',...
            'Position', [0.65, 0, 0.3, 1],'Title','Across sessions','FontWeight','bold','BackgroundColor','w');
        
%---
%The parameters box
       
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,10/11,0.35,1/11],'style', 'text', 'String','Animal ID:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,9/11,0.35,1/11],'style', 'text', 'String','Experiment phase:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,8/11,0.35,1/11],'style', 'text', 'String','Body weight:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,7/11,0.35,1/11],'style', 'text', 'String','Number of trials:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,6/11,0.35,1/11],'style', 'text', 'String','Early withdrawal rate:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,5/11,0.35,1/11],'style', 'text', 'String','No choice made:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,4/11,0.35,1/11],'style', 'text', 'String','Median trial duration (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,3/11,0.35,1/11],'style', 'text', 'String','Reward size (ml):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,2/11,0.35,1/11],'style', 'text', 'String','Earned on left (ml):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,1/11,0.35,1/11],'style', 'text', 'String','Earned on right (ml):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0,0,0.35,1/11],'style', 'text', 'String','Total reward (ml):', 'HorizontalAlignment','left','BackgroundColor','w');

 %uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,10/11,0.12,1/11],'style', 'text', 'String',SessionData.TrialSettings(end).SubjectName, 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,9/11,0.12,1/11],'style', 'text', 'String',experimentName{end}, 'HorizontalAlignment','left','BackgroundColor','w');
 %uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,8/11,0.12,1/11],'style', 'text', 'String',num2str(bodyWeight(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,7/11,0.12,1/11],'style', 'text', 'String',num2str(numTrials(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,6/11,0.12,1/11],'style', 'text', 'String',num2str(earlyWithdrawalSessions(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,5/11,0.12,1/11],'style', 'text', 'String',num2str(noChoice(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,4/11,0.12,1/11],'style', 'text', 'String',num2str(nanmedian(trialDuration{end})), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,3/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).leftRewardVolume), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,2/11,0.12,1/11],'style', 'text', 'String',num2str(totalLeftRewards/1000), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,1/11,0.12,1/11],'style', 'text', 'String',num2str(totalRightRewards/1000), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,0,0.12,1/11],'style', 'text', 'String',num2str(totalReward(end)/1000), 'HorizontalAlignment','left','BackgroundColor','w');
 if newChipmunk
     uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,10/11,0.12,1/11],'style', 'text', 'String',SessionData.TrialSettings(end).demonID, 'HorizontalAlignment','left','BackgroundColor','w');
     uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,8/11,0.12,1/11],'style', 'text', 'String',num2str(bodyWeight(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 else
     uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,10/11,0.12,1/11],'style', 'text', 'String',SessionData.TrialSettings(end).SubjectName, 'HorizontalAlignment','left','BackgroundColor','w');
     uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.36,8/11,0.12,1/11],'style', 'text', 'String',num2str(bodyWeight(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 end

 
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,10/11,0.35,1/11],'style', 'text', 'String','Early timeout (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,9/11,0.35,1/11],'style', 'text', 'String','Wrong timeout (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,8/11,0.35,1/11],'style', 'text', 'String','Fraction left side trials:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,7/11,0.35,1/11],'style', 'text', 'String','Mean extended stim (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,6/11,0.35,1/11],'style', 'text', 'String','Median wait time (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,5/11,0.35,1/11],'style', 'text', 'String','Median gap time (s):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,4/11,0.35,1/11],'style', 'text', 'String','Left side preference:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,3/11,0.35,1/11],'style', 'text', 'String','Performance on easiest:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,2/11,0.35,1/11],'style', 'text', 'String','d'' on easiest:', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,1/11,0.35,1/11],'style', 'text', 'String','Perceptual bias (Hz):', 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.52,0/11,0.35,1/11],'style', 'text', 'String','Perceptual sensitivity:', 'HorizontalAlignment','left','BackgroundColor','w');


 %uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,10/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).EarlyPunishDuration), 'HorizontalAlignment','left','BackgroundColor','w');
 %uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,9/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).PunishDuration), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,8/11,0.12,1/11],'style', 'text', 'String',num2str(fractionLeftAssigned(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,7/11,0.12,1/11],'style', 'text', 'String',num2str(extendedStim(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,6/11,0.12,1/11],'style', 'text', 'String',num2str(nanmedian(absoluteWaitTime{end})), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,5/11,0.12,1/11],'style', 'text', 'String',num2str(gapTimeSessions(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,4/11,0.12,1/11],'style', 'text', 'String',num2str(leftwardBias(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,3/11,0.12,1/11],'style', 'text', 'String',num2str(performanceEasiest(end)), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,2/11,0.12,1/11],'style', 'text', 'String',num2str(dPrime), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,1/11,0.12,1/11],'style', 'text', 'String',num2str(freqBias), 'HorizontalAlignment','left','BackgroundColor','w');
 uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,0/11,0.12,1/11],'style', 'text', 'String',num2str(freqSensitivity), 'HorizontalAlignment','left','BackgroundColor','w');
if newChipmunk 
  uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,10/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).earlyPunishTimeout), 'HorizontalAlignment','left','BackgroundColor','w');
  uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,9/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).wrongPunishTimeout), 'HorizontalAlignment','left','BackgroundColor','w');
else
  uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,10/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).EarlyPunishDuration), 'HorizontalAlignment','left','BackgroundColor','w');
  uicontrol('Parent', figPanels.parameters,'Units', 'normal', 'Position',[0.88,9/11,0.12,1/11],'style', 'text', 'String',num2str(SessionData.TrialSettings(end).PunishDuration), 'HorizontalAlignment','left','BackgroundColor','w');
end

%--------------------------------------------------------------------------
%The read-outs across time
acrossTrials = struct();

acrossTrials.performance = axes('Parent', figPanels.acrossTrials,...
            'Units', 'Normal', 'Position', [0.1, 0.74, 0.8, 0.2],'tickdir','out');
 
acrossTrials.waitTimes = axes('Parent', figPanels.acrossTrials,...
            'Units', 'Normal', 'Position', [0.1, 0.42, 0.8, 0.2],'tickdir','out');      

acrossTrials.trialHist = axes('Parent', figPanels.acrossTrials,...
            'Units', 'Normal', 'Position', [0.1, 0.1, 0.8, 0.2],'tickdir','out');
        
%now the plotting
%performance
axes(acrossTrials.performance)
errorPlottingOptions = [];
errorPlottingOptions.handle = acrossTrials.performance; %setting the plotting parameters
errorPlottingOptions.error = 'sem';
errorPlottingOptions.color_area = performancePlotsColor;
errorPlottingOptions.alpha = 0.3;
errorPlottingOptions.color_line = performancePlotsColor;
errorPlottingOptions.line_width = 1.3;
errorPlottingOptions.x_axis = trialVect;
plot_areaerrorbar(withinPerformance, errorPlottingOptions); %commented figure command in this function (line 59)
xlabel('Trial number')
ylabel(sprintf('Fraction correct / %d trials',round(window)))
title('Performance over trials')
ylim([0 1]); xlim([0 SessionData.nTrials])
set(gca,'TickDir','out')
line(xlim,[0.5 0.5], 'color', 'k', 'LineStyle','--')
box off
%
hold on
yyaxis right
plot(trialVect, withinEarlyWithdrawals,'LineWidth',1.3,'Color',[0 0.4470 0.7410]);
ylim([0 1])
ylabel(sprintf('Early withdrawal rate / %d trials',round(window)))
set(gca,'YColor',[0 0.4470 0.7410], 'TickDir','out');
box off

%-----
% Wait time difference
axes(acrossTrials.waitTimes)
hold on
dotSize  = 36 * 100/SessionData.nTrials;
sucH = scatter(find(earlyWithdrawals{s} == 0), waitTimeDiff(find(earlyWithdrawals{end} == 0)),dotSize, 'filled','MarkerFaceColor',[0.1 0.4 1]);
if isfield(SessionData, 'ObsOutcomeRecord')
    sucM = scatter(find(SessionData.ObsOutcomeRecord==1 & SessionData.OutcomeRecord > -1), waitTimeDiff(find(SessionData.ObsOutcomeRecord==1 & SessionData.OutcomeRecord > -1)), dotSize,'filled','MarkerFaceColor',[0, 0.8, 0.4]);
end
eWitH = scatter(find(earlyWithdrawals{s} == 1), waitTimeDiff(find(earlyWithdrawals{end} == 1)),dotSize,'MarkerEdgeColor','r');
%include info about the wait time progresion, very clunky...
if ischar(SessionData.TrialSettings(2).minWaitTime) && ischar(SessionData.TrialSettings(end).minWaitTime)
    report_string = 'Wait time start: %s, end: %s s';
elseif ischar(SessionData.TrialSettings(2).minWaitTime) && isnumeric(SessionData.TrialSettings(end).minWaitTime)
     report_string = 'Wait time start: %s, end: %0.3f s'; %0.3f specifies 3 display digits after the zero with f focing decimal display
elseif isnumeric(SessionData.TrialSettings(2).minWaitTime) && ischar(SessionData.TrialSettings(end).minWaitTime)
    report_string = 'Wait time start: %0.3f, end: %s s';
elseif isnumeric(SessionData.TrialSettings(2).minWaitTime) && isnumeric(SessionData.TrialSettings(end).minWaitTime)
     report_string = 'Wait time start: %0.3f, end: %0.3f s';
end
time_disp = plot(NaN, 'w');
xlim([0 length(waitTimeDiff)+1]); ylim([-1.5 1]);
line(xlim,[0 0],'LineStyle','--','Color','k')
set(gca,'tickDir','out'); box off;
xlabel('Trial number')
ylabel('Time difference (s)')
title('Difference between target and recorded wait time')
if isfield(SessionData, 'ObsOutcomeRecord')
    lg = legend([sucM, sucH, eWitH, time_disp],'Mutual success','Demon success','Demon early withdrawal',...
    sprintf(report_string, SessionData.TrialSettings(2).minWaitTime, SessionData.TrialSettings(end).minWaitTime),'Location','Best');
else
    lg = legend([sucH, eWitH, time_disp],'Success','Early withdrawal',...
        sprintf(report_string, SessionData.TrialSettings(2).minWaitTime, SessionData.TrialSettings(end).minWaitTime),'Location','Best');
end
lg.BoxFace.ColorType='truecoloralpha'; %Set semi-transparent background for legend
lg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
%-----------------
% The trial count histogram now
axes(acrossTrials.trialHist)
bar(binEdges, trialDist,'FaceColor',[0.4 0.4 0.4],'FaceAlpha',0.5)
xlabel('Time within session')
ylabel('Trial count')
title(sprintf('Trial counts per %d min', round(trialTimeBinSize)))
box off; set(gca,'tickDir','out')


%-----------------------------
%The read-outs from the entire session
withinSession = struct();

withinSession.psychometric = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.1, 0.8, 0.45, 0.16],'tickdir','out');

withinSession.modelBiases = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.65, 0.8, 0.25, 0.16],'tickdir','out');

withinSession.waitTime = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.1, 0.55, 0.45, 0.16],'tickdir','out');

withinSession.reactionTime = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.65, 0.55, 0.25, 0.16],'tickdir','out');

withinSession.motionTime = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.1, 0.3, 0.45, 0.16],'tickdir','out');

withinSession.gapTime = axes('Parent', figPanels.withinSession,...
    'Units','Normal','Position',[0.65, 0.3, 0.25, 0.16],'tickdir','out');

withinSession.interTrialInterval = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.1, 0.05, 0.45, 0.16],'tickdir','out');

withinSession.initiationDelay = axes('Parent', figPanels.withinSession,...
    'Units', 'Normal', 'Position', [0.65, 0.05, 0.25, 0.16],'tickdir','out');

%-------------
axes(withinSession.psychometric)
for n=1:3
    if ~isnan(avRightSideChoice(n,1))
        errorbar(uniqueFreqVals',avRightSideChoice(n,:),avRightSideChoice(n,:) - lowerBound(n,:),upperBound(n,:) -avRightSideChoice(n,:),'o','Color',psychometricPlotsColor(n,:),'MarkerSize',3,'MarkerFaceColor',psychometricPlotsColor(n,:))
        hold on
        if fitPsychometric
            plot(PFxdata(n,:),PFydata(n,:),'Color',psychometricPlotsColor(n,:))
        end
    end
end
    %plot(uniqueFreqVals,avRightSideChoice,'-o','Color',performancePlotsColor)
box off; set(gca,'TickDir','out');
ylim([0 1]); yticks(0:0.2:1); xlim([0 24]); xticks(4:4:20) 
grid on
xlabel('Stimulation frequency (Hz)')
ylabel('Fraction of right side choices')
title('Psychometric plot')

axes(withinSession.modelBiases)
scatter(choiceCoef(2,end), choiceCoef(3,end),'filled','MarkerFaceColor',performancePlotsColor)
box off; grid on;
xlim([-5 10]); ylim([-10 10]);
title('Stim & choice history')
cuAx = gca;
cuAx.XAxisLocation = 'origin'; cuAx.YAxisLocation = 'origin';
Ylm=ylim;                          % get x, y axis limits 
Xlm=xlim;                          % so can position relative instead of absolute
Xlb=mean(Xlm);                    % set horizontally at midpoint
Ylb=0.99*Ylm(1);                  % and just 1% below minimum y value
xlabel('Stimulus weight','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 
Xlb = 1.01*Xlm(1);
Ylb = mean(Ylm);
ylabel('Switch <> Stay','Position',[Xlb,Ylb],'VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90)

axes(withinSession.reactionTime)
violinplot(goCueReactionTimes); %represent the reaction time distribution as a violin plot
set(gca,'xTick',[]); set(gca,'TickDir','out'); ylim([0 1]);
yticks(0:0.2:1)
ylabel('Reaction time (s)')
title('Go cue reaction times')

axes(withinSession.waitTime)
plot(waitTimeXdata,waitTimeYdata)
hold on
line([nanmedian(targetWaitTime{end}) nanmedian(targetWaitTime{end})],ylim,'LineStyle','--','Color','k')
xlim([0 3]);% ylim([0 sum(waitTimeYdata)/length(waitTimeYdata)]) %Y limit to not exceed max if every observation was inside one xvalue bin
set(gca,'tickDir','out'); box off
xlabel('Wait time (s)')
ylabel('Density')
title('Wait time distribution')
lg = legend({'Observed','Median target'},'Location','Best');
lg.BoxFace.ColorType='truecoloralpha'; %Set semi-transparent background for legend
lg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

axes(withinSession.motionTime)
plot(motionTimeXdata,motionTimeYdata,'color',[0.1 0.76 0.3])
xlim([0 3]);%
set(gca,'tickDir','out'); box off
xlabel('Motion time (s)')
ylabel('Density')
title('Time to move to side port')

axes(withinSession.gapTime)
hold on
spacingSuccess = (0.5 - rand(length(find(performance{end} == 1)),1))*0.6;
spacingFail = (0.5 - rand(length(find(performance{end} == 0)),1))*0.6;
scatter(1+spacingSuccess, gapTime{end}(find(performance{end} == 1)),15,'filled','MarkerFaceColor',[0.35 0.8 0.1])
scatter(2+spacingFail,gapTime{end}(find(performance{end} == 0)),15,'filled','MarkerFaceColor',[0.79 0.2 0.1])
%line([1-0.3 1+0.3],[mean(reportingSpeed(find(OutcomeRecord == 1))) mean(reportingSpeed(find(OutcomeRecord == 1)))],'color','k','LineWidth',1.6)
%line([2-0.3 2+0.3],[mean(reportingSpeed(find(OutcomeRecord == 0))) mean(reportingSpeed(find(OutcomeRecord == 0)))],'color','k','LineWidth',1.6)
set(gca,'tickDir','out'); box off; grid on;
xlim([0 3]); xticks([1 2]); xticklabels({'Correct', 'Incorrect'}); xtickangle(45)
ylim([0 1.5])
ylabel('Gap time (s)')
title('Gap time')

axes(withinSession.interTrialInterval)
plot(itiXdata,itiYdata,'color',[0.4 0.4 0.4])
xlim([0 30]);%
set(gca,'tickDir','out'); box off
xlabel('Inter-trial interval (s)')
ylabel('Density')
title('Inter-trial interval')

axes(withinSession.initiationDelay)
bar(initiationDelayBinCenters,initiationDelayCounts,'FaceColor',[0.6 0.6 0.6],'LineWidth',0.05)
set(gca,'tickDir','out'); box off;
ylim([0 100]) %Possibly many trials are not initiated...
xlabel('Initiation delay (s)')
ylabel('Counts')
title('Initiation delays')

%----------------------------------
%Get the axes objects for the across sessions panel
acrossSessions = struct();

acrossSessions.performance = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.9, 0.8, 0.08],'tickdir','out');

acrossSessions.trialHistoryWeights = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.78, 0.8, 0.08],'tickdir','out');

acrossSessions.earlyWithdrawalRate = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.66, 0.8, 0.08],'tickdir','out');

acrossSessions.motionTime = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.54, 0.8, 0.08],'tickdir','out');

acrossSessions.gapTime = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.42, 0.8, 0.08],'tickdir','out');

acrossSessions.trialNum = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.3, 0.8, 0.08],'tickdir','out');

acrossSessions.waterRewards = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.18, 0.8, 0.08],'tickdir','out');
title('Water rewards earned')

acrossSessions.bodyWeight = axes('Parent', figPanels.acrossSessions,...
    'Units', 'Normal', 'Position', [0.1, 0.06, 0.8, 0.08],'tickdir','out');

%---------------
axes(acrossSessions.performance)
plot(performanceEasiest,'-o','color',performancePlotsColor)
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([0:sessionsBack+1]); xticklabels([]);
yticks([0.5 0.75 1])
ylim([0.4 1])
box off; grid on
set(gca,'TickDir','out')
ylabel('Fract. correct')
title('Performance on easiest')

axes(acrossSessions.trialHistoryWeights)
plot(choiceCoef(2,:),'-o','color',performancePlotsColor)
hold on
plot((1:sessionsBack+1),choiceCoef(3,:),'-o','color','k')
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
ylim([-10 10]);
box off; grid on
set(gca,'TickDir','out')
ylabel('log(odds)')
title('Choice history weights')
%legend({'Stimulus','Prior choice'},'Location','NorthWest')

axes(acrossSessions.earlyWithdrawalRate)
plot(earlyWithdrawalSessions,'-o')
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([1:sessionsBack+1]); xticklabels([]);
box off; grid on
ylim([0 1])
set(gca,'TickDir','out')
ylabel('Fraction')
title('Early withdrawal rate over sessions')

axes(acrossSessions.motionTime)
plot(movementSessions,'-o','color',[0.1 0.76 0.3])
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([1:sessionsBack+1]); xticklabels([]);
box off; grid on
ylim([0 3])
set(gca,'TickDir','out')
ylabel('Time (s)')
title('Median motion time over session')

axes(acrossSessions.gapTime)
plot(gapTimeSessions,'-o','color',[0.1 0.76 0.3])
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([1:sessionsBack+1]); xticklabels([]);
box off; grid on
ylim([0 1.5])
set(gca,'TickDir','out')
ylabel('Time (s)')
title('Median gap between stim offset and reporting')

axes(acrossSessions.trialNum)
plot(numTrials,'-o','color',[0.4 0.4 0.4])
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([1:sessionsBack+1]); xticklabels([]);
box off; grid on
set(gca,'TickDir','out')
ylabel('#')
title('Trial number over sessions')

axes(acrossSessions.waterRewards)
plot(totalReward/1000,'-o','color',[0.4 0.4 0.4])
xlim([0 sessionsBack+1])
xticks([0 sessionsBack+1]); xticklabels([]);
%xticks([1:sessionsBack+1]); xticklabels([]);
ylim([0 1.6])
box off; grid on
set(gca,'TickDir','out')
ylabel('Volume (ml)')
title('Water reward earned over sessions')

axes(acrossSessions.bodyWeight)
plot(bodyWeight,'-o','color',[0.4 0.4 0.4])
xlim([0 sessionsBack+1])
% xticks([0 sessionsBack+1]); xticklabels([]);
xticks([0:sessionsBack+1]);  xticklabels([{''},dates]); xtickangle(45);
ylim([12 30])
box off; grid on
set(gca,'TickDir','out')
ylabel('Weight (g)')
title('Body weight over sessions')

saveas(fi,fullfile(outputLoc,[fileName(1:end-4) '_summary.png']))
figureHandle = fi;

end





