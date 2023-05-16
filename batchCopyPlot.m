function batchCopyPlot(animals, inputLoc, labdataLoc, outputLoc, sessionsBack, datestring);
% batchCopyPlot(animals, inputLoc, copyTo, labdataLoc, sessionsBack, datestring);
%
% Function to plot the performance metrics in the chipmunk task for a set
% of animals. The function copies data from a remote location to a local
% directory and runs the performanceOverviewDecisionTask function.
% Please note the following: This function is built to maintain the
% Churchland lab data structure conventions in the file destination folder.
% Currently, this function doesn't handle instances with multiple sessions
% per mouse very well...
% Also, this function doesn't automatically down-load any other sessions
% other than the ones recorded on the specified date. Make sure to have the
% other sessions local already if you use sessionsBack > 0.
%
% Inputs:
% animals (mandatory): cell array of animal identifiers.
% inputLoc (mandatory): string, path to where the data are stored
% labdataLoc (mandatory): string, path of the local labdata folder
% outputLoc (mandatory): string, path to save figures to
% sessionsBack (mandatory): int, look at the last x sessions
% datestring (optional): string, specify a custom date to plot the data for
%
% LO, 05/16/2023
%--------------------------------------------------------------------------
%Input checks and assignment of defaults
if ~exist('sessionsBack') || isempty(sessionsBack)
    sessionsBack = 10;
end

if ~exist('datestring') || isempty(datestring)
    tmp = datetime("now"); %Use the current date as the reference if nothing is specified
    datestring = datestr(tmp, 'yyyymmdd'); 
end


%Start looping through the specified animals
for k =1:length(animals)
    
    %First find the directory of the specified date
    folderPath = fullfile(inputLoc, animals{k});
    dirList = struct2cell(dir(folderPath)); %Where rows are: name, folder, date, bytes, isdir and datenum
    directories = dirList(1, cell2mat(dirList(5,:)) == 1); %Find the entries that are directories and retrieve folder names for these only
    for n=1:length(directories)
        if ~isempty(strfind(directories{n}, datestring))
            matched = n;
        end
    end
    
    %Find the only .mat file in the corresponding session and chipmunk
    %folder
    cFileStruct = dir(fullfile(folderPath, directories{matched},'chipmunk', '*.mat'));
    
    %Create the analog directory structure locally, pretty clunky...
    mkdir(labdataLoc, animals{k});
    mkdir(fullfile(labdataLoc, animals{k}), directories{matched});
    mkdir(fullfile(labdataLoc, animals{k}, directories{matched}),'chipmunk');
    copyPath = fullfile(labdataLoc, animals{k},directories{matched},'chipmunk');
    
    %Now the copying
    [copyStatus, copyErrorMsg] = copyfile(fullfile(cFileStruct.folder, cFileStruct.name),copyPath, 'f');
    
    %Define the session file
    sessionFile = fullfile(copyPath, cFileStruct.name);
    
    %Now call the performance plot function
    performanceOverviewDecisionTask(sessionFile,outputLoc, sessionsBack);
    %No returns needed here
    
    %And that's it already...
end
    