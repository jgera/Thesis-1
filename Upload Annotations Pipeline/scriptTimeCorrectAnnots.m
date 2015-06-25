%% Script to clip interictal, preictal, ictal segments from Neocortical Epilepsy Mayo Datasets from IEEG Portal

%% Open session and queue study snapshots
clear; clc;
snapList = {'Study 004-2', ...
    'Study 006', ...
    'Study 010', ...
    'Study 011', ...
    'Study 012-2', ...
    'Study 014', ...
    'Study 016', ...
    'Study 019', ...
    'Study 020', ...
    'Study 023', ...
    'Study 026', ...
    'Study 028', ...
    'Study 030', ...
    'Study 031', ...
    'Study 033', ...
    'Study 037'};

for i = 1:length(snapList)
    try 
        if i == 1
            session = IEEGSession(snapList{i}, 'tblevins', 'tbl_ieeglogin');
        else
            session.openDataSet(snapList{i});
        end
        disp(sprintf('Loading %s', snapList{i}));    
    catch
        disp(sprintf('Error loading %s', snapList{i}));
    end
end

%% Run through each snapshot and parse annotation events
pad_secs = 300;

for ptID = 1:length(session.data)
    lclStudyName = sprintf('MAYO_NEO_%0.3d', ptID);
    disp(sprintf('Local Study: %s', lclStudyName));
    
    % Dataset info
    chanNum = length(session.data(ptID).channels);
    
    % Standard datafile contents
    Fs = round(session.data(ptID).sampleRate);
    channels = cell(chanNum,1);
    for c = 1:chanNum
        channels{c} = session.data(ptID).channels(c).label;
    end
    
    % Process Annotations
    eventNum = 1;
    try 
        eventList = session.data(ptID).annLayer.getEvents(0);
    catch
        disp(sprintf('Error getting annotations for %s', snapList{ptID}));
    end    
    for eventID = 1:length(eventList)
        if strcmp(eventList(eventID).type, 'Seizure')
            lclStudyFile = [lclStudyName, sprintf('-Sz%0.3d', eventNum)];
            disp(sprintf('   Clip: %s', lclStudyFile));
            
            % Clip sample period           
            szTimeStart = eventList(eventID).start;
            clipTimeLength = 2*pad_secs*1e6;
            clipTimeStart = szTimeStart - pad_secs*1e6;
            clipTimeStop = clipTimeStart + clipTimeLength - 1;
            
            clipIdx = Fs*clipTimeStart/1e6 : Fs*clipTimeStop/1e6;
            clipSampLength = length(clipIdx);
            
            % Allocate memory and retrieve ieeg data channel-by-channel

                       
            % Placeholder for marking seizure type
            szType = 'NA'; %'2ndGen', 'CPS', 'SPS'
            
            % Save clipped data
            
            evTimeStart = clipTimeStart;
            evTimeStop = clipTimeStop;
            save([lclStudyFile], 'Fs', 'evTimeStart', 'evTimeStop', 'szType', '-v7.3'); 

            clear evData szClip preszClip;
            
            eventNum = eventNum + 1;
        end
    end
end