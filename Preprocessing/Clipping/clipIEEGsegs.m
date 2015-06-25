%% Script to pull and organize data into 1 second clips for Kaggle algorithms
clear
clc

% Add ieeg toolbox to path
addpath(genpath('ieeg-matlab-1.8.3'))

% Mayo Neo Data Sets
study = {'TB Study 004-2', ...
    'TB Study 006', ...
    'TB Study 010', ...
    'TB Study 011', ...
    'TB Study 012-2', ...
    'TB Study 014', ...
    'TB Study 016', ...
    'TB Study 019'};

% Select study
idx = 7;

% Connect to portal
session = IEEGSession(study{idx},'tblevins','tbl_ieeglogin'); % change to your portal username and pw file

fprintf('\nConnected to the Portal!\n')
% Establish which layer is EEC
layerNames = {session.data.annLayer.name};
ind=find(ismember(layerNames,'KD_BO_EEC'));

% Get EEC start times
EECstart = session.data.annLayer(ind).getEvents(0);
EECtimes = {EECstart.start};

% Print number of seizures in study
fprintf('\n%d seizures detected in data set',length(EECtimes));

% Establish seizure base layer
szind = find(ismember(layerNames,strcat(study{idx}(4:end),' base')));

% Get Seizure end times
szEnd = session.data.annLayer(szind).getEvents(0);
szTimes = {szEnd.stop};

% Only use seizure end times accompanied by an EEC
if idx == 1
    szTimes = {szTimes{1:3}};
elseif idx == 2
    szTimes = {szTimes{[2,3]}};
elseif idx == 3
    szTimes = {szTimes{1:3}};
elseif idx == 4
    szTimes = {szTimes{[1,3]}};
elseif idx == 6
    szTimes = {szTimes{[1,4]}};
elseif idx == 7
    szTimes = {szTimes{[1,2,3,5,6,7]}};
end
    


fprintf('Pulling ictal data...')
% Pull ictal data
ictalSeg = cell(1,length(EECtimes));

for i = 1:length(EECtimes)
    % Pull data between EEC and seizure end
    ictalSeg{i} = session.data.getvalues(EECtimes{i},szTimes{i}-EECtimes{i},1:length(session.data.channels));    
    
    % Remove NaNs
    if sum(sum(isnan(ictalSeg{i}))) > 0 || sum(sum(isinf(ictalSeg{i}))) > 0
        [row,col] = find(isnan(ictalSeg{i}));
        ictalSeg{i}(row,:) = [];
    end
    
end
fprintf('Done!\n')
fprintf('Pulling interictal data...')

% Pull interictal data the same length as ictal data
interSeg = cell(1,length(EECtimes));

% Average ictal time
avgIctalTime = mean(cell2mat(szTimes(1:length(EECtimes)))-cell2mat(EECtimes));

% Sampling rate
fs = round(session.data.sampleRate);


for i = 1:length(EECtimes)
    i
    % training interictal segments
    if i <= round(.4*length(interSeg))
        
        % pull interictal clip 5 times the length of average ictal time
        % interictal clip must not overlap with seizure
        if i==1 || 5*avgIctalTime < EECtimes{i}-szTimes{i-1}
        interSeg{i} = session.data.getvalues(round(EECtimes{i}/2)-2.5*avgIctalTime,avgIctalTime*5,1:length(session.data.channels));
        else
            error('Ictal/interictal overlap detected')
        end
            
    else
        % Testing interictal segments
        % pull interictal segments that lead up to seizure
        if 5*avgIctalTime < EECtimes{i}-szTimes{i-1}
            interSeg{i} = session.data.getvalues(EECtimes{i}-5*avgIctalTime,5*avgIctalTime,1:length(session.data.channels));

        else
            % Remove seizure
            interSeg{i} = [];
            ictalSeg{i} = [];
            fprintf('seizure %d has been removed',i)
        end
    end
    

    
    % Remove NaNs
    if sum(sum(isnan(interSeg{i}))) > 0 || sum(sum(isinf(interSeg{i}))) > 0
        [row,col] = find(isnan(interSeg{i}));
        interSeg{i}(row,:) = [];
    end
    
end
fprintf('Done!\n')


%% Divide segments into one second clips

freq = round(session.data.sampleRate);
channels = struct(session.data.channels.label);

fprintf('Clipping ictal data into 1-second clips...')

% ictal clips
testcount=0;
traincount=0;
for i = 1:length(ictalSeg);
    nclips = floor(length(ictalSeg{i})/fs);
    if i <= round(.4*length(ictalSeg))
        for j = 1:nclips
            traincount = traincount +1;
            data = ictalSeg{i}(fs*(j-1)+1:fs*j,:)';
            latency = j-1;
            str = sprintf('Patient_%d_ictal_segment_%d',idx,traincount);
            %save(str,'channels','data','freq','latency')
        end
    else
        for j = 1:nclips
            testcount = testcount+1;
            data = ictalSeg{i}(fs*(j-1)+1:fs*j,:)';
            latency = j-1;
            str = sprintf('Patient_%d_test_segment_%d',idx,testcount);
            %save(str,'channels','data','freq','latency')
        end
        szendclipt(i-round(.4*length(ictalSeg))) = testcount;
    end
end
fprintf('Done!\n')

% interictal clips  
fprintf('Clipping interictal data into 1-second clips...')

cutoff = testcount;
cutstr = sprintf('ictalCutoff_%d',idx);
save(cutstr,'cutoff')

traincount=0;
for i = 1:length(interSeg);
    
    nclips = floor(length(interSeg{i})/fs);
    ran = randperm(nclips);
    if i <= round(.4*length(interSeg))
        for j = 1:round(nclips/2)
            
            traincount = traincount + 1;
            data = interSeg{i}(fs*(ran(j)-1)+1:fs*ran(j),:)';
            latency = j-1;
            str = sprintf('Patient_%d_interictal_segment_%d',idx,traincount);
            %save(str,'channels','data','freq')
        end
    else
        for j = 1:nclips
            testcount = testcount+1;
            data = interSeg{i}(fs*(j-1)+1:fs*j,:)';
            latency = j-1;
            str = sprintf('Patient_%d_test_segment_%d',idx,testcount);
            %save(str,'channels','data','freq')
        end
        nszendclipt(i - round(.4*length(interSeg))) = testcount;
    end
end
fprintf('Done!\n')

testszend = [szendclipt nszendclipt]; 
szendstr = sprintf('testszend_%d',idx);
save(szendstr,'testszend')

% clear variables to free RAM
clear