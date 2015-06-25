%% Upload Annotations to the Portal
function uploadAnnots(csvfile,stop,dataset)

	
%	csvfile     -	string of name of csvfile with annotations
%   stop        -   final seizure number for the patient
%	dataset		-	IEEGDataset object

% Sampling rate
fs = dataset.sampleRate;

% Extract annotations from CSV file
[channels,times] = csv2annot(csvfile,stop,fs);
[EECloc,UEOloc] = chanstr2num(channels,dataset);

%% EEC
layerName = 'KD_BO_EEC';
eventTimesUSec = times(:,1);
eventChannels = EECloc;
label = 'EEC';
uploadAnnotations_v2(dataset,layerName,eventTimesUSec,eventChannels,label)

%% TRN
layerName = 'KD_BO_TRN';
eventTimesUSec = times(:,2);
eventChannels = cell(1,length(times(:,2)));
for i = 1:length(eventChannels)
    eventChannels{i} = 1:length(dataset.channels);
end
label = 'TRN';
uploadAnnotations_v2(dataset,layerName,eventTimesUSec,eventChannels,label)

%% UEO
layerName = 'KD_BO_UEO';
eventTimesUSec = times(:,3);
eventChannels = UEOloc;
label = 'UEO';
uploadAnnotations_v2(dataset,layerName,eventTimesUSec,eventChannels,label)



