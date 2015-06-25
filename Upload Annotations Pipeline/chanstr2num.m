%% Function to convert channel strings to (location) number

function [EECloc,UEOloc] = chanstr2num(channels,dataset)

comparison = cell(1,1);

% All channels in the dataset
for i = 1:length(dataset.channels)
    comparison{1}{i} = dataset.channels(i).label;
end

EECloc = cell(1,length(channels));
UEOloc = cell(1,length(channels));

% Find matches
for j = 1:length(channels)
     
    for k = 1:length(channels{j}{1})
        
        try    
            EECloc{j}(k) = find(strcmp(channels{j}{1}(k),comparison{1,1}));
        catch
        end
        
    end
    for k = 1:length(channels{j}{2})
        try
            UEOloc{j}(k) = find(strcmp(channels{j}{2}(k),comparison{1,1}));
        catch
        end
    end
end



