%% Script to create data structure for global seizure start times
clear; clc;

globStartTimes = cell(16,1);
directory = dir('SzTimesFix');
filenames = extractfield(directory,'name');

j=0;
k=1;
for i = 1:length(filenames)-2
    j=j+1;
    file = filenames{i+2};
    load(file)
    globStartTimes{k}(j) = evTimeStart;
    clear evTimeStart
    try
        if file(12) ~= filenames{i+3}(12)
            k = k+1;
            j = 0;
        end
    catch
    end
end
    
    