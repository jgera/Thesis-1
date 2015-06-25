%% Extract Annotations from csv files
function [channels,times] = csv2annot(csvfile,stop,fs)

split = strsplit(csvfile,'Sz');
load('gstartTimes')
var = str2num(csvfile(10:12));
globalStartTimes = globStartTimes{var};
start = str2num(csvfile(16:18));

for i = start:stop
    if i < 10
        sz = sprintf('00%d.xlsx',i);
    else
        sz = sprintf('0%d.xlsx',i);
    end
    
    filename = (strcat(char(split(1)),'Sz', sz));
     
    C = importdata(filename);
    

    chanEEC = strsplit(C{2,1},{' ',', '},'CollapseDelimiters',true);
    chanUEO = strsplit(C{4,1},{' ',', '},'CollapseDelimiters',true);
    chanEEC(1) = [];
    chanUEO(1) = [];
    channels{i-start+1,:} = {chanEEC,chanUEO};

    % Start time of EEC in microseconds
    startTimeEEC = C{2,2}(4:end);
    [Y,M,D,H,MN,S] = datevec(startTimeEEC);
    EECtime = (S + MN * 60 + H * 60^2) * 10^6 + globalStartTimes(i);

    % Start time of UEO in microseconds
    startTimeUEO = C{4,2}(4:end);
    [Y,M,D,H,MN,S] = datevec(startTimeUEO);
    UEOtime = (S + MN * 60 + H * 60^2) * 10^6 + globalStartTimes(i);

    % Start time of TRN in microseconds
    startTimeTRN = C{3,2}(4:end);
    [Y,M,D,H,MN,S] = datevec(startTimeTRN);
    TRNtime = (S + MN * 60 + H * 60^2) * 10^6 + globalStartTimes(i);

    times(i-start+1,:) = [EECtime TRNtime UEOtime];
end

