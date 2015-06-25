%% Standard Seizure Detection Algorithm for variable training set analysis
% <latex>
% \title{Litt Lab - Standard Seizure Detector\\{\normalsize Spring 2015}}
% \author{Tyler Blevins}
% \date{\today}
% \maketitle
% </latex>

% clear the workspace and console
clear all; close all; clc;
warning('off')

% This seizure detection algorithm is based off of the detector described
% in "One-Class Novelty Detection for Seizure Analysis from Intracranial
% EEG" by Gardner et al.

% Patient #
pn = 8; 

% number of seizures in training set
n = 3;

fprintf('Loading Library...')
% add path to clips
addpath(genpath('lib'))


pathstr = sprintf('Patient_%d',pn);
addpath(genpath(pathstr))
fprintf('Done! \n')

%% Features
% Mean curve length feature function
CL = @(x) log(sum(abs(diff(x)))/length(x));

% Mean energy
E = @(x) log(sum(x.^2)/length(x));


% Mean Teager Energy
% see function file TE


%% Process Training Clips
% Process patient training clips

% load filter
load('hum500')

fprintf('Processing Patient training clips...')

i=0;
while 1 == 1
    i=i+1;
    try
        % load ictal data
        str = sprintf('Patient_%d_ictal_segment_%d.mat',pn,i);
        data = load(str);
        latency = data.latency;
        data = data.data;
    catch
        break
    end
    % filter data and downsample to 400hz
    data = filtfilt(hum500,1,data')';
    data = downsample(data',round(size(data,2)/200))';
    cldat = [];
    edat = [];
    TEdat = [];
    
    % Feature Extraction
    for c = 1:size(data,1)
        cldat = [cldat CL(data(c,:))];
        edat = [edat E(data(c,:))];
        TEdat = [TEdat TE(data(c,:))];
    end
    traindata(i,:) = [cldat edat TEdat];
    
    % Construct Labels
    trainlabelssz(i) = 1;
    if latency <= 15
        trainlabelsearly(i) = 1;
    else
        trainlabelsearly(i) = 0;
    end

end


                       
% interictal segments

i=0;
while 1 == 1
    i=i+1;
    try
        % load interictal data
        str = sprintf('Patient_%d_interictal_segment_%d.mat',pn,i);
        data = load(str);
    catch
        break
    end
    
    data = data.data;
    data = filtfilt(hum500,1,data')';
    data = downsample(data',round(size(data,2)/200))';
    cldat = [];
    edat = [];
    TEdat = [];
    
    % Feature Extraction
    for c = 1:size(data,1)
        cldat = [cldat CL(data(c,:))];
        edat = [edat E(data(c,:))];
        TEdat = [TEdat TE(data(c,:))];
    end
    traindata(end+1,:) = [cldat edat TEdat];
    
    % Construct Labels
    trainlabelssz(end+1) = -1;
    trainlabelsearly(end+1) = 0;
end

fprintf('Done! \n')


%% Process Testing clips

% Process patient testing clips
fprintf('Processing Patient testing clips...')

% load ictal/interictal changeover point
cutpoint = sprintf('ictalCutoff_%d',pn);
load(cutpoint)


i=0;
while 1 == 1
    i=i+1;
    try
        % load test segments
        str = sprintf('Patient_%d_test_segment_%d.mat',pn,i);
        data = load(str);
    catch
        break
    end
    try
        latency = data.latency;
    catch
        latency = 16;
    end
    data = data.data;  
    data = filtfilt(hum500,1,data')';
    data = downsample(data',round(size(data,2)/200))';
   
    cldat = [];
    edat = [];
    TEdat = [];
    
    % Feature Extraction
    for c = 1:size(data,1)
        cldat = [cldat CL(data(c,:))];
        edat = [edat E(data(c,:))];
        TEdat = [TEdat TE(data(c,:))];
    end
    testdata(i,:) = [cldat edat TEdat];
    if i<= cutoff
        testlabelssz(i) = 1;
    else
        testlabelssz(i) = -1;
    end
    try
        if latency <= 15
            testlabelsearly(i) = 1;
        else
            testlabelsearly(i) = 0;
        end
    catch
        testlabelsearly(i) = 0;
    end

end

fprintf('Done! \n')

% savestr = sprintf('Thesis_Patient_%d',pn);
% save(savestr,'traindata','trainlabelssz','trainlabelsearly','testdata','testlabelssz','testlabelsearly')
%% Train SVM models
% Construct SVM models for each subject


% Train model for sz labeling
svmModelsz = svmtrain(trainlabelssz',traindata,'-b 1');

% Train model for early labeling
svmModelearly = svmtrain(trainlabelsearly',traindata,'-b 1');

%%
% figure(1)
% scatter3(traindata(trainlabelssz==1,1),traindata(trainlabelssz==1,97),traindata(trainlabelssz==1,198))
% hold on
% scatter3(traindata(trainlabelssz==-1,1),traindata(trainlabelssz==-1,97),traindata(trainlabelssz==-1,193))
% hold off
% 
% [coeff, score, latent, tsquared, explained] = pca([traindata;testdata]);
% 
% figure(2)
% scatter3(score(trainlabelssz==1,1),score(trainlabelssz==1,2),score(trainlabelssz==1,3))
% hold on
% scatter3(score(trainlabelssz==-1,1),score(trainlabelssz==-1,2),score(trainlabelssz==-1,3))
% hold off


% % Train model for sz labeling
% svmModelsz = svmtrain(trainlabelssz',score(1:length(trainlabelssz),1:3),'-b 1');
% 
% % Train model for early labeling
% svmModelearly = svmtrain(trainlabelsearly',score(1:length(trainlabelssz),1:3),'-b 1');

%% Test SVM models

[testDataSVMsz,testAccSVMsz,probsz] = svmpredict(testlabelssz',testdata,svmModelsz,'-b 1');
[testDataSVMearly,testAccSVMearly,probearly] = svmpredict(testlabelsearly',testdata,svmModelearly, '-b 1');

% testDataknnsz = knnclassify(testdata,traindata,trainlabelssz');
% testDataknnearly = knnclassify(testdata,traindata,trainlabelsearly');  

%% Calculate AUC
% load('hills')
[Xsvmsz,Ysvmsz,Tsvmsz,AUCsvmsz] = perfcurve(testlabelssz,probsz(:,1),1);
[Xsvmear,Ysvmear,Tsvmear,AUCsvmear] = perfcurve(testlabelsearly,probearly(:,1),1);

str = sprintf('stan_n%d',n);
save(str,'Xsvmsz','Ysvmsz','AUCsvmsz')
% [Xhillssz,Yhillssz,Thillssz,AUChillssz] = perfcurve(testlabelssz,hillssz,1);
% [Xhillsear,Yhillsear,Thillsear,AUChillsear] = perfcurve(testlabelsearly,hillsear,1);
% 
% hillsszlabels = round(hillssz);
% hillsearlabels = round(hillsear);
% 
% figure(1)
% plot(Xsvmsz,Ysvmsz,Xhillssz,Yhillssz),legend('Standard','Michael Hills','Location','Best')
% xlabel('false positive rate')
% ylabel('true positive rate')
% title('ROC for Ictal vs NonIctal Classification')
% figure(2)
% plot(Xsvmear,Ysvmear,Xhillsear,Yhillsear),legend('Standard','Michael Hills','Location','Best')
% title('ROC for Early Seizure Classification')
% ScoreSVM = .5 * (AUCsvmsz + AUCsvmear);
% ScoreHills = .5 * (AUChillssz + AUChillsear);
% 
% load('testszend')
% 
% szstart = testDataSVMsz(testszend(16)+1:testszend(17))';
% szend = testDataSVMsz(1:testszend(1))';
% 
% szbarlocs = [szstart szend];
% 
% nszstart = testDataSVMearly(testszend(16)+1:testszend(17))';
% nszend = testDataSVMearly(1:testszend(1))';
% 
% earbarlocs = [nszstart nszend];
% 
% [szB,szI] = find(diff(szbarlocs));
% 
% 
% session = IEEGSession('TB Study 019','tblevins','tbl_IEEGlogin');
% 
% % Establish which layer is UEO
% layerNames = {session.data.annLayer.name};
% ind=find(ismember(layerNames,'KD_BO_UEO'));
% 
% % Get UEO start times
% UEOstart = session.data.annLayer(ind).getEvents(0);
% UEOtimes = {UEOstart.start};
% 
% % Establish which layer is EEC
% layerNames = {session.data.annLayer.name};
% ind=find(ismember(layerNames,'KD_BO_EEC'));
% 
% % Get EEC start times
% EECstart = session.data.annLayer(ind).getEvents(0);
% EECtimes = {EECstart.start};
% 
% UEOplottime = (UEOtimes{1}-EECtimes{1})/1e6;
% 
% figure(3)
% axis([-(length(earbarlocs)/2) length(szbarlocs)/2 0 1])
% vline(find(szbarlocs==1)-(length(szbarlocs)/2),'r')
% try
%     vline(find(earbarlocs==1)-(length(earbarlocs)/2),'y')
% catch
% end
% vline(0,'k')
% vline(UEOplottime,'b')
% set(gca,'YTickLabel',[]);
% xlabel('time from EEC (sec)')
% title('Standard')
% 
% 
% 
% szstart = hillsszlabels(testszend(16)+1:testszend(17))';
% szend = hillsszlabels(1:testszend(1))';
% 
% szbarlocs = [szstart szend];
% 
% nszstart = hillsearlabels(testszend(16)+1:testszend(17))';
% nszend = hillsearlabels(1:testszend(1))';
% 
% earbarlocs = [nszstart nszend];
% 
% figure(4)
% axis([-(length(earbarlocs)/2) length(szbarlocs)/2 0 1])
% vline(find(szbarlocs==1)-(length(szbarlocs)/2),'r')
% try
%     vline(find(earbarlocs==1)-(length(earbarlocs)/2),'y')
% catch
% end
% vline(0,'k')
% vline(UEOplottime,'b')
% set(gca,'YTickLabel',[]);
% xlabel('time from EEC (sec)')
% title('Hills')
% 
% [hillsB,hillsI] = find(diff(szbarlocs));
% 
% figure(5)
% axis([-(length(earbarlocs)/2) length(szbarlocs)/2 0 1])
% vline(hillsI(3)-(length(earbarlocs)/2),'m','Hills')
% vline(szI(1)-(length(earbarlocs)/2),'g')
% vline(0,'k')
% vline(UEOplottime,'b','UEO')
% set(gca,'YTickLabel',[]);
% xlabel('time from EEC (sec)')
