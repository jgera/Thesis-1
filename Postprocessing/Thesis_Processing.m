%%%%%%%%%%%%%%%%%%%%%%%%%%% THESIS PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% Patients in Study
Patients = {'Thesis_Patient_7';...
    'Thesis_Patient_8'};
 

addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('lib'))

%% SVM models
% Construct SVM models for each subject

szLabels = [];
earLabels = [];
szProb = [];
earProb = [];
testlabsz = cell(1,2);
testlabear = cell(1,2);
for pat = 1:length(Patients)
    
    % Load Training and Testing Feature Sets and Labels
    load(Patients{pat})
   
    % Train model for sz labeling
    svmModelsz = svmtrain(trainlabelssz',traindata,'-b 1');

    % Train model for early labeling
    svmModelearly = svmtrain(trainlabelsearly',traindata,'-b 1');
    
    % Test SVM models
    [testDataSVMsz{pat},testAccSVMsz{pat},probsz{pat}] = svmpredict(testlabelssz',testdata,svmModelsz,'-b 1');
    [testDataSVMearly{pat},testAccSVMearly{pat},probearly{pat}] = svmpredict(testlabelsearly',testdata,svmModelearly, '-b 1');
    
    % Concatenate Test Labels for all patients
    szLabels = [szLabels testlabelssz];
    earLabels = [earLabels testlabelsearly];
    
    % Concatenate Probabilities for all patients
    szProb = [szProb; probsz{pat}(:,1)];
    earProb = [earProb; probearly{pat}(:,1)];
    
    testlabsz{pat} = testlabelssz;
    testlabear{pat} = testlabelsearly;

end




%% Calculate Metrics

[Xsz1,Ysz1,Tsz1,AUCsz1] = perfcurve(testlabsz{1},probsz{1}(:,1),1);
[Xsz2,Ysz2,Tsz2,AUCsz2] = perfcurve(testlabsz{2},probsz{2}(:,1),1);
[C,I] = unique(Xsz1);
[C2,I2] = unique(Xsz2);
Ysz1 = spline(C,Ysz1(I),C2);
Ysz = mean([Ysz1';Ysz2(I2)']);
Xsz = unique(Xsz2);
AUCsz = mean([AUCsz1 AUCsz2]);

[Xear1,Year1,Tear1,AUCear1] = perfcurve(testlabear{1},probearly{1}(:,1),1);
[Xear2,Year2,Tear2,AUCear2] = perfcurve(testlabear{2},probearly{2}(:,1),1);
[C,I] = unique(Xear1);
[C2,I2] = unique(Xear2);
Year1 = spline(C,Year1(I),C2);
Year = mean([Year1';Year2(I2)']);
Xear = unique(Xear2);
AUCear = mean([AUCear1 AUCear2]);    

% ROC Plot
figure(1)
plot(Xsz,Ysz)
hold on
plot(Xear,Year)
plot(linspace(0,1,length(Xear)),linspace(0,1,length(Xear)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Standard Detector')
axis([0 1 0 1])
hold off
    
[Xsz7,Ysz7,Tsz7,AUCsz7] = perfcurve(testlabsz{1},szProb(1:3004),1);
[Xear7,Year7,Tear7,AUCear7] = perfcurve(testlabear{1},earProb(1:3004),1);
% ROC Plot
figure(23)
plot(Xsz7,Ysz7)
hold on
plot(Xear7,Year7)
plot(linspace(0,1,length(Xear7)),linspace(0,1,length(Xear7)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Standard Detector Pat 7')
hold off
%% HILLS ALGORITHM

sub = 'Hills_2.csv';


    % load patient data
    Hillsprob = csvread(sub,1,1);

    % Seizure Probabilities
    HszProb = Hillsprob(:,1);

    % Early Probabilities
    HearProb = Hillsprob(:,2);
    
hillsszlabels = round(HszProb);
hillsearlabels = round(HearProb);

Hszlabels_7 = hillsszlabels(1:3004);
Hearlabels_7 = hillsearlabels(1:3004);

Hszlabels_8 = hillsszlabels(3005:end);
Hearlabels_8 = hillsearlabels(3005:end);

% Calculate Metrics

[HXsz1,HYsz1,HTsz1,HAUCsz1] = perfcurve(testlabsz{1},HszProb(1:3004),1);
[HXsz2,HYsz2,HTsz2,HAUCsz2] = perfcurve(testlabsz{2},HszProb(3005:end),1);
[C,I] = unique(HXsz1);
[C2,I2] = unique(HXsz2);
HYsz1 = spline(C,HYsz1(I),C2);
HYsz = mean([HYsz1';HYsz2(I2)']);
HXsz = unique(HXsz2);
HAUCsz = mean([HAUCsz1 HAUCsz2]);

[HXear1,HYear1,HTear1,HAUCear1] = perfcurve(testlabear{1},HearProb(1:3004),1);
[HXear2,HYear2,HTear2,HAUCear2] = perfcurve(testlabear{2},HearProb(3005:end),1);    
[C,I] = unique(HXear1);
[C2,I2] = unique(HXear2);
HYear1 = spline(C,HYear1(I),C2);
HYear = mean([HYear1';HYear2(I2)']);
HXear = unique(HXear2);
HAUCear = mean([HAUCear1 HAUCear2]);  


% ROC Plot
figure(2)
plot(HXsz,HYsz)
hold on
plot(HXear,HYear)
plot(linspace(0,1,length(HXear)),linspace(0,1,length(HXear)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Hills Detector')
axis([0 1 0 1])
hold off

[HXsz7,HYsz7,HTsz7,HAUCsz7] = perfcurve(testlabsz{1},HszProb(1:3004),1);
[HXear7,HYear7,HTear7,HAUCear7] = perfcurve(testlabear{1},HearProb(1:3004),1);
% ROC Plot
figure(20)
plot(HXsz7,HYsz7)
hold on
plot(HXear7,HYear7)
plot(linspace(0,1,length(HXear7)),linspace(0,1,length(HXear7)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Hills Detector Pat 7')
hold off
%% OLSON & MINGLE ALGORITHM

sub = 'olson_submission_p7.csv';


    % load patient data
    OMprob = csvread(sub,1,1);

    % Seizure Probabilities
    OMszProb = OMprob(:,1);

    % Early Probabilities
    OMearProb = OMprob(:,2);

% Create labels
OMszlabels = round(OMszProb);
OMearlabels = round(OMearProb);

% Patient 7 labels
OMszlabels_7 = OMszlabels(1:3004);
OMearlabels_7 = OMearlabels(1:3004);


% Calculate Metrics

[OMXsz,OMYsz,OMTsz,OMAUCsz] = perfcurve(testlabsz{1},OMszProb,1);
[OMXear,OMYear,OMTear,OMAUCear] = perfcurve(testlabear{1},OMearProb,1);
    

% ROC Plot
figure(3)
plot(OMXsz,OMYsz)
hold on
plot(OMXear,OMYear)
plot(linspace(0,1,length(OMXear)),linspace(0,1,length(OMXear)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Olson & Mingle Detector')
hold off


[OMXsz7,OMYsz7,OMTsz7,OMAUCsz7] = perfcurve(testlabsz{1},OMszProb(1:3004),1);
[OMXear7,OMYear7,OMTear7,OMAUCear7] = perfcurve(testlabear{1},OMearProb(1:3004),1);
% ROC Plot
figure(21)
plot(OMXsz7,OMYsz7)
hold on
plot(OMXear7,OMYear7)
plot(linspace(0,1,length(OMXear7)),linspace(0,1,length(OMXear7)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Olson & Mingle Detector Pat 7')
hold off

%% CDIPSTER ALGORITHM

sub = 'cdip_submission.csv';


    % load patient data
    Cdipprob = csvread(sub,1,1);

    % Seizure Probabilities
    CszProb = Cdipprob(:,1);

    % Early Probabilities
    CearProb = Cdipprob(:,2);
    
Cdipszlabels = round(CszProb);
Cdipearlabels = round(CearProb);

Cszlabels_7 = Cdipszlabels(1:3004);
Cearlabels_7 = Cdipearlabels(1:3004);

Cszlabels_8 = Cdipszlabels(3005:end);
Cearlabels_8 = Cdipearlabels(3005:end);

% Calculate Metrics

[CXsz1,CYsz1,CTsz1,CAUCsz1] = perfcurve(testlabsz{1},CszProb(1:3004),1);
[CXsz2,CYsz2,CTsz2,CAUCsz2] = perfcurve(testlabsz{2},CszProb(3005:end),1);
[C,I] = unique(CXsz1);
[C2,I2] = unique(CXsz2);
CYsz1 = spline(C,CYsz1(I),C2);
CYsz = mean([CYsz1';CYsz2(I2)']);
CXsz = unique(CXsz2);
CAUCsz = mean([CAUCsz1 CAUCsz2]);

[CXear1,CYear1,CTear1,CAUCear1] = perfcurve(testlabear{1},CearProb(1:3004),1);
[CXear2,CYear2,CTear2,CAUCear2] = perfcurve(testlabear{2},CearProb(3005:end),1); 
[C,I] = unique(CXear1);
[C2,I2] = unique(CXear2);
CYear1 = spline(C,CYear1(I),C2);
CYear = mean([CYear1';CYear2(I2)']);
CXear = unique(CXear2);
CAUCear = mean([CAUCear1 CAUCear2]);  

% ROC Plot
figure(4)
plot(CXsz,CYsz)
hold on
plot(CXear,CYear)
plot(linspace(0,1,length(CXear)),linspace(0,1,length(CXear)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Cdipster Detector')
axis([0 1 0 1])
hold off

[CXsz7,CYsz7,CTsz7,CAUCsz7] = perfcurve(testlabsz{1},CszProb(1:3004),1);
[CXear7,CYear7,CTear7,CAUCear7] = perfcurve(testlabear{1},CearProb(1:3004),1);
% ROC Plot
figure(22)
plot(CXsz7,CYsz7)
hold on
plot(CXear7,CYear7)
plot(linspace(0,1,length(CXear7)),linspace(0,1,length(CXear7)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Ictal/Nonictal','Early','Chance','Location','Best')
title('ROC Curves Cdipster Detector Pat 7')
hold off
%% Comparing Results (ROC)

% Seizure ROC Plot
figure(5)
plot(Xsz,Ysz)
hold on
plot(HXsz,HYsz)
plot(CXsz,CYsz)
plot(linspace(0,1,length(CXsz)),linspace(0,1,length(CXsz)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Standard','1st Place','3rd Place','Chance','Location','Best')
title('ROC Curves for Seizure Detection')
axis([0 1 0 1])
hold off

% Patient 7 ROC Sz
figure(24)
plot(Xsz7,Ysz7)
hold on
plot(HXsz7,HYsz7)
plot(OMXsz7,OMYsz7)
plot(CXsz7,CYsz7)
plot(linspace(0,1,length(CXsz)),linspace(0,1,length(CXsz)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Standard','1st Place','2nd Place','3rd Place','Chance','Location','Best')
title('ROC Curves for Seizure Detection - Validation Data Set 1')
axis([0 1 0 1])
hold off

% Early Seizure ROC Plot
figure(6)
plot(Xear,Year)
hold on
plot(HXear,HYear)
plot(CXear,CYear)
plot(linspace(0,1,length(CXear)),linspace(0,1,length(CXear)),'--k')
xlabel('FPR')
ylabel('TPR')
legend('Standard','Michael Hills','Cdipsters','Chance','Location','Best')
title('ROC Curves for Early Seizure Detection')
axis([0 1 0 1])
hold off

%% Comparing Results (Detection)

%                             PATIENT 8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--Standard--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('EEG_testsz1')
load('testszend_8')
load('ictalCutoff_8')

% Process channels to overlay on detection timing graphs
ictalsz1 = ictalsz1';
intersz1 = intersz1';

sz1 = [intersz1 ictalsz1];

% sz 19
sz1chan = [51 52 58 63 70];

% % sz 30
% sz1chan = [51 52 57 58 59];

% % sz 15
% sz1chan = [67 75];

% sz 16
% sz1chan = [52 53];

sz1 = sz1(sz1chan,:);

sz1 = bsxfun(@minus, sz1, mean(sz1,2));


a = -.1;
b = .1;

for i = 1:length(sz1chan)
    sz1(i,:) = a + (sz1(i,:) - min(sz1(i,:)))*(b-a)/(max(sz1(i,:))-min(sz1(i,:)));
end


% interictal start
%iiloc = find(testszend == cutoff);
iiloc = 26;
% Identify seizure start and end times
szstart = testDataSVMsz{2}(testszend(iiloc)+1:testszend(iiloc+1))'; 
szend = testDataSVMsz{2}(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Vector with classifications of sz or nonsz
szbarlocs = [szstart szend];
stanszbarlocs = szbarlocs;

% Identify nonseizure start and end times
nszstart = testDataSVMearly{2}(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = testDataSVMearly{2}(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Vector with classifications of early or nonearly
earbarlocs = [nszstart nszend];
stanearbarlocs = earbarlocs;

% Identify first seizure flag
[szB,szI] = find(diff(szbarlocs));

% % Connect to IEEG portal for Patient 8 data
% session = IEEGSession('TB Study 019','tblevins','tbl_ieeglogin');
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
% % UEO location in seizure plot
% UEOplottime = (UEOtimes{15}-EECtimes{15})/1e6;

% plot classification
figure(7)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Standard')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%--Michael Hills--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seizure start and end times
szstart = Hszlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
szend = Hszlabels_8(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Seizure locations
szbarlocs = [szstart szend];
Hszbarlocs = szbarlocs;

% Nonseizure start and end times
nszstart = Hearlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = Hearlabels_8(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Early bar locations
earbarlocs = [nszstart nszend];
Hearbarlocs = earbarlocs;

% Plot
figure(8)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Hills - Patient 8')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')
hold off
[hillsB,hillsI] = find(diff(szbarlocs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%--Olson & Mingle--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Seizure start and end times
% szstart = OMszlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
% szend = OMszlabels_8(1:testszend(1))';
% 
% % Seizure locations
% szbarlocs = [szstart szend];
% 
% % Nonseizure start and end times
% nszstart = OMearlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
% nszend = OMearlabels_8(1:testszend(1))';
% 
% % Early bar locations
% earbarlocs = [nszstart nszend];
% 
% % Plot
% figure(9)
% axis([-(length(szstart)) length(szend) 0 1])
% vline(find(szbarlocs==1)-(length(szstart)),'r')
% try
%     vline(find(earbarlocs==1)-(length(szstart)),'y')
% catch
% end
% vline(0,'k')
% vline(UEOplottime,'b')
% hold on
% plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(1,:) + .25,'k')
% plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(2,:) + .75,'k')
% set(gca,'YTickLabel',[]);
% xlabel('time from EEC (sec)')
% title('Olson & Mingle - Patient 8')
% hold off
% 
% [omB,omI] = find(diff(szbarlocs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--Cdipster--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Seizure start and end times
szstart = Cszlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
szend = Cszlabels_8(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Seizure locations
szbarlocs = [szstart szend];
cdipszbarlocs = szbarlocs;

% Nonseizure start and end times
nszstart = Cearlabels_8(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = Cearlabels_8(testszend(iiloc-22)+1:testszend(iiloc-22+1))';

% Early bar locations
earbarlocs = [nszstart nszend];
cdipearbarlocs = earbarlocs;
% Plot
figure(10)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Cdipster - Patient 8')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')
hold off

[cdipB,cdipI] = find(diff(szbarlocs));


%% Combined graph

% Plot
figure(11)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz_hills(find(Hszbarlocs==1)-(length(szstart)),'k')
vline_sz_standard(find(stanszbarlocs==1)-(length(szstart)),'k')
vline_sz_cdip(find(cdipszbarlocs==1)-(length(szstart)),'k')

vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz1)), sz1(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Algorithm Comparison - Patient 8')
text(-(length(szstart)),.975,'Hills')
text(-(length(szstart)),.925,'Standard')
text(-(length(szstart)),.875,'Cdipster')
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               PATIENT 7
load('EEG_testsz2')
load('testszend_7')
load('ictalCutoff_7')


% Process channels to overlay on detection timing graphs
ictalsz2 = ictalsz2';
intersz2 = intersz2';

sz2 = [intersz2 ictalsz2];

sz2chan = [1 2 3 4 13 17 36];

sz2 = sz2(sz2chan,:);

sz2 = bsxfun(@minus, sz2, mean(sz2,2));


a = -.05;
b = .05;

for i = 1:length(sz2chan)
    sz2(i,:) = a + (sz2(i,:) - min(sz2(i,:)))*(b-a)/(max(sz2(i,:))-min(sz2(i,:)));
end



%                                Standard

% interictal start
iiloc = 4;

% Identify seizure start and end times
szstart = testDataSVMsz{1}(testszend(iiloc)+1:testszend(iiloc+1))'; 
szend = testDataSVMsz{1}(1:testszend(iiloc-4+1))';

% Vector with classifications of sz or nonsz
szbarlocs = [szstart szend];
stanszbarlocs = szbarlocs;

% Identify nonseizure start and end times
nszstart = testDataSVMearly{1}(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = testDataSVMearly{1}(1:testszend(iiloc-4+1))';

% Vector with classifications of early or nonearly
earbarlocs = [nszstart nszend];

% Identify first seizure flag
[szB,szI] = find(diff(szbarlocs));

% % Connect to IEEG portal for Patient 7 data
% session = IEEGSession('TB Study 016','tblevins','tbl_ieeglogin');
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
% % UEO location in seizure plot
% UEOplottime = (UEOtimes{3}-EECtimes{3})/1e6;

% plot classification
figure(12)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Standard - Patient 7')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')




hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%--Michael Hills--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seizure start and end times
szstart = Hszlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
szend = Hszlabels_7(1:testszend(iiloc-4+1))';

% Seizure locations
szbarlocs = [szstart szend];
Hszbarlocs = szbarlocs;

% Nonseizure start and end times
nszstart = Hearlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = Hearlabels_7(1:testszend(iiloc-4+1))';

% Early bar locations
earbarlocs = [nszstart nszend];

% Plot
figure(13)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Hills - Patient 7')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')


hold off
[hillsB,hillsI] = find(diff(szbarlocs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%--Olson & Mingle--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Seizure start and end times
szstart = OMszlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
szend = OMszlabels_7(1:testszend(iiloc-4+1))';

% Seizure locations
szbarlocs = [szstart szend];
olsonszbarlocs = szbarlocs;

% Nonseizure start and end times
nszstart = OMearlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = OMearlabels_7(1:testszend(iiloc-4+1))';

% Early bar locations
earbarlocs = [nszstart nszend];

% Plot
figure(14)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Olson & Mingle - Patient 7')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')




hold off
[olsonB,olsonI] = find(diff(szbarlocs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--Cdipster--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Seizure start and end times
szstart = Cszlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
szend = Cszlabels_7(1:testszend(iiloc-4+1))';

% Seizure locations
szbarlocs = [szstart szend];
cdipszbarlocs = szbarlocs;
% Nonseizure start and end times
nszstart = Cearlabels_7(testszend(iiloc)+1:testszend(iiloc+1))';
nszend = Cearlabels_7(1:testszend(iiloc-4+1))';

% Early bar locations
earbarlocs = [nszstart nszend];

% Plot
figure(15)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz(find(szbarlocs==1)-(length(szstart)),'k')
try
    vline_ear(find(earbarlocs==1)-(length(szstart)),'k')
catch
end
vline(0,'k')
% vline(UEOplottime,'b')
hold on
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(2,:) + .25,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(3,:) + .40,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(4,:) + .55,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(5,:) + .70,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Cdipster - Patient 7')
text(-(length(szstart)),.95,'Seizure')
text(-(length(szstart)),.85,'Early')


hold off

%% Combined graph

% Plot
figure(16)
axis([-(length(szstart)) length(szend) 0 1])
vline_sz_hills(find(Hszbarlocs==1)-(length(szstart)),'k')
vline_sz_standard(find(stanszbarlocs==1)-(length(szstart)),'k')
vline_sz_cdip(find(cdipszbarlocs==1)-(length(szstart)),'k')
vline_sz_olson(find(olsonszbarlocs==1)-(length(szstart)),'k')
ylim=get(gca,'YLim');
xlim=get(gca,'XLim');
% vline(UEOplottime,'b')
hold on
rectangle('Position',[xlim(1),0.05,xlim(2)-xlim(1),ylim(2)-.6],...
          'Curvature',0.4,...
          'LineWidth',2,...
          'FaceColor',[.9,.9,.9])
vline_eec(0,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(1,:) + .1,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(2,:) + .2,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(3,:) + .3,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(4,:) + .4,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(5,:) + .5,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(6,:) + .575,'k')
plot(linspace(-(length(szstart)), length(szend), length(sz2)), sz2(7,:) + .675,'k')
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
xlabel('time from EEC (sec)')
title('Seizure Detection Precision')
% text(-(length(szstart)),.975,'Hills')
% text(-(length(szstart)),.925,'Standard')
% text(-(length(szstart)),.875,'Cdipster')
% text(-(length(szstart)),.825,'Olson')



% label algos
% text(xlim(1)-40,ylim(2)-.038,'Hills',...
%    'VerticalAlignment','bottom',...
%    'HorizontalAlignment','left')
text(xlim(1)-80,ylim(2)-.188,'Standard',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
% text(xlim(1)-75,ylim(2)-.138,'Cdipster',...
%    'VerticalAlignment','bottom',...
%    'HorizontalAlignment','left')
% text(xlim(1)-55,ylim(2)-.188,'Olson',...
%    'VerticalAlignment','bottom',...
%    'HorizontalAlignment','left')

text(xlim(1)-81,ylim(2)-.038,'1st Place',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')

text(xlim(1)-86,ylim(2)-.088,'2nd Place',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')

text(xlim(1)-82,ylim(2)-.138,'3rd Place',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')

% label channels
text(xlim(1)-55,ylim(2)-.917,'RAS1',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-55,ylim(2)-.816,'RAS2',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-55,ylim(2)-.717,'RAS3',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-55,ylim(2)-.617,'RAS4',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-65,ylim(2)-.523,'RFG17',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-65,ylim(2)-.423,'RFG20',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
text(xlim(1)-55,ylim(2)-.337,'ROF4',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')

rectangle('Position',[xlim(1),0,xlim(2)-xlim(1),1],...
          'Curvature',0,...
          'LineWidth',1)
rectangle('Position',[xlim(1),0.9,xlim(2)-xlim(1),.05],...
          'Curvature',0,...
          'LineWidth',1)
rectangle('Position',[xlim(1),0.85,xlim(2)-xlim(1),.05],...
          'Curvature',0,...
          'LineWidth',1)
rectangle('Position',[xlim(1),0.8,xlim(2)-xlim(1),.05],...
          'Curvature',0,...
          'LineWidth',1)

hold off



figure(17)
bar([AUCsz,HAUCsz,OMAUCsz,CAUCsz])
ylabel('Seizure AUC')
Labels = {'Standard', '1st Place', '2nd Place', '3rd Place'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);

Hills = [HAUCsz1,HAUCsz2];
Cdip = [CAUCsz1,CAUCsz2];
Stand = [AUCsz1,AUCsz2];

group = [Stand',Hills',Cdip'];

[h,p] = ttest(Hills,Stand);

anovap = anova1(group);