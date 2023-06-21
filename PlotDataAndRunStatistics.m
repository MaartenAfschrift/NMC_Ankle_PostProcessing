%% StatisticsResuls_CSVBased.m
%-----------------------------

% This script does statistics on the outcome variables to compare the
% subject and exoskeleton data (EMG, motion capture, exoskeleton sensors)
% between the three different controllers
%   - NMC_Default = default neuromuscular controller
%   - NMC_COM = neuromuscular controller with COM velocity feedback
%   - MinImp = minimal impedance controller

% Note that this script uses the generated .csv files and therefore the
% results might slightly deviate from the analysis done on the raw data 
% (as reported in the paper).
% The raw data was interpolated to 100Hz to share the data in a reasonably
% sized format. Note that EMG, forceplate and exoskeleton sensor data was
% sampled at 1000Hz. 

% clear all variables
clear all;

% add current path and subfolder to your matlab path
addpath(genpath(pwd));

% settings
Settings.UseMean = false; % use mean of median of the 10-12 repetitions of each perturbation per subject
Settings.NormToMinImpPert = false; % normalize muscle response to muscle resposne to perturbation in zero impedance mode
Settings.SaveFigures = false; % save the figures

% index and name controller conditions
ControlConditions = [3 2 4]; % (3= NMC_COM, 2 = NMC_Default, 4 = minimal impedance)
NamesConditions = {'GeyerBal','GeyerCOM','GeyerD','MinImp'};
NamesConditions_header = {' ','NMC COM','NMC default', 'Minimal impedance'};


%% Load the main data matrix

% path with all data
DatPath = GetDataAfschrift2022();

% get subject information from the .yaml files
[SubjStruct, SubjFolders, SubjPreFix, SubjID_Exo, mass, height, age ] = ...
    GetSubjInfo(fullfile(DatPath,'SubjectInformation'));

% data matrix perturbed walking (generated with CreateDataMatrix.m)
load(fullfile(DatPath,'ResultsFiles','DataMatrix_PerResponse_csv.mat'),'data','headers');

% data matrix steady-state walking (generated with CreateDataMatrix.m)
Adapt = load(fullfile(DatPath,'ResultsFiles','EMG_Unperturbed_csv.mat'));

%% write all results statistical tests in a .txt file

DiaryFile = fullfile(fullfile(DatPath,'Statistics.txt'));
if exist(DiaryFile,'file')
    delete(DiaryFile)
end
diary(DiaryFile)

% create a folder to save figures
if ~isfolder(fullfile(DatPath,'ResultsFigures'))
    mkdir(fullfile(DatPath,'ResultsFigures'));
end

%% Display Subject statistics

disp('subject statistics');
disp([ 'mass: ' num2str(nanmean(mass))  '(' num2str(std(mass)) ') kg'  ]);
disp([ 'height: ' num2str(nanmean(height))  '(' num2str(std(height)) ') m'  ]);
disp([ 'age: ' num2str(nanmean(age))  '(' num2str(std(age)) ') years'  ]);

%% Average muscle activity after perturbation (Figure 6 in paper)

% we run the analysis for
Muscles = {'Soleus R','Gastroc R','Tibialis R','Vastus R','Soleus L','Gastroc L','Tibialis L'};
h = figure('Name','MuscleResponse');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[150   127   710   594]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% select control conditions
ContrNam = NamesConditions(ControlConditions);

% colors for each control condition
Cols = [28/255 117/255 188/255; % NMC_Def
    114/255 191/255 111/255; % NMC_COM
    0.6 0.6 0.6]; % minimal impedance
mk = 4;

% pre-allocate matrix with subject average responses
DatAvStore = nan(7,3,2);

% print to txt file
disp('Statistics Muscle Response to perturbation')
disp(' ');
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over soleus, gastrocnemius and tibialis anterior muscle on the
    % right leg
    for iM = 1:3 % we run this analys
        % display muscle name
        disp(['  ' Muscles{iM}]);
        
        % get matrix with average/median of nperturbations for each subject
        Muscle = data(:,strcmp(headers,Muscles{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                % select perturabtion trials for specific subject in a
                % specific controller condition and perturbation direction
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(Muscle(iSel));
                else
                    DatComp(s,ctr) = nanmedian(Muscle(iSel));
                end
            end
        end
        DatAvStore(iM,:,iPertDir) = nanmean(DatComp);

        % norm to response minimal impedance
        DatCompNormMin = DatComp./DatComp(:,3);
        if length(ControlConditions) == 3 && Settings.NormToMinImpPert
            DatComp = DatComp./DatComp(:,3);
        end

        % create table for anova test
        within = table(ContrNam','VariableNames',{'Controller'});
        within.Controller = categorical(within.Controller);
        table_DatComp = array2table(DatComp,'VariableNames',{'V1','V2','V3'});
        rm = fitrm(table_DatComp,'V1-V3 ~ 1','WithinDesign',within);
        ranova_tbl = ranova(rm,'WithinModel','Controller');

        % Run the multiple comparison
        mc = multcompare(rm,'Controller','ComparisonType','tukey-kramer');

        % Mauchly’s test for sphericity
        tbl_sph = mauchly(rm);
        rm_eps = epsilon(rm);

        % display statistic results
        PercChange_NMC_COM_MinImp = nanmean(DatCompNormMin(:,2))-1;
        PercChange_NMC_Default_MinImp =  nanmean(DatCompNormMin(:,1))-1;
        disp('    Relative change Muscle activity');
        disp(['      NMC-Default / MinImp ', num2str(round(PercChange_NMC_Default_MinImp*100,2))]);
        disp(['      NMC-COM / MinImp ', num2str(round(PercChange_NMC_COM_MinImp*100,2))]);
        disp('    Stats');
        disp(['     F(' num2str(ranova_tbl.DF(3)) ',' num2str(ranova_tbl.DF(4)) ') = ' num2str(ranova_tbl.F(3))   ', p =' num2str(ranova_tbl.pValue(3))])
        disp(['       P(NMC-Default / NMC-COM) = ' num2str(round(mc.pValue(1),3)) ]);
        disp(['       P(NMC-COM / MinImpedance)  = ' num2str(round(mc.pValue(2),3)) ]);
        disp(['       P(NMC-Default / MinImpedance) = ' num2str(round(mc.pValue(4),3)) ]);

        % Plot data and run statistics
        subplot(2,3,iM+(iPertDir-1)*3)
        for i=1:length(ControlConditions)
            PlotBar(i*2-1,DatComp(:,i),Cols(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;

        set(gca,'box','off')
        set(gca,'XTick',[1 3 5]);
        set(gca,'XTickLabel',NamesConditions(ControlConditions));
        set(gca,'XTickLabelRotation',60);
        title(Muscles(iM))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5)

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_MuscleResponse_.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_MuscleResponse_.svg'),'svg');
end

% display descriptive statics (stored in the log file)
disp(' ');
disp('NMC-COM w.r.t. NMC-Default  ');
PercChange = (DatAvStore(:,2,:)-DatAvStore(:,1,:))./DatAvStore(:,1,:);
PertDir = {'Push', 'Pull'};
for ip = 1:2
    disp(['Percentage change in muscle activity after  ' PertDir{ip}])
    for m=1:length(PercChange)
        dsel = squeeze(PercChange(m,ip));
        disp(['   ' Muscles{m} '  ' num2str(round(dsel*100,1))]);
    end
end

disp(' ');
disp('NMC-COM w.r.t. Minimal impedance  ');
PercChange = (DatAvStore(:,2,:)-DatAvStore(:,3,:))./DatAvStore(:,3,:);
PertDir = {'Push', 'Pull'};
for ip = 1:2
    disp(['Percentage change in muscle activity after  ' PertDir{ip}])
    for m=1:length(PercChange)
        dsel = squeeze(PercChange(m,ip));
        disp(['   ' Muscles{m} '  ' num2str(round(dsel*100,1))]);
    end
end
disp(' ');


%% Adaptation to Geyer assistance: Figure 3 paper

% post adaptation values
h = figure('Name','EMG adaptation');
set(h,'color',[1, 1, 1]);
set(h,'Position',[318   476   763   303]);
mk = 4;

% get average zero impedance
ZeroImp_meangc= squeeze(nanmean(Adapt.ZeroImpEMGStore,1));
ZeroImp_meantrials= squeeze(nanmean(ZeroImp_meangc,3));

% get average post adaptation
NMC_meangc = squeeze(nanmean(Adapt.NMCEMGStore,1));
NMC_meantrials = squeeze(nanmean(NMC_meangc,4));

% exclude EMG data of some particular subjects due to sensor problems. See
% AnalyzeOutliers for more details.
% 1. Subject 7: Gastrocnemius - L:
ZeroImp_meantrials(7,6) = NaN;
NMC_meantrials(7,:,6) = NaN;


CGeyer = [143 25 28; 195 31 72; 240 52 36]./255;
DatAvStore = nan(3,4);
for iM = 1:3 % loop over the muscles
    % get the average muscle activity for the right leg
    MinImpRight = ZeroImp_meantrials(:,iM);
    dRight =  squeeze(NMC_meantrials(:,:,iM));
    % get the average muscle activity for the right leg
    MinImpLeft = ZeroImp_meantrials(:,iM+4);
    dLeft =  squeeze(NMC_meantrials(:,:,iM+4));
    % get average of left and right leg
    DatComp(:,1) = (MinImpRight + MinImpLeft)/2;
    DatComp(:,2:4) = (dRight+dLeft)/2;
    % run Shapiro-Wilk test to see if data is normally distributed
    [Wilk(1).H, Wilk(1).pValue, Wilk(1).W] = swtest(DatComp(:,1), 0.05);
    [Wilk(2).H, Wilk(2).pValue, Wilk(2).W] = swtest(DatComp(:,4), 0.05);
    % run statistics
    if Wilk(1).H == 0 && Wilk(2).H == 0
        [pairedttest,p,ci,stats] = ttest(DatComp(:,1),DatComp(:,4),0.05);
        StatHeader = ['t(' num2str(stats.df), ') = ' num2str(round(stats.tstat,4)) ', p = ' num2str(round(p,3))];
    else
        [p, hstat, stats] = signrank(DatComp(:,1),DatComp(:,4),'method','approximate');
        StatHeader = ['z = ' num2str(round(stats.zval,4)) ', p = ' num2str(round(p,4))];
    end

    % Plot data and run statistics
    subplot(1,3,iM)
    PlotBar(8,DatComp(:,1),[0.4 0.4 0.4],mk);
    PlotBar(1,DatComp(:,2),CGeyer(1,:),mk);
    PlotBar(3,DatComp(:,3),CGeyer(2,:),mk);
    PlotBar(5,DatComp(:,4),CGeyer(3,:),mk);
    yLineInput = [squeeze(DatComp(:,2)) squeeze(DatComp(:,3)) squeeze(DatComp(:,4)) squeeze(DatComp(:,1))];
    xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5 ones(nsubj,1)*8];
    line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;

    set(gca,'box','off')
    set(gca,'XTick',[1 3 5 8]);
    set(gca,'XTickLabel',{'Assist-start','Assist-mid','Assist-end','Zero-Imp'});
    set(gca,'XTickLabelRotation',60);
    title({Muscles{iM},StatHeader})
    set(gca,'FontSize',11);
    set(gca,'Box','off');
    set(gca,'LineWidth',1.5)
    set(gca,'YLim', [0 1.6]);
    DatAvStore(iM,:) = nanmean(DatComp);
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_EMGadaptation_v3.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_EMGadaptation_v3.svg'),'svg');
end

PercChange = (DatAvStore(:,4)-DatAvStore(:,1))./DatAvStore(:,1);
disp('Unperturbed walking - Percentage change in muscle activity: ')
for m=1:length(PercChange)
    disp(['   ' Muscles{m} '  ' num2str(PercChange(m))]);
end

%% Figure average exoskeleton torque after push and pull perturbations (Figure 4 paper)

% open figure
h = figure('Name','AvExoTorque');
set(h,'Position',[159   104   321   577]);
set(h,'color',[1, 1, 1]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
mk = 4;
% we only don't plot the minimal impedance response to the perturbation
ControlConditionsS =ControlConditions(1:2); 

% we want to plot the assistance during unperturbed walking
UnpAv_Assist = squeeze(nanmean(Adapt.ExoAdaptation_NMC,3));
UnpMean = squeeze(nanmean(UnpAv_Assist,2));

for iPertDir = 1:2 % loop over the two perturbation direciton
    % get matrix with average of nperturbations for each subject in
    % controller condtion GeyerD and Geyer COM
    Dsel= data(:,strcmp(headers,'avExoTorque'));
    DatComp= nan(length(subj_unique),2); ctr = 0;
    for iContr = ControlConditions
        ctr = ctr+1;
        for s=1:length(subj_unique)
            iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
            if Settings.UseMean
                DatComp(s,ctr) = nanmean(Dsel(iSel))./mass(s);
            else
                DatComp(s,ctr) = nanmedian(Dsel(iSel))./mass(s);
            end
        end
    end
    % we want geyer assistance during unperturbed walking
    DatComp(:,3) = UnpMean(:,1)./mass';

    % one way onava statistics
    DatCompNormMin = DatComp./DatComp(:,3);

    % create table for anova test
    within = table(ContrNam','VariableNames',{'Controller'});
    within.Controller = categorical(within.Controller);
    table_DatComp = array2table(DatComp,'VariableNames',{'V1','V2','V3'});
    rm = fitrm(table_DatComp,'V1-V3 ~ 1','WithinDesign',within);
    ranova_tbl = ranova(rm,'WithinModel','Controller');

    % Run the multiple comparison
    mc = multcompare(rm,'Controller','ComparisonType','tukey-kramer');

    % display statistic results
    PercChange_COM_MinImp = nanmean(DatCompNormMin(:,2))-1;
    PercChange_Def_MinImp =  nanmean(DatCompNormMin(:,1))-1;

    disp('    Relative change exoskeleton torque');
    disp(['      GeyerD / GeyerUnp ', num2str(round(PercChange_Def_MinImp*100,2))]);
    disp(['      NMC-COM / GeyerUnp ', num2str(round(PercChange_COM_MinImp*100,2))]);
    disp('    Stats');
    disp(['     F(' num2str(ranova_tbl.DF(3)) ',' num2str(ranova_tbl.DF(4)), ...
        ') = ' num2str(ranova_tbl.F(3))   ', p =' num2str(ranova_tbl.pValue(3))])
    disp(['       P(NMC-Default / NMC-COM) = ' num2str(round(mc.pValue(1),3)) ]);
    disp(['       P(NMC-COM / GeyerUnp)  = ' num2str(round(mc.pValue(2),3)) ]);
    disp(['       P(NMC-Default / GeyerUnp) = ' num2str(round(mc.pValue(4),3)) ]);

    % Plot data
    subplot(2,1,iPertDir)
    for i=1:length(ControlConditions)
        PlotBar(i*2-1,DatComp(:,i),Cols(i,:),mk);
    end
    if length(ControlConditions) == 2
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3];
    elseif length(ControlConditions) == 3
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
    end
    line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]);

    % label axis and so on
    set(gca,'box','off')
    if length(ControlConditionsS) == 2
        set(gca,'XTick',[1 3]);
    elseif length(ControlConditionsS) == 3
        set(gca,'XTick',[1 3 5]);
    end
    set(gca,'XTickLabel',NamesConditions(ControlConditionsS));
    set(gca,'XTickLabelRotation',60);
    title('avExoTorque')
    if iM ==1
        if iPertDir ==1
            ylabel('PUSH');
        elseif iPertDir ==2
            ylabel('Pull');
        end
    end
    set(gca,'FontSize',11);
    set(gca,'Box','off');
    set(gca,'LineWidth',1.5)
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','AvExoMoment.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','AvExoMoment.svg'),'svg');
end

%% COM displacement after perturbation (Figure 5 paper)

OutStats = {'COMstance'};
OutStats_title = {'COMstance'};

% loop over all possible outcomes
h = figure('Name','COMstance');
set(h,'Position',[1166         502         258         622]);
set(h,'color',[1, 1, 1]);

Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
mk = 4;

disp(' ')
disp(' Statistics COM movement')
for iPertDir = 1:2 % loop over the two perturbation direciton
    disp([PertDirLab{iPertDir} ' perturbations']);
    for iM = 1:length(OutStats)

        % get matrix with average of nperturbations for each subject in
        % controller condtion GeyerD and Geyer COM
        Dsel= data(:,strcmp(headers,OutStats{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(Dsel(iSel));
                else
                    DatComp(s,ctr) = nanmedian(Dsel(iSel));
                end
            end
        end

        % one way onava statistics
        DatCompNormMin = DatComp./DatComp(:,3);

        % create table for anova test
        within = table(ContrNam','VariableNames',{'Controller'});
        within.Controller = categorical(within.Controller);
        table_DatComp = array2table(DatComp,'VariableNames',{'V1','V2','V3'});
        rm = fitrm(table_DatComp,'V1-V3 ~ 1','WithinDesign',within);
        ranova_tbl = ranova(rm,'WithinModel','Controller');

        % Run the multiple comparison
        mc = multcompare(rm,'Controller','ComparisonType','tukey-kramer');

        % Mauchly’s test for sphericity
        tbl_sph = mauchly(rm);
        rm_eps = epsilon(rm);

        % display statistic results
        PercChange_COM_MinImp = nanmean(DatCompNormMin(:,2))-1;
        PercChange_Default_MinImp =  nanmean(DatCompNormMin(:,1))-1;

        disp('    Relative change COM movement');
        disp(['      NMC-Default / MinImp ', num2str(round(PercChange_Default_MinImp*100,2))]);
        disp(['      NMC-COM / MinImp ', num2str(round(PercChange_COM_MinImp*100,2))]);
        disp('    Absolute change COM movement');
        disp(['      NMC-Default / MinImp ', num2str(nanmean(DatComp(:,1)-DatComp(:,3) ))]);
        disp(['      NMC-COM / MinImp ', num2str(nanmean(DatComp(:,2)-DatComp(:,3) ))]);
        disp('    Stats');
        disp(['     F(' num2str(ranova_tbl.DF(3)) ',' num2str(ranova_tbl.DF(4)) ') = ' num2str(ranova_tbl.F(3))   ', p =' num2str(ranova_tbl.pValue(3))])
        disp(['       P(NMC-Default / NMC-COM) = ' num2str(round(mc.pValue(1),3)) ]);
        disp(['       P(NMC-COM / MinImpedance)  = ' num2str(round(mc.pValue(2),3)) ]);
        disp(['       P(NMC-Default / MinImpedance) = ' num2str(round(mc.pValue(4),3)) ]);


        % Plot data 
        subplot(2,length(OutStats),iM+(iPertDir-1)*length(OutStats))
        for i=1:length(ControlConditions)
            PlotBar(i*2-1,DatComp(:,i),Cols(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]);

        % label axis and so on
        set(gca,'box','off')
        set(gca,'XTick',[1 3 5]);
        set(gca,'XTickLabel',NamesConditions(ControlConditions));
        set(gca,'XTickLabelRotation',60);
        title(OutStats_title{iM})
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5);
    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','COMstance.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','COMstance.svg'),'svg');
end


%% Adaptation to Geyer assistance (old version of figure 2)
% this old version of figure 2 runs a separate statistical test on the left
% and right leg (which fortunatly does not change the results as we do not
% expect left-right differences during the steady-state walking)

% Figure 
h = figure('Name','EMG adaptation - vOld');
set(h,'color',[1, 1, 1]);
set(h,'Position',[33 124 2010 303]);
mk = 4;

% get average zero impedance
ZeroImp_meangc= squeeze(nanmean(Adapt.ZeroImpEMGStore,1));
ZeroImp_meantrials= squeeze(nanmean(ZeroImp_meangc,3));

% get average post adaptation
Geyer_meangc = squeeze(nanmean(Adapt.NMCEMGStore,1));
Geyer_meantrials = squeeze(nanmean(Geyer_meangc,4));

% output header
OutStats ={'Adaptation start','Adaptation mid','Adaptation end'};

% exclude EMG data of some particular subjects due to sensor problems. See
% AnalyzeOutliers for more details.
% 1. Subject 7: Gastrocnemius - L:
ZeroImp_meantrials(7,6) = NaN;
Geyer_meantrials(7,:,6) = NaN;
% problem with minimal impedance values in subject 2 (reason is unclear). W
% e had to we had to re-solder the cable of the reference electrode in this
% test subject.
% (I'm probably using the wrong Qualisys file for the zero impedance walking)
% ZeroImp_meantrials(2,:) = NaN;
% Geyer_meantrials(2,:,:) = NaN;
ZeroImp_meantrials(2,4) = NaN;
Geyer_meantrials(2,:,4) = NaN;


CGeyer = copper(5);
DatAvStore = nan(7,4);
for iM = 1:7
    % get the data
    DatComp(:,1) = ZeroImp_meantrials(:,iM);
    DatComp(:,2:4) = squeeze(Geyer_meantrials(:,:,iM));

    [Wilk(1).H, Wilk(1).pValue, Wilk(1).W] = swtest(DatComp(:,1), 0.05);
    [Wilk(2).H, Wilk(2).pValue, Wilk(2).W] = swtest(DatComp(:,4), 0.05);

    % run statistics
    if Wilk(1).H == 0 && Wilk(2).H == 0
        [pairedttest,p,ci,stats] = ttest(DatComp(:,1),DatComp(:,4),0.05);
        StatHeader = ['t(' num2str(stats.df), ') = ' num2str(round(stats.tstat,4)) ', p = ' num2str(round(p,3))];
    else
        [p, hstat, stats] = signrank(DatComp(:,1),DatComp(:,4),'method','approximate');
        StatHeader = ['z = ' num2str(round(stats.zval,4)) ', p = ' num2str(round(p,4))];
    end

    % Plot data and run statistics
    subplot(1,7,iM)
    PlotBar(8,DatComp(:,1),[0.4 0.4 0.4],mk);
    PlotBar(1,DatComp(:,2),CGeyer(2,:),mk);
    PlotBar(3,DatComp(:,3),CGeyer(3,:),mk);
    PlotBar(5,DatComp(:,4),CGeyer(4,:),mk);
    yLineInput = [squeeze(DatComp(:,2)) squeeze(DatComp(:,3)) squeeze(DatComp(:,4)) squeeze(DatComp(:,1))];
    xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5 ones(nsubj,1)*8];
    line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;

    set(gca,'box','off')
    set(gca,'XTick',[1 3 5 8]);
    set(gca,'XTickLabel',{'Assist-start','Assist-mid','Assist-end','Zero-Imp'});
    set(gca,'XTickLabelRotation',60);
    title({Muscles{iM},StatHeader})
    set(gca,'FontSize',11);
    set(gca,'Box','off');
    set(gca,'LineWidth',1.5)
    set(gca,'YLim', [0 1.6]);
    DatAvStore(iM,:) = nanmean(DatComp);

end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_EMGadaptation_v2.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_EMGadaptation_v2.svg'),'svg');
end

PercChange = (DatAvStore(:,4)-DatAvStore(:,1))./DatAvStore(:,1);
disp('Unperturbed walking - Percentage change in muscle activity: ')
for m=1:length(PercChange)
    disp(['   ' Muscles{m} '  ' num2str(PercChange(m))]);
end

%% Plot exploratory analysis -- "balance" (related?) outcomes reviewer 1

% question reviewer 1: plot additional balance related outcomes.

% so we explored several balance related outecomes. not sure if this really
% adds something but you can decide yourself :).

% we run the analysis for (you can use this script to plot any outcome from the data array)
PlotVars = {'StrideLength_x','dPelvisFoot_x','Foot_Xcom_x'};
PlotVars_header = {'\Delta StrideLength-x [m]','\Delta Pelvis-Foot-x [m]','\Delta Xcom-Foot-x [m]'};
iRef = [6 8 10];
BoolRelativeToUnp = true; % changes in outcomes w.r.t. unperturbed walking 
nVar = length(PlotVars);
h = figure('Name','Balance outcomes');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[181.0000  149.0000  628.0000  708.6667]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% select control conditions
ContrNam = NamesConditions(ControlConditions);

% colors for each control condition
mk = 4;

% pre-allocate matrix with subject average responses
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over all outcomes
    for iM = 1:nVar       
        % get matrix with average/median of nperturbations for each subject
        DatSel = data(:,strcmp(headers,PlotVars{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                    % select perturabtion trials for specific subject in a
                    % specific controller condition and perturbation direction
                    iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                    if Settings.UseMean
                        DatComp(s,ctr) = nanmean(DatSel(iSel));
                    else
                        DatComp(s,ctr) = nanmedian(DatSel(iSel));
                    end
            end
        end
        % there is an outlier in the foot position (wrong label marker in subject 5)
        DatComp(5,:) = NaN; % To Do check this in dataset.       

        % deviation from unperturbed values
        if BoolRelativeToUnp
            for s=1:length(subj_unique)
                DRef = squeeze(Adapt.ZeroImpSpatioTempStore(s,iRef(iM),:));
                DRef_mean = nanmean(DRef);
                DatComp(s,:) = DatComp(s,:)-DRef_mean;
            end
        end

        % Plot data and run statistics
        subplot(2,nVar,iM+(iPertDir-1)*nVar)
        for i=1:length(ControlConditions)
            PlotBar(i*2-1,DatComp(:,i),Cols(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;

        set(gca,'box','off')
        set(gca,'XTick',[1 3 5]);
        set(gca,'XTickLabel',NamesConditions_header(ControlConditions));
        set(gca,'XTickLabelRotation',60);
        title(PlotVars_header(iM))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5);
        set(gca,'YLim',[-0.2 0.2]);

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_StabilityOutcomes.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_StabilityOutcomes.svg'),'svg');
end

disp(' ');


%% Plot exploratory analysis -- exoskeleton work reviewer 1

% Reviewer 1 asked (the excelled) question if the exoskeleton compensated
% for the energetic effect of the perturbation.

% we run the analysis for (you can use this script to plot any outcome from the data array)
PlotVars = {'ExoWork R Rstance'};
nVar = length(PlotVars);
h = figure('Name','Exo work');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[150   366   932   355]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% colors for each control condition
mk = 4;

% pre-allocate matrix with subject average responses
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over all outcomes
    for iM = 1:nVar       
        % get matrix with average/median of nperturbations for each subject
        DatSel = data(:,strcmp(headers,PlotVars{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                % select perturabtion trials for specific subject in a
                % specific controller condition and perturbation direction
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(DatSel(iSel));
                else
                    DatComp(s,ctr) = nanmedian(DatSel(iSel));
                end
            end
        end    
        % add unperturbed values
        Exo_mean = squeeze(nanmedian(Adapt.ExoAdaptation_MinImp,2));
        Exo_mean_Work_MinImp = Exo_mean(:,2);
        Exo_mean = squeeze(nanmedian(Adapt.ExoAdaptation_NMC,3));
        Exo_mean_Work_NMC = Exo_mean(:,3,2);
        DatComp = [DatComp Exo_mean_Work_MinImp Exo_mean_Work_NMC];
        ColsSel = [Cols; Cols(3,:); Cols(1,:)];

        % Plot data
        subplot(1,2,iPertDir);
        xLoc = [1 3 5 9 11];
        for i=1:length(ControlConditions)+2
            PlotBar(xLoc(i),DatComp(:,i),ColsSel(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3)), ...
            squeeze(DatComp(:,4)) squeeze(DatComp(:,5))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5 ones(nsubj,1)*9 ones(nsubj,1)*11];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;        
        NamesConditions_temp = [NamesConditions_header(ControlConditions) {'Minimal impdance','NMC default'}];

        set(gca,'box','off')
        set(gca,'XTick',xLoc);
        set(gca,'XTickLabel',NamesConditions_temp);
        set(gca,'XTickLabelRotation',60);
        title(PlotVars(iM))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5)

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Exo_WorkRStance.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Exo_WorkRStance.svg'),'svg');
end

disp(' ');

%% Plot exploratory analysis -- pusher work work reviewer 1

% plot work done by pusher

% we run the analysis for (you can use this script to plot any outcome from the data array)
PlotVars = {'WorkPusher'};
nVar = length(PlotVars);
h = figure('Name','Exo work');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[150   366   932   355]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% colors for each control condition
mk = 4;

% pre-allocate matrix with subject average responses
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over all outcomes
    for iM = 1:nVar       
        % get matrix with average/median of nperturbations for each subject
        DatSel = data(:,strcmp(headers,PlotVars{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                % select perturabtion trials for specific subject in a
                % specific controller condition and perturbation direction
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(DatSel(iSel));
                else
                    DatComp(s,ctr) = nanmedian(DatSel(iSel));
                end
            end
        end        

        % Plot data
        subplot(1,2,iPertDir);
        xLoc = [1 3 5];
        for i=1:length(ControlConditions)
            PlotBar(xLoc(i),DatComp(:,i),ColsSel(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;
        NamesConditions_temp = [NamesConditions_header(ControlConditions)];

        set(gca,'box','off')
        set(gca,'XTick',xLoc);
        set(gca,'XTickLabel',NamesConditions_temp);
        set(gca,'XTickLabelRotation',60);
        title(PlotVars(iM))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5)

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Pusher_WorkRStance.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Pusher_WorkRStance.svg'),'svg');
end
disp(' ');




%% Average muscle activity during left stance phase after perturbation

% this was a question of reviewer 1: is muscle activity different in the
% contra-lateral limb ?

% we run the analysis for
Muscles = {'Soleus R LStance','Gastroc R LStance','Tibialis R LStance',...
    'Vastus R LStance','Soleus L LStance','Gastroc L LStance','Tibialis L LStance'};
h = figure('Name','MuscleResponse L Stance');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[150   127   710   594]);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% select control conditions
ContrNam = NamesConditions(ControlConditions);

% colors for each control condition
Cols = [28/255 117/255 188/255; % NMC_Def
    114/255 191/255 111/255; % NMC_COM
    0.6 0.6 0.6]; % minimal impedance
mk = 4;

% pre-allocate matrix with subject average responses
DatAvStore = nan(7,3,2);

% print to txt file
disp('Statistics Muscle Response to perturbation')
disp(' ');
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over soleus, gastrocnemius and tibialis anterior muscle on the
    % right leg
    for iM = 1:3 % we run this analys
        % display muscle name
        disp(['  ' Muscles{iM+4}]);
        
        % get matrix with average/median of nperturbations for each subject
        Muscle = data(:,strcmp(headers,Muscles{iM+4}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                % select perturabtion trials for specific subject in a
                % specific controller condition and perturbation direction
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(Muscle(iSel));
                else
                    DatComp(s,ctr) = nanmedian(Muscle(iSel));
                end
            end
        end
        DatAvStore(iM,:,iPertDir) = nanmean(DatComp);

        % norm to response minimal impedance
        DatCompNormMin = DatComp./DatComp(:,3);
        if length(ControlConditions) == 3 && Settings.NormToMinImpPert
            DatComp = DatComp./DatComp(:,3);
        end


        % create table for anova test
        within = table(ContrNam','VariableNames',{'Controller'});
        within.Controller = categorical(within.Controller);
        table_DatComp = array2table(DatComp,'VariableNames',{'V1','V2','V3'});
        rm = fitrm(table_DatComp,'V1-V3 ~ 1','WithinDesign',within);
        ranova_tbl = ranova(rm,'WithinModel','Controller');

        % Run the multiple comparison
        mc = multcompare(rm,'Controller','ComparisonType','tukey-kramer');

        % Mauchly’s test for sphericity
        tbl_sph = mauchly(rm);
        rm_eps = epsilon(rm);

        % display statistic results
        PercChange_NMC_COM_MinImp = nanmean(DatCompNormMin(:,2))-1;
        PercChange_NMC_Default_MinImp =  nanmean(DatCompNormMin(:,1))-1;
        disp('    Relative change Muscle activity');
        disp(['      NMC-Default / MinImp ', num2str(round(PercChange_NMC_Default_MinImp*100,2))]);
        disp(['      NMC-COM / MinImp ', num2str(round(PercChange_NMC_COM_MinImp*100,2))]);
        disp('    Stats');
        disp(['     F(' num2str(ranova_tbl.DF(3)) ',' num2str(ranova_tbl.DF(4)) ') = ' num2str(ranova_tbl.F(3))   ', p =' num2str(ranova_tbl.pValue(3))])
        disp(['       P(NMC-Default / NMC-COM) = ' num2str(round(mc.pValue(1),3)) ]);
        disp(['       P(NMC-COM / MinImpedance)  = ' num2str(round(mc.pValue(2),3)) ]);
        disp(['       P(NMC-Default / MinImpedance) = ' num2str(round(mc.pValue(4),3)) ]);

        % Plot data and run statistics
        subplot(2,3,iM+(iPertDir-1)*3)
        for i=1:length(ControlConditions)
            PlotBar(i*2-1,DatComp(:,i),Cols(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;

        set(gca,'box','off')
        set(gca,'XTick',[1 3 5]);
        set(gca,'XTickLabel',NamesConditions(ControlConditions));
        set(gca,'XTickLabelRotation',60);
        title(Muscles(iM+4))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5)

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_MuscleResponse_.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','Stats_MuscleResponse_.svg'),'svg');
end

disp(' ');
disp('NMC-COM w.r.t. NMC-Default  ');
PercChange = (DatAvStore(:,2,:)-DatAvStore(:,1,:))./DatAvStore(:,1,:);
PertDir = {'Push', 'Pull'};
for ip = 1:2
    disp(['Percentage change in muscle activity after  ' PertDir{ip}])
    for m=1:length(PercChange)
        dsel = squeeze(PercChange(m,ip));
        disp(['   ' Muscles{m} '  ' num2str(round(dsel*100,1))]);
    end
end

disp(' ');
disp('NMC-COM w.r.t. Minimal impedance  ');
PercChange = (DatAvStore(:,2,:)-DatAvStore(:,3,:))./DatAvStore(:,3,:);
PertDir = {'Push', 'Pull'};
for ip = 1:2
    disp(['Percentage change in muscle activity after  ' PertDir{ip}])
    for m=1:length(PercChange)
        dsel = squeeze(PercChange(m,ip));
        disp(['   ' Muscles{m} '  ' num2str(round(dsel*100,1))]);
    end
end
disp(' ');

%% Plot COM position at first 4 heelstrikes after perturbation onset -- reviewers question

% the reviewer asked if subjects need multiple steps to return to the
% treadmill origin.

% we run the analysis for (you can use this script to plot any outcome from the data array)
PlotVars = {'COM_hs1','COM_hs2','COM_hs3','COM_hs4'};
PlotVars_headers = {'heel strike 1','heel strike 2','heel strike 3','heel strike 4'};
nVar = length(PlotVars);
h = figure('Name','COM nhs');

% open figure
set(h,'color',[1, 1, 1]);
set(h,'Position',[0.0790    0.3660    1.0030    0.5637]*1000);

% select some specific columns in the data matrix
Subject = data(:,strcmp(headers,'Subject-id'));
Controller = data(:,strcmp(headers,'Controller-id'));
PertDir = data(:,strcmp(headers,'Perturbation-direction'));
subj_unique = unique(Subject);
nsubj = length(subj_unique);

% select control conditions
ContrNam = NamesConditions(ControlConditions);

% colors for each control condition
mk = 4;

% pre-allocate matrix with subject average responses
PertDirLab = {'Push', 'Pull'};

% loop over push and pull perturbations
for iPertDir = 1:2 % push and pull perturbations
    disp(' ');
    disp([PertDirLab{iPertDir} ' perturbations']);
    % loop over all outcomes
    for iM = 1:nVar       
        % get matrix with average/median of nperturbations for each subject
        DatSel = data(:,strcmp(headers,PlotVars{iM}));
        DatComp= nan(length(subj_unique),2); ctr = 0;
        for iContr = ControlConditions
            ctr = ctr+1;
            for s=1:length(subj_unique)
                % select perturabtion trials for specific subject in a
                % specific controller condition and perturbation direction
                iSel = Subject == subj_unique(s) & Controller == iContr & PertDir == iPertDir;
                if Settings.UseMean
                    DatComp(s,ctr) = nanmean(DatSel(iSel));
                else
                    DatComp(s,ctr) = nanmedian(DatSel(iSel));
                end
            end
        end         

        % Plot data
        subplot(2,nVar,(iPertDir-1)*nVar+iM);
        xLoc = [1 3 5];
        for i=1:length(ControlConditions)
            PlotBar(xLoc(i),DatComp(:,i),ColsSel(i,:),mk);
        end
        yLineInput = [squeeze(DatComp(:,1)) squeeze(DatComp(:,2)) squeeze(DatComp(:,3))];
        xLineInput = [ones(nsubj,1) ones(nsubj,1)*3 ones(nsubj,1)*5];
        line(xLineInput', yLineInput','Color',[0.5 0.5 0.5]); hold on;
        NamesConditions_temp = [NamesConditions_header(ControlConditions)];

        set(gca,'box','off')
        set(gca,'XTick',xLoc);
        set(gca,'XTickLabel',NamesConditions_temp);
        set(gca,'XTickLabelRotation',60);
        title(PlotVars_headers(iM))
        if iM ==1
            if iPertDir ==1
                ylabel('PUSH');
            elseif iPertDir ==2
                ylabel('Pull');
            end
        end
        set(gca,'FontSize',11);
        set(gca,'Box','off');
        set(gca,'LineWidth',1.5);
        if iPertDir == 1
            set(gca,'Ylim',[0, 0.25]);
        else
            set(gca,'Ylim',[-0.25, 0]);
        end

    end
end
if Settings.SaveFigures
    saveas(h,fullfile(DatPath,'ResultsFigures','COM_movement_nhs.fig'),'fig');
    saveas(h,fullfile(DatPath,'ResultsFigures','COM_movement_nhs.svg'),'svg');
end

disp(' ');


%% end diary
diary off;
