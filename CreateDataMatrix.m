%% Create Data Table
%-------------------


% this scripts create a large data table with the outcomes after each
% perturabtion.
clear all;

%% Path to the datafiles

% add main folder and subfolders to repository
addpath(genpath(pwd));

% download the data from dropbox
DatPath = GetDataAfschrift2022();

% Read the Yaml files with subject information
[SubjStruct, SubjFolders, SubjPreFix, SubjID_Exo, mass, height, age ] = ...
    GetSubjInfo(fullfile(DatPath,'SubjectInformation'));


%% Flow control

% post processing unperturbed data
Set.CreateDatMatUnp = true;

% settings for post processing perturbed data
Set.CreateDatMatPert = true;

% time window for EMG analysis and average encoder angle
Set.ComputePertResponse_TimeWindow = [0 0.5];
Set.ComputeEncoderAngle_TimeWindow = [0 0.5];

% how to norm EMG data
Set.EMGNorm_Max = false; % max of gait cycle average EMG data
Set.EMGNorm_Mean = false; % mean of gait cycle average EMG data
Set.MeanEMG = true; % average of low pass filtered EMG data during zero impedance walking (this was used in the paper)

%% Event detection settings
treshold = 30;  % 30 N as (vertical) force treshold for event detection
dtOffPlate = 0.2; % swing duration is minimal 200 ms (to avoid wrong detection events)
treshold_dtStride = 3; % excluding strides longer than 3s

%% Initialize data matrix

% additional information on headers:
% - MuscleName: average data first n-seconds after perturbation onset
%             (time window based on Set.ComputePertResponse_TimeWindow)
% - MuscleName LStance: Average muscle activity during first Lefst stance
%              phase after perturbation onset
% - COMstance  displacement in walking direction of the pelvis marker from
%             right heelstrike (before perturb) till left heelstrike

% controller conditions
NamesConditions = {'NMC_Bal','NMC_Default','NMC_COM','MinImp'};

% set warnings off
warning off;

%% PostProcessing Unperturbed walking

if Set.CreateDatMatUnp
    disp('Postprocessing steady-state walking data');
    % pre allocatie data matrices
    nsubj = length(SubjStruct);
    NMCEMGStore = nan(100,nsubj,3,7,200); 
    ZeroImpEMGStore = nan(100,nsubj,7,200);
    NMCSpatioTempStore = nan(nsubj,3,5,200); % (1) = stride duraction, (2)dt_stance, (3) percStance, (4) percSwing, (5) percDS
    ZeroImpSpatioTempStore = nan(nsubj,5,200); % (1) = stride duraction, (2)dt_stance, (3) percStance, (4) percSwing, (5) percDS
    ExoAdaptation_NMC = nan(nsubj,3,200);
    ExoAdaptation_MinImp = nan(nsubj,200);
    for s=1:length(SubjStruct)
        sID = SubjStruct{s}.subject.id;
        FolderSteadyState = fullfile(DatPath,['subject_' num2str(sID)],'SteadyState');
        for f = 1:4
            if f<4
                % walking with NMC controller
                FileSteadyState = fullfile(FolderSteadyState,['NMC_Default' num2str(f) '.csv']);
            else
                FileSteadyState = fullfile(FolderSteadyState,'MinImp.csv');
            end

            if exist(FileSteadyState,'file')
                % read the .csv file
                T = readtable(FileSteadyState);

                % get EMG data
                EMG_norm = [T.Soleus_R T.Gastroc_R T.Tibialis_R T.Vastus_R, ...
                    T.Soleus_L T.Gastroc_L T.Tibialis_L];
                
                % detect gait cycle events
                [Event,Step] = GetSpatioTemporalParam(T.time,T.Fz_R,...
                    T.Fz_R,treshold,dtOffPlate);
                
                % get the gait cycle averages between minute 2 en 4 of the
                % trial
                lowEMG_gcL = NormCycle(T.time,Event.ths_l,EMG_norm,...
                    treshold_dtStride,true);
                [~,~,ntrials] = size(lowEMG_gcL);
                lowEMG_gcR = NormCycle(T.time,Event.ths_r,EMG_norm,...
                    treshold_dtStride,true);
                [~,~,ntrials] = size(lowEMG_gcR);                

                % get foot placement information
                Pelvis = [T.COMPelvis_x T.COMPelvis_y T.COMPelvis_z];
                FootL = [T.LTOE_x T.LTOE_y T.LTOE_z];
                FootR = [T.RTOE_x T.RTOE_y T.RTOE_z];
                StabInfo = getFootPlacementInfo(Event,T.time,Pelvis,FootL,FootR);

                % exoskeleton torque during adaptation
                tau_r = T.RA_JointTorque_Nm;
                tPert_Qual = Event.ths_r;
                [TauExoAv,~] = getAvergePertRStance(Event,tau_r,T.time,Event.ths_r(2:end-1));

                % save results
                ntrials = length(Step.StrideTime);
                iSel = 1:ntrials;
                if f<4
                    % muscle information
                    NMCEMGStore(:,s,f,5:7,1:ntrials) =  lowEMG_gcL(:,5:7,:);
                    NMCEMGStore(:,s,f,1:4,1:ntrials) =  lowEMG_gcR(:,1:4,:);
                    % average torques
                    ExoAdaptation_NMC(s,f,1:length(TauExoAv)) = TauExoAv;
                    % spatio-temporal information
                    NMCSpatioTempStore(s,f,1,1:ntrials) = Step.StrideTime(1:ntrials);
                    NMCSpatioTempStore(s,f,2,1:ntrials) = Step.dtStance_n(1:ntrials);
                    NMCSpatioTempStore(s,f,3,1:ntrials) = Step.PercDS_n(1:ntrials);
                    NMCSpatioTempStore(s,f,4,1:ntrials) = Step.PercStance_n(1:ntrials);
                    NMCSpatioTempStore(s,f,5,1:ntrials) = Step.PercSwing_n(1:ntrials);
                    NMCSpatioTempStore(s,f,6,1:ntrials) = StabInfo.ths_l.dFeet(1:ntrials);
                    NMCSpatioTempStore(s,f,7,1:ntrials) = StabInfo.ths_l.dFootPelvis(1:ntrials);
                else
                    % muscle information
                    ZeroImpEMGStore(:,s,5:7,1:ntrials) =  lowEMG_gcL(:,5:7,:);
                    ZeroImpEMGStore(:,s,1:4,1:ntrials) =  lowEMG_gcR(:,1:4,:);
                    % average torques
                    ExoAdaptation_MinImp(s,1:length(TauExoAv)) = TauExoAv;
                    % spatio-temporal information
                    ZeroImpSpatioTempStore(s,1,1:ntrials) = Step.StrideTime(iSel);
                    ZeroImpSpatioTempStore(s,2,1:ntrials) = Step.dtStance_n(iSel);
                    ZeroImpSpatioTempStore(s,3,1:ntrials) = Step.PercDS_n(iSel);
                    ZeroImpSpatioTempStore(s,4,1:ntrials) = Step.PercStance_n(iSel);
                    ZeroImpSpatioTempStore(s,5,1:ntrials) = Step.PercSwing_n(iSel);
                    ZeroImpSpatioTempStore(s,6,1:ntrials) = StabInfo.ths_l.dFeet(iSel);
                    ZeroImpSpatioTempStore(s,7,1:ntrials) = StabInfo.ths_l.dFootPelvis(iSel);
                end
            end
        end       
        disp([' Subject ' num2str(s) ' finished']);
    end

    % structure with information on multi-dim matrices
    Info_NMC_EMG.dim1 = 'norm gait cycle';
    Info_NMC_EMG.dim2 = SubjPreFix;
    Info_NMC_EMG.dim3= 'Qual File Number';
    Info_NMC_EMG.dim4 = {'Soleus R','Gastroc R','Tibialis R','Vastus R','Soleus L','Gastroc L','Tibialis L'};
    Info_NMC_EMG.dim5 = 'Trials';
    Info_NMC_EMG.script = 'Batch_AdaptationGeyerModel.m';

    Info_NMC_SpatioTemp.dim1 = SubjPreFix;
    Info_NMC_SpatioTemp.dim2= 'Qual File Number';
    Info_NMC_SpatioTemp.dim4 = 'Trials';
    Info_NMC_SpatioTemp.dim3 = '(1) = stride duraction, (2)dt_stance, (3) percStance, (4) percSwing, (5) percDS, (6) stride length, (7) distance foot pelvis';
    Info_NMC_SpatioTemp.script = 'Batch_AdaptationGeyerModel.m';

    Info_ZeroImp_SpatioTemp.dim1 = SubjPreFix;
    Info_ZeroImp_SpatioTemp.dim3 = 'Trials';
    Info_ZeroImp_SpatioTemp.dim2 = '(1) = stride duraction, (2)dt_stance, (3) percStance, (4) percSwing, (5) percDS, (6) stride length, (7) distance foot pelvis';
    Info_ZeroImp_SpatioTemp.script = 'Batch_AdaptationGeyerModel.m';

    InfoZeroImp.dim1 = 'norm gait cycle';
    InfoZeroImp.dim2 = SubjPreFix;
    InfoZeroImp.dim3 = {'Soleus R','Gastroc R','Tibialis R','Vastus R','Soleus L','Gastroc L','Tibialis L'};
    InfoZeroImp.dim4 = 'Trials';
    InfoZeroImp.script = 'Batch_AdaptationGeyerModel.m';

    % save as .mat file
    if ~isfolder(fullfile(DatPath,'ResultsFiles'))
        mkdir(fullfile(DatPath,'ResultsFiles'));
    end
    save(fullfile(DatPath,'ResultsFiles','EMG_Unperturbed_csv.mat'),'ZeroImpEMGStore',...
        'NMCEMGStore','Info_NMC_EMG','NMCSpatioTempStore',...
        'Info_NMC_SpatioTemp','Info_ZeroImp_SpatioTemp','ZeroImpSpatioTempStore',...
        'ExoAdaptation_MinImp','ExoAdaptation_NMC');
end

%% Postprocessing perturbed data
if Set.CreateDatMatPert
    disp('Postprocessing perturbed data');
    % headers for output matrix
    headers = {'Subject-id','Controller-id','Perturbation-direction',...
        'Soleus R','Gastroc R','Tibialis R','Vastus R','Soleus L','Gastroc L','Tibialis L',...
        'SwingTime','StrideLength','dPelvisFoot','COMstance','avExoTorque','av_AnkleAngleR', ...
        'Soleus R LStance','Gastroc R LStance','Tibialis R LStance','Vastus R LStance',...
        'Soleus L LStance','Gastroc L LStance','Tibialis L LStance','LstanceTime','WorkExoRstance'};
    headers_Unp = {'Subject-id','Controller-id','Soleus R','Gastroc R','Tibialis R',...
        'Vastus R','Soleus L','Gastroc L','Tibialis L'};

    % pre-allocate large 2D matrix with all outputs
    data = nan(10000,length(headers));
    ctr = 1; % row counter (each row is a perturbation trial)
    for s=1:length(SubjFolders)
        sID = SubjStruct{s}.subject.id;
        for i = 1:length(NamesConditions)
            % path to folder with perturbation trials
            FolderPertFiles = fullfile(DatPath,['subject_' num2str(sID)],NamesConditions{i});
            for PertDir = 1:2
                if PertDir == 1
                    PertName = 'Push';
                else
                    PertName = 'Pull';
                end
                % loop over all 20 iterations of perturbations
                for ip = 1:20
                    % test if the file exsists
                    DatFile = fullfile(FolderPertFiles,[PertName '_' num2str(ip) '.csv']);
                    if exist(DatFile,'file')
                        % read the .csv file with data
                        T = readtable(DatFile);

                        % perturbation is always at t=0s
                        tPert = 0; 

                        % get matrix with EMG data in time window
                        EMG_norm = [T.Soleus_R T.Gastroc_R T.Tibialis_R T.Vastus_R, ...
                            T.Soleus_L T.Gastroc_L T.Tibialis_L];
                        tWindow = Set.ComputePertResponse_TimeWindow;
                        iSel = T.time>(tWindow(1)+tPert) & T.time<(tWindow(2)+tPert); % perturbation onset is at T.time == 0
                        EMGResponse= nanmean(EMG_norm(iSel,:));

                        % detect gait cycle events
                        [Event,Step] = GetSpatioTemporalParam(T.time,T.Fz_L,...
                            T.Fz_R,treshold,dtOffPlate);

                        % movement pelvis on treadmill
                        Pelvis = [T.COMPelvis_x  T.COMPelvis_y T.COMPelvis_z];
                        PelvisPos_SingleStance= getPelvisMovementTreadmill(Event,T.time,Pelvis,tPert);

                        % get foot placement at first event after perturbation
                        Pelvis = [T.COMPelvis_x T.COMPelvis_y T.COMPelvis_z];
                        FootL = [T.LTOE_x T.LTOE_y T.LTOE_z];
                        FootR = [T.RTOE_x T.RTOE_y T.RTOE_z];
                        StabInfo = getFootPlacementInfo(Event,T.time,Pelvis,FootL,FootR,tPert);

                        % Get average exo torque stance phase
                        [TauExoAv,~] = getAvergePertRStance(Event,T.RA_JointTorque_Nm,T.time,tPert);

                        % get the average encoder angle
                        tWindow = tPert + Set.ComputeEncoderAngle_TimeWindow;
                        iSel = T.time>tWindow(1) & T.time<tWindow(2);
                        qa_mean = nanmean(T.RA_JointAngle_rad(iSel));

                        % swing duration
                        % toe-off around perturbation onset
                        [~,iMin] = min(abs(Event.tto_l-tPert));
                        t_toe_l = Event.tto_l(iMin);
                        % first left heelstrike after toe-off
                        t_hs_l = Event.ths_l(find(Event.ths_l>t_toe_l,1,'first'));
                        % store duration swing phase
                        SwingDuration = t_hs_l-t_toe_l;

                        % Lstance duration
                        % toe-off around perturbation onset
                        t0 = Event.ths_l(find(Event.ths_l>tPert,1,'first'));
                        % first left heelstrike after toe-off
                        tend = Event.tto_l(find(Event.tto_l>t0,1,'first'));
                        % store duration swing phase
                        LstanceDuration = tend-t0;

                        % Muscle activity during first left stance phase after perturbation onset
                        [~, nEMG] = size(EMG_norm);
                        EMG_LeftStance = nan(nEMG,1);
                        for iM = 1:nEMG
                            EMG_LeftStance(iM) = getAvergeLeftStanceAfterPert(Event,EMG_norm(:,iM),T.time,tPert);
                        end

                        % store outcomes in the large 2d matrix
                        %---------------------------------------

                        % general information
                        data(ctr,1) = s; % subject id
                        data(ctr,2) = i; % control condition
                        data(ctr,3) = PertDir; % perturbation direction (1 = push, 2 = pull)

                        % EMG response
                        iEMGs = find(strcmp(headers,'Soleus R'));
                        data(ctr,iEMGs:iEMGs+6) = EMGResponse;

                        % 'StrideLength'
                        iStrideLength = find(strcmp(headers,'StrideLength'));
                        if isstruct(StabInfo) && length(StabInfo.ths_l.dFeet)>1
                            data(ctr,iStrideLength) = StabInfo.ths_l.dFeet(1);
                        else
                            data(ctr,iStrideLength) = NaN;
                        end

                        % position foot w.r.t. pelvis
                        iSel = find(strcmp(headers,'dPelvisFoot'));
                        if isstruct(StabInfo) && length(StabInfo.ths_l.dFootPelvis)>1
                            data(ctr,iSel) = StabInfo.ths_l.dFootPelvis(1);
                        else
                            data(ctr,iSel) = NaN;
                        end

                        % movement COM stance phase
                        iSel = find(strcmp(headers,'COMstance'));
                        if length(PelvisPos_SingleStance)>1
                            data(ctr,iSel) = PelvisPos_SingleStance(1);
                        else
                            data(ctr,iSel) = NaN;
                        end

                        % average exoskeleton torque
                        iSel = strcmp(headers,'avExoTorque');
                        data(ctr,iSel) = TauExoAv;

                        % duration stride
                        iSel = strcmp(headers,'SwingTime');
                        data(ctr,iSel) = SwingDuration;

                        % average encoder angle
                        iSel = strcmp(headers,'av_AnkleAngleR');
                        data(ctr,iSel) = qa_mean;

                        % EMG left stance
                        iEMGs = find(strcmp(headers,'Soleus R LStance'));
                        data(ctr,iEMGs:iEMGs+6) = EMG_LeftStance;
                        ctr = ctr+1;
                    end
                end
            end
        end
        disp([' Subject ' num2str(s) ' finished']);
    end
    data(ctr:end,:) = [];

    % save the data
    if ~isfolder(fullfile(DatPath,'ResultsFiles'))
        mkdir(fullfile(DatPath,'ResultsFiles'));
    end
    save(fullfile(DatPath,'ResultsFiles','DataMatrix_PerResponse_csv.mat'),'data',...
        'headers');
end