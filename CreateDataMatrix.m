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
% NMC_Bal is the 
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
    NMCSpatioTempStore = nan(nsubj,3,11,200);
    ZeroImpSpatioTempStore = nan(nsubj,11,200); 
    ExoAdaptation_NMC = nan(nsubj,3,200,2);
    ExoAdaptation_MinImp = nan(nsubj,200,2);
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
                [Event,Step] = GetSpatioTemporalParam(T.time,T.Fy_L,...
                    T.Fy_R,treshold,dtOffPlate);

                % get the gait cycle averages between minute 2 en 4 of the
                % trial
                lowEMG_gcL = NormCycle(T.time,Event.ths_l,EMG_norm,...
                    treshold_dtStride,true);
                [~,~,ntrials] = size(lowEMG_gcL);
                lowEMG_gcR = NormCycle(T.time,Event.ths_r,EMG_norm,...
                    treshold_dtStride,true);
                [~,~,ntrials] = size(lowEMG_gcR);

                % get foot placement information (and Xcom for reviewer 1)
                Pelvis = [T.COMPelvis_x T.COMPelvis_y T.COMPelvis_z]./1000;
                FootL = [T.LTOE_x T.LTOE_y T.LTOE_z]./1000;
                FootR = [T.RTOE_x T.RTOE_y T.RTOE_z]./1000;
                StabInfo = getFootPlacementInfo(Event,T.time,Pelvis,FootL,FootR);

                % exoskeleton torque during adaptation
                tau_r = T.RA_JointTorque_Nm;
                tPert_Qual = Event.ths_r;
                [TauExoAv,~] = getAvergePertRStance(Event,tau_r,T.time,...
                    Event.ths_r(2:end-1));

                % exoskeleton work during perturbed Rstance phase «reviewer 1»
                tau_r = T.RA_JointTorque_Nm;
                time = T.time;
                q = T.RA_JointAngle_rad;
                qd = [diff(q)./diff(time); 0];
                ExoPower = tau_r.*qd;
                [PowerExoAv, ~, dtRstance] = getAvergePertRStance(Event,...
                    ExoPower,time, Event.ths_r(2:end-1));
                WorkExoAv = PowerExoAv.* dtRstance;

                % save results
                ntrials = length(Event.ths_l)-1;
                ntrialsR = length(Event.ths_r)-1;
                iSel = 1:ntrials;
                if f<4
                    % muscle information
                    NMCEMGStore(:,s,f,5:7,1:ntrials) =  lowEMG_gcL(:,5:7,:);
                    NMCEMGStore(:,s,f,1:4,1:ntrialsR) =  lowEMG_gcR(:,1:4,:);
                    % average torques
                    ExoAdaptation_NMC(s,f,1:length(TauExoAv),1) = TauExoAv;
                    ExoAdaptation_NMC(s,f,1:length( WorkExoAv),2) =  WorkExoAv;
                    % spatio-temporal information
                    NMCSpatioTempStore(s,f,1,1:ntrials) = Step.StrideTime(1:ntrials);
                    NMCSpatioTempStore(s,f,2,1:ntrials) = Step.dtStance_n(1:ntrials);
                    NMCSpatioTempStore(s,f,3,1:ntrials) = Step.PercDS_n(1:ntrials);
                    NMCSpatioTempStore(s,f,4,1:ntrials) = Step.PercStance_n(1:ntrials);
                    NMCSpatioTempStore(s,f,5,1:ntrials) = Step.PercSwing_n(1:ntrials);
                    NMCSpatioTempStore(s,f,6,1:ntrials) = StabInfo.ths_l.dFeet(1:ntrials,1);
                    NMCSpatioTempStore(s,f,7,1:ntrials) = StabInfo.ths_l.dFeet(1:ntrials,3);
                    NMCSpatioTempStore(s,f,8,1:ntrials) = StabInfo.ths_l.dFootPelvis(1:ntrials,1);
                    NMCSpatioTempStore(s,f,9,1:ntrials) = StabInfo.ths_l.dFootPelvis(1:ntrials,3);
                    NMCSpatioTempStore(s,f,10,1:ntrials) = StabInfo.ths_l.dFootXcom(1:ntrials,1);
                    NMCSpatioTempStore(s,f,11,1:ntrials) = StabInfo.ths_l.dFootXcom(1:ntrials,3);
                else
                    % muscle information
                    ZeroImpEMGStore(:,s,5:7,1:ntrials) =  lowEMG_gcL(:,5:7,:);
                    ZeroImpEMGStore(:,s,1:4,1:ntrialsR) =  lowEMG_gcR(:,1:4,:);
                    % average torques
                    ExoAdaptation_MinImp(s,1:length(TauExoAv),1) = TauExoAv;
                    ExoAdaptation_MinImp(s,1:length( WorkExoAv),2) = WorkExoAv;
                    % spatio-temporal information
                    ZeroImpSpatioTempStore(s,1,1:ntrials) = Step.StrideTime(iSel);
                    ZeroImpSpatioTempStore(s,2,1:ntrials) = Step.dtStance_n(iSel);
                    ZeroImpSpatioTempStore(s,3,1:ntrials) = Step.PercDS_n(iSel);
                    ZeroImpSpatioTempStore(s,4,1:ntrials) = Step.PercStance_n(iSel);
                    ZeroImpSpatioTempStore(s,5,1:ntrials) = Step.PercSwing_n(iSel);
                    ZeroImpSpatioTempStore(s,6,1:ntrials) = StabInfo.ths_l.dFeet(iSel,1);
                    ZeroImpSpatioTempStore(s,7,1:ntrials) = StabInfo.ths_l.dFeet(iSel,3);
                    ZeroImpSpatioTempStore(s,8,1:ntrials) = StabInfo.ths_l.dFootPelvis(iSel,1);
                    ZeroImpSpatioTempStore(s,9,1:ntrials) = StabInfo.ths_l.dFootPelvis(iSel,3);
                    ZeroImpSpatioTempStore(s,10,1:ntrials) = StabInfo.ths_l.dFootXcom(iSel,1);
                    ZeroImpSpatioTempStore(s,11,1:ntrials) = StabInfo.ths_l.dFootXcom(iSel,3);
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
    Info_NMC_EMG.script = 'CreateDataMatrix.m';

    Info_NMC_SpatioTemp.dim1 = SubjPreFix;
    Info_NMC_SpatioTemp.dim2= 'Qual File Number';
    Info_NMC_SpatioTemp.dim4 = 'Trials';
    Info_NMC_SpatioTemp.dim3 = ['(1) = stride duraction, (2)dt_stance, (3) percStance,' ...
        '(4) percSwing, (5) percDS, (6) stride length-x, (7) stride length-z, (8) position foot w.r.t com -x',...
        '(9) position foot w.r.t com -z','(10) position foot w.r.t Xcom -x','(11) position foot w.r.t Xcom -z'];
    Info_NMC_SpatioTemp.script = 'CreateDataMatrix.m';

    Info_ZeroImp_SpatioTemp.dim1 = SubjPreFix;
    Info_ZeroImp_SpatioTemp.dim3 = 'Trials';
    Info_ZeroImp_SpatioTemp.dim2 = ['(1) = stride duraction, (2)dt_stance, (3) percStance,' ...
        '(4) percSwing, (5) percDS, (6) stride length-x, (7) stride length-z, (8) position foot w.r.t com -x',...
        '(9) position foot w.r.t com -z','(10) position foot w.r.t Xcom -x','(11) position foot w.r.t Xcom -z'];
    Info_ZeroImp_SpatioTemp.script = 'CreateDataMatrix.m';

    InfoZeroImp.dim1 = 'norm gait cycle';
    InfoZeroImp.dim2 = SubjPreFix;
    InfoZeroImp.dim3 = {'Soleus R','Gastroc R','Tibialis R','Vastus R','Soleus L','Gastroc L','Tibialis L'};
    InfoZeroImp.dim4 = 'Trials';
    InfoZeroImp.script = 'CreateDataMatrix.m';

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
        'SwingTime','StrideLength_x','StrideLength_z','dPelvisFoot_x','dPelvisFoot_z','COMstance','avExoTorque','av_AnkleAngleR', ...
        'Soleus R LStance','Gastroc R LStance','Tibialis R LStance','Vastus R LStance',...
        'Soleus L LStance','Gastroc L LStance','Tibialis L LStance', 'Foot_Xcom_x', 'Foot_Xcom_z', ...
        'ExoWork R Rstance','WorkPusher','COM_hs1','COM_hs2','COM_hs3','COM_hs4'};
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
                        [Event,Step] = GetSpatioTemporalParam(T.time,T.Fy_L,...
                            T.Fy_R,treshold,dtOffPlate);

                        % movement pelvis on treadmill
                        Pelvis = [T.COMPelvis_x  T.COMPelvis_y T.COMPelvis_z]./1000;
                        [PelvisPos_SingleStance, dPelvis_nhs]= getPelvisMovementTreadmill(Event,T.time,Pelvis,tPert);

                        % get foot placement at first event after perturbation
                        FootL = [T.LTOE_x T.LTOE_y T.LTOE_z]./1000;
                        FootR = [T.RTOE_x T.RTOE_y T.RTOE_z]./1000;
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

                        % exoskeleton work during perturbed Rstance phase «reviewer 1»
                        tau_r = T.RA_JointTorque_Nm;
                        time = T.time;
                        q = T.RA_JointAngle_rad;
                        qd = [diff(q)./diff(time); 0];
                        ExoPower = tau_r.*qd;
                        [PowerExoAv, ~, dtRstance] = getAvergePertRStance(Event,...
                            ExoPower,time, tPert);
                        WorkExoAv = PowerExoAv.* dtRstance;

                        % work done by pusher
                        r_pusher = T.SACR_x/1000 ;
                        r_pusher_dot = [diff(r_pusher)./diff(T.time); 0] + 2/3.6;
                        F_pusher_mag = mass(s)*9.81*0.12;
                        F_pusher = zeros(size(r_pusher));
                        i_pert = find(T.time>=0 & T.time<=0.2);
                        if PertDir == 1 
                            F_pusher(i_pert) = F_pusher_mag;
                        else
                            F_pusher(i_pert) = -F_pusher_mag;
                        end
                        Power_pusher = F_pusher.*r_pusher_dot;
                        [Power_pusherAv, ~, dtRstance] = getAvergePertRStance(Event,...
                            Power_pusher,time, tPert);
                        WorkPusher = Power_pusherAv.* dtRstance;

                        
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
                        iStrideLength = find(strcmp(headers,'StrideLength_x'));
                        iStrideLength_z = find(strcmp(headers,'StrideLength_z'));
                        if isstruct(StabInfo) && length(StabInfo.ths_l.dFeet)>1
                            data(ctr,iStrideLength) = StabInfo.ths_l.dFeet(1);
                             data(ctr,iStrideLength_z) = StabInfo.ths_l.dFeet(3);
                        else
                            data(ctr,iStrideLength) = NaN;
                            data(ctr,iStrideLength_z) = NaN;
                        end

                        % position foot w.r.t. pelvis
                        iSelx = find(strcmp(headers,'dPelvisFoot_x'));
                        iSelz = find(strcmp(headers,'dPelvisFoot_z'));
                        if isstruct(StabInfo) && length(StabInfo.ths_l.dFootPelvis)>1
                            data(ctr,iSelx) = StabInfo.ths_l.dFootPelvis(1);
                            data(ctr,iSelz) = StabInfo.ths_l.dFootPelvis(3);
                        else
                            data(ctr,iSelx) = NaN;
                            data(ctr,iSelz) = NaN;
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

                        % position foot w.r.t. Xcom
                        if isstruct(StabInfo) && length(StabInfo.ths_l.dFootPelvis)>1
                            iSel = find(strcmp(headers,'Foot_Xcom_x'));
                            data(ctr,iSel) = StabInfo.ths_l.dFootXcom(1);
                            iSel = find(strcmp(headers,'Foot_Xcom_z'));
                            data(ctr,iSel) = StabInfo.ths_l.dFootXcom(3);
                        end

                        % exoskeleton work done
                        iSel = find(strcmp(headers,'ExoWork R Rstance'));
                        data(ctr,iSel) = WorkExoAv;   

                        % work done by pusher
                        iSel = find(strcmp(headers,'WorkPusher'));
                        data(ctr,iSel) = WorkPusher;  

                        % question reviewer: COM motion at heelstrikes afte
                        % r perturbations
                        for id = 1:4
                            if length(dPelvis_nhs(:,1))>=id
                                iSel = find(strcmp(headers,['COM_hs' num2str(id)]));
                                data(ctr,iSel) = dPelvis_nhs(id,1);  
                            end
                        end
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