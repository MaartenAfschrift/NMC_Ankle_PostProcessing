function [StabInfo] = getFootPlacementInfo(Event,time,Pelvis,FootL,FootR,varargin)
%getFootPlacementInfo Computes distance between pelvis and foot at
% gait events based on raw Marker coordinates
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) time = time vector
%       (3) Pelvis = array with pelvis xyz coordinates
%       (4) FootL = array with FootL coordinates
%       (5) FootR = array with FootR coordinates
%       (6) varargin:
%           (1): timing perturbation
%           (2): treadmill speed
%           (3): LegLength

tPert = [];
if ~isempty(varargin)
    tPert = varargin{1};
end
vTreadmill = 2./3.6;
if length(varargin)>1
    vTreadmill = varargin{2};
end
LegLength = 0.9;
if length(varargin)>2
    LegLength = varargin{3};
end


% compute pelvis velocity
Pelvis_dot = diff(Pelvis)./diff(time);
Pelvis_dot = [Pelvis_dot; zeros(1,3)];

% compute Xcom position
Xcom = zeros(size(Pelvis));
Xcom(:,1) = Pelvis(:,1) + (Pelvis_dot(:,1) + vTreadmill)./(sqrt(9.81/LegLength));
Xcom(:,3) =  Pelvis(:,3) + (Pelvis_dot(:,3))./(sqrt(9.81/LegLength));

% we want to loop over all event names
eventNames = fieldnames(Event);

% Get Marker information
for i=1:length(eventNames)
    % get timing of events
    tEvent = Event.(eventNames{i});
    if ~isempty(tPert)
        tEventO = tEvent;
        tEvent = nan(1,length(tPert));
        for ip = 1:length(tPert)
            % first first event after preturbation onset
            tEvent(ip) = tEventO(find(tEventO>tPert(ip),1,'first'));
        end
    end
    iEvent = nan(size(tEvent));
    for ie=1:length(iEvent)
        iEvent(ie) = find(time>=tEvent(ie),1,'first');
    end

    % compute distance between marker coordinates
    if length(FootR)>1 && length(FootL)>1 && length(Pelvis)>1
        if strcmp(eventNames{i}(end-1:end),'_l')
            dFootPelvis = FootL(iEvent,:)-Pelvis(iEvent,:);
            dFeet = FootL(iEvent,:)-FootR(iEvent,:);
            dFoot_Xcom = FootL(iEvent,:)-Xcom(iEvent,:);
        else
            dFootPelvis = FootR(iEvent,:)-Pelvis(iEvent,:);
            dFeet = FootR(iEvent,:)-FootL(iEvent,:);
            dFoot_Xcom = FootR(iEvent,:)-Xcom(iEvent,:);
        end
    else
        dFootPelvis = NaN;
        dFeet = NaN;
        dFoot_Xcom = NaN;
    end
    StabInfo.(eventNames{i}).dFeet = dFeet;
    StabInfo.(eventNames{i}).dFootPelvis = dFootPelvis;
    StabInfo.(eventNames{i}).dFootXcom = dFoot_Xcom;
end
end