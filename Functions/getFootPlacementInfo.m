function [StabInfo] = getFootPlacementInfo(Event,time,Pelvis,FootL,FootR,varargin)
%getFootPlacementInfo Computes distance between pelvis and foot at
% gait events based on raw Marker coordinates
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) Qualdat = data structure from qualisys
%       (3) tPert = only report outcomes at first event after tPert

tPert = [];
if ~isempty(varargin)
    tPert = varargin{1};
end

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
        else
            dFootPelvis = FootR(iEvent,:)-Pelvis(iEvent,:);
            dFeet = FootR(iEvent,:)-FootL(iEvent,:);
        end
    else
        dFootPelvis = NaN;
        dFeet = NaN;
    end
    StabInfo.(eventNames{i}).dFeet = dFeet;
    StabInfo.(eventNames{i}).dFootPelvis = dFootPelvis;
end
end