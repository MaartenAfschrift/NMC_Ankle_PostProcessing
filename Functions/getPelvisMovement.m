function [dPelvis] = getPelvisMovement(Event,Qualdat,varargin)
%getFootPlacement_MarkerBased Computes distance between pelvis and foot at
% gait events based on raw Marker coordinates
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) Qualdat = data structure from qualisys
%       (3) tPert = only report outcomes at first event after tPert

tPert = [];
if ~isempty(varargin)
    tPert = varargin{1};
end

NamePlevisMarker = 'SACR';
if length(varargin)>1
    NamePlevisMarker = varargin{2};
end

% Get Marker information
if isfield(Qualdat,'Trajectories')
    Labels = Qualdat.Trajectories.Labeled.Labels;
    Markers = Qualdat.Trajectories.Labeled.Data;
    [nc,ndim,nframes] = size(Markers);

    if any(strcmp(Labels,NamePlevisMarker))
        Pelvis = squeeze(Markers(strcmp(Labels,NamePlevisMarker),1:3,:))';
        t_Marker = (0:nframes-1)./Qualdat.FrameRate;

        % we want to compute the pelvis movement from right heelstrike untill left
        % heelstrike
        ths_r = Event.ths_r;

        % only first right heelstrike before perturbation onset
        if ~isempty(tPert)
            ths_rO = ths_r;
            ths_r = nan(1,length(tPert));
            for ip = 1:length(tPert)
                % first first event after preturbation onset
                ths_r(ip) = ths_rO(find(ths_rO<tPert(ip),1,'last'));
            end
        end


        dPelvis = nan(length(ths_r),3);
        for i=1:length(ths_r)
            % find next lef theelstrike
            t0 = ths_r(i);
            tend = Event.ths_l(find(Event.ths_l>ths_r(i),1,'first'));
            if ~isempty(tend)
                % get index associated with event
                i0 = find(t_Marker>=t0,1,'first');
                iend = find(t_Marker>=tend,1,'first');
                % compute distance between marker coordinates
                dPelvis(i,:) = Pelvis(iend,:) - Pelvis(i0,:);
            end
        end
    else
        dPelvis = NaN;
    end
else
    dPelvis = NaN;
end