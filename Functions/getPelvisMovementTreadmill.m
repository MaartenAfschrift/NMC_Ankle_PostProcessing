function [dPelvis, dPelvis_nhs] = getPelvisMovementTreadmill(Event,time,Pelvis,varargin)
%getPelvisMovementTreadmill Computes the displacement of the pelvis between
%right heelstrike and the subsequent right heelstrike
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) Qualdat = data structure from qualisys
%       (3) tPert = only report outcomes at first event after tPert

tPert = [];
if ~isempty(varargin)
    tPert = varargin{1};
end


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
        i0 = find(time>=t0,1,'first');
        iend = find(time>=tend,1,'first');
        % compute distance between marker coordinates
        dPelvis(i,:) = Pelvis(iend,:) - Pelvis(i0,:);
    end
end

% addition based on question reviewer 1: also compute pelvis position at
% first n-heelstrikes after perturbation onset:
ths = [Event.ths_r; Event.ths_l]; % vector with all heelstrikes
ths = sort(ths);
ihs = find(ths>0); % all heelstrikes after perturbation onset
dPelvis_nhs = nan(length(ihs),3);
for i=1:length(ihs)
    % get pelvis position before perturbation onset
    t0 = ths_r;
    i0 = find(time>=t0,1,'first');
    x0 = Pelvis(i0,:);
    % get pelvis position at ith heelstrike
    tend = ths(ihs(i));
    iend = find(time>=tend,1,'first');
    xs = Pelvis(iend,:);
    % store pelvis displacement
    dPelvis_nhs(i,:) = xs-x0;
end





