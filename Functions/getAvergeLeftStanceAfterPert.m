function [datmean] = getAvergeLeftStanceAfterPert(Event,dat,t,tPert)
%getAvergeLeftStanceAfterPer Computes the average value (of input argument dat)
% during the first left stance phase after perturbation onset
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) data = vector with data
%       (3) t = time vector
%       (4) tPert = perturnation onset


% we want to compute the average exoskeleton moment from right heelstrike
% before perturbation onset until right toe-off
ths = Event.ths_l;

% only first right heelstrike before perturbation onset

thsO = ths;
ths = nan(1,length(tPert));
for ip = 1:length(tPert)
    % first first event after preturbation onset
    ths(ip) = thsO(find(thsO>tPert(ip),1,'first'));
end


datmean= nan(length(ths),1);
for i=1:length(ths)
    % find next lef theelstrike
    t0 = ths(i);
    tend = Event.tto_l(find(Event.tto_l>ths(i),1,'first'));
    if ~isempty(tend)
        % indices between events
        iSel = t>=t0 & t<=tend;
        datmean(i) = nanmean(dat(iSel));
    end
end