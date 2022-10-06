function [datAv,datNorm,dt_Rstance] = getAvergePertRStance(Event,dat,t,tPert)
%getAvergePertRStance Computes the average value of input argument dat during Rstance phase
%   % input arguments:
%       (1) Event = structure with gait cycle events
%       (2) dat = input data
%       (3) t = time vector
%       (4) tPert = timing perturbation


% we want to compute the average exoskeleton moment from right heelstrike
% before perturbation onset until right toe-off
ths_r = Event.ths_r;

% only first right heelstrike before perturbation onset
ths_rO = ths_r;
ths_r = nan(1,length(tPert));
for ip = 1:length(tPert)
    % first first event after preturbation onset
    ths_r(ip) = ths_rO(find(ths_rO<tPert(ip),1,'last'));
end


datAv = nan(length(ths_r),1);
datNorm = nan(100,length(ths_r));
dt_Rstance = nan(length(ths_r),1);
for i=1:length(ths_r)
    % find next left theelstrike
    t0 = ths_r(i);
    tend = Event.tto_r(find(Event.tto_r>ths_r(i),1,'first'));
    if ~isempty(tend)
        % indices between events
        iSel = t>=t0 & t<=tend;
        dsel = dat(iSel);
        datAv(i) = nanmean(dsel);
        datNorm(:,i) = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
        dt_Rstance(i) = tend-t0;
    end
end