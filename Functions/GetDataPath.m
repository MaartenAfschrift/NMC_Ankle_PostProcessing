function [Datapath] = GetDataPath()
%GetDataPath selects the data folder (depending on the computer)
%   Detailed explanation goes here

% get the name of the current computer
CompName = getenv('COMPUTERNAME');

% datapath on desktop @Home
if strcmp(CompName,'DESKTOP-D77FCSJ')
    Datapath = 'D:\Twente_AnkleExo';
elseif strcmp(CompName,'VU-MWP-28313857')
    Datapath = 'D:\Data\Twente_AnkleExo';
else
    disp('Warning: No datapath specified for this computer');
    disp('Using the default datapath');
    Datapath = 'D:\Twente_AnkleExo';
end


end