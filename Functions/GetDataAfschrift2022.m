function [DataPath] = GetDataAfschrift2022()
% GetDataAfschrift2022 Downloads the data of steady-state and perturbed
% walking with ankle-foot exoskeleton controlled with NMC controller (with
% or without COM feedback)

% Publication [under review]
% Bio-inspired control of ankle-foot exoskeleton reduces ankle muscle activity
% in steady-state and perturbed walking

% get path to DataFolder
[FuncPath,~,~] = fileparts(mfilename('fullpath'));
[MainPath,~] = fileparts(FuncPath);
DataPath = fullfile(MainPath,'Data');

if isfolder(DataPath) && isfile(fullfile(DataPath,'subject_1','MinImp','Pull_1.csv'))
    disp('Data was already downloaded');
else
    if ~isfolder(DataPath)
        mkdir(DataPath);
    end
    % download the executable from google drive
    disp('Hit enter if you want to download Afschrift 2022 data from https://osf.io/ekrqj (500 MB)');
    pause();
    disp(' download started')
    % download from dropbox
    %LinkZipFile = 'https://www.dropbox.com/s/nikwfsov5ch4nb4/Data_NMC_Balance_Ankle.zip?dl=1';
%     LinkZipFile = 'https://www.dropbox.com/s/l7qr5m4ik7ar757/Data_NMC_Balance_Ankle.zip?dl=1';
    % test download from osf (test if this works....)
    LinkZipFile = 'https://osf.io/download/6496e6506c0981071b7d1b7a/';    
    websave(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'),LinkZipFile);
    disp(' unzip DataFile');
    unzip(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'),DataPath);
    delete(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'));
    disp(' Downloading Afschrift 2022 data finished');
    disp([' You can find the data in: ' DataPath]);
    disp(' You might want to read the README.md file in this folder')
end

end
