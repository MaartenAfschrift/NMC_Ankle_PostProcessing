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
    disp('Hit enter if you want to download Afschrift 2022 data (500 MB)');
    pause();
    disp(' download started')
    LinkZipFile = 'https://www.dropbox.com/s/pwe0hkrwexhpqej/Data_NMC_Balance_Ankle.zip?dl=1';    
    websave(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'),LinkZipFile);
    disp(' unzip DataFile');
    unzip(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'),DataPath);
    delete(fullfile(DataPath,'Data_NMC_Balance_Ankle.zip'));
    disp(' Downloading Afschrift 2022 data finished');
end

end