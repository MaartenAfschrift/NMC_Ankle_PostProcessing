function [Subj, SubjFolders, SubjPreFix, SubjID_Exo, mass, height, ...
    age, sID ] = GetSubjInfo(SubjFileDatPath)
%GetSubjInfo Reads the .yml files of all subjects

ct = 1;
for s = 1:30
    % yml file
    ymlFile = fullfile(SubjFileDatPath, ['subject' num2str(s) '.yml']);
    % read the .yml file
    if exist(ymlFile,'file')
        Subject = yaml.loadFile(ymlFile);
        
        % add information
        Subj{ct} = Subject;

        % path information
        SubjFolders{ct} = char(Subject.subject.id_folder);
        SubjPreFix{ct} = char(Subject.subject.id_Qualisys);
        temp = char(Subject.subject.id_exo);
        SubjID_Exo{ct} = temp(1:end-1);

        % subject statistics
        mass(ct) = Subject.subject.mass;
        height(ct) = Subject.subject.height;
        age(ct) = Subject.subject.age;

        % subject ID
        sID(ct) = Subject.subject.id;
        ct = ct + 1;
    end
end