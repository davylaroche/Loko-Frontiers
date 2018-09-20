% Permet de lister tous les .c3d de tous les patients
% Patients concern�s : inscrits dans le fichier excel
% Protocole : Locox 2
clc; clear all; 

FILENAME_DBPATIENTS = 'D:\projets termin�s\LOCOX_2\liste_volontaires.xlsx';
folder = 'D:\projets termin�s\LOCOX_2\Volontaires\';

[~,~,DB] = xlsread(FILENAME_DBPATIENTS);

id = char(DB(1:end,1));
    
%% Listing C3D par Volontaires
for c1 = 1 : length(id)
    c1
    id(c1, :)
    
    C3Dfolder = [folder,id(c1, :),'\Marche\'];
    listec3d = dir(fullfile([C3Dfolder,'*.c3d']));
    
    file_list={listec3d.name};
  
    %% Param�tres ST
    paramST = calculAllParamST_MG(C3Dfolder,file_list);
    
    %% Cin�matique
    [kinematic] = calculAllKinematic_MG(C3Dfolder,file_list);
    
    %% Cin�tique
    [kinetic] = calculAllKinetic_MG(C3Dfolder,file_list);
       
end

