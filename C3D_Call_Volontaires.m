% Permet de lister tous les .c3d de tous les patients
% Patients concernés : inscrits dans le fichier excel
% Protocole : Locox 2
clc; clear all; 

FILENAME_DBPATIENTS = 'D:\projets terminés\LOCOX_2\liste_volontaires.xlsx';
folder = 'D:\projets terminés\LOCOX_2\Volontaires\';

[~,~,DB] = xlsread(FILENAME_DBPATIENTS);

id = char(DB(1:end,1));
    
%% Listing C3D par Volontaires
for c1 = 1 : length(id)
    c1
    id(c1, :)
    
    C3Dfolder = [folder,id(c1, :),'\Marche\'];
    listec3d = dir(fullfile([C3Dfolder,'*.c3d']));
    
    file_list={listec3d.name};
  
    %% Paramètres ST
    paramST = calculAllParamST_MG(C3Dfolder,file_list);
    
    %% Cinématique
    [kinematic] = calculAllKinematic_MG(C3Dfolder,file_list);
    
    %% Cinétique
    [kinetic] = calculAllKinetic_MG(C3Dfolder,file_list);
       
end

