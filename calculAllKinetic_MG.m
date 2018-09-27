% FONCTION calculAllParamST
% Fonction retournant une structure contanant pour le coté Gauche et Droit
% la cinématique et cinétique calculer sur les fichiers c3d données
% en Entrées
% Entrées : - C3Dfolder : chemin des fichiers .c3d
%           - listeC3D liste des fichiers de posture
% Sortie : kinetic : structure contenant la cinématique moyenne, paramètres spatio temporels à
% gauche et à droite


function [kinetic] = calculAllKinetic_MG(C3Dfolder,listec3d)

% Définition des valeurs que l'on souhaite extraire 
pointsLabel = {'Ptilt', 'Pobli','Prota',...
    'Hflex','Habd','Hrota',...
    'Kflex','Kabd','Krota',...
    'Aflex','Frota','Fvert',...
    'PELA','PELO','PELL','PELP',...
    'FEA','FEO','FEL','FEP',...
    'TIA','TIO','TIL','TIP',...
    'TOA','TOO','TOL','TOP',...
    'PlanCovIndex','PlanCov_u1', 'PlanCov_u3', ...
    'HipExtMoment','HipAbdMoment','HipRotMoment',...
    'KneeExtMoment','KneeAbdMoment','KneeRotMoment',...
    'AnklePlantMoment',...
    'HipPower','KneePower','AnklePower',...
    'LatMedGRF','AntPostGRF','VerticalGRF','CoM_vert','SpeedWalkingCoM'};

%% Initialisation
num = length(pointsLabel);
for k = 1:num
    kinetic.Left.(pointsLabel{k}) = [];
    kinetic.Right.(pointsLabel{k}) = [];
end

%% Calcul Cinématique
for c = 1:length(listec3d)
    fileName = char(strcat(C3Dfolder,listec3d(c)));
    [meanAngleLeft] = moyenneAllCycleKinetic_MG(fileName,'Left',C3Dfolder);
    [meanAngleRight] = moyenneAllCycleKinetic_MG(fileName,'Right',C3Dfolder);
    for k = 1:num
        kinetic.Left.(pointsLabel{k}) = [kinetic.Left.(pointsLabel{k}),meanAngleLeft.(pointsLabel{k})];
        kinetic.Right.(pointsLabel{k}) = [kinetic.Right.(pointsLabel{k}),meanAngleRight.(pointsLabel{k})];
    end 

end