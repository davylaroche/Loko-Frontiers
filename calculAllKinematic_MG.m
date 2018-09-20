

function [kinematic] = calculAllKinematic_MG(C3Dfolder,listec3d)

%acq = btkReadAcquisition(char(strcat(C3Dfolder,listec3d(1))));

pointsLabel = {'Ptilt', 'Pobli','Prota',...
    'Hflex','Habd','Hrota',...
    'Kflex','Kabd','Krota',...
    'Aflex','Frota','Fvert',...
    'PELA','PELO','PELL','PELP',...
    'FEA','FEO','FEL','FEP',...
    'TIA','TIO','TIL','TIP',...
    'TOA','TOO','TOL','TOP',...
    'PlanCovIndex','PlanCov_u1', 'PlanCov_u3'};

%% Initialisation
num = length(pointsLabel);
for k = 1:num
    kinematic.Left.(pointsLabel{k}) = [];
    kinematic.Right.(pointsLabel{k}) = [];
end

%% Calcul Cinématique
for c = 1:length(listec3d)
    fileName = char(strcat(C3Dfolder,listec3d(c)));
    
    [meanAngleLeft] = moyenneAllCycleKinematic_MG(fileName,'Left');
    [meanAngleRight] = moyenneAllCycleKinematic_MG(fileName,'Right');
    
    for k = 1:num
        kinematic.Left.(pointsLabel{k}) = [kinematic.Left.(pointsLabel{k}),meanAngleLeft.(pointsLabel{k})];
        kinematic.Right.(pointsLabel{k}) = [kinematic.Right.(pointsLabel{k}),meanAngleRight.(pointsLabel{k})];
    end 
end

                    
                    
                    