
% FONCTION calculAllParamST
% Fonction retournant une structure contanant pour le coté Gauche et Droit
% l'ensemble des paramètres spatio temporel calculer sur les fichiers c3d
% données en Entrées
% Entrées : - C3Dfolder : chemin des fichiers .c3d
%           - listeC3D liste des fichiers de posture
% Sortie : paramST : structure contenant les paramètres spatio temporels à
% gauche et à droite

function [paramST] = calculAllParamST_MG(C3Dfolder,listec3d)

%% Initialisation de paramST
% lecture du premier fichier par mokka
acq = btkReadAcquisition(char(strcat(C3Dfolder,listec3d(1))));
% Calcul des paramètres spatiotemporel pour ce fichier 
paramST = c3d2ParamST_MG(acq);
% Détermination des paramètres spatio temporel contenu dans paramST
paramSTSideField = fieldnames(paramST);
paramSTName = fieldnames(paramST.(paramSTSideField{1}));

%% Calcul pour chaque fichier
% Vérification qu'il y a plus d'un ficheir
if length(listec3d)>1
    % On commence la boucle for a 2 car on a deja initialisé avec 1
    for c = 2:length(listec3d)
        %Lecture fichier c3d et calcul des paramètres spatioTemp
        acq = btkReadAcquisition(char(strcat(C3Dfolder,listec3d(c))));
        paramSTtemp = c3d2ParamST_MG(acq);
        % Pour chaque coté ...
        for c1 = 1:size(paramSTSideField,1)
            % ... et chaque paramètres spatio temporel...
            for c2 = 1:size(paramSTName,1)
                % ...on enregistre dans paramST global 
                paramST.(paramSTSideField{c1}).(paramSTName{c2}) = [...
                    paramST.(paramSTSideField{c1}).(paramSTName{c2}),...
                    paramSTtemp.(paramSTSideField{c1}).(paramSTName{c2})];
            end
        end
    end
end
