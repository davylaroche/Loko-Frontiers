% Extraction à partir des .enf des informations plateforme
% modification du 12/09/2017 : integration du cas special des nouveaux
% fichiers enf 
function [plateformDetect] = plateformDetection(filename,side)

fileID = fopen(filename,'r');
fileCharac = fscanf(fileID,'%s',Inf);
fclose(fileID);

% on teste si la chaine de charactère "ECL_" est dans le filename
if not(isempty(strfind(fileCharac,'ECL_')))
    
    % On detecte la position des endroits ou la chaine de charactere FP
    % apparait (c'est la ou il y a les informations plateformes)
    FPposition = strfind(fileCharac,'FP');
    
    % La premiere occurence de FP decrit le numéro de plateforme la deuxième si
    % le pied gauche/droit/invalide a été dessus. On parcour l'ensemble des
    % valeur de FPposition de 2 en 2 et on enregistre la premiere comme le
    % numero de plateforme et la deuxième comme le pied
    
    plateformDetect = [];
    for i = 1:2:size(FPposition,2)
        if strcmpi(side,'right')&& strcmpi(fileCharac(FPposition(i+1)+2),'r')
            plateformDetect = [plateformDetect,str2num(fileCharac(FPposition(i)+2))];
            
        elseif strcmpi(side,'left')&& strcmpi(fileCharac(FPposition(i+1)+2),'l')
            plateformDetect = [plateformDetect,str2num(fileCharac(FPposition(i)+2))];
        end
    end
    
else
    
    FPposition = strfind(fileCharac,'FP');
    plateformDetect = [];
    for i = 1:size(FPposition,2)
        if strcmpi(side,'right')&& strcmpi(fileCharac(FPposition(i)+4),'r')
            
            plateformDetect = [plateformDetect,str2num(fileCharac(FPposition(i)+2))];
            
        elseif strcmpi(side,'left')&& strcmpi(fileCharac(FPposition(i)+4),'l')
            
            plateformDetect = [plateformDetect,str2num(fileCharac(FPposition(i)+2))];
        
        end
    end
end