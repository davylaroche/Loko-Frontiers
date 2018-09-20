function cycles_cut(filename)
% Permet de découper les cycles d'un fichier .c3d
% Paramètre - filename : chemin complet du fichier

filename

H = btkReadAcquisition(filename); % Lecture du fichier
frequ = btkGetPointFrequency(H); % Fréquence d'acquisition
FF = btkGetFirstFrame(H); % 1ère frame
data = btkGetMarkers(H); % Position de tous les marqueurs au cours du temps
btkClearEvents(H); % Supprime les événements du fichier

% Détermination des paramètres des événements à enregistrer
FootStrike = {1, 'Foot Strike',  'The instant the heel strikes the ground'};
FootOff =  {2, 'Foot Off',  'The instant the toe leaves the ground'};

for c5 = 1:2 % Pour différencier côté droit/gauche
    clear Xheel XTO Xsacrum tHS tTO HS TO TOFF METADATA
    if c5 == 1
        side = 'Left';
        Yheel = data.LHEE(:,2); %Position en Y du marqueur du talon            
    else
        side = 'Right';
        Yheel = data.RHEE(:,2); 
    end

    % Calcul de la position en Y du sacrum
    Ysacrum = (data.RPSI(:,2) + data.LPSI(:,2))/2;
    
    % Calcul de l'écart entre le sacrum et le talon au cours du temps
    % Suivant le sens de marche (+Y ou -Y)
    if Ysacrum(end, 1)-Ysacrum(1, 1)>0
        tHS = -Ysacrum + Yheel;
    else
        tHS = Ysacrum - Yheel;
    end
    % Sélection du maximum de cet écart pour déterminer les événements
    [HS, ~, TO]= Peaks_Select(tHS, 'semi-automatic');
    % Recalage à 0
    HS = HS + FF; % heel strike
    TOFF = HS(1:end-1) + TO; % toe off

    %%% Inscription des événements
    METADATA = btkGetMetaData(H, 'SUBJECTS', 'NAMES');
    for c3 = 1:length(HS)
         btkAppendEvent(H, FootStrike{2}, HS(c3)/frequ, side, char(METADATA.info.values), FootStrike{3}, FootStrike{1});
    end
    for c4 = 1:length(TOFF)
         btkAppendEvent(H, FootOff{2}, TOFF(c4)/frequ, side, char(METADATA.info.values), FootOff{3}, FootOff{1});
    end
end

% Enregistrement du .c3d
btkWriteAcquisition(H, filename);

end