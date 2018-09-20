% FONCTION c3d2ParamST
% Fonction retournant pour 1 fichier c3d l'ensemble des valeurs des 
% paramètres spatio temporel non moyenné (1 valeur pour chaque cycle) de
% chaque coté
% Entrées : - acq : fichier c3d transformé par btkReadAcquisition(nomfichierC3d);
% Sortie : paramST : structure contenant les paramètres spatio temporels à
% gauche et à droite pour chaque cycle (ex: paramST.Left.SpeedCoM = [1.1,1.2,1.15,..,1.01]

% Alexandre Naaim : 19/09/2017 (Mise au propre) 

function [paramST] = c3d2ParamST_MG(acq)
%% Transformation du referentiel temps :
% Dans les c3d les événements peuvent être mis par rapport au moment ou
% l'acquisition a débuté et non le moment ou elle a été coupé. Alors que
% les données type point (ex position de marqueurs etc..) sont référencé
% par rapport à l'endroit ou le fichier est coupé. Grace à cette partie ou
% recoordonne les deux pour pouvoir faire les calculs.
[Events ~] = btkGetEvents(acq);
FF = btkGetFirstFrame(acq);
freq = btkGetPointFrequency(acq);
EventType = fieldnames(Events);
Point = btkGetPoints(acq);

% Event are integrated in the number of frame of the acqusition:
for i = 1:size(EventType,1)
    Events.(EventType{i}) = round(Events.(EventType{i})*freq)-FF+1;
end

%% Initialisation
% Left
% distance and time
paramST.Left.stepLenght = [];
paramST.Left.stepTime = [];

paramST.Left.strideLenght = [];
paramST.Left.strideWidth = [];
paramST.Left.strideTime = [];

paramST.Left.speedCoM = [];
paramST.Left.cadence = [];
paramST.Left.speedCadenceStrideLenght = [];

%percentage
paramST.Left.stancePhase = [];
paramST.Left.swingPhase = [];
paramST.Left.totalDoubleSupport = [];
paramST.Left.totalSingleSupport = [];
paramST.Left.DS1 = [];
paramST.Left.DS2 = [];

paramST.Left.SingleSupportTime = [];
paramST.Left.DoubleSupportTime = [];
paramST.Left.stanceTime = [];
paramST.Left.swingTime = [];

% Similar parameter on the left and the right
paramST.Right = paramST.Left;


%% Analyse de chaque coté
% Définition des coté et des coté controlatéral associées
cote = {'Right','Left'};
coteAbr = {'R','L'};
controLateral = {'Left','Right'};
controLateralAbr = {'L','R'};

% Pour chaque coté :
for i = 1:2
    
    % Pour chaque pas on calcul
    
    % Foot Strike
    FSstr = strcat(cote(i),'_Foot_Strike');
    % Cell2Str before we have a cell with strcat using str{} allow to get a
    % string that can be used in dynamic structure calling
    FS = Events.(FSstr{1});
    
    % Controlateral FS
    C_FSstr = strcat(controLateral(i),'_Foot_Strike');
    C_FS = Events.(C_FSstr{1});
    
    % Foot off
    FOstr = strcat(cote(i),'_Foot_Off');
    FO = Events.(FOstr{1});
    
    % Controlateral Foot strike
    C_FOstr =  strcat(controLateral(i),'_Foot_Off');
    C_FO = Events.(C_FOstr{1});
    
    % Point used
    Heestr =  strcat(coteAbr(i),'HEE');
    C_Heestr = strcat(controLateralAbr(i),'HEE');
    
    Hee = Point.(Heestr{1});
    C_Hee = Point.(C_Heestr{1});
    % COM barycent of LASI RASI LPSI RPSI
    CoM = Point.LASI+Point.RASI+Point.LPSI+Point.RPSI;
    CoM = CoM/4;
    
    % Calcul de la direction dans laquelle le patient ce déplace
    [~,indDeplMax] = max(abs(CoM(end,:)-CoM(1,:)));
    
    % We remove from all the events the one before the first ispsolateral
    % foot strike    
    FO = FO(FO>FS(1));
    C_FS = C_FS(C_FS>FS(1));
    C_FO = C_FO(C_FO>FS(1));
    
    for c1 = 1:(size(FS,2)-1)
        % We need to test that between 2 foot strike there is a foot off and
        % a controlateral foot strike and off
        testFO = isempty(FO(FO>FS(c1) & FO<FS(c1+1)));
        testC_FO = isempty(C_FO(C_FO>FS(c1) & C_FO<FS(c1+1)));
        testC_FS = isempty(C_FS(C_FS>FS(c1) & C_FS<FS(c1+1)));
        % If one of this condition is not respected we add at the beginning
        % of the Event vector a 0 value in order to keep the correspondence
        % between the foot strike event and 
        if testFO
            FO = [0,FO];
        end
        if testC_FO
            C_FO = [0,C_FO];
        end
        if testC_FS
            C_FS = [0,C_FS];
        end
        % If no test is empty (which mean that there is the event between 
        % the 2 foot strike) it is possible to compute all the different
        % spation temporal parameter.
        
        if not(testFO | testC_FO | testC_FS)
            %% Step lenght and time
            %   Lenght m
            Steptemp = Hee(FS(c1+1),:) - C_Hee(C_FS(c1),:);
            StepLenght = abs(Steptemp(indDeplMax)/1000);
            
            %   Time = swing phase?(we remove one frame because we consider that
            %   the cycle is finished before the foot strike
            StepTime = (FS(c1+1)-C_FS(c1))/freq;
            
            paramST.(cote{i}).stepLenght = [paramST.(cote{i}).stepLenght,StepLenght];
            paramST.(cote{i}).stepTime = [paramST.(cote{i}).stepTime,StepTime];
            
            
            %% Stride lenght and time
            %   Lenght
            stridetemp = Hee(FS(c1+1),:) - Hee(FS(c1),:);
            strideLenght = abs(stridetemp(indDeplMax)/1000);
            
            % Distance between a droite and a point
            % https://fr.wikipedia.org/wiki/Distance_d%27un_point_%C3%A0_une_droite
            u = Hee(FS(c1+1),:) - Hee(FS(c1),:);
            BA = C_Hee(C_FS(c1),:) - Hee(FS(c1),:);
            strideWidth = norm(cross(BA,u))/norm(u)/1000;
            
            % Time
            strideTime = (FS(c1+1)-FS(c1))/freq;
            
            paramST.(cote{i}).strideLenght = [paramST.(cote{i}).strideLenght,strideLenght];
            paramST.(cote{i}).strideWidth = [paramST.(cote{i}).strideWidth,strideWidth];
            paramST.(cote{i}).strideTime = [paramST.(cote{i}).strideTime,strideTime];
            
            %% Cadence (step/min)
            cadence = 2*60/strideTime;
            
            % Vitesse : 2 méthodes
            
            % Méthode : cadence*Stride_lenght/120
            speed = cadence*strideLenght/120;
            
            % Speed CoM : The horizontal difference between the middle of the pelvis
            % markers at the beginning and the end of the acqusition divided by the
            % time of the acqusition:
            CoMtemp = CoM(FS(c1+1),:) - CoM(FS(c1),:);
            CoMtemp = abs(CoMtemp(indDeplMax)/1000);
            speedCoM = CoMtemp/strideTime;
            
            paramST.(cote{i}).speedCoM = [paramST.(cote{i}).speedCoM,speedCoM];
            paramST.(cote{i}).cadence = [paramST.(cote{i}).cadence,cadence];
            paramST.(cote{i}).speedCadenceStrideLenght = [paramST.(cote{i}).speedCadenceStrideLenght,speed];
            
            
            %% Paramètre de cycle
            strideTot = FS(c1+1)-FS(c1);
            
            % End and Begin of the first and second double support
            firstDoubleSupport = 100*(C_FO(c1)-FS(c1))/strideTot;
            secondeDoubleSupport = 100*(C_FS(c1)-FS(c1))/strideTot;
                        
            % Pourcentage/time Stance et Swing Phase
            stanceFrame = FO(c1)-FS(c1);
            swingFrame = FS(c1+1)-FO(c1);
            
            stancePerc = 100*stanceFrame/strideTot;
            stanceTime = stanceFrame/freq;
            swingPerc = 100*swingFrame/strideTot;
            swingTime = swingFrame/freq;         
            
            % Percentage end 1st and 2nd double Support :
            firstDoubleSupportNbrFrame = C_FO(c1)-FS(c1);
            secondeDoubleSupportNbrFrame = FO(c1)-C_FS(c1);
                       
            firstDoubleSupportPerc = 100*firstDoubleSupportNbrFrame/strideTot;
            secondeDoubleSupportPerc = 100*secondeDoubleSupportNbrFrame/strideTot;
            
            doubleSupportPerc = firstDoubleSupportPerc+secondeDoubleSupportPerc;
            doubleSupportTime = (firstDoubleSupportNbrFrame+secondeDoubleSupportNbrFrame)/freq;
            
            % Total single and total double support :
            singleSupportNbrFrame = C_FS(c1)-C_FO(c1);
            singleSupportPerc = 100*singleSupportNbrFrame/strideTot;
            singleSupportTime = singleSupportNbrFrame/freq;
            
            %singleSupport = stancePerc-doubleSupport;
            
            paramST.(cote{i}).stancePhase = [paramST.(cote{i}).stancePhase,stancePerc];
            paramST.(cote{i}).swingPhase = [paramST.(cote{i}).swingPhase,swingPerc];
            
            paramST.(cote{i}).DS1 = [paramST.(cote{i}).DS1,firstDoubleSupport];
            paramST.(cote{i}).DS2 = [paramST.(cote{i}).DS2,secondeDoubleSupport];
            
            paramST.(cote{i}).totalDoubleSupport = [paramST.(cote{i}).totalDoubleSupport,doubleSupportPerc];
            paramST.(cote{i}).totalSingleSupport = [paramST.(cote{i}).totalSingleSupport,singleSupportPerc];
            
            paramST.(cote{i}).SingleSupportTime = [paramST.(cote{i}).SingleSupportTime,singleSupportTime];
            paramST.(cote{i}).DoubleSupportTime = [paramST.(cote{i}).DoubleSupportTime,doubleSupportTime];
            paramST.(cote{i}).stanceTime = [paramST.(cote{i}).stanceTime,stanceTime];
            paramST.(cote{i}).swingTime = [paramST.(cote{i}).swingTime,swingTime];
            
        end       
    end
end