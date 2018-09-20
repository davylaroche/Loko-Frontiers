function [meanAngle] = moyenneAllCycleKinetic_MG(fileName,side,pathname)

pointsLabel = {'RPelvisAngles', 'LPelvisAngles',...
    'RHipAngles', 'LHipAngles',...
    'RKneeAngles', 'LKneeAngles',...
    'RAnkleAngles', 'LAnkleAngles',...
    'LFootProgressAngles','RFootProgressAngles',...
    'RHipMoment','LHipMoment',...
    'RKneeMoment','LKneeMoment',...
    'RAnkleMoment','LAnkleMoment',...
    'RHipPower','LHipPower',...
    'RKneePower','LKneePower',...
    'RAnklePower','LAnklePower',...
    'RNormalisedGRF','LNormalisedGRF'};

num = length(pointsLabel);
indglob = 1;


% pointsLabelStruct = {'Ptilt', 'Pobli','Prota',...
%     'Hflex','Habd','Hrota',...
%     'Kflex','Aflex','Frota'};

pointsLabelStruct = {'Ptilt', 'Pobli','Prota',...
    'Hflex','Habd','Hrota',...
    'Kflex','Kabd','Krota',...
    'Aflex','Frota','Fvert',...
    'HipExtMoment','HipAbdMoment','HipRotMoment',...
    'KneeExtMoment','KneeAbdMoment','KneeRotMoment',...
    'AnklePlantMoment',...
    'HipPower','KneePower','AnklePower',...
    'LatMedGRF','AntPostGRF','VerticalGRF'};


numStruct = length(pointsLabelStruct);

for k = 1:numStruct
    %meanAngleLeft.(pointsLabelStruct{k}) = [];
    meanAngle.(pointsLabelStruct{k}) = [];
end


if strcmp(side,'Right')
    C_side = 'Left';
    forceName = 'RGroundReactionForce';
else
    C_side ='Right';
    forceName = 'LGroundReactionForce';
end


meanAngleRight = [];
meanAngleLeft = [];

acq = btkReadAcquisition(fileName);
points = btkGetPoints(acq);
Events = btkGetEvents(acq);

%% version avec PROCESSING
metaData = btkGetMetaData(acq);
bodyMass = metaData.children.PROCESSING.children.Bodymass.info.values;


%% version sans PROCESSING (utilisation fichier sujet.mp pour obtenir les données anthropométriques)
% % test = cd;
% 
% filefile = dir(fullfile(pathname, '*.mp'));
% 
% if length(filefile)>1
%     cd(pathname)
%     filefile = uigetfile('*.mp');
% end
% 
% filefilefile = [pathname,filefile(1).name];
% 
% fid = fopen(filefilefile, 'rt');
% datacell = textscan( fid, '%s=%f' );
% fclose(fid);
% textdata = datacell{1};
% data = datacell{2};
% 
% for fff = 1:size(textdata,1);
%     if strcmp(textdata(fff),'$Bodymass')
%         bodyMass = data(fff)
%     elseif strcmp(textdata(fff),'$LLegLength')
%         LLegLength = data(fff)
%     elseif strcmp(textdata(fff),'$RLegLength')
%         RLegLength = data(fff)
%     end
% end
 
%%

fc = 6;%fréquence de coupure
frq = btkGetPointFrequency(acq);
[b,a] = butter(4,fc/(frq/2));

% On va normalisé ca sur un cycle de Right Heel Strike a Right Heel Strike

firstTime = btkGetFirstFrame(acq)/btkGetPointFrequency(acq);
lastTime = btkGetLastFrame(acq)/btkGetPointFrequency(acq);

if isfield(points,'RPelvisAngles') && isfield(points,'LPelvisAngles') ...
        && isfield(points,'RHipAngles') && isfield(points,'LHipAngles')...
        && isfield(points,'RKneeAngles') && isfield(points,'LKneeAngles')...
        && isfield(points,'RAnkleAngles') && isfield(points,'LAnkleAngles')...
        && isfield(points,'RFootProgressAngles') && isfield(points,'LFootProgressAngles')...
        && isfield(points,'RHipMoment') && isfield(points,'LHipMoment')...
        && isfield(points,'RKneeMoment') && isfield(points,'LKneeMoment')...
        && isfield(points,'RAnkleMoment') && isfield(points,'LAnkleMoment')...
        && isfield(points,'RHipPower') && isfield(points,'LHipPower')...
        && isfield(points,'RKneePower') && isfield(points,'LKneePower')...
        && isfield(points,'RAnklePower') && isfield(points,'LAnklePower')...
        && (isfield(points,'RNormalisedGRF') || isfield(points,'LNormalisedGRF'))
    
    %     if all(any(points.LSHO,2))&& all(any(points.RSHO,2))...
    %             && all(any(points.LPSI,2))&& all(any(points.RPSI,2))...
    %             && all(any(points.LFHD,2))&& all(any(points.RFHD,2))...
    %             && all(any(points.LASI,2))&& all(any(points.RASI,2))
    
    if isfield(Events,'Right_Foot_Strike') && isfield(Events,'Left_Foot_Strike') ...
            && isfield(Events,'Right_Foot_Off') && isfield(Events,'Left_Foot_Off')
        
        
        % On enleve les evenement avant la premiere frame
        Events.Right_Foot_Off = Events.Right_Foot_Off (round(Events.Right_Foot_Off * frq) >= round(firstTime * frq) & round(Events.Right_Foot_Off * frq) <= round(lastTime * frq));
        Events.Right_Foot_Strike = Events.Right_Foot_Strike (round(Events.Right_Foot_Strike * frq) >= round(firstTime * frq) & round(Events.Right_Foot_Strike * frq) <= round(lastTime * frq));
        Events.Left_Foot_Strike = Events.Left_Foot_Strike (round(Events.Left_Foot_Strike * frq) >= round(firstTime * frq) &  round(Events.Left_Foot_Strike * frq) <= round(lastTime * frq));
        Events.Left_Foot_Off = Events.Left_Foot_Off (round(Events.Left_Foot_Off * frq) >= round(firstTime * frq) & round(Events.Left_Foot_Off * frq) <= round(lastTime * frq));
        
        % On reteste si c'est pas vide
        if isfield(Events,'Right_Foot_Strike') && isfield(Events,'Left_Foot_Strike') ...
                && isfield(Events,'Right_Foot_Off') && isfield(Events,'Left_Foot_Off')
            
            
            FSstr = [side,'_Foot_Strike'];
            % Cell2Str before we have a cell with strcat using str{} allow to get a
            % string that can be used in dynamic structure calling
            FS = Events.(FSstr);
            
            % Controlateral FS
            %C_FSstr = 'Left_Foot_Strike';
            C_FSstr = [C_side,'_Foot_Strike'];
            C_FS = Events.(C_FSstr);
            
            % Foot off
            %FOstr = 'Right_Foot_Off';
            FOstr = [side,'_Foot_Off'];
            FO = Events.(FOstr);
            
            % Controlateral Foot strike
            %C_FOstr =  'Left_Foot_Off';
            C_FOstr = [C_side,'_Foot_Off'];
            C_FO = Events.(C_FOstr);
            
            FO = FO(FO>FS(1));
            C_FS = C_FS(C_FS>FS(1));
            C_FO = C_FO(C_FO>FS(1));
            
            for ind4 = 1:(size(FS,2)-1)
                % We need to test that between 2 foot strike there is a foot of and
                % a controlateral foot strike and off
                testFO = isempty(FO(FO>FS(ind4) & FO<FS(ind4+1)));
                % If one of this condition is not respected we add at the beginning
                % of the Event vector a 0 value in order to keep the correspondence
                % between the foot strike event and teh other
                if testFO
                    FO = [0,FO];
                end
                
                % If no test is empty (which mean that there is the event between
                % the 2 foot strike) it is possible to compute all the different
                % spation temporal parameter.
                
                if not(testFO)
                    % on transforme les evenemetns qui sont en temps en
                    % frame pour pouvoir les utiliser comme réference dans
                    % les points extrait de c3d
                    RHS1 = round(FS(ind4)*frq);
                    RHS2 = round(FS(ind4+1)*frq);
                    
                    RHS1 = RHS1-btkGetFirstFrame(acq)+1;
                    RHS2 = RHS2-btkGetFirstFrame(acq)+1;
                    
                    for k = 1:num
                        if isfield(points,pointsLabel{k})
                            %n = size(points.(pointsLabel{k})(RHS1:RHS2,:),1);
                            %dataExtracted.(pointsLabel{k}) = points.(pointsLabel{k})(RHS1:RHS2,:);
                            tempdata = points.(pointsLabel{k})(RHS1:RHS2,:);
                            tempdataButterFilt = filtfilt(b,a,points.(pointsLabel{k})(RHS1:RHS2,:));
                            %dataExtracted.(pointsLabel{k}) = interp1((1:n)',points.(pointsLabel{k})(RHS1:RHS2,:),...
                            %    linspace(1,n,101)','spline');
                            
                            n = size(tempdataButterFilt,1);
                            
                            dataExtracted.(pointsLabel{k}) = interp1((1:n)',tempdataButterFilt,...
                                linspace(1,n,101)','spline');
                        end
                        
                        
                    end
                    
                    forcePlateforme = btkGetGroundReactionWrenches(acq);
                    factor = btkGetAnalogFrequency(acq)/btkGetPointFrequency(acq);
                    % fichier = dir([fileName(1:end-4),'*Trial*.enf']);==>
                    % changer pour la version du dessous pour éviter les
                    % confusion entre 1 et 10
                    fichier = dir([fileName(1:end-4),'.Trial*.enf']);
                    fileNameEnf = fullfile(pathname,fichier.name);
                    %fileNameEnf = strcat(fileName(1:end-4),'.','Trial',fileName(end-5:end-4),'.enf');
                    
                    numeroPlateform = plateformDetection_MG(fileNameEnf,side);
                    
                    for ind5 = numeroPlateform
                        n0plateforme = ind5;
                        
                        nbrFrameWthPlat = sum(any(abs(forcePlateforme(n0plateforme).F(RHS1*factor:RHS2*factor,:))>5,2));
                        perctFrameWthPlat = nbrFrameWthPlat/(RHS2*factor-RHS1*factor);
                        
                        if  perctFrameWthPlat>0.05

                            % Ancien bout de code avant que l'extraction
                            % des plateforme ne viennent pas du c3d
%                             for ind6=1:3
%                                 tempForcePlateform = forcePlateforme(n0plateforme).F(RHS1*factor:RHS2*factor,:);
%                                 n = size(forcePlateforme(n0plateforme).F(RHS1*factor:RHS2*factor,:),1);
%                                 forcePlateformExtracted = interp1((1:n)',tempForcePlateform,...
%                                     linspace(1,n,101)','spline');
%                                 forcePlateformNormalize = forcePlateformExtracted*(1/(bodyMass*9.81));
%                             end
                            
%                             if strcmpi(side,'right')
%                                 positionFileName = strfind(fileName,'\')
%                                 positionFileName = positionFileName(end)+1
%                                 titlefigure = fileName(positionFileName:end-4)
%                                 titlefigure = strcat(titlefigure,'_',side)
%                                 figure('Name',titlefigure)
%                                 subplot(1,2,1)
%                                 hold on
%                                 plot(points.RAnkleMoment(:,1))
%                                 plot([RHS1,RHS1],[-100,100])
%                                 plot([RHS2,RHS2],[-100,100])
%                                 hold off
%                                 subplot(1,2,2)
%                                 hold on
%                                 plot(points.RKneeAngles(:,1))
%                                 plot([RHS1,RHS1],[-10,70])
%                                 plot([RHS2,RHS2],[-10,70])
%                                 hold off
%                                 
%                                 GarderAdroite = questdlg('Garder cette courbe?');
                                
%                                 if strcmp(GarderAdroite,'Yes')
                                    % avec PROCESSING
                                    
                                    if strcmpi(side,'right')
                                        LegLenght = metaData.children.PROCESSING.children.RLegLength.info.values/1000;
                                        % sans PROCESSING
                                        %                                     LegLenght = RLegLength/1000;
                                        
                                        % Cinématique
                                        meanAngle.Ptilt(:,indglob) = dataExtracted.RPelvisAngles(:,1);
                                        meanAngle.Pobli(:,indglob) = dataExtracted.RPelvisAngles(:,2);
                                        meanAngle.Prota(:,indglob) = dataExtracted.RPelvisAngles(:,3);
                                        
                                        meanAngle.Hflex(:,indglob) = dataExtracted.RHipAngles(:,1);
                                        meanAngle.Habd(:,indglob) = dataExtracted.RHipAngles(:,2);
                                        meanAngle.Hrota(:,indglob) = dataExtracted.RHipAngles(:,3);
                                        
                                        meanAngle.Kflex(:,indglob) = dataExtracted.RKneeAngles(:,1);
                                        meanAngle.Kabd(:,indglob) = dataExtracted.RKneeAngles(:,2);
                                        meanAngle.Krota(:,indglob) = dataExtracted.RKneeAngles(:,3);
                                        
                                        meanAngle.Aflex(:,indglob) = dataExtracted.RAnkleAngles(:,1);
                                        
                                        meanAngle.Frota(:,indglob) = dataExtracted.RFootProgressAngles(:,3);
                                        meanAngle.Fvert(:,indglob) = -dataExtracted.RFootProgressAngles(:,1)-90;
                                        
                                        NormMoment = bodyMass*9.81*LegLenght;
                                        
                                        % Moment
                                        meanAngle.HipExtMoment = dataExtracted.RHipMoment(:,1)/NormMoment;
                                        meanAngle.HipAbdMoment = dataExtracted.RHipMoment(:,2)/NormMoment;
                                        meanAngle.HipRotMoment = dataExtracted.RHipMoment(:,3)/NormMoment;
                                        
                                        meanAngle.KneeExtMoment = dataExtracted.RKneeMoment(:,1)/NormMoment;
                                        meanAngle.KneeAbdMoment = dataExtracted.RKneeMoment(:,2)/NormMoment;
                                        meanAngle.KneeRotMoment = dataExtracted.RKneeMoment(:,3)/NormMoment;
                                        
                                        meanAngle.AnklePlantMoment = dataExtracted.RAnkleMoment(:,1)/NormMoment;
                                        %                                 meanAngle.KneeRotMoment = dataExtracted.RAnkleMoment(:,2);
                                        %                                 meanAngle.KneeEvertMoment = dataExtracted.RAnkleMoment(:,3);
                                        
                                        % Power
                                        meanAngle.HipPower(:,indglob) = dataExtracted.RHipPower(:,3);
                                        meanAngle.KneePower(:,indglob) = dataExtracted.RKneePower(:,3);
                                        meanAngle.AnklePower(:,indglob) = dataExtracted.RAnklePower(:,3);
                                        
                                        % Plateforme
                                        %                                     meanAngle.LatMedGRF(:,indglob) = forcePlateformNormalize(:,1);
                                        %                                     meanAngle.AntPostGRF(:,indglob) = forcePlateformNormalize(:,2);
                                        %                                     meanAngle.VerticalGRF(:,indglob) = forcePlateformNormalize(:,3);
                                        
                                        meanAngle.LatMedGRF(:,indglob) = dataExtracted.RNormalisedGRF(:,1)/100;
                                        meanAngle.AntPostGRF(:,indglob) = dataExtracted.RNormalisedGRF(:,2)/100;
                                        meanAngle.VerticalGRF(:,indglob) = dataExtracted.RNormalisedGRF(:,3)/100;
                                        
                                        %                                 end
                                        
                                        %                             elseif strcmpi(side,'left')
                                        %                                 positionFileName = strfind(fileName,'\')
                                        %                                 positionFileName = positionFileName(end)+1
                                        %                                 titlefigure = fileName(positionFileName:end-4)
                                        %                                 titlefigure = strcat(titlefigure,'_',side)
                                        %                                 figure('Name',titlefigure)
                                        %                                 subplot(1,2,1)
                                        %                                 hold on
                                        %                                 plot(points.LAnkleMoment(:,1))
                                        %                                 plot([RHS1,RHS1],[-10,70])
                                        %                                 plot([RHS2,RHS2],[-10,70])
                                        %                                 hold off
                                        %                                 subplot(1,2,2)
                                        %                                 hold on
                                        %                                 plot(points.LKneeAngles(:,1))
                                        %                                 plot([RHS1,RHS1],[-10,70])
                                        %                                 plot([RHS2,RHS2],[-10,70])
                                        %                                 hold off
                                        %                                 GarderAgauche = questdlg('Garder cette courbe?');
                                        %
                                        %                                 if strcmp(GarderAgauche,'Yes')
                                        % avec PROCESSING
                                    elseif strcmpi(side,'left')
                                        LegLenght = metaData.children.PROCESSING.children.LLegLength.info.values/1000;
                                        % sans PROCESSING
                                        %                                     LegLenght = LLegLength/1000;
                                        
                                        meanAngle.Ptilt(:,indglob) = dataExtracted.LPelvisAngles(:,1);
                                        meanAngle.Pobli(:,indglob) = dataExtracted.LPelvisAngles(:,2);
                                        meanAngle.Prota(:,indglob) = dataExtracted.LPelvisAngles(:,3);
                                        
                                        meanAngle.Hflex(:,indglob) = dataExtracted.LHipAngles(:,1);
                                        meanAngle.Habd(:,indglob) = dataExtracted.LHipAngles(:,2);
                                        meanAngle.Hrota(:,indglob) = dataExtracted.LHipAngles(:,3);
                                        
                                        meanAngle.Kflex(:,indglob) = dataExtracted.LKneeAngles(:,1);
                                        meanAngle.Kabd(:,indglob) = dataExtracted.LKneeAngles(:,2);
                                        meanAngle.Krota(:,indglob) = dataExtracted.LKneeAngles(:,3);
                                        
                                        meanAngle.Aflex(:,indglob) = dataExtracted.LAnkleAngles(:,1);
                                        
                                        meanAngle.Frota(:,indglob) = dataExtracted.LFootProgressAngles(:,3);
                                        meanAngle.Fvert(:,indglob) = -dataExtracted.LFootProgressAngles(:,1)-90;
                                        
                                        NormMoment = bodyMass*9.81*LegLenght;
                                        % Moment
                                        meanAngle.HipExtMoment = dataExtracted.LHipMoment(:,1)/NormMoment;
                                        meanAngle.HipAbdMoment = dataExtracted.LHipMoment(:,2)/NormMoment;
                                        meanAngle.HipRotMoment = dataExtracted.LHipMoment(:,3)/NormMoment;
                                        
                                        meanAngle.KneeExtMoment = dataExtracted.LKneeMoment(:,1)/NormMoment;
                                        meanAngle.KneeAbdMoment = dataExtracted.LKneeMoment(:,2)/NormMoment;
                                        meanAngle.KneeRotMoment = dataExtracted.LKneeMoment(:,3)/NormMoment;
                                        
                                        meanAngle.AnklePlantMoment = dataExtracted.LAnkleMoment(:,1)/NormMoment;
                                        %                                 meanAngle.KneeRotMoment = dataExtracted.RAnkleMoment(:,2);
                                        %                                 meanAngle.KneeEvertMoment = dataExtracted.RAnkleMoment(:,3);
                                        
                                        % Power
                                        meanAngle.HipPower(:,indglob) = dataExtracted.LHipPower(:,3);
                                        meanAngle.KneePower(:,indglob) = dataExtracted.LKneePower(:,3);
                                        meanAngle.AnklePower(:,indglob) = dataExtracted.LAnklePower(:,3);
                                        
                                        % Plateforme
                                        %                                     meanAngle.LatMedGRF(:,indglob) = forcePlateformNormalize(:,1);
                                        %                                     meanAngle.AntPostGRF(:,indglob) = forcePlateformNormalize(:,2);
                                        %                                     meanAngle.VerticalGRF(:,indglob) = forcePlateformNormalize(:,3);
                                        
                                        meanAngle.LatMedGRF(:,indglob) = dataExtracted.LNormalisedGRF(:,1)/100;
                                        meanAngle.AntPostGRF(:,indglob) = dataExtracted.LNormalisedGRF(:,2)/100;
                                        meanAngle.VerticalGRF(:,indglob) = dataExtracted.LNormalisedGRF(:,3)/100;
                                    end
%                                 end
                                
%                             end
                            
                            indglob = indglob+1;
                        end
                    end
                end
            end
        end
    end
    %end
end
