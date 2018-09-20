function [meanAngle] = moyenneAllCycleKinematic_MG(fileName,side)

pointsLabel = {'RPelvisAngles', 'LPelvisAngles',...
    'RHipAngles', 'LHipAngles',...
    'RKneeAngles', 'LKneeAngles',...
    'RAnkleAngles', 'LAnkleAngles',...
    'LFootProgressAngles','RFootProgressAngles',...
    'PELA','PELO','PELL','PELP',...
    'LFEA','LFEO','LFEL','LFEP',...
    'LTIA','LTIO','LTIL','LTIP',...
    'LTOA','LTOO','LTOL','LTOP',...
    'RFEA','RFEO','RFEL','RFEP',...
    'RTIA','RTIO','RTIL','RTIP',...
    'RTOA','RTOO','RTOL','RTOP'};
num = length(pointsLabel);
indglob = 1;


pointsLabelStruct = {'Ptilt', 'Pobli','Prota',...
    'Hflex','Habd','Hrota',...
    'Kflex','Kabd','Krota',...
    'Aflex','Frota','Fvert',...
    'PELA','PELO','PELL','PELP',...
    'FEA','FEO','FEL','FEP',...
    'TIA','TIO','TIL','TIP',...
    'TOA','TOO','TOL','TOP',...
    'PlanCovIndex','PlanCovAngle'};

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
frq = btkGetPointFrequency(acq);
% On va normalisé ca sur un cycle de Right Heel Strike a Right Heel Strike

fc = 6;%fréquence de coupure
[b,a] = butter(4,fc/(frq/2));

firstTime = btkGetFirstFrame(acq)/btkGetPointFrequency(acq);
lastTime = btkGetLastFrame(acq)/btkGetPointFrequency(acq);

if isfield(points,'RPelvisAngles') && isfield(points,'LPelvisAngles') ...
        && isfield(points,'RHipAngles') && isfield(points,'LHipAngles')...
        && isfield(points,'RKneeAngles') && isfield(points,'LKneeAngles')...
        && isfield(points,'RAnkleAngles') && isfield(points,'LAnkleAngles')...
        && isfield(points,'RFootProgressAngles') && isfield(points,'LFootProgressAngles')
    
    
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
                % between the foot strike event and
                if testFO
                    FO = [0,FO];
                end
                
                % If no test is empty (which mean that there is the event between
                % the 2 foot strike) it is possible to compute all the different
                % spation temporal parameter.
                
                if not(testFO)
                    
                    RHS1 = round(FS(ind4)*frq);
                    RHS2 = round(FS(ind4+1)*frq);
                    
                    RHS1 = RHS1-btkGetFirstFrame(acq)+1;
                    RHS2 = RHS2-btkGetFirstFrame(acq)+1;
                    
                    for k = 1:num
                        
                        tempdata = points.(pointsLabel{k})(RHS1:RHS2,:);
                        tempdataButterFilt = filtfilt(b,a,tempdata);
                        n = size(tempdataButterFilt,1);
                        %dataExtracted.(pointsLabel{k}) = points.(pointsLabel{k})(RHS1:RHS2,:);
                        
                        
                        dataExtracted.(pointsLabel{k}) = interp1((1:n)',points.(pointsLabel{k})(RHS1:RHS2,:),...
                            linspace(1,n,101)','spline');
                    end
                    
                    
                    if strcmp(side,'Right')
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
                        
                                                
                    else
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
                    end
                    
                    if strcmp(side,'Right')
                        meanAngle.PELA(:,indglob*3-2:indglob*3) = dataExtracted.PELA;
                        meanAngle.PELO(:,indglob*3-2:indglob*3) = dataExtracted.PELO;
                        meanAngle.PELL(:,indglob*3-2:indglob*3) = dataExtracted.PELL;
                        meanAngle.PELP(:,indglob*3-2:indglob*3) = dataExtracted.PELP;
                                                
                        meanAngle.FEA(:,indglob*3-2:indglob*3) = dataExtracted.RFEA;
                        meanAngle.FEO(:,indglob*3-2:indglob*3) = dataExtracted.RFEO;
                        meanAngle.FEL(:,indglob*3-2:indglob*3) = dataExtracted.RFEL;
                        meanAngle.FEP(:,indglob*3-2:indglob*3) = dataExtracted.RFEP;
                        
                        meanAngle.TIA(:,indglob*3-2:indglob*3) = dataExtracted.RTIA;
                        meanAngle.TIO(:,indglob*3-2:indglob*3) = dataExtracted.RTIO;
                        meanAngle.TIL(:,indglob*3-2:indglob*3) = dataExtracted.RTIL;
                        meanAngle.TIP(:,indglob*3-2:indglob*3) = dataExtracted.RTIP;
                        
                        meanAngle.TOA(:,indglob*3-2:indglob*3) = dataExtracted.RTOA;
                        meanAngle.TOO(:,indglob*3-2:indglob*3) = dataExtracted.RTOO;
                        meanAngle.TOL(:,indglob*3-2:indglob*3) = dataExtracted.RTOL;
                        meanAngle.TOP(:,indglob*3-2:indglob*3) = dataExtracted.RTOP;
                        
                    else 
                        meanAngle.PELA(:,indglob*3-2:indglob*3) = dataExtracted.PELA;
                        meanAngle.PELO(:,indglob*3-2:indglob*3) = dataExtracted.PELO;
                        meanAngle.PELL(:,indglob*3-2:indglob*3) = dataExtracted.PELL;
                        meanAngle.PELP(:,indglob*3-2:indglob*3) = dataExtracted.PELP;
                                                
                        meanAngle.FEA(:,indglob*3-2:indglob*3) = dataExtracted.LFEA;
                        meanAngle.FEO(:,indglob*3-2:indglob*3) = dataExtracted.LFEO;
                        meanAngle.FEL(:,indglob*3-2:indglob*3) = dataExtracted.LFEL;
                        meanAngle.FEP(:,indglob*3-2:indglob*3) = dataExtracted.LFEP;
                        
                        meanAngle.TIA(:,indglob*3-2:indglob*3) = dataExtracted.LTIA;
                        meanAngle.TIO(:,indglob*3-2:indglob*3) = dataExtracted.LTIO;
                        meanAngle.TIL(:,indglob*3-2:indglob*3) = dataExtracted.LTIL;
                        meanAngle.TIP(:,indglob*3-2:indglob*3) = dataExtracted.LTIP;
                        
                        meanAngle.TOA(:,indglob*3-2:indglob*3) = dataExtracted.LTOA;
                        meanAngle.TOO(:,indglob*3-2:indglob*3) = dataExtracted.LTOO;
                        meanAngle.TOL(:,indglob*3-2:indglob*3) = dataExtracted.LTOL;
                        meanAngle.TOP(:,indglob*3-2:indglob*3) = dataExtracted.LTOP;
                        
                    end
                    
                    %% traitement Plan de Covariation
                    [COEFF, EXPLAINED]=CovariationPlane(meanAngle,indglob); 
                    meanAngle.PlanCovIndex(:,indglob)=EXPLAINED;
                    meanAngle.PlanCovAngle(:,indglob*3-2:indglob*3)=COEFF;     
                                       
                    %%
                    indglob = indglob+1;
                    
                end
            end
        end
    end
    %end
end
