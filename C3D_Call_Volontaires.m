% Permet de lister tous les .c3d de tous les patients
% Patients concernés : inscrits dans le fichier excel
% Protocole : Locox 2
clc; clear all; 

%% Chemin ordi VICON
% FILENAME_DBPATIENTS = 'D:\projets terminés\LOCOX_2\liste_volontaires.xlsx';
% folder = 'D:\projets terminés\LOCOX_2\Volontaires\';

%% Chemin ordi MG 
FILENAME_DBPATIENTS = 'Y:\Backup PIT\Backup Vicon\projets terminés\LOCOX_2\liste_volontaires.xlsx';
folder = 'Y:\Backup PIT\Backup Vicon\projets terminés\LOCOX_2\Volontaires\';
repertorySave = 'B:\Documents\4. Valorisation\2. Articles\4. Puissance & vitesse marche\Results_Locox_Frontiers\';

% %% Chemin ordi DL 
% FILENAME_DBPATIENTS = 'Z:\Backup PIT\Backup Vicon\projets terminés\LOCOX_2\liste_volontaires.xlsx';
% folder = 'Z:\Backup PIT\Backup Vicon\projets terminés\LOCOX_2\Volontaires\';
% repertorySave = 'C:\Users\davy.laroche\Documents\';

%% listing subject
[~,~,DB] = xlsread(FILENAME_DBPATIENTS);

id = char(DB(1:end,1));

%% Listing des variables de sorties sur le fichier EXCEL
tableauExcelTotal = {'Sujet','vitesseMarche',...
    'AmplitudeCoMVertical','RoMSagittalHip','RoMSagittalKnee','RoMSagittalAnkle',...
    'PlanCovIndex','PlanCovIndexu1','PlanCovIndexu3','StiffnessKnee','StiffnessAnkle','PmaxHip','PmaxKnee','PmaxAnkle',...
    'IntegralePHipPositive','IntegralePHipNegative','IntegralePKneePositive','IntegralePKneeNegative','IntegralePAnklePositive','IntegralePAnkleNegative'};

%% Listing C3D par Volontaires
for c1 = 1 : length(id)
    c1
    id(c1, :)
    
    clear kinetic C3Dfolder listec3d file_list
    
    C3Dfolder = [folder,id(c1, :),'\Marche\'];
    listec3d = dir(fullfile([C3Dfolder,'*.c3d']));
    
    file_list={listec3d.name};
       
    %% CoM / Cinématique & Cinétique
    [kinetic] = calculAllKinetic_MG(C3Dfolder,file_list);
    
    %% Traitement et extraction variables cibles
    % Sélection gauche/droite
    clear KINsel
    
    if mod(c1,2) == 0; % nombre pair
        KINsel = kinetic.Left;
    elseif mod(c1,2) == 1; % nombre impair
        KINsel = kinetic.Right;
    end
    
    clear var2extract
    
    % Vitesse de marche
    var2extract.GaitSpeedCoM =  KINsel.SpeedWalkingCoM';
    
    % Amplitude max des mouvements verticaux du CoM pdt cycle de marche
    var2extract.CoM_vert =  KINsel.CoM_vert';
    
    % RoM Hanche / Genou et Cheville 
    clear i
    for i = 1:size(KINsel.Hflex,2)
        if max(KINsel.Hflex(:,i))>0 && min(KINsel.Hflex(:,i))<0
           var2extract.RoMSagittalHip(i,1) =  max(KINsel.Hflex(:,i))+abs(min(KINsel.Hflex(:,i)));
        else
           var2extract.RoMSagittalHip(i,1) =  max(KINsel.Hflex(:,i))-min(KINsel.Hflex(:,i));
        end
        
        if max(KINsel.Kflex(:,i))>0 && min(KINsel.Kflex(:,i))<0
           var2extract.RoMSagittalKnee(i,1) =  max(KINsel.Kflex(:,i))+abs(min(KINsel.Kflex(:,i)));
        else
           var2extract.RoMSagittalKnee(i,1) =  max(KINsel.Kflex(:,i))-min(KINsel.Kflex(:,i));
        end
        
        if max(KINsel.Aflex(:,i))>0 && min(KINsel.Aflex(:,i))<0
           var2extract.RoMSagittalAnkle(i,1) =  max(KINsel.Aflex(:,i))+abs(min(KINsel.Aflex(:,i)));
        else
           var2extract.RoMSagittalAnkle(i,1) =  max(KINsel.Aflex(:,i))-min(KINsel.Aflex(:,i));
        end
    end
    
    % Coordination : PI, u1 et u3
    clear i
    for i = 1:size(KINsel.PlanCovIndex,2)
        var2extract.PlanCovIndex(i,1) =  sum(KINsel.PlanCovIndex(1:2,i));
    end
    
    var2extract.PlanCovIndex_u1 =  KINsel.PlanCov_u1';
    var2extract.PlanCovIndex_u3 =  KINsel.PlanCov_u3';
    
    % Stiffness Genou et Cheville entre HeelStrike et MidStance (moment du
    % pic) - ne prendre que KINsel même pour la cinématique afin de ne
    % calculer la stiffness uniquement sur les pas plateformes (lorsque les
    % moments sont calculables)
    clear i
%     [maxKneeMoment,index_maxKneeMoment] = max(KINsel.KneeExtMoment,[],1);
%     [maxAnkleMoment,index_maxAnkleMoment] = max(KINsel.AnklePlantMoment,[],1);

    [maxCoMvert,index_maxCoMvert] = max(KINsel.CoM,[],1);
    
    for i = 1:size(KINsel.KneeExtMoment,2)
        clear deltaKneeMoment deltaAnkleMoment deltaKneeAngle deltaAnkleAngle K_Knee K_Ankle
        deltaKneeMoment = KINsel.KneeExtMoment(index_maxCoMvert(i),i)-KINsel.KneeExtMoment(1,i);
        deltaAnkleMoment = KINsel.AnklePlantMoment(index_maxCoMvert(i),i)-KINsel.AnklePlantMoment(1,i);
        
        deltaKneeMoment = deltaKneeMoment/1000;
        deltaAnkleMoment = deltaAnkleMoment/1000;
        
        deltaKneeAngle = KINsel.Kflex(index_maxCoMvert(i),i)-KINsel.Kflex(1,i);
        deltaAnkleAngle = KINsel.Aflex(index_maxCoMvert(i),i)-KINsel.Aflex(1,i);
        
        K_Knee = deltaKneeMoment/deltaKneeAngle;
        K_Ankle = deltaAnkleMoment/deltaAnkleAngle;
        
        var2extract.StiffnessKnee(i,1) =  K_Knee;
        var2extract.StiffnessAnkle(i,1) =  K_Ankle;
    end  
    
    % Puissance max Hanche / Genou et Cheville
    clear i
    for i = 1:size(KINsel.HipPower,2)
        var2extract.PmaxHip(i,1) =  max(KINsel.HipPower(:,i));
        var2extract.PmaxKnee(i,1) =  max(KINsel.KneePower(:,i));
        var2extract.PmaxAnkle(i,1) =  max(KINsel.AnklePower(:,i));
    end
    
    % Intégrale de la courbe Puissance Absorbée et Puissance Générée Hanche / Genou et Cheville 
    clear i 
    for i = 1:size(KINsel.HipPower,2)
        % Hips
        clear ii HipP Neg_HipP Pos_HipP        
        HipP = KINsel.HipPower(:,i);
        Neg_HipP = HipP;
        Pos_HipP = HipP;
        
        for ii = 1: size(HipP,1)
            if Neg_HipP(ii)>0;
               Neg_HipP(ii)=0;
            end
            if Pos_HipP(ii)<0;
               Pos_HipP(ii)=0;
            end
        end
       
        var2extract.Integrale_PHip_Positive(i,1) = sum(Pos_HipP);
        var2extract.Integrale_PHip_Negative(i,1) = sum(abs(Neg_HipP));
        
        % Knee
        clear ii KneeP Neg_KneeP Pos_KneeP
        KneeP = KINsel.KneePower(:,i);
        Neg_KneeP = KneeP;
        Pos_KneeP = KneeP;
        
        for ii = 1: size(KneeP,1)
            if Neg_KneeP(ii)>0;
               Neg_KneeP(ii)=0;
            end
            if Pos_KneeP(ii)<0;
               Pos_KneeP(ii)=0;
            end
        end
       
        var2extract.Integrale_PKnee_Positive(i,1) = sum(Pos_KneeP);
        var2extract.Integrale_PKnee_Negative(i,1) = sum(abs(Neg_KneeP));
       
        % Ankle
        clear ii AnkleP Neg_AnkleP Pos_AnkleP 
        AnkleP = KINsel.AnklePower(:,i);
        Neg_AnkleP = AnkleP;
        Pos_AnkleP = AnkleP;
        
        for ii = 1: size(AnkleP,1)
            if Neg_AnkleP(ii)>0;
               Neg_AnkleP(ii)=0;
            end
            if Pos_AnkleP(ii)<0;
               Pos_AnkleP(ii)=0;
            end
        end
       
        var2extract.Integrale_PAnkle_Positive(i,1) = sum(Pos_AnkleP);
        var2extract.Integrale_PAnkle_Negative(i,1) = sum(abs(Neg_AnkleP)); 
        
    end
    
    %% Copier données dans fichier Excel
    for iii = 1:size(var2extract.RoMSagittalHip,1);
        Recapitulatifbis = {char(id(c1,:)),var2extract.GaitSpeedCoM(iii)...
            var2extract.CoM_vert(iii),var2extract.RoMSagittalHip(iii), var2extract.RoMSagittalKnee(iii), var2extract.RoMSagittalAnkle(iii), ...
            var2extract.PlanCovIndex(iii), var2extract.PlanCovIndex_u1(iii), var2extract.PlanCovIndex_u3(iii), var2extract.StiffnessKnee(iii), ...
            var2extract.StiffnessAnkle(iii), var2extract.PmaxHip(iii), var2extract.PmaxKnee(iii), var2extract.PmaxAnkle(iii),...
            var2extract.Integrale_PHip_Positive(iii), var2extract.Integrale_PHip_Negative(iii),var2extract.Integrale_PKnee_Positive(iii), var2extract.Integrale_PKnee_Negative(iii),...
            var2extract.Integrale_PAnkle_Positive(iii), var2extract.Integrale_PAnkle_Negative(iii)};
        
        tableauExcelTotal  = [tableauExcelTotal;Recapitulatifbis];
        
        clear Recapitulatifbis 
    end
end    

% Enregistrement des fichiers excel de données - TOU SUJETS
cd(repertorySave)
xlswrite('LOCOX_Frontiers.xls',tableauExcelTotal)



