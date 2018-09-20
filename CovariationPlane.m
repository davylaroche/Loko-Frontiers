function [COEFF, EXPLAINED]=CovariationPlane(kinematic,indglob)

%%% Construction des matrices 
Othigh = [kinematic.PELA(:,indglob*3-2:indglob*3)-kinematic.PELO(:,indglob*3-2:indglob*3) kinematic.PELL(:,indglob*3-2:indglob*3)-kinematic.PELO(:,indglob*3-2:indglob*3) kinematic.PELP(:,indglob*3-2:indglob*3)-kinematic.PELO(:,indglob*3-2:indglob*3)];
Oshank = [kinematic.FEA(:,indglob*3-2:indglob*3)-kinematic.FEO(:,indglob*3-2:indglob*3) kinematic.FEL(:,indglob*3-2:indglob*3)-kinematic.FEO(:,indglob*3-2:indglob*3) kinematic.FEP(:,indglob*3-2:indglob*3)-kinematic.FEO(:,indglob*3-2:indglob*3)];
Ofoot = [kinematic.TIA(:,indglob*3-2:indglob*3)-kinematic.TIO(:,indglob*3-2:indglob*3) kinematic.TIL(:,indglob*3-2:indglob*3)-kinematic.TIO(:,indglob*3-2:indglob*3) kinematic.TIP(:,indglob*3-2:indglob*3)-kinematic.TIO(:,indglob*3-2:indglob*3)];
distaFoot = [kinematic.TOA(:,indglob*3-2:indglob*3)-kinematic.TOO(:,indglob*3-2:indglob*3) kinematic.TOL(:,indglob*3-2:indglob*3)-kinematic.TOO(:,indglob*3-2:indglob*3) kinematic.TOP(:,indglob*3-2:indglob*3)-kinematic.TOO(:,indglob*3-2:indglob*3)];

%%% define vertical matrix
n=size(Othigh, 1);
mz = [0 0 1];
for c1=1:n
    clear mx my mx2
    mx = kinematic.PELA(c1,indglob*3-2:indglob*3)-kinematic.PELO(c1, indglob*3-2:indglob*3);
    my = cross(mz, mx);
    mx2 = cross(my, mz);
    myp = my./(sum(my.^2)^0.5);
    mxp = mx2./(sum(mx2.^2)^0.5);
    O(c1, :) = [mxp myp mz];
end

% Othigh = [points.PELA(:, [2 1 3])-points.PELO(:, [2 1 3]) points.PELL(:, [2 1 3])-points.PELO(:, [2 1 3]) points.PELP(:, [2 1 3])-points.PELO(:, [2 1 3])];
% Oshank = [points.LFEA(:, [2 1 3])-points.LFEO(:, [2 1 3]) points.LFEL(:, [2 1 3])-points.LFEO(:, [2 1 3]) points.LFEP(:, [2 1 3])-points.LFEO(:, [2 1 3])];
% Ofoot = [points.LTIA(:, [2 1 3])-points.LTIO(:, [2 1 3]) points.LTIL(:, [2 1 3])-points.LTIO(:, [2 1 3]) points.LTIP(:, [2 1 3])-points.LTIO(:, [2 1 3])];
% distaFoot = [points.LTOA(:, [2 1 3])-points.LTOO(:, [2 1 3]) points.LTOL(:, [2 1 3])-points.LTOO(:, [2 1 3]) points.LTOP(:, [2 1 3])-points.LTOO(:, [2 1 3])];

%%% résolution des angles matrix to matrix
% [thigh, T1] = angles_solver(Othigh, Oshank, 'YXZ');
% [shank, T2] = angles_solver(Oshank, Ofoot, 'YXZ');
% [foot, T3] = angles_solver(Ofoot, distaFoot, 'YXZ');

[thigh, T1] = angles_solver(O, Oshank, 'YXZ');
[shank, T2] = angles_solver(O, Ofoot, 'YXZ');
[foot, T3] = angles_solver(O, distaFoot, 'XYZ');

%% Centrer / réduire données avant PCA
M = [thigh(:,1) shank(:,1) foot(:,1)];
Mc = M - (ones(size(M, 1), 1) * nanmean(M));
Mcr = Mc/diag(nanstd(Mc, 1));

%% Effectuer la PCA sur la matrice de covariance
V = cov(Mcr);
[COEFF,~,EXPLAINED] = pcacov(V);
