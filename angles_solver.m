function [angles, T] = angles_solver(O, M, rotation)
% This function compute 3D angles by succcessive rotation of rotation input
% rotation input could be XYZ XZY YXZ YZX ZXY ZYX
% O input matrix is the reference matrix, Angles will be calculated fonction of it must be 1*9 or m*9 matrix ordered as X Y Z
% M input matrix must be m*9 matrix ordered as X Y Z
% T output is the 4*4 rotation matrix [R, d;0 0 0 1] with R rotation matrix [Rx Ry Rz]
% and d translation matrix [dx dy dz]

if nargin ~= 3
    error 'wrong number of input arguments'
end

if size(O, 1) ~= 1 && size(O, 1) ~= size(M, 1)
    error 'reference matrix is bad calibrated'
end

if size(O, 2) ~= 9 || size(M, 2) ~= 9 
    error 'reference matrix or destination maxtrix are bad calibrated'
end   
    
    
for yy = 1:size(M, 1)
    clear T datain
    if size(O, 1) > 1, datain = [O(yy, :); M(yy, :)];else datain = [O; M(yy, :)];end
    [T] = SODER(datain);
    
    if ~isnan(det(T))
        switch rotation
            case 'XYZ'
                angles(yy, :) = RXYZSOLV(T);
            case 'XZY'
                angles(yy, :) = RXZYSOLV(T);
            case 'YXZ'
                angles(yy, :) = RYXZSOLV(T);
            case 'YZX'
                angles(yy, :) = RYZXSOLV(T);
            case 'ZXY'
                angles(yy, :) = RZXYSOLV(T);
            case 'ZYX'
                angles(yy, :) = RZYXSOLV(T);   
            otherwise
                error 'such type of rotation does not exist'
        end
    else
        angles(yy, 1:3) = NaN;   
    end
end

angles = angles(:, 1:3);

end

