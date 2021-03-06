function [out] = rxzysolv(T)
% Description:	Solves for alpha,beta,gama,Hx,Hy,Hz of the transformation
%               matrix with the order Rxzy ([ry] [rz] [rx])
%
% Input:	- T      Transformation matrix
% Output:	- out    [alpha,beta,gama,Hx,Hy,Hz]
%                        Note that the angles are given in the range 
%                        from -180 to +180 deg.
%                        gama is not a typo: gamma could not be used because
%                        it is an existing matlab function 
% Author:	Christoph Reinschmidt, HPL, The University of Calgary
% Date:		October, 1994
% Last changes: February 13, 1995
% Version:	1.0

if size(T)~=[4,4]; 
   disp('Error: transformation matrix has to be a 4x4 matrix')
   return;
end;

if sum(isnan(T(:)))~=0, out=[NaN,NaN,NaN,NaN,NaN,NaN]; return; end

gama = asin(T(2,1));   %'assumption' that cos(gama)>0

alphasin = asin(-T(2,3)/cos(gama)); 
alphacos = acos(T(2,2)/cos(gama)); 
if (alphacos>pi/2 & alphasin>0); alpha=pi-alphasin; end;
if (alphacos>pi/2 & alphasin<0); alpha=-pi-alphasin; end; 
if (alphacos<=pi/2); alpha=alphasin; end;

betasin = asin(-T(3,1)/cos(gama));
betacos = acos(T(1,1)/cos(gama));
if (betacos>pi/2 & betasin>0); beta=pi-betasin; end;
if (betacos>pi/2 & betasin<0); beta=-pi-betasin; end; 
if (betacos<=pi/2); beta=betasin; end;



% Calculation of Hx,Hy,Hz

A=[cos(beta)*cos(gama),0, sin(beta);...
   sin(gama),1,0;...
   -sin(beta)*cos(gama),0,cos(beta)];
H=A\T(1:3,4); H=H';


out=[rad2deg(alpha),rad2deg(beta),rad2deg(gama),H];

