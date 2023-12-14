% ------------------------------------------------------------------------------
%
%                            function angl
%
%  this function calculates the angle between two vectors.  the output is
%    set to 999999.1 to indicate an undefined value.  be sure to check for
%    this at the output phase.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - fix tolerances                                 5 sep 2002
%
%  inputs          description                    range / units
%    vec1        - vector number 1
%    vec2        - vector number 2
%
%  outputs       :
%    theta       - angle between the two vectors  -pi to pi
%
%  locals        :
%    temp        - temporary real variable
%
%  coupling      :
%
% [theta] = angl ( vec1,vec2 );
% ----------------------------------------------------------------------------- }

function [theta] = angl_vec ( vec1,vec2 )

%     small     = 0.00000001;
    undefined = 999999.1;

    magv1 = sqrt(sum(vec1.^2,2));
    magv2 = sqrt(sum(vec2.^2,2));
    
    theta = zeros(size(magv1));
    mag_prod = magv1.*magv2;
    check_mag = mag_prod > 1e-16;
    
    temp = dot(vec1(check_mag,:),vec2(check_mag,:),2) ./ mag_prod(check_mag);
    check_temp = abs( temp ) > 1.0;
    temp(check_temp) = sign(temp(check_temp)) * 1.0;
    theta(check_mag) = acos( temp );
    theta(~check_mag) = undefined;