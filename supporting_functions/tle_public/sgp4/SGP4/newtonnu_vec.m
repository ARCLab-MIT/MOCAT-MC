% ------------------------------------------------------------------------------
%
%                           function newtonnu
%
%  this function solves keplers equation when the true anomaly is known.
%    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
%    the parabolic limit at 168 is arbitrary. the hyperbolic anomaly is also
%    limited. the hyperbolic sine is used because it's not double valued.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - fix small                                     24 sep 2002
%
%  inputs          description                    range / units
%    ecc         - eccentricity                   0.0  to
%    nu          - true anomaly                   -2pi to 2pi rad
%
%  outputs       :
%    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 deg
%    m           - mean anomaly                   0.0  to 2pi rad       151.7425 deg
%
%  locals        :
%    e1          - eccentric anomaly, next value  rad
%    sine        - sine of e
%    cose        - cosine of e
%    ktr         - index
%
%  coupling      :
%    arcsinh     - arc hyperbolic sine
%    sinh        - hyperbolic sine
%
%  references    :
%    vallado       2007, 85, alg 5
%
% [e0,m] = newtonnu ( ecc,nu );
% ------------------------------------------------------------------------------

function [e0,m] = newtonnu_vec ( ecc,nu )

    % ---------------------  implementation   ---------------------
    e0_error = 999999.9;
    m_error = 999999.9;
    small = 0.00000001;
    
    check_ecc = abs( ecc ) > small;
    idx_ecc = find(check_ecc);
    idx_notecc = find(~check_ecc);
    m = zeros(size(nu));
    e0 = zeros(size(nu));

    % --------------------------- circular ------------------------
    m(~check_ecc) = nu(~check_ecc);
    e0(~check_ecc) = nu(~check_ecc);
    % ---------------------- elliptical -----------------------
    check_e = ecc(idx_ecc)<(1.0-small);
    idx_e = idx_ecc(check_e);
    ecci = ecc(idx_e);
    nui = nu(idx_e);
    c_nui = cos(nui);
    one_p_eccicos = ( 1.0 +ecci.*c_nui );
    sine = ( sqrt( 1.0 -ecci.*ecci ) .* sin(nui) ) ./ one_p_eccicos;
    cose = ( ecci + c_nui ) ./ one_p_eccicos;
    e0_temp  = atan2( sine,cose );
    e0(idx_e) = e0_temp;
    m(idx_e)   = e0_temp - ecci.*sin(e0_temp);
    % -------------------- hyperbolic  --------------------
    check_h = ecc(idx_notecc)>(1.0-small);
    idx_h = idx_notecc(check_h);
    idx_p = idx_notecc(~check_h);
    ecch = ecc(idx_h);
    nuh = nu(idx_h);
    
    check_h2 = (ecch > 1.0 ) & (abs(nuh)+0.00001 < pi-acos(1.0./ecch));
    idx_h2 = idx_h(check_h2);
    ecch2 = ecch(check_h2);
    nuh2 = nuh(check_h2);
    
    sine = ( sqrt( ecch2.*ecch2-1.0  ) .* sin(nuh2) ) ./ ( 1.0  + ecch2.*cos(nuh2) );
    e0_temp  = asinh( sine );
    e0(idx_h2) = e0_temp;
    m(idx_h2)   = ecch2.*sinh(e0_temp) - e0_temp;
    e0(idx_h(~check_h2)) = e0_error;
    m(idx_h(~check_h2))  = m_error;
    % ----------------- parabolic ---------------------
    check_nu = ( abs(nu(idx_p)) < 168.0*pi/180.0  );
    e0_temp = tan( nu(idx_p(check_nu))*0.5  );
    e0(idx_p(check_nu)) = e0_temp;
    m(idx_p(check_nu)) = e0_temp + (e0_temp.^3)./3.0;
    e0(idx_p(~check_nu)) = e0_error;
    m(idx_p(~check_nu)) = m_error;
    
    check_ecc = ecc < 1.0;

    m_temp = rem( m((check_ecc)),2.0 *pi );
    check_m = m_temp < 0.0;
    m_temp(check_m) = m(check_m) + 2.0 *pi;
    m(check_ecc) = m_temp;
    e0(check_ecc) = rem( e0(check_ecc),2.0 *pi );