function [r,v]=oe2rv_vec(input_all,big_e,param)

mu = param.mu;

a = input_all(:,1);
ecc = input_all(:,2);
inc = input_all(:,3);
big_omega = input_all(:,4);
w = input_all(:,5);

if isempty(big_e)
    % Solve Kepler's Equation
    big_m = input_all(:,6);
    big_e = big_m;
    check_mag = 1:numel(big_e);
    for l=1:1000
        Epw_check = big_e(check_mag);
        e_check = ecc(check_mag);
        DeltaEpw = - (big_m(check_mag) - Epw_check + e_check.*sin(Epw_check))./(-1+e_check.*cos(Epw_check));
        big_e(check_mag) = Epw_check+DeltaEpw;
        check_mag = check_mag(DeltaEpw > 10e-14);
        if isempty(check_mag) %if all elements have gone below threshold
            break
        end
    end
end

% Get Initial Position and Velocity
c_bige = cos(big_e);
s_bige = sin(big_e);
one_m_e_sq = 1-ecc.^2;
c_bigO = cos(big_omega);
s_bigO = sin(big_omega);
c_i = cos(inc);
s_i = sin(inc);
c_w = cos(w);
s_w = sin(w);

rmag=a.*(1-ecc.*c_bige);
rp=a.*(c_bige-ecc);
rq=a.*sqrt(one_m_e_sq).*s_bige;
vp=-sqrt(mu*a)./rmag.*s_bige;
vq=sqrt(mu*a.*one_m_e_sq)./rmag.*c_bige;
c11=c_bigO.*c_w-s_bigO.*s_w.*c_i;
c12=-c_bigO.*s_w-s_bigO.*c_w.*c_i;
c21=s_bigO.*c_w+c_bigO.*s_w.*c_i;
c22=-s_bigO.*s_w+c_bigO.*c_w.*c_i;
c31=s_w.*s_i;
c32=c_w.*s_i;
r1=c11.*rp+c12.*rq;
r2=c21.*rp+c22.*rq;
r3=c31.*rp+c32.*rq;
v1=c11.*vp+c12.*vq;
v2=c21.*vp+c22.*vq;
v3=c31.*vp+c32.*vq;
r=[r1,r2,r3];
v=[v1,v2,v3];

