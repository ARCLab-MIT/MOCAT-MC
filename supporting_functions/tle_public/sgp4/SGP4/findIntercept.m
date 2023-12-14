function llh_test=findIntercept(h_intercept,u_sat,rng_sat,origin_llh,varargin)
%        llh_test=findIntercept(h_intercept,u_sat,origin_llh)
%  Find llh coordinates of h_intercept along ray to satellite defined by u_sat
%  
%  h_intercept  = desired intercept height
%  u_sat        = unit vector in tcs system pointing to satellite
%  origin_llh   = tcs origin (llh)
%
if isempty(varargin)
    err=100;
else
    varargin{1}=err;
end
h_test=0;
rng_test_min=0;
rng_test_max=rng_sat;
while abs(h_test-h_intercept)>err
    rng_test=(rng_test_min+rng_test_max)/2;
    tcs_test=u_sat*rng_test;
    llh_test=tcs2llhT(tcs_test,origin_llh);
    h_test=llh_test(3);
    if h_test>=h_intercept;
        rng_test_max=rng_test;
    else
        rng_test_min=rng_test;
    end
end
return