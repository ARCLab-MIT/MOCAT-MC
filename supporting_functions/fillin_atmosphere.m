% ATMOSPHERIC MODEL (pre-computed JB2008)
% valid from March 2008 - Feb 2224 for density file 'dens_jb2008_032020_022224.mat'
density_profile = 'JB2008'; % The options are 'static' or 'JB2008'
cfgMC.density_profile = density_profile;
if strcmpi(cfgMC.density_profile,'JB2008')
    cfgMC = initJB2008(cfgMC);
end

function cfgMCout = initJB2008(cfgMC)
    
    fn = which('dens_jb2008_032020_022224.mat');
    load(fn);
    
    dens_times=zeros(length(dens_highvar.month),1);
    
    for k=1:length(dens_highvar.month)
        dens_times(k,1)= juliandate(datetime(dens_highvar.year(k),dens_highvar.month(k),0));
    end
    
    [dens_times2,alt2] = meshgrid(dens_times,dens_highvar.alt);

    cfgMCout = cfgMC;
    cfgMCout.param.dens_times = dens_times2;
    cfgMCout.param.alt = alt2;
    cfgMCout.param.dens_value = dens_highvar.dens;
end