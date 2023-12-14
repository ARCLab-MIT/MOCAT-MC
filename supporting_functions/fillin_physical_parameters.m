% Fill in missing satellite physical parameters
cfgMC.fillMassRadius = 2;     % 0: don't fill in missing DISCOS data (many objects with 0 radius and/or mass)
                              % 1: ESA's method -- assume spherical aluminum depending on RCS size (S/M/L) 
                              % 2: resampling method
cfgMC.initpopMultiplier = 1;  % multiplier for initial population 
cfgMC.physicalBstar = 1;        % [0/1] recalculate B* as Bstar = 1/2 * Cd * A/m * 0.157e6
                                % calc in initSim; from setupTLEssem_uniform.m