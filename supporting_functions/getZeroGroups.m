function [g1,g2,g3] = getZeroGroups(inmatsat)
idx_mass = 8; idx_radius = 9; idx_objectclass = 23;

% mass vs radius
g1.allclass = []; % group 1: payloads; logical index
g2.allclass = []; % group 2: RBs
g3.allclass = []; % group 3: all debris

for ii = 1:12
    msinds = inmatsat(:,idx_objectclass) == ii; % all obj w/ objclass
    if ii == 1
        g1.allclass = find(msinds);                                 % all payload entries
        g1.zr = g1.allclass(inmatsat(g1.allclass,idx_radius) == 0); % index of zero radius
        g1.zm = g1.allclass(inmatsat(g1.allclass,idx_mass) == 0);   % index of zero mass
        g1.nz = setdiff(g1.allclass, union(g1.zr,g1.zm));           % index of non-zero radius and mass
        g1.nzno = g1.nz(~isoutlier(inmatsat(g1.nz,idx_radius)) ...  % non-zero, non-outlier
            & ~isoutlier(inmatsat(g1.nz,idx_mass)) );
    elseif ii == 5
        g2.allclass = find(msinds);
        g2.zr = g2.allclass(inmatsat(g2.allclass,idx_radius) == 0); % all RB entries
        g2.zm = g2.allclass(inmatsat(g2.allclass,idx_mass) == 0);
        g2.nz = setdiff(g2.allclass, union(g2.zr,g2.zm));
        g2.nzno = g2.nz(~isoutlier(inmatsat(g2.nz,idx_radius)) ...
            & ~isoutlier(inmatsat(g2.nz,idx_mass)) );
    else
        g3.allclass = union(g3.allclass, find(msinds));             % all debris entries
    end
end
g3.zr = g3.allclass(inmatsat(g3.allclass,idx_radius) == 0);
g3.zm = g3.allclass(inmatsat(g3.allclass,idx_mass) == 0);
g3.nz = setdiff(g3.allclass, union(g3.zr, g3.zm));
g3.nzno = g3.nz(~isoutlier(inmatsat(g3.nz,idx_radius)) ...
    & ~isoutlier(inmatsat(g3.nz,idx_mass)) );
end
