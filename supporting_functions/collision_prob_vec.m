function  pr=collision_prob_vec(p1_radius, p1_v, p2_radius, p2_v, CUBE_RES)
%p1_v and p2_v should be [n_sats x 3]

sigma = (p1_radius + p2_radius).^2*(pi/1e6); % crossectional collision area
dU = CUBE_RES^3; % volume of the cube
Vimp = sqrt(sum((p1_v - p2_v).^2,2)); % rel. velocity
pr = Vimp/dU .* sigma;