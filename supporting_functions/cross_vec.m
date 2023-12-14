function c = cross_vec(a,b)
% This function returns the cross product of the input vectors, observing
% the input ordering. Since the code does not have the error checking in
% the MATLAB function of the same name it is much faster
c = [a(:,2).*b(:,3)-b(:,2).*a(:,3),...
     a(:,3).*b(:,1)-b(:,3).*a(:,1),...
     a(:,1).*b(:,2)-b(:,1).*a(:,2)];
