function R = rot2D(theta)
% Standard 2D rotation using 3D function
R = [cos(theta) -sin(theta) 0; ...
     sin(theta)  cos(theta) 0; ...
             0           0  1];
R = R(1:2,1:2);

