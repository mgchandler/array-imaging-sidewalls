p1 = [1, 0, 1];   % Point 1
p2 = [-3, 0, 2];  % Point 2
r = 6;            % Radius
N = 100;          % Number of points to discretise arc into
go_left = false;  % At point 1 looking towards point 2, which way does the
                  % arc go?
                  
geom = fn_make_arc("C", p1, p2, r, N, go_left, 0);
                  
% %% Arc
% 
% if r < norm(p2-p1)/2
%     error("Radius too small to construct an arc.")
% end
% 
% % Points 1 and 2 sit on the arc by definition. Find the centre point using
% % the radius.
% v_2 = (p2 - p1)/2;   % Half of the chord which joins point 1 to point 2.
% % Vector from mid-point between p1 and p2 to the centre of the circle. Must
% % be normal to v_2. Will be on LHS of v_2 if go_left == true.
% n = (-1)^go_left * sqrt(r^2 - norm(v_2)^2) * [-v_2(3)/norm(v_2), 0, v_2(1)/norm(v_2)];
% 
% % Absolute location of arc centre point. n1 defined wrt v_2; v_2 defined
% % wrt p1.
% d = p1 + v_2 + n;
% 
% % Defining phi = 0 along the +ve x-axis, what is the angle made by p1 from
% % the centre of the circle.
% phi1 = atan2(p1(3)-d(3), p1(1)-d(1));
% % Angle made by p2 from circle centre.
% phi2 = atan2(p2(3)-d(3), p2(1)-d(1));
% phi = linspace(phi1, phi2, N+1);
% 
% % Make the wall.
% wall = [real(r*exp(1i*phi))+d(1); zeros(1, N+1); imag(r*exp(1i*phi))+d(3)].';
% basis = zeros(N+1, 3, 3);
% basis(:, 1, 1) = sin(phi);
% basis(:, 1, 3) = cos(phi);
% basis(:, 2, 2) =  1;
% basis(:, 3, 1) = -basis(:, 1, 3);
% basis(:, 3, 3) =  basis(:, 1, 1);
% 
% %% Linear
% 
% geom = fn_make_wall('a', p1, p2, N, 1);
% 
% %% Plotting
% 
figure, grid on, hold on, axis equal
plot(p1(1), p1(3), '*', "MarkerSize", 10)
plot(p2(1), p2(3), '*', "MarkerSize", 10)
plot(geom.coords(:, 1), geom.coords(:, 3))

M = 60;
plot([geom.coords(M, 1), geom.coords(M, 1)+geom.basis(M, 1, 1)], [geom.coords(M, 3), geom.coords(M, 3)+geom.basis(M, 3, 1)])
plot([geom.coords(M, 1), geom.coords(M, 1)+geom.basis(M, 1, 3)], [geom.coords(M, 3), geom.coords(M, 3)+geom.basis(M, 3, 3)])