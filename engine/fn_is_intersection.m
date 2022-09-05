function is_intersect = fn_is_intersection(x1, x2, y1, y2)
% Tests whether the line segments x and y (with end points x1, x2 and y1,
% y2 respectively) intersect anywhere. The proceedure is as follows:
% START
% (1) : do bounding boxes intersect?
%       y -> 2;     n -> 5
% (2) : does line x intersect line segment y?
%       y -> 3;     n -> 5
% (3) : does line y intersect line segment x?
%       y -> 4;     n -> 5
% (4) : lines do intersect, return 1
% (5) : lines do not intersect, return 0
%
% INPUTS:
%   - x1 : array (1, 3)
%       The cartesian coordinates which describe one end of the line
%       segment labelled x.
%   - x2 : array (1, 3)
%       The cartesian coordinates which describe the other end of the line
%       segment labelled x.
%   - y1 : array (1, 3)
%       The cartesian coordinates which describe one end of the line
%       segment labelled y.
%   - y2 : array (1, 3)
%       The cartesian coordinates which describe the other end of the line
%       segment labelled y.
%
% OUTPUTS:
%   - is_intersect : logical
%       Logical value for whether there is an intersection between line
%       segments x and y (i.e. 1 if there is an intersection, 0 otherwise).

[no_scats, ~] = size(y1); % one scatterer:3x1

% Check bounding boxes. General line L has corners
%   [[L1_x, L1_z];
%    [L1_x, L2_z];
%    [L2_x, L1_z];
%    [L2_x, L2_z]]
box_x_TL = [min(x1(1), x2(1)), ...
            0.0, ...
            min(x1(3), x2(3))];
box_x_BR = [max(x1(1), x2(1)), ...
            0.0, ...
            max(x1(3), x2(3))];
box_y_TL = [min(y1(:, 1), y2(:, 1)), ...
            zeros(no_scats, 1), ...
            min(y1(:, 3), y2(:, 3))];
box_y_BR = [max(y1(:, 1), y2(:, 1)), ...
            zeros(no_scats, 1), ...
            max(y1(:, 3), y2(:, 3))];
is_overlap = and(and(box_x_TL(1) <= box_y_BR(:, 1), ...
                     box_x_BR(1) >= box_y_TL(:, 1)), ...
                 and(box_x_TL(3) <= box_y_BR(:, 3), ...
                     box_x_BR(3) >= box_y_TL(:, 3)) ...
             );
if all(~is_overlap)
    is_intersect = is_overlap;
    return
end

where_is_y1_wrt_x = point_loc_wrt_line(x1, x2, y1);
where_is_y2_wrt_x = point_loc_wrt_line(x1, x2, y2);
where_is_x1_wrt_y = point_loc_wrt_line(y1, y2, x1);
where_is_x2_wrt_y = point_loc_wrt_line(y1, y2, x2);

is_y1_on_x = abs(where_is_y1_wrt_x) < eps;
is_y2_on_x = abs(where_is_y2_wrt_x) < eps;
is_x1_on_y = abs(where_is_x1_wrt_y) < eps;
is_x2_on_y = abs(where_is_x2_wrt_y) < eps;

is_y1_right_of_x = where_is_y1_wrt_x < 0;
is_y2_right_of_x = where_is_y2_wrt_x < 0;
is_x1_right_of_y = where_is_x1_wrt_y < 0;
is_x2_right_of_y = where_is_x2_wrt_y < 0;

touches_or_crosses_x_y = or(xor(is_y1_right_of_x, is_y2_right_of_x), ...
                             or(is_y1_on_x, is_y2_on_x));
touches_or_crosses_y_x = or(xor(is_x1_right_of_y, is_x2_right_of_y), ...
                             or(is_x1_on_y, is_x2_on_y));
                         
is_intersect = and(and(touches_or_crosses_x_y, touches_or_crosses_y_x), is_overlap);



end

function point_loc = point_loc_wrt_line(a1, a2, y)
% Works out where point y is with respect to line a described by start and
% end location a1 and a2 respectively.
%
% INPUTS:
%   - a1 : array (1, 3)
%       Start of line segment a.
%   - a2 : array (1, 3)
%       End of line segment a.
%   - y : array (1, 3)
%       A point which is to be determined whether it is on line a.
%
% OUTPUTS:
%   - r : double
%       Cross product value which indicates where y is wrt a.
%       r < 0  =>  y right of a
%       r = 0  =>  y on a
%       r > 0  =>  y left of a

[no_scats_y, ~] = size(y);
[no_scats_a, ~] = size(a1);

atemp = [a2(:, 1) - a1(:, 1), ...
         zeros(no_scats_a, 1), ...
         a2(:, 3) - a1(:, 3)];
ytemp = [y(:, 1) - a1(:, 1), ...
         zeros(no_scats_y * no_scats_a, 1), ...
         y(:, 3) - a1(:, 3)];
     
point_loc = atemp(:, 1) .* ytemp(:, 3) - ytemp(:, 1) .* atemp(:, 3);
end