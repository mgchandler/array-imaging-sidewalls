function fn_print_time(process, time)
% Simply prints the time provided to the screen, in format:
%       "{Process} completed in {time}"
% INPUTS:
% - process : string
%       Name of the process which the time is associated with.
% - time : double
%       Time taken to run the process. Expected to be provided in seconds.

if time > 3600
    time = time/3600;
    units = 'hrs';
elseif time > 60
    time = time/60;
    units = 'mins';
else
    units = 'secs';
end

fprintf('%s in %.2f %s\n', process, time, units)

end