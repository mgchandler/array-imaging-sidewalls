function [new_data, new_tx, new_rx] = fn_swap_FMC_HMC(data, tx, rx)
% Converts between FMC data and HMC data, depending on which was supplied.
% INPUTS:
%   - data : array of size (time_pts, N)
%       Array of timetraces. Assumes that first axis is time, second is
%       tx-rx pairs. For FMC, this should equal (probe_els^2); for HMC this
%       should equal (Σ_(n=1)^(probe_els) (n)).
%   - tx : array of size (1, N)
%       Specifies the transmitting element index in the phased array for
%       each timetrace in data.
%   - rx : array of size (1, N)
%       Specifies the receiving element index in the phased array for each
%       timetrace in data.
%
% OUTPUTS:
%   - new_data : array of size(time_pts, N_new)
%       Array of timetraces, equivalent to data but now has the other form
%       (i.e. if data is FMC, new_data is HMC and vice versa).
%   - new_tx : array of size (1, N_new)
%       Specifies transmitting element index for each timetrace in
%       new_data.
%   - new_rx : array of size (1, N_new)
%       Specifies receiving element index for each timetrace in new_data.

probe_els = max(tx);
assert(probe_els == max(rx), "fn_swap_FMC_HMC: invalid tx/rx indices.")
assert(size(tx, 2) == size(rx, 2), "fn_swap_FMC_HMC: tx, rx are not the same length.")
N_HMC = sum(1:probe_els);
N_FMC = probe_els^2;

if size(data, 2) == N_FMC
    output_FMC = false;
elseif size(data, 2) == N_HMC
    output_FMC = true;
else
    error("fn_swap_FMC_HMC: data is not FMC or HMC. Check tx/rx indices.")
end

if output_FMC
    %% Input is HMC; output FMC.
    new_data = zeros(size(data, 1), N_FMC);
    new_tx   = zeros(1, N_FMC);
    new_rx   = zeros(1, N_FMC);
    % Treat diagonals and off-diagonals separately.
    diags = tx==rx;
    new_data(:, 1:probe_els) = data(:, diags);
    new_tx(:, 1:probe_els)   = tx(:, diags);
    new_rx(:, 1:probe_els)   = rx(:, diags);
    
    % Copy off-diags to both upper and lower sections of FMC.
    new_data(:, probe_els+1:N_HMC) = data(:, ~diags);
    new_tx(:, probe_els+1:N_HMC)   = tx(:, ~diags);
    new_rx(:, probe_els+1:N_HMC)   = rx(:, ~diags);
    new_data(:, N_HMC+1:end) = data(:, ~diags);
    new_tx(:, N_HMC+1:end)   = rx(:, ~diags);
    new_rx(:, N_HMC+1:end)   = tx(:, ~diags);
    
    % Sort data so it looks as expected - first by rx then by tx.
    [~,I] = sort(new_rx);
    new_data = new_data(:, I);
    new_tx   = new_tx(:, I);
    new_rx   = new_rx(:, I);
    [~,I] = sort(new_tx);
    new_data = new_data(:, I);
    new_tx   = new_tx(:, I);
    new_rx   = new_rx(:, I);
    
else 
    %% Input is FMC; output HMC.
    new_data = zeros(size(data, 1), N_HMC);
    new_tx   = zeros(1, N_HMC);
    new_rx   = zeros(1, N_HMC);
    % Can either outright delete one entry, or take average of two entries.
    % Take average here.
    diags = tx==rx;
    lower = tx<rx;
    upper = tx>rx;
    data_upper = data(:, upper);
    tx_upper   = tx(:, upper);
    rx_upper   = rx(:, upper);
    data_lower = data(:, lower);
    tx_lower   = tx(:, lower);
    rx_lower   = rx(:, lower);
    
    new_data(:, 1:probe_els) = data(:, diags);
    new_tx(:, 1:probe_els)   = tx(:, diags);
    new_rx(:, 1:probe_els)   = rx(:, diags);
    
    % Arrange the upper and lower off-diag sections in such a way that
    % timetrace (tx, rx) in the upper section is equivalent to timetrace
    % (rx, tx) in the lower section (i.e. first timetrace has (tx=1,rx=2)
    % in lower, so first timetrace in upper should be (tx=2,rx=1)). Can
    % then average these two matrices directly.
    [~,I] = sort(tx_upper);
    data_upper = data_upper(:, I);
    tx_upper   = tx_upper(:, I);
    rx_upper   = rx_upper(:, I);
    [~,I] = sort(rx_upper);
    data_upper = data_upper(:, I);
    tx_upper   = tx_upper(:, I);    % Not required as lower used later,
    rx_upper   = rx_upper(:, I);    % but computed for completeness.
    [~,I] = sort(rx_lower);
    data_lower = data_lower(:, I);
    tx_lower   = tx_lower(:, I);
    rx_lower   = rx_lower(:, I);
    [~,I] = sort(tx_lower);
    data_lower = data_lower(:, I);
    tx_lower   = tx_lower(:, I);
    rx_lower   = rx_lower(:, I);
    
    % Do averaging, take lower (tx,rx) values for consistency with BRAIN.
    new_data(:, probe_els+1:end) = (data_upper + data_lower) / 2;
    new_tx(:, probe_els+1:end)   = tx_lower;
    new_rx(:, probe_els+1:end)   = rx_lower;
    
    % Sort.
    [~,I] = sort(new_rx);
    new_data = new_data(:, I);
    new_tx   = new_tx(:, I);
    new_rx   = new_rx(:, I);
    [~,I] = sort(new_tx);
    new_data = new_data(:, I);
    new_tx   = new_tx(:, I);
    new_rx   = new_rx(:, I);
end

end