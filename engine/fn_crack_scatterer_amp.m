function [scat_amp, extras] = fn_crack_scatterer_amp(inc_angle, out_angle, dens, vel_L, vel_T, freq, inc_mode, out_mode, crack_length, nodes_per_wavelength)
% Computes the scattering amplitude amplitude for a crack inside the block
% (i.e. no interaction with the surface) for some incident or scatterer
% amplitude.

% INPUTS:
% - inc_angle : double
%       Incident angle on the scatterer.
% - out_angle : double
%       Outgoing angle from the scatterer.
% - dens : double
%       Density of the material.
% - vel_L : double
%       Longitudinal wave speed in the material.
% - vel_T : double
%       Shear wave speed in the material.
% - freq : double
%       Frequency of the wave.
% - inc_mode : logical
%       Logical test of whether the incident mode is shear.
% - out_mode : logical
%       Logical test of whether the outgoing mode is shear.
% - scat_info : struct (1, 1)
%       Contains all information pertaining to the scatterer.
% - crack_length : double
%       Length of the crack from tip to tip.
% - nodes_per_length : integer
%       Defaults to 20.
%
% OUTPUTS:
% - scat_amp : complex
%       Scattering amplitude for a crack.
% - extras : idk
%       I can't remember why I included this. Trace it back and check.
%       Might not be needed.

lame_lambda = dens * (vel_L^2 - 2 * vel_T^2);
lame_mu = dens * vel_T^2;

omega = 2 * pi * freq;
xi1 = omega / vel_L;
xi2 = omega / vel_T;
xi = vel_T / vel_L;

lambda_L = vel_L / freq;
lambda_T = vel_T / freq;
extras.a = 1;

num_nodes = ceil(crack_length / lambda_L * nodes_per_wavelength);
p = 0.1133407986; % Magic constant.
h_nodes = crack_length / (num_nodes + 2*p);
x_nodes = zeros(num_nodes, 1);
x_nodes(1) = h_nodes * (1/2 + p);
x_nodes(2:num_nodes) = x_nodes(1) + [1:num_nodes-1] * h_nodes;
x_nodes = x_nodes - crack_length / 2;
x_nodes = x_nodes.';

mx_nodes = x_nodes .* ones(1,num_nodes) - ones(num_nodes,1) .* x_nodes.';    
A_x = fn_Qx_matrix(xi,xi2,h_nodes,num_nodes,mx_nodes);
A_z = fn_Qz_matrix(xi,xi2,h_nodes,num_nodes,mx_nodes);
       
a_L = -1i * xi1 * pi / xi2^2;
a_T = -1i * xi2 * pi / xi2^2;

scat_amp = zeros(length(inc_angle), length(out_angle));

for ii = 1:length(inc_angle)

    nv = [0; 1];
    sv = [-sin(inc_angle(ii)); -cos(inc_angle(ii))];
    tv = [sv(2); -sv(1)];

    if ~inc_mode
        b_L = exp(1i * xi1 * x_nodes * sv(1)) * basis_func(-xi1 * h_nodes * sv(1)); % Define basis function.
        b_x = -2 * sv(1) * sv(2) * b_L;
        b_z = -(1 / xi^2 - 2 * sv(1)^2) * b_L;
        b_x = b_x.';
        b_z = b_z.';
        vxL = linsolve(A_x, b_x); % Solve A_x, b_x
        vzL = linsolve(A_z, b_z); % Solve A_z, b_z
    
    else
        b_T = exp(1i * xi2 * x_nodes * sv(1)) * basis_func(-xi2 * h_nodes * sv(1));
        b_x = -(tv(1) * sv(2) + tv(2) * sv(1)) * b_T;
        b_z = -2 * tv(2) * sv(2) * b_T;
        b_x = b_x.';
        b_z = b_z.';
        vxT = linsolve(A_x, b_x); % Solve A_x, b_x
        vzT = linsolve(A_z, b_z); % Solve A_z, b_z  
        
    end
    
    for jj = 1:length(out_angle)

        ev = [sin(out_angle(jj)); cos(out_angle(jj))];
        tv = [ev(2); -ev(1)];
        c_L = basis_func(xi1 * h_nodes * ev(1)) * exp(-1i * xi1 * ev(1) * x_nodes);
        c_L = c_L.';
        c_T = basis_func(xi2 * h_nodes * ev(1)) * exp(-1i * xi2 * ev(1) * x_nodes);
        c_T = c_T.';

        if ~inc_mode % If incident mode is longitudinal.
            v_L = [a_L * (vxL.' * c_L); a_L * (vzL.' * c_L)];
            v_T = [a_L * (vxL.' * c_T); a_L * (vzL.' * c_T)];
            if ~out_mode % If outgoing mode is longitudinal.
                scat_amp(ii, jj) = ( ...
                    1/4 * sqrt(2/pi) * exp(-1i * pi/4) * xi1^(5/2) * ( ...
                    lame_lambda / (dens * omega^2) * dot(v_L, nv) + ...
                    2 * lame_mu / (dens * omega^2) * dot(v_L, ev) * dot(ev.', nv) ...
                    ) / sqrt(lambda_L) ...
                );

            else % If incident mode is longitudinal and outgoing mode is shear.
                scat_amp(ii, jj) = ( ...
                    1/4 * sqrt(2/pi) * exp(-1i * pi/4) * xi2^(5/2) * lame_mu / (dens * omega^2) * ( ...
                    dot(v_T, tv) * dot(ev.', nv) + dot(v_T, ev) * dot(tv.', nv) ...
                    ) / sqrt(lambda_T) ...
                );

            end

        else % If incident mode is shear.
            v_L = [a_T * (vxT.' * c_L); a_T * (vzT.' * c_L)];
            v_T = [a_T * (vxT.' * c_T); a_T * (vzT.' * c_T)];
            
            if ~out_mode % If outgoing mode is longitudinal.
                scat_amp(ii, jj) = ( ...
                    1/4 * sqrt(2/pi) * exp(-1i * pi/4) * xi1^(5/2) * ( ...
                    lame_lambda / (dens * omega^2) * dot(v_L, nv) + ...
                    2 * lame_mu / (dens * omega^2) * dot(v_L, ev) * dot(ev.', nv) ...
                    ) / sqrt(lambda_L) ...
                );

            else % If incident mode is shear and outgoing mode is shear.
                scat_amp(ii, jj) = ( ...
                    1/4 * sqrt(2/pi) * exp(-1i * pi/4) * xi2^(5/2) * lame_mu / (dens * omega^2) * ( ...
                    (v_T.' * tv) * (ev.' * nv) + (v_T.' * ev) * (tv.' * nv) ...
                    ) / sqrt(lambda_T) ...
                );

            end

        end
    end

end

end



function F = basis_func(k)
k00 = 1e-1;
ind0 = find(abs(k)<=k00); k0 = k(ind0);
F(ind0) = 32/35-16/315*k0.^2+4/3465*k0.^4;
ind_0 = find(abs(k)>k00); k_0 = k(ind_0);
F(ind_0) = 96./k_0.^7 .* ( k_0 .* (k_0.^2 - 15) .* cos(k_0) - (6 * k_0.^2-15) .* sin(k_0) );
F = 35/32 * F;
return;
end



function A_x = fn_Qx_matrix(xi,xi2,h_nodes,N_nodes,mx_nodes)

 %integral around branch point
 z1 = xi2*h_nodes*[0:N_nodes-1];
   
  for ii=1:length(z1)
   int_F1 = @(x) fn_F(xi,xi2,h_nodes,1-x.^2) ./ sqrt(2-x.^2) .* cos((1-x.^2)*z1(ii));
   
   int_F2 = @(x) fn_F(xi,xi2,h_nodes,1+x.^2) ./ sqrt(2+x.^2) .* cos((1+x.^2)*z1(ii));
   
   int_F = @(x) ( -((1+1i*x.^2).^2 - 0.5).^2 .* fn_P(xi2*h_nodes*(1+1i*x.^2)) ./ sqrt(2+1i*x.^2) * exp(-1i*pi/4) * exp(1i*(z1(ii)-2*xi2*h_nodes)) + ...
                 (xi+1i*x.^2).^2 .* fn_P(xi2*h_nodes*(xi+1i*x.^2)) .* sqrt(2*xi+1i*x.^2) .* x.^2 *  exp(1i*pi/4) * exp(1i*xi*(z1(ii)-2*xi2*h_nodes)) ) .*...
                  exp(-(z1(ii)-2*xi2*h_nodes)*x.^2);
   
  
   if ii<3
      I12(ii) =  4 * 1i * (quad(int_F1,0,1)) +  4 * quad(int_F2,0,50); 
   else
      I12(ii) = 4 * 1i * quad(int_F,0,70);  
   end
   
  end
  
  I12 = [I12(N_nodes:-1:2),I12]; 
  v_ind = [1:N_nodes]';
  m_ind = N_nodes*ones(size(mx_nodes)) + v_ind*ones(1,N_nodes) - ones(N_nodes,1)*v_ind';
  A_x = I12(m_ind);
return;
end



function A_z = fn_Qz_matrix(xi,xi2,h_nodes,N_nodes,mx_nodes)

  %integral around branch point
  z1 = xi2*h_nodes*[0:N_nodes-1];
  
  for ii=1:length(z1)
   int_F1 = @(x) fn_F(xi,xi2,h_nodes,xi-x.^2) ./ sqrt(2*xi-x.^2) .* cos((xi-x.^2)*z1(ii));
      
   int_F2 = @(x) fn_F(xi,xi2,h_nodes,xi+x.^2) ./ sqrt(2*xi+x.^2) .* cos((xi+x.^2)*z1(ii));
   
   int_F = @(x) ( -((xi+1i*x.^2).^2 - 0.5).^2 .* fn_P(xi2*h_nodes*(xi+1i*x.^2)) ./ sqrt(2*xi+1i*x.^2) * exp(-1i*pi/4) * exp(1i*xi*(z1(ii)-2*xi2*h_nodes)) + ...
                 (1+1i*x.^2).^2 .* fn_P(xi2*h_nodes*(1+1i*x.^2)) .* sqrt(2+1i*x.^2) .* x.^2 *  exp(1i*pi/4) * exp(1i*(z1(ii)-2*xi2*h_nodes)) ) .*...
                  exp(-(z1(ii)-2*xi2*h_nodes)*x.^2);
   
   
   if ii<3
      I12(ii) =  4 * 1i * (quad(int_F1,0,sqrt(xi))) +  4 * quad(int_F2,0,50); 
   else
      I12(ii) = 4 * 1i * quad(int_F,0,70);  
   end
   
   
  end

  I12 = [I12(N_nodes:-1:2),I12]; 
  v_ind = [1:N_nodes]';
  m_ind = N_nodes*ones(size(mx_nodes)) + v_ind*ones(1,N_nodes) - ones(N_nodes,1)*v_ind';
  A_z = I12(m_ind);
return;
end



function f = fn_P(k)
  k00 = 1e-1;
  ind0 = find(abs(k)<=k00); k0 = k(ind0);
  F(ind0) = 32/35 + 32/35*1i*k0 - 32/63*k0.^2 - 64/315*1i*k0.^3;    
  ind1 = find(abs(k)>k00); k1 = k(ind1);
  sk = (exp(2*1i*k1)-1) / (2*1i); ck = (exp(2*1i*k1)+1) / 2;
  F(ind1) = 96./k1.^7 .* ( k1 .* (k1.^2 - 15) .* ck - (6 * k1.^2-15) .* sk );
  F = 35/32 * F;
  
  f = F.^2;
return;
end



function f = fn_F(xi,xi2,h,betta)
  k = xi2*h*betta;
  
  k00 = 1e-1;
  ind0 = find(abs(k)<=k00); k0 = k(ind0);
  F(ind0) = 32/35-16/315*k0.^2+4/3465*k0.^4;
  ind_0 = find(abs(k)>k00); k_0 = k(ind_0);
  F(ind_0) = 96./k_0.^7 .* ( k_0 .* (k_0.^2 - 15) .* cos(k_0) - (6 * k_0.^2-15) .* sin(k_0) );
  F = 35/32 * F;
  
  sigma1 = fn_sigma(betta, xi); sigma2 = fn_sigma(betta, 1);
  L2 = ( -(betta.^2 - 0.5).^2 + betta.^2 .* sigma1 .* sigma2 );
  
  f = L2 .* F.^2;
return;
end



function sigma = fn_sigma(k, k0)
 s = k.^2 - k0.^2;
 ind_pos = find(s>=0); ind_neg = find(s<0);
 sigma = zeros(size(k));
 sigma(ind_pos) = sqrt(s(ind_pos)); sigma(ind_neg) = -1i*sqrt(-s(ind_neg)); 
return;
end