% Idea flat lens encircled power calculation from:
% Sang, Di, et al. "Toward high‐efficiency ultrahigh numerical aperture freeform metalens: from vector diffraction theory to topology optimization." Laser & Photonics Reviews 16.10 (2022): 2200265.
function [Eout, total_transmitted_power, x_out, y_out] =  plot_ideal_flat_lens_eff_function(lambda0, NA, D, to_plot)

% define coordinates in lens plane
num_x_samples = 4096;
x1 = linspace(-10*D/2,10*D/2,num_x_samples);
y1 = x1;
x = meshgrid(x1,y1);
y = x';
f0 = D/2/tan(asin(NA)); % focal length

% phase function of ideal lens
r = sqrt(x.^2+y.^2);  % y, x are the cartesian coordinates. 
ph_lens = -2*pi/lambda0*(sqrt(r.^2 +f0^2)-f0); % r is the radial coordinate. 

theta = atan(r/f0);
% complex lens transmission function, including apodization. 
t_lens = 1./sqrt(cos(theta)).*exp(j*ph_lens); % note that there is a j missing in the equation in the paper. 
index = find(r>D/2); t_lens(index) = 0;

psi = atan(y./x);
% field transmitted after lens
Ex = t_lens.*(cos(theta).*cos(psi).^2 + sin(psi).^2);

Ey = t_lens.*(cos(theta)-1).*cos(psi).*sin(psi);

Ez = t_lens.*sin(theta).*cos(psi);

if to_plot == 1
figure(1); subplot(4,1,1); imagesc(x1, y1, abs(Ex)); axis equal; caxis([0 2]); colormap('cool'); colorbar; axis([-D/2 D/2 -D/2 D/2]);
subplot(4,1,2); imagesc(x1, y1, abs(Ey)); caxis([0 2]);  axis equal; colormap('cool'); colorbar;  axis([-D/2 D/2 -D/2 D/2]);
subplot(4,1,3); caxis([0 2]);  imagesc(x1, y1, abs(Ez)); axis equal;  caxis([0 2]); colormap('cool'); colorbar;  axis([-D/2 D/2 -D/2 D/2]);
subplot(4,1,4); imagesc(x1, y1, sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2)); axis equal;  caxis([0 2]); colormap('cool'); colorbar;  axis([-D/2 D/2 -D/2 D/2]);
end

[Ex_f, x_out, y_out] = angularSpectrumMethod(Ex, x1, y1, lambda0, f0); 
[Ey_f, x_out, y_out] = angularSpectrumMethod(Ey, x1, y1, lambda0, f0); 
Ez_f = calculateEz(Ex_f, Ey_f, Ez, x_out, y_out, f0); % based on divergence of E = 0. 
% normlization to make sure power flow in each plane is the same. 
A_const = sqrt( (sum(sum(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2)) - sum(sum(abs(Ex_f).^2 + abs(Ey_f).^2)))/sum(sum(abs(Ez_f).^2)) );
Ez_f = Ez_f*A_const;

Eout = sqrt(abs(Ex_f).^2 + abs(Ey_f).^2 + abs(Ez_f).^2);

if to_plot == 1
figure(2); subplot(4,1,1); imagesc( x_out/um, y_out/um, abs(Ex_f)); axis equal; colormap('cool'); colorbar; caxis([0 max(max(Eout))]); title('|Ex|');  axis([-2 2 -2 2]); % NA=0.8 % axis([-3/2 3/2 -3/2 3/2]); %NA=0.95  %axis([-D/2 D/2 -D/2 D/2]/um); 
subplot(4,1,2); imagesc( x_out/um, y_out/um, abs(Ey_f));  axis equal; colormap('cool'); colorbar;  caxis([0 max(max(Eout))]); title('|Ey|');  axis([-2 2 -2 2]); % NA=0.8; % axis([-3/2 3/2 -3/2 3/2]);  %NA=0.95  %axis([-D/2 D/2 -D/2 D/2]/um); 
subplot(4,1,3);  imagesc( x_out/um, y_out/um, abs(Ez_f)); axis equal;   colormap('cool'); colorbar;   caxis([0 max(max(Eout))]); title('|Ez|'); axis([-2 2 -2 2]); % NA=0.8 % axis([-3/2 3/2 -3/2 3/2]); %NA=0.95  %axis([-D/2 D/2 -D/2 D/2]/um); 
subplot(4,1,4); imagesc(x_out/um, y_out/um, Eout); axis equal;  colormap('cool'); colorbar;  caxis([0 max(max(Eout))]); title('|E|'); axis([-2 2 -2 2]); % NA=0.8 % axis([-3/2 3/2 -3/2 3/2]); %NA=0.95  %axis([-D/2 D/2 -D/2 D/2]/um); 
end

total_transmitted_power =  sum(sum(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2)); % this should be denominator for encircled power calculation. 

return;
