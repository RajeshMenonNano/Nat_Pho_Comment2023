% main file to generate encircled power

clear all; close all;

% constants
nm = 1e-9; 
um = 1e-6; 

% lens parameters:
lambda0 = 633*nm; 
D = 32*lambda0; % 20*um; 
NA = [0.7:0.01:0.99]; % change focal length inside function. 
F = 0.5*lambda0./NA; % FWHM
% Use line below if using radius in units of FWHM. 
% R = [0:0.01:10]; % radius of iris in unts of FWHM. 
% Use line below if using radius in units of lambda. 
R = [0:0.02:3]*lambda0; 
Encirc = zeros(length(NA), length(R)); 

for cnt1 = 1:length(NA)
    [Eout, total_transmitted_power, x_out, y_out] = plot_ideal_flat_lens_eff_function(lambda0, NA(cnt1), D, 0);
    d = sqrt(x_out.^2 + (y_out').^2); % Transpose y_out for correct dimension matching
    for cnt2 = 1:length(R)
        % Use line below if using radius in units of FWHM. 
        % Encirc(cnt1, cnt2) = sum(sum(abs(Eout.*(d  <= R(cnt2)*F(cnt1))).^2))/total_transmitted_power;  % if using units of FWHM. 
        % Use line below if using radius in units of lambda. 
        Encirc(cnt1, cnt2) = sum(sum(abs(Eout.*(d  <= R(cnt2))).^2))/total_transmitted_power;  % if using units of lambda0. 
    end
end

figure; imagesc(NA, R/lambda0, Encirc'); colorbar; xlabel('NA'); 
ylabel('Radius of iris (units of \lambda)'); 
% ylabel('Radius of iris (units of \lambda/2NA)'); 
axis xy; 
fontsize(gcf, 15,"points"); title('Encircled power'); 
%hold on; [C, h] = contour(Encirc, [0.9 0.9], 'k','LineWidth',2); clabel(C, h);
h_ax_c = axes('position', get(gca, 'position'), 'Color', 'none');
contour(h_ax_c, Encirc', [0.5 0.75 0.9], 'ShowText', 'on','LineWidth',2,'LineColor','k');
h_ax_c.Color = 'none'; h_ax_c.YTick = []; h_ax_c.XTick = [];
