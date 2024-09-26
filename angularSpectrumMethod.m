function [Uout, x_out, y_out] = angularSpectrumMethod(Uin, x, y, lambda, z)
    % Uin: Input complex field distribution
    % x, y: Spatial cartesian coordinates for the input field
    % lambda: Wavelength of the light
    % z: Propagation distance
    % Uout: Output complex field distribution
    % x_out, y_out: Spatial cartesian coordinates for the output field
    
    % Calculate spatial frequencies
    [Nx, Ny] = size(Uin);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    fx = (-Nx/2:Nx/2-1) / (Nx*dx);
    fy = (-Ny/2:Ny/2-1) / (Ny*dy);
    [FX, FY] = meshgrid(fx, fy);
    
    % Angular spectrum
    k = 2*pi/lambda;  % Wavenumber
    FS = fftshift(fft2(Uin));  % Fourier transform of the input field
    H = exp(1i * k * z * sqrt(1 - (lambda * FX).^2 - (lambda * FY).^2));  % Transfer function
    Uout = ifft2(ifftshift(FS .* H));  % Inverse Fourier transform
    
    % Output coordinates
    x_out = x(1) + (0:Nx-1) * dx;
    y_out = y(1) + (0:Ny-1) * dy;
end
