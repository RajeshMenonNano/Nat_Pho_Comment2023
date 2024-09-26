function Ez_out = calculateEz(Ex, Ey, Ez, x, y, f)
    % calculateEz - Calculate the z-component of the electric field
    %   Ex: x-component of the electric field (2D matrix)
    %   Ey: y-component of the electric field (2D matrix)
    %   x: x-coordinates (vector)
    %   y: y-coordinates (vector)
    
    % Check if dimensions of Ex and Ey match
    if size(Ex) ~= size(Ey)
        error('Dimensions of Ex and Ey must match.');
    end
    
    % Calculate the spatial derivatives of Ex and Ey
    % We provide two spacing arguments for the gradient function,
    % corresponding to the x and y dimensions.
    [dEx_dx, ~] = gradient(Ex, x(2)-x(1), y(2)-y(1)); % derivative of Ex with respect to x
    [~, dEy_dy] = gradient(Ey, x(2)-x(1), y(2)-y(1)); % derivative of Ey with respect to y

    % Calculate Ez
    Ez_out = Ez - f*(dEx_dx + dEy_dy);
end
