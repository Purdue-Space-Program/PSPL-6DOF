lift = computeLift(rocket, mach, alpha);

% Lift is pretty easy (approximately).

% mach greater than 1, use the linear supersonic theory:
if mach > 1
    lift = cL; % Assign the calculated lift for supersonic conditions
end
cL = 4*alpha / sqrt(mach^2 - 1);
