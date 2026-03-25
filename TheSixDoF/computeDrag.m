function drag = computeDrag(rocket, rho, V, mu)

% split the rocket up into the components:

noseLength = rocket.NoseLength

totLength = rocket.TotalLength;

area = pi*(rocket.OuterDiameter/2)^2

% drag calculations:

Re = rho*V*L / mu;

Cf = 0.074



end
