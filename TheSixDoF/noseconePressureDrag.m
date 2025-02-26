%% Caleb Rice
% creating piecewise drag models for nose cone pressure drag

function dragFunc = noseconePressureDrag(radius, length, param, type)

    % Initialization
    a = sym("a");
    b = sym("b");
    M = sym("M");
    x = sym("x");

    fineness = length / (2 * radius);
    gamma = 1.4;

    if type == 2

        % Subsonic interpolation
        phi = acot(length / radius);
        C0 = 0.8 * (sin(phi))^2;
        C1 = sin(phi);
        Cprime1 = (4 / (gamma + 1)) * (1 - 0.5 * C1);

        sub_eq1 = a * (1^b) + C0 == C1;
        sub_eq2 = diff((a * (M^b) + C0), M) == Cprime1;
        sub_eq2 = subs(sub_eq2, M, 1);

        SubsonicCoeff = solve([sub_eq1 sub_eq2], [a, b]);

        SubCoeff(M) = SubsonicCoeff.a * (M^SubsonicCoeff.b) + C0;

        % supersonic equation
        SupCoeff(M) = (2.1 / (sqrt(1 + 4 * fineness^2))) + (0.5 / (sqrt((1 + 4 * fineness^2) * (M^2 - 1))));

        % transonic linear interpolation
        b0 = SubCoeff(1) - ((SupCoeff(1.3) - SubCoeff(1)) / 0.3);
        transCoeff(M) = ((SupCoeff(1.3) - SubCoeff(1)) / 0.3) * M + b0;

        % final function output
        dragFunc(M) = piecewise((M >= 0) & (M <= 1), SubCoeff, (M > 1) & (M < 1.3), transCoeff, (M >= 1.3), SupCoeff);
        
    elseif type == 3
        p = ((length^2 + radius^2) * (((2 - param) * length)^2 + (param * radius)^2)) / (4 * (param * radius)^2);
        slope(x) = ((length / param) - x) / sqrt(p - ((length / param) - x)^2);
        phi = acos(1 / (sqrt(1 + slope(length)^2)));

        C0 = 0.8 * (sin(phi))^2;
        C1 = sin(phi);
        Cprime1 = (4 / (gamma + 1)) * (1 - 0.5 * C1);

        sub_eq1 = a * (1^b) + C0 == C1;
        sub_eq2 = diff((a * (M^b) + C0), M) == Cprime1;
        sub_eq2 = subs(sub_eq2, M, 1);

        SubsonicCoeff = solve([sub_eq1 sub_eq2], [a, b]);

        SubCoeff(M) = SubsonicCoeff.a * (M^SubsonicCoeff.b) + C0;

        % supersonic equation
        SupCoeff(M) = (0.72 * (param-0.5)^2 + 0.82) * ((2.1 / (sqrt(1 + 4 * fineness^2))) + (0.5 / (sqrt((1 + 4 * fineness^2) * (M^2 - 1)))));

        % transonic linear interpolation
        b0 = SubCoeff(1) - ((SupCoeff(1.3) - SubCoeff(1)) / 0.3);
        transCoeff(M) = (0.72 * (param-0.5)^2 + 0.82) * (((SupCoeff(1.3) - SubCoeff(1)) / 0.3) * M + b0);

        % final function output
        dragFunc(M) = piecewise((M >= 0) & (M <= 1), SubCoeff, (M > 1) & (M < 1.3), transCoeff, (M >= 1.3), SupCoeff);

    elseif type == 4
    elseif type == 5
    elseif type == 6
    elseif type == 1
        fprintf("WARNING: Haack series pressure calculations are based on data for a Von Karman nose cone, which may be inaccurate for parameters other than 0")
    else
        fprintf("ERROR: value of nc_type must be an integer between 1 and 6")
        dragFunc = "ERROR";
    end
end