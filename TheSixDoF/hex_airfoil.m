function [xy, params] = hex_airfoil(varargin)
% HEX_AIRFOIL  Generate, plot, and export a hexagonal airfoil.
%
%   [XY, PARAMS] = HEX_AIRFOIL('Chord', 1.0, 'TOverC', 0.12, ...
%       'FlatRatio', 0.40, 'FlatCenter', 0.50, 'Npts', 201, ...
%       'NormalizeToChord', false, 'Plot', true)
%
% Parameters (Name-Value pairs):
%   'Chord'           : Chord length c (default 1.0).
%   'TOverC'          : Maximum thickness as a fraction of chord (default 0.12).
%                       Max thickness occurs on the flat region; total thickness = TOverC * Chord.
%   'FlatRatio'       : Flat length as a fraction of chord (0 < FlatRatio < 1). (default 0.40)
%                       This is the length of the flat segment on both top and bottom.
%   'FlatCenter'      : Location of the *center* of the flat as a fraction of chord (0..1). (default 0.50)
%                       0 = leading edge, 1 = trailing edge. Flat is centered here.
%   'Npts'            : Total number of output points including the repeated TE point (default 201).
%                       Points start at TE, go over the upper surface to LE, then lower surface back to TE.
%   'NormalizeToChord': If true, writes x/c and y/c to the .dat file (default false -> absolute units).
%   'Plot'            : If true, plots the airfoil (default true).
%
% Returns:
%   XY     : [N x 2] array of (x,y) coordinates ordered for XFOIL (TE->upper->LE->lower->TE).
%   PARAMS : Struct of resolved parameters.
%
% Geometry model:
%   - Sharp LE (x=0) and TE (x=c).
%   - Three segments per surface:
%       LE wedge: linear rise from thickness 0 at x=0 to +/-(t/2) at x = x1
%       Flat    : constant +/-(t/2) between x1 and x2
%       TE wedge: linear drop from +/-(t/2) at x = x2 to 0 at x = c
%   - x1 and x2 are set by FlatRatio and FlatCenter.
%
% Notes:
%   - XFOIL ordering: starts at TE (x=c, y=0), counter-clockwise, ends at TE.
%   - If FlatRatio ~ 0, the shape becomes a diamond / double wedge (no flat). If ~1, it becomes a rectangle
%     with vanishing wedges—both are handled with guards.
%
% Example calls:
%   hex_airfoil('Chord',1,'TOverC',0.12,'FlatRatio',0.40,'FlatCenter',0.50,...
%               'Npts',241);
%   hex_airfoil('Chord',0.5,'TOverC',0.18,'FlatRatio',0.25,'FlatCenter',0.45,...
%               'NormalizeToChord',true);
%   [xy, P] = hex_airfoil('TOverC',0.10,'FlatRatio',0.60,'Npts',301,'Plot',false);
%

% -------- Parse inputs --------
p = inputParser;
addParameter(p, 'Chord', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'TOverC', 0.12, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p, 'FlatRatio', 0.40, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p, 'FlatCenter', 0.50, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p, 'Npts', 201, @(x)isnumeric(x)&&isscalar(x)&&x>=8&&mod(round(x),1)==0);
addParameter(p, 'NormalizeToChord', false, @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p, 'Plot', true, @(x)islogical(x)||ismember(x,[0 1]));
parse(p, varargin{:});
P = p.Results;

c     = P.Chord;
t     = P.TOverC * c;
FR    = min(max(P.FlatRatio,0),1);
xctr  = P.FlatCenter * c;
tol   = 1e-9;

% -------- Point budgeting (per-surface count) --------
N_total = round(P.Npts);
N_side  = max(floor((N_total - 1)/2) + 1, 6);   % points along TE->LE (upper), incl both ends

% -------- Build per-case geometry --------
if FR <= tol
    % Logic for DOUBLE-WEDGE (diamond)
    xpk = max(0, min(c, xctr));           % apex location
    le_len = xpk;                          % 0 -> xpk
    te_len = c - xpk;                      % xpk -> c
    prop   = [te_len, le_len] / max(te_len + le_len, eps);  % [TE wedge, LE wedge]
    counts = max([2 2], round(prop * (N_side + 1)));        % +1 slack
    target = N_side + 1;
    while sum(counts) ~= target
        diffCt = sum(counts) - target;
        if diffCt > 0
            [~, idx] = max(counts); if counts(idx) > 2, counts(idx) = counts(idx) - 1; end
        else
            [~, idx] = max(prop); counts(idx) = counts(idx) + 1;
        end
    end
    n_te = counts(1); n_le = counts(2);

    % Upper (TE -> apex -> LE)
    xu_te = linspace(c, xpk, n_te);
    yu_te = (t/2) * (c - xu_te) / max(c - xpk, eps);

    xu_le = linspace(xpk, 0, n_le);
    yu_le = (t/2) * (xu_le) / max(xpk, eps);

    xu = [xu_te, xu_le(2:end)];
    yu = [yu_te, yu_le(2:end)];

    % Lower (LE -> apex -> TE)
    xl_le = linspace(0, xpk, n_le);
    yl_le = -(t/2) * (xl_le) / max(xpk, eps);

    xl_te = linspace(xpk, c, n_te);
    yl_te = -(t/2) * (c - xl_te) / max(c - xpk, eps);

    xl = [xl_le, xl_te(2:end)];
    yl = [yl_le, yl_te(2:end)];

    x1 = xpk; x2 = xpk; flatL = 0;

elseif FR >= 1 - tol
    % Logic for RECTANGULAR
    % Upper: straight line y = +t/2 from TE->LE
    xu = linspace(c, 0, N_side);
    yu = (t/2) * ones(size(xu));

    % Lower: straight line y = -t/2 from LE->TE
    xl = linspace(0, c, N_side);
    yl = -(t/2) * ones(size(xl));

    x1 = 0; x2 = c; flatL = c;

else
    % Logic for HEX (flat + wedges)
    flatL = FR * c;
    x1 = max(0, min(c, xctr - flatL/2));
    x2 = max(0, min(c, xctr + flatL/2));

    le_len   = x1 - 0;
    flat_len = x2 - x1;
    te_len   = c - x2;
    total_side = le_len + flat_len + te_len;
    if total_side <= 0, error('Degenerate geometry.'); end

    prop   = [te_len, flat_len, le_len] / total_side;
    counts = max([2 2 2], round(prop * (N_side + 2)));
    target = N_side + 2;
    while sum(counts) ~= target
        diffCt = sum(counts) - target;
        if diffCt > 0
            [~, idx] = max(counts);
            if counts(idx) > 2, counts(idx) = counts(idx) - 1; end
        else
            [~, idx] = max(prop); counts(idx) = counts(idx) + 1;
        end
    end
    n_te = counts(1); n_fl = counts(2); n_le = counts(3);

    % Upper (TE -> flat -> LE)
    xu_te = linspace(c, x2, n_te);
    yu_te = (t/2) * (c - xu_te) / max(c - x2, eps);

    xu_fl = linspace(x2, x1, n_fl);
    yu_fl = (t/2) * ones(size(xu_fl));

    xu_le = linspace(x1, 0, n_le);
    yu_le = (t/2) * (xu_le / max(x1, eps));

    xu = [xu_te,  xu_fl(2:end),  xu_le(2:end)];
    yu = [yu_te,  yu_fl(2:end),  yu_le(2:end)];

    % Lower (LE -> flat -> TE)
    xl_le = linspace(0, x1, n_le);
    yl_le = -(t/2) * (xl_le / max(x1, eps));

    xl_fl = linspace(x1, x2, n_fl);
    yl_fl = -(t/2) * ones(size(xl_fl));

    xl_te = linspace(x2, c, n_te);
    yl_te = -(t/2) * (c - xl_te) / max(c - x2, eps);

    xl = [xl_le,  xl_fl(2:end),  xl_te(2:end)];
    yl = [yl_le,  yl_fl(2:end),  yl_te(2:end)];
end

% -------------------- Assemble full XFOIL ordering --------------------
% Start at TE (upper path already starts at TE) -> ... -> LE,
% then lower path from LE to TE (skip duplicate LE), end at TE (duplicate allowed).
if FR >= 1 - tol
    % Rectangular: include both LE and TE vertical edges in the file
    % Upper runs: (c,+t/2) -> (0,+t/2)
    % Lower runs: (0,-t/2) -> (c,-t/2)
    x = [xu, xl];                 % KEEP xl(1) so we have (0,-t/2)
    y = [yu, yl];

    % close at TE by returning to the first point (c,+t/2)
    x(end+1) = c; 
    y(end+1) = +t/2;
else
    % Hex / double-wedge: skip duplicate LE point as usual
    x = [xu, xl(2:end)];
    y = [yu, yl(2:end)];
    if x(end) ~= x(1) || y(end) ~= y(1)
        x(end+1) = x(1); 
        y(end+1) = y(1);
    end
end
xy = [x(:), y(:)];

% Optional normalization for file export
xy_out = xy;
if P.NormalizeToChord
    xy_out = xy / c;
end

% -------------------- Plot --------------------
if P.Plot
    figure('Color',[1 1 1]); hold on;
    plot(xy(:,1), xy(:,2), '-', 'LineWidth', 1.25);
    plot([0 c],[0 0],'--');
    axis equal; grid on;
    xlabel('x'); ylabel('y');
    if FR <= tol
        ttl = sprintf('Double-Wedge  (c=%.3g, t/c=%.3g, apex=%.0f%% chord)', c, P.TOverC, (xctr/c)*100);
    elseif FR >= 1 - tol
        ttl = sprintf('Rectangular  (c=%.3g, t/c=%.3g)', c, P.TOverC);
    else
        ttl = sprintf('Hex  (c=%.3g, t/c=%.3g, flat=%.0f%% @ %.0f%% chord)', c, P.TOverC, (flatL/c)*100, (xctr/c)*100);
    end
    title(ttl);
end

% -------------------- Write .dat file --------------------
fid = fopen('HEXAirfoil.dat', 'w');
if fid < 0, error('Could not open HEXAirfoil.dat for writing.'); end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%.8f %.8f\n', xy_out.');

% -------------------- Output params --------------------
params = struct();
params.Chord            = c;
params.TOverC           = P.TOverC;
params.MaxThickness     = t;
params.FlatRatio        = FR;
params.FlatLength       = flatL;
params.FlatCenter       = P.FlatCenter;
params.Npts             = size(xy,1);
params.NormalizeToChord = logical(P.NormalizeToChord);
params.Case             = ternary(FR<=tol,'double-wedge', ternary(FR>=1-tol,'rectangular','hex'));
params.FlatX1           = x1;    % for double-wedge x1=x2=xpk; for rect x1=0, x2=c
params.FlatX2           = x2;
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end