% https://en.wikipedia.org/wiki/NACA_airfoil
% This function takes input of a 4-digit NACA airfoil identifier and
% and outputs a plot of the airfoil and .dat file in Selig format.
% The airfoil datapoints are spaced using a cosine distribution to group
% points around the LE and TE where pressure gradients are highest.
%_________________________________________________________________________
% FORMAT
% n4dig('ID',1 or 0)
%
% Example: n4dig('2412',1) plots a Naca 2412 airfoil and writes the points
% to a .dat file named NACA4Airfoil.dat. If the second argument == 0, then no file
% is written.
%_________________________________________________________________________
function [s,output] = n4dig(ID,writefile);
    
    % INPUTS -------------------------------------------------------------  
    c = 1; %non-dimensional chord length
    s = ID; %numeric (MUST BE CHARACTER ARRAY) ID of 4-digit section (i.e. '2412')
    res = 200; %resolution of plot (200 default)
    % --------------------------------------------------------------------
    
    %% nargin statements
    if ischar(s) == 0
        error("Airfoil ID must be input as character array (Ex: '2412')");
    end
    
    if nargin == 0
        error('At least one argument is required (airfoil ID)');
    end
    if nargin == 1
        output = 0;
    end
    if nargin == 2 %if 2 input arguments use input data for export yes or no
        output = writefile;
    end
    
  
    %% Main Calcs
    d1 = str2double(s(1)); %pulls the first digit
    d2 = str2double(s(2)); %pulls the second digit
    d34 = str2double(s(3:4)); %pulls the third and fourth digits
    m = d1/100; %max camber (as percentage of chord length)
    p = d2/10; %max camber loc. from LE (as percentage of 1/10 of chord)
    t = d34/100; %max thickness ratio (as percentage of chord)
    
    % EXAMPLE AIRFOIL
    % NACA 2412
    % First Digit: 2/100 = 2% max camber
    % Second Digit: 4/10 = 40% from LE is max camber
    % Third/Fourth Digit: 12/100 = 12% max thickness
    
    
    % this portion adapted from: https://www.mathworks.com/matlabcentral/answers/276773-need-help-with-4-digit-airfoil
    % ----------------------------------------------------------------
    
    %x = linspace(0, c, res);
    x = nonLinspace(0, c, res, 'cos'); %cosine dist. function by Connor Ott
    % Cosine distribution of function y = 1/2(1-cos(x))
    
    if d1 == 0 %if symmetric 4-digit airfoil
        yt = 5*t*(.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 +...
            0.2843*x.^3 - 0.1015*x.^4);
        xu = x;
        xl = x;
        yu = yt;
        yl = -yt;
        
    else
        yt = 5*t*c*(.2969*(sqrt(x/c))+-.1260*(x/c)+...
        -.3516*(x/c).^2+.2843*(x/c).^3+-.1015*(x/c).^4);
        for i = 1:length(x)
              if x(i) <= p*c
                  yc(i) = m*(x(i)/p^2)*(2*p-(x(i)/c));
                  dx(i) = (2*m)/p^2*(p-(x(i)/c));
              elseif x(i) > p*c
                  yc(i) = m*((c-x(i))/(1-p)^2)*(1+(x(i)/c)-(2*p));
                  dx(i) = ((2*m)/(1-p)^2)*(p-(x(i)/c));
              end
              %upper and lower limits of the airfoil (xu,yu) ; (xl,yl)
              theta = atan(dx(i));
              xu(i) = x(i)-yt(i)*sin(theta);
              yu(i) = yc(i)+yt(i)*cos(theta);
              xl(i) = x(i)+yt(i)*sin(theta);
              yl(i) = yc(i)-yt(i)*cos(theta);
        end    
    end
    
    %% Plot Airfoil Shape
    plot(xu,yu, 'Color', '#0072BD', 'LineWidth', 1.25) %upper surface
    hold on
    plot(xl,yl, 'Color', '#0072BD', 'LineWidth', 1.25) %bottom surface
    if d1 ~= 0
        plot(x,yc,'r') %camber line
    end
    plot([0 c],[0 yt(1,res)],'--') %chord line
    axis equal; grid on;
    title(ID)
    xlabel('Chord Length')
    
    hold off
    % ----------------------------------------------------------------
    
    %% Output file
    if output == 1  
        write_ID = 'NACA4Airfoil.dat';
        
        if exist(write_ID, 'file') == 2 %delete .dat file if same exists
            delete(write_ID);
        end
        
        %write header
        header = ['NACA', convertCharsToStrings(s)];
        writematrix(header, write_ID,'Delimiter',' ');
        
        %selig_out = [1 0]; %manually write tail point
        % CAUSES EXCESSIVE PANEL ANGLE ON TE ^^^
        selig_out = [];
        for i = res-1:-1:2 %write upper surface, tail to nose
            selig_out = [selig_out; xu(i), yu(i)];
        end
        for i = 1:res %write lower surface
            selig_out = [selig_out; xl(i), yl(i)];
        end
        %selig_out = [selig_out; 1,0]; %manually write last 1,0 point
        % CAUSES EXCESSIVE PANEL ANGLE ON TE ^^^        
        writematrix(selig_out, write_ID,'Delimiter',' ','WriteMode','append');        
    end
       
end
function [ nonLinVec ] = nonLinspace( mn, mx, num, spacetype )

%(https://www.mathworks.com/matlabcentral/fileexchange/64831-non-linearly-spaced-vector-generator), 
% MATLAB Central File Exchange. Retrieved April 6, 2021
% -------------------------------------------------------------------------
% nonLinspace(mn, mx, num, spacetype) returns a vector of non-linearly 
% spaced elements based on spacing specified by spacetype. 
%
% nonLinVec = nonLinspace(mn, mx, num, 'exp10') returns a vector of
% elements with smaller spacing at the beginning of the vector and greater
% spacing at the end of the vector based on the curve y = 10^x.
%
% nonLinVec = nonLinspace(mn, mx, num, 'cos') returns a vector of elements
% with smaller spacing at the beginning and end of the vector, and greater
% spacing in the middle based on the curve y = 1/2(1-cos(x)).
%
% nonLinVec = nonLinspace(mn, mx, num, 'log10') returns a vector of
% elements with greater spacing at the beginning of the vector and smaller
% spacing at the end of the vector. 
% 
%   Inputs: 
%       mn        - The minimum value in the vector. 
%       mx        - The maximum value in the vector.
%       num       - The number of elements in the vector. 
%       spacetype - Specifies the type of spacing needed. 
%
%   Outputs:
%       nonLinVec - A vector consisting of elements with spacing specified 
%                   by spacetype.
% -------------------------------------------------------------------------
if strcmpi(spacetype, 'exp10')
    % exponentialliness is the upper bound of the original 10^x curve
    % before it is scaled to fit the limits requested by the user. Since
    % the concavity of 10^x changes in different parts of its domain,
    % different spacing is seen when using different bounds. After some
    % basic qualitative analysis, an exponentialliness of 20 seemed to be a
    % good fit for my purposes. Increasing this value will increase the
    % spacing towards the end of the vector and decrease it towards the
    % beginning. 
    exponentialliness = 20;
    nonLinVec = (mx-mn)/exponentialliness * ...
                (10.^(linspace(0, log10(exponentialliness+1), num)) - 1)...
                + mn;
            
elseif strcmpi(spacetype, 'cos')
    nonLinVec = (mx - mn)*(0.5*(1-cos(linspace(0, pi, num)))) + mn;
    
elseif strcmpi(spacetype, 'log10')
    % As with exponentialliness, this defines the bounds on the log10(x)
    % curve. Increasing loginess will decreasing the spacing towards the
    % end of the vector and increase it towards the beginning. 
    loginess = 1.5;
    nonLinVec = (mx - mn)/loginess* ...
                log10((linspace(0, 10^(loginess) - 1, num)+ 1)) + mn;
            
end
    
end