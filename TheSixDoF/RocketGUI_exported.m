classdef RocketGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RocketDesignTab                matlab.ui.container.Tab
        GridLayout9                    matlab.ui.container.GridLayout
        GUIHelpButton                  matlab.ui.control.Button
        SetBasePathButton              matlab.ui.control.Button
        PropGrid2                      matlab.ui.container.GridLayout
        Tree                           matlab.ui.container.Tree
        ListBox                        matlab.ui.control.ListBox
        Switchto3D                     matlab.ui.control.Button
        UpdatePlotButton               matlab.ui.control.Button
        Switchto2D                     matlab.ui.control.Button
        Panel                          matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        TabGroup2                      matlab.ui.container.TabGroup
        RocketTab                      matlab.ui.container.Tab
        RocketGrid                     matlab.ui.container.GridLayout
        WetMasskgEditField             matlab.ui.control.NumericEditField
        WetMasskgEditFieldLabel        matlab.ui.control.Label
        RocketColor                    matlab.ui.control.ColorPicker
        ColorColorPickerLabel          matlab.ui.control.Label
        ButtonPanel                    matlab.ui.container.Panel
        ButtonGrid                     matlab.ui.container.GridLayout
        SaveRocketButton               matlab.ui.control.Button
        LoadRocketButton               matlab.ui.control.Button
        RASAeroDataLabel               matlab.ui.control.Label
        AeroDataButton                 matlab.ui.control.Button
        RocketNameEditField            matlab.ui.control.EditField
        rocket_real_name               matlab.ui.control.EditField % idk if i did ts right - david
        launch_site_name               matlab.ui.control.EditField % idk if i did ts right - david
        NameEditField_2Label           matlab.ui.control.Label
        NoseConeLengthmEditField       matlab.ui.control.NumericEditField
        NoseConeLengthmEditFieldLabel  matlab.ui.control.Label
        NoseConeGeometryDropDown       matlab.ui.control.DropDown
        NoseConeGeometryDropDownLabel  matlab.ui.control.Label
        RocketDiameterEditField        matlab.ui.control.NumericEditField
        RocketOuterDiametermLabel      matlab.ui.control.Label
        RocketLengthEditField          matlab.ui.control.NumericEditField
        TotalLengthmLabel              matlab.ui.control.Label
        ComponentsTab                  matlab.ui.container.Tab
        ComponentsGrid                 matlab.ui.container.GridLayout
        PropertyPanel                  matlab.ui.container.Panel
        PropertyGrid                   matlab.ui.container.GridLayout
        AddComponentButton             matlab.ui.control.Button
        ComponentSelectionDropDown     matlab.ui.control.DropDown
        ComponentSelectionDropDownLabel  matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
        SimulationTab                  matlab.ui.container.Tab
        LaunchLocationPanel            matlab.ui.container.Panel
        Panel_2                        matlab.ui.container.Panel
        GridLayout2                    matlab.ui.container.GridLayout
        GetWeatherConditionsButton     matlab.ui.control.Button
        SimulateLaunchButton           matlab.ui.control.Button
        LongitudedegEditField          matlab.ui.control.NumericEditField
        LongitudedegEditFieldLabel     matlab.ui.control.Label
        LatitudedegEditField           matlab.ui.control.NumericEditField
        LatitudedegEditFieldLabel      matlab.ui.control.Label
        DateSelectionDatePicker        matlab.ui.control.DatePicker
        DateSelectionDatePickerLabel   matlab.ui.control.Label
    end


    properties (Access = private)
        lineColor % line color for plots
        ThreeDPlot=0; % flag to turn 3d plotting on and off
        %ComponentList = "" % a list of the components on the rocket
        ComponentDetails
        rocket
        PropertyEditFields
        PropertyEditLabels
        autoRefresh = 0;
        rootNode;
        env
        Settings
        Date datetime
        AeroLoc
    end

    methods (Access = private)

        function RocketPlotter(app)
            % create general parameters for the rocket:
            leng = app.RocketLengthEditField.Value;
            noseLeng = app.NoseConeLengthmEditField.Value;
            dia = app.RocketDiameterEditField.Value;
            R = dia/2;

            % if the plot is in 3d
            if app.ThreeDPlot


                view(app.UIAxes, 3) %set view to 3D
                cla(app.UIAxes) %reset axis
                cameratoolbar("show");
                cameratoolbar("SetMode","orbit");
                cameratoolbar("SetCoordSys", "none")

                % body
                [Z, Y, X] = cylinder(R,100); %make unit cyliner along x axis
                X_body = X*(leng-noseLeng) + noseLeng;

                surf(app.UIAxes, X_body,Y,Z, "FaceColor",app.RocketColor.Value,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                axis(app.UIAxes, "equal")
                axis(app.UIAxes, 'auto')
                hold(app.UIAxes, "on")

                % nose
                resolution = 100;
                x_res_nose = 0:noseLeng/resolution:noseLeng;
                switch app.NoseConeGeometryDropDown.Value

                    case 'Conic'
                        nose_radius_func_ish = R.*(x_res_nose./noseLeng);
                    case 'Tangent Ogive'
                        L = noseLeng;
                        rho = (R^2 + L^2) / (2*R);
                        nose_radius_func_ish = sqrt(rho^2-(L-x_res_nose).^2) + R - rho;
                    case 'Von Karman'
                        theta = linspace(0,pi,resolution);
                        L = noseLeng;
                        xNose = L/2 * (1-cos(theta));
                        y_nose = R/sqrt(pi) * sqrt(theta-sin(2*theta)/2);
                        nose_radius_func_ish = interp1(xNose,y_nose,x_res_nose);
                    case 'Elliptical'
                        nose_radius_func_ish = sqrt((R^2) -(R^2).*((x_res_nose-noseLeng).^2)./(noseLeng^2));
                end
                [Z, Y, X] = cylinder(nose_radius_func_ish, 100);
                X_nose = X*noseLeng;

                surf(app.UIAxes, X_nose,Y,Z, "FaceColor",app.RocketColor.Value,'FaceAlpha', 0.7, 'EdgeAlpha',0);

                % wireframe
                thet = 0:2*pi/resolution:2*pi;
                y = R*sin(thet);
                z = R*cos(thet);
                x = noseLeng.*thet.^0;
                plot3(app.UIAxes, x,y,z, 'w')
                plot3(app.UIAxes, x+leng-noseLeng, y, z, 'w')

                % plot the fins:
                app.PlotFins();

                % components:
                app.PlotComponents();


            else % plot is in 2D

                view(app.UIAxes, 2)
                cla(app.UIAxes)
                cameratoolbar("SetMode","dollyhv");
                %cameratoolbar("hide");

                % define the geometry over the nose cone:
                xNose = linspace(0,noseLeng, 50);

                % change the y profile based on the selection.
                switch app.NoseConeGeometryDropDown.Value

                    case 'Conic'
                        yNose = xNose.*dia./(noseLeng*2);
                    case 'Tangent Ogive'
                        L = noseLeng;
                        rho = (R^2 + L^2) / (2*R);
                        yNose = sqrt(rho^2-(L-xNose).^2) + R - rho;

                    case 'Von Karman'
                        theta = linspace(0,pi,50);
                        L = noseLeng;
                        xNose = L/2 * (1-cos(theta));
                        yNose = R/sqrt(pi) * sqrt(theta-sin(2*theta)/2);

                    case 'Elliptical'
                        L = noseLeng;
                        yNose = R*sqrt(1-((xNose-L).^2./L^2));

                end

                x = [noseLeng,leng];
                y = dia* ones(1,numel(x));
                x = [xNose,x];
                y = [2*yNose,y];

                % plot the base body of the rocket:
                plot(app.UIAxes, x,y/2, app.lineColor)
                hold(app.UIAxes, "on")
                plot(app.UIAxes, x,-y/2, app.lineColor)
                plot(app.UIAxes, [x(end),x(end)], [dia/2,-dia/2], app.lineColor)


                % plot the fins of the rocket
                app.PlotFins()

                app.PlotComponents();

                % define the standard limits for the plot
                hold(app.UIAxes, 'off')
                xlim(app.UIAxes, [-0.02*leng,leng*1.15])
                axis(app.UIAxes, "equal")

            end

        end

        function Geoplotter(app)

            delete(findall(app.LaunchLocationPanel, 'Type', 'axes'))

            lat = app.LatitudedegEditField.Value;
            long = app.LongitudedegEditField.Value;

            g = geoaxes(app.LaunchLocationPanel);

            set(g, 'Basemap', 'satellite')

            geoplot(g, lat, long, 'ro')

            size = 0.1;

            geolimits(g, [lat-size, lat+size], [long-size,long+size]);
        end

        function [finPtsX, finPtsY] = FinPlotter(app, nFields)

            % This function is only called when a fin plot is present, pull
            % the parameters from here:

            delete(findall(app.PropertyGrid, 'Type', 'axes'))
            % Create an axes for the plot (this will take the last row)
            plotAxes = axes(app.PropertyGrid);
            plotAxes.Layout.Row = nFields+1;  % Put plot in the last row (nFields + 1)
            plotAxes.Layout.Column = [1, 2];    % Span across both columns

            % plot the fins:
            span = app.PropertyEditFields(3).Value;
            rootChord = app.PropertyEditFields(4).Value;
            tipChord = app.PropertyEditFields(5).Value;
            sweep = app.PropertyEditFields(6).Value;

            % the fin will always be defined by four points, with the first
            % point begin [0,0]. These are defined in the order:
            % [0,0]
            % [sweep,span]
            % [sweep+tipChord, span]
            % [rootChord, 0]

            finPtsX = [0, sweep, sweep+tipChord, rootChord, 0];
            finPtsY = [0, span, span, 0, 0];

            plot(plotAxes, finPtsX, finPtsY, app.lineColor);
            axis(plotAxes, "equal");
        end

        function PlotFins(app)

            % look at the rocket object and see if there are fins
            % associated with it:

            % first, check if the rocket is empty. If not, proceed:
            if ~isempty(app.rocket)
                compList = app.rocket.ComponentList;

                % in the event that the rocket has components:
                if numEntries(compList) > 0
                    exists = 0;
                    count = 0;
                    len = numEntries(compList);
                    values = compList.values;

                    % go through each and check for fins
                    for idx = 1:len
                        if isa(values{idx}, 'Fins')
                            exists = 1;
                            count = count + 1;
                            finObject = values{idx};
                            numFins(count) = double(finObject.Count);
                            rootChord(count) = finObject.RootChord;
                            tipChord(count) = finObject.TipChord;
                            span(count) = finObject.Span;
                            sweep(count) = finObject.Sweep;
                            finOffset(count) = finObject.Position(1);
                            finColor = finObject.Color;
                        end
                    end

                    if exists == 0
                        return
                    end
                else
                    return
                end
            else
                return
            end

            % 2d plot

            for idx = 1:count

                if ~app.ThreeDPlot

                    % replace this (get rid of catch statement)
                    try
                        finPtsX = [0, sweep(idx), sweep(idx)+tipChord(idx), rootChord(idx), 0];
                    catch
                        return
                    end
                    finPtsY = [0, span(idx), span(idx), 0, 0];

                    % get parameters from user input:
                    rearDist = finOffset(idx);
                    dia = app.rocket.OuterDiameter;
                    leng = app.rocket.TotalLength;

                    % first, check which fins should be plotted based on occlusion
                    % (manually for now lmao, don't know how to write this
                    % programmatically)

                    switch numFins(idx)
                        case 1
                            plotFin = 1;
                        case 2
                            plotFin = [1,2];
                        case 3
                            plotFin = [1,2,3];
                        case 4
                            plotFin = [1,2,4];
                        case 5
                            plotFin = [1,2,5];
                        case 6
                            plotFin = [1,2,6];
                    end


                    % run a for loop for each of the fins, calculated the projected
                    % view
                    for idx2 = plotFin

                        % the first fin always points up towards us, so use
                        % that as the baseline reference (theta = 0):
                        theta = (2*pi)/numFins(idx) * (idx2-1);

                        % generate an array of matrices for the projected fins.
                        xFinShifted = finPtsX + rearDist;
                        yFinProjection = finPtsY*sin(theta);

                        % add the radial component to the y-values:
                        rad = dia/2 * sin(theta);

                        yFinProjection = yFinProjection+rad;

                        % if the y-component is less than the rocket body
                        % radius, it is occluded from the view. Calculate how much
                        % of the fin to show in this case:

                        % start by calculating the slope of the front and back
                        % lines:
                        slopeFront = (yFinProjection(2)-yFinProjection(1)) ...
                            / (xFinShifted(2) - xFinShifted(1));

                        slopeBack = (yFinProjection(2)-yFinProjection(1)) ...
                            / (xFinShifted(3) - xFinShifted(4));

                        % in the region from the top to the back of the view:
                        if theta > pi/2 && theta <= pi
                            % figure out the distance for each of the y-points:
                            dist = dia/2 - yFinProjection;

                            % update the x-position based on the y-points if the
                            % distance is less than zero (inside the body), shift
                            % the points by the slope:

                            yFinProjection(yFinProjection<dia/2) = dia/2;

                            % in the region from the back to the bottom side:
                        elseif theta > pi && theta < 3*pi/2
                            yFinProjection(yFinProjection>-dia/2) = -dia/2;
                        end

                        plot(app.UIAxes, xFinShifted, yFinProjection, app.lineColor);

                    end

                else

                    if numFins(idx) ~= 0

                        % fins

                        % create rectangle until introduce naca airfoils for fins
                        xOut = [0,0,1,1,0];
                        yOut = [-0.01,0.01,0.01,-0.01,-0.01];
                        %assumes y comes pre-scaled

                        X_fin = zeros(length(xOut),2);
                        Y_fin = zeros(length(xOut),2);
                        Z_fin = zeros(length(xOut),2);
                        X_fin_top = zeros(length(xOut));
                        Y_fin_top = zeros(length(xOut));
                        Z_fin_top = zeros(length(xOut));

                        dia = app.RocketDiameterEditField.Value;
                        R = dia/2;

                        for n = 1:length(xOut)
                            X_fin(n,1) = (finOffset(idx)) + rootChord(idx).*xOut(n);
                            X_fin(n,2) = (finOffset(idx)) + tipChord(idx).*xOut(n) + sweep(idx);
                            Y_fin(n,1) = yOut(n);
                            Y_fin(n,2) = yOut(n);
                            Z_fin(n,1) = R;
                            Z_fin(n,2) = R + span(idx);
                        end

                        for n = 1:length(xOut)
                            X_fin_top(n) = (finOffset(idx)) + tipChord(idx).*xOut(n) + sweep(idx);
                            Y_fin_top(n) = yOut(n);
                            Z_fin_top(n) = R+span(idx);
                        end

                        scopy = zeros(numFins(idx));
                        stcopy = zeros(length(xOut), numFins(idx));

                        for i = 1:numFins(idx)
                            scopy(i) = surf(app.UIAxes,X_fin,Y_fin,Z_fin, "FaceColor",finColor,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                            stcopy(:,i) = fill3(app.UIAxes,X_fin_top,Y_fin_top,Z_fin_top, [1,1,1], 'FaceColor',finColor, 'EdgeAlpha', 0);
                            direction = [1 0 0];
                            origin = [0 0 0];
                            rotate(scopy(i),direction,rad2deg((i-1)*(2*pi)/numFins(idx)),origin);
                            rotate(stcopy(:,i), direction,rad2deg((i-1)*(2*pi)/numFins(idx)),origin)

                        end
                    end
                end
            end
        end

        function PlotComponents(app)

            % get the components from the rocket

            if ~isempty(app.rocket)
                compList = app.rocket.ComponentList;

                % in the event that the rocket has components:
                if numEntries(compList) > 0
                    len = numEntries(compList);
                    values = compList.values;

                    % go through each and check for fins
                    for idx = 1:len
                        % first case is a tank
                        if isa(values{idx}, 'Tank')

                            tankObj = values{idx};

                            leng = tankObj.Length;
                            rad = tankObj.TankDia/2;
                            FuelOx = tankObj.FuelOx;
                            dist = tankObj.Position(1);
                            if strcmp(FuelOx, 'Fuel')
                                color = 'r';
                            else
                                color = 'b';
                            end

                            if app.ThreeDPlot

                                isSpherica = false;
                                if isSpherica
                                    [Z, Y, X] = cylinder(rad,100); %make unit cyliner along x axis
                                    X_body = X*(leng-2*rad)+dist+rad;
                                    surf(app.UIAxes, X_body,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);

                                    x_res_nose = 0:rad/50:rad;
                                    nose_radius_func_ish = sqrt((rad^2) -(rad^2).*((x_res_nose-rad).^2)./(rad^2));
                                    [Z, Y, X] = cylinder(nose_radius_func_ish, 100);
                                    X_nose = X*rad+dist;
                                    surf(app.UIAxes, X_nose,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);

                                    x_res_nose = rad:rad/50:2*rad;
                                    nose_radius_func_ish = sqrt((rad^2) -(rad^2).*((x_res_nose-rad).^2)./(rad^2));
                                    [Z, Y, X] = cylinder(nose_radius_func_ish, 100);
                                    X_nose = X*rad+dist+leng-rad;
                                    surf(app.UIAxes, X_nose,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                                else
                                    l_a = sqrt(2)*rad;
                                    [Z, Y, X] = cylinder(rad,100); %make unit cyliner along x axis
                                    X_body = X*(leng-2*l_a)+dist+l_a;
                                    surf(app.UIAxes, X_body,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);

                                    x_res_nose = 0:l_a/50:l_a;
                                    nose_radius_func_ish = sqrt((rad^2) -(rad^2).*((x_res_nose-l_a).^2)./(l_a^2));
                                    [Z, Y, X] = cylinder(nose_radius_func_ish, 100);
                                    X_nose = X*l_a+dist;
                                    surf(app.UIAxes, X_nose,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);

                                    x_res_nose = l_a:l_a/50:2*l_a;
                                    nose_radius_func_ish = sqrt((rad^2) -(rad^2).*((x_res_nose-l_a).^2)./(l_a^2));
                                    [Z, Y, X] = cylinder(nose_radius_func_ish, 100);
                                    X_nose = X*l_a+dist+leng-l_a;
                                    surf(app.UIAxes, X_nose,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                                end

                            else % 2d plot
                                isSpehical = false;
                                if isSpehical
                                    xTank = [dist+l_a,  dist+leng-l_a];
                                    yTank = [rad, rad];
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    xTank = [dist+l_a,  dist+leng-l_a];
                                    yTank = [-rad, -rad];
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)

                                    x_res_nose = 0:rad/50:rad;
                                    yTank = sqrt((rad^2) -(rad^2).*((x_res_nose-rad).^2)./(rad^2));
                                    xTank = x_res_nose+dist;
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    plot(app.UIAxes, xTank, -1*yTank, 'LineStyle','-', 'Color', color)

                                    x_res_nose = rad:rad/50:2*rad;
                                    yTank = sqrt((rad^2) -(rad^2).*((x_res_nose-rad).^2)./(rad^2));
                                    xTank = x_res_nose+dist+leng-2*rad;
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    plot(app.UIAxes, xTank, -1*yTank, 'LineStyle','-', 'Color', color)
                                else
                                    l_a = sqrt(2)*rad;
                                    xTank = [dist+l_a,  dist+leng-l_a];
                                    yTank = [rad, rad];
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    xTank = [dist+l_a,  dist+leng-l_a];
                                    yTank = [-rad, -rad];
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)

                                    x_res_nose = 0:l_a/50:l_a;
                                    yTank = sqrt((rad^2) -(rad^2).*((x_res_nose-l_a).^2)./(l_a^2));
                                    xTank = x_res_nose+dist;
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    plot(app.UIAxes, xTank, -1*yTank, 'LineStyle','-', 'Color', color)

                                    x_res_nose = l_a:l_a/50:2*l_a;
                                    yTank = sqrt((rad^2) -(rad^2).*((x_res_nose-l_a).^2)./(l_a^2));
                                    xTank = x_res_nose+dist+leng-2*l_a;
                                    plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)
                                    plot(app.UIAxes, xTank, -1*yTank, 'LineStyle','-', 'Color', color)
                                end
                            end

                        elseif isa(values{idx}, 'PointMass')

                            ptObj = values{idx};

                            xPos = ptObj.Position(1);
                            yPos = ptObj.Position(2);
                            zPos = ptObj.Position(3);
                            color = ptObj.Color;


                            if app.ThreeDPlot
                                plot3(app.UIAxes, xPos, yPos, zPos, '.', 'MarkerSize', 30, 'Color', color)
                            else
                                plot(app.UIAxes, xPos, yPos, '.', MarkerSize = 30, Color = color);
                            end

                        elseif isa(values{idx}, 'PropulsionSystem')

                            propSys = values{idx};

                            dist = propSys.Position(1);
                            dia = 2*sqrt(propSys.ExitArea / pi);
                            color = propSys.Color;

                            % create the contour:

                            xProp = linspace(0,3*dia);

                            yProp = .2/.25 * (dia-0.3*dia*sin(1.5*pi*(xProp-.06)/(3*dia)));

                            xProp = xProp + dist;

                            endcapX = [xProp(end), xProp(end)];
                            endcapY = [yProp(end), -yProp(end)];

                            if app.ThreeDPlot
                                [Z, Y, X] = cylinder(yProp, 100);

                                X = X*3*dia + dist;

                                plot3(app.UIAxes, X, Y, Z, 'LineStyle','-', 'Color', color)

                            else
                                plot(app.UIAxes, xProp,yProp, "Color", app.lineColor)
                                plot(app.UIAxes, xProp, -yProp, "Color", app.lineColor)
                                plot(app.UIAxes, endcapX, endcapY, app.lineColor)
                            end





                        end
                    end
                end

            else
                return
            end

            leng = numel(compList);

            % for idx = 1:leng
            %
            %     itemData = app.ComponentDetails(idx,:);
            %
            %     componentType = itemData(1);
            %
            %     switch componentType
            %         case 5 % point mass
            %
            %             % for the point mass, the data is in the order:
            %             % Component Type (5)
            %             % Dry Mass
            %             % Dist from Nose
            %             % y dist
            %             % z dist
            %
            %             dist = itemData(2);
            %             y = itemData(3);
            %             z = itemData(4);
            %
            %             x = dist;
            %             y = y;
            %
            %             plot(app.UIAxes, x, y, 'MarkerSize', 10, 'Marker','.','Color','k')
            %
            %
            %
            %     end
            % end
        end

        function plotCoM()
            % plot the center of mass of the rocket
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % get the color scheme of the app, and determine the plotting
            % color from that

            t = theme;
            close(gcf)

            if strcmp(t.BaseColorStyle, 'dark')
                app.lineColor = 'w';
            else
                app.lineColor = 'k';
            end

            app.Geoplotter();

        end

        % Value changed function: RocketLengthEditField
        function RocketLengthChanged(app, event)

            rocketLeng = app.RocketLengthEditField.Value;

            if rocketLeng < app.NoseConeLengthmEditField.Value
                uialert(app.UIFigure, "The rocket cannot be " + ...
                    "shorter than the nose cone length!", "Input Error", "Icon","error");
                return

            end

            if app.autoRefresh
                app.RocketPlotter();
            end


        end

        % Value changed function: LatitudedegEditField
        function LatitudeChanged(app, event)
            value = app.LatitudedegEditField.Value;
            app.Geoplotter();
        end

        % Value changed function: RocketDiameterEditField
        function RocketDiaChanged(app, event)

            if app.autoRefresh
                app.RocketPlotter();
            end

        end

        % Value changed function: NoseConeLengthmEditField
        function NoseCoseLengthChanged(app, event)

            noseLeng = app.NoseConeLengthmEditField.Value;

            if noseLeng >= app.RocketLengthEditField.Value
                uialert(app.UIFigure, "Nose cone cannot be longer than rocket body!", ...
                    "Dimension Error")
                return
            end

            if app.autoRefresh
                app.RocketPlotter()
            end

        end

        % Value changed function: NoseConeGeometryDropDown
        function NoseConeTypeChanged(app, event)

            if app.autoRefresh
                app.RocketPlotter();
            end
        end

        % Button pushed function: Switchto3D
        function ConvertToThreeD(app, event)
            app.Switchto2D.Enable = "on";
            app.Switchto2D.Visible = 'on';
            app.Switchto3D.Enable = 'off';
            app.Switchto3D.Visible = "off";

            app.ThreeDPlot = 1;

            app.RocketPlotter();
        end

        % Button pushed function: Switchto2D
        function SwitchToTwoD(app, event)
            app.Switchto3D.Enable = "on";
            app.Switchto3D.Visible = 'on';
            app.Switchto2D.Enable = 'off';
            app.Switchto2D.Visible = 'off';

            app.ThreeDPlot = 0;
            app.RocketPlotter();
        end

        % Value changed function: LongitudedegEditField
        function longitudeChanged(app, event)
            value = app.LongitudedegEditField.Value;

            app.Geoplotter();
        end

        % Button pushed function: AddComponentButton
        function AddComponent(app, event)
            if isempty(app.rocket)
                uialert(app.UIFigure, 'No Rocket Object Found', 'Please create or load a rocket first!')
                return
            else
                nSuperProperties = length(properties('RocketComponent'))-1;
                nProperties = length(app.PropertyEditFields);
                propertyArray = app.PropertyEditFields;

                componentVal = string(propertyArray(nProperties-(nSuperProperties-1)).Value);
                componentClass = string(app.ComponentSelectionDropDown.Value);
                component = feval(componentClass, componentVal);

                for idx = 1:nProperties
                    property = string(app.PropertyEditLabels(idx).Text);
                    value = app.PropertyEditFields(idx).Value;

                    if property == 'Position'
                        value = str2num(app.PropertyEditFields(idx).Value);
                    end

                    if isempty(value)
                        uialert(app.UIFigure, sprintf("Component not added, missing property: %s", property), "Error");
                        return

                    elseif property == 'Position'
                        component.(property) = str2num(app.PropertyEditFields(idx).Value);
                    else
                        component.(property) = app.PropertyEditFields(idx).Value;
                    end
                end

                % add a warning if that name already exists:
                if numEntries(app.rocket.ComponentList) ~= 0
                    if any(app.rocket.ComponentList.keys == component.Name)
                        uialert(app.UIFigure, "A component with this name already exists", "Edit object name")
                        return
                    end
                end

                % add the component to the rocket object
                app.rocket.addComponent(component);

                % save the file with the new components:
                name = string(app.RocketNameEditField.Value);
                subpath = filesep + "Inputs" + filesep + "Saved Rockets" + filesep + name + ".mat";

                fullpath = fullfile(pwd, subpath);

                rocketObj = app.rocket;
                save(fullpath, "rocketObj")

                % add the component to the tree:
                uitreenode(app.rootNode, 'Text', component.Name);
                expand(app.rootNode)

                uiconfirm(app.UIFigure, 'Component Addition Successful', 'Component Addition')

                app.RocketPlotter();
            end
        end

        % Selection change function: TabGroup2
        function TabGroupChanged(app, event)

        end

        % Button pushed function: SimulateLaunchButton
        function SimulateLaunchClicked(app, event)
            if isempty(app.rocket)
                uialert(app.UIFigure, "The rocket is empty!", "Setup error!")
            end

            if isempty(app.env)
                uialert(app.UIFigure, "The environment is empty!", "Setup error!")
            end


            app.Settings = IntegratorSettings("apogee", 0.1, "medium");

            Main(app.rocket, app.env, app.Settings, app.rocket_real_name.Value);
        end

        % Drop down opening function: ComponentSelectionDropDown
        function ComponentSelectionDropDownOpening(app, event)
            pathy = ['Classes' filesep 'Components'];

            fileStruct = dir(pathy);

            % Convert to string array
            for idx = 3:length(fileStruct)
                files(idx-2) = string(fileStruct(idx).name);
            end

            % Remove 'RocketComponent.m', it is abstract
            files(files == "RocketComponent.m") = [];

            componentNames = erase(files, ".m");

            app.ComponentSelectionDropDown.Items = componentNames;
        end

        % Button pushed function: UpdatePlotButton
        function UpdatePlotButtonPushed(app, event)
            app.RocketPlotter();
        end

        % Button pushed function: AeroDataButton
        function AeroDataButtonPushed(app, event)
            [file, inPath] = uigetfile('*.csv', 'Select a RASAero csv File');
            if file == 0
                uialert(app.UIFigure, 'No File Selected!', 'File Selection Error')
                return
            end
            name = string(app.RocketNameEditField.Value);
            if isempty(name) || name == "" || name == "Rocket Name"
                uialert(app.UIFigure, ...
                    "Please enter the Rocket Name first.", ...
                    "Input Error");
                return
            end
            src = [inPath, file];
            app.AeroLoc = src;
            % Update UI text
            app.AeroDataButton.Text = file;
            % Update rocket if it exists
            if ~isempty(app.rocket)
                app.rocket.AeroData = src;
                app.rocket.AeroData = setAeroData(app.rocket, app.rocket.AeroData);
            end
        end

        % Button pushed function: LoadRocketButton
        function LoadRocketButtonPushed(app, event)
            subpath = [filesep 'Inputs' filesep 'Saved Rockets'];

            pathy = fullfile(pwd,subpath);

            [file, path] = uigetfile('*.mat', 'Select a Stored Rocket File', pathy);

            if file ~= 0

                filepath = fullfile(path, file);
                app.rocket = load(filepath, "rocketObj").rocketObj;

                % set the properties in the rocket, components are auto
                % added

                app.AeroDataButton.Text = app.rocket.Name + "_aero.csv";
                app.RocketNameEditField.Value = app.rocket.Name;

                if ~isempty(app.rocket.TotalLength)
                    app.RocketLengthEditField.Value = app.rocket.TotalLength;
                end

                if ~isempty(app.rocket.OuterDiameter)
                    app.RocketDiameterEditField.Value = app.rocket.OuterDiameter;
                end

                app.NoseConeLengthmEditField.Value = app.rocket.NoseLength;
                app.NoseConeGeometryDropDown.Value = app.rocket.NoseGeometry;
                app.WetMasskgEditField.Value = app.rocket.TotalMass;

                app.autoRefresh = 1;

                % create the root node:
                rootNode = uitreenode(app.Tree, 'Text', app.rocket.Name);
                expand(rootNode);

                app.rootNode = rootNode;

                % create the nodes for each of the components:
                numVals = numEntries(app.rocket.ComponentList);

                if numVals ~= 0
                    names = app.rocket.ComponentList.keys;

                    for idx = 1:numVals
                        uitreenode(app.rootNode, 'Text', names(idx));
                    end
                    expand(app.rootNode)
                end

                app.RocketPlotter();
            end


        end

        % Button pushed function: SaveRocketButton
        function SaveRocketButtonPushed(app, event)
            name = string(app.RocketNameEditField.Value);
            matfilePath = "Inputs" + filesep + "Saved Rockets" + filesep + name + ".mat";
            filepath = "Inputs" + filesep + "RASAero" + filesep + name + ".csv";
            if isempty(app.rootNode)
                rootNode = uitreenode(app.Tree, 'Text', name);
                expand(rootNode);
                app.rootNode = rootNode;
                uialert(app.UIFigure, "Rocket saved.", "Congratulations!", "Icon","success")
            else
                uialert(app.UIFigure, "Rocket Object already exists." + ...
                    " Components have been saved.", "Rocket Saved",Icon="success")
            end

            % only create a new object if it is empty
            if isempty(app.rocket)
                app.rocket = Rocket(name);
            end

            app.rocket.TotalLength = app.RocketLengthEditField.Value;
            app.rocket.OuterDiameter = app.RocketDiameterEditField.Value;

            % only add the aero data the first time, otherwise don't need
            % it

            if isempty(app.rocket.AeroData)
                if ~strcmp(app.AeroDataButton.Text, "Select File")
                    app.rocket.AeroData = setAeroData(app.rocket, app.AeroLoc);
                end

            else
                app.rocket.AeroData = app.rocket.AeroData;
            end
            app.rocket.NoseLength = app.NoseConeLengthmEditField.Value;
            app.rocket.NoseGeometry = app.NoseConeGeometryDropDown.Value;
            app.rocket.TotalMass = app.WetMasskgEditField.Value;
            rocketObj = app.rocket;
            save(fullfile(pwd, matfilePath), "rocketObj")
            % update the tree node with the rocket object:
            % after creating a rocket object, auto refresh the plot with
            % changes:
            app.autoRefresh = 1;
            app.RocketPlotter();
        end

        % Value changed function: RocketNameEditField
        function RocketNameChanged(app, event)
            value = app.RocketNameEditField.Value;

            app.UIAxes.Title.String = [value, ' Layout'];

        end

        % Value changed function: ComponentSelectionDropDown
        function ComponentSelectionDropDownValueChanged(app, event)
            type = string(app.ComponentSelectionDropDown.Value);

            propertyList = string(properties(type));
            propertyArray = matlab.metadata.Class.fromName(type).PropertyList;

            % set the names to the property description, if it exists:

            for idx = 1:numel(propertyList)
                if ~isempty(propertyArray(idx).Description)
                    propertyName(idx) = string(propertyArray(idx).Description);

                else
                    % if there is no associated description, default to the
                    % name
                    propertyName(idx) = propertyArray(idx).Name;
                end
            end

            delete(app.PropertyGrid.Children);

            nFields = length(propertyList);
            app.PropertyGrid.RowHeight = repmat({'fit'}, 1, nFields+1);


            for idx = 1:nFields
                propertyEditLabels(idx) = uilabel(app.PropertyGrid, 'Text', propertyList(idx));
                propertyEditLabels(idx).Layout.Row = idx;
                propertyEditLabels(idx).Layout.Column = 1;
                app.PropertyEditLabels = propertyEditLabels;

                % Check if the property is 'Position' and handle it differently
                if propertyList(idx) == "Position"
                    % Create a custom input field for Position (3x1 vector)
                    propertyEditFields(idx) = uieditfield(app.PropertyGrid, 'text');
                    propertyEditFields(idx).Layout.Row = idx;
                    propertyEditFields(idx).Layout.Column = 2;
                    propertyEditFields(idx).Value = '[0, 0, 0]'; % Default value for Position
                    %propertyEditFields(idx).UserData = 'Position'; % Mark as Position for validation

                elseif propertyList(idx) == "Color"
                    % Create a color picker for the Color property
                    propertyEditFields(idx) = uicolorpicker(app.PropertyGrid);
                    propertyEditFields(idx).Layout.Row = idx;
                    propertyEditFields(idx).Layout.Column = 2;  % Input field on the left side
                    propertyEditFields(idx).Value = [0.7, 0.7, 0.7];  % Default to white color (RGB)
                    propertyEditFields(idx).UserData = 'Color';  % Mark as Color for validation

                elseif contains(string(propertyArray(idx).Validation.Class.Name), ["int", "double", "single"])
                    fieldType = 'numeric';
                    propertyEditFields(idx) = uieditfield(app.PropertyGrid, fieldType);
                    propertyEditFields(idx).Layout.Row = idx;
                    propertyEditFields(idx).Layout.Column = 2;
                else
                    fieldType = 'text';
                    propertyEditFields(idx) = uieditfield(app.PropertyGrid, fieldType);
                    propertyEditFields(idx).Layout.Row = idx;
                    propertyEditFields(idx).Layout.Column = 2;
                end

                app.PropertyEditFields = propertyEditFields;
                if type == "Fins"
                    propertyEditFields(idx).ValueChangedFcn = @(src, event) FinPlotter(app, nFields);
                end
            end

            % Add a plot for "Fins" if it is selected from the component selection
            if type == "Fins"
                app.FinPlotter(nFields)
            end
        end

        % Button pushed function: GetWeatherConditionsButton
        function getWeather(app, event)
            % get the weather if the appropriate fields are filled out:

            if isempty(app.Date)
                uialert(app.UIFigure, "Date field is empty!", "Input Error")
            end

            if isempty(app.LatitudedegEditField.Value) || isempty(app.LongitudedegEditField.Value)
                uialert(app.UIFigure, "Location Field is empty!", "Input Error")
            end

            lat = app.LatitudedegEditField.Value;
            lon = app.LongitudedegEditField.Value;


            switch app.launch_site_name.Value
                case "FAR"
                    rail_height = 18.29; % Rail at FAR [m]
                case "Launch Trailer"
                    rail_height = 6.096; % Rail on the Launch Trailer [m]
                otherwise
                    error("Invalid Launch Site Name!")
            end

            app.env = Environment(lat, lon, app.Date, rail_height);

            % get the weather for that environment
            app.env = getLocalWeather(app.env);

            uialert(app.UIFigure, "Weather Succesfully Downloaded!", "Success", "Icon","success");


        end

        % Callback function: RocketColor
        function BaseColorChanged(app, event)
            value = app.RocketColor.Value;

            if app.autoRefresh
                app.RocketPlotter();
            end
        end

        % Clicked callback: Tree
        function NodeClicked(app, event)
            node = event.InteractionInformation.Node;
            
            
            if ~isempty(node)
            % if ~isempty(node) || node
                print(node)

                selected = node.Text;

                % parse which node was clicked. This will have the same name as
                % the dictionary item:

                compList = app.rocket.ComponentList;

                % return the cell array with the correct values:
                fprintf("selected: %s\n", selected)
                componentCell = compList(selected);

                comp = componentCell{1};
            end

        end

        % Value changed function: DateSelectionDatePicker
        function dateSelected(app, event)
            value = app.DateSelectionDatePicker.Value;

            % illusion of free will, set the date to now:

            app.Date = datetime("now", "TimeZone","UTC");

        end


        % Double-clicked callback: Tree
        function NodeDoubleClick(app, event)
            node = event.InteractionInformation.Node;
            if ~isempty(node)
                selected = node.Text;

                compList = app.rocket.ComponentList;

                % return the cell array with the correct values:
                componentCell = compList(selected);

                % confirm that you want to delete this component:
                selection = uiconfirm(app.UIFigure, "Are you sure you want to delete this component?" ...
                    , "Deletion Confirmation");

                if selection == 'OK'
                    % remove that component from the rocket
                    app.rocket.removeComponent(selected);

                    % remove it from the tree:
                    node.delete;
                end

                % save the file with the components removed:
                name = string(app.RocketNameEditField.Value);
                path = filesep + "Inputs" + filesep + "Saved Rockets" + filesep + name + ".mat";

                fullpath = fullfile(pwd, path);

                rocketObj = app.rocket;
                save(fullpath, "rocketObj")

                app.RocketPlotter();
            end


        end

        % Value changed function: WetMasskgEditField
        function TotalMassChanged(app, event)
            value = app.WetMasskgEditField.Value;

            app.rocket.TotalMass = value;

            if app.autoRefresh
                app.RocketPlotter();
            end

        end

        % Button pushed function: SetBasePathButton
        function setBasePath(app, event)
            % set the base
            selpath = uigetdir(pwd,"Set the root directory to your local 'PSPL-6DOF\TheSixDoF' folder");

            % make sure that this path includes the 'TheSixDoF folder'

            if ~contains(selpath, 'TheSixDoF')

                uialert(app.UIFigure, "The path does not contain the correct folder!", "Path selection error")

                return
            end


            % with the correct selection, change to the correct folder and
            % change and add the correct subfolders:
            cd(selpath)

            addpath(genpath('Classes'))

            % enable the main user input now that this has been done:
            app.Panel.Enable = "on";
            app.UpdatePlotButton.Enable = "on";
            app.Switchto3D.Enable = "on";
        end

        % Button pushed function: GUIHelpButton
        function OpenHelp(app, event)

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % disp("SEXSEXSEXSEXSEX")

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 481];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 2 640 480];

            % Create RocketDesignTab
            app.RocketDesignTab = uitab(app.TabGroup);
            app.RocketDesignTab.Title = 'Rocket Design';

            % Create GridLayout9
            % I think this is the rocket design grid? - david
            app.GridLayout9 = uigridlayout(app.RocketDesignTab);
            app.GridLayout9.ColumnWidth = {175, 275, '2x', '2x'};
            app.GridLayout9.RowHeight = {25, '1.3x', 22, '1x'};
            app.GridLayout9.ColumnSpacing = 5.04001007080078;
            app.GridLayout9.RowSpacing = 5.07499361038208;
            app.GridLayout9.Padding = [5.04001007080078 5.07499361038208 5.04001007080078 5.07499361038208];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout9);
            title(app.UIAxes, 'Rocket Layout', 'Interpreter', 'latex')
            xlabel(app.UIAxes, 'X', 'Interpreter', 'latex')
            ylabel(app.UIAxes, 'Y', 'Interpreter', 'latex')
            zlabel(app.UIAxes, 'Z', 'Interpreter', 'latex')
            app.UIAxes.Layout.Row = [1 2];
            app.UIAxes.Layout.Column = [3 4];

            % Create Panel
            app.Panel = uipanel(app.GridLayout9);
            app.Panel.Enable = 'off';
            app.Panel.Layout.Row = [2 4];
            app.Panel.Layout.Column = [1 2];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Panel);
            app.GridLayout.ColumnWidth = {165, 156};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.GridLayout);
            app.TabGroup2.SelectionChangedFcn = createCallbackFcn(app, @TabGroupChanged, true);
            app.TabGroup2.Layout.Row = [1 7];
            app.TabGroup2.Layout.Column = [1 2];

            % Create RocketTab
            app.RocketTab = uitab(app.TabGroup2);
            app.RocketTab.Title = 'Rocket';

            % Create RocketGrid
            app.RocketGrid = uigridlayout(app.RocketTab);
            app.RocketGrid.ColumnWidth = {'1x', '2x'};
            app.RocketGrid.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % the name of the variable should be the same as the name in
            % the GUI but im not fixing allat - david
            % Create TotalLengthmLabel
            app.TotalLengthmLabel = uilabel(app.RocketGrid);
            app.TotalLengthmLabel.HorizontalAlignment = 'center';
            app.TotalLengthmLabel.WordWrap = 'on';
            app.TotalLengthmLabel.Layout.Row = 2;
            app.TotalLengthmLabel.Layout.Column = 1;
            app.TotalLengthmLabel.Text = 'Airframe Length [m]';

            % Create RocketLengthEditField
            app.RocketLengthEditField = uieditfield(app.RocketGrid, 'numeric');
            app.RocketLengthEditField.Limits = [0 Inf];
            app.RocketLengthEditField.AllowEmpty = 'on';
            app.RocketLengthEditField.ValueChangedFcn = createCallbackFcn(app, @RocketLengthChanged, true);
            app.RocketLengthEditField.HorizontalAlignment = 'center';
            app.RocketLengthEditField.Layout.Row = 2;
            app.RocketLengthEditField.Layout.Column = 2;
            app.RocketLengthEditField.Value = [];

            % Create RocketOuterDiametermLabel
            app.RocketOuterDiametermLabel = uilabel(app.RocketGrid);
            app.RocketOuterDiametermLabel.HorizontalAlignment = 'center';
            app.RocketOuterDiametermLabel.WordWrap = 'on';
            app.RocketOuterDiametermLabel.Layout.Row = 3;
            app.RocketOuterDiametermLabel.Layout.Column = 1;
            app.RocketOuterDiametermLabel.Text = 'Airframe Outer Diameter [m]';

            % Create RocketDiameterEditField
            app.RocketDiameterEditField = uieditfield(app.RocketGrid, 'numeric');
            app.RocketDiameterEditField.Limits = [0 Inf];
            app.RocketDiameterEditField.AllowEmpty = 'on';
            app.RocketDiameterEditField.ValueChangedFcn = createCallbackFcn(app, @RocketDiaChanged, true);
            app.RocketDiameterEditField.HorizontalAlignment = 'center';
            app.RocketDiameterEditField.Layout.Row = 3;
            app.RocketDiameterEditField.Layout.Column = 2;
            app.RocketDiameterEditField.Value = [];

            % Create NoseConeGeometryDropDownLabel
            app.NoseConeGeometryDropDownLabel = uilabel(app.RocketGrid);
            app.NoseConeGeometryDropDownLabel.HorizontalAlignment = 'center';
            app.NoseConeGeometryDropDownLabel.WordWrap = 'on';
            app.NoseConeGeometryDropDownLabel.Layout.Row = 5;
            app.NoseConeGeometryDropDownLabel.Layout.Column = 1;
            app.NoseConeGeometryDropDownLabel.Text = 'Nose Cone Geometry';

            % Create NoseConeGeometryDropDown
            app.NoseConeGeometryDropDown = uidropdown(app.RocketGrid);
            app.NoseConeGeometryDropDown.Items = {'Von Karman', 'Tangent Ogive', 'Conic', 'Elliptical'};
            app.NoseConeGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @NoseConeTypeChanged, true);
            app.NoseConeGeometryDropDown.Layout.Row = 5;
            app.NoseConeGeometryDropDown.Layout.Column = 2;
            app.NoseConeGeometryDropDown.Value = 'Von Karman';

            % Create NoseConeLengthmEditFieldLabel
            app.NoseConeLengthmEditFieldLabel = uilabel(app.RocketGrid);
            app.NoseConeLengthmEditFieldLabel.HorizontalAlignment = 'center';
            app.NoseConeLengthmEditFieldLabel.WordWrap = 'on';
            app.NoseConeLengthmEditFieldLabel.Layout.Row = 6;
            app.NoseConeLengthmEditFieldLabel.Layout.Column = 1;
            app.NoseConeLengthmEditFieldLabel.Text = 'Nose Cone Length [m]';

            % Create NoseConeLengthmEditField
            app.NoseConeLengthmEditField = uieditfield(app.RocketGrid, 'numeric');
            app.NoseConeLengthmEditField.Limits = [0 Inf];
            app.NoseConeLengthmEditField.AllowEmpty = 'on';
            app.NoseConeLengthmEditField.ValueChangedFcn = createCallbackFcn(app, @NoseCoseLengthChanged, true);
            app.NoseConeLengthmEditField.HorizontalAlignment = 'center';
            app.NoseConeLengthmEditField.Layout.Row = 6;
            app.NoseConeLengthmEditField.Layout.Column = 2;
            app.NoseConeLengthmEditField.Value = [];

            % Create NameEditField_2Label
            app.NameEditField_2Label = uilabel(app.RocketGrid);
            app.NameEditField_2Label.HorizontalAlignment = 'center';
            app.NameEditField_2Label.Layout.Row = 1;
            app.NameEditField_2Label.Layout.Column = 1;
            app.NameEditField_2Label.Text = 'Name';

            % Create RocketNameEditField
            app.RocketNameEditField = uieditfield(app.RocketGrid, 'text');
            app.RocketNameEditField.ValueChangedFcn = createCallbackFcn(app, @RocketNameChanged, true);
            app.RocketNameEditField.HorizontalAlignment = 'center';
            app.RocketNameEditField.Layout.Row = 1;
            app.RocketNameEditField.Layout.Column = 2;
            app.RocketNameEditField.Value = 'Rocket Name';


            % fuck this - david

            rocket_real_name = "Rocket A";
            % rocket_real_name  = "CMS";
            % rocket_real_name  = "Copperhead";

            switch rocket_real_name
                case "Rocket A"
                    launch_site_name = "Launch Trailer";
                case {"CMS", "Copperhead"}
                    launch_site_name = "FAR";
                otherwise
                    error("Invalid Rocket Name!")
            end

            fprintf("--------- Rocket Name: %s ---------\n", rocket_real_name);
            fprintf("--------- Launch Site Name: %s ---------\n", launch_site_name);


            % what the fuck am i doing - david
            app.rocket_real_name = uieditfield(app.RocketGrid);
            app.rocket_real_name.Value = rocket_real_name;

            app.launch_site_name = uieditfield(app.RocketGrid);
            app.launch_site_name.Value = launch_site_name;



            % Create AeroDataButton
            app.AeroDataButton = uibutton(app.RocketGrid, 'push');
            app.AeroDataButton.ButtonPushedFcn = createCallbackFcn(app, @AeroDataButtonPushed, true);
            app.AeroDataButton.Layout.Row = 4;
            app.AeroDataButton.Layout.Column = 2;
            app.AeroDataButton.Text = 'Select File';

            % Create RASAeroDataLabel
            app.RASAeroDataLabel = uilabel(app.RocketGrid);
            app.RASAeroDataLabel.HorizontalAlignment = 'center';
            app.RASAeroDataLabel.WordWrap = 'on';
            app.RASAeroDataLabel.Layout.Row = 4;
            app.RASAeroDataLabel.Layout.Column = 1;
            app.RASAeroDataLabel.Text = 'RASAero Data';

            % Create ButtonPanel
            app.ButtonPanel = uipanel(app.RocketGrid);
            app.ButtonPanel.BorderType = 'none';
            app.ButtonPanel.Layout.Row = [9 10];
            app.ButtonPanel.Layout.Column = [1 2];

            % Create ButtonGrid
            app.ButtonGrid = uigridlayout(app.ButtonPanel);

            % Create LoadRocketButton
            app.LoadRocketButton = uibutton(app.ButtonGrid, 'push');
            app.LoadRocketButton.ButtonPushedFcn = createCallbackFcn(app, @LoadRocketButtonPushed, true);
            app.LoadRocketButton.Layout.Row = 1;
            app.LoadRocketButton.Layout.Column = [1 2];
            app.LoadRocketButton.Text = 'Load Rocket';

            % Create SaveRocketButton
            app.SaveRocketButton = uibutton(app.ButtonGrid, 'push');
            app.SaveRocketButton.ButtonPushedFcn = createCallbackFcn(app, @SaveRocketButtonPushed, true);
            app.SaveRocketButton.Layout.Row = 2;
            app.SaveRocketButton.Layout.Column = [1 2];
            app.SaveRocketButton.Text = 'Save Rocket';

            % Create ColorColorPickerLabel
            app.ColorColorPickerLabel = uilabel(app.RocketGrid);
            app.ColorColorPickerLabel.HorizontalAlignment = 'center';
            app.ColorColorPickerLabel.Layout.Row = 7;
            app.ColorColorPickerLabel.Layout.Column = 1;
            app.ColorColorPickerLabel.Text = 'Color';

            % Create RocketColor
            app.RocketColor = uicolorpicker(app.RocketGrid);
            app.RocketColor.Value = [0.8 0.8 0.8];
            app.RocketColor.ValueChangedFcn = createCallbackFcn(app, @BaseColorChanged, true);
            app.RocketColor.Layout.Row = 7;
            app.RocketColor.Layout.Column = 2;

            % Create WetMasskgEditFieldLabel
            app.WetMasskgEditFieldLabel = uilabel(app.RocketGrid);
            app.WetMasskgEditFieldLabel.HorizontalAlignment = 'center';
            app.WetMasskgEditFieldLabel.Layout.Row = 8;
            app.WetMasskgEditFieldLabel.Layout.Column = 1;
            app.WetMasskgEditFieldLabel.Text = 'Wet Mass [kg]';

            % Create WetMasskgEditField
            app.WetMasskgEditField = uieditfield(app.RocketGrid, 'numeric');
            app.WetMasskgEditField.ValueChangedFcn = createCallbackFcn(app, @TotalMassChanged, true);
            app.WetMasskgEditField.Layout.Row = 8;
            app.WetMasskgEditField.Layout.Column = 2;

            % Create ComponentsTab
            app.ComponentsTab = uitab(app.TabGroup2);
            app.ComponentsTab.Title = 'Components';

            % Create ComponentsGrid
            app.ComponentsGrid = uigridlayout(app.ComponentsTab);
            app.ComponentsGrid.ColumnWidth = {'1x', '2x'};
            app.ComponentsGrid.RowHeight = {50, '1x', '1x', '1x', '1x', '1x', 55};

            % Create ComponentSelectionDropDownLabel
            app.ComponentSelectionDropDownLabel = uilabel(app.ComponentsGrid);
            app.ComponentSelectionDropDownLabel.HorizontalAlignment = 'center';
            app.ComponentSelectionDropDownLabel.WordWrap = 'on';
            app.ComponentSelectionDropDownLabel.Layout.Row = 1;
            app.ComponentSelectionDropDownLabel.Layout.Column = 1;
            app.ComponentSelectionDropDownLabel.Text = 'Component Selection';

            % Create ComponentSelectionDropDown
            app.ComponentSelectionDropDown = uidropdown(app.ComponentsGrid);
            app.ComponentSelectionDropDown.Items = {};
            app.ComponentSelectionDropDown.DropDownOpeningFcn = createCallbackFcn(app, @ComponentSelectionDropDownOpening, true);
            app.ComponentSelectionDropDown.ValueChangedFcn = createCallbackFcn(app, @ComponentSelectionDropDownValueChanged, true);
            app.ComponentSelectionDropDown.Layout.Row = 1;
            app.ComponentSelectionDropDown.Layout.Column = 2;
            app.ComponentSelectionDropDown.Value = {};

            % Create AddComponentButton
            app.AddComponentButton = uibutton(app.ComponentsGrid, 'push');
            app.AddComponentButton.ButtonPushedFcn = createCallbackFcn(app, @AddComponent, true);
            app.AddComponentButton.Layout.Row = 7;
            app.AddComponentButton.Layout.Column = [1 2];
            app.AddComponentButton.Text = 'Add Component';

            % Create PropertyPanel
            app.PropertyPanel = uipanel(app.ComponentsGrid);
            app.PropertyPanel.BorderType = 'none';
            app.PropertyPanel.Layout.Row = [2 6];
            app.PropertyPanel.Layout.Column = [1 2];

            % Create PropertyGrid
            app.PropertyGrid = uigridlayout(app.PropertyPanel);
            app.PropertyGrid.ColumnWidth = {'1x', '2x'};
            app.PropertyGrid.RowHeight = {'1x', '1x', '1x', '1x', '1x'};
            app.PropertyGrid.Scrollable = 'on';

            % Create Switchto2D
            app.Switchto2D = uibutton(app.GridLayout9, 'push');
            app.Switchto2D.ButtonPushedFcn = createCallbackFcn(app, @SwitchToTwoD, true);
            app.Switchto2D.Enable = 'off';
            app.Switchto2D.Visible = 'off';
            app.Switchto2D.Layout.Row = 3;
            app.Switchto2D.Layout.Column = 4;
            app.Switchto2D.Text = '2D View';

            % Create UpdatePlotButton
            app.UpdatePlotButton = uibutton(app.GridLayout9, 'push');
            app.UpdatePlotButton.ButtonPushedFcn = createCallbackFcn(app, @UpdatePlotButtonPushed, true);
            app.UpdatePlotButton.Enable = 'off';
            app.UpdatePlotButton.Layout.Row = 3;
            app.UpdatePlotButton.Layout.Column = 3;
            app.UpdatePlotButton.Text = 'Update Plot';

            % Create Switchto3D
            app.Switchto3D = uibutton(app.GridLayout9, 'push');
            app.Switchto3D.ButtonPushedFcn = createCallbackFcn(app, @ConvertToThreeD, true);
            app.Switchto3D.Enable = 'off';
            app.Switchto3D.Layout.Row = 3;
            app.Switchto3D.Layout.Column = 4;
            app.Switchto3D.Text = '3D View';

            % Create ListBox
            app.ListBox = uilistbox(app.GridLayout9);
            app.ListBox.Items = {};
            app.ListBox.Enable = 'off';
            app.ListBox.Visible = 'off';
            app.ListBox.Layout.Row = 4;
            app.ListBox.Layout.Column = 3;
            app.ListBox.Value = {};

            % Create Tree
            app.Tree = uitree(app.GridLayout9);
            app.Tree.Layout.Row = 4;
            app.Tree.Layout.Column = 3;
            app.Tree.ClickedFcn = createCallbackFcn(app, @NodeClicked, true);
            app.Tree.DoubleClickedFcn = createCallbackFcn(app, @NodeDoubleClick, true);

            % Create PropGrid2
            app.PropGrid2 = uigridlayout(app.GridLayout9);
            app.PropGrid2.ColumnWidth = {'1x', '2x'};
            app.PropGrid2.RowHeight = {'1x', '1x', '1x'};
            app.PropGrid2.Layout.Row = 4;
            app.PropGrid2.Layout.Column = 4;

            % Create SetBasePathButton
            app.SetBasePathButton = uibutton(app.GridLayout9, 'push');
            app.SetBasePathButton.ButtonPushedFcn = createCallbackFcn(app, @setBasePath, true);
            app.SetBasePathButton.Layout.Row = 1;
            app.SetBasePathButton.Layout.Column = 1;
            app.SetBasePathButton.Text = 'Set Base Path';

            % Create GUIHelpButton
            app.GUIHelpButton = uibutton(app.GridLayout9, 'push');
            app.GUIHelpButton.ButtonPushedFcn = createCallbackFcn(app, @OpenHelp, true);
            app.GUIHelpButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.GUIHelpButton.FontWeight = 'bold';
            app.GUIHelpButton.FontColor = [1 1 1];
            app.GUIHelpButton.Layout.Row = 1;
            app.GUIHelpButton.Layout.Column = 2;
            app.GUIHelpButton.Text = 'GUI Help';

            % Create SimulationTab
            app.SimulationTab = uitab(app.TabGroup);
            app.SimulationTab.Title = 'Simulation';

            % Create Panel_2
            app.Panel_2 = uipanel(app.SimulationTab);
            app.Panel_2.Position = [1 0 260 457];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel_2);
            app.GridLayout2.ColumnWidth = {114, 200, "1x"};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create DateSelectionDatePickerLabel
            app.DateSelectionDatePickerLabel = uilabel(app.GridLayout2);
            app.DateSelectionDatePickerLabel.HorizontalAlignment = 'right';
            app.DateSelectionDatePickerLabel.Layout.Row = 1;
            app.DateSelectionDatePickerLabel.Layout.Column = 1;
            app.DateSelectionDatePickerLabel.Text = 'Date Selection';

            % Create DateSelectionDatePicker
            app.DateSelectionDatePicker = uidatepicker(app.GridLayout2);
            app.DateSelectionDatePicker.Limits = [datetime([1940 1 1]) datetime([9999 12 31])];
            app.DateSelectionDatePicker.ValueChangedFcn = createCallbackFcn(app, @dateSelected, true);
            app.DateSelectionDatePicker.Layout.Row = 1;
            app.DateSelectionDatePicker.Layout.Column = 2;
            app.LongitudedegEditField.Value 

  
            FAR_latitude = 35.3474; % [degrees]
            FAR_longitude = -117.8091; % [degrees]


            dairy_farm_latitude = 40.509936707682144; % [degrees]
            dairy_farm_longitude = -87.02196060569756; % [degrees]
            
            switch app.rocket_real_name.Value
                case "Rocket A"
                    default_latitude_value = dairy_farm_latitude;
                    default_longitude_value = dairy_farm_longitude;
                case {"CMS", "Copperhead"}
                    default_latitude_value = FAR_latitude;
                    default_longitude_value = FAR_longitude;
                otherwise
                    error("Invalid Rocket Name!")
            end


            % Create LatitudedegEditFieldLabel
            app.LatitudedegEditFieldLabel = uilabel(app.GridLayout2);
            app.LatitudedegEditFieldLabel.HorizontalAlignment = 'right';
            app.LatitudedegEditFieldLabel.Layout.Row = 2;
            app.LatitudedegEditFieldLabel.Layout.Column = 1;
            app.LatitudedegEditFieldLabel.Text = 'Latitude (deg)';
            
            % Create LatitudedegEditField
            app.LatitudedegEditField = uieditfield(app.GridLayout2, 'numeric');
            app.LatitudedegEditField.Limits = [-90 90];
            app.LatitudedegEditField.ValueChangedFcn = createCallbackFcn(app, @LatitudeChanged, true);
            app.LatitudedegEditField.Layout.Row = 2;
            app.LatitudedegEditField.Layout.Column = 2;
            app.LatitudedegEditField.Value = default_latitude_value;

            % Create LongitudedegEditFieldLabel
            app.LongitudedegEditFieldLabel = uilabel(app.GridLayout2);
            app.LongitudedegEditFieldLabel.HorizontalAlignment = 'right';
            app.LongitudedegEditFieldLabel.Layout.Row = 3;
            app.LongitudedegEditFieldLabel.Layout.Column = 1;
            app.LongitudedegEditFieldLabel.Text = 'Longitude (deg)';


            % Create LongitudedegEditField
            app.LongitudedegEditField = uieditfield(app.GridLayout2, 'numeric');
            app.LongitudedegEditField.Limits = [-180 180];
            app.LongitudedegEditField.ValueChangedFcn = createCallbackFcn(app, @longitudeChanged, true);
            app.LongitudedegEditField.Layout.Row = 3;
            app.LongitudedegEditField.Layout.Column = 2;
            app.LongitudedegEditField.Value = default_longitude_value;

            % Create SimulateLaunchButton
            app.SimulateLaunchButton = uibutton(app.GridLayout2, 'push');
            app.SimulateLaunchButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateLaunchClicked, true);
            app.SimulateLaunchButton.Layout.Row = 7;
            app.SimulateLaunchButton.Layout.Column = [1 2];
            app.SimulateLaunchButton.Text = 'Simulate Launch';

            % Create GetWeatherConditionsButton
            app.GetWeatherConditionsButton = uibutton(app.GridLayout2, 'push');
            app.GetWeatherConditionsButton.ButtonPushedFcn = createCallbackFcn(app, @getWeather, true);
            app.GetWeatherConditionsButton.Layout.Row = 4;
            app.GetWeatherConditionsButton.Layout.Column = [1 2];
            app.GetWeatherConditionsButton.Text = 'Get Weather Conditions';

            % Create LaunchLocationPanel
            app.LaunchLocationPanel = uipanel(app.SimulationTab);
            app.LaunchLocationPanel.TitlePosition = 'centertop';
            app.LaunchLocationPanel.Title = 'Launch Location';
            app.LaunchLocationPanel.Position = [275 237 352 210];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = RocketGUI_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end