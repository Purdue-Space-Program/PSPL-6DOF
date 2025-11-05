classdef RocketGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RocketDesignTab                matlab.ui.container.Tab
        GridLayout9                    matlab.ui.container.GridLayout
        Tree                           matlab.ui.container.Tree
        UITable                        matlab.ui.control.Table
        ListBox                        matlab.ui.control.ListBox
        Switchto3D                     matlab.ui.control.Button
        UpdatePlotButton               matlab.ui.control.Button
        Switchto2D                     matlab.ui.control.Button
        Panel                          matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        TabGroup2                      matlab.ui.container.TabGroup
        RocketTab                      matlab.ui.container.Tab
        RocketGrid                     matlab.ui.container.GridLayout
        RocketColor                    matlab.ui.control.ColorPicker
        ColorColorPickerLabel          matlab.ui.control.Label
        ButtonPanel                    matlab.ui.container.Panel
        ButtonGrid                     matlab.ui.container.GridLayout
        RevertButton                   matlab.ui.control.Button
        CreateNewRocketButton          matlab.ui.control.Button
        SaveRocketButton               matlab.ui.control.Button
        LoadRocketButton               matlab.ui.control.Button
        RASAeroDataLabel               matlab.ui.control.Label
        AeroDataButton                 matlab.ui.control.Button
        RocketNameEditField            matlab.ui.control.EditField
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
        ThreeDPlot % flag to turn 3d plotting on and off
        ComponentList = "" % a list of the components on the rocket
        ComponentDetails
        rocket Rocket
        PropertyEditFields
        PropertyEditLabels
        autoRefresh = 0;
    end

    methods (Access = private)

        function RocketPlotter(app)

            % create general parameters for the rocket:
            leng = app.RocketLengthEditField.Value;
            noseLeng = app.NoseConeLengthmEditField.Value;
            dia = app.RocketDiameterEditField.Value;
            R = dia/2;
            %fins

            % look at the rocket object and see if there are fins
            % associated with it:

            % first, check if the rocket is empty. If not, proceed:
            if ~isempty(app.rocket)
                compList = app.rocket.componentList;

                % in the event that the rocket has components:
                if numEntries(compList) > 0
                    len = numEntries(compList);
                    values = compList.values;

                    % go through each and check for fins
                    for idx = 1:len
                        if isa(values{idx}, 'Fins')
                            finObject = values{idx};
                            num_Fins = double(finObject.Count);
                            rootChord = finObject.RootChord;
                            tipChord = finObject.TipChord;
                            span = finObject.Span;
                            sweep = finObject.Sweep;
                            fin_offset = finObject.Position(1);
                            finColor = finObject.Color;
                        end
                    end

                else
                    num_Fins = 0;
                end


            else
                num_Fins = 0;
            end


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

                % components:
                app.PlotComponents();

                if num_Fins ~= 0

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

                    for n = 1:length(xOut)
                        X_fin(n,1) = (fin_offset) + rootChord.*xOut(n);
                        X_fin(n,2) = (fin_offset) + tipChord.*xOut(n) + sweep;
                        Y_fin(n,1) = yOut(n);
                        Y_fin(n,2) = yOut(n);
                        Z_fin(n,1) = R;
                        Z_fin(n,2) = R + span;
                    end

                    for n = 1:length(xOut)
                        X_fin_top(n) = (fin_offset) + tipChord.*xOut(n) + sweep;
                        Y_fin_top(n) = yOut(n);
                        Z_fin_top(n) = R+span;
                    end

                    scopy = zeros(num_Fins);
                    stcopy = zeros(length(xOut), num_Fins);

                    for i = 1:num_Fins
                        scopy(i) = surf(app.UIAxes,X_fin,Y_fin,Z_fin, "FaceColor",finColor,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                        stcopy(:,i) = fill3(app.UIAxes,X_fin_top,Y_fin_top,Z_fin_top, [1,1,1], 'FaceColor',finColor, 'EdgeAlpha', 0);
                        direction = [1 0 0];
                        origin = [0 0 0];
                        rotate(scopy(i),direction,rad2deg((i-1)*(2*pi)/num_Fins),origin);
                        rotate(stcopy(:,i), direction,rad2deg((i-1)*(2*pi)/num_Fins),origin)

                    end
                end


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
                xlim(app.UIAxes, [-0.02*leng,leng*1.02])
                axis(app.UIAxes, "equal")

            end

        end

        function Geoplotter(app)

            lat = app.LatitudedegEditField.Value;
            long = app.LongitudedegEditField.Value;

            g = geoaxes(app.LaunchLocationPanel);

            set(g, 'Basemap', 'satellite')

            geoplot(g, lat, long, 'ro')

            size = 0.1;

            geolimits(g, [lat-size, lat+size], [long-size,long+size]);
        end

        function [finPtsX, finPtsY] = FinPlotter(app)

            % look at the rocket object and see if there are fins
            % associated with it:

            % first, check if the rocket is empty. If not, proceed:
            if ~isempty(app.rocket)
                compList = app.rocket.componentList;

                % in the event that the rocket has components:
                if numEntries(compList) > 0
                    len = numEntries(compList);
                    values = compList.values;

                    % go through each and check for fins
                    for idx = 1:len
                        if isa(values{idx}, 'Fins')
                            finObject = values{idx};
                            num_Fins = double(finObject.Count);
                            rootChord = finObject.RootChord;
                            tipChord = finObject.TipChord;
                            span = finObject.Span;
                            sweep = finObject.Sweep;
                            fin_offset = finObject.Position(1);
                        end
                    end
                end

            else
                num_Fins = 0;
                return
            end


            % get the basic parameters of the fin:
            rootChord = app.RootChordEditField.Value;
            tipChord = app.TipChordEditField.Value;
            span = app.SpanEditField.Value;
            sweep = app.SweepEditField.Value;

            % the fin will always be defined by four points, with the first
            % point begin [0,0]. These are defined in the order:
            % [0,0]
            % [sweep,span]
            % [sweep+tipChord, span]
            % [rootChord, 0]

            finPtsX = [0, sweep, sweep+tipChord, rootChord, 0];
            finPtsY = [0, span, span, 0, 0];

            plot(app.FinGraph, finPtsX, finPtsY, app.lineColor);
            axis(app.FinGraph, "equal");
        end

        function PlotFins(app)

            % get the default fin geometry from the fin plotter function

            if ~isempty(app.rocket)
                compList = app.rocket.componentList;

                % in the event that the rocket has components:
                if numEntries(compList) > 0
                    len = numEntries(compList);
                    values = compList.values;

                    % go through each and check for fins
                    for idx = 1:len
                        if isa(values{idx}, 'Fins')
                            finObject = values{idx};
                            numFins = double(finObject.Count);
                            rootChord = finObject.RootChord;
                            tipChord = finObject.TipChord;
                            span = finObject.Span;
                            sweep = finObject.Sweep;
                            finOffset = finObject.Position(1);
                        end
                    end

                else
                    return
                end

            else
                numFins = 0;
                return
            end

            % replace this

            try
                finPtsX = [0, sweep, sweep+tipChord, rootChord, 0];
            catch
                return
            end
            finPtsY = [0, span, span, 0, 0];

            % get parameters from user input:
            rearDist = finOffset;
            dia = app.rocket.outerDiameter;
            leng = app.rocket.totalLength;

            % first, check which fins should be plotted based on occlusion
            % (manually for now lmao, don't know how to write this
            % programmatically)

            switch numFins
                case 1
                    % why?
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
            for idx = plotFin

                % the first fin always points up towards us, so use
                % that as the baseline reference (theta = 0):
                theta = (2*pi)/numFins * (idx-1);

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
        end

        function PlotComponents(app)

            % get the components from the rocket

            if ~isempty(app.rocket)
                compList = app.rocket.componentList;

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

                            xTank = [dist, dist, dist+leng, dist+leng, dist];
                            yTank = [-rad, rad, rad, -rad, -rad];

                            if app.ThreeDPlot

                                [Z, Y, X] = cylinder(rad,100); %make unit cyliner along x axis
                                X_body = X*(leng)+dist;
                
                                surf(app.UIAxes, X_body,Y,Z, "FaceColor",color,'FaceAlpha', 0.7, 'EdgeAlpha',0);
                            else

                            plot(app.UIAxes, xTank, yTank, 'LineStyle','-', 'Color', color)

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

            app.RevertButton.Enable = 'off';

            % set up the filepath:
            path

        end

        % Value changed function: RocketLengthEditField
        function RocketLengthChanged(app, event)
            app.RevertButton.Enable = 'on';

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
            app.RevertButton.Enable = 'on';

            if app.autoRefresh
                app.RocketPlotter();
            end

        end

        % Value changed function: NoseConeLengthmEditField
        function NoseCoseLengthChanged(app, event)
            app.RevertButton.Enable = 'on';

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
            app.RevertButton.Enable = 'on';

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

                componentName = string(propertyArray(nProperties-(nSuperProperties-1)).Value);
                componentClass = string(app.ComponentSelectionDropDown.Value);
                component = feval(componentClass, componentName);

                for idx = 1:nProperties
                    property = string(app.PropertyEditLabels(idx).Text);
                    value = app.PropertyEditFields(idx).Value;

                    if property == 'Position'
                        value = str2num(app.PropertyEditFields(idx).Value);
                    end

                    if isempty(value)
                        return
                    elseif property == 'Position'
                        component.(property) = str2num(app.PropertyEditFields(idx).Value);
                    else
                        component.(property) = app.PropertyEditFields(idx).Value;
                    end
                end

                % add the component to the rocket object
                app.rocket.addComponent(component);

                % add the component to the tree:

                uiconfirm(app.UIFigure, 'Component Addition Successful', 'Component Addition')

                app.RocketPlotter();
            end
        end

        % Selection change function: TabGroup2
        function TabGroupChanged(app, event)

        end

        % Button pushed function: SimulateLaunchButton
        function SimulateLaunchClicked(app, event)
            Main();
        end

        % Drop down opening function: ComponentSelectionDropDown
        function ComponentSelectionDropDownOpening(app, event)
            fileStruct = dir("TheSixDoF\Classes\Components");
            for idx = 3:length(fileStruct)
                files(idx-2) = string(fileStruct(idx).name);
            end
            componentNames = erase(files, ".m");

            app.ComponentSelectionDropDown.Items = componentNames;
        end

        % Button pushed function: UpdatePlotButton
        function UpdatePlotButtonPushed(app, event)
            app.RocketPlotter();
        end

        % Button pushed function: AeroDataButton
        function AeroDataButtonPushed(app, event)
            [file, path] = uigetfile('*.csv', 'Select a RASAero csv File');

            if file == 0
                uialert(app.UIFigure, 'No File Selected!', 'File Selection Error')
                return
            else
                name = app.RocketNameEditField.Value;
                file = string(file);
                path = string(path);
                if isempty(name)
                    error('Please enter the Rocket Name first')
                else
                    name = string(name);
                    filename = name + "_aero.csv";
                    filepath = fullfile(path, file);
                    savepath = "TheSixDoF" + filesep + "Inputs" + filesep + "RASAero" + filesep + name + ".csv";
                    path = pwd;
                    savepath = fullfile(path, savepath);

                    if savepath ~= filepath
                        % maybe change this
                        movefile(filepath, savepath);
                    end

                    app.AeroDataButton.Text = filename;
                end

            end
        end

        % Button pushed function: LoadRocketButton
        function LoadRocketButtonPushed(app, event)
            [file, path] = uigetfile('*.mat', 'Select a Stored Rocket File');

            if file ~= 0
    
                filepath = fullfile(path, file);
                app.rocket = load(filepath, "rocketObj").rocketObj;
    
                app.AeroDataButton.Text = app.rocket.name + "_aero.csv";
                app.RocketNameEditField.Value = app.rocket.name;
    
                if ~isempty(app.rocket.totalLength)
                    app.RocketLengthEditField.Value = app.rocket.totalLength;
                end
    
                if ~isempty(app.rocket.outerDiameter)
                    app.RocketDiameterEditField.Value = app.rocket.outerDiameter;
                end
            end


        end

        % Button pushed function: SaveRocketButton
        function SaveRocketButtonPushed(app, event)
            name = string(app.RocketNameEditField.Value);
            path = "TheSixDoF" + filesep + "Inputs" + filesep + "Saved Rockets" + filesep + name + ".mat";

            app.rocket = Rocket(name);

            app.rocket.totalLength = app.RocketLengthEditField.Value;
            app.rocket.outerDiameter = app.RocketDiameterEditField.Value;
            app.rocket.aeroData = name;

            rocketObj = app.rocket;

            save(path, "rocketObj")

            % update the tree node with the rocket object:
            rootNode = uitreenode(app.Tree, 'Text', name);
            expand(rootNode);

            uiconfirm(app.UIFigure, "Rocket saved.", "Congratulations!")

            % after creating a rocket object, auto refresh the plot with
            % changes:
            app.autoRefresh = 1;
        end

        % Value changed function: RocketNameEditField
        function RocketNameChanged(app, event)
            app.RevertButton.Enable = 'on';
            value = app.RocketNameEditField.Value;

            app.UIAxes.Title.String = [value, ' Layout'];

        end

        % Value changed function: ComponentSelectionDropDown
        function ComponentSelectionDropDownValueChanged(app, event)
            type = string(app.ComponentSelectionDropDown.Value);

            propertyList = string(properties(type));
            propertyArray = matlab.metadata.Class.fromName(type).PropertyList;

            delete(app.PropertyGrid.Children);

            nFields = length(propertyList);
            app.PropertyGrid.RowHeight = repmat({'fit'}, 1, nFields);

            % Set up grid layout for PropertyGrid with extra row for the plot
            % gridLayout = uigridlayout(app.PropertyGrid, [nFields + 1, 2]);  % Extra row for the plot
            % gridLayout.RowHeight = [repmat({'fit'}, 1, nFields), '1x'];  % Make the last row flexible for the plot

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

                % Add a plot for "Fins" if it is selected from the component selection
                % if type == "Fins"
                %     % Create an axes for the plot (this will take the last row)
                %     plotAxes = axes(gridLayout);
                %     plotAxes.Layout.Row = nFields + 1;  % Put plot in the last row (nFields + 1)
                %     plotAxes.Layout.Column = [1, 2];    % Span across both columns
                %     %plotAxes.Position = [0.1, 0.1, 0.8, 0.8];  % Adjust plot size within the cell
                %     plot(plotAxes, 1:10, rand(1, 10));  % Example plot (replace with actual data for Fins)
                % end
            end
        end

        % Button pushed function: RevertButton
        function RevertButtonPushed(app, event)
            app.RocketNameEditField.Value = app.rocket.name;
            app.AeroDataButton.Text = app.rocket.name + "_aero.csv";
            app.RocketNameEditField.Value = app.rocket.name;

            if ~isempty(app.rocket.totalLength)
                app.RocketLengthEditField.Value = app.rocket.totalLength;
            end

            if ~isempty(app.rocket.outerDiameter)
                app.RocketDiameterEditField.Value = app.rocket.outerDiameter;
            end
        end

        % Button pushed function: GetWeatherConditionsButton
        function getWeather(app, event)
            % get the weather if the appropriate fields are filled out:

            date = app.DateSelectionDatePicker.Value;

            if isempty(date)
                uialert(app.UIFigure, "Date field is empty!", "Input Error")
            end

        end

        % Callback function: RocketColor
        function BaseColorChanged(app, event)
            value = app.RocketColor.Value;

            if app.autoRefresh
                app.RocketPlotter();
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

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
            app.GridLayout9 = uigridlayout(app.RocketDesignTab);
            app.GridLayout9.ColumnWidth = {'2x', '2x', '2x'};
            app.GridLayout9.RowHeight = {'1.3x', 22, '1x'};
            app.GridLayout9.ColumnSpacing = 5.04001007080078;
            app.GridLayout9.RowSpacing = 5.07499361038208;
            app.GridLayout9.Padding = [5.04001007080078 5.07499361038208 5.04001007080078 5.07499361038208];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout9);
            title(app.UIAxes, 'Rocket Layout', 'Interpreter', 'latex')
            xlabel(app.UIAxes, 'X', 'Interpreter', 'latex')
            ylabel(app.UIAxes, 'Y', 'Interpreter', 'latex')
            zlabel(app.UIAxes, 'Z', 'Interpreter', 'latex')
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = [2 3];

            % Create Panel
            app.Panel = uipanel(app.GridLayout9);
            app.Panel.Layout.Row = [1 3];
            app.Panel.Layout.Column = 1;

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Panel);
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
            app.RocketGrid.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create TotalLengthmLabel
            app.TotalLengthmLabel = uilabel(app.RocketGrid);
            app.TotalLengthmLabel.HorizontalAlignment = 'center';
            app.TotalLengthmLabel.WordWrap = 'on';
            app.TotalLengthmLabel.Layout.Row = 2;
            app.TotalLengthmLabel.Layout.Column = 1;
            app.TotalLengthmLabel.Text = 'Total Length [m]';

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
            app.ButtonPanel.Layout.Row = [8 9];
            app.ButtonPanel.Layout.Column = [1 2];

            % Create ButtonGrid
            app.ButtonGrid = uigridlayout(app.ButtonPanel);

            % Create LoadRocketButton
            app.LoadRocketButton = uibutton(app.ButtonGrid, 'push');
            app.LoadRocketButton.ButtonPushedFcn = createCallbackFcn(app, @LoadRocketButtonPushed, true);
            app.LoadRocketButton.Layout.Row = 1;
            app.LoadRocketButton.Layout.Column = 1;
            app.LoadRocketButton.Text = 'Load';

            % Create SaveRocketButton
            app.SaveRocketButton = uibutton(app.ButtonGrid, 'push');
            app.SaveRocketButton.ButtonPushedFcn = createCallbackFcn(app, @SaveRocketButtonPushed, true);
            app.SaveRocketButton.Layout.Row = 2;
            app.SaveRocketButton.Layout.Column = 1;
            app.SaveRocketButton.Text = 'Save';

            % Create CreateNewRocketButton
            app.CreateNewRocketButton = uibutton(app.ButtonGrid, 'push');
            app.CreateNewRocketButton.Layout.Row = 2;
            app.CreateNewRocketButton.Layout.Column = 2;
            app.CreateNewRocketButton.Text = 'Create New';

            % Create RevertButton
            app.RevertButton = uibutton(app.ButtonGrid, 'push');
            app.RevertButton.ButtonPushedFcn = createCallbackFcn(app, @RevertButtonPushed, true);
            app.RevertButton.Layout.Row = 1;
            app.RevertButton.Layout.Column = 2;
            app.RevertButton.Text = 'Revert';

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
            app.Switchto2D.Layout.Row = 2;
            app.Switchto2D.Layout.Column = 3;
            app.Switchto2D.Text = '2D View';

            % Create UpdatePlotButton
            app.UpdatePlotButton = uibutton(app.GridLayout9, 'push');
            app.UpdatePlotButton.ButtonPushedFcn = createCallbackFcn(app, @UpdatePlotButtonPushed, true);
            app.UpdatePlotButton.Layout.Row = 2;
            app.UpdatePlotButton.Layout.Column = 2;
            app.UpdatePlotButton.Text = 'Update Plot';

            % Create Switchto3D
            app.Switchto3D = uibutton(app.GridLayout9, 'push');
            app.Switchto3D.ButtonPushedFcn = createCallbackFcn(app, @ConvertToThreeD, true);
            app.Switchto3D.Layout.Row = 2;
            app.Switchto3D.Layout.Column = 3;
            app.Switchto3D.Text = '3D View';

            % Create ListBox
            app.ListBox = uilistbox(app.GridLayout9);
            app.ListBox.Items = {};
            app.ListBox.Enable = 'off';
            app.ListBox.Visible = 'off';
            app.ListBox.Layout.Row = 3;
            app.ListBox.Layout.Column = 2;
            app.ListBox.Value = {};

            % Create UITable
            app.UITable = uitable(app.GridLayout9);
            app.UITable.ColumnName = {'Property'; 'Value'};
            app.UITable.RowName = {'Property1, Property2, Property3'};
            app.UITable.Layout.Row = 3;
            app.UITable.Layout.Column = 3;

            % Create Tree
            app.Tree = uitree(app.GridLayout9);
            app.Tree.Layout.Row = 3;
            app.Tree.Layout.Column = 2;

            % Create SimulationTab
            app.SimulationTab = uitab(app.TabGroup);
            app.SimulationTab.Title = 'Simulation';

            % Create Panel_2
            app.Panel_2 = uipanel(app.SimulationTab);
            app.Panel_2.Position = [1 0 260 457];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel_2);
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
            app.DateSelectionDatePicker.Layout.Row = 1;
            app.DateSelectionDatePicker.Layout.Column = 2;

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
            app.LatitudedegEditField.Value = 35.3474;

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
            app.LongitudedegEditField.Value = -117.8091;

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