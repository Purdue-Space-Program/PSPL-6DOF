classdef RocketGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RocketDesignTab                matlab.ui.container.Tab
        GridLayout9                    matlab.ui.container.GridLayout
        ComponentBrowserPanel          matlab.ui.container.Panel
        GridLayout10                   matlab.ui.container.GridLayout
        EditComponentListBox           matlab.ui.control.ListBox
        EditComponentListBoxLabel      matlab.ui.control.Label
        FinDesignPanel                 matlab.ui.container.Panel
        GridLayout5                    matlab.ui.container.GridLayout
        TabGroup3                      matlab.ui.container.TabGroup
        AirfoilTab                     matlab.ui.container.Tab
        GridLayout7                    matlab.ui.container.GridLayout
        DoubleAngleLabel               matlab.ui.control.Label
        DoubleAngleEditField           matlab.ui.control.EditField
        NACAEditField                  matlab.ui.control.NumericEditField
        NACAEditFieldLabel             matlab.ui.control.Label
        AirfoilGeometryDropDown        matlab.ui.control.DropDown
        AirfoilGeometryDropDownLabel   matlab.ui.control.Label
        PlanformTab                    matlab.ui.container.Tab
        GridLayout6                    matlab.ui.container.GridLayout
        SweepEditField                 matlab.ui.control.NumericEditField
        SweepEditFieldLabel            matlab.ui.control.Label
        SpanEditField                  matlab.ui.control.NumericEditField
        SpanEditFieldLabel             matlab.ui.control.Label
        TipChordEditField              matlab.ui.control.NumericEditField
        TipChordEditFieldLabel         matlab.ui.control.Label
        RootChordEditField             matlab.ui.control.NumericEditField
        RootChordEditFieldLabel        matlab.ui.control.Label
        FinGraph                       matlab.ui.control.UIAxes
        Switchto2D                     matlab.ui.control.Button
        Switchto3DButton               matlab.ui.control.Button
        Panel                          matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        TabGroup2                      matlab.ui.container.TabGroup
        ComponentTab                   matlab.ui.container.Tab
        GridLayout4                    matlab.ui.container.GridLayout
        ComponentOptionsLabel          matlab.ui.control.Label
        DryMasskgEditField             matlab.ui.control.NumericEditField
        DryMasskgEditFieldLabel        matlab.ui.control.Label
        AddComponentButton             matlab.ui.control.Button
        TankDiametermEditField         matlab.ui.control.NumericEditField
        TankDiametermEditFieldLabel    matlab.ui.control.Label
        TankLengthmEditField           matlab.ui.control.NumericEditField
        TankLengthmEditFieldLabel      matlab.ui.control.Label
        DistancefromNosemEditField     matlab.ui.control.NumericEditField
        DistancefromNosemEditFieldLabel  matlab.ui.control.Label
        ComponentSelectionDropDown     matlab.ui.control.DropDown
        ComponentSelectionDropDownLabel  matlab.ui.control.Label
        GeometryTab                    matlab.ui.container.Tab
        GridLayout3                    matlab.ui.container.GridLayout
        RocketDryMass                  matlab.ui.control.NumericEditField
        DryMasskgEditField_2Label      matlab.ui.control.Label
        Switch                         matlab.ui.control.Switch
        DistancefromRear               matlab.ui.control.NumericEditField
        DistancefromRearmEditFieldLabel  matlab.ui.control.Label
        FinDesignLabel                 matlab.ui.control.Label
        FinNumberSpinner               matlab.ui.control.Spinner
        FinNumberSpinnerLabel          matlab.ui.control.Label
        NoseConeLengthmEditField       matlab.ui.control.NumericEditField
        NoseConeLengthmEditFieldLabel  matlab.ui.control.Label
        NoseConeGeometryDropDown       matlab.ui.control.DropDown
        NoseConeGeometryDropDownLabel  matlab.ui.control.Label
        RocketDiameterEditField        matlab.ui.control.NumericEditField
        RocketDiametermEditFieldLabel  matlab.ui.control.Label
        RocketLengthEditField          matlab.ui.control.NumericEditField
        RocketLengthmEditFieldLabel    matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
        SimulationTab                  matlab.ui.container.Tab
        LaunchLocationPanel            matlab.ui.container.Panel
        Panel_2                        matlab.ui.container.Panel
        GridLayout2                    matlab.ui.container.GridLayout
        LongitudedegEditField          matlab.ui.control.NumericEditField
        LongitudedegEditFieldLabel     matlab.ui.control.Label
        LatitudedegEditField           matlab.ui.control.NumericEditField
        LatitudedegEditFieldLabel      matlab.ui.control.Label
        DateSelectionDatePicker        matlab.ui.control.DatePicker
        DateSelectionDatePickerLabel   matlab.ui.control.Label
    end


    properties (Access = private)
        lineColor; % line color for plots
        ThreeDPlot; % flag to turn 3d plotting on and off
        ComponentList = ""; % a list of the components on the rocket
        ComponentDetails;
    end

    methods (Access = private)

        function RocketPlotter(app)

            % create general parameters for the rocket:
            leng = app.RocketLengthEditField.Value;
            noseLeng = app.NoseConeLengthmEditField.Value;
            dia = app.RocketDiameterEditField.Value;

            % define the geometry over the nose cone:
            xNose = linspace(0,noseLeng, 50);

            % change the y profile based on the selection.

            R = dia/2;

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

            % if the plot is in 3d, modify it
            if app.ThreeDPlot

                % [X, Y, Z] = cylinder(y);
                %
                % surf(app.UIAxes, X,Y,Z, "FaceColor",'r', 'EdgeAlpha',0);
                % axis(app.UIAxes, "equal")
            else

                % plot the base body of the rocket:
                plot(app.UIAxes, x,y/2, app.lineColor)
                hold(app.UIAxes, "on")
                plot(app.UIAxes, x,-y/2, app.lineColor)
                plot(app.UIAxes, [x(end),x(end)], [dia/2,-dia/2], app.lineColor)


                % plot the fins of the rocket, if they exist

                numFin = app.FinNumberSpinner.Value;

                if numFin ~= 0
                    app.PlotFins(numFin)
                end

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

        function PlotFins(app, numFin)

            % get the default fin geometry from the fin plotter function
            [xFin, yFin] = app.FinPlotter;

            % get parameters from user input:
            rearDist = app.DistancefromRear.Value;
            dia = app.RocketDiameterEditField.Value;
            leng = app.RocketLengthEditField.Value;

            % first, check which fins should be plotted based on occlusion
            % (manually for now lmao, don't know how to write this
            % programmatically)

            switch numFin
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
                theta = (2*pi)/numFin * (idx-1);

                % generate an array of matrices for the projected fins.
                xFinShifted = xFin + leng - rearDist;
                yFinProjection = yFin*sin(theta);

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

            leng = numel(app.EditComponentListBox.Items);

            for idx = 1:leng

                itemData = app.ComponentDetails(idx,:);

                componentType = itemData(1);

                switch componentType
                    case 1 % tank

                        % for the tank, the data is in the order:
                        % Component Type
                        % Length
                        % Diameter
                        % Dist from Nose
                        % Dry Mass

                        length = itemData(2);
                        rad = itemData(3)/2;
                        dist = itemData(4);

                        xTank = [dist, dist, dist+length, dist+length, dist];
                        yTank = [-rad, rad, rad, -rad, -rad];

                        plot(app.UIAxes, xTank, yTank, 'b-')
                end
            end
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

            app.RocketPlotter();
            app.Geoplotter();
            app.FinPlotter();
            app.AirfoilGeoChanged();
        end

        % Value changed function: RocketLengthEditField
        function RocketLengthChanged(app, event)
            rocketLeng = app.RocketLengthEditField.Value;

            if rocketLeng < app.NoseConeLengthmEditField.Value
                uialert(app.UIFigure, "The rocket cannot be " + ...
                    "shorter than the nose cone length!", "Input Error", "Icon","error");
                return

            end

            app.RocketPlotter();

        end

        % Value changed function: LatitudedegEditField
        function LatitudeChanged(app, event)
            value = app.LatitudedegEditField.Value;
            app.Geoplotter();
        end

        % Value changed function: RocketDiameterEditField
        function RocketDiaChanged(app, event)
            value = app.RocketDiameterEditField.Value;
            app.RocketPlotter();
        end

        % Value changed function: NoseConeLengthmEditField
        function NoseCoseLengthChanged(app, event)
            noseLeng = app.NoseConeLengthmEditField.Value;

            if noseLeng >= app.RocketLengthEditField.Value
                uialert(app.UIFigure, "Nose cone cannot be longer than rocket body!", ...
                    "Dimension Error")
                return
            end

            app.RocketPlotter();
        end

        % Value changed function: NoseConeGeometryDropDown
        function NoseConeTypeChanged(app, event)
            app.RocketPlotter();
        end

        % Value changed function: RootChordEditField
        function RootChordChanged(app, event)
            app.FinPlotter();
            app.RocketPlotter();
        end

        % Value changed function: TipChordEditField
        function TipChordChanged(app, event)
            app.FinPlotter();
            app.RocketPlotter();
        end

        % Value changed function: SpanEditField
        function SpanChanged(app, event)
            app.FinPlotter();
            app.RocketPlotter();
        end

        % Value changed function: SweepEditField
        function SweepChanged(app, event)
            app.FinPlotter();
            app.RocketPlotter();
        end

        % Value changed function: FinNumberSpinner
        function FinNumberChanged(app, event)
            value = app.FinNumberSpinner.Value;
            app.RocketPlotter();
        end

        % Value changed function: DistancefromRear
        function DistFromRearChanged(app, event)
            value = app.DistancefromRear.Value;
            app.RocketPlotter();
        end

        % Button pushed function: Switchto3DButton
        function ConvertToThreeD(app, event)
            app.Switchto2D.Enable = "on";
            app.Switchto2D.Visible = 'on';
            app.Switchto3DButton.Enable = 'off';
            app.Switchto3DButton.Visible = "off";

            app.ThreeDPlot = 1;

            app.RocketPlotter();
        end

        % Button pushed function: Switchto2D
        function SwitchToTwoD(app, event)
            app.Switchto3DButton.Enable = "on";
            app.Switchto3DButton.Visible = 'on';
            app.Switchto2D.Enable = 'off';
            app.Switchto2D.Visible = 'off';

            app.ThreeDPlot = 0;
            app.RocketPlotter();
        end

        % Value changed function: Switch
        function SiImperialSwitch(app, event)
            value = app.Switch.Value;

            if strcmp(value, "Imperial")
            app.RocketDiametermEditFieldLabel.Text = 'Rocket Diameter [in]';
            app.RocketLengthmEditFieldLabel.Text = 'Rocket Length [in]';
            app.NoseConeLengthmEditFieldLabel.Text = 'Nose Cone Length [in]';
            app.DistancefromRearmEditFieldLabel.Text = 'Distance from Rear [in]';
            else
            app.RocketDiametermEditFieldLabel.Text = 'Rocket Diameter [m]';
            app.RocketLengthmEditFieldLabel.Text = 'Rocket Length [m]';
            app.NoseConeLengthmEditFieldLabel.Text = 'Nose Cone Length [m]'; 
            app.DistancefromRearmEditFieldLabel.Text = 'Distance from Rear [m]';

            end
        end

        % Value changed function: AirfoilGeometryDropDown
        function AirfoilGeoChanged(app, event)
            value = app.AirfoilGeometryDropDown.Value;
            
        switch value

            case 'NACA 4 Series'

                % switch on the NACA number input
                app.NACAEditField.Enable = "on";
                app.NACAEditField.Visible = "on";
                app.NACAEditFieldLabel.Enable = "on";
                app.NACAEditFieldLabel.Visible = "on";

                app.DoubleAngleEditField.Enable = 'off';
                app.DoubleAngleEditField.Visible = 'off';
                app.DoubleAngleLabel.Visible = 'off';
                app.DoubleAngleLabel.Enable = 'off';

            case 'Double Wedge'
                app.DoubleAngleEditField.Enable = 'on';
                app.DoubleAngleEditField.Visible = 'on';
                app.DoubleAngleLabel.Visible = 'on';
                app.DoubleAngleLabel.Enable = 'on';


                app.NACAEditField.Enable = "off";
                app.NACAEditField.Visible = "off";
                app.NACAEditFieldLabel.Enable = "off";
                app.NACAEditFieldLabel.Visible = "off";

        end
        end

        % Callback function
        function GeoplotZoomChanged(app, event)
            value = app.ZoomLevelSlider.Value;
            
            app.Geoplotter();
        end

        % Callback function
        function GeoplotClicked(app, event)
                clickedPoint = ax.CurrentPoint; % [x, y] coordinates where clicked

                % Convert the clicked point to geographical coordinates (latitude, longitude)
                lat = clickedPoint(1,2);
                lon = clickedPoint(1,1);
            
                % Display the coordinates
                disp(['Latitude: ', num2str(lat), ', Longitude: ', num2str(lon)]);
        end

        % Value changed function: LongitudedegEditField
        function longitudeChanged(app, event)
            value = app.LongitudedegEditField.Value;
            
            app.Geoplotter();
        end

        % Button pushed function: AddComponentButton
        function AddComponent(app, event)
            name = string(app.ComponentSelectionDropDown.Value);

            % check if the identifier is unique
            
            if contains(app.ComponentList, name)
                % do something if an element already exists
            end

            % if it is empty
            if (app.ComponentList == "")
                app.ComponentList = name;
            else
                app.ComponentList = [app.ComponentList, name];
            end

            app.EditComponentListBox.Items = app.ComponentList;

            leng = numel(app.EditComponentListBox.Items);

            switch name
                case "Tank"
                    app.ComponentDetails(leng,:) = [1,app.TankLengthmEditField.Value, ...
                    app.TankDiametermEditField.Value, ...
                    app.DistancefromNosemEditField.Value, app.DryMasskgEditField.Value];

                otherwise
                    app.EditComponentListBox.ItemsData = {};
            end 

            app.RocketPlotter();



        end

        % Selection change function: TabGroup2
        function TabGroupChanged(app, event)
            selectedTab = app.TabGroup2.SelectedTab.Title;

            if strcmp(selectedTab, "Geometry")
                app.ComponentBrowserPanel.Visible = 'off';
                app.ComponentBrowserPanel.Enable = 'off';
            else
                app.ComponentBrowserPanel.Visible = 'on';
                app.ComponentBrowserPanel.Enable = 'on';
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
            app.TabGroup.Position = [2 2 640 480];

            % Create RocketDesignTab
            app.RocketDesignTab = uitab(app.TabGroup);
            app.RocketDesignTab.Title = 'Rocket Design';

            % Create GridLayout9
            app.GridLayout9 = uigridlayout(app.RocketDesignTab);
            app.GridLayout9.ColumnWidth = {'fit', '1.08x', 86, '1.08x'};
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
            app.UIAxes.Layout.Column = [2 4];

            % Create FinDesignPanel
            app.FinDesignPanel = uipanel(app.GridLayout9);
            app.FinDesignPanel.Title = 'Fin Design';
            app.FinDesignPanel.Layout.Row = 3;
            app.FinDesignPanel.Layout.Column = [2 4];

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.FinDesignPanel);
            app.GridLayout5.ColumnWidth = {100, 110, '2x'};
            app.GridLayout5.RowHeight = {'1x', '1x', '1x', '1x', '1x'};

            % Create FinGraph
            app.FinGraph = uiaxes(app.GridLayout5);
            title(app.FinGraph, 'Fin Design', 'Interpreter', 'latex')
            xlabel(app.FinGraph, 'X', 'Interpreter', 'latex')
            ylabel(app.FinGraph, 'Y', 'Interpreter', 'latex')
            zlabel(app.FinGraph, 'Z', 'Interpreter', 'latex')
            app.FinGraph.Layout.Row = [1 5];
            app.FinGraph.Layout.Column = 3;

            % Create TabGroup3
            app.TabGroup3 = uitabgroup(app.GridLayout5);
            app.TabGroup3.Layout.Row = [1 5];
            app.TabGroup3.Layout.Column = [1 2];

            % Create PlanformTab
            app.PlanformTab = uitab(app.TabGroup3);
            app.PlanformTab.Title = 'Planform';

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.PlanformTab);
            app.GridLayout6.ColumnWidth = {90, 90};
            app.GridLayout6.RowHeight = {'1x', '1x', '1x', '1x'};

            % Create RootChordEditFieldLabel
            app.RootChordEditFieldLabel = uilabel(app.GridLayout6);
            app.RootChordEditFieldLabel.HorizontalAlignment = 'center';
            app.RootChordEditFieldLabel.Layout.Row = 1;
            app.RootChordEditFieldLabel.Layout.Column = 1;
            app.RootChordEditFieldLabel.Text = 'Root Chord';

            % Create RootChordEditField
            app.RootChordEditField = uieditfield(app.GridLayout6, 'numeric');
            app.RootChordEditField.Limits = [0 Inf];
            app.RootChordEditField.ValueChangedFcn = createCallbackFcn(app, @RootChordChanged, true);
            app.RootChordEditField.Layout.Row = 1;
            app.RootChordEditField.Layout.Column = 2;
            app.RootChordEditField.Value = 0.2;

            % Create TipChordEditFieldLabel
            app.TipChordEditFieldLabel = uilabel(app.GridLayout6);
            app.TipChordEditFieldLabel.HorizontalAlignment = 'center';
            app.TipChordEditFieldLabel.Layout.Row = 2;
            app.TipChordEditFieldLabel.Layout.Column = 1;
            app.TipChordEditFieldLabel.Text = 'Tip Chord';

            % Create TipChordEditField
            app.TipChordEditField = uieditfield(app.GridLayout6, 'numeric');
            app.TipChordEditField.Limits = [0 Inf];
            app.TipChordEditField.ValueChangedFcn = createCallbackFcn(app, @TipChordChanged, true);
            app.TipChordEditField.Layout.Row = 2;
            app.TipChordEditField.Layout.Column = 2;
            app.TipChordEditField.Value = 0.1;

            % Create SpanEditFieldLabel
            app.SpanEditFieldLabel = uilabel(app.GridLayout6);
            app.SpanEditFieldLabel.HorizontalAlignment = 'center';
            app.SpanEditFieldLabel.Layout.Row = 3;
            app.SpanEditFieldLabel.Layout.Column = 1;
            app.SpanEditFieldLabel.Text = 'Span';

            % Create SpanEditField
            app.SpanEditField = uieditfield(app.GridLayout6, 'numeric');
            app.SpanEditField.Limits = [0 Inf];
            app.SpanEditField.ValueChangedFcn = createCallbackFcn(app, @SpanChanged, true);
            app.SpanEditField.Layout.Row = 3;
            app.SpanEditField.Layout.Column = 2;
            app.SpanEditField.Value = 0.1;

            % Create SweepEditFieldLabel
            app.SweepEditFieldLabel = uilabel(app.GridLayout6);
            app.SweepEditFieldLabel.HorizontalAlignment = 'center';
            app.SweepEditFieldLabel.Layout.Row = 4;
            app.SweepEditFieldLabel.Layout.Column = 1;
            app.SweepEditFieldLabel.Text = 'Sweep';

            % Create SweepEditField
            app.SweepEditField = uieditfield(app.GridLayout6, 'numeric');
            app.SweepEditField.ValueChangedFcn = createCallbackFcn(app, @SweepChanged, true);
            app.SweepEditField.Layout.Row = 4;
            app.SweepEditField.Layout.Column = 2;
            app.SweepEditField.Value = 0.15;

            % Create AirfoilTab
            app.AirfoilTab = uitab(app.TabGroup3);
            app.AirfoilTab.Title = 'Airfoil';

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.AirfoilTab);
            app.GridLayout7.RowHeight = {'1x', '1x', '1x'};

            % Create AirfoilGeometryDropDownLabel
            app.AirfoilGeometryDropDownLabel = uilabel(app.GridLayout7);
            app.AirfoilGeometryDropDownLabel.HorizontalAlignment = 'right';
            app.AirfoilGeometryDropDownLabel.Layout.Row = 1;
            app.AirfoilGeometryDropDownLabel.Layout.Column = 1;
            app.AirfoilGeometryDropDownLabel.Text = 'Airfoil Geometry';

            % Create AirfoilGeometryDropDown
            app.AirfoilGeometryDropDown = uidropdown(app.GridLayout7);
            app.AirfoilGeometryDropDown.Items = {'NACA 4 Series', 'Hexagonal', 'Double Wedge'};
            app.AirfoilGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @AirfoilGeoChanged, true);
            app.AirfoilGeometryDropDown.Layout.Row = 1;
            app.AirfoilGeometryDropDown.Layout.Column = 2;
            app.AirfoilGeometryDropDown.Value = 'NACA 4 Series';

            % Create NACAEditFieldLabel
            app.NACAEditFieldLabel = uilabel(app.GridLayout7);
            app.NACAEditFieldLabel.HorizontalAlignment = 'right';
            app.NACAEditFieldLabel.Enable = 'off';
            app.NACAEditFieldLabel.Visible = 'off';
            app.NACAEditFieldLabel.Layout.Row = 2;
            app.NACAEditFieldLabel.Layout.Column = 1;
            app.NACAEditFieldLabel.Text = 'NACA';

            % Create NACAEditField
            app.NACAEditField = uieditfield(app.GridLayout7, 'numeric');
            app.NACAEditField.ValueDisplayFormat = '%04d';
            app.NACAEditField.Editable = 'off';
            app.NACAEditField.Enable = 'off';
            app.NACAEditField.Visible = 'off';
            app.NACAEditField.Layout.Row = 2;
            app.NACAEditField.Layout.Column = 2;

            % Create DoubleAngleEditField
            app.DoubleAngleEditField = uieditfield(app.GridLayout7, 'text');
            app.DoubleAngleEditField.Editable = 'off';
            app.DoubleAngleEditField.Enable = 'off';
            app.DoubleAngleEditField.Visible = 'off';
            app.DoubleAngleEditField.Layout.Row = 2;
            app.DoubleAngleEditField.Layout.Column = 2;

            % Create DoubleAngleLabel
            app.DoubleAngleLabel = uilabel(app.GridLayout7);
            app.DoubleAngleLabel.HorizontalAlignment = 'center';
            app.DoubleAngleLabel.Enable = 'off';
            app.DoubleAngleLabel.Visible = 'off';
            app.DoubleAngleLabel.Layout.Row = 2;
            app.DoubleAngleLabel.Layout.Column = 1;
            app.DoubleAngleLabel.Text = 'Double Angle';

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

            % Create GeometryTab
            app.GeometryTab = uitab(app.TabGroup2);
            app.GeometryTab.Title = 'Geometry';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GeometryTab);
            app.GridLayout3.RowHeight = {20, '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create RocketLengthmEditFieldLabel
            app.RocketLengthmEditFieldLabel = uilabel(app.GridLayout3);
            app.RocketLengthmEditFieldLabel.HorizontalAlignment = 'center';
            app.RocketLengthmEditFieldLabel.WordWrap = 'on';
            app.RocketLengthmEditFieldLabel.Layout.Row = 2;
            app.RocketLengthmEditFieldLabel.Layout.Column = 1;
            app.RocketLengthmEditFieldLabel.Text = 'Rocket Length [m]';

            % Create RocketLengthEditField
            app.RocketLengthEditField = uieditfield(app.GridLayout3, 'numeric');
            app.RocketLengthEditField.ValueChangedFcn = createCallbackFcn(app, @RocketLengthChanged, true);
            app.RocketLengthEditField.HorizontalAlignment = 'left';
            app.RocketLengthEditField.Layout.Row = 2;
            app.RocketLengthEditField.Layout.Column = 2;
            app.RocketLengthEditField.Value = 5;

            % Create RocketDiametermEditFieldLabel
            app.RocketDiametermEditFieldLabel = uilabel(app.GridLayout3);
            app.RocketDiametermEditFieldLabel.HorizontalAlignment = 'center';
            app.RocketDiametermEditFieldLabel.WordWrap = 'on';
            app.RocketDiametermEditFieldLabel.Layout.Row = 3;
            app.RocketDiametermEditFieldLabel.Layout.Column = 1;
            app.RocketDiametermEditFieldLabel.Text = 'Rocket Diameter [m]';

            % Create RocketDiameterEditField
            app.RocketDiameterEditField = uieditfield(app.GridLayout3, 'numeric');
            app.RocketDiameterEditField.ValueChangedFcn = createCallbackFcn(app, @RocketDiaChanged, true);
            app.RocketDiameterEditField.HorizontalAlignment = 'left';
            app.RocketDiameterEditField.Layout.Row = 3;
            app.RocketDiameterEditField.Layout.Column = 2;
            app.RocketDiameterEditField.Value = 0.2;

            % Create NoseConeGeometryDropDownLabel
            app.NoseConeGeometryDropDownLabel = uilabel(app.GridLayout3);
            app.NoseConeGeometryDropDownLabel.HorizontalAlignment = 'center';
            app.NoseConeGeometryDropDownLabel.WordWrap = 'on';
            app.NoseConeGeometryDropDownLabel.Layout.Row = 4;
            app.NoseConeGeometryDropDownLabel.Layout.Column = 1;
            app.NoseConeGeometryDropDownLabel.Text = 'Nose Cone Geometry';

            % Create NoseConeGeometryDropDown
            app.NoseConeGeometryDropDown = uidropdown(app.GridLayout3);
            app.NoseConeGeometryDropDown.Items = {'Von Karman', 'Tangent Ogive', 'Conic', 'Elliptical'};
            app.NoseConeGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @NoseConeTypeChanged, true);
            app.NoseConeGeometryDropDown.Layout.Row = 4;
            app.NoseConeGeometryDropDown.Layout.Column = 2;
            app.NoseConeGeometryDropDown.Value = 'Von Karman';

            % Create NoseConeLengthmEditFieldLabel
            app.NoseConeLengthmEditFieldLabel = uilabel(app.GridLayout3);
            app.NoseConeLengthmEditFieldLabel.HorizontalAlignment = 'center';
            app.NoseConeLengthmEditFieldLabel.WordWrap = 'on';
            app.NoseConeLengthmEditFieldLabel.Layout.Row = 5;
            app.NoseConeLengthmEditFieldLabel.Layout.Column = 1;
            app.NoseConeLengthmEditFieldLabel.Text = 'Nose Cone Length [m]';

            % Create NoseConeLengthmEditField
            app.NoseConeLengthmEditField = uieditfield(app.GridLayout3, 'numeric');
            app.NoseConeLengthmEditField.Limits = [0 Inf];
            app.NoseConeLengthmEditField.ValueChangedFcn = createCallbackFcn(app, @NoseCoseLengthChanged, true);
            app.NoseConeLengthmEditField.Layout.Row = 5;
            app.NoseConeLengthmEditField.Layout.Column = 2;
            app.NoseConeLengthmEditField.Value = 0.5;

            % Create FinNumberSpinnerLabel
            app.FinNumberSpinnerLabel = uilabel(app.GridLayout3);
            app.FinNumberSpinnerLabel.HorizontalAlignment = 'center';
            app.FinNumberSpinnerLabel.Layout.Row = 8;
            app.FinNumberSpinnerLabel.Layout.Column = 1;
            app.FinNumberSpinnerLabel.Text = 'Fin Number';

            % Create FinNumberSpinner
            app.FinNumberSpinner = uispinner(app.GridLayout3);
            app.FinNumberSpinner.Limits = [0 6];
            app.FinNumberSpinner.ValueChangedFcn = createCallbackFcn(app, @FinNumberChanged, true);
            app.FinNumberSpinner.Layout.Row = 8;
            app.FinNumberSpinner.Layout.Column = 2;
            app.FinNumberSpinner.Value = 4;

            % Create FinDesignLabel
            app.FinDesignLabel = uilabel(app.GridLayout3);
            app.FinDesignLabel.HorizontalAlignment = 'center';
            app.FinDesignLabel.FontWeight = 'bold';
            app.FinDesignLabel.Layout.Row = 7;
            app.FinDesignLabel.Layout.Column = [1 2];
            app.FinDesignLabel.Text = 'Fin Design';

            % Create DistancefromRearmEditFieldLabel
            app.DistancefromRearmEditFieldLabel = uilabel(app.GridLayout3);
            app.DistancefromRearmEditFieldLabel.HorizontalAlignment = 'center';
            app.DistancefromRearmEditFieldLabel.WordWrap = 'on';
            app.DistancefromRearmEditFieldLabel.Layout.Row = 9;
            app.DistancefromRearmEditFieldLabel.Layout.Column = 1;
            app.DistancefromRearmEditFieldLabel.Text = 'Distance from Rear [m]';

            % Create DistancefromRear
            app.DistancefromRear = uieditfield(app.GridLayout3, 'numeric');
            app.DistancefromRear.ValueChangedFcn = createCallbackFcn(app, @DistFromRearChanged, true);
            app.DistancefromRear.Layout.Row = 9;
            app.DistancefromRear.Layout.Column = 2;
            app.DistancefromRear.Value = 0.25;

            % Create Switch
            app.Switch = uiswitch(app.GridLayout3, 'slider');
            app.Switch.Items = {'SI', 'Imperial'};
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SiImperialSwitch, true);
            app.Switch.Layout.Row = 1;
            app.Switch.Layout.Column = [1 2];
            app.Switch.Value = 'SI';

            % Create DryMasskgEditField_2Label
            app.DryMasskgEditField_2Label = uilabel(app.GridLayout3);
            app.DryMasskgEditField_2Label.HorizontalAlignment = 'right';
            app.DryMasskgEditField_2Label.Layout.Row = 6;
            app.DryMasskgEditField_2Label.Layout.Column = 1;
            app.DryMasskgEditField_2Label.Text = 'Dry Mass [kg]';

            % Create RocketDryMass
            app.RocketDryMass = uieditfield(app.GridLayout3, 'numeric');
            app.RocketDryMass.Layout.Row = 6;
            app.RocketDryMass.Layout.Column = 2;

            % Create ComponentTab
            app.ComponentTab = uitab(app.TabGroup2);
            app.ComponentTab.Title = 'Component';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.ComponentTab);
            app.GridLayout4.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create ComponentSelectionDropDownLabel
            app.ComponentSelectionDropDownLabel = uilabel(app.GridLayout4);
            app.ComponentSelectionDropDownLabel.HorizontalAlignment = 'center';
            app.ComponentSelectionDropDownLabel.WordWrap = 'on';
            app.ComponentSelectionDropDownLabel.Layout.Row = 1;
            app.ComponentSelectionDropDownLabel.Layout.Column = 1;
            app.ComponentSelectionDropDownLabel.Text = 'Component Selection';

            % Create ComponentSelectionDropDown
            app.ComponentSelectionDropDown = uidropdown(app.GridLayout4);
            app.ComponentSelectionDropDown.Items = {'Tank', 'Engine', 'Point Mass', 'Sensor'};
            app.ComponentSelectionDropDown.Layout.Row = 1;
            app.ComponentSelectionDropDown.Layout.Column = 2;
            app.ComponentSelectionDropDown.Value = 'Tank';

            % Create DistancefromNosemEditFieldLabel
            app.DistancefromNosemEditFieldLabel = uilabel(app.GridLayout4);
            app.DistancefromNosemEditFieldLabel.HorizontalAlignment = 'center';
            app.DistancefromNosemEditFieldLabel.WordWrap = 'on';
            app.DistancefromNosemEditFieldLabel.Layout.Row = 5;
            app.DistancefromNosemEditFieldLabel.Layout.Column = 1;
            app.DistancefromNosemEditFieldLabel.Text = 'Distance from Nose [m]';

            % Create DistancefromNosemEditField
            app.DistancefromNosemEditField = uieditfield(app.GridLayout4, 'numeric');
            app.DistancefromNosemEditField.Layout.Row = 5;
            app.DistancefromNosemEditField.Layout.Column = 2;

            % Create TankLengthmEditFieldLabel
            app.TankLengthmEditFieldLabel = uilabel(app.GridLayout4);
            app.TankLengthmEditFieldLabel.HorizontalAlignment = 'center';
            app.TankLengthmEditFieldLabel.WordWrap = 'on';
            app.TankLengthmEditFieldLabel.Layout.Row = 3;
            app.TankLengthmEditFieldLabel.Layout.Column = 1;
            app.TankLengthmEditFieldLabel.Text = 'Tank Length [m]';

            % Create TankLengthmEditField
            app.TankLengthmEditField = uieditfield(app.GridLayout4, 'numeric');
            app.TankLengthmEditField.Limits = [0 Inf];
            app.TankLengthmEditField.Layout.Row = 3;
            app.TankLengthmEditField.Layout.Column = 2;

            % Create TankDiametermEditFieldLabel
            app.TankDiametermEditFieldLabel = uilabel(app.GridLayout4);
            app.TankDiametermEditFieldLabel.HorizontalAlignment = 'center';
            app.TankDiametermEditFieldLabel.WordWrap = 'on';
            app.TankDiametermEditFieldLabel.Layout.Row = 4;
            app.TankDiametermEditFieldLabel.Layout.Column = 1;
            app.TankDiametermEditFieldLabel.Text = 'Tank Diameter [m]';

            % Create TankDiametermEditField
            app.TankDiametermEditField = uieditfield(app.GridLayout4, 'numeric');
            app.TankDiametermEditField.Limits = [0 Inf];
            app.TankDiametermEditField.Layout.Row = 4;
            app.TankDiametermEditField.Layout.Column = 2;

            % Create AddComponentButton
            app.AddComponentButton = uibutton(app.GridLayout4, 'push');
            app.AddComponentButton.ButtonPushedFcn = createCallbackFcn(app, @AddComponent, true);
            app.AddComponentButton.Layout.Row = 7;
            app.AddComponentButton.Layout.Column = [1 2];
            app.AddComponentButton.Text = 'Add Component';

            % Create DryMasskgEditFieldLabel
            app.DryMasskgEditFieldLabel = uilabel(app.GridLayout4);
            app.DryMasskgEditFieldLabel.HorizontalAlignment = 'center';
            app.DryMasskgEditFieldLabel.WordWrap = 'on';
            app.DryMasskgEditFieldLabel.Layout.Row = 6;
            app.DryMasskgEditFieldLabel.Layout.Column = 1;
            app.DryMasskgEditFieldLabel.Text = 'Dry Mass [kg]';

            % Create DryMasskgEditField
            app.DryMasskgEditField = uieditfield(app.GridLayout4, 'numeric');
            app.DryMasskgEditField.Limits = [0 Inf];
            app.DryMasskgEditField.Layout.Row = 6;
            app.DryMasskgEditField.Layout.Column = 2;

            % Create ComponentOptionsLabel
            app.ComponentOptionsLabel = uilabel(app.GridLayout4);
            app.ComponentOptionsLabel.HorizontalAlignment = 'center';
            app.ComponentOptionsLabel.FontWeight = 'bold';
            app.ComponentOptionsLabel.Layout.Row = 2;
            app.ComponentOptionsLabel.Layout.Column = [1 2];
            app.ComponentOptionsLabel.Text = 'Component Options';

            % Create Switchto3DButton
            app.Switchto3DButton = uibutton(app.GridLayout9, 'push');
            app.Switchto3DButton.ButtonPushedFcn = createCallbackFcn(app, @ConvertToThreeD, true);
            app.Switchto3DButton.Layout.Row = 2;
            app.Switchto3DButton.Layout.Column = 3;
            app.Switchto3DButton.Text = 'Switch to 3D';

            % Create Switchto2D
            app.Switchto2D = uibutton(app.GridLayout9, 'push');
            app.Switchto2D.ButtonPushedFcn = createCallbackFcn(app, @SwitchToTwoD, true);
            app.Switchto2D.Enable = 'off';
            app.Switchto2D.Visible = 'off';
            app.Switchto2D.Layout.Row = 2;
            app.Switchto2D.Layout.Column = 3;
            app.Switchto2D.Text = 'Switch to 2D';

            % Create ComponentBrowserPanel
            app.ComponentBrowserPanel = uipanel(app.GridLayout9);
            app.ComponentBrowserPanel.Enable = 'off';
            app.ComponentBrowserPanel.TitlePosition = 'centertop';
            app.ComponentBrowserPanel.Title = 'Component Browser';
            app.ComponentBrowserPanel.Visible = 'off';
            app.ComponentBrowserPanel.Layout.Row = 3;
            app.ComponentBrowserPanel.Layout.Column = [2 4];
            app.ComponentBrowserPanel.FontWeight = 'bold';
            app.ComponentBrowserPanel.FontSize = 14;

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.ComponentBrowserPanel);
            app.GridLayout10.ColumnWidth = {'0.4x', '1x'};
            app.GridLayout10.RowHeight = {28, '1x'};

            % Create EditComponentListBoxLabel
            app.EditComponentListBoxLabel = uilabel(app.GridLayout10);
            app.EditComponentListBoxLabel.HorizontalAlignment = 'right';
            app.EditComponentListBoxLabel.Layout.Row = 1;
            app.EditComponentListBoxLabel.Layout.Column = 1;
            app.EditComponentListBoxLabel.Text = 'Edit Component';

            % Create EditComponentListBox
            app.EditComponentListBox = uilistbox(app.GridLayout10);
            app.EditComponentListBox.Items = {};
            app.EditComponentListBox.Layout.Row = 2;
            app.EditComponentListBox.Layout.Column = 1;
            app.EditComponentListBox.Value = {};

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