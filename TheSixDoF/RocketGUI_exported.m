classdef RocketGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RocketDesignTab                matlab.ui.container.Tab
        Switchto2D                     matlab.ui.control.Button
        Switchto3DButton               matlab.ui.control.Button
        FinDesignPanel                 matlab.ui.container.Panel
        GridLayout5                    matlab.ui.container.GridLayout
        TabGroup3                      matlab.ui.container.TabGroup
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
        AirfoilTab                     matlab.ui.container.Tab
        GridLayout7                    matlab.ui.container.GridLayout
        AirfoilGeometryDropDown        matlab.ui.control.DropDown
        AirfoilGeometryDropDownLabel   matlab.ui.control.Label
        FinGraph                       matlab.ui.control.UIAxes
        Panel                          matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        TabGroup2                      matlab.ui.container.TabGroup
        GeometryTab                    matlab.ui.container.Tab
        GridLayout3                    matlab.ui.container.GridLayout
        Switch                         matlab.ui.control.Switch
        DistancefromRearmEditField     matlab.ui.control.NumericEditField
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
        ComponentTab                   matlab.ui.container.Tab
        GridLayout4                    matlab.ui.container.GridLayout
        PlaceComponentButton           matlab.ui.control.Button
        ComponentLocationEditField     matlab.ui.control.NumericEditField
        ComponentLocationEditFieldLabel  matlab.ui.control.Label
        ComponentSelectionDropDown     matlab.ui.control.DropDown
        ComponentSelectionDropDownLabel  matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
        SimulationTab                  matlab.ui.container.Tab
        LocationPlotPanel              matlab.ui.container.Panel
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
        ThreeDPlot;
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

            switch app.NoseConeGeometryDropDown.Value

                case 'Conic'
                    yNose = xNose.*dia./(noseLeng*2);
                case 'Tangent Ogive'
                    R = dia/2;
                    L = noseLeng;
                    rho = (R^2 + L^2) / (2*R);
                    yNose = sqrt(rho^2-(L-xNose).^2) + R - rho;

                case 'Von Karman'
                    theta = linspace(0,pi,50);
                    R = dia/2;
                    L = noseLeng;
                    xNose = L/2 * (1-cos(theta));
                    yNose = R/sqrt(pi) * sqrt(theta-sin(2*theta)/2);
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

                    app.plotFins(numFin)

                end

                % define the standard limits for the plot
                hold(app.UIAxes, 'off')
                xlim(app.UIAxes, [-0.02*leng,leng*1.02])
                axis(app.UIAxes, "equal")

            end

        end

        function Geoplotter(app)

            lat = app.LatitudedegEditField.Value;
            long = app.LongitudedegEditField.Value;

            g = geoaxes(app.LocationPlotPanel);

            geoplot(g, lat, long, 'ro')
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

        function plotFins(app, numFin)

            % get the default fin geometry from the fin plotter function
            [xFin, yFin] = app.FinPlotter;

            % get parameters from user input:
            rearDist = app.DistancefromRearmEditField.Value;
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
            value = app.NoseConeLengthmEditField.Value;
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

        % Value changed function: DistancefromRearmEditField
        function DistFromRearChanged(app, event)
            value = app.DistancefromRearmEditField.Value;
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
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 640 480];

            % Create RocketDesignTab
            app.RocketDesignTab = uitab(app.TabGroup);
            app.RocketDesignTab.Title = 'Rocket Design';

            % Create UIAxes
            app.UIAxes = uiaxes(app.RocketDesignTab);
            title(app.UIAxes, 'Rocket Layout', 'Interpreter', 'latex')
            xlabel(app.UIAxes, 'X', 'Interpreter', 'latex')
            ylabel(app.UIAxes, 'Y', 'Interpreter', 'latex')
            zlabel(app.UIAxes, 'Z', 'Interpreter', 'latex')
            app.UIAxes.Position = [226 264 403 185];

            % Create Panel
            app.Panel = uipanel(app.RocketDesignTab);
            app.Panel.Position = [1 1 213 456];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Panel);
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.GridLayout);
            app.TabGroup2.Layout.Row = [1 7];
            app.TabGroup2.Layout.Column = [1 2];

            % Create GeometryTab
            app.GeometryTab = uitab(app.TabGroup2);
            app.GeometryTab.Title = 'Geometry';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GeometryTab);
            app.GridLayout3.RowHeight = {20, '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

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
            app.NoseConeGeometryDropDown.Items = {'Von Karman', 'Tangent Ogive', 'Conic'};
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
            app.FinNumberSpinnerLabel.Layout.Row = 7;
            app.FinNumberSpinnerLabel.Layout.Column = 1;
            app.FinNumberSpinnerLabel.Text = 'Fin Number';

            % Create FinNumberSpinner
            app.FinNumberSpinner = uispinner(app.GridLayout3);
            app.FinNumberSpinner.Limits = [0 6];
            app.FinNumberSpinner.ValueChangedFcn = createCallbackFcn(app, @FinNumberChanged, true);
            app.FinNumberSpinner.Layout.Row = 7;
            app.FinNumberSpinner.Layout.Column = 2;
            app.FinNumberSpinner.Value = 4;

            % Create FinDesignLabel
            app.FinDesignLabel = uilabel(app.GridLayout3);
            app.FinDesignLabel.HorizontalAlignment = 'center';
            app.FinDesignLabel.FontWeight = 'bold';
            app.FinDesignLabel.Layout.Row = 6;
            app.FinDesignLabel.Layout.Column = [1 2];
            app.FinDesignLabel.Text = 'Fin Design';

            % Create DistancefromRearmEditFieldLabel
            app.DistancefromRearmEditFieldLabel = uilabel(app.GridLayout3);
            app.DistancefromRearmEditFieldLabel.HorizontalAlignment = 'center';
            app.DistancefromRearmEditFieldLabel.WordWrap = 'on';
            app.DistancefromRearmEditFieldLabel.Layout.Row = 8;
            app.DistancefromRearmEditFieldLabel.Layout.Column = 1;
            app.DistancefromRearmEditFieldLabel.Text = 'Distance from Rear [m]';

            % Create DistancefromRearmEditField
            app.DistancefromRearmEditField = uieditfield(app.GridLayout3, 'numeric');
            app.DistancefromRearmEditField.ValueChangedFcn = createCallbackFcn(app, @DistFromRearChanged, true);
            app.DistancefromRearmEditField.Layout.Row = 8;
            app.DistancefromRearmEditField.Layout.Column = 2;
            app.DistancefromRearmEditField.Value = 0.25;

            % Create Switch
            app.Switch = uiswitch(app.GridLayout3, 'slider');
            app.Switch.Items = {'SI', 'Imperial'};
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SiImperialSwitch, true);
            app.Switch.Layout.Row = 1;
            app.Switch.Layout.Column = [1 2];
            app.Switch.Value = 'SI';

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

            % Create ComponentLocationEditFieldLabel
            app.ComponentLocationEditFieldLabel = uilabel(app.GridLayout4);
            app.ComponentLocationEditFieldLabel.HorizontalAlignment = 'center';
            app.ComponentLocationEditFieldLabel.WordWrap = 'on';
            app.ComponentLocationEditFieldLabel.Layout.Row = 2;
            app.ComponentLocationEditFieldLabel.Layout.Column = 1;
            app.ComponentLocationEditFieldLabel.Text = 'Component Location';

            % Create ComponentLocationEditField
            app.ComponentLocationEditField = uieditfield(app.GridLayout4, 'numeric');
            app.ComponentLocationEditField.Layout.Row = 2;
            app.ComponentLocationEditField.Layout.Column = 2;

            % Create PlaceComponentButton
            app.PlaceComponentButton = uibutton(app.GridLayout4, 'push');
            app.PlaceComponentButton.Layout.Row = 7;
            app.PlaceComponentButton.Layout.Column = [1 2];
            app.PlaceComponentButton.Text = 'Place Component';

            % Create FinDesignPanel
            app.FinDesignPanel = uipanel(app.RocketDesignTab);
            app.FinDesignPanel.Title = 'Fin Design';
            app.FinDesignPanel.Position = [227 12 400 218];

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
            app.AirfoilGeometryDropDown.Items = {'NACA 4 Series', 'Hexagonal', 'Wedge'};
            app.AirfoilGeometryDropDown.Layout.Row = 1;
            app.AirfoilGeometryDropDown.Layout.Column = 2;
            app.AirfoilGeometryDropDown.Value = 'NACA 4 Series';

            % Create Switchto3DButton
            app.Switchto3DButton = uibutton(app.RocketDesignTab, 'push');
            app.Switchto3DButton.ButtonPushedFcn = createCallbackFcn(app, @ConvertToThreeD, true);
            app.Switchto3DButton.Position = [378 237 100 22];
            app.Switchto3DButton.Text = 'Switch to 3D';

            % Create Switchto2D
            app.Switchto2D = uibutton(app.RocketDesignTab, 'push');
            app.Switchto2D.ButtonPushedFcn = createCallbackFcn(app, @SwitchToTwoD, true);
            app.Switchto2D.Enable = 'off';
            app.Switchto2D.Visible = 'off';
            app.Switchto2D.Position = [378 237 100 22];
            app.Switchto2D.Text = 'Switch to 2D';

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
            app.LongitudedegEditField.Layout.Row = 3;
            app.LongitudedegEditField.Layout.Column = 2;
            app.LongitudedegEditField.Value = -117.8091;

            % Create LocationPlotPanel
            app.LocationPlotPanel = uipanel(app.SimulationTab);
            app.LocationPlotPanel.TitlePosition = 'centertop';
            app.LocationPlotPanel.Title = 'Location Plot';
            app.LocationPlotPanel.Position = [293 213 319 221];

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