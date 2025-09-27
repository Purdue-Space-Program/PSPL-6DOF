classdef RocketGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RocketDesignTab                matlab.ui.container.Tab
        Panel                          matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        RocketDiametermEditField       matlab.ui.control.NumericEditField
        RocketDiametermEditFieldLabel  matlab.ui.control.Label
        RocketLengthmEditField         matlab.ui.control.NumericEditField
        RocketLengthmEditFieldLabel    matlab.ui.control.Label
        ComponentLocationEditField     matlab.ui.control.NumericEditField
        ComponentLocationEditFieldLabel  matlab.ui.control.Label
        ComponentSelectionDropDown     matlab.ui.control.DropDown
        ComponentSelectionDropDownLabel  matlab.ui.control.Label
        PlaceComponentButton           matlab.ui.control.Button
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

    
    methods (Access = private)
        
        function RocketPlotter(app)

            x = linspace(0,1);
            y = x.^2;

            plot(app.UIAxes, x,y)
        end
        
        function Geoplotter(app)

            lat = app.LatitudedegEditField.Value;
            long = app.LongitudedegEditField.Value;

            g = geoaxes(app.LocationPlotPanel);

            geoplot(g, lat, long, 'ro')

            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.RocketPlotter
            app.Geoplotter
        end

        % Value changed function: RocketLengthmEditField
        function RocketLengthChanged(app, event)
            value = app.RocketLengthmEditField.Value;
            
        end

        % Value changed function: LatitudedegEditField
        function LatitudeChanged(app, event)
            value = app.LatitudedegEditField.Value;
            app.Geoplotter();
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

            % Create PlaceComponentButton
            app.PlaceComponentButton = uibutton(app.GridLayout, 'push');
            app.PlaceComponentButton.Layout.Row = 7;
            app.PlaceComponentButton.Layout.Column = [1 2];
            app.PlaceComponentButton.Text = 'Place Component';

            % Create ComponentSelectionDropDownLabel
            app.ComponentSelectionDropDownLabel = uilabel(app.GridLayout);
            app.ComponentSelectionDropDownLabel.HorizontalAlignment = 'center';
            app.ComponentSelectionDropDownLabel.WordWrap = 'on';
            app.ComponentSelectionDropDownLabel.Layout.Row = 3;
            app.ComponentSelectionDropDownLabel.Layout.Column = 1;
            app.ComponentSelectionDropDownLabel.Text = 'Component Selection';

            % Create ComponentSelectionDropDown
            app.ComponentSelectionDropDown = uidropdown(app.GridLayout);
            app.ComponentSelectionDropDown.Items = {'Tank', 'Engine', 'Point Mass', 'Sensor'};
            app.ComponentSelectionDropDown.Layout.Row = 3;
            app.ComponentSelectionDropDown.Layout.Column = 2;
            app.ComponentSelectionDropDown.Value = 'Tank';

            % Create ComponentLocationEditFieldLabel
            app.ComponentLocationEditFieldLabel = uilabel(app.GridLayout);
            app.ComponentLocationEditFieldLabel.HorizontalAlignment = 'center';
            app.ComponentLocationEditFieldLabel.WordWrap = 'on';
            app.ComponentLocationEditFieldLabel.Layout.Row = 4;
            app.ComponentLocationEditFieldLabel.Layout.Column = 1;
            app.ComponentLocationEditFieldLabel.Text = 'Component Location';

            % Create ComponentLocationEditField
            app.ComponentLocationEditField = uieditfield(app.GridLayout, 'numeric');
            app.ComponentLocationEditField.Layout.Row = 4;
            app.ComponentLocationEditField.Layout.Column = 2;

            % Create RocketLengthmEditFieldLabel
            app.RocketLengthmEditFieldLabel = uilabel(app.GridLayout);
            app.RocketLengthmEditFieldLabel.HorizontalAlignment = 'center';
            app.RocketLengthmEditFieldLabel.WordWrap = 'on';
            app.RocketLengthmEditFieldLabel.Layout.Row = 1;
            app.RocketLengthmEditFieldLabel.Layout.Column = 1;
            app.RocketLengthmEditFieldLabel.Text = 'Rocket Length [m]';

            % Create RocketLengthmEditField
            app.RocketLengthmEditField = uieditfield(app.GridLayout, 'numeric');
            app.RocketLengthmEditField.ValueChangedFcn = createCallbackFcn(app, @RocketLengthChanged, true);
            app.RocketLengthmEditField.HorizontalAlignment = 'left';
            app.RocketLengthmEditField.Layout.Row = 1;
            app.RocketLengthmEditField.Layout.Column = 2;

            % Create RocketDiametermEditFieldLabel
            app.RocketDiametermEditFieldLabel = uilabel(app.GridLayout);
            app.RocketDiametermEditFieldLabel.HorizontalAlignment = 'center';
            app.RocketDiametermEditFieldLabel.WordWrap = 'on';
            app.RocketDiametermEditFieldLabel.Layout.Row = 2;
            app.RocketDiametermEditFieldLabel.Layout.Column = 1;
            app.RocketDiametermEditFieldLabel.Text = 'Rocket Diameter [m]';

            % Create RocketDiametermEditField
            app.RocketDiametermEditField = uieditfield(app.GridLayout, 'numeric');
            app.RocketDiametermEditField.HorizontalAlignment = 'left';
            app.RocketDiametermEditField.Layout.Row = 2;
            app.RocketDiametermEditField.Layout.Column = 2;

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