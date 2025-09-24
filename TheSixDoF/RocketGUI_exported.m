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
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: RocketLengthmEditField
        function RocketLengthChanged(app, event)
            value = app.RocketLengthmEditField.Value;

            
            
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