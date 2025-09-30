function Enhanced_CSTR_Fluid_Simulation()
    % Enhanced CSTR Saponification Reactor with Continuous Feed System
    % Features: Two inlet pipes at top, one outlet pipe at bottom
    % Continuous particle flow through inlet/outlet system
    % Advanced Controls: Density control, viscosity control, complete dark/bright mode
    % Clean Design: Multiple Rushton turbines, motor, baffle support ring, cooling coil
    % Reaction: CH3COOC2H5 + NaOH → CH3COONa + C2H5OH
    
    fprintf('=== ENHANCED CSTR WITH CONTINUOUS FEED SYSTEM ===\n');
    fprintf('✓ Continuous Feed System: 2 Inlet + 1 Outlet Pipes\n');
    fprintf('✓ Comprehensive Dark/Bright Mode for ALL UI Elements\n');
    fprintf('✓ Pilot Plant Scale (22,000L capacity)\n');
    fprintf('✓ Triple 6-Blade Rushton Turbine Impellers (0-600 RPM)\n');
    fprintf('✓ Enhanced Fluid Physics & Wave Animation\n');
    fprintf('✓ Chemical Species Tracking & Saponification Reaction\n');
    fprintf('✓ Real-time Particle Feed/Exit System\n');
    fprintf('\nInitializing Enhanced Pilot Plant Reactor...\n');
    
    % Initialize enhanced global variables
    initializeEnhancedGlobals();
    
    % Create enhanced UI with improved dark mode support
    [fig, ax] = createAdvancedUIWithImprovedDarkMode();
    
    % Store axes handle globally
    global ax_handle;
    ax_handle = ax;
    
    % Force immediate drawing
    drawnow;
    
    fprintf('Advanced UI with continuous feed system created successfully\n');
    fprintf('Use the Dark Mode toggle to switch between themes!\n');
    fprintf('Ready for simulation - Turn on MAIN POWER and press START!\n');
end

% Enhanced global variables initialization
function initializeEnhancedGlobals()
    % Original reactor control variables
    global isPaused isStopped isRestart impeller_speed simulation_time;
    global isAnimating speedValue isDescending isRunning;
    global turbulence_level fluid_viscosity;
    global rpm_value current_time ax_handle;
    global tank_diameter tank_height tank_volume;
    global gravity_effect;
    
    % Reactor data structure
    global reactor_data;
    
    % Control components
    global knobControl rpmGauge statusLamp powerSwitch;
    global startButton pauseButton stopButton resetButton clearParticlesButton;
    global rpmValueLabel;
    
    % Advanced control variables
    global flowToggle particleToggle;
    global densitySlider viscositySlider darkModeToggle;
    global densityLabel viscosityLabel modeLabel;
    global show_flow_lines show_particles;
    global particle_density fluid_viscosity_factor dark_mode;
    global motor_graphics baffle_ring_graphics;
    
    % Feed system controls
    global feedRateSlider feedRateLabel feedToggle;
    global feed_rate_factor continuous_feeding;
    
    % UI Components for comprehensive dark mode
    global leftPanel centerPanel powerPanel rpmPanel buttonPanel;
    global vizPanel advPanel modePanel feedPanel mainFig;
    global allUIComponents;
    
    % Graphics and simulation handles
    global frame_count last_fps_time fps_display;
    global boundary_layer_thickness wall_shear_rate wall_friction_factor;
    global wall_heat_transfer_coeff reynolds_number;
    global performance_metrics simTimer;
    
    % Feed system variables
    global inlet_pipes outlet_pipe pipe_graphics;
    global feed_particle_queue exit_particle_queue;
    global total_particles_fed total_particles_exited;
    
    % Residence time controls
global residenceTimeSlider residenceTimeLabel residence_time;
residence_time = 10;  % Default 120 seconds
    
    % === PILOT PLANT REACTOR SPECIFICATIONS ===
    tank_volume = 22000;      % 22,000L pilot plant
    tank_diameter = 3.0;      % 3m diameter  
    tank_height = 4.0;        % 4m height
    
    % Simulation control flags
    isPaused = false;
    isStopped = false;
    isRestart = false;
    isAnimating = false;
    isDescending = false;
    isRunning = false;
    
    % Simulation parameters
    impeller_speed = 150.0;
    rpm_value = 300;
    simulation_time = 300;
    speedValue = 500;
    current_time = 0;
    ax_handle = [];
    
    % Enhanced fluid properties
    turbulence_level = 0.8;
    fluid_viscosity = 0.3;
    gravity_effect = 0.005;
    
    % Advanced control states
    show_flow_lines = false;
    show_particles = true;
    particle_density = 1.0;
    fluid_viscosity_factor = 1.0;
    dark_mode = false;  % Start in bright mode
    
    % Feed system parameters
    feed_rate_factor = 1.0;
    continuous_feeding = true;
    
    % Initialize feed system variables
    inlet_pipes = [];
    outlet_pipe = [];
    pipe_graphics = struct();
    feed_particle_queue = [];
    exit_particle_queue = [];
    total_particles_fed = 0;
    total_particles_exited = 0;
    
    % Timer and simulation state
    simTimer = [];
    
    % Performance monitoring
    frame_count = 0;
    last_fps_time = 0;
    fps_display = 0;
    
    % Graphics handles
    motor_graphics = [];
    baffle_ring_graphics = [];
    
    % Wall effects parameters
    boundary_layer_thickness = 0.002;
    wall_shear_rate = 50;
    wall_friction_factor = 0.02;
    wall_heat_transfer_coeff = 500;
    reynolds_number = 1000;
    
    % Add these lines in the global variables section
global targetConversionSlider targetConversionLabel target_conversion;
target_conversion = 80;  % Default 80% conversion
    
    performance_metrics = struct('fps', 0, 'conversion', 0, 'mixing_efficiency', 0);
    
    % Initialize UI components storage
    allUIComponents = struct();
    
    % Initialize reactor data
    reactor_data = initialize_pilot_plant_reactor_parameters();
    
    fprintf('Enhanced pilot plant reactor with feed system initialized successfully\n');
end

% Initialize pilot plant reactor parameters
function reactor = initialize_pilot_plant_reactor_parameters()
    global tank_diameter tank_height particle_density fluid_viscosity_factor;
    
    reactor = struct();
    
    % Pilot Plant Reactor Geometry
    reactor.tank_radius = tank_diameter / 2;
    reactor.tank_height = tank_height;
    reactor.fluid_level = tank_height * 0.8;
    
    % Multiple High-Speed Impeller Parameters
    reactor.impeller_count = 3;
    reactor.impeller_speed = 150.0;
    reactor.max_rpm = 600.0;
    reactor.impeller_diameter = 1.0;
    
    % Impeller positions
    reactor.impeller_heights = [0.8, 2.0, 3.2];
    reactor.blade_count = 6;
    reactor.impeller_running = true;
    
    % Feed System Parameters
    reactor.inlet_pipe_count = 2;
    reactor.outlet_pipe_count = 1;
    reactor.feed_rate = 10.0;  % particles per second
    reactor.exit_rate = 8.0;   % particles per second
    reactor.residence_time = 120.0;  % seconds
    
    % Inlet positions (at top of tank)
    reactor.inlet_positions = [
        [reactor.tank_radius * 0.7, 0, reactor.fluid_level * 0.95];  % Inlet 1: EtOAc
        [-reactor.tank_radius * 0.7, 0, reactor.fluid_level * 0.95]  % Inlet 2: NaOH
    ];
    
    % Outlet position (at bottom of tank)
    reactor.outlet_position = [0, reactor.tank_radius * 0.8, 0.1];
    
    % Fluid Properties
    base_particle_count = 1000;
    reactor.particle_count = round(base_particle_count * particle_density);
    reactor.density = 950;
    reactor.viscosity = 1.2 * fluid_viscosity_factor;
    reactor.surface_tension = 0.025;
    
    % Chemical Reaction Parameters
    reactor.ethyl_acetate_conc = 1.0;
    reactor.naoh_conc = 1.2;
    reactor.conversion = 0;
    reactor.reaction_rate = 0.05;
    
    % Enhanced Physics Parameters
    reactor.dt = 0.015;
    reactor.gravity = 9.81;
    reactor.mixing_efficiency = 0;
    reactor.reynolds_number = 0;
    reactor.power_number = 6.3;
    
    % Visualization parameters
    reactor.show_flow_lines = true;
    reactor.show_particles = true;
    reactor.show_reaction_zones = true;
    reactor.show_waves = true;
    reactor.show_feed_system = true;
    reactor.camera_angle = [35, 20];
    
    % Animation state
    reactor.time = 0;
    reactor.frame_count = 0;
    reactor.impeller_angles = zeros(1, reactor.impeller_count);
    reactor.start_time = tic;
    
    % Feed system tracking
    reactor.particles_in_system = 0;
    reactor.total_fed = 0;
    reactor.total_exited = 0;
    reactor.current_holdup = 0;
    
    % Initialize chemical species storage
    reactor.species = struct();
    reactor.species.ethyl_acetate = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
    reactor.species.naoh = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
    reactor.species.products = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
    
    % Initialize blade data storage
    reactor.blade_original_vertices = cell(reactor.impeller_count, reactor.blade_count);
    for imp = 1:reactor.impeller_count
        for blade = 1:reactor.blade_count
            reactor.blade_original_vertices{imp, blade} = [];
        end
    end
    
    % Graphics handles initialization
    reactor.fig = [];
    reactor.ax = [];
    reactor.tank_walls = [];
    reactor.tank_bottom = [];
    reactor.fluid_surface = [];
    reactor.impeller_shafts = [];
    reactor.impeller_blades = cell(reactor.impeller_count, 1);
    reactor.impeller_hubs = [];
    reactor.motor = [];
    reactor.species_plots = struct();
    reactor.flow_lines = [];
    reactor.feed_pipes = [];
    reactor.outlet_pipes = [];
    
    fprintf('Pilot plant reactor with feed system initialized: %.1fm tank, %d impellers, %d inlets, %d outlet\n', ...
            tank_diameter, reactor.impeller_count, reactor.inlet_pipe_count, reactor.outlet_pipe_count);
end

% Get current color scheme
function colors = getCurrentColorScheme()
    global dark_mode;
    
    if dark_mode
        % Complete Dark color scheme - ALL elements dark
        colors.bg_color = [0.08, 0.08, 0.12];           % Very dark main background
        colors.panel_color = [0.12, 0.12, 0.16];        % Dark panels
        colors.text_color = [0.9, 0.9, 0.95];           % Light text
        colors.center_panel_color = [0.08, 0.08, 0.12]; % Dark center panel (same as main)
        colors.center_text_color = [0.8, 0.85, 0.9];    % Light center text
        colors.button_bg = [0.2, 0.2, 0.25];            % Dark buttons
        colors.button_text = [0.9, 0.9, 0.95];          % Light button text
        colors.accent_color = [0.4, 0.7, 1.0];          % Bright blue accent
        colors.warning_color = [1.0, 0.7, 0.3];         % Bright orange accent
        colors.success_color = [0.4, 0.9, 0.5];         % Bright green accent
        colors.error_color = [1.0, 0.4, 0.4];           % Bright red accent
        colors.grid_color = [0.3, 0.3, 0.35];           % Dark grid
        colors.axes_bg = [0.05, 0.05, 0.08];            % Very dark axes background
    else
        % Bright color scheme
        colors.bg_color = [0.95, 0.95, 0.98];           % Light background
        colors.panel_color = [0.88, 0.92, 0.98];        % Light panels
        colors.text_color = [0.1, 0.2, 0.5];            % Dark text
        colors.center_panel_color = [0.98, 0.98, 1.0];  % Very light center panel
        colors.center_text_color = [0.2, 0.1, 0.4];     % Dark center text
        colors.button_bg = [0.9, 0.9, 0.95];            % Light buttons
        colors.button_text = [0.2, 0.2, 0.4];           % Dark button text
        colors.accent_color = [0.2, 0.4, 0.8];          % Blue accent
        colors.warning_color = [0.8, 0.5, 0.1];         % Orange accent
        colors.success_color = [0.2, 0.6, 0.3];         % Green accent
        colors.error_color = [0.8, 0.2, 0.2];           % Red accent
        colors.grid_color = [0.3, 0.4, 0.6];            % Light grid
        colors.axes_bg = [0.96, 0.98, 1.0];             % Very light axes background
    end
end

% Create advanced UI with improved dark mode support
function [fig, ax] = createAdvancedUIWithImprovedDarkMode()
    global mainFig allUIComponents;
    
    % Get current color scheme
    colors = getCurrentColorScheme();
    
    % Create main figure
    fig = uifigure('Position', [30, 30, 1400, 1100], ...
                   'Name', 'Enhanced CSTR - Continuous Feed System', ...
                   'Color', colors.bg_color, ...
                   'AutoResizeChildren', 'off');
    
    mainFig = fig;
    allUIComponents.mainFig = fig;
    
    % Create main grid layout
    mainGrid = uigridlayout(fig, [1, 2]);
    mainGrid.ColumnWidth = {350, '1x'};
    mainGrid.RowHeight = {'1x'};
    mainGrid.Padding = [15, 15, 15, 15];
    mainGrid.ColumnSpacing = 15;
    mainGrid.BackgroundColor = colors.bg_color;
    allUIComponents.mainGrid = mainGrid;
    
    % === LEFT PANEL ===
    global leftPanel;
    leftPanel = uipanel(mainGrid, 'Title', 'PILOT PLANT REACTOR CONTROL SYSTEM', ...
                       'FontWeight', 'bold', 'FontSize', 14, ...
                       'BackgroundColor', colors.panel_color, ...
                       'ForegroundColor', colors.text_color);
    leftPanel.Layout.Row = 1;
    leftPanel.Layout.Column = 1;
    allUIComponents.leftPanel = leftPanel;
    
    % === CENTER PANEL ===
    global centerPanel;
    centerPanel = uipanel(mainGrid, 'Title', 'CONTINUOUS FEED SAPONIFICATION REACTOR', ...
                         'FontWeight', 'bold', 'FontSize', 16, ...
                         'BackgroundColor', colors.center_panel_color, ...
                         'ForegroundColor', colors.center_text_color);
    centerPanel.Layout.Row = 1;
    centerPanel.Layout.Column = 2;
    allUIComponents.centerPanel = centerPanel;
    
    % Create center grid
    centerGrid = uigridlayout(centerPanel, [1, 1]);
    centerGrid.RowHeight = {'1x'};
    centerGrid.ColumnWidth = {'1x'};
    centerGrid.Padding = [10, 10, 10, 10];
    centerGrid.BackgroundColor = colors.center_panel_color;
    allUIComponents.centerGrid = centerGrid;
    
    % 3D visualization axes
    ax = uiaxes(centerGrid);
    ax.Layout.Row = 1;
    ax.Layout.Column = 1;
    allUIComponents.mainAxes = ax;
    
    % Setup axes properties
    setupAxesProperties(ax);
    
    % Create controls
    createAdvancedControlsWithDarkMode(leftPanel);
    
    % Draw reactor geometry
    drawEnhancedReactorGeometry(ax);
    
    % Set close callback
    fig.CloseRequestFcn = @(src, event) figureCloseCallback(src, event);
    
    fprintf('Enhanced UI with continuous feed system created successfully\n');
end

% Setup axes properties with complete dark mode support
function setupAxesProperties(ax)
    global tank_diameter tank_height dark_mode;
    
    % Get current colors
    colors = getCurrentColorScheme();
    
    hold(ax, 'on');
    axis(ax, 'equal');
    CSTR_radius = tank_diameter / 2;
    xlim(ax, [-CSTR_radius*1.3 CSTR_radius*1.3]);
    ylim(ax, [-CSTR_radius*1.3 CSTR_radius*1.3]);
    zlim(ax, [-0.2 tank_height*1.2]);
    view(ax, 35, 25);
    grid(ax, 'on');
    
    if dark_mode
        % Complete dark theme for axes
        ax.GridAlpha = 0.2;
        ax.GridColor = colors.grid_color;
        ax.Color = colors.axes_bg;               % Very dark axes background
        ax.XColor = [0.7, 0.8, 0.9];            % Light axis colors
        ax.YColor = [0.7, 0.8, 0.9];
        ax.ZColor = [0.7, 0.8, 0.9];
        title_color = [0.8, 0.9, 1.0];          % Bright title
        
        % Set parent figure color to ensure complete dark background
        fig = ancestor(ax, 'figure');
        if ~isempty(fig)
            fig.Color = colors.bg_color;
        end
    else
        % Bright theme for axes
        ax.GridAlpha = 0.6;
        ax.GridColor = colors.grid_color;
        ax.Color = colors.axes_bg;               % Very light axes background
        ax.XColor = [0.2, 0.3, 0.6];
        ax.YColor = [0.2, 0.3, 0.6];
        ax.ZColor = [0.2, 0.3, 0.6];
        title_color = [0.1, 0.3, 0.7];
        
        % Set parent figure color
        fig = ancestor(ax, 'figure');
        if ~isempty(fig)
            fig.Color = colors.bg_color;
        end
    end
    
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    xlabel(ax, 'X (m)', 'Color', ax.XColor, 'FontWeight', 'bold');
    ylabel(ax, 'Y (m)', 'Color', ax.YColor, 'FontWeight', 'bold');
    zlabel(ax, 'Z (m)', 'Color', ax.ZColor, 'FontWeight', 'bold');
    
    % Enable 3D navigation
    rotate3d(ax, 'on');
    zoom(ax, 'on');
    pan(ax, 'on');
    
    % Add enhanced lighting for better visualization
    lighting(ax, 'gouraud');
    if dark_mode
        light(ax, 'Position', [1, 1, 2], 'Color', [0.8, 0.8, 0.9]);
        light(ax, 'Position', [-1, -1, 0.5], 'Color', [0.6, 0.7, 0.8]);
        light(ax, 'Position', [0, 1, -1], 'Color', [0.5, 0.6, 0.7]);
    else
        light(ax, 'Position', [1, 1, 2], 'Color', [1, 1, 1]);
        light(ax, 'Position', [-1, -1, 0.5], 'Color', [0.8, 0.8, 1]);
        light(ax, 'Position', [0, 1, -1], 'Color', [0.6, 0.8, 0.8]);
    end
end

% Create advanced controls with comprehensive dark mode support
function createAdvancedControlsWithDarkMode(leftPanel)
    % Create grid layout with immediate color application
   leftGrid = uigridlayout(leftPanel, [6, 1]);
leftGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};  % All 'fit' to prevent overflow
    leftGrid.ColumnWidth = {'1x'};
    leftGrid.Padding = [10, 10, 10, 10];
    leftGrid.RowSpacing = 8;
    
    % Apply current colors immediately
    colors = getCurrentColorScheme();
    leftGrid.BackgroundColor = colors.panel_color;
    
    % Store grid for theme updates
    global allUIComponents;
    allUIComponents.leftGrid = leftGrid;
    
    % Create control sections in order
    createPowerAndModeControls(leftGrid);     % Row 1: Combined Power + Mode
    createRPMControls(leftGrid);              % Row 2: RPM Controls  
    createFeedSystemControls(leftGrid);       % Row 3: Feed System
    createVisualizationControlsWithDarkMode(leftGrid); % Row 4: Visualization
    createResidenceTimeControls(leftGrid); 
    createAdvancedSlidersWithDarkMode(leftGrid);       % Row 5: Advanced
    createOperationControls(leftGrid);        % Row 6: Operation Buttons
end

% Create combined power and mode controls side by side
function createPowerAndModeControls(parent)
    global powerSwitch statusLamp darkModeToggle modeLabel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Combined Panel for Power and Mode
    combinedPanel = uipanel(parent, 'Title', 'SYSTEM POWER & DISPLAY MODE', ...
                           'BackgroundColor', colors.panel_color, ...
                           'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    combinedPanel.Layout.Row = 1;
    allUIComponents.combinedPanel = combinedPanel;
    
    % Grid with 2 columns for side-by-side layout
    combinedGrid = uigridlayout(combinedPanel, [2, 2]);
    combinedGrid.ColumnWidth = {'1x', '1x'};
    combinedGrid.RowHeight = {'fit', 'fit'};
    combinedGrid.Padding = [10, 10, 10, 10];
    combinedGrid.RowSpacing = 8;
    combinedGrid.ColumnSpacing = 15;
    combinedGrid.BackgroundColor = colors.panel_color;
    allUIComponents.combinedGrid = combinedGrid;
    
    % === LEFT SIDE: POWER CONTROLS ===
    powerSubGrid = uigridlayout(combinedGrid, [2, 2]);
    powerSubGrid.Layout.Row = [1, 2];
    powerSubGrid.Layout.Column = 1;
    powerSubGrid.ColumnWidth = {'fit', '1x'};
    powerSubGrid.RowHeight = {'fit', 'fit'};
    powerSubGrid.Padding = [5, 5, 5, 5];
    powerSubGrid.RowSpacing = 8;
    powerSubGrid.BackgroundColor = colors.panel_color;
    allUIComponents.powerSubGrid = powerSubGrid;
    
    % Power switch
    powerSwitch = uiswitch(powerSubGrid, 'toggle');
    powerSwitch.Layout.Row = 1;
    powerSwitch.Layout.Column = 1;
    powerSwitch.ValueChangedFcn = @switchValueChanged;
    allUIComponents.powerSwitch = powerSwitch;
    
    powerLabel = uilabel(powerSubGrid, 'Text', 'MAIN POWER', ...
                        'FontWeight', 'bold', 'FontSize', 11, ...
                        'FontColor', colors.text_color);
    powerLabel.Layout.Row = 1;
    powerLabel.Layout.Column = 2;
    allUIComponents.powerLabel = powerLabel;
    
    % Status lamp
    statusLamp = uilamp(powerSubGrid, 'Color', colors.error_color);
    statusLamp.Layout.Row = 2;
    statusLamp.Layout.Column = 1;
    allUIComponents.statusLamp = statusLamp;
    
    statusLabel = uilabel(powerSubGrid, 'Text', 'SYSTEM STATUS', ...
                         'FontWeight', 'bold', 'FontSize', 11, ...
                         'FontColor', colors.text_color);
    statusLabel.Layout.Row = 2;
    statusLabel.Layout.Column = 2;
    allUIComponents.statusLabel = statusLabel;
    
    % === RIGHT SIDE: DISPLAY MODE CONTROLS ===
    modeSubGrid = uigridlayout(combinedGrid, [2, 3]);
    modeSubGrid.Layout.Row = [1, 2];
    modeSubGrid.Layout.Column = 2;
    modeSubGrid.ColumnWidth = {'1x', 'fit', 'fit'};
    modeSubGrid.RowHeight = {'fit', 'fit'};
    modeSubGrid.Padding = [5, 5, 5, 5];
    modeSubGrid.RowSpacing = 8;
    modeSubGrid.BackgroundColor = colors.panel_color;
    allUIComponents.modeSubGrid = modeSubGrid;
    
    % Current mode label
    modeLabel = uilabel(modeSubGrid, 'Text', 'Bright Mode', ...
                       'FontWeight', 'bold', 'FontSize', 11, ...
                       'FontColor', colors.text_color);
    modeLabel.Layout.Row = 1;
    modeLabel.Layout.Column = [1, 3];
    allUIComponents.modeLabel = modeLabel;
    
    % Mode switch label
    switchLabel = uilabel(modeSubGrid, 'Text', 'Dark:', ...
                         'FontWeight', 'bold', 'FontSize', 11, ...
                         'FontColor', colors.text_color, ...
                         'HorizontalAlignment', 'right');
    switchLabel.Layout.Row = 2;
    switchLabel.Layout.Column = 2;
    allUIComponents.switchLabel = switchLabel;
    
    % Dark/Bright mode toggle
    darkModeToggle = uiswitch(modeSubGrid, 'toggle', 'Value', 'Off', ...
                             'ValueChangedFcn', @enhancedDarkModeToggleCallback);
    darkModeToggle.Layout.Row = 2;
    darkModeToggle.Layout.Column = 3;
    allUIComponents.darkModeToggle = darkModeToggle;
end

% Create RPM controls only (no power controls)
function createRPMControls(parent)
    global knobControl rpmGauge rpmValueLabel rpmPanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % RPM Control Section
    rpmPanel = uipanel(parent, 'Title', 'TRIPLE RUSHTON TURBINE CONTROL (0-600 RPM)', ...
                      'BackgroundColor', colors.panel_color, ...
                      'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    rpmPanel.Layout.Row = 2;
    allUIComponents.rpmPanel = rpmPanel;
    
    rpmGrid = uigridlayout(rpmPanel, [3, 2]);
    rpmGrid.ColumnWidth = {'1x', '1x'};
    rpmGrid.RowHeight = {100, 'fit', 'fit'};
    rpmGrid.Padding = [15, 15, 15, 15];
    rpmGrid.RowSpacing = 8;
    rpmGrid.BackgroundColor = colors.panel_color;
    allUIComponents.rpmGrid = rpmGrid;
    
    % RPM Gauge
    rpmGauge = uigauge(rpmGrid, 'circular', 'Limits', [0, 600], ...
                      'MajorTicks', [0, 100, 200, 300, 400, 500, 600]);
    rpmGauge.Layout.Row = 1;
    rpmGauge.Layout.Column = 1;
    rpmGauge.BackgroundColor = colors.button_bg;
    rpmGauge.ScaleColors = colors.accent_color;
    allUIComponents.rpmGauge = rpmGauge;
    
    % RPM Control Knob
    knobControl = uiknob(rpmGrid, 'continuous', 'Limits', [0, 600], 'Value', 300);
    knobControl.Layout.Row = 1;
    knobControl.Layout.Column = 2;
    knobControl.ValueChangedFcn = @knobValueChanged;
    allUIComponents.knobControl = knobControl;
    
    % Labels
    gaugeLabel = uilabel(rpmGrid, 'Text', 'RPM METER', ...
                        'FontSize', 10, 'HorizontalAlignment', 'center', ...
                        'FontColor', colors.text_color, 'FontWeight', 'bold');
    gaugeLabel.Layout.Row = 2;
    gaugeLabel.Layout.Column = 1;
    allUIComponents.gaugeLabel = gaugeLabel;
    
    knobLabel = uilabel(rpmGrid, 'Text', 'RPM CONTROL', ...
                       'FontSize', 10, 'HorizontalAlignment', 'center', ...
                       'FontColor', colors.text_color, 'FontWeight', 'bold');
    knobLabel.Layout.Row = 2;
    knobLabel.Layout.Column = 2;
    allUIComponents.knobLabel = knobLabel;
    
    % RPM Display
    rpmValueLabel = uilabel(rpmGrid, 'Text', '300 RPM', ...
                           'FontWeight', 'bold', 'FontSize', 14, ...
                           'HorizontalAlignment', 'center', ...
                           'FontColor', colors.accent_color);
    rpmValueLabel.Layout.Row = 3;
    rpmValueLabel.Layout.Column = [1, 2];
    allUIComponents.rpmValueLabel = rpmValueLabel;
end

% Create separate operation controls panel
% Create separate operation controls panel
function createOperationControls(parent)
    global startButton pauseButton stopButton resetButton clearParticlesButton;
    global buttonPanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Operation Buttons Panel
    buttonPanel = uipanel(parent, 'Title', 'OPERATION CONTROL', ...
                         'BackgroundColor', colors.panel_color, ...
                         'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    buttonPanel.Layout.Row = 6;
    allUIComponents.buttonPanel = buttonPanel;
    
    buttonGrid = uigridlayout(buttonPanel, [3, 3]);
    buttonGrid.ColumnWidth = {'1x', '1x', '1x'};
    buttonGrid.RowHeight = {35, 35, 35};
    buttonGrid.Padding = [8, 8, 8, 8];
    buttonGrid.RowSpacing = 6;
    buttonGrid.ColumnSpacing = 6;
    buttonGrid.BackgroundColor = colors.panel_color;
    allUIComponents.buttonGrid = buttonGrid;
    
    % Control buttons
    startButton = uibutton(buttonGrid, 'push', 'Text', 'START', ...
        'BackgroundColor', colors.success_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @startFluidSimulation);
    startButton.Layout.Row = 1;
    startButton.Layout.Column = 1;
    resetButtonState(startButton);  % ADDED: Fix button state
    allUIComponents.startButton = startButton;
    
    pauseButton = uibutton(buttonGrid, 'push', 'Text', 'PAUSE', ...
        'BackgroundColor', colors.warning_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @pauseButtonCallback);
    pauseButton.Layout.Row = 1;
    pauseButton.Layout.Column = 2;
    resetButtonState(pauseButton);  % ADDED: Fix button state
    allUIComponents.pauseButton = pauseButton;
    
    stopButton = uibutton(buttonGrid, 'push', 'Text', 'STOP', ...
        'BackgroundColor', colors.error_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @stopButtonCallback);
    stopButton.Layout.Row = 1;
    stopButton.Layout.Column = 3;
    resetButtonState(stopButton);  % ADDED: Fix button state
    allUIComponents.stopButton = stopButton;
    
    resetButton = uibutton(buttonGrid, 'push', 'Text', 'RESET', ...
        'BackgroundColor', colors.accent_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @resetButtonCallback);
    resetButton.Layout.Row = 2;
    resetButton.Layout.Column = 1;
    resetButtonState(resetButton);  % ADDED: Fix button state
    allUIComponents.resetButton = resetButton;
    
    % Start/Stop Feed button
    startFeedButton = uibutton(buttonGrid, 'push', 'Text', 'START FEED', ...
        'BackgroundColor', colors.accent_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 10, ...
        'ButtonPushedFcn', @startFeedCallback);
    startFeedButton.Layout.Row = 2;
    startFeedButton.Layout.Column = 2;
    resetButtonState(startFeedButton);  % ADDED: Fix button state
    allUIComponents.startFeedButton = startFeedButton;
    
    % Emergency stop
    emergencyButton = uibutton(buttonGrid, 'push', 'Text', 'E-STOP', ...
        'BackgroundColor', colors.error_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @emergencyStopCallback);
    emergencyButton.Layout.Row = 2;
    emergencyButton.Layout.Column = 3;
    resetButtonState(emergencyButton);  % ADDED: Fix button state
    allUIComponents.emergencyButton = emergencyButton;
    
    % Clear All Particles button
    clearParticlesButton = uibutton(buttonGrid, 'push', 'Text', 'CLEAR REACTOR', ...
        'BackgroundColor', colors.warning_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 10, ...
        'ButtonPushedFcn', @clearAllParticlesCallback);
    clearParticlesButton.Layout.Row = 3;
    clearParticlesButton.Layout.Column = [1, 3];
    resetButtonState(clearParticlesButton);  % ADDED: Fix button state
    allUIComponents.clearParticlesButton = clearParticlesButton;
    
    % Initialize disabled state
    disableAllControls();
end


% Create feed system controls
function createFeedSystemControls(parent)
    global feedRateSlider feedRateLabel feedToggle feedPanel allUIComponents;
    global feed_rate_factor continuous_feeding;
    
    colors = getCurrentColorScheme();
    
    % Feed System Control Panel
    feedPanel = uipanel(parent, 'Title', 'CONTINUOUS FEED SYSTEM', ...
                       'BackgroundColor', colors.panel_color, ...
                       'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    feedPanel.Layout.Row = 3;
    allUIComponents.feedPanel = feedPanel;
    
    feedGrid = uigridlayout(feedPanel, [4, 2]);
    feedGrid.ColumnWidth = {'fit', '1x'};
    feedGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    feedGrid.Padding = [10, 10, 10, 10];
    feedGrid.RowSpacing = 8;
    feedGrid.BackgroundColor = colors.panel_color;
    allUIComponents.feedGrid = feedGrid;
    
    % Feed toggle
    feedLabel = uilabel(feedGrid, 'Text', 'Continuous Feed:', ...
                       'FontWeight', 'bold', 'FontSize', 10, ...
                       'FontColor', colors.text_color);
    feedLabel.Layout.Row = 1;
    feedLabel.Layout.Column = 1;
    allUIComponents.feedLabel = feedLabel;
    
    feedToggle = uicheckbox(feedGrid, 'Text', 'Enable', 'Value', true, ...
                           'FontColor', colors.text_color, ...
                           'ValueChangedFcn', @feedToggleCallback);
    feedToggle.Layout.Row = 1;
    feedToggle.Layout.Column = 2;
    allUIComponents.feedToggle = feedToggle;
    
    % Feed rate slider
    feedRateTitleLabel = uilabel(feedGrid, 'Text', 'Feed Rate:', ...
                                'FontWeight', 'bold', 'FontSize', 10, ...
                                'FontColor', colors.text_color);
    feedRateTitleLabel.Layout.Row = 2;
    feedRateTitleLabel.Layout.Column = 1;
    allUIComponents.feedRateTitleLabel = feedRateTitleLabel;
    
    feedRateSlider = uislider(feedGrid, 'Limits', [0.1, 3.0], 'Value', 1.0, ...
                             'MajorTicks', [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], ...
                             'ValueChangedFcn', @feedRateSliderCallback);
    feedRateSlider.Layout.Row = 3;
    feedRateSlider.Layout.Column = [1, 2];
    allUIComponents.feedRateSlider = feedRateSlider;
    
    feedRateLabel = uilabel(feedGrid, 'Text', '1.00x', ...
                           'FontWeight', 'bold', 'FontSize', 10, ...
                           'FontColor', colors.accent_color, ...
                           'HorizontalAlignment', 'center');
    feedRateLabel.Layout.Row = 2;
    feedRateLabel.Layout.Column = 2;
    allUIComponents.feedRateLabel = feedRateLabel;
    
    % Feed status
    feedStatusLabel = uilabel(feedGrid, 'Text', 'Status: Ready for continuous operation', ...
                             'FontSize', 9, 'FontColor', colors.success_color, ...
                             'HorizontalAlignment', 'left');
    feedStatusLabel.Layout.Row = 4;
    feedStatusLabel.Layout.Column = [1, 2];
    allUIComponents.feedStatusLabel = feedStatusLabel;
    
    % Initialize feed system
    continuous_feeding = true;
    feed_rate_factor = 1.0;
end

% Create residence time controls
function createResidenceTimeControls(parent)
    global residenceTimeSlider residenceTimeLabel residence_time allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Residence Time Control Panel
    residencePanel = uipanel(parent, 'Title', 'RESIDENCE TIME CONTROL', ...
                            'BackgroundColor', colors.panel_color, ...
                            'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    residencePanel.Layout.Row = 4;
    allUIComponents.residencePanel = residencePanel;
    
    residenceGrid = uigridlayout(residencePanel, [3, 2]);
    residenceGrid.ColumnWidth = {'fit', '1x'};
    residenceGrid.RowHeight = {'fit', 'fit', 'fit'};
    residenceGrid.Padding = [10, 10, 10, 10];
    residenceGrid.RowSpacing = 8;
    residenceGrid.BackgroundColor = colors.panel_color;
    allUIComponents.residenceGrid = residenceGrid;
    
    % Residence time title
    residenceTitleLabel = uilabel(residenceGrid, 'Text', 'Reaction Time:', ...
                                 'FontWeight', 'bold', 'FontSize', 10, ...
                                 'FontColor', colors.text_color);
    residenceTitleLabel.Layout.Row = 1;
    residenceTitleLabel.Layout.Column = 1;
    allUIComponents.residenceTitleLabel = residenceTitleLabel;
    
    % Residence time value display
    residenceTimeLabel = uilabel(residenceGrid, 'Text', '10s', ...  % CHANGED: Default to 10s
                                'FontWeight', 'bold', 'FontSize', 10, ...
                                'FontColor', colors.accent_color, ...
                                'HorizontalAlignment', 'center');
    residenceTimeLabel.Layout.Row = 1;
    residenceTimeLabel.Layout.Column = 2;
    allUIComponents.residenceTimeLabel = residenceTimeLabel;
    
    % Residence time slider - CHANGED: Range from 5-120s, default 10s
    residenceTimeSlider = uislider(residenceGrid, 'Limits', [5, 120], 'Value', 10, ...
                                  'MajorTicks', [5, 10, 20, 30, 60, 120], ...
                                  'ValueChangedFcn', @residenceTimeSliderCallback);
    residenceTimeSlider.Layout.Row = 2;
    residenceTimeSlider.Layout.Column = [1, 2];
    allUIComponents.residenceTimeSlider = residenceTimeSlider;
    
    % Info label
    infoLabel = uilabel(residenceGrid, 'Text', 'Time for particles to turn GREEN (5-120s)', ...
                       'FontSize', 8, 'FontColor', colors.success_color, ...
                       'HorizontalAlignment', 'center');
    infoLabel.Layout.Row = 3;
    infoLabel.Layout.Column = [1, 2];
    allUIComponents.residenceInfoLabel = infoLabel;
    
    % Initialize residence time - CHANGED: Default to 10 seconds
    residence_time = 10;  % Much faster conversion!
end

% Create main controls with dark mode support
function createMainControlsWithDarkMode(parent)
    global knobControl rpmGauge statusLamp powerSwitch;
    global startButton pauseButton stopButton resetButton clearParticlesButton;
    global rpmValueLabel powerPanel rpmPanel buttonPanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Power & Status Section
    powerPanel = uipanel(parent, 'Title', 'SYSTEM POWER', ...
                        'BackgroundColor', colors.panel_color, ...
                        'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    powerPanel.Layout.Row = 2;  % Changed from 1 to 2
    allUIComponents.powerPanel = powerPanel;
    
    powerGrid = uigridlayout(powerPanel, [2, 2]);
    powerGrid.ColumnWidth = {'fit', '1x'};
    powerGrid.RowHeight = {'fit', 'fit'};
    powerGrid.Padding = [15, 15, 15, 15];
    powerGrid.RowSpacing = 10;
    powerGrid.BackgroundColor = colors.panel_color;
    
    % Power switch
    powerSwitch = uiswitch(powerGrid, 'toggle');
    powerSwitch.Layout.Row = 1;
    powerSwitch.Layout.Column = 1;
    powerSwitch.ValueChangedFcn = @switchValueChanged;
    allUIComponents.powerSwitch = powerSwitch;
    
    powerLabel = uilabel(powerGrid, 'Text', 'MAIN POWER', ...
                        'FontWeight', 'bold', 'FontSize', 12, ...
                        'FontColor', colors.text_color);
    powerLabel.Layout.Row = 1;
    powerLabel.Layout.Column = 2;
    allUIComponents.powerLabel = powerLabel;
    
    % Status lamp
    statusLamp = uilamp(powerGrid, 'Color', colors.error_color);
    statusLamp.Layout.Row = 2;
    statusLamp.Layout.Column = 1;
    allUIComponents.statusLamp = statusLamp;
    
    statusLabel = uilabel(powerGrid, 'Text', 'SYSTEM STATUS', ...
                         'FontWeight', 'bold', 'FontSize', 12, ...
                         'FontColor', colors.text_color);
    statusLabel.Layout.Row = 2;
    statusLabel.Layout.Column = 2;
    allUIComponents.statusLabel = statusLabel;
    
    % RPM Control Section
    rpmPanel = uipanel(parent, 'Title', 'TRIPLE RUSHTON TURBINE CONTROL (0-600 RPM)', ...
                      'BackgroundColor', colors.panel_color, ...
                      'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    rpmPanel.Layout.Row = 3;  % Changed from 2 to 3
    allUIComponents.rpmPanel = rpmPanel;
    
    rpmGrid = uigridlayout(rpmPanel, [3, 2]);
    rpmGrid.ColumnWidth = {'1x', '1x'};
    rpmGrid.RowHeight = {100, 'fit', 'fit'};
    rpmGrid.Padding = [15, 15, 15, 15];
    rpmGrid.RowSpacing = 8;
    rpmGrid.BackgroundColor = colors.panel_color;
    allUIComponents.rpmGrid = rpmGrid;
    
    % RPM Gauge
    rpmGauge = uigauge(rpmGrid, 'circular', 'Limits', [0, 600], ...
                      'MajorTicks', [0, 100, 200, 300, 400, 500, 600]);
    rpmGauge.Layout.Row = 1;
    rpmGauge.Layout.Column = 1;
    rpmGauge.BackgroundColor = colors.button_bg;
    rpmGauge.ScaleColors = colors.accent_color;
    allUIComponents.rpmGauge = rpmGauge;
    
    % RPM Control Knob
    knobControl = uiknob(rpmGrid, 'continuous', 'Limits', [0, 600], 'Value', 300);
    knobControl.Layout.Row = 1;
    knobControl.Layout.Column = 2;
    knobControl.ValueChangedFcn = @knobValueChanged;
    allUIComponents.knobControl = knobControl;
    
    % Labels
    gaugeLabel = uilabel(rpmGrid, 'Text', 'RPM METER', ...
                        'FontSize', 10, 'HorizontalAlignment', 'center', ...
                        'FontColor', colors.text_color, 'FontWeight', 'bold');
    gaugeLabel.Layout.Row = 2;
    gaugeLabel.Layout.Column = 1;
    allUIComponents.gaugeLabel = gaugeLabel;
    
    knobLabel = uilabel(rpmGrid, 'Text', 'RPM CONTROL', ...
                       'FontSize', 10, 'HorizontalAlignment', 'center', ...
                       'FontColor', colors.text_color, 'FontWeight', 'bold');
    knobLabel.Layout.Row = 2;
    knobLabel.Layout.Column = 2;
    allUIComponents.knobLabel = knobLabel;
    
    % RPM Display
    rpmValueLabel = uilabel(rpmGrid, 'Text', '300 RPM', ...
                           'FontWeight', 'bold', 'FontSize', 14, ...
                           'HorizontalAlignment', 'center', ...
                           'FontColor', colors.accent_color);
    rpmValueLabel.Layout.Row = 3;
    rpmValueLabel.Layout.Column = [1, 2];
    allUIComponents.rpmValueLabel = rpmValueLabel;
    
    % Operation Buttons
    buttonPanel = uipanel(parent, 'Title', 'OPERATION CONTROL', ...
                         'BackgroundColor', colors.panel_color, ...
                         'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    buttonPanel.Layout.Row = 7;  % Changed from 6 to 7
    allUIComponents.buttonPanel = buttonPanel;
    
    buttonGrid = uigridlayout(buttonPanel, [3, 3]);
    buttonGrid.ColumnWidth = {'1x', '1x', '1x'};
    buttonGrid.RowHeight = {'1x', '1x', '1x'};
    buttonGrid.Padding = [8, 8, 8, 8];
    buttonGrid.RowSpacing = 6;
    buttonGrid.ColumnSpacing = 6;
    buttonGrid.BackgroundColor = colors.panel_color;
    
    % Control buttons
    startButton = uibutton(buttonGrid, 'push', 'Text', 'START', ...
        'BackgroundColor', colors.success_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @startFluidSimulation);
    startButton.Layout.Row = 1;
    startButton.Layout.Column = 1;
    allUIComponents.startButton = startButton;
    
    pauseButton = uibutton(buttonGrid, 'push', 'Text', 'PAUSE', ...
        'BackgroundColor', colors.warning_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @pauseButtonCallback);
    pauseButton.Layout.Row = 1;
    pauseButton.Layout.Column = 2;
    allUIComponents.pauseButton = pauseButton;
    
    stopButton = uibutton(buttonGrid, 'push', 'Text', 'STOP', ...
        'BackgroundColor', colors.error_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @stopButtonCallback);
    stopButton.Layout.Row = 1;
    stopButton.Layout.Column = 3;
    allUIComponents.stopButton = stopButton;
    
    resetButton = uibutton(buttonGrid, 'push', 'Text', 'RESET', ...
        'BackgroundColor', colors.accent_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @resetButtonCallback);
    resetButton.Layout.Row = 2;
    resetButton.Layout.Column = 1;
    allUIComponents.resetButton = resetButton;
    
    % Start/Stop Feed button
    startFeedButton = uibutton(buttonGrid, 'push', 'Text', 'START FEED', ...
        'BackgroundColor', colors.accent_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 10, ...
        'ButtonPushedFcn', @startFeedCallback);
    startFeedButton.Layout.Row = 2;
    startFeedButton.Layout.Column = 2;
    allUIComponents.startFeedButton = startFeedButton;
    
    % Emergency stop
    emergencyButton = uibutton(buttonGrid, 'push', 'Text', 'E-STOP', ...
        'BackgroundColor', colors.error_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 11, ...
        'ButtonPushedFcn', @emergencyStopCallback);
    emergencyButton.Layout.Row = 2;
    emergencyButton.Layout.Column = 3;
    allUIComponents.emergencyButton = emergencyButton;
    
    % Clear All Particles button
    clearParticlesButton = uibutton(buttonGrid, 'push', 'Text', 'CLEAR REACTOR', ...
        'BackgroundColor', colors.warning_color, 'FontWeight', 'bold', ...
        'FontColor', [1, 1, 1], 'FontSize', 10, ...
        'ButtonPushedFcn', @clearAllParticlesCallback);
    clearParticlesButton.Layout.Row = 3;
    clearParticlesButton.Layout.Column = [1, 3];
    allUIComponents.clearParticlesButton = clearParticlesButton;
    
    % Initialize disabled state
    disableAllControls();
end

% Create visualization controls with comprehensive dark mode support
function createVisualizationControlsWithDarkMode(parent)
    global flowToggle particleToggle vizPanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Visualization Controls Panel
    vizPanel = uipanel(parent, 'Title', 'VISUALIZATION CONTROLS', ...
                      'BackgroundColor', colors.panel_color, ...
                      'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    vizPanel.Layout.Row = 4;
    allUIComponents.vizPanel = vizPanel;
    
    vizGrid = uigridlayout(vizPanel, [3, 2]);
    vizGrid.ColumnWidth = {'fit', '1x'};
    vizGrid.RowHeight = {'fit', 'fit', 'fit'};
    vizGrid.Padding = [10, 10, 10, 10];
    vizGrid.RowSpacing = 8;
    vizGrid.BackgroundColor = colors.panel_color;
    allUIComponents.vizGrid = vizGrid;
    
    % Flow lines toggle
    flowLabel = uilabel(vizGrid, 'Text', 'Flow Lines:', ...
                       'FontWeight', 'bold', 'FontSize', 10, ...
                       'FontColor', colors.text_color);
    flowLabel.Layout.Row = 1;
    flowLabel.Layout.Column = 1;
    allUIComponents.flowLabel = flowLabel;
    
    flowToggle = uicheckbox(vizGrid, 'Text', 'Show', 'Value', true, ...
                           'FontColor', colors.text_color, ...
                           'ValueChangedFcn', @flowToggleCallback);
    flowToggle.Layout.Row = 1;
    flowToggle.Layout.Column = 2;
    allUIComponents.flowToggle = flowToggle;
    
    % Particle species toggle
    particleLabel = uilabel(vizGrid, 'Text', 'Particles:', ...
                           'FontWeight', 'bold', 'FontSize', 10, ...
                           'FontColor', colors.text_color);
    particleLabel.Layout.Row = 2;
    particleLabel.Layout.Column = 1;
    allUIComponents.particleLabel = particleLabel;
    
    particleToggle = uicheckbox(vizGrid, 'Text', 'Show', 'Value', true, ...
                               'FontColor', colors.text_color, ...
                               'ValueChangedFcn', @particleToggleCallback);
    particleToggle.Layout.Row = 2;
    particleToggle.Layout.Column = 2;
    allUIComponents.particleToggle = particleToggle;
    
    % Feed system toggle
    feedSystemLabel = uilabel(vizGrid, 'Text', 'Feed System:', ...
                             'FontWeight', 'bold', 'FontSize', 10, ...
                             'FontColor', colors.text_color);
    feedSystemLabel.Layout.Row = 3;
    feedSystemLabel.Layout.Column = 1;
    allUIComponents.feedSystemLabel = feedSystemLabel;
    
    feedSystemToggle = uicheckbox(vizGrid, 'Text', 'Show Pipes', 'Value', true, ...
                                 'FontColor', colors.text_color, ...
                                 'ValueChangedFcn', @feedSystemToggleCallback);
    feedSystemToggle.Layout.Row = 3;
    feedSystemToggle.Layout.Column = 2;
    allUIComponents.feedSystemToggle = feedSystemToggle;
end

% Create advanced sliders with comprehensive dark mode support
function createAdvancedSlidersWithDarkMode(parent)
    global densitySlider viscositySlider densityLabel viscosityLabel advPanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Advanced Parameters Panel
    advPanel = uipanel(parent, 'Title', 'ADVANCED PARAMETERS', ...
                      'BackgroundColor', colors.panel_color, ...
                      'ForegroundColor', colors.text_color, 'FontWeight', 'bold');
    advPanel.Layout.Row = 5;
    allUIComponents.advPanel = advPanel;
    
    advGrid = uigridlayout(advPanel, [4, 2]);
    advGrid.ColumnWidth = {'fit', '1x'};
    advGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    advGrid.Padding = [10, 10, 10, 10];
    advGrid.RowSpacing = 2;
    advGrid.BackgroundColor = colors.panel_color;
    allUIComponents.advGrid = advGrid;
    
    % Particle density slider
    densityTitleLabel = uilabel(advGrid, 'Text', 'Particle Density:', ...
                               'FontWeight', 'bold', 'FontSize', 10, ...
                               'FontColor', colors.text_color);
    densityTitleLabel.Layout.Row = 1;
    densityTitleLabel.Layout.Column = 1;
    allUIComponents.densityTitleLabel = densityTitleLabel;
    
    densitySlider = uislider(advGrid, 'Limits', [0.1, 2.0], 'Value', 1.0, ...
                            'MajorTicks', [0.1, 0.5, 1.0, 1.5, 2.0], ...
                            'ValueChangedFcn', @densitySliderCallback);
    densitySlider.Layout.Row = 2;
    densitySlider.Layout.Column = [1, 2];
    allUIComponents.densitySlider = densitySlider;
    
    densityLabel = uilabel(advGrid, 'Text', '1.00x', ...
                          'FontWeight', 'bold', 'FontSize', 10, ...
                          'FontColor', colors.accent_color, ...
                          'HorizontalAlignment', 'center');
    densityLabel.Layout.Row = 1;
    densityLabel.Layout.Column = 2;
    allUIComponents.densityLabel = densityLabel;
    
    % Fluid viscosity slider
    viscosityTitleLabel = uilabel(advGrid, 'Text', 'Fluid Viscosity:', ...
                                 'FontWeight', 'bold', 'FontSize', 10, ...
                                 'FontColor', colors.text_color);
    viscosityTitleLabel.Layout.Row = 3;
    viscosityTitleLabel.Layout.Column = 1;
    allUIComponents.viscosityTitleLabel = viscosityTitleLabel;
    
    viscositySlider = uislider(advGrid, 'Limits', [0.1, 3.0], 'Value', 1.0, ...
                              'MajorTicks', [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], ...
                              'ValueChangedFcn', @viscositySliderCallback);
    viscositySlider.Layout.Row = 4;
    viscositySlider.Layout.Column = [1, 2];
    allUIComponents.viscositySlider = viscositySlider;
    
    viscosityLabel = uilabel(advGrid, 'Text', '1.00x', ...
                            'FontWeight', 'bold', 'FontSize', 10, ...
                            'FontColor', colors.accent_color, ...
                            'HorizontalAlignment', 'center');
    viscosityLabel.Layout.Row = 3;
    viscosityLabel.Layout.Column = 2;
    allUIComponents.viscosityLabel = viscosityLabel;
    
    % Target conversion slider
conversionTitleLabel = uilabel(advGrid, 'Text', 'Target Conversion:', ...
                               'FontWeight', 'bold', 'FontSize', 10, ...
                               'FontColor', colors.text_color);
conversionTitleLabel.Layout.Row = 5;
conversionTitleLabel.Layout.Column = 1;
allUIComponents.conversionTitleLabel = conversionTitleLabel;

global targetConversionSlider targetConversionLabel target_conversion;
targetConversionSlider = uislider(advGrid, 'Limits', [0, 100], 'Value', 80, ...
                                 'MajorTicks', [0, 20, 40, 60, 80, 100], ...
                                 'ValueChangedFcn', @targetConversionSliderCallback);
targetConversionSlider.Layout.Row = 6;
targetConversionSlider.Layout.Column = [1, 2];
allUIComponents.targetConversionSlider = targetConversionSlider;

targetConversionLabel = uilabel(advGrid, 'Text', '80%', ...
                               'FontWeight', 'bold', 'FontSize', 10, ...
                               'FontColor', colors.accent_color, ...
                               'HorizontalAlignment', 'center');
targetConversionLabel.Layout.Row = 5;
targetConversionLabel.Layout.Column = 2;
allUIComponents.targetConversionLabel = targetConversionLabel;

% Initialize target conversion
target_conversion = 80;
end

% Create enhanced mode controls with comprehensive dark mode support
function createModeControlsEnhanced(parent)
    global darkModeToggle modeLabel modePanel allUIComponents;
    
    colors = getCurrentColorScheme();
    
    % Mode Control Panel - NOW AT TOP FOR VISIBILITY
    modePanel = uipanel(parent, 'Title', '🌓 DISPLAY MODE', ...
                       'BackgroundColor', colors.panel_color, ...
                       'ForegroundColor', colors.text_color, 'FontWeight', 'bold', ...
                       'FontSize', 12);
    modePanel.Layout.Row = 1;  % NOW AT TOP!
    allUIComponents.modePanel = modePanel;
    
    modeGrid = uigridlayout(modePanel, [1, 3]);
    modeGrid.ColumnWidth = {'1x', 'fit', 'fit'};
    modeGrid.RowHeight = {'fit'};
    modeGrid.Padding = [15, 15, 15, 15];
    modeGrid.RowSpacing = 8;
    modeGrid.BackgroundColor = colors.panel_color;
    allUIComponents.modeGrid = modeGrid;
    
    % Current mode label
    modeLabel = uilabel(modeGrid, 'Text', 'Bright Mode', ...
                       'FontWeight', 'bold', 'FontSize', 12, ...
                       'FontColor', colors.text_color);
    modeLabel.Layout.Row = 1;
    modeLabel.Layout.Column = 1;
    allUIComponents.modeLabel = modeLabel;
    
    % Mode switch label
    switchLabel = uilabel(modeGrid, 'Text', 'Dark:', ...
                         'FontWeight', 'bold', 'FontSize', 11, ...
                         'FontColor', colors.text_color, ...
                         'HorizontalAlignment', 'right');
    switchLabel.Layout.Row = 1;
    switchLabel.Layout.Column = 2;
    allUIComponents.switchLabel = switchLabel;
    
    % Dark/Bright mode toggle - LARGER AND MORE VISIBLE
    darkModeToggle = uiswitch(modeGrid, 'toggle', 'Value', 'Off', ...
                             'ValueChangedFcn', @enhancedDarkModeToggleCallback);
    darkModeToggle.Layout.Row = 1;
    darkModeToggle.Layout.Column = 3;
    allUIComponents.darkModeToggle = darkModeToggle;
end

% Draw enhanced reactor geometry with feed system
function drawEnhancedReactorGeometry(ax)
    global tank_diameter tank_height;
    global motor_graphics baffle_ring_graphics reactor_data;
    
    % Clear existing content
    cla(ax);
    
    % Create reactor geometry
    create_pilot_plant_reactor_geometry(ax);
    create_multiple_impeller_systems(ax);
    
    % Draw additional components
    motor_graphics = drawMotor(ax, tank_height);
    [baffle_ring_graphics] = drawBafflesWithRing(ax, tank_diameter/2, tank_height);
    
    % Draw enhanced feed system
    drawContinuousFeedSystem(ax);
    drawBasicCoolingCoil(ax, tank_diameter/2, tank_height);
end

% Draw continuous feed system (2 inlets + 1 outlet)
function drawContinuousFeedSystem(ax)
    global reactor_data dark_mode pipe_graphics;
    
    if isempty(reactor_data)
        return;
    end
    
    if dark_mode
        etac_pipe_color = [1.0, 0.8, 0.2];    % Bright yellow for EtOAc
        naoh_pipe_color = [0.4, 0.7, 1.0];    % Bright blue for NaOH  
        outlet_pipe_color = [0.8, 0.3, 1.0];  % Bright purple for outlet
        text_color = [0.9, 0.9, 0.95];
    else
        etac_pipe_color = [1.0, 0.7, 0.2];    % Orange for EtOAc
        naoh_pipe_color = [0.2, 0.4, 1.0];    % Blue for NaOH
        outlet_pipe_color = [0.6, 0.1, 0.9];  % Purple for outlet
        text_color = [0.2, 0.4, 1.0];
    end
    
    pipe_radius = 0.03;
    pipe_length = 0.4;
    
    % === INLET PIPE 1: Ethyl Acetate (Top Right) ===
    inlet1_pos = reactor_data.inlet_positions(1, :);
    
    % Horizontal section
    [Xp1, Yp1, Zp1] = cylinder(pipe_radius, 12);
    Yp1 = Yp1 * pipe_length;
    Zp1 = Zp1 * 0.05 + inlet1_pos(3) - 0.025;
    
    for i = 1:size(Xp1, 2)
        Xp1(:, i) = Xp1(:, i) + inlet1_pos(1) + (i-1) * pipe_length / size(Xp1, 2);
    end
    
    pipe_graphics.inlet1_horizontal = surf(ax, Xp1, Yp1, Zp1, ...
                                          'FaceColor', etac_pipe_color, ...
                                          'EdgeColor', 'none', ...
                                          'SpecularStrength', 0.8);
    
    % Vertical drop section
    [Xp1v, Yp1v, Zp1v] = cylinder(pipe_radius, 12);
    Xp1v = Xp1v + inlet1_pos(1);
    Yp1v = Yp1v + inlet1_pos(2);
    Zp1v = inlet1_pos(3) - Zp1v * 0.3;
    
    pipe_graphics.inlet1_vertical = surf(ax, Xp1v, Yp1v, Zp1v, ...
                                        'FaceColor', etac_pipe_color, ...
                                        'EdgeColor', 'none', ...
                                        'SpecularStrength', 0.8);
    
    % Pipe label
    text(ax, inlet1_pos(1) + 0.2, inlet1_pos(2), inlet1_pos(3) + 0.1, ...
         'EtOAc INLET →', 'Color', text_color, 'FontWeight', 'bold', 'FontSize', 10);
    
    % === INLET PIPE 2: NaOH (Top Left) ===
    inlet2_pos = reactor_data.inlet_positions(2, :);
    
    % Horizontal section
    [Xp2, Yp2, Zp2] = cylinder(pipe_radius, 12);
    Yp2 = Yp2 * pipe_length;
    Zp2 = Zp2 * 0.05 + inlet2_pos(3) - 0.025;
    
    for i = 1:size(Xp2, 2)
        Xp2(:, i) = Xp2(:, i) + inlet2_pos(1) - (i-1) * pipe_length / size(Xp2, 2);
    end
    
    pipe_graphics.inlet2_horizontal = surf(ax, Xp2, Yp2, Zp2, ...
                                          'FaceColor', naoh_pipe_color, ...
                                          'EdgeColor', 'none', ...
                                          'SpecularStrength', 0.8);
    
    % Vertical drop section
    [Xp2v, Yp2v, Zp2v] = cylinder(pipe_radius, 12);
    Xp2v = Xp2v + inlet2_pos(1);
    Yp2v = Yp2v + inlet2_pos(2);
    Zp2v = inlet2_pos(3) - Zp2v * 0.3;
    
    pipe_graphics.inlet2_vertical = surf(ax, Xp2v, Yp2v, Zp2v, ...
                                        'FaceColor', naoh_pipe_color, ...
                                        'EdgeColor', 'none', ...
                                        'SpecularStrength', 0.8);
    
    % Pipe label
    text(ax, inlet2_pos(1) - 0.2, inlet2_pos(2), inlet2_pos(3) + 0.1, ...
         '← NaOH INLET', 'Color', text_color, 'FontWeight', 'bold', 'FontSize', 10);
    
    % === OUTLET PIPE: Products (Bottom) ===
    outlet_pos = reactor_data.outlet_position;
    
    % Vertical section from tank bottom
    [Xpo, Ypo, Zpo] = cylinder(pipe_radius * 1.2, 12);
    Xpo = Xpo + outlet_pos(1);
    Ypo = Ypo + outlet_pos(2);
    Zpo = Zpo * 0.15 + outlet_pos(3);
    
    pipe_graphics.outlet_vertical = surf(ax, Xpo, Ypo, Zpo, ...
                                        'FaceColor', outlet_pipe_color, ...
                                        'EdgeColor', 'none', ...
                                        'SpecularStrength', 0.8);
    
    % Horizontal outlet section
    [Xpoh, Ypoh, Zpoh] = cylinder(pipe_radius * 1.2, 12);
    Xpoh = outlet_pos(1) + Xpoh * pipe_length;
    Ypoh = Ypoh + outlet_pos(2);
    Zpoh = Zpoh * 0.05 + outlet_pos(3);
    
    pipe_graphics.outlet_horizontal = surf(ax, Xpoh, Ypoh, Zpoh, ...
                                          'FaceColor', outlet_pipe_color, ...
                                          'EdgeColor', 'none', ...
                                          'SpecularStrength', 0.8);
    
    % Outlet label
    text(ax, outlet_pos(1) + 0.1, outlet_pos(2) + 0.1, outlet_pos(3) + 0.2, ...
         'PRODUCTS OUT ↓', 'Color', text_color, 'FontWeight', 'bold', 'FontSize', 10);
    
    % Flow rate indicators
    createFlowRateIndicators(ax);
    
   end

% Create flow rate indicators for pipes
function createFlowRateIndicators(ax)
    global reactor_data dark_mode;
    
    if dark_mode
        indicator_color = [0.8, 0.9, 1.0];
    else
        indicator_color = [0.3, 0.5, 0.8];
    end
    
    % Add flow direction arrows for inlet pipes
    for i = 1:reactor_data.inlet_pipe_count
        inlet_pos = reactor_data.inlet_positions(i, :);
        
        % Create arrow indicating flow direction (downward)
        arrow_x = [inlet_pos(1), inlet_pos(1)];
        arrow_y = [inlet_pos(2), inlet_pos(2)];
        arrow_z = [inlet_pos(3) + 0.05, inlet_pos(3) - 0.1];
        
        plot3(ax, arrow_x, arrow_y, arrow_z, ...
              'Color', indicator_color, 'LineWidth', 4, ...
              'Marker', 'v', 'MarkerSize', 8, 'MarkerFaceColor', indicator_color);
    end
    
    % Add flow direction arrow for outlet pipe
    outlet_pos = reactor_data.outlet_position;
    arrow_x = [outlet_pos(1), outlet_pos(1) + 0.15];
    arrow_y = [outlet_pos(2), outlet_pos(2)];
    arrow_z = [outlet_pos(3), outlet_pos(3)];
    
    plot3(ax, arrow_x, arrow_y, arrow_z, ...
          'Color', indicator_color, 'LineWidth', 4, ...
          'Marker', '>', 'MarkerSize', 8, 'MarkerFaceColor', indicator_color);
end

% Enhanced dark mode toggle callback with COMPLETE UI updates
function enhancedDarkModeToggleCallback(src, ~)
    global dark_mode modeLabel ax_handle allUIComponents;
    
    % Update dark mode state
    if strcmp(src.Value, 'On')
        dark_mode = true;
        fprintf('=== ACTIVATING COMPLETE DARK MODE ===\n');
    else
        dark_mode = false;
        fprintf('=== ACTIVATING COMPLETE BRIGHT MODE ===\n');
    end
    
    % Get new color scheme
    colors = getCurrentColorScheme();
    
    % Update mode label immediately
    if dark_mode
        modeLabel.Text = 'Dark Mode';
    else
        modeLabel.Text = 'Bright Mode';
    end
    modeLabel.FontColor = colors.text_color;
    
    % Update all UI components
    fprintf('Updating all UI components...\n');
    updateAllUIComponentsWithColors(colors);
    
    % ADD THIS LINE HERE:
refreshAllButtonStates();  % ADDED: Fix button states after theme change
    
    % Update figure background
    if isfield(allUIComponents, 'mainFig') && isvalid(allUIComponents.mainFig)
        allUIComponents.mainFig.Color = colors.bg_color;
    end
    
    % Update axes
    if ~isempty(ax_handle) && isvalid(ax_handle)
        setupAxesProperties(ax_handle);
    end
    
    % Update center panel
    if isfield(allUIComponents, 'centerPanel') && isvalid(allUIComponents.centerPanel)
        allUIComponents.centerPanel.BackgroundColor = colors.center_panel_color;
        allUIComponents.centerPanel.ForegroundColor = colors.center_text_color;
    end
    
    if isfield(allUIComponents, 'centerGrid') && isvalid(allUIComponents.centerGrid)
        allUIComponents.centerGrid.BackgroundColor = colors.center_panel_color;
    end
    
    % Force update ALL label colors including panel titles
    fprintf('Updating ALL label colors...\n');
    forceUpdateAllLabelColors();
    
    % Redraw reactor
    if ~isempty(ax_handle) && isvalid(ax_handle)
        redrawReactorWithPreservedParticles(ax_handle);
    end
    
    % Final verification and cleanup
    fprintf('Final verification...\n');
    if isfield(allUIComponents, 'mainFig') && isvalid(allUIComponents.mainFig)
        figure(allUIComponents.mainFig);  % Bring to front and refresh
        allUIComponents.mainFig.Color = colors.bg_color;
    end
    
    drawnow;
    pause(0.1);
    drawnow;
    
    if dark_mode
        fprintf('✓ DARK MODE ACTIVATED!\n');
    else
        fprintf('✓ BRIGHT MODE ACTIVATED!\n');
    end
end

% Update all UI components with new colors - COMPLETE COVERAGE
function updateAllUIComponentsWithColors(colors)
    global allUIComponents;
    
    % Update main figure and ALL grids
    if isfield(allUIComponents, 'mainFig') && isvalid(allUIComponents.mainFig)
        allUIComponents.mainFig.Color = colors.bg_color;
    end
    
    if isfield(allUIComponents, 'mainGrid') && isvalid(allUIComponents.mainGrid)
        allUIComponents.mainGrid.BackgroundColor = colors.bg_color;
    end
    
    if isfield(allUIComponents, 'centerGrid') && isvalid(allUIComponents.centerGrid)
        allUIComponents.centerGrid.BackgroundColor = colors.center_panel_color;
    end
    
panel_fields = {'leftPanel', 'centerPanel', 'combinedPanel', 'rpmPanel', ...
               'buttonPanel', 'vizPanel', 'advPanel', 'feedPanel', 'residencePanel'};
    
    for i = 1:length(panel_fields)
        field = panel_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            if strcmp(field, 'centerPanel')
                allUIComponents.(field).BackgroundColor = colors.center_panel_color;
                allUIComponents.(field).ForegroundColor = colors.center_text_color;
            else
                allUIComponents.(field).BackgroundColor = colors.panel_color;
                allUIComponents.(field).ForegroundColor = colors.text_color;
            end
        end
    end
    
    % Update ALL internal grids within panels
grid_fields = {'leftGrid', 'centerGrid', 'combinedGrid', 'powerSubGrid', 'modeSubGrid', ...
               'rpmGrid', 'buttonGrid', 'vizGrid', 'advGrid', 'feedGrid', 'residenceGrid'};
    for i = 1:length(grid_fields)
        field = grid_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                if contains(field, 'center') || contains(field, 'Center')
                    allUIComponents.(field).BackgroundColor = colors.center_panel_color;
                else
                    allUIComponents.(field).BackgroundColor = colors.panel_color;
                end
            catch
                % Some grids may not support BackgroundColor
            end
        end
    end
    
    % Update ALL labels with comprehensive coverage
label_fields = {'powerLabel', 'statusLabel', 'gaugeLabel', 'knobLabel', ...
               'flowLabel', 'particleLabel', 'densityTitleLabel', ...
               'viscosityTitleLabel', 'feedLabel', 'residenceTitleLabel', ...
               'feedRateTitleLabel', 'feedStatusLabel', 'feedSystemLabel', ...
               'switchLabel', 'modeLabel'};
    
    for i = 1:length(label_fields)
        field = label_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.text_color;
            catch
                % Continue if error
            end
        end
    end
    
    % Update special labels with accent colors
   accent_label_fields = {'rpmValueLabel', 'densityLabel', 'viscosityLabel', 'feedRateLabel', 'residenceTimeLabel'};
    for i = 1:length(accent_label_fields)
        field = accent_label_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.accent_color;
            catch
                % Continue if error
            end
        end
    end
    
    % Update ALL checkboxes and toggles
    checkbox_fields = {'flowToggle', 'particleToggle', 'feedToggle', 'feedSystemToggle'};
    for i = 1:length(checkbox_fields)
        field = checkbox_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.text_color;
            catch
                % Continue if error
            end
        end
    end
    
    % Update ALL sliders
    slider_fields = {'densitySlider', 'viscositySlider', 'feedRateSlider'};
    for i = 1:length(slider_fields)
        field = slider_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.text_color;
            catch
                % Ignore if not supported
            end
        end
    end
    
    % Update gauges and knobs
    if isfield(allUIComponents, 'rpmGauge') && isvalid(allUIComponents.rpmGauge)
        allUIComponents.rpmGauge.BackgroundColor = colors.button_bg;
        allUIComponents.rpmGauge.ScaleColors = colors.accent_color;
    end
    
    % Update ALL buttons with theme-appropriate colors
    button_configs = {
        {'startButton', colors.success_color}, 
        {'pauseButton', colors.warning_color}, 
        {'stopButton', colors.error_color},
        {'resetButton', colors.accent_color},
        {'startFeedButton', colors.accent_color},
        {'emergencyButton', colors.error_color},
        {'clearParticlesButton', colors.warning_color}
    };
    
    for i = 1:length(button_configs)
        field = button_configs{i}{1};
        button_color = button_configs{i}{2};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            allUIComponents.(field).BackgroundColor = button_color;
            allUIComponents.(field).FontColor = [1, 1, 1];
        end
    end
    
    fprintf('✓ Updated ALL UI components with %s theme colors\n', ...
            iif(colors.bg_color(1) < 0.5, 'DARK', 'BRIGHT'));
end

% Update ALL font colors comprehensively 
function updateAllFontColors(colors)
    global allUIComponents;
    
    % Update ALL regular labels with text color
    regular_label_fields = {'powerLabel', 'statusLabel', 'gaugeLabel', 'knobLabel', ...
                       'flowLabel', 'particleLabel', 'feedLabel', 'switchLabel', ...
                       'densityTitleLabel', 'viscosityTitleLabel', 'feedRateTitleLabel', ...
                       'feedStatusLabel', 'feedSystemLabel'};
    
 for i = 1:length(regular_label_fields)
    field = regular_label_fields{i};
    if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
        try
            allUIComponents.(field).FontColor = colors.text_color;
        catch ME
            fprintf('Warning: Could not update %s FontColor: %s\n', field, ME.message);
        end
    else
        fprintf('Warning: Label %s not found or invalid\n', field);
    end
end
    
    % Update accent color labels (value displays)
    accent_label_fields = {'rpmValueLabel', 'densityLabel', 'viscosityLabel', 'feedRateLabel'};
    for i = 1:length(accent_label_fields)
        field = accent_label_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.accent_color;
            catch
                % Continue if error
            end
        end
    end
    
    % Update status labels with appropriate colors
    if isfield(allUIComponents, 'feedStatusLabel') && isvalid(allUIComponents.feedStatusLabel)
        try
            allUIComponents.feedStatusLabel.FontColor = colors.success_color;
        catch
            % Continue if error
        end
    end
    
    % Update checkboxes and toggle text
    checkbox_fields = {'flowToggle', 'particleToggle', 'feedToggle', 'feedSystemToggle'};
    for i = 1:length(checkbox_fields)
        field = checkbox_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            try
                allUIComponents.(field).FontColor = colors.text_color;
            catch
                % Continue if error
            end
        end
    end
    
    % Update mode label specifically
    if isfield(allUIComponents, 'modeLabel') && isvalid(allUIComponents.modeLabel)
        try
            allUIComponents.modeLabel.FontColor = colors.text_color;
        catch
            % Continue if error
        end
    end
    
    % Update any remaining text elements by finding all labels in the figure
    if isfield(allUIComponents, 'mainFig') && isvalid(allUIComponents.mainFig)
        try
            % Find all label objects in the figure
            all_labels = findall(allUIComponents.mainFig, 'Type', 'uilabel');
            for i = 1:length(all_labels)
                try
                    % Check if this is an accent label (contains numbers or specific keywords)
                    label_text = all_labels(i).Text;
                    if contains(label_text, 'RPM') || contains(label_text, 'x') || ...
                       contains(label_text, '%') || contains(label_text, 'Mode')
                        % Use accent color for value displays and mode
                        if contains(label_text, 'Mode')
                            all_labels(i).FontColor = colors.text_color;
                        else
                            all_labels(i).FontColor = colors.accent_color;
                        end
                    else
                        % Use regular text color for other labels
                        all_labels(i).FontColor = colors.text_color;
                    end
                catch
                    % Continue if individual label fails
                end
            end
            
            % Find all checkbox objects and update their text
            all_checkboxes = findall(allUIComponents.mainFig, 'Type', 'uicheckbox');
            for i = 1:length(all_checkboxes)
                try
                    all_checkboxes(i).FontColor = colors.text_color;
                catch
                    % Continue if individual checkbox fails
                end
            end
            
        catch
            % Continue if figure search fails
        end
    end
    
    fprintf('✓ Updated ALL font colors for %s theme\n', ...
            iif(colors.bg_color(1) < 0.5, 'DARK', 'BRIGHT'));
end

% Utility function for conditional values
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

% Create pilot plant reactor geometry
function create_pilot_plant_reactor_geometry(ax)
    global reactor_data dark_mode;
    
    if dark_mode
        wall_color = [0.25, 0.25, 0.35];
        bottom_color = [0.2, 0.2, 0.3];
        fluid_color = [0.2, 0.5, 0.8];
    else
        wall_color = [0.7, 0.7, 0.8];
        bottom_color = [0.6, 0.6, 0.7];
        fluid_color = [0.3, 0.5, 0.8];
    end
    
    % Reactor vessel walls
    [theta, z] = meshgrid(linspace(0, 2*pi, 60), linspace(0, reactor_data.tank_height, 40));
    x_wall = reactor_data.tank_radius * cos(theta);
    y_wall = reactor_data.tank_radius * sin(theta);
    
    reactor_data.tank_walls = surf(ax, x_wall, y_wall, z, ...
                                  'FaceColor', wall_color, ...
                                  'FaceAlpha', 0.3, ...
                                  'EdgeColor', 'none', ...
                                  'AmbientStrength', 0.4, ...
                                  'DiffuseStrength', 0.6, ...
                                  'SpecularStrength', 0.8);
    
    % Reactor bottom
    [r, theta_bottom] = meshgrid(linspace(0, reactor_data.tank_radius, 25), ...
                                linspace(0, 2*pi, 50));
    x_bottom = r .* cos(theta_bottom);
    y_bottom = r .* sin(theta_bottom);
    z_bottom = zeros(size(x_bottom));
    
    reactor_data.tank_bottom = surf(ax, x_bottom, y_bottom, z_bottom, ...
                                   'FaceColor', bottom_color, ...
                                   'FaceAlpha', 0.9, ...
                                   'EdgeColor', 'none', ...
                                   'SpecularStrength', 0.7);
    
    % Liquid surface
    reactor_data.fluid_surface = surf(ax, x_bottom*0.98, y_bottom*0.98, ...
                                     ones(size(x_bottom))*reactor_data.fluid_level, ...
                                     'FaceColor', fluid_color, ...
                                     'FaceAlpha', 0.7, ...
                                     'EdgeColor', 'none', ...
                                     'AmbientStrength', 0.3, ...
                                     'SpecularStrength', 0.9);
    
    % Industrial supports
    create_industrial_supports(ax);
    create_reactor_fittings(ax);
end

% Create industrial support structure
function create_industrial_supports(ax)
    global reactor_data dark_mode;
    
    if dark_mode
        support_color = [0.35, 0.35, 0.45];
        platform_color = [0.3, 0.3, 0.4];
    else
        support_color = [0.4, 0.4, 0.5];
        platform_color = [0.3, 0.3, 0.4];
    end
    
    support_height = reactor_data.tank_height + 1.0;
    support_radius = 0.08;
    
    % Main support columns
    for i = 1:6
        angle = (i-1) * pi/3;
        x_pos = (reactor_data.tank_radius + 0.4) * cos(angle);
        y_pos = (reactor_data.tank_radius + 0.4) * sin(angle);
        
        [theta_s, z_s] = meshgrid(linspace(0, 2*pi, 12), ...
                                 linspace(0, support_height, 25));
        x_support = x_pos + support_radius * cos(theta_s);
        y_support = y_pos + support_radius * sin(theta_s);
        
        surf(ax, x_support, y_support, z_s, ...
             'FaceColor', support_color, ...
             'EdgeColor', 'none', ...
             'SpecularStrength', 0.6);
    end
    
    % Platform at working level
    platform_height = reactor_data.fluid_level + 0.5;
    platform_radius = reactor_data.tank_radius + 0.6;
    
    [r_plat, theta_plat] = meshgrid(linspace(reactor_data.tank_radius + 0.2, platform_radius, 10), ...
                                   linspace(0, 2*pi, 40));
    x_platform = r_plat .* cos(theta_plat);
    y_platform = r_plat .* sin(theta_plat);
    z_platform = platform_height * ones(size(x_platform));
    
    surf(ax, x_platform, y_platform, z_platform, ...
         'FaceColor', platform_color, ...
         'FaceAlpha', 0.8, ...
         'EdgeColor', 'none');
end

% Add reactor nozzles and fittings
function create_reactor_fittings(ax)
    global reactor_data dark_mode;
    
    if dark_mode
        text_color = [0.9, 0.9, 0.95];
        fitting_color = [0.45, 0.45, 0.55];
    else
        text_color = [0.2, 0.4, 1.0];
        fitting_color = [0.5, 0.5, 0.6];
    end
    
    % Additional fittings for continuous operation
    feed_height = reactor_data.fluid_level * 0.7;
    [theta_nozzle, z_nozzle] = meshgrid(linspace(0, 2*pi, 12), ...
                                       linspace(0, 0.2, 8));
    x_nozzle = reactor_data.tank_radius + 0.05 + 0.06 * cos(theta_nozzle);
    y_nozzle = 0.06 * sin(theta_nozzle);
    z_nozzle = feed_height + z_nozzle;
    
    surf(ax, x_nozzle, y_nozzle, z_nozzle, ...
         'FaceColor', fitting_color, ...
         'EdgeColor', 'none');
    
    text(ax, reactor_data.tank_radius + 0.15, 0, feed_height, 'CONT. FEED →', ...
         'Color', text_color, 'FontWeight', 'bold', 'FontSize', 9);
end

% Create multiple impeller systems
function create_multiple_impeller_systems(ax)
    global reactor_data dark_mode;
    
    if dark_mode
        shaft_color = [0.2, 0.2, 0.3];
    else
        shaft_color = [0.2, 0.2, 0.3];
    end
    
    % Impeller shaft
    shaft_radius = 0.05;
    shaft_height = reactor_data.fluid_level + 1.5;
    
    [theta_shaft, z_shaft] = meshgrid(linspace(0, 2*pi, 20), ...
                                     linspace(0, shaft_height, 40));
    x_shaft = shaft_radius * cos(theta_shaft);
    y_shaft = shaft_radius * sin(theta_shaft);
    
    reactor_data.impeller_shafts = surf(ax, x_shaft, y_shaft, z_shaft, ...
                                       'FaceColor', shaft_color, ...
                                       'EdgeColor', 'none', ...
                                       'AmbientStrength', 0.3, ...
                                       'SpecularStrength', 0.9);
    
    % Create Rushton turbines
    for imp_idx = 1:reactor_data.impeller_count
        impeller_height = reactor_data.impeller_heights(imp_idx);
        create_rushton_turbine_at_height(ax, imp_idx, impeller_height);
        create_impeller_hub_at_height(ax, impeller_height);
    end
end

% Create Rushton turbine blades
function create_rushton_turbine_at_height(ax, impeller_index, height)
    global reactor_data dark_mode;
    
    blade_length = reactor_data.impeller_diameter / 2 * 0.8;
    blade_height = reactor_data.impeller_diameter / 5;
    
    if dark_mode
        blade_colors = {[0.3, 0.8, 1.0], [0.8, 0.3, 1.0], [1.0, 0.8, 0.3]};
    else
        blade_colors = {[0.1, 0.6, 0.9], [0.6, 0.1, 0.9], [0.9, 0.6, 0.1]};
    end
    
    blade_color = blade_colors{min(impeller_index, 3)};
    
    % Initialize blade array
    reactor_data.impeller_blades{impeller_index} = gobjects(reactor_data.blade_count, 1);
    
    for i = 1:reactor_data.blade_count
        blade_angle = (i-1) * 2*pi / reactor_data.blade_count;
        
        % Create blade geometry
        blade_r_inner = 0.08;
        blade_r_outer = blade_r_inner + blade_length;
        
        [r_blade, z_blade] = meshgrid(linspace(blade_r_inner, blade_r_outer, 15), ...
                                     linspace(height - blade_height/2, ...
                                             height + blade_height/2, 12));
        
        theta_blade = blade_angle * ones(size(r_blade));
        x_blade = r_blade .* cos(theta_blade);
        y_blade = r_blade .* sin(theta_blade);
        
        reactor_data.impeller_blades{impeller_index}(i) = surf(ax, x_blade, y_blade, z_blade, ...
                                                              'FaceColor', blade_color, ...
                                                              'EdgeColor', 'none', ...
                                                              'AmbientStrength', 0.4, ...
                                                              'SpecularStrength', 0.8);
        
        % Store original data
        reactor_data.blade_original_vertices{impeller_index, i} = struct('r', r_blade, 'theta', blade_angle, 'z', z_blade);
    end
    
    % Add impeller label
    if dark_mode
        label_color = blade_color;
    else
        label_color = blade_color * 0.8;
    end
    
    label_text = sprintf('Impeller %d', impeller_index);
    text(ax, reactor_data.tank_radius * 0.7, 0, height + 0.2, label_text, ...
         'FontWeight', 'bold', 'FontSize', 8, 'Color', label_color);
end

% Create impeller hub
function create_impeller_hub_at_height(ax, height)
    global dark_mode;
    
    if dark_mode
        hub_color = [0.4, 0.4, 0.5];
    else
        hub_color = [0.4, 0.4, 0.5];
    end
    
    hub_radius = 0.15;
    [theta_hub, phi_hub] = meshgrid(linspace(0, 2*pi, 24), linspace(0, pi, 18));
    x_hub = hub_radius * sin(phi_hub) .* cos(theta_hub);
    y_hub = hub_radius * sin(phi_hub) .* sin(theta_hub);
    z_hub = height + hub_radius * cos(phi_hub);
    
    surf(ax, x_hub, y_hub, z_hub, ...
         'FaceColor', hub_color, ...
         'EdgeColor', 'none', ...
         'SpecularStrength', 0.8);
end

% Draw motor
function motor_handle = drawMotor(ax, tankHeight)
    global dark_mode;
    
    if dark_mode
        motor_color = [0.4, 0.4, 0.5];
        fin_color = [0.35, 0.35, 0.45];
        text_color = [0.7, 0.9, 1.0];
    else
        motor_color = [0.2, 0.3, 0.4];
        fin_color = [0.3, 0.3, 0.4];
        text_color = [0.6, 0.8, 1.0];
    end
    
    % Motor housing
    motor_radius = 0.4;
    motor_height = 0.8;
    motor_z = tankHeight + 0.2;
    
    [theta_motor, z_motor] = meshgrid(linspace(0, 2*pi, 24), ...
                                     linspace(motor_z, motor_z + motor_height, 15));
    x_motor = motor_radius * cos(theta_motor);
    y_motor = motor_radius * sin(theta_motor);
    
    motor_handle = surf(ax, x_motor, y_motor, z_motor, ...
         'FaceColor', motor_color, ...
         'EdgeColor', 'none', ...
         'SpecularStrength', 0.7);
    
    % Motor cooling fins
    fin_count = 12;
    for i = 1:fin_count
        angle = (i-1) * 2*pi / fin_count;
        fin_length = 0.12;
        fin_thickness = 0.015;
        
        [theta_fin, z_fin] = meshgrid(linspace(-fin_thickness, fin_thickness, 3), ...
                                     linspace(motor_z + 0.15, motor_z + 0.65, 6));
        x_fin = (motor_radius + fin_length * linspace(0, 1, size(theta_fin, 2))) .* cos(angle + theta_fin);
        y_fin = (motor_radius + fin_length * linspace(0, 1, size(theta_fin, 2))) .* sin(angle + theta_fin);
        
        surf(ax, x_fin, y_fin, z_fin, ...
             'FaceColor', fin_color, ...
             'EdgeColor', 'none');
    end
    
    % Motor label
    text(ax, 0.5, 0, motor_z + motor_height/2, 'HIGH-SPEED MOTOR', ...
         'FontWeight', 'bold', 'FontSize', 8, 'Color', text_color);
end

% Draw baffles with support ring
function ring_handle = drawBafflesWithRing(ax, tankRadius, tankHeight)
    global dark_mode;
    
    if dark_mode
        baffle_color = [0.35, 0.35, 0.45];
        ring_color = [0.35, 0.35, 0.45];
        text_color = [0.6, 0.6, 0.8];
    else
        baffle_color = [0.2, 0.6, 0.9];
        ring_color = [0.3, 0.3, 0.4];
        text_color = [0.4, 0.4, 0.6];
    end
    
    numBaffles = 4;
    baffleWidth = tankRadius / 10;
    ring_height = tankHeight / 2;
    ring_thickness = 0.05;
    ring_width = 0.08;
    
    % Draw baffles
    for i = 1:numBaffles
        angle = (i-1) * (2*pi/numBaffles);
        x_outer = tankRadius * cos(angle);
        y_outer = tankRadius * sin(angle);
        x_inner = (tankRadius - baffleWidth) * cos(angle);
        y_inner = (tankRadius - baffleWidth) * sin(angle);
        
        try
            patch(ax, ...
                [x_outer, x_inner, x_inner, x_outer], ...
                [y_outer, y_inner, y_inner, y_outer], ...
                [0, 0, tankHeight, tankHeight], ...
                baffle_color, 'EdgeColor', baffle_color * 0.8, ...
                'FaceAlpha', 0.9, 'LineWidth', 2);
        catch ME
            fprintf('Warning: Could not create baffle %d: %s\n', i, ME.message);
        end
    end
    
    % Draw support ring
    ring_outer_radius = tankRadius - baffleWidth/2;
    ring_inner_radius = ring_outer_radius - ring_width;
    
    [theta_ring, z_ring] = meshgrid(linspace(0, 2*pi, 50), ...
                                   linspace(ring_height - ring_thickness/2, ...
                                           ring_height + ring_thickness/2, 10));
    x_ring_outer = ring_outer_radius * cos(theta_ring);
    y_ring_outer = ring_outer_radius * sin(theta_ring);
    
    ring_handle = surf(ax, x_ring_outer, y_ring_outer, z_ring, ...
         'FaceColor', ring_color, ...
         'FaceAlpha', 0.95, ...
         'EdgeColor', 'none', ...
         'SpecularStrength', 0.8);
    
    % Ring inner surface
    x_ring_inner = ring_inner_radius * cos(theta_ring);
    y_ring_inner = ring_inner_radius * sin(theta_ring);
    
    surf(ax, x_ring_inner, y_ring_inner, z_ring, ...
         'FaceColor', ring_color, ...
         'FaceAlpha', 0.95, ...
         'EdgeColor', 'none', ...
         'SpecularStrength', 0.8);
    
    % Add label
    text(ax, ring_outer_radius + 0.1, 0, ring_height, 'Support Ring', ...
         'FontWeight', 'bold', 'FontSize', 8, 'Color', text_color);
end

% Draw cooling coil
function drawBasicCoolingCoil(ax, tankRadius, tankHeight)
    global dark_mode;
    
    if dark_mode
        coil_color = [0.9, 0.6, 0.3];
    else
        coil_color = [0.8, 0.35, 0.1];
    end
    
    coil_radius = tankRadius * 0.50;
    num_turns = 10;
    points_per_turn = 60;
    
    total_points = num_turns * points_per_turn;
    theta = linspace(0, 2*pi*num_turns, total_points);
    
    x_coil = coil_radius * cos(theta);
    y_coil = coil_radius * sin(theta);
    z_coil = linspace(tankHeight * 0.92, tankHeight * 0.08, total_points);
    
    try
        plot3(ax, x_coil, y_coil, z_coil, ...
              'Color', coil_color, ...
              'LineWidth', 5, ...
              'LineStyle', '-');
    catch ME
        fprintf('Warning: Could not create cooling coil: %s\n', ME.message);
    end
end

% Redraw reactor while preserving particles
function redrawReactorWithPreservedParticles(ax)
    global reactor_data;
    
    % Store current particle data and feed system state
    if ~isempty(reactor_data) && isfield(reactor_data, 'species')
        stored_species = reactor_data.species;
        stored_conversion = reactor_data.conversion;
        stored_time = reactor_data.time;
        stored_frame_count = reactor_data.frame_count;
        stored_impeller_angles = reactor_data.impeller_angles;
        stored_total_fed = reactor_data.total_fed;
        stored_total_exited = reactor_data.total_exited;
        stored_current_holdup = reactor_data.current_holdup;
    else
        stored_species = [];
        stored_conversion = 0;
        stored_time = 0;
        stored_frame_count = 0;
        stored_impeller_angles = [];
        stored_total_fed = 0;
        stored_total_exited = 0;
        stored_current_holdup = 0;
    end
    
    % Redraw reactor geometry
    drawEnhancedReactorGeometry(ax);
    
    % Restore particle data and feed system state
    if ~isempty(stored_species)
        reactor_data.species = stored_species;
        reactor_data.conversion = stored_conversion;
        reactor_data.time = stored_time;
        reactor_data.frame_count = stored_frame_count;
        reactor_data.total_fed = stored_total_fed;
        reactor_data.total_exited = stored_total_exited;
        reactor_data.current_holdup = stored_current_holdup;
        
        if ~isempty(stored_impeller_angles)
            reactor_data.impeller_angles = stored_impeller_angles;
        end
        
        % Recreate visualizations
        createSpeciesVisualization(ax);
        createFlowVisualization(ax);
    end
end

% === SIMULATION AND CALLBACK FUNCTIONS ===

% Start fluid simulation with continuous feed
function startFluidSimulation(~, ~)
    global isPaused isStopped powerSwitch ax_handle isRunning reactor_data;
    global continuous_feeding;
    
    fprintf('=== STARTING CONTINUOUS FEED CSTR SIMULATION ===\n');
    
    if strcmp(powerSwitch.Value, 'Off')
        powerSwitch.Value = 'On';
        enableAllControls();
    end
    
    if ~isempty(ax_handle)
        isPaused = false;
        isStopped = false;
        isRunning = false;
        
        try
            % Initialize reactor
            initializeReactorParameters();
            initializeContinuousFeedSystem();
            
            % Create visualizations
            createSpeciesVisualization(ax_handle);
            createFlowVisualization(ax_handle);
            
            % Run simulation
            runFluidSimulation(ax_handle);
            
            fprintf('Continuous feed simulation started successfully!\n');
        catch ME
            fprintf('Error starting simulation: %s\n', ME.message);
        end
    end
end

% Initialize continuous feed system
function initializeContinuousFeedSystem()
    global reactor_data continuous_feeding feed_rate_factor;
    global total_particles_fed total_particles_exited;
    
    if isempty(reactor_data)
        initializeReactorParameters();
    end
    
    % Reset feed system counters
    total_particles_fed = 0;
    total_particles_exited = 0;
    reactor_data.total_fed = 0;
    reactor_data.total_exited = 0;
    reactor_data.current_holdup = 0;
    reactor_data.particles_in_system = 0;
    
    % Initialize feed rates based on slider settings
    base_feed_rate = 10.0;  % particles per second
    reactor_data.feed_rate = base_feed_rate * feed_rate_factor;
    reactor_data.exit_rate = base_feed_rate * 0.8 * feed_rate_factor;  % Slightly less than feed
    
    % Start with some initial particles if continuous feeding is enabled
    if continuous_feeding
        initializeInitialParticles();
    end
    
    fprintf('Continuous feed system initialized: Feed=%.1f p/s, Exit=%.1f p/s\n', ...
            reactor_data.feed_rate, reactor_data.exit_rate);
end

% Initialize some initial particles in the reactor
function initializeInitialParticles()
    global reactor_data;
    
    % Start with moderate number of initial particles
    initial_etac_particles = 100;
    initial_naoh_particles = 100;
    
    % Initialize Ethyl Acetate particles
    angles_ea = 2*pi * rand(initial_etac_particles, 1);
    radii_ea = reactor_data.tank_radius * 0.8 * sqrt(rand(initial_etac_particles, 1));
    heights_ea = reactor_data.fluid_level * 0.6 + reactor_data.fluid_level * 0.3 * rand(initial_etac_particles, 1);
    
    % FIXED: Initialize as column vectors
    reactor_data.species.ethyl_acetate.x = radii_ea .* cos(angles_ea);
    reactor_data.species.ethyl_acetate.y = radii_ea .* sin(angles_ea);
    reactor_data.species.ethyl_acetate.z = heights_ea;
    reactor_data.species.ethyl_acetate.vx = zeros(initial_etac_particles, 1);
    reactor_data.species.ethyl_acetate.vy = zeros(initial_etac_particles, 1);
    reactor_data.species.ethyl_acetate.vz = zeros(initial_etac_particles, 1);
    reactor_data.species.ethyl_acetate.age = zeros(initial_etac_particles, 1);
    
    % Initialize NaOH particles
    angles_naoh = 2*pi * rand(initial_naoh_particles, 1);
    radii_naoh = reactor_data.tank_radius * 0.8 * sqrt(rand(initial_naoh_particles, 1));
    heights_naoh = reactor_data.fluid_level * 0.5 * rand(initial_naoh_particles, 1) + 0.2;
    
    % FIXED: Initialize as column vectors
    reactor_data.species.naoh.x = radii_naoh .* cos(angles_naoh);
    reactor_data.species.naoh.y = radii_naoh .* sin(angles_naoh);
    reactor_data.species.naoh.z = heights_naoh;
    reactor_data.species.naoh.vx = zeros(initial_naoh_particles, 1);
    reactor_data.species.naoh.vy = zeros(initial_naoh_particles, 1);
    reactor_data.species.naoh.vz = zeros(initial_naoh_particles, 1);
    reactor_data.species.naoh.age = zeros(initial_naoh_particles, 1);
    
    % FIXED: Initialize empty products with proper structure
    reactor_data.species.products = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', [], 'age', []);
    
    % Update particle count
    reactor_data.current_holdup = initial_etac_particles + initial_naoh_particles;
    reactor_data.particles_in_system = reactor_data.current_holdup;
    
    fprintf('Initial particles loaded: %d EtOAc, %d NaOH (all as column vectors)\n', ...
            initial_etac_particles, initial_naoh_particles);
end

% Initialize reactor parameters
function initializeReactorParameters()
    global reactor_data particle_density fluid_viscosity_factor;
    
    if isempty(reactor_data)
        reactor_data = initialize_pilot_plant_reactor_parameters();
    end
    
    fprintf('Reactor initialized: %.1fm tank, %d impellers, continuous feed system\n', ...
            reactor_data.tank_radius*2, reactor_data.impeller_count);
end

% Create species visualization
function createSpeciesVisualization(ax)
    global reactor_data show_particles;
    
    if isempty(reactor_data) || isempty(reactor_data.species)
        return;
    end
    
    % Ethyl Acetate - Yellow particles (representing EtOAc feed)
    if ~isempty(reactor_data.species.ethyl_acetate.x)
        reactor_data.species_plots.ethyl_acetate = scatter3(ax, ...
            reactor_data.species.ethyl_acetate.x, ...
            reactor_data.species.ethyl_acetate.y, ...
            reactor_data.species.ethyl_acetate.z, ...
            25, [1.0, 0.8, 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
    else
        % Initialize empty plot for ethyl acetate
        reactor_data.species_plots.ethyl_acetate = scatter3(ax, [], [], [], ...
            25, [1.0, 0.8, 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
    end
    
    if ~show_particles
        reactor_data.species_plots.ethyl_acetate.Visible = 'off';
    end
    
    % NaOH - Blue particles (representing NaOH feed)
    if ~isempty(reactor_data.species.naoh.x)
        reactor_data.species_plots.naoh = scatter3(ax, ...
            reactor_data.species.naoh.x, ...
            reactor_data.species.naoh.y, ...
            reactor_data.species.naoh.z, ...
            25, [0.2, 0.4, 1.0], 'filled', 'MarkerFaceAlpha', 0.8);
    else
        % Initialize empty plot for NaOH
        reactor_data.species_plots.naoh = scatter3(ax, [], [], [], ...
            25, [0.2, 0.4, 1.0], 'filled', 'MarkerFaceAlpha', 0.8);
    end
    
    if ~show_particles
        reactor_data.species_plots.naoh.Visible = 'off';
    end
    
    % Products - Green particles (representing product formation)
    % FIX: Always initialize the products plot, even when empty
    if ~isempty(reactor_data.species.products.x)
        reactor_data.species_plots.products = scatter3(ax, ...
            reactor_data.species.products.x, ...
            reactor_data.species.products.y, ...
            reactor_data.species.products.z, ...
            30, [0.2, 0.8, 0.3], 'filled', 'MarkerFaceAlpha', 0.9);
    else
        % Initialize empty plot for products - THIS IS THE KEY FIX
        reactor_data.species_plots.products = scatter3(ax, [], [], [], ...
            30, [0.2, 0.8, 0.3], 'filled', 'MarkerFaceAlpha', 0.9);
        fprintf('Products plot initialized (empty) with GREEN color\n');
    end
    
    if ~show_particles
        reactor_data.species_plots.products.Visible = 'off';
    end
end

% Create flow visualization
function createFlowVisualization(ax)
    global reactor_data show_flow_lines dark_mode;
    
    if isempty(reactor_data)
        return;
    end
    
    if dark_mode
        flow_colors = {[0.3, 0.9, 1.0], [0.9, 0.3, 1.0], [1.0, 0.9, 0.3]};
    else
        flow_colors = {[0, 0.8, 1], [0.8, 0, 1], [1, 0.8, 0]};
    end
    
    % Enhanced flow lines for multiple impellers
    flow_line_count = 20;
    total_flow_lines = flow_line_count * reactor_data.impeller_count;
    reactor_data.flow_lines = gobjects(total_flow_lines, 1);
    
    line_idx = 1;
    for imp = 1:reactor_data.impeller_count
        impeller_height = reactor_data.impeller_heights(imp);
        flow_color = flow_colors{min(imp, 3)};
        
        for i = 1:flow_line_count
            angle_offset = (i-1) * 2*pi / flow_line_count;
            
            % Create 3D flow patterns
            t = linspace(0, 3*pi, 60);
            radius_variation = 0.8 + 0.3 * sin(t * 2);
            base_radius = reactor_data.impeller_diameter * 0.6;
            
            x_flow = base_radius * radius_variation .* cos(angle_offset + t) .* exp(-t*0.05);
            y_flow = base_radius * radius_variation .* sin(angle_offset + t) .* exp(-t*0.05);
            z_flow = impeller_height + 0.4 * sin(t) .* exp(-t*0.1);
            
            try
                reactor_data.flow_lines(line_idx) = plot3(ax, x_flow, y_flow, z_flow, ...
                                                          'Color', flow_color, 'LineWidth', 2, ...
                                                          'LineStyle', '-');
                
                if ~show_flow_lines
                    reactor_data.flow_lines(line_idx).Visible = 'off';
                end
            catch ME
                fprintf('Warning: Could not create flow line %d: %s\n', line_idx, ME.message);
            end
            
            line_idx = line_idx + 1;
        end
    end
end

% Run fluid simulation with continuous feed
function runFluidSimulation(ax)
    global simTimer isRunning isStopped isPaused statusLamp;
    global current_time reactor_data;
    
    if isRunning
        fprintf('Fluid simulation already running!\n');
        return;
    end
    
    % Set simulation state
    isRunning = true;
    isStopped = false;
    isPaused = false;
    current_time = 0;
    
    if ~isempty(reactor_data)
        reactor_data.time = 0;
        reactor_data.frame_count = 0;
        reactor_data.conversion = 0;
        reactor_data.start_time = tic;
    end
    
    % Update status
    if ~isempty(statusLamp) && isvalid(statusLamp)
        statusLamp.Color = [0.2, 1, 0.2];
    end
    
    % Create simulation timer
    dt = 0.015;
    update_interval = 0.03;
    
    simTimer = timer('ExecutionMode', 'fixedRate', ...
                     'Period', update_interval, ...
                     'TimerFcn', @(~,~) continuousFeedSimulationStep(ax, dt), ...
                     'StartDelay', 0.1);
    
    try
        start(simTimer);
        fprintf('Continuous feed pilot plant simulation started\n');
    catch ME
        fprintf('Error starting simulation timer: %s\n', ME.message);
        isRunning = false;
        if ~isempty(statusLamp) && isvalid(statusLamp)
            statusLamp.Color = [1, 0.2, 0.2];
        end
    end
end

% Main continuous feed simulation step
function continuousFeedSimulationStep(ax, dt)
    global current_time isRunning isPaused isStopped;
    global reactor_data impeller_speed continuous_feeding;
    
    if isStopped || ~isRunning || isPaused
        return;
    end
    
    try
        current_time = current_time + dt;
        
        if ~isempty(reactor_data)
            reactor_data.time = current_time;
            reactor_data.frame_count = reactor_data.frame_count + 1;
            reactor_data.impeller_speed = impeller_speed;
            
            % Update impeller rotations
            update_multiple_impeller_rotations();
            
            % Continuous feed system operations
            if continuous_feeding
                processContinuousFeed(dt);
                processParticleExit(dt);
            end
            
            % Update species physics
            update_chemical_species_physics_multiple_impellers();
            
            % Simulate reaction
            simulate_saponification_reaction();
            
            % Update flow visualization
            update_enhanced_flow_lines_multiple_impellers();
            
            % Update fluid surface
            update_fluid_surface_with_waves();
            
            % Update UI displays
            updateEnhancedUIDisplaysWithFeedSystem();
        end
        
        drawnow limitrate;
        
    catch ME
        fprintf('Error in continuous feed simulation step: %s\n', ME.message);
        stopSimulation();
    end
end

% Process continuous feed from inlet pipes
function processContinuousFeed(dt)
    global reactor_data feed_rate_factor;
    
    if isempty(reactor_data)
        return;
    end
    
    % Calculate number of particles to add this timestep
    base_feed_rate = reactor_data.feed_rate * feed_rate_factor;
    particles_to_add = poissrnd(base_feed_rate * dt);
    
    if particles_to_add > 0
        % Split between two inlets (EtOAc and NaOH)
        etac_particles = ceil(particles_to_add / 2);
        naoh_particles = floor(particles_to_add / 2);
        
        % FIXED: Ensure all arrays are initialized as column vectors
        if ~isfield(reactor_data.species.ethyl_acetate, 'age') || isempty(reactor_data.species.ethyl_acetate.age)
            reactor_data.species.ethyl_acetate.age = [];
        end
        if ~isfield(reactor_data.species.naoh, 'age') || isempty(reactor_data.species.naoh.age)
            reactor_data.species.naoh.age = [];
        end
        
        % Add EtOAc particles from Inlet 1
        if etac_particles > 0
            inlet1_pos = reactor_data.inlet_positions(1, :);
            
            % Pre-allocate arrays for new particles
            new_x = zeros(etac_particles, 1);  % FIXED: Force column vector
            new_y = zeros(etac_particles, 1);
            new_z = zeros(etac_particles, 1);
            new_vx = zeros(etac_particles, 1);
            new_vy = zeros(etac_particles, 1);
            new_vz = zeros(etac_particles, 1);
            new_age = zeros(etac_particles, 1);
            
            valid_particles = 0;
            for i = 1:etac_particles
                x_new = inlet1_pos(1) + (rand() - 0.5) * 0.1;
                y_new = inlet1_pos(2) + (rand() - 0.5) * 0.1;
                z_new = inlet1_pos(3) - rand() * 0.2;
                
                if sqrt(x_new^2 + y_new^2) < reactor_data.tank_radius * 0.95
                    valid_particles = valid_particles + 1;
                    new_x(valid_particles) = x_new;
                    new_y(valid_particles) = y_new;
                    new_z(valid_particles) = z_new;
                    new_vx(valid_particles) = (rand() - 0.5) * 0.05;
                    new_vy(valid_particles) = (rand() - 0.5) * 0.05;
                    new_vz(valid_particles) = -0.1 - rand() * 0.05;
                    new_age(valid_particles) = 0;
                end
            end
            
            % Trim arrays to actual valid particles
            if valid_particles > 0
                new_x = new_x(1:valid_particles);
                new_y = new_y(1:valid_particles);
                new_z = new_z(1:valid_particles);
                new_vx = new_vx(1:valid_particles);
                new_vy = new_vy(1:valid_particles);
                new_vz = new_vz(1:valid_particles);
                new_age = new_age(1:valid_particles);
                
                % FIXED: Ensure consistent dimensions before concatenation
                reactor_data.species.ethyl_acetate.x = ensureColumnVector(reactor_data.species.ethyl_acetate.x);
                reactor_data.species.ethyl_acetate.y = ensureColumnVector(reactor_data.species.ethyl_acetate.y);
                reactor_data.species.ethyl_acetate.z = ensureColumnVector(reactor_data.species.ethyl_acetate.z);
                reactor_data.species.ethyl_acetate.vx = ensureColumnVector(reactor_data.species.ethyl_acetate.vx);
                reactor_data.species.ethyl_acetate.vy = ensureColumnVector(reactor_data.species.ethyl_acetate.vy);
                reactor_data.species.ethyl_acetate.vz = ensureColumnVector(reactor_data.species.ethyl_acetate.vz);
                reactor_data.species.ethyl_acetate.age = ensureColumnVector(reactor_data.species.ethyl_acetate.age);
                
                % Concatenate new particles
                reactor_data.species.ethyl_acetate.x = [reactor_data.species.ethyl_acetate.x; new_x];
                reactor_data.species.ethyl_acetate.y = [reactor_data.species.ethyl_acetate.y; new_y];
                reactor_data.species.ethyl_acetate.z = [reactor_data.species.ethyl_acetate.z; new_z];
                reactor_data.species.ethyl_acetate.vx = [reactor_data.species.ethyl_acetate.vx; new_vx];
                reactor_data.species.ethyl_acetate.vy = [reactor_data.species.ethyl_acetate.vy; new_vy];
                reactor_data.species.ethyl_acetate.vz = [reactor_data.species.ethyl_acetate.vz; new_vz];
                reactor_data.species.ethyl_acetate.age = [reactor_data.species.ethyl_acetate.age; new_age];
            end
        end
        
        % Add NaOH particles from Inlet 2
        if naoh_particles > 0
            inlet2_pos = reactor_data.inlet_positions(2, :);
            
            % Pre-allocate arrays for new particles
            new_x = zeros(naoh_particles, 1);  % FIXED: Force column vector
            new_y = zeros(naoh_particles, 1);
            new_z = zeros(naoh_particles, 1);
            new_vx = zeros(naoh_particles, 1);
            new_vy = zeros(naoh_particles, 1);
            new_vz = zeros(naoh_particles, 1);
            new_age = zeros(naoh_particles, 1);
            
            valid_particles = 0;
            for i = 1:naoh_particles
                x_new = inlet2_pos(1) + (rand() - 0.5) * 0.1;
                y_new = inlet2_pos(2) + (rand() - 0.5) * 0.1;
                z_new = inlet2_pos(3) - rand() * 0.2;
                
                if sqrt(x_new^2 + y_new^2) < reactor_data.tank_radius * 0.95
                    valid_particles = valid_particles + 1;
                    new_x(valid_particles) = x_new;
                    new_y(valid_particles) = y_new;
                    new_z(valid_particles) = z_new;
                    new_vx(valid_particles) = (rand() - 0.5) * 0.05;
                    new_vy(valid_particles) = (rand() - 0.5) * 0.05;
                    new_vz(valid_particles) = -0.1 - rand() * 0.05;
                    new_age(valid_particles) = 0;
                end
            end
            
            % Trim arrays to actual valid particles
            if valid_particles > 0
                new_x = new_x(1:valid_particles);
                new_y = new_y(1:valid_particles);
                new_z = new_z(1:valid_particles);
                new_vx = new_vx(1:valid_particles);
                new_vy = new_vy(1:valid_particles);
                new_vz = new_vz(1:valid_particles);
                new_age = new_age(1:valid_particles);
                
                % FIXED: Ensure consistent dimensions before concatenation
                reactor_data.species.naoh.x = ensureColumnVector(reactor_data.species.naoh.x);
                reactor_data.species.naoh.y = ensureColumnVector(reactor_data.species.naoh.y);
                reactor_data.species.naoh.z = ensureColumnVector(reactor_data.species.naoh.z);
                reactor_data.species.naoh.vx = ensureColumnVector(reactor_data.species.naoh.vx);
                reactor_data.species.naoh.vy = ensureColumnVector(reactor_data.species.naoh.vy);
                reactor_data.species.naoh.vz = ensureColumnVector(reactor_data.species.naoh.vz);
                reactor_data.species.naoh.age = ensureColumnVector(reactor_data.species.naoh.age);
                
                % Concatenate new particles
                reactor_data.species.naoh.x = [reactor_data.species.naoh.x; new_x];
                reactor_data.species.naoh.y = [reactor_data.species.naoh.y; new_y];
                reactor_data.species.naoh.z = [reactor_data.species.naoh.z; new_z];
                reactor_data.species.naoh.vx = [reactor_data.species.naoh.vx; new_vx];
                reactor_data.species.naoh.vy = [reactor_data.species.naoh.vy; new_vy];
                reactor_data.species.naoh.vz = [reactor_data.species.naoh.vz; new_vz];
                reactor_data.species.naoh.age = [reactor_data.species.naoh.age; new_age];
            end
        end
        
        % Update feed counters
        reactor_data.total_fed = reactor_data.total_fed + particles_to_add;
        reactor_data.particles_in_system = reactor_data.particles_in_system + particles_to_add;
    end
end

% Process particle exit from outlet pipe
function processParticleExit(dt)
    global reactor_data feed_rate_factor;
    
    if isempty(reactor_data)
        return;
    end
    
    % Calculate number of particles to remove this timestep
    base_exit_rate = reactor_data.exit_rate * feed_rate_factor;
    particles_to_remove = poissrnd(base_exit_rate * dt);
    
    if particles_to_remove > 0
        outlet_pos = reactor_data.outlet_position;
        exit_radius = 0.3;  % Radius around outlet where particles can exit
        
        % Collect particles near outlet
        all_species = {'ethyl_acetate', 'naoh', 'products'};
        particles_removed = 0;
        
        for s = 1:length(all_species)
            species_name = all_species{s};
            species = reactor_data.species.(species_name);
            
            if isempty(species.x) || particles_removed >= particles_to_remove
                continue;
            end
            
            % FIXED: Ensure consistent dimensions for distance calculation
            species.x = ensureColumnVector(species.x);
            species.y = ensureColumnVector(species.y);
            species.z = ensureColumnVector(species.z);
            
            % Find particles near outlet
            distances = sqrt((species.x - outlet_pos(1)).^2 + ...
                           (species.y - outlet_pos(2)).^2 + ...
                           (species.z - outlet_pos(3)).^2);
            
            % Also consider particles in lower region of tank (gravity effect)
            low_particles = find(species.z < reactor_data.fluid_level * 0.3);
            near_outlet = find(distances < exit_radius);
            
            % Combine both criteria (near outlet OR in lower region)
            candidates = unique([low_particles; near_outlet]);
            
            if ~isempty(candidates)
                % Remove particles (prioritize those closest to outlet)
                [~, sort_idx] = sort(distances(candidates));
                sorted_candidates = candidates(sort_idx);
                
                to_remove = min(length(sorted_candidates), particles_to_remove - particles_removed);
                remove_indices = sorted_candidates(1:to_remove);
                
                % FIXED: Remove particles using logical indexing (safer)
                keep_mask = true(length(species.x), 1);
                keep_mask(remove_indices) = false;
                
                species.x = species.x(keep_mask);
                species.y = species.y(keep_mask);
                species.z = species.z(keep_mask);
                species.vx = ensureColumnVector(species.vx);
                species.vy = ensureColumnVector(species.vy);
                species.vz = ensureColumnVector(species.vz);
                species.vx = species.vx(keep_mask);
                species.vy = species.vy(keep_mask);
                species.vz = species.vz(keep_mask);
                
                if isfield(species, 'age')
                    species.age = ensureColumnVector(species.age);
                    species.age = species.age(keep_mask);
                end
                
                % Update species data
                reactor_data.species.(species_name) = species;
                
                particles_removed = particles_removed + to_remove;
            end
        end
        
        % Update exit counters
        reactor_data.total_exited = reactor_data.total_exited + particles_removed;
        reactor_data.particles_in_system = reactor_data.particles_in_system - particles_removed;
        
        % Ensure particle count doesn't go negative
        reactor_data.particles_in_system = max(0, reactor_data.particles_in_system);
    end
end

% Update multiple impeller rotations
function update_multiple_impeller_rotations()
    global reactor_data;
    
    if isempty(reactor_data) || isempty(reactor_data.impeller_blades)
        return;
    end
    
    % Calculate angular velocity ONLY if running and speed > 0
    if reactor_data.impeller_running && reactor_data.impeller_speed > 1.0
        angular_velocity = reactor_data.impeller_speed * 2 * pi / 60;
        
        % Update all impellers
        for imp = 1:reactor_data.impeller_count
            phase_shift = (imp - 1) * pi / 3;
            reactor_data.impeller_angles(imp) = reactor_data.impeller_angles(imp) + angular_velocity * reactor_data.dt;
            
            % Update blades
            for blade = 1:reactor_data.blade_count
                if ~isempty(reactor_data.blade_original_vertices{imp, blade}) && ...
                   isvalid(reactor_data.impeller_blades{imp}(blade))
                    
                    original_data = reactor_data.blade_original_vertices{imp, blade};
                    new_angle = original_data.theta + reactor_data.impeller_angles(imp) + phase_shift;
                    
                    new_x = original_data.r .* cos(new_angle);
                    new_y = original_data.r .* sin(new_angle);
                    
                    try
                        set(reactor_data.impeller_blades{imp}(blade), 'XData', new_x, 'YData', new_y);
                    catch ME
                        fprintf('Warning: Could not update impeller %d blade %d: %s\n', imp, blade, ME.message);
                    end
                end
            end
        end
    end
    % If not running or RPM is 0, impellers stay in current position (no rotation)
end

% Update chemical species physics with enhanced continuous flow effects
function update_chemical_species_physics_multiple_impellers()
    global reactor_data current_time fluid_viscosity_factor;
    
    if isempty(reactor_data) || isempty(reactor_data.species)
        return;
    end
    
    species_names = {'ethyl_acetate', 'naoh', 'products'};
    
    for s = 1:length(species_names)
        species_name = species_names{s};
        species = reactor_data.species.(species_name);
        
        if isempty(species.x), continue; end
        
        n_particles = length(species.x);
        
                  % Update particle ages
        if isfield(species, 'age')
            species.age = species.age + reactor_data.dt;
        else
            species.age = zeros(size(species.x));  % Initialize if missing
            fprintf('Initializing age field for %s (%d particles)\n', species_name, n_particles);
        end
        
        for i = 1:n_particles
            x = species.x(i);
            y = species.y(i);
            z = species.z(i);
            
            % Initialize force accumulation
            total_vx = 0;
            total_vy = 0;
            total_vz = 0;
            
            % ONLY APPLY IMPELLER FORCES IF RUNNING AND RPM > 0
            if reactor_data.impeller_running && reactor_data.impeller_speed > 1.0
                for imp = 1:reactor_data.impeller_count
                    impeller_height = reactor_data.impeller_heights(imp);
                    
                    % Distance from impeller
                    dx = x;
                    dy = y;
                    dz = z - impeller_height;
                    dist_from_impeller = sqrt(dx^2 + dy^2 + dz^2);
                    
                    % Impeller influence
                    influence = max(0, 1 - dist_from_impeller / 4.0);
                    
                    if influence > 0
                        radius = sqrt(x^2 + y^2);
                        if radius > 0.05
                            angle = atan2(y, x);
                            
                            % Tangential velocity proportional to RPM
                            tip_speed = reactor_data.impeller_speed / 60 * reactor_data.impeller_diameter * pi / 2;
                            tangential_speed = tip_speed * influence * (radius / (reactor_data.impeller_diameter/2));
                            
                            % Radial flow
                            impeller_angle = reactor_data.impeller_angles(imp);
                            radial_component = tangential_speed * 0.3 * cos(impeller_angle + angle * reactor_data.blade_count);
                            
                            % Blade effects
                            blade_effect = 1 + 0.5 * sin(impeller_angle * reactor_data.blade_count + angle * reactor_data.blade_count);
                            
                            % Scale all forces by RPM (so at 0 RPM = 0 force)
                            rpm_factor = reactor_data.impeller_speed / 600;  % Scale 0-1
                            
                            % Accumulate forces
                            total_vx = total_vx + (-sin(angle) * tangential_speed + cos(angle) * radial_component) * blade_effect * influence * 0.15 * rpm_factor;
                            total_vy = total_vy + (cos(angle) * tangential_speed + sin(angle) * radial_component) * blade_effect * influence * 0.15 * rpm_factor;
                            
                            % Vertical mixing - only when running
                            turbulent_intensity = rpm_factor;
                            vertical_pumping = sin(current_time * 8 + angle * 2) * influence * turbulent_intensity * 0.4;
                            total_vz = total_vz + vertical_pumping * 0.1 * rpm_factor;
                            
                            % Turbulence - only when running fast enough
                            if reactor_data.impeller_speed > 50
                                turbulent_vx = (rand() - 0.5) * turbulent_intensity * 0.08 * influence;
                                turbulent_vy = (rand() - 0.5) * turbulent_intensity * 0.08 * influence;
                                turbulent_vz = (rand() - 0.5) * turbulent_intensity * 0.04 * influence;
                                
                                total_vx = total_vx + turbulent_vx * rpm_factor;
                                total_vy = total_vy + turbulent_vy * rpm_factor;
                                total_vz = total_vz + turbulent_vz * rpm_factor;
                            end
                        end
                    end
                end
            end
            
            % BUOYANCY FORCE - Particles should be neutrally buoyant in fluid
            % Instead of strong gravity, add small buoyancy to keep particles floating
            buoyancy_force = 0.001;  % Very small upward force
            particle_density_factor = 0.98;  % Slightly less dense than fluid
            
            if z < reactor_data.fluid_level * 0.3
                % Add stronger buoyancy near bottom to prevent settling
                total_vz = total_vz + buoyancy_force * (1 - particle_density_factor) * 2;
            else
                % Small buoyancy throughout fluid
                total_vz = total_vz + buoyancy_force * (1 - particle_density_factor);
            end
            
            % Very weak gravity for realism, but not enough to overcome buoyancy
            weak_gravity = -0.002;  % Much weaker than before
            total_vz = total_vz + weak_gravity;
            
            % Add flow towards outlet only if particles are very close to bottom
            if z < reactor_data.fluid_level * 0.15
                outlet_pos = reactor_data.outlet_position;
                dx_outlet = outlet_pos(1) - x;
                dy_outlet = outlet_pos(2) - y;
                dz_outlet = outlet_pos(3) - z;
                dist_to_outlet = sqrt(dx_outlet^2 + dy_outlet^2 + dz_outlet^2);
                
                if dist_to_outlet > 0
                    drainage_strength = 0.01 / (1 + dist_to_outlet);
                    total_vx = total_vx + dx_outlet * drainage_strength;
                    total_vy = total_vy + dy_outlet * drainage_strength;
                    total_vz = total_vz + dz_outlet * drainage_strength;
                end
            end
            
            % Apply forces to velocity
            species.vx(i) = species.vx(i) + total_vx;
            species.vy(i) = species.vy(i) + total_vy;
            species.vz(i) = species.vz(i) + total_vz;
            
            % Enhanced viscosity damping
            base_damping = 0.98;  % Base damping when running
            if reactor_data.impeller_speed < 1.0
                % Strong damping when not running - particles should stop
                base_damping = 0.85;  
            end
            
            damping = base_damping - reactor_data.viscosity * 0.01 * fluid_viscosity_factor;
            damping = max(damping, 0.7);  % Minimum damping
            
            species.vx(i) = species.vx(i) * damping;
            species.vy(i) = species.vy(i) * damping;
            species.vz(i) = species.vz(i) * damping;
            
            % Update positions
            species.x(i) = species.x(i) + species.vx(i) * reactor_data.dt;
            species.y(i) = species.y(i) + species.vy(i) * reactor_data.dt;
            species.z(i) = species.z(i) + species.vz(i) * reactor_data.dt;
            
            % Boundary conditions
            r = sqrt(species.x(i)^2 + species.y(i)^2);
            if r > reactor_data.tank_radius * 0.95
                norm_x = species.x(i) / r;
                norm_y = species.y(i) / r;
                
                species.x(i) = norm_x * reactor_data.tank_radius * 0.95;
                species.y(i) = norm_y * reactor_data.tank_radius * 0.95;
                
                species.vx(i) = species.vx(i) * -0.3;
                species.vy(i) = species.vy(i) * -0.3;
            end
            
            % Vertical boundaries - particles should float in fluid
            if species.z(i) < 0.1
                species.z(i) = 0.1;
                species.vz(i) = abs(species.vz(i)) * 0.5;  % Bounce up gently
            elseif species.z(i) > reactor_data.fluid_level * 0.98
                species.z(i) = reactor_data.fluid_level * 0.98;
                species.vz(i) = -abs(species.vz(i)) * 0.3;  % Gentle bounce down
            end
        end
             
       
        % Update species data
        reactor_data.species.(species_name) = species;
    end
       
    % Update visualizations
    update_species_visualizations();
end

% Enhanced saponification reaction with better product creation
function simulate_saponification_reaction()
    global reactor_data residence_time;
    
    if isempty(reactor_data) || isempty(reactor_data.species)
        return;
    end
    
    % Calculate mixing efficiency based on RPM
    rpm = reactor_data.impeller_speed;
    max_rpm = 600;
    
    if rpm < 10
        mixing_efficiency = 0.1;
    else
        mixing_efficiency = 0.1 + 2.9 * (log10(rpm) - 1) / (log10(max_rpm) - 1);
        mixing_efficiency = max(0.1, min(3.0, mixing_efficiency));
    end
    
    effective_residence_time = residence_time / mixing_efficiency;
    
    % Convert aged particles to products
    species_names = {'ethyl_acetate', 'naoh'};
    max_conversions_per_frame = round(5 * mixing_efficiency);
    total_converted = 0;
    
    for s = 1:length(species_names)
        species_name = species_names{s};
        species = reactor_data.species.(species_name);
        
        if isempty(species.x), continue; end
        
        % FIXED: Ensure all arrays are column vectors
        species.x = ensureColumnVector(species.x);
        species.y = ensureColumnVector(species.y);
        species.z = ensureColumnVector(species.z);
        species.vx = ensureColumnVector(species.vx);
        species.vy = ensureColumnVector(species.vy);
        species.vz = ensureColumnVector(species.vz);
        
        % Initialize age field if missing
        if ~isfield(species, 'age') || isempty(species.age)
            species.age = zeros(size(species.x));
        else
            species.age = ensureColumnVector(species.age);
        end
        
        % Find particles that have reached effective residence time
        aged_particles = find(species.age >= effective_residence_time);
        
        if ~isempty(aged_particles)
            % Convert particles based on mixing efficiency
            num_to_convert = min(length(aged_particles), max_conversions_per_frame - total_converted);
            particles_to_convert = aged_particles(1:num_to_convert);
            
            % FIXED: Ensure products arrays are properly initialized as column vectors
            if ~isfield(reactor_data.species.products, 'age') || isempty(reactor_data.species.products.age)
                reactor_data.species.products.age = [];
            end
            
            % Prepare new product particles
            new_product_x = species.x(particles_to_convert);
            new_product_y = species.y(particles_to_convert);
            new_product_z = species.z(particles_to_convert);
            new_product_vx = species.vx(particles_to_convert);
            new_product_vy = species.vy(particles_to_convert);
            new_product_vz = species.vz(particles_to_convert);
            new_product_age = zeros(length(particles_to_convert), 1);
            
            % FIXED: Ensure products arrays are column vectors before concatenation
            reactor_data.species.products.x = ensureColumnVector(reactor_data.species.products.x);
            reactor_data.species.products.y = ensureColumnVector(reactor_data.species.products.y);
            reactor_data.species.products.z = ensureColumnVector(reactor_data.species.products.z);
            reactor_data.species.products.vx = ensureColumnVector(reactor_data.species.products.vx);
            reactor_data.species.products.vy = ensureColumnVector(reactor_data.species.products.vy);
            reactor_data.species.products.vz = ensureColumnVector(reactor_data.species.products.vz);
            reactor_data.species.products.age = ensureColumnVector(reactor_data.species.products.age);
            
            % Add to products (GREEN)
            reactor_data.species.products.x = [reactor_data.species.products.x; new_product_x];
            reactor_data.species.products.y = [reactor_data.species.products.y; new_product_y];
            reactor_data.species.products.z = [reactor_data.species.products.z; new_product_z];
            reactor_data.species.products.vx = [reactor_data.species.products.vx; new_product_vx];
            reactor_data.species.products.vy = [reactor_data.species.products.vy; new_product_vy];
            reactor_data.species.products.vz = [reactor_data.species.products.vz; new_product_vz];
            reactor_data.species.products.age = [reactor_data.species.products.age; new_product_age];
            
            % FIXED: Remove converted particles using logical indexing
            keep_mask = true(length(species.x), 1);
            keep_mask(particles_to_convert) = false;
            
            species.x = species.x(keep_mask);
            species.y = species.y(keep_mask);
            species.z = species.z(keep_mask);
            species.vx = species.vx(keep_mask);
            species.vy = species.vy(keep_mask);
            species.vz = species.vz(keep_mask);
            species.age = species.age(keep_mask);
            
            % Update species data
            reactor_data.species.(species_name) = species;
            total_converted = total_converted + length(particles_to_convert);
            
            % Print conversion info
            if mod(reactor_data.frame_count, 120) == 0
                fprintf('✓ Converted %d %s particles to GREEN! (RPM: %.0f, Eff Time: %.1fs)\n', ...
                        length(particles_to_convert), species_name, rpm, effective_residence_time);
            end
        end
        
        if total_converted >= max_conversions_per_frame
            break;
        end
    end
    
    % Store mixing efficiency for display
    reactor_data.mixing_efficiency = mixing_efficiency;
    
    % Calculate conversion
    total_particles = length(reactor_data.species.ethyl_acetate.x) + ...
                     length(reactor_data.species.naoh.x) + ...
                     length(reactor_data.species.products.x);
    
    if total_particles > 0
        reactor_data.conversion = (length(reactor_data.species.products.x) / total_particles) * 100;
    else
        reactor_data.conversion = 0;
    end
end

% Update species visualizations
function update_species_visualizations()
    global reactor_data show_particles;
    
    if isempty(reactor_data) || isempty(reactor_data.species) || ~show_particles
        return;
    end
    
    % Update Ethyl Acetate particles
    if isfield(reactor_data.species_plots, 'ethyl_acetate') && isvalid(reactor_data.species_plots.ethyl_acetate)
        if ~isempty(reactor_data.species.ethyl_acetate.x)
            set(reactor_data.species_plots.ethyl_acetate, ...
                'XData', reactor_data.species.ethyl_acetate.x, ...
                'YData', reactor_data.species.ethyl_acetate.y, ...
                'ZData', reactor_data.species.ethyl_acetate.z, ...
                'CData', [1.0, 0.8, 0.2], ...  % Ensure yellow color
                'Visible', 'on');
        else
            set(reactor_data.species_plots.ethyl_acetate, ...
                'XData', [], 'YData', [], 'ZData', [], 'Visible', 'off');
        end
    end
    
    % Update NaOH particles
    if isfield(reactor_data.species_plots, 'naoh') && isvalid(reactor_data.species_plots.naoh)
        if ~isempty(reactor_data.species.naoh.x)
            set(reactor_data.species_plots.naoh, ...
                'XData', reactor_data.species.naoh.x, ...
                'YData', reactor_data.species.naoh.y, ...
                'ZData', reactor_data.species.naoh.z, ...
                'CData', [0.2, 0.4, 1.0], ...  % Ensure blue color
                'Visible', 'on');
        else
            set(reactor_data.species_plots.naoh, ...
                'XData', [], 'YData', [], 'ZData', [], 'Visible', 'off');
        end
    end
    
    % Update Product particles - ENHANCED WITH DEBUG INFO
    if isfield(reactor_data.species_plots, 'products') && isvalid(reactor_data.species_plots.products)
        if ~isempty(reactor_data.species.products.x)
            num_products = length(reactor_data.species.products.x);
            fprintf('Updating %d GREEN product particles\n', num_products);
            
            set(reactor_data.species_plots.products, ...
                'XData', reactor_data.species.products.x, ...
                'YData', reactor_data.species.products.y, ...
                'ZData', reactor_data.species.products.z, ...
                'CData', [0.2, 0.8, 0.3], ...  % Force GREEN color
                'SizeData', 35, ...             % Slightly larger for visibility
                'MarkerFaceAlpha', 0.9, ...    % Make more opaque
                'Visible', 'on');
        else
            set(reactor_data.species_plots.products, ...
                'XData', [], 'YData', [], 'ZData', [], 'Visible', 'off');
        end
    else
        fprintf('Warning: Products plot not found or invalid\n');
    end
end

% Update enhanced flow lines
function update_enhanced_flow_lines_multiple_impellers()
    global reactor_data show_flow_lines current_time dark_mode;
    
    if isempty(reactor_data) || ~show_flow_lines || isempty(reactor_data.flow_lines)
        return;
    end
    
    if dark_mode
        flow_colors = {[0.3, 0.9, 1.0], [0.9, 0.3, 1.0], [1.0, 0.9, 0.3]};
        alpha_value = 0.8;
    else
        flow_colors = {[0, 0.8, 1], [0.8, 0, 1], [1, 0.8, 0]};
        alpha_value = 0.7;
    end
    
    % Update flow lines with time-dependent patterns
    flow_line_count = 20;
    line_idx = 1;
    
    for imp = 1:reactor_data.impeller_count
        impeller_height = reactor_data.impeller_heights(imp);
        flow_color = flow_colors{min(imp, 3)};
        
        for i = 1:flow_line_count
            if line_idx <= length(reactor_data.flow_lines) && isvalid(reactor_data.flow_lines(line_idx))
                angle_offset = (i-1) * 2*pi / flow_line_count;
                
                % Create time-varying flow patterns
                t = linspace(0, 3*pi, 60);
                time_phase = current_time * reactor_data.impeller_speed / 60 * 0.5;
                impeller_phase = reactor_data.impeller_angles(imp);
                
                % Enhanced turbulent flow
                turbulence_intensity = reactor_data.impeller_speed / 600;
                radius_variation = 0.8 + 0.3 * sin(t * 2 + time_phase + impeller_phase) * turbulence_intensity;
                base_radius = reactor_data.impeller_diameter * 1.5;
                
                x_flow = base_radius * radius_variation .* cos(angle_offset + t + time_phase + impeller_phase) .* exp(-t*0.05);
                y_flow = base_radius * radius_variation .* sin(angle_offset + t + time_phase + impeller_phase) .* exp(-t*0.05);
                z_flow = impeller_height + 0.4 * sin(t + time_phase) .* exp(-t*0.1);
                
                % Apply wall reflection
                for j = 1:length(x_flow)
                    r = sqrt(x_flow(j)^2 + y_flow(j)^2);
                    if r > reactor_data.tank_radius * 0.9
                        scale_factor = reactor_data.tank_radius * 0.9 / r;
                        x_flow(j) = x_flow(j) * scale_factor;
                        y_flow(j) = y_flow(j) * scale_factor;
                    end
                end
                
                try
                    set(reactor_data.flow_lines(line_idx), 'XData', x_flow, 'YData', y_flow, 'ZData', z_flow, ...
                        'Color', flow_color, 'LineWidth', 2 + turbulence_intensity);
                catch ME
                    fprintf('Warning: Could not update flow line %d: %s\n', line_idx, ME.message);
                end
            end
            line_idx = line_idx + 1;
        end
    end
end

% Update fluid surface with waves
function update_fluid_surface_with_waves()
    global reactor_data current_time dark_mode;
    
    if isempty(reactor_data) || ~isfield(reactor_data, 'fluid_surface') || ~isvalid(reactor_data.fluid_surface)
        return;
    end
    
    if ~reactor_data.show_waves
        return;
    end
    
    % Create wave patterns based on impeller speed and continuous flow
    wave_amplitude = reactor_data.impeller_speed / 600 * 0.05;
    wave_frequency = reactor_data.impeller_speed / 60 * 2;
    
    % Multiple impellers create more complex waves
    impeller_boost = 1 + (reactor_data.impeller_count - 1) * 0.5;
    wave_amplitude = wave_amplitude * impeller_boost;
    
    % Continuous feeding creates additional surface disturbance
    feed_disturbance = 1.2;
    wave_amplitude = wave_amplitude * feed_disturbance;
    
    % Get surface geometry
    surface_data = get(reactor_data.fluid_surface, 'ZData');
    base_height = reactor_data.fluid_level;
    
    [rows, cols] = size(surface_data);
    [X, Y] = meshgrid(linspace(-reactor_data.tank_radius*0.98, reactor_data.tank_radius*0.98, cols), ...
                      linspace(-reactor_data.tank_radius*0.98, reactor_data.tank_radius*0.98, rows));
    
    % Create wave pattern following particle circulation
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X);
    
    % === KEY MODIFICATION: CIRCULATION-FOLLOWING WAVES ===
    
    % Calculate the tangential flow direction (FORCED ANTICLOCKWISE)
    % Always rotate counterclockwise regardless of impeller direction
    rotation_direction = 1; % FIXED: Always +1 for anticlockwise circulation
    
    % Angular velocity of fluid circulation (always positive for anticlockwise)
    angular_velocity = abs(reactor_data.impeller_speed) * 2 * pi / 60; % rad/s
    
    % Combine waves from all impellers with circulation consideration
    total_waves = zeros(size(R));
    for imp = 1:reactor_data.impeller_count
        impeller_height = reactor_data.impeller_heights(imp);
        height_factor = exp(-abs(base_height - impeller_height) / reactor_data.fluid_level);
        
        % === CIRCULATION-FOLLOWING RADIAL WAVES ===
        % Waves that propagate outward but are "bent" by anticlockwise circulation
        circulation_effect = angular_velocity * current_time; % Always anticlockwise
        
        % Modified radial waves that follow circulation
        % FIXED: Use element-wise multiplication
        radial_phase = R .* 3 - wave_frequency * current_time;
        circulation_phase = Theta + circulation_effect * 0.3; % Circulation bends the wave fronts
        
        radial_waves = wave_amplitude * height_factor .* ...
                      sin(radial_phase + reactor_data.impeller_angles(imp) + circulation_phase .* 0.5) .* ...
                      exp(-R ./ reactor_data.tank_radius);
        
        % === CIRCULATION-FOLLOWING CIRCUMFERENTIAL WAVES ===
        % Waves that travel in the tangential direction (ANTICLOCKWISE)
        tangential_wave_speed = angular_velocity * 0.8; % Waves travel slightly slower than fluid
        tangential_phase = (Theta .* reactor_data.blade_count - tangential_wave_speed * current_time); % Always anticlockwise
        
        % Add radial variation to make waves more realistic
        % FIXED: Use element-wise multiplication
        radial_modulation = 1 + 0.3 * sin(R .* 4 + current_time * 2);
        
        circumferential_waves = wave_amplitude * height_factor * 0.7 .* ...
                               sin(tangential_phase + reactor_data.impeller_angles(imp)) .* ...
                               radial_modulation .* exp(-R ./ reactor_data.tank_radius);
        
        % === SPIRAL WAVES (New Feature) ===
        % Waves that spiral outward following the anticlockwise circulation pattern
        spiral_pitch = 2.0; % Controls how tight the spiral is
        spiral_phase = (Theta - R .* spiral_pitch - tangential_wave_speed * current_time * 0.5); % Always anticlockwise
        
        spiral_waves = wave_amplitude * height_factor * 0.4 .* ...
                      sin(spiral_phase + reactor_data.impeller_angles(imp) * 0.5) .* ...
                      exp(-R ./ reactor_data.tank_radius) .* ...
                      (R ./ reactor_data.tank_radius); % Stronger at outer radius
        
        % === BLADE PASSAGE WAVES ===
        % Discrete waves created by each blade passing
        blade_frequency = reactor_data.blade_count * angular_velocity;
        blade_waves = zeros(size(R));
        
        for blade = 1:reactor_data.blade_count
            blade_angle = (blade - 1) * 2*pi / reactor_data.blade_count + reactor_data.impeller_angles(imp);
            
            % Angular distance from current blade position
            angular_distance = mod(Theta - blade_angle + pi, 2*pi) - pi;
            
            % Create wake pattern behind each blade
            wake_strength = exp(-abs(angular_distance) ./ 0.5) .* exp(-abs(R - reactor_data.impeller_diameter/4) ./ 0.3);
            % FIXED: Use element-wise multiplication
            wake_phase = -blade_frequency * current_time + R .* 2;
            
            blade_waves = blade_waves + wave_amplitude * height_factor * 0.3 .* ...
                         sin(wake_phase) .* wake_strength;
        end
        
        % Combine all wave types
        total_waves = total_waves + radial_waves + circumferential_waves + spiral_waves + blade_waves;
    end
    
    % === INLET JET WAVES (Enhanced) ===
    % Add waves from inlet jets that also follow anticlockwise circulation
    for i = 1:reactor_data.inlet_pipe_count
        inlet_pos = reactor_data.inlet_positions(i, :);
        inlet_distance = sqrt((X - inlet_pos(1)).^2 + (Y - inlet_pos(2)).^2);
        
        % Jet waves that get caught up in anticlockwise circulation
        jet_angle = atan2(inlet_pos(2), inlet_pos(1));
        circulation_effect = angular_velocity * current_time; % Always anticlockwise
        jet_circulation_phase = Theta - jet_angle + circulation_effect * 0.2;
        
        inlet_waves = wave_amplitude * 0.4 .* ...
                     sin(wave_frequency * current_time * 3 + jet_circulation_phase) .* ...
                     exp(-inlet_distance ./ 0.4);
        total_waves = total_waves + inlet_waves;
    end
    
    % Surface tension damping (enhanced)
    edge_distance = reactor_data.tank_radius - R;
    damping_factor = 1 - exp(-edge_distance ./ 0.1);
    
    % Add circulation-based damping near walls (always anticlockwise)
    wall_circulation_damping = 1 - 0.3 .* exp(-edge_distance ./ 0.05) .* (angular_velocity / 10);
    damping_factor = damping_factor .* wall_circulation_damping;
    
    % Combine effects
    wave_height = base_height + total_waves .* damping_factor;
    
    % Apply constraints
    wave_height = max(wave_height, base_height - wave_amplitude * 2);
    wave_height = min(wave_height, base_height + wave_amplitude * 2);
    
    try
        set(reactor_data.fluid_surface, 'ZData', wave_height);
        
        % Update surface color based on circulation intensity
        if dark_mode
            base_color = [0.2, 0.4, 0.7];
        else
            base_color = [0.3, 0.5, 0.8];
        end
        
        % Color varies with circulation strength
        circulation_intensity = abs(angular_velocity) / 60 * impeller_boost; % 0-1 scale
        surface_color = base_color + [0.1, 0.1, 0.2] * circulation_intensity;
        surface_color = min(surface_color, [1, 1, 1]);
        
        set(reactor_data.fluid_surface, 'FaceColor', surface_color);
    catch ME
        fprintf('Warning: Could not update fluid surface: %s\n', ME.message);
    end
end

% === OPTIONAL: Add wave visualization function ===
function visualize_circulation_vectors(ax)
    global reactor_data;
    
    % Add arrows showing circulation direction on surface
    if isempty(reactor_data)
        return;
    end
    
    % Create circulation vector field on surface
    [X_vec, Y_vec] = meshgrid(linspace(-reactor_data.tank_radius*0.8, reactor_data.tank_radius*0.8, 12), ...
                              linspace(-reactor_data.tank_radius*0.8, reactor_data.tank_radius*0.8, 12));
    
    R_vec = sqrt(X_vec.^2 + Y_vec.^2);
    Theta_vec = atan2(Y_vec, X_vec);
    
    % Calculate circulation velocity at each point
    angular_velocity = reactor_data.impeller_speed * 2 * pi / 60;
    tangential_velocity = angular_velocity * R_vec * 0.3; % Scale down for visualization
    
    % Convert to Cartesian components
    U_circulation = -tangential_velocity .* sin(Theta_vec);
    V_circulation = tangential_velocity .* cos(Theta_vec);
    
    % Only show vectors inside tank
    mask = R_vec < reactor_data.tank_radius * 0.9;
    
    % Plot circulation vectors on fluid surface
    Z_surface = reactor_data.fluid_level * ones(size(X_vec));
    
    quiver3(ax, X_vec(mask), Y_vec(mask), Z_surface(mask), ...
            U_circulation(mask), V_circulation(mask), zeros(sum(mask(:)), 1), ...
            'Color', [1, 1, 0], 'LineWidth', 1.5, 'MaxHeadSize', 0.3, ...
            'AutoScale', 'on', 'AutoScaleFactor', 0.8);
end


% Update enhanced UI displays with feed system information
function updateEnhancedUIDisplaysWithFeedSystem()
    global reactor_data current_time;
    global rpmGauge rpmValueLabel statusLamp;
    global performance_metrics frame_count last_fps_time;
    global allUIComponents;
    
      % Update RPM gauge and display
    if ~isempty(rpmGauge) && isvalid(rpmGauge) && ~isempty(reactor_data)
        try
            rpmGauge.Value = reactor_data.impeller_speed;
        catch
        end
    end
    
if ~isempty(rpmValueLabel) && isvalid(rpmValueLabel) && ~isempty(reactor_data)
    try
        % Show RPM with mixing efficiency
        mixing_eff = 1.0;
        if isfield(reactor_data, 'mixing_efficiency')
            mixing_eff = reactor_data.mixing_efficiency;
        end
        rpmValueLabel.Text = sprintf('%.0f RPM (%.1fx)', reactor_data.impeller_speed, mixing_eff);
    catch
        % Ignore errors
    end
end
    
    % Update feed system status with RPM effect
    if isfield(allUIComponents, 'feedStatusLabel') && isvalid(allUIComponents.feedStatusLabel)
        try
            colors = getCurrentColorScheme();
            if ~isempty(reactor_data)
                mixing_eff = 1.0;
                if isfield(reactor_data, 'mixing_efficiency')
                    mixing_eff = reactor_data.mixing_efficiency;
                end
                
                status_text = sprintf('RPM: %.0f (%.1fx) | Conv: %.1f%% | Holdup: %d', ...
                                    reactor_data.impeller_speed, mixing_eff, ...
                                    reactor_data.conversion, reactor_data.particles_in_system);
                allUIComponents.feedStatusLabel.Text = status_text;
                
                % Color based on mixing efficiency
                if mixing_eff > 2.0
                    allUIComponents.feedStatusLabel.FontColor = colors.success_color;  % Green - excellent mixing
                elseif mixing_eff > 1.5
                    allUIComponents.feedStatusLabel.FontColor = colors.accent_color;   % Blue - good mixing
                elseif mixing_eff > 1.0
                    allUIComponents.feedStatusLabel.FontColor = colors.warning_color;  % Orange - fair mixing
                else
                    allUIComponents.feedStatusLabel.FontColor = colors.error_color;    % Red - poor mixing
                end
            end
        catch
        end
    end
    
    % Update feed system status
if isfield(allUIComponents, 'feedStatusLabel') && isvalid(allUIComponents.feedStatusLabel)
    try
        colors = getCurrentColorScheme();
        if ~isempty(reactor_data)
            mixing_eff = 1.0;
            if isfield(reactor_data, 'mixing_efficiency')
                mixing_eff = reactor_data.mixing_efficiency;
            end
            
            status_text = sprintf('RPM: %.0f (%.1fx) | Conv: %.1f%% | Holdup: %d', ...
                                reactor_data.impeller_speed, mixing_eff, ...
                                reactor_data.conversion, reactor_data.particles_in_system);
            allUIComponents.feedStatusLabel.Text = status_text;
            
            % Color based on mixing efficiency
            if mixing_eff > 2.0
                allUIComponents.feedStatusLabel.FontColor = colors.success_color;  % Green
            elseif mixing_eff > 1.5
                allUIComponents.feedStatusLabel.FontColor = colors.accent_color;   % Blue
            elseif mixing_eff > 1.0
                allUIComponents.feedStatusLabel.FontColor = colors.warning_color;  % Orange
            else
                allUIComponents.feedStatusLabel.FontColor = colors.error_color;    % Red
            end
        end
    catch
        % Ignore errors
    end
end
    
    % Calculate performance metrics
    frame_count = frame_count + 1;
    current_fps_time = toc(reactor_data.start_time);
    
    if current_fps_time - last_fps_time > 1.0  % Update every second
        fps = frame_count / current_fps_time;
        performance_metrics.fps = fps;
        
        if ~isempty(reactor_data)
            performance_metrics.conversion = reactor_data.conversion;
            performance_metrics.mixing_efficiency = reactor_data.mixing_efficiency;
        end
        
        last_fps_time = current_fps_time;
        
         end
end

% === CALLBACK FUNCTIONS ===

% Switch value changed
function switchValueChanged(src, ~)
    if strcmp(src.Value, 'On')
        enableAllControls();
        fprintf('Main power ON - Continuous feed system ready\n');
    else
        disableAllControls();
        stopSimulation();
        fprintf('Main power OFF - Feed system stopped\n');
    end
end

% Knob value changed
function knobValueChanged(src, ~)
    global impeller_speed rpmGauge rpmValueLabel reactor_data;
    
    impeller_speed = src.Value;
    
    if ~isempty(reactor_data)
        reactor_data.impeller_speed = impeller_speed;
        
        % Set impeller running status based on RPM
        if impeller_speed < 1.0
            reactor_data.impeller_running = false;
            fprintf('Impellers STOPPED (RPM: %.0f)\n', impeller_speed);
        else
            reactor_data.impeller_running = true;
           % fprintf('Impellers RUNNING (RPM: %.0f)\n', impeller_speed);
        end
    end
    
    % Update displays
    if ~isempty(rpmGauge) && isvalid(rpmGauge)
        rpmGauge.Value = impeller_speed;
    end
    
    if ~isempty(rpmValueLabel) && isvalid(rpmValueLabel)
        rpmValueLabel.Text = sprintf('%.0f RPM', impeller_speed);
    end
end

% Feed toggle callback
function feedToggleCallback(src, ~)
    global continuous_feeding;
    
    continuous_feeding = src.Value;
    
    if continuous_feeding
        fprintf('Continuous feeding ENABLED\n');
    else
        fprintf('Continuous feeding DISABLED\n');
    end
end


% Residence time slider callback
function residenceTimeSliderCallback(src, ~)
    global residence_time residenceTimeLabel;
    
    residence_time = src.Value;
    residenceTimeLabel.Text = sprintf('%.0fs', residence_time);
    
    fprintf('Residence time set to %.0f seconds - Particles will turn GREEN after %.0fs\n', ...
            residence_time, residence_time);
end


% Target conversion slider callback
function targetConversionSliderCallback(src, ~)
    global target_conversion targetConversionLabel reactor_data;
    target_conversion = src.Value;
    targetConversionLabel.Text = sprintf('%.0f%%', target_conversion);
    
    if ~isempty(reactor_data)
        fprintf('Target conversion set to %.0f%%\n', target_conversion);
    end
end

% Feed rate slider callback
function feedRateSliderCallback(src, ~)
    global feed_rate_factor feedRateLabel reactor_data;
    
    feed_rate_factor = src.Value;
    feedRateLabel.Text = sprintf('%.2fx', feed_rate_factor);
    
    % Update feed and exit rates
    if ~isempty(reactor_data)
        base_feed_rate = 10.0;
        reactor_data.feed_rate = base_feed_rate * feed_rate_factor;
        reactor_data.exit_rate = base_feed_rate * 0.8 * feed_rate_factor;
      %  fprintf('Feed rate set to %.2fx (Feed: %.1f p/s, Exit: %.1f p/s)\n', ...
       %         feed_rate_factor, reactor_data.feed_rate, reactor_data.exit_rate);
    end
end

% Start feed callback
function startFeedCallback(~, ~)
    global continuous_feeding allUIComponents;
    
    continuous_feeding = true;
    
    % Update button text
    if isfield(allUIComponents, 'startFeedButton') && isvalid(allUIComponents.startFeedButton)
        if continuous_feeding
            allUIComponents.startFeedButton.Text = 'STOP FEED';
            allUIComponents.startFeedButton.ButtonPushedFcn = @stopFeedCallback;
        end
    end
    
    % Update feed toggle
    if isfield(allUIComponents, 'feedToggle') && isvalid(allUIComponents.feedToggle)
        allUIComponents.feedToggle.Value = true;
    end
    
    fprintf('Continuous feed system STARTED\n');
end

% Stop feed callback
function stopFeedCallback(~, ~)
    global continuous_feeding allUIComponents;
    
    continuous_feeding = false;
    
    % Update button text
    if isfield(allUIComponents, 'startFeedButton') && isvalid(allUIComponents.startFeedButton)
        if ~continuous_feeding
            allUIComponents.startFeedButton.Text = 'START FEED';
            allUIComponents.startFeedButton.ButtonPushedFcn = @startFeedCallback;
        end
    end
    
    % Update feed toggle
    if isfield(allUIComponents, 'feedToggle') && isvalid(allUIComponents.feedToggle)
        allUIComponents.feedToggle.Value = false;
    end
    
    fprintf('Continuous feed system STOPPED\n');
end

% Feed system toggle callback
function feedSystemToggleCallback(src, ~)
    global pipe_graphics;
    
    show_pipes = src.Value;
    
    if ~isempty(pipe_graphics)
        pipe_fields = fieldnames(pipe_graphics);
        for i = 1:length(pipe_fields)
            if isfield(pipe_graphics, pipe_fields{i}) && isvalid(pipe_graphics.(pipe_fields{i}))
                if show_pipes
                    pipe_graphics.(pipe_fields{i}).Visible = 'on';
                else
                    pipe_graphics.(pipe_fields{i}).Visible = 'off';
                end
            end
        end
    end
    
    if show_pipes
        fprintf('Feed system pipes shown\n');
    else
        fprintf('Feed system pipes hidden\n');
    end
end

% Pause button callback
function pauseButtonCallback(~, ~)
    global isPaused statusLamp;
    
    isPaused = ~isPaused;
    
    if isPaused
        fprintf('Continuous feed simulation PAUSED\n');
        if ~isempty(statusLamp) && isvalid(statusLamp)
            statusLamp.Color = [1, 1, 0.2];  % Yellow
        end
    else
        fprintf('Continuous feed simulation RESUMED\n');
        if ~isempty(statusLamp) && isvalid(statusLamp)
            statusLamp.Color = [0.2, 1, 0.2];  % Green
        end
    end
end

% Stop button callback
function stopButtonCallback(~, ~)
    stopSimulation();
end

% Reset button callback
function resetButtonCallback(~, ~)
    global reactor_data current_time continuous_feeding;
    
    stopSimulation();
    
    current_time = 0;
    
    if ~isempty(reactor_data)
        reactor_data.time = 0;
        reactor_data.frame_count = 0;
        reactor_data.conversion = 0;
        reactor_data.impeller_angles = zeros(1, reactor_data.impeller_count);
        
        % Reset feed system counters
        reactor_data.total_fed = 0;
        reactor_data.total_exited = 0;
        reactor_data.particles_in_system = 0;
        reactor_data.current_holdup = 0;
        
        % Clear species data
        reactor_data.species.ethyl_acetate = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
        reactor_data.species.naoh = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
        reactor_data.species.products = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
    end
    
    % Reset feed system
    continuous_feeding = true;
    
    fprintf('Continuous feed simulation RESET\n');
end

% Emergency stop callback
function emergencyStopCallback(~, ~)
    global powerSwitch statusLamp continuous_feeding;
    
    stopSimulation();
    
    % Stop feed system
    continuous_feeding = false;
    
    % Turn off main power
    if ~isempty(powerSwitch) && isvalid(powerSwitch)
        powerSwitch.Value = 'Off';
    end
    
    disableAllControls();
    
    if ~isempty(statusLamp) && isvalid(statusLamp)
        statusLamp.Color = [1, 0.2, 0.2];  % Red
    end
    
    fprintf('*** EMERGENCY STOP ACTIVATED - FEED SYSTEM STOPPED ***\n');
end

% Clear all particles callback
function clearAllParticlesCallback(~, ~)
    global reactor_data;
    
    if isempty(reactor_data)
        fprintf('Reactor not initialized\n');
        return;
    end
    
    % Clear all species particles including age arrays
    reactor_data.species.ethyl_acetate = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', [], 'age', []);
    reactor_data.species.naoh = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', [], 'age', []);
    reactor_data.species.products = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', [], 'age', []);
    
    % Reset conversion and counters
    reactor_data.conversion = 0;
    reactor_data.particles_in_system = 0;
    reactor_data.current_holdup = 0;
    
    % Update visualizations
    if isfield(reactor_data, 'species_plots')
        if isfield(reactor_data.species_plots, 'ethyl_acetate') && isvalid(reactor_data.species_plots.ethyl_acetate)
            set(reactor_data.species_plots.ethyl_acetate, 'XData', [], 'YData', [], 'ZData', []);
        end
        if isfield(reactor_data.species_plots, 'naoh') && isvalid(reactor_data.species_plots.naoh)
            set(reactor_data.species_plots.naoh, 'XData', [], 'YData', [], 'ZData', []);
        end
        if isfield(reactor_data.species_plots, 'products') && isvalid(reactor_data.species_plots.products)
            set(reactor_data.species_plots.products, 'XData', [], 'YData', [], 'ZData', []);
        end
    end
    
    fprintf('All particles cleared from reactor tank\n');
    fprintf('Feed system counters reset\n');
    fprintf('Conversion reset to 0.0%%\n');
end

% Flow toggle callback
function flowToggleCallback(src, ~)
    global show_flow_lines reactor_data;
    
    if isa(src, 'matlab.ui.control.CheckBox')
        show_flow_lines = src.Value;
    else
        show_flow_lines = ~show_flow_lines;
    end
    
    if ~isempty(reactor_data) && isfield(reactor_data, 'flow_lines')
        for i = 1:length(reactor_data.flow_lines)
            if isvalid(reactor_data.flow_lines(i))
                if show_flow_lines
                    reactor_data.flow_lines(i).Visible = 'on';
                else
                    reactor_data.flow_lines(i).Visible = 'off';
                end
            end
        end
    end
    
    if show_flow_lines
        fprintf('Flow lines enabled\n');
    else
        fprintf('Flow lines disabled\n');
    end
end

% Particle toggle callback
function particleToggleCallback(src, ~)
    global show_particles reactor_data;
    
    if isa(src, 'matlab.ui.control.CheckBox')
        show_particles = src.Value;
    else
        show_particles = ~show_particles;
    end
    
    if ~isempty(reactor_data) && isfield(reactor_data, 'species_plots')
        species_fields = fieldnames(reactor_data.species_plots);
        for i = 1:length(species_fields)
            if isfield(reactor_data.species_plots, species_fields{i}) && ...
               isvalid(reactor_data.species_plots.(species_fields{i}))
                if show_particles
                    reactor_data.species_plots.(species_fields{i}).Visible = 'on';
                else
                    reactor_data.species_plots.(species_fields{i}).Visible = 'off';
                end
            end
        end
    end
    
    if show_particles
        fprintf('Particle species enabled\n');
    else
        fprintf('Particle species disabled\n');
    end
end

% Density slider callback
function densitySliderCallback(src, ~)
    global particle_density densityLabel reactor_data;
    particle_density = src.Value;
    densityLabel.Text = sprintf('%.2fx', particle_density);
    
    % Update particle count
    if ~isempty(reactor_data)
        base_particle_count = 1000;
        reactor_data.particle_count = round(base_particle_count * particle_density);
        fprintf('Particle density set to %.2fx (Count: %d)\n', particle_density, reactor_data.particle_count);
    end
end

% Viscosity slider callback
function viscositySliderCallback(src, ~)
    global fluid_viscosity_factor viscosityLabel reactor_data;
    fluid_viscosity_factor = src.Value;
    viscosityLabel.Text = sprintf('%.2fx', fluid_viscosity_factor);
    
    % Update viscosity
    if ~isempty(reactor_data)
        reactor_data.viscosity = 1.2 * fluid_viscosity_factor;
        fprintf('Fluid viscosity set to %.2fx\n', fluid_viscosity_factor);
    end
end

% Stop simulation helper
function stopSimulation()
    global simTimer isRunning isStopped statusLamp;
    
    isStopped = true;
    isRunning = false;
    
    if ~isempty(simTimer) && isvalid(simTimer)
        stop(simTimer);
        delete(simTimer);
        simTimer = [];
    end
    
    if ~isempty(statusLamp) && isvalid(statusLamp)
        statusLamp.Color = [1, 0.2, 0.2];  % Red
    end
    
    fprintf('Continuous feed simulation STOPPED\n');
end

% Enable all controls
function enableAllControls()
    global knobControl startButton pauseButton stopButton resetButton;
    global flowToggle particleToggle densitySlider viscositySlider darkModeToggle;
    global clearParticlesButton feedToggle feedRateSlider allUIComponents;
    
    controls = {knobControl, startButton, pauseButton, stopButton, resetButton, ...
               flowToggle, particleToggle, densitySlider, viscositySlider, darkModeToggle, ...
               clearParticlesButton, feedToggle, feedRateSlider};
    
    for i = 1:length(controls)
        if ~isempty(controls{i}) && isvalid(controls{i})
            try
                controls{i}.Enable = 'on';
            catch
                % Ignore if control doesn't support Enable property
            end
        end
    end
    
    % Enable feed system buttons
    if isfield(allUIComponents, 'startFeedButton') && isvalid(allUIComponents.startFeedButton)
        allUIComponents.startFeedButton.Enable = 'on';
    end
    
    if isfield(allUIComponents, 'feedSystemToggle') && isvalid(allUIComponents.feedSystemToggle)
        allUIComponents.feedSystemToggle.Enable = 'on';
    end
     refreshAllButtonStates();
    
end

% Disable all controls
function disableAllControls()
    global knobControl startButton pauseButton stopButton resetButton;
    global flowToggle particleToggle densitySlider viscositySlider darkModeToggle;
    global clearParticlesButton feedToggle feedRateSlider allUIComponents;
    
    controls = {knobControl, pauseButton, stopButton, resetButton, ...
               flowToggle, particleToggle, densitySlider, viscositySlider, darkModeToggle, ...
               clearParticlesButton, feedToggle, feedRateSlider};
    
    for i = 1:length(controls)
        if ~isempty(controls{i}) && isvalid(controls{i})
            try
                controls{i}.Enable = 'off';
            catch
                % Ignore if control doesn't support Enable property
            end
        end
    end
    
    % Disable feed system buttons
    if isfield(allUIComponents, 'startFeedButton') && isvalid(allUIComponents.startFeedButton)
        allUIComponents.startFeedButton.Enable = 'off';
    end
    
    if isfield(allUIComponents, 'feedSystemToggle') && isvalid(allUIComponents.feedSystemToggle)
        allUIComponents.feedSystemToggle.Enable = 'off';
    end
    
    % Keep start button enabled
    if ~isempty(startButton) && isvalid(startButton)
        startButton.Enable = 'on';
    end
    
     % ADD THIS LINE HERE:
    refreshAllButtonStates();  % ADDED: Fix button states after disabling
    
end

% Figure close callback
function figureCloseCallback(src, ~)
    global simTimer;
    
    % Stop simulation timer
    if ~isempty(simTimer) && isvalid(simTimer)
        stop(simTimer);
        delete(simTimer);
    end
    
    % Close figure
    delete(src);
    
    fprintf('Enhanced CSTR with continuous feed system closed\n');
end

% Add this debug function to track product creation
function debugProductCreation()
    global reactor_data;
    
    if isempty(reactor_data) || isempty(reactor_data.species)
        fprintf('DEBUG: No reactor data\n');
        return;
    end
    
    ea_count = length(reactor_data.species.ethyl_acetate.x);
    naoh_count = length(reactor_data.species.naoh.x);
    product_count = length(reactor_data.species.products.x);
    
    fprintf('DEBUG: EtOAc=%d, NaOH=%d, Products=%d, RPM=%.0f, Conversion=%.1f%%\n', ...
            ea_count, naoh_count, product_count, reactor_data.impeller_speed, reactor_data.conversion);
    
    % Check if products plot exists and is valid
    if isfield(reactor_data, 'species_plots') && isfield(reactor_data.species_plots, 'products')
        if isvalid(reactor_data.species_plots.products)
            vis_status = reactor_data.species_plots.products.Visible;
            color_data = reactor_data.species_plots.products.CData;
            fprintf('DEBUG: Products plot valid, Visible=%s, Color=[%.1f,%.1f,%.1f]\n', ...
                    vis_status, color_data(1), color_data(2), color_data(3));
        else
            fprintf('DEBUG: Products plot invalid\n');
        end
    else
        fprintf('DEBUG: Products plot not found\n');
    end
end
% Force update all label colors - call this after dark mode toggle
function forceUpdateAllLabelColors()
    global dark_mode allUIComponents;
    
    if ~isfield(allUIComponents, 'mainFig') || ~isvalid(allUIComponents.mainFig)
        return;
    end
    
    % Get all labels and panel titles in the figure
    all_labels = findall(allUIComponents.mainFig, 'Type', 'uilabel');
    all_panels = findall(allUIComponents.mainFig, 'Type', 'uipanel');
    
    % Define colors
    if dark_mode
        regular_text_color = [0.9, 0.9, 0.95];  % Light for dark mode
        accent_color = [0.4, 0.7, 1.0];         % Bright blue for values
        panel_title_color = [0.9, 0.9, 0.95];    % White for panel titles
    else
        regular_text_color = [0.1, 0.2, 0.5];   % Dark for bright mode
        accent_color = [0.2, 0.4, 0.8];         % Blue for values
        panel_title_color = [0.1, 0.2, 0.5];     % Dark for panel titles
    end
    
    % Update each label
    for i = 1:length(all_labels)
        try
            label_text = all_labels(i).Text;
            
            % Check if this is a value display or regular label
            if contains(label_text, 'RPM') && ~isempty(regexp(label_text, '\d', 'once'))
                % This is an RPM value (like "300 RPM")
                all_labels(i).FontColor = accent_color;
            elseif contains(label_text, 'x') && ~isempty(regexp(label_text, '\d', 'once'))
                % This is a multiplier value (like "1.00x")
                all_labels(i).FontColor = accent_color;
            else
                % This is a regular label
                all_labels(i).FontColor = regular_text_color;
            end
            
        catch
            % If there's an error, just use regular text color
            try
                all_labels(i).FontColor = regular_text_color;
            catch
                % Skip this label if it can't be updated
            end
        end
    end
    
    % Update panel titles
    for i = 1:length(all_panels)
        try
            if ~isempty(all_panels(i).Title)
                all_panels(i).ForegroundColor = panel_title_color;
            end
        catch
            % Skip if there's an error
        end
    end
    
    drawnow;
end


% DIAGNOSTIC FUNCTION: Check array dimensions
function checkArrayDimensions()
    global reactor_data;
    
    if isempty(reactor_data) || isempty(reactor_data.species)
        fprintf('No reactor data to check\n');
        return;
    end
    
    species_names = {'ethyl_acetate', 'naoh', 'products'};
    field_names = {'x', 'y', 'z', 'vx', 'vy', 'vz', 'age'};
    
    fprintf('=== ARRAY DIMENSION CHECK ===\n');
    for s = 1:length(species_names)
        species_name = species_names{s};
        species = reactor_data.species.(species_name);
        
        fprintf('%s:\n', upper(species_name));
        for f = 1:length(field_names)
            field_name = field_names{f};
            if isfield(species, field_name) && ~isempty(species.(field_name))
                dims = size(species.(field_name));
                fprintf('  %s: [%dx%d] %s\n', field_name, dims(1), dims(2), ...
                        iif(dims(1) > dims(2), 'COLUMN', 'ROW'));
            else
                fprintf('  %s: EMPTY\n', field_name);
            end
        end
        fprintf('\n');
    end
end


% HELPER FUNCTION: Ensure array is a column vector
function vec = ensureColumnVector(vec)
    if isempty(vec)
        vec = [];
        return;
    end
    
    % Convert to column vector if it's a row vector
    if size(vec, 1) == 1 && size(vec, 2) > 1
        vec = vec';
    end
    
    % If it's already a column vector, leave it as is
    % If it's a scalar, convert to 1x1
    if isscalar(vec)
        vec = vec;  % Keep as scalar
    end
end


% Helper function to calculate advanced mixing efficiency
function mixing_eff = calculateMixingEfficiency(rpm)
    global reactor_data;
    
    if isempty(reactor_data)
        mixing_eff = 1.0;
        return;
    end
    
    % Advanced mixing model based on chemical engineering principles
    if rpm < 5
        mixing_eff = 0.1;  % Almost no mixing
    elseif rpm < 50
        mixing_eff = 0.1 + 0.4 * (rpm - 5) / 45;  % Linear increase to 0.5
    elseif rpm < 200
        mixing_eff = 0.5 + 1.0 * (rpm - 50) / 150;  % Linear to 1.5
    elseif rpm < 400
        mixing_eff = 1.5 + 1.0 * (rpm - 200) / 200;  % Linear to 2.5
    else
        mixing_eff = 2.5 + 0.5 * (rpm - 400) / 200;  % Diminishing returns to 3.0
        mixing_eff = min(mixing_eff, 3.0);  % Cap at 3.0x
    end
    
    % Store Reynolds number for display (optional)
    if isfield(reactor_data, 'impeller_diameter') && isfield(reactor_data, 'viscosity')
        Re = (rpm/60) * reactor_data.impeller_diameter^2 * reactor_data.density / reactor_data.viscosity;
        reactor_data.reynolds_number = Re;
    end
end
% === ADD THESE NEW FUNCTIONS AT THE END OF YOUR CODE ===

% HELPER FUNCTION: Reset button state to prevent sinking appearance
function resetButtonState(button)
    if ~isempty(button) && isvalid(button)
        try
            % Force button to normal (unpressed) state
            if isprop(button, 'Value')
                button.Value = false;  % For toggle buttons
            end
            
            % Ensure button is enabled and visible
            button.Enable = 'on';
            button.Visible = 'on';
            
            % Force refresh of button appearance
            drawnow;
            
        catch ME
            % Ignore errors - some properties might not exist
            fprintf('Warning: Could not reset button state: %s\n', ME.message);
        end
    end
end

% FUNCTION: Refresh all button states
function refreshAllButtonStates()
    global startButton pauseButton stopButton resetButton clearParticlesButton;
    global allUIComponents;
    
    fprintf('Refreshing all button states...\n');
    
    % Reset main operation buttons
    buttons_to_reset = {startButton, pauseButton, stopButton, resetButton, clearParticlesButton};
    
    for i = 1:length(buttons_to_reset)
        if ~isempty(buttons_to_reset{i}) && isvalid(buttons_to_reset{i})
            resetButtonState(buttons_to_reset{i});
        end
    end
    
    % Reset additional buttons stored in allUIComponents
    button_fields = {'startFeedButton', 'emergencyButton'};
    
    for i = 1:length(button_fields)
        field = button_fields{i};
        if isfield(allUIComponents, field) && isvalid(allUIComponents.(field))
            resetButtonState(allUIComponents.(field));
        end
    end
    
    % Force complete UI refresh
    drawnow;
    pause(0.05);
    drawnow;
    
    fprintf('All button states refreshed\n');
end

% DIAGNOSTIC FUNCTION: Check button states
function checkButtonStates()
    global startButton pauseButton stopButton resetButton clearParticlesButton;
    global allUIComponents;
    
    fprintf('=== BUTTON STATE DIAGNOSTIC ===\n');
    
    buttons = {
        {'startButton', startButton}, 
        {'pauseButton', pauseButton}, 
        {'stopButton', stopButton}, 
        {'resetButton', resetButton}, 
        {'clearParticlesButton', clearParticlesButton}
    };
    
    for i = 1:length(buttons)
        button_name = buttons{i}{1};
        button_obj = buttons{i}{2};
        
        if ~isempty(button_obj) && isvalid(button_obj)
            fprintf('%s:\n', button_name);
            fprintf('  Enable: %s\n', button_obj.Enable);
            fprintf('  Visible: %s\n', button_obj.Visible);
            
            if isprop(button_obj, 'Value')
                fprintf('  Value: %s\n', mat2str(button_obj.Value));
            end
            
            bg_color = button_obj.BackgroundColor;
            fprintf('  BackgroundColor: [%.2f, %.2f, %.2f]\n', bg_color(1), bg_color(2), bg_color(3));
        else
            fprintf('%s: INVALID OR EMPTY\n', button_name);
        end
        fprintf('\n');
    end
end
