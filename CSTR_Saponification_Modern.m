classdef CSTR_Saponification_Modern < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        MainGridLayout              matlab.ui.container.GridLayout
        ConversionInfo               matlab.ui.control.HTML
        
        % Left Panel Components
        LeftPanel                   matlab.ui.container.Panel
        LeftGridLayout              matlab.ui.container.GridLayout
        
        % Logo
        LogoHTML                    matlab.ui.control.HTML
        
        % Data Section
        DataPanel                   matlab.ui.container.Panel
        DataGrid                    matlab.ui.container.GridLayout
        UploadButton                matlab.ui.control.Button
        FileLabel                   matlab.ui.control.Label
        
        % Parameters Section
        ParamsPanel                 matlab.ui.container.Panel
        ParamsGrid                  matlab.ui.container.GridLayout
        TempDropDown                matlab.ui.control.DropDown
        TempLabel                   matlab.ui.control.Label
        VolumeField                 matlab.ui.control.NumericEditField
        VolumeLabel                 matlab.ui.control.Label
        FlowrateField               matlab.ui.control.NumericEditField
        FlowrateLabel               matlab.ui.control.Label
        CinField                    matlab.ui.control.NumericEditField
        CinLabel                    matlab.ui.control.Label
        SimTimeField                matlab.ui.control.NumericEditField
        SimTimeLabel                matlab.ui.control.Label
        
        % Actions Section
        ActionsPanel                matlab.ui.container.Panel
        ActionsGrid                 matlab.ui.container.GridLayout
        ManualKButton               matlab.ui.control.Button
        OptimizeButton              matlab.ui.control.Button
        OptimizeUltraButton         matlab.ui.control.Button
        ClearAllButton              matlab.ui.control.Button
        
        % Display Options
        DisplayPanel                matlab.ui.container.Panel
        DisplayGrid                 matlab.ui.container.GridLayout
        JacketSwitch                matlab.ui.control.Switch
        JacketLabel                 matlab.ui.control.Label
        ProductsSwitch              matlab.ui.control.Switch
        ProductsLabel               matlab.ui.control.Label
        
        % Export Section
        ExportPanel                 matlab.ui.container.Panel
        ExportGrid                  matlab.ui.container.GridLayout
        SaveButton                  matlab.ui.control.Button
        ExportDataButton            matlab.ui.control.Button
        ExportPlotsButton           matlab.ui.control.Button
        
        % Status Section
        StatusPanel                 matlab.ui.container.Panel
        StatusGrid                  matlab.ui.container.GridLayout
        StatusLabel                 matlab.ui.control.Label
        ProgressLabel               matlab.ui.control.Label
        
        % Right Panel - Visualization
        RightPanel                  matlab.ui.container.Panel
        RightGrid                   matlab.ui.container.GridLayout
        TabGroup                    matlab.ui.container.TabGroup
        
        % Tab 1: Concentrations
        ConcentrationTab            matlab.ui.container.Tab
        ConcentrationGrid           matlab.ui.container.GridLayout
        UIAxes1                     matlab.ui.control.UIAxes
        ConcentrationInfo           matlab.ui.control.HTML
        
        % Tab 2: Temperature
        TemperatureTab              matlab.ui.container.Tab
        TemperatureGrid             matlab.ui.container.GridLayout
        UIAxes2                     matlab.ui.control.UIAxes
        TemperatureInfo             matlab.ui.control.HTML
        
        % Tab 3: Products
        ProductsTab                 matlab.ui.container.Tab
        ProductsGrid                matlab.ui.container.GridLayout
        UIAxes3                     matlab.ui.control.UIAxes
        ProductsInfo                matlab.ui.control.HTML
        
        % Tab 4: Analysis
        AnalysisTab                 matlab.ui.container.Tab
        AnalysisGrid                matlab.ui.container.GridLayout
        ConvergencePanel            matlab.ui.container.Panel
        ConvergenceText             matlab.ui.control.TextArea
        ComparisonTable             matlab.ui.control.Table
        ConcentrationTablePanel     matlab.ui.container.Panel
        ConcentrationTable          matlab.ui.control.Table
        TemperatureTablePanel       matlab.ui.container.Panel
        ComparisonTablePanel        matlab.ui.container.Panel 
        TemperatureTable            matlab.ui.control.Table
        
        % NEW TAB 5: Quick Sim
        QuickSimTab                 matlab.ui.container.Tab
        QuickSimGrid                matlab.ui.container.GridLayout
        CinQuickField               matlab.ui.control.NumericEditField
        FlowrateQuickField          matlab.ui.control.NumericEditField
        SimTimeQuickField           matlab.ui.control.NumericEditField
        RateConstQuickField         matlab.ui.control.NumericEditField
        TempQuickDropDown           matlab.ui.control.DropDown
        RunQuickSimButton           matlab.ui.control.Button
        UIAxesQuickConc             matlab.ui.control.UIAxes
        UIAxesQuickTemp             matlab.ui.control.UIAxes
        UIAxesQuickProducts         matlab.ui.control.UIAxes
        MaxErrorLabel               matlab.ui.control.Label
        MinErrorLabel               matlab.ui.control.Label
        MeanErrorLabel              matlab.ui.control.Label
        QuickSimInfo                matlab.ui.control.HTML
        % Tab 6: Conversion Analysis
ConversionTab               matlab.ui.container.Tab
ConversionGrid              matlab.ui.container.GridLayout
UIAxesConversionProduct     matlab.ui.control.UIAxes
UIAxesConversionReactant    matlab.ui.control.UIAxes

% Context Menus
PlotContextMenu             matlab.ui.container.ContextMenu
        
        CopyDataMenuItem            matlab.ui.container.Menu
        SaveImageMenuItem           matlab.ui.container.Menu
        
        % Loading UI
        LoadingDialog               matlab.ui.dialog.ProgressDialog
        
         % NEW: Enhanced buttons (
    ValidateDataButton          matlab.ui.control.Button
    ExportEnhancedButton        matlab.ui.control.Button

    end
    
    % Properties for app functionality
    properties (Access = private)
        CoolingFlowrateHistory    double    % [L/min] adaptive coolant flow vs time
CoolingTimeHistory        double    % [min] time points for cooling data
CoolingCapacityHistory    double    % [J/min] required cooling capacity vs time
CoolingAdequacyHistory    double    % [-] adequacy factor vs time
TemperatureDeviationHistory double  % [°C] temperature deviation from setpoint vs time
TempDifferenceHistory       double    % [°C] actual temp difference vs time
    TotalCoolingRequiredHistory double    % [J/min] total cooling required vs time
          % ADD THESE NEW PROPERTIES:
    % NEW: Adaptive Cooling System Properties
    Ea = 30000;                 % Activation energy [J/mol]
    A_preexp = 3.187e5;         % Pre-exponential factor [L/mol·min]
    R_gas = 8.314;             % Gas constant [J/(mol·K)]
    
    Fj_min = 0.1;                   % [L/min] minimum coolant flow
    Fj_max = 50;                    % [L/min] maximum coolant flow  
    Fj_adaptive                     % [L/min] current adaptive coolant flow
    isothermal_tolerance = 0.1;     % [°C] temperature control tolerance
    cooling_control_gain = 2.0;     % PI controller gain for cooling
    integral_error = 0;             % Accumulated temperature error
    IsothermalMode = true;          % Enable/disable isothermal control
    
    % Cooling system capacity tracking
    max_cooling_capacity            % [J/min] maximum theoretical cooling
    required_cooling_capacity       % [J/min] currently required cooling
    cooling_adequacy_factor         % Ratio of available/required cooling
    
    % UI Components for isothermal control
    IsothermalLabel                 matlab.ui.control.Label
    IsothermalSwitch                matlab.ui.control.Switch
        
        Data                        table
        SelectedTemp                double = 30
        V = 2.5;                    % Reactor volume [L]
        V0 = 0.18;                  % Total volumetric flow rate [L/min] (180 ml/min = 0.18 L/min)
        V0_EtOAc = 0.09;            % EtOAc flowrate [L/min]
        V1_NaOH = 0.09;             % NaOH flowrate [L/min]
        Cin = 0.05;                 % Inlet concentration [mol/L]
        SimTime = 30;               % Simulation time [min]
        k_fit                       double
        k_manual                    double
        SimResults                  struct
        ManualSimResults            struct
        QuickSimResults             struct  % NEW: Quick simulation results
        ShowJacket                  logical = false
        ShowProducts                logical = true
        RSquared                    double = NaN
        RSquaredManual              double = NaN
        OptimizationResults         struct
        ConvergenceData             struct
        ManualConvergenceData       struct
        QuickSimConvergenceData     struct  % NEW: Quick sim convergence data
        ComparisonTableData         table
        ComparisonTableDataManual   table
        
        
        % ENHANCED CONVERGENCE TRACKING with different criteria
        IsConverged                 logical = false
        IsManualConverged           logical = false
        IsQuickSimConverged         logical = false  % NEW: Quick sim convergence status
        

        % Manual simulation convergence thresholds (more relaxed)
        ManualReactantThreshold     double  = 0.5  % Higher threshold for manual reactants
        ManualProductThreshold      double  = 0.5 % Higher threshold for manual temperature
        ManualTemperatureThreshold  double  = 0.5  % Higher threshold for manual temperature   % <-- ADD THIS LINE

        % Optimized simulation convergence thresholds (stricter)
        OptimizedReactantThreshold  double  = 0.5  % Stricter threshold for optimized reactants
        OptimizedProductThreshold   double  = 0.5  % Stricter threshold for optimized products
        OptimizedTemperatureThreshold double = 0.5  % Stricter threshold for optimized temperature
        
        % Quick simulation convergence thresholds (medium)
        QuickSimReactantThreshold   double = 0.5   % Medium threshold for quick sim reactants
        QuickSimProductThreshold    double  = 0.5  % Medium threshold for quick sim products
        QuickSimTemperatureThreshold double  = 0.5 % Medium threshold for quick sim temperature
        
        % SEPARATE CONVERGENCE STATUS TRACKING
        % Manual convergence status
        ManualReactantsConverged    logical = false
        ManualProductsConverged     logical = false
        ManualTemperatureConverged  logical = false
        
        % Optimized convergence status
        OptimizedReactantsConverged logical = false
        OptimizedProductsConverged  logical = false
        OptimizedTemperatureConverged logical = false
        
        % Quick sim convergence status
        QuickSimReactantsConverged  logical = false
        QuickSimProductsConverged   logical = false
        QuickSimTemperatureConverged logical = false
        
        OriginalSimTime             double = NaN
OriginalEndConcentration    double = NaN  
ExtensionMode               logical = false
        
        % CONSISTENT REACTOR PARAMETERS (all in consistent units)
        Fj = 1;                     % [L/min] coolant flow - converted from m³/h to L/min
        Vj = 0.086;                  % [L] jacket volume - converted from m³ to L
        pj = 999.7;                 % [kg/m³] coolant density
        Cpj = 4030;                 % [J/(kg·K)] coolant heat capacity
        Tj0 = 25;                   % [°C] coolant inlet temperature
        
        % Heat transfer coefficients
        h_reactor = 100;            % [W/(m²·K)] 
        h_coolant = 1200;           % [W/(m²·K)]
        delta_wall = 0.00001;         % [m] wall thickness
        k_wall = 16.3;                % [W/(m·K)] thermal conductivity
        
        % Surface area (estimated based on volume)
        a = 0.102;                   % [m²] reactor surface area
        UA                          % Overall heat transfer coefficient × area [J/(min·K)] - CONSISTENT UNITS
        
        % Design System
        AccentColor = [0 122 255]/255;      % iOS Blue
        SuccessColor = [52 199 89]/255;     % iOS Green
        WarningColor = [255 149 0]/255;     % iOS Orange
        ErrorColor = [255 59 48]/255;       % iOS Red
        BackgroundColor = [242 242 247]/255; % iOS System Gray 6
        SurfaceColor = [255 255 255]/255;   % White
        TextColor = [0 0 0]/255;            % Black
        SecondaryTextColor = [142 142 147]/255; % iOS Secondary Label
        
        % Animation timers
        LoadingTimer                timer
        IsAnimating                 logical = false
        
        UltraPrecisionResults       struct
ComparisonTableDataUltra    table
RSquaredUltra              double = NaN
k_ultra                    double
UltraConvergenceData       struct
IsUltraConverged           logical = false
UltraReactantsConverged    logical = false
UltraProductsConverged     logical = false
UltraTemperatureConverged  logical = false

    end
    
    methods (Access = private)
        % ADD THIS NEW METHOD:
    function k_arrhenius = calculateArrheniusRateConstant(app, T_celsius)
        % Calculate rate constant using Arrhenius equation
        % T_celsius: Temperature in Celsius
        % Returns: k in L/mol·min
        
        T_kelvin = T_celsius + 273.15;  % Convert to Kelvin
        k_arrhenius = app.A_preexp * exp(-app.Ea / (app.R_gas * T_kelvin));
    end
    
    % ADD THIS NEW METHOD:
    function dydt = odefun_arrhenius_temperature(app, t, y, V0_EtOAc, V1_NaOH, V, Cin, To)
    try
        % Extract variables
        CA = y(1);      % NaOH concentration [mol/L]
        CB = y(2);      % EtOAc concentration [mol/L]
        CC = y(3);      % NaOAc concentration [mol/L]
        CD = y(4);      % EtOH concentration [mol/L]
        T_reactor = y(5);  % Reactor temperature [°C]
        Tj = y(6);      % Jacket temperature [°C]
        
        % Calculate rate constant using Arrhenius equation
        k = app.calculateArrheniusRateConstant(T_reactor);
        
        % Total flowrate [L/min]
        V0_total = V0_EtOAc + V1_NaOH;
        
        % Calculate mixture properties
        [p_avg, Cp_avg] = app.calculate_mixture_properties(CA, CB, CC, CD);
        
        % Reaction rate (second-order) [mol/(L·min)]
        r = k * CA * CB;
        
        % Temperature-dependent heat of reaction [J/mol]
        lambda_T = app.heat_of_reaction(T_reactor + 273.15);
        
        % Heat generation from reaction [J/min]
        Q_rxn = -lambda_T * r * V;
        
        % MODIFIED: Different cooling behavior based on isothermal mode
        if app.IsothermalMode
            % ISOTHERMAL MODE: Active temperature control with Arrhenius kinetics
            % Calculate temperature error
            temp_error = T_reactor - To;
            
            % PI controller for cooling flowrate adjustment
            app.integral_error = app.integral_error + temp_error;
            
            % Calculate required cooling to maintain isothermal conditions
            Q_base_required = abs(Q_rxn);  % Base cooling for reaction heat
            
            % Additional cooling for temperature correction (PI control)
            Q_correction = app.cooling_control_gain * temp_error + ...
                          0.1 * app.cooling_control_gain * app.integral_error;
            
            % Total required cooling
            Q_total_required = Q_base_required + Q_correction;
            
            % Calculate required coolant flowrate using ACTUAL reactor temperature
            temp_difference = max(T_reactor - app.Tj0, 0.1);
            Fj_required = abs(Q_total_required) / ...
                         (app.pj/1000 * app.Cpj * temp_difference);
            
            % Apply limits and update adaptive flowrate
            app.Fj_adaptive = max(app.Fj_min, min(app.Fj_max, Fj_required));
            
            % Store cooling system status
            app.required_cooling_capacity = abs(Q_total_required);
            app.max_cooling_capacity = app.Fj_max * app.pj/1000 * app.Cpj * temp_difference;
            app.cooling_adequacy_factor = app.max_cooling_capacity / max(app.required_cooling_capacity, 1);
            
            % Effective heat transfer with adaptive cooling
            heat_transfer_efficiency = 1.0;  % Full efficiency in isothermal mode
            
        else
            % NON-ISOTHERMAL MODE WITH ARRHENIUS: NO COOLING JACKET EXISTS
            % Zero coolant flowrate - jacket is completely disconnected
            app.Fj_adaptive = 0;  % NO cooling - jacket doesn't exist
            
            % No heat transfer efficiency - no jacket system
            heat_transfer_efficiency = 0;  % 0% efficient - no cooling jacket
            
            % Store cooling system status for non-isothermal mode (no jacket)
            app.required_cooling_capacity = abs(Q_rxn);  % Heat that needs removal but can't be
            app.max_cooling_capacity = 0;  % NO cooling capacity - no jacket
            app.cooling_adequacy_factor = 0;  % 0% adequacy - no cooling available
        end
        
        % Mass balance equations [mol/L·min]
        dCA_dt = (V1_NaOH/V) * (Cin - CA) - r;
        dCB_dt = (V0_EtOAc/V) * (Cin - CB) - r;
        dCC_dt = (V0_total/V) * (0 - CC) + r;
        dCD_dt = (V0_total/V) * (0 - CD) + r;
        
        % MODIFIED: Heat transfer calculation with efficiency factor
        if app.IsothermalMode
            Q_jacket_theoretical = app.UA * (T_reactor - Tj);
            Q_jacket = Q_jacket_theoretical * heat_transfer_efficiency;
        else
            % NON-ISOTHERMAL MODE: NO JACKET EXISTS
            Q_jacket_theoretical = 0;  % No jacket to transfer heat to
            Q_jacket = 0;  % No heat transfer to jacket
        end
        
        % Convert mixture density from [kg/m³] to [kg/L]
        p_avg_kg_L = p_avg / 1000;
        
        % MODIFIED: Reactor energy balance with enhanced temperature effects
        if app.IsothermalMode
            % For isothermal operation with Arrhenius kinetics, strong control
            dT_reactor_dt = (V0_total/V) * (To - T_reactor) + ...
                           (Q_rxn) / (p_avg_kg_L * Cp_avg * V) - ...
                           (Q_jacket) / (p_avg_kg_L * Cp_avg * V) - ...
                           5 * (T_reactor - To);  % Strong temperature control for Arrhenius
        else
            % Non-isothermal with Arrhenius: Allow significant temperature rise
            % This creates a positive feedback loop where higher T -> higher k -> more reaction -> more heat -> higher T
            dT_reactor_dt = (V0_total/V) * (To - T_reactor) + ...
                           (Q_rxn) / (p_avg_kg_L * Cp_avg * V);
            % NO cooling term - let Arrhenius feedback drive temperature up with NO heat removal!
        end
        
        % MODIFIED: Jacket energy balance
        if app.IsothermalMode
            % Isothermal mode: Effective cooling with adaptive flowrate
            pj_kg_L = app.pj / 1000;
            dTj_dt = (app.Fj_adaptive / app.Vj) * (app.Tj0 - Tj) + ...
                    (Q_jacket) / (pj_kg_L * app.Cpj * app.Vj);
        else
            % Non-isothermal mode: NO JACKET EXISTS
            dTj_dt = 0;  % No jacket temperature change - jacket doesn't exist
            % Keep Tj at inlet temperature since jacket is disconnected
        end
        
        % Store cooling history data for Arrhenius simulation
        temp_deviation = abs(T_reactor - To);
        app.storeCoolingHistory(t, app.Fj_adaptive, app.required_cooling_capacity, ...
                               app.cooling_adequacy_factor, temp_deviation);
        
        % Return derivatives
        dydt = [dCA_dt; dCB_dt; dCC_dt; dCD_dt; dT_reactor_dt; dTj_dt];
        
    catch
        dydt = zeros(6,1);
    end
end
        % NEW METHOD: Update Conversion Info Panel
        function exportCombinedConcentrationPlots(app, folder, timestamp)
    % Export combined concentration profile plots (reactants + products)
    % Creates two versions: 15 min and full simulation time
    
    try
        exportedCount = 0;
        
        % Check if we have simulation data to export
        hasOptimized = ~isempty(app.SimResults) && isfield(app.SimResults, 't');
        hasManual = ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't');
        hasUltra = ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't');
        hasQuickSim = ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't');
        
        if ~(hasOptimized || hasManual || hasUltra || hasQuickSim)
            return;
        end
        
        % 1. COMBINED CONCENTRATION PROFILES - 15 MINUTES
        app.updateLoadingDialog(0.15, 'Creating 15-minute combined concentration plot...');
        
        try
            fig1 = figure('Visible', 'off', 'Position', [100, 100, 1400, 900]);
            ax1 = axes(fig1, 'Position', [0.08, 0.12, 0.85, 0.78]);
            
            hold(ax1, 'on');
            
            % Set time limit to 15 minutes
            time_limit = 15;
            
            % Plot experimental data if available
            if ~isempty(app.Data)
                exp_mask = app.Data.Time <= time_limit;
                if sum(exp_mask) > 0
                    scatter(ax1, app.Data.Time(exp_mask), app.Data.Concentration(exp_mask), 120, ...
                           'MarkerEdgeColor', [0.8 0.2 0.2], 'MarkerFaceColor', [1 0.8 0.8], ...
                           'LineWidth', 3, 'DisplayName', 'Experimental NaOH', 'Marker', 'o');
                end
            end
            
            % Plot Manual Simulation (if available)
            if hasManual
                manual_mask = app.ManualSimResults.t <= time_limit;
                if sum(manual_mask) > 0
                    % Reactants (Manual)
                    plot(ax1, app.ManualSimResults.t(manual_mask), app.ManualSimResults.C_NaOH(manual_mask), ':', ...
                         'LineWidth', 4, 'Color', [0.9 0.3 0.3], 'DisplayName', 'Manual NaOH');
                    plot(ax1, app.ManualSimResults.t(manual_mask), app.ManualSimResults.C_EtOAc(manual_mask), ':', ...
                         'LineWidth', 4, 'Color', [0.9 0.6 0.2], 'DisplayName', 'Manual EtOAc');
                    
                    % Products (Manual)
                    plot(ax1, app.ManualSimResults.t(manual_mask), app.ManualSimResults.C_NaOAc(manual_mask), ':', ...
                         'LineWidth', 4, 'Color', [0.8 0.7 0.2], 'DisplayName', 'Manual NaOAc');
                    plot(ax1, app.ManualSimResults.t(manual_mask), app.ManualSimResults.C_EtOH(manual_mask), ':', ...
                         'LineWidth', 4, 'Color', [0.6 0.4 0.8], 'DisplayName', 'Manual EtOH');
                end
            end
            
            % Plot Optimized Simulation (if available)
            if hasOptimized
                opt_mask = app.SimResults.t <= time_limit;
                if sum(opt_mask) > 0
                    % Reactants (Optimized)
                    plot(ax1, app.SimResults.t(opt_mask), app.SimResults.C_NaOH(opt_mask), '-', ...
                         'LineWidth', 4, 'Color', [0.2 0.4 0.8], 'DisplayName', 'Optimized NaOH');
                    plot(ax1, app.SimResults.t(opt_mask), app.SimResults.C_EtOAc(opt_mask), '-', ...
                         'LineWidth', 4, 'Color', [0.3 0.7 0.3], 'DisplayName', 'Optimized EtOAc');
                    
                    % Products (Optimized)
                    plot(ax1, app.SimResults.t(opt_mask), app.SimResults.C_NaOAc(opt_mask), '--', ...
                         'LineWidth', 4, 'Color', [0.9 0.6 0.0], 'DisplayName', 'Optimized NaOAc');
                    plot(ax1, app.SimResults.t(opt_mask), app.SimResults.C_EtOH(opt_mask), '--', ...
                         'LineWidth', 4, 'Color', [0.5 0.3 0.7], 'DisplayName', 'Optimized EtOH');
                end
            end
            
            % Plot Ultra-Precision Simulation (if available)
            if hasUltra
                ultra_mask = app.UltraPrecisionResults.t <= time_limit;
                if sum(ultra_mask) > 0
                    % Reactants (Ultra-Precision)
                    plot(ax1, app.UltraPrecisionResults.t(ultra_mask), app.UltraPrecisionResults.C_NaOH(ultra_mask), '-', ...
                         'LineWidth', 5, 'Color', [1 0.08 0.58], 'DisplayName', 'Ultra-Precision NaOH');
                    plot(ax1, app.UltraPrecisionResults.t(ultra_mask), app.UltraPrecisionResults.C_EtOAc(ultra_mask), '-', ...
                         'LineWidth', 5, 'Color', [0.54 0 0.54], 'DisplayName', 'Ultra-Precision EtOAc');
                    
                    % Products (Ultra-Precision)
                    plot(ax1, app.UltraPrecisionResults.t(ultra_mask), app.UltraPrecisionResults.C_NaOAc(ultra_mask), '-', ...
                         'LineWidth', 5, 'Color', [1 0.27 0], 'DisplayName', 'Ultra-Precision NaOAc');
                    plot(ax1, app.UltraPrecisionResults.t(ultra_mask), app.UltraPrecisionResults.C_EtOH(ultra_mask), '-', ...
                         'LineWidth', 5, 'Color', [0.13 0.55 0.13], 'DisplayName', 'Ultra-Precision EtOH');
                end
            end
            
            % Plot Quick Simulation (if available)
            if hasQuickSim
                quick_mask = app.QuickSimResults.t <= time_limit;
                if sum(quick_mask) > 0
                    % Reactants (Quick Sim)
                    plot(ax1, app.QuickSimResults.t(quick_mask), app.QuickSimResults.C_NaOH(quick_mask), '-', ...
                         'LineWidth', 3, 'Color', [0.0 0.5 0.0], 'DisplayName', 'Quick Sim NaOH');
                    plot(ax1, app.QuickSimResults.t(quick_mask), app.QuickSimResults.C_EtOAc(quick_mask), '-', ...
                         'LineWidth', 3, 'Color', [1.0 0.65 0], 'DisplayName', 'Quick Sim EtOAc');
                    
                    % Products (Quick Sim)
                    plot(ax1, app.QuickSimResults.t(quick_mask), app.QuickSimResults.C_NaOAc(quick_mask), '--', ...
                         'LineWidth', 3, 'Color', [0.18 0.31 0.31], 'DisplayName', 'Quick Sim NaOAc');
                    plot(ax1, app.QuickSimResults.t(quick_mask), app.QuickSimResults.C_EtOH(quick_mask), '--', ...
                         'LineWidth', 3, 'Color', [0.5 0.0 0.5], 'DisplayName', 'Quick Sim EtOH');
                end
            end
            
            % Formatting for 15-minute plot
            xlim(ax1, [0 time_limit]);
            ylim(ax1, [0 max(app.Cin * 1.1, 0.055)]);
            
            % Create title with parameter information
            if ~isnan(app.k_fit)
                titleStr = sprintf('Complete Concentration Profiles (15 min) - k_{opt} = %.4f L/mol·min, V₀ = %.3f L/min', ...
                                  app.k_fit, app.FlowrateField.Value);
            elseif ~isnan(app.k_manual)
                titleStr = sprintf('Complete Concentration Profiles (15 min) - k_{manual} = %.4f L/mol·min', app.k_manual);
            else
                titleStr = 'Complete Concentration Profiles (15 min) - Reactants & Products';
            end
            
            title(ax1, titleStr, 'FontWeight', 'bold', 'FontSize', 16);
            xlabel(ax1, 'Time (min)', 'FontWeight', 'bold', 'FontSize', 14);
            ylabel(ax1, 'Concentration (mol/L)', 'FontWeight', 'bold', 'FontSize', 14);
            
            % Enhanced legend
            leg1 = legend(ax1, 'Location', 'eastoutside', 'FontSize', 11);
            leg1.Title.String = 'Species & Methods';
            leg1.Title.FontWeight = 'bold';
            
            grid(ax1, 'off');
            ax1.GridAlpha = 0.3;
            ax1.Box = 'on';
            ax1.FontSize = 12;
            
            % Add conversion information text box
            if hasOptimized
                final_conversion = (app.Cin - app.SimResults.C_NaOH(end)) / app.Cin * 100;
                text(ax1, time_limit * 0.7, app.Cin * 0.95, ...
                     sprintf('Final Conversion: %.1f%%', final_conversion), ...
                     'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white', ...
                     'EdgeColor', 'black', 'Margin', 5);
            end
            
            filename1 = sprintf('Combined_Concentration_15min_%s.png', timestamp);
            if exist('exportgraphics', 'file')
                exportgraphics(ax1, fullfile(folder, filename1), 'Resolution', 300);
            else
                print(fig1, fullfile(folder, filename1), '-dpng', '-r300');
            end
            close(fig1);
            exportedCount = exportedCount + 1;
            
        catch ME
            fprintf('15-minute combined plot export failed: %s\n', ME.message);
        end
        
        % 2. COMBINED CONCENTRATION PROFILES - FULL SIMULATION TIME
        app.updateLoadingDialog(0.85, 'Creating full-time combined concentration plot...');
        
        % Only create full-time plot if simulation time > 15 minutes
        full_sim_time = app.SimTimeField.Value;
        if full_sim_time > 15
            try
                fig2 = figure('Visible', 'off', 'Position', [100, 100, 1400, 900]);
                ax2 = axes(fig2, 'Position', [0.08, 0.12, 0.85, 0.78]);
                
                hold(ax2, 'on');
                
                % Plot experimental data if available
                if ~isempty(app.Data)
                    scatter(ax2, app.Data.Time, app.Data.Concentration, 120, ...
                           'MarkerEdgeColor', [0.8 0.2 0.2], 'MarkerFaceColor', [1 0.8 0.8], ...
                           'LineWidth', 3, 'DisplayName', 'Experimental NaOH', 'Marker', 'o');
                end
                
                % Plot Manual Simulation (if available)
                if hasManual
                    % Reactants (Manual)
                    plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_NaOH, ':', ...
                         'LineWidth', 4, 'Color', [0.9 0.3 0.3], 'DisplayName', 'Manual NaOH');
                    plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_EtOAc, ':', ...
                         'LineWidth', 4, 'Color', [0.9 0.6 0.2], 'DisplayName', 'Manual EtOAc');
                    
                    % Products (Manual)
                    plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_NaOAc, ':', ...
                         'LineWidth', 4, 'Color', [0.8 0.7 0.2], 'DisplayName', 'Manual NaOAc');
                    plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_EtOH, ':', ...
                         'LineWidth', 4, 'Color', [0.6 0.4 0.8], 'DisplayName', 'Manual EtOH');
                end
                
                % Plot Optimized Simulation (if available)
                if hasOptimized
                    % Reactants (Optimized)
                    plot(ax2, app.SimResults.t, app.SimResults.C_NaOH, '-', ...
                         'LineWidth', 4, 'Color', [0.2 0.4 0.8], 'DisplayName', 'Optimized NaOH');
                    plot(ax2, app.SimResults.t, app.SimResults.C_EtOAc, '-', ...
                         'LineWidth', 4, 'Color', [0.3 0.7 0.3], 'DisplayName', 'Optimized EtOAc');
                    
                    % Products (Optimized)
                    plot(ax2, app.SimResults.t, app.SimResults.C_NaOAc, '--', ...
                         'LineWidth', 4, 'Color', [0.9 0.6 0.0], 'DisplayName', 'Optimized NaOAc');
                    plot(ax2, app.SimResults.t, app.SimResults.C_EtOH, '--', ...
                         'LineWidth', 4, 'Color', [0.5 0.3 0.7], 'DisplayName', 'Optimized EtOH');
                end
                
                % Plot Ultra-Precision Simulation (if available)
                if hasUltra
                    % Reactants (Ultra-Precision)
                    plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, '-', ...
                         'LineWidth', 5, 'Color', [1 0.08 0.58], 'DisplayName', 'Ultra-Precision NaOH');
                    plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOAc, '-', ...
                         'LineWidth', 5, 'Color', [0.54 0 0.54], 'DisplayName', 'Ultra-Precision EtOAc');
                    
                    % Products (Ultra-Precision)
                    plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
                         'LineWidth', 5, 'Color', [1 0.27 0], 'DisplayName', 'Ultra-Precision NaOAc');
                    plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
                         'LineWidth', 5, 'Color', [0.13 0.55 0.13], 'DisplayName', 'Ultra-Precision EtOH');
                end
                
                % Plot Quick Simulation (if available)
                if hasQuickSim
                    % Reactants (Quick Sim)
                    plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_NaOH, '-', ...
                         'LineWidth', 3, 'Color', [0.0 0.5 0.0], 'DisplayName', 'Quick Sim NaOH');
                    plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_EtOAc, '-', ...
                         'LineWidth', 3, 'Color', [1.0 0.65 0], 'DisplayName', 'Quick Sim EtOAc');
                    
                    % Products (Quick Sim)
                    plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_NaOAc, '--', ...
                         'LineWidth', 3, 'Color', [0.18 0.31 0.31], 'DisplayName', 'Quick Sim NaOAc');
                    plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_EtOH, '--', ...
                         'LineWidth', 3, 'Color', [0.5 0.0 0.5], 'DisplayName', 'Quick Sim EtOH');
                end
                
                % Formatting for full-time plot
                xlim(ax2, [0 full_sim_time]);
                ylim(ax2, [0 max(app.Cin * 1.1, 0.055)]);
                
                % Create title with parameter information
                if ~isnan(app.k_fit)
                    titleStr = sprintf('Complete Concentration Profiles (%.0f min) - k_{opt} = %.4f L/mol·min, V₀ = %.3f L/min', ...
                                      full_sim_time, app.k_fit, app.FlowrateField.Value);
                elseif ~isnan(app.k_manual)
                    titleStr = sprintf('Complete Concentration Profiles (%.0f min) - k_{manual} = %.4f L/mol·min', ...
                                      full_sim_time, app.k_manual);
                else
                    titleStr = sprintf('Complete Concentration Profiles (%.0f min) - Reactants & Products', full_sim_time);
                end
                
                title(ax2, titleStr, 'FontWeight', 'bold', 'FontSize', 16);
                xlabel(ax2, 'Time (min)', 'FontWeight', 'bold', 'FontSize', 14);
                ylabel(ax2, 'Concentration (mol/L)', 'FontWeight', 'bold', 'FontSize', 14);
                
                % Enhanced legend
                leg2 = legend(ax2, 'Location', 'eastoutside', 'FontSize', 11);
                leg2.Title.String = 'Species & Methods';
                leg2.Title.FontWeight = 'bold';
                
                grid(ax2, 'off');
                ax2.GridAlpha = 0.3;
                ax2.Box = 'on';
                ax2.FontSize = 12;
                
                % Add conversion and flowrate information text box
                text_y_pos = app.Cin * 0.95;
                if hasOptimized
                    final_conversion = (app.Cin - app.SimResults.C_NaOH(end)) / app.Cin * 100;
                    total_flowrate = app.FlowrateField.Value;
                    each_flowrate = total_flowrate / 2;
                    
                    info_text = sprintf(['Final Conversion: %.1f%%\n' ...
                                        'Total Flowrate: %.3f L/min\n' ...
                                        'Each Reactant: %.3f L/min'], ...
                                        final_conversion, total_flowrate, each_flowrate);
                    
                    text(ax2, full_sim_time * 0.65, text_y_pos, info_text, ...
                         'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', 'white', ...
                         'EdgeColor', 'black', 'Margin', 5, 'VerticalAlignment', 'top');
                end
                
                filename2 = sprintf('Combined_Concentration_FullTime_%dmin_%s.png', round(full_sim_time), timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax2, fullfile(folder, filename2), 'Resolution', 300);
                else
                    print(fig2, fullfile(folder, filename2), '-dpng', '-r300');
                end
                close(fig2);
                exportedCount = exportedCount + 1;
                
            catch ME
                fprintf('Full-time combined plot export failed: %s\n', ME.message);
            end
        end
        
        % Return number of successfully exported plots
        fprintf('Successfully exported %d combined concentration plots\n', exportedCount);
        
    catch ME
        fprintf('Combined concentration plots export failed: %s\n', ME.message);
    end
end
function updateConversionInfo(app)
    try
        infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: linear-gradient(135deg, #F2F2F7 0%, #E5E5EA 100%); border-radius: 8px; border: 2px solid #D1D1D6;">';
        infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #1D1D1F; text-align: center;">🎯 Conversion Steady-State Analysis</h3>'];
        
        % Create two-column layout using CSS
        infoHTML = [infoHTML '<div style="display: flex; gap: 20px;">'];
        
        % Left column: Product Formation Steady States
        infoHTML = [infoHTML '<div style="flex: 1; background: white; padding: 8px; border-radius: 6px; border: 1px solid #D1D1D6;">'];
        infoHTML = [infoHTML '<h4 style="margin: 0 0 8px 0; color: #007AFF; font-size: 12px;">📈 PRODUCT FORMATION STEADY-STATE</h4>'];
        
        % Check optimized product steady states
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            conv_NaOAc_opt = (app.SimResults.C_NaOAc / app.Cin) * 100;
            conv_EtOH_opt = (app.SimResults.C_EtOH / app.Cin) * 100;
            
            [ss_time_NaOAc_opt, ss_conv_NaOAc_opt] = app.detectSteadyStateHelper(app.SimResults.t, conv_NaOAc_opt);
            [ss_time_EtOH_opt, ss_conv_EtOH_opt] = app.detectSteadyStateHelper(app.SimResults.t, conv_EtOH_opt);
            
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 10px;"><strong style="color: #FF9500;">🔶 Optimized Results:</strong></p>'];
            if ~isnan(ss_time_NaOAc_opt)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ NaOAc: t=%.1f min, Conv=%.1f%%</p>', ss_time_NaOAc_opt, ss_conv_NaOAc_opt)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ NaOAc: No steady-state detected</p>'];
            end
            
            if ~isnan(ss_time_EtOH_opt)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ EtOH: t=%.1f min, Conv=%.1f%%</p>', ss_time_EtOH_opt, ss_conv_EtOH_opt)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ EtOH: No steady-state detected</p>'];
            end
        end
        
        % Check manual product steady states
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            conv_NaOAc_manual = (app.ManualSimResults.C_NaOAc / app.Cin) * 100;
            conv_EtOH_manual = (app.ManualSimResults.C_EtOH / app.Cin) * 100;
            
            [ss_time_NaOAc_manual, ss_conv_NaOAc_manual] = app.detectSteadyStateHelper(app.ManualSimResults.t, conv_NaOAc_manual);
            [ss_time_EtOH_manual, ss_conv_EtOH_manual] = app.detectSteadyStateHelper(app.ManualSimResults.t, conv_EtOH_manual);
            
            infoHTML = [infoHTML '<p style="margin: 4px 0 2px 0; font-size: 10px;"><strong style="color: #9366D6;">🔷 Manual Results:</strong></p>'];
            if ~isnan(ss_time_NaOAc_manual)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ NaOAc: t=%.1f min, Conv=%.1f%%</p>', ss_time_NaOAc_manual, ss_conv_NaOAc_manual)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ NaOAc: No steady-state detected</p>'];
            end
            
            if ~isnan(ss_time_EtOH_manual)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ EtOH: t=%.1f min, Conv=%.1f%%</p>', ss_time_EtOH_manual, ss_conv_EtOH_manual)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ EtOH: No steady-state detected</p>'];
            end
        end
        
        infoHTML = [infoHTML '</div>'];
        
        % Right column: Unreacted Reactant Steady States
        infoHTML = [infoHTML '<div style="flex: 1; background: white; padding: 8px; border-radius: 6px; border: 1px solid #D1D1D6;">'];
        infoHTML = [infoHTML '<h4 style="margin: 0 0 8px 0; color: #FF3B30; font-size: 12px;">📉 UNREACTED REACTANT STEADY-STATE</h4>'];
        
        % Check optimized reactant steady states
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            conv_NaOH_opt = ((app.Cin - app.SimResults.C_NaOH) / app.Cin) * 100;
            conv_EtOAc_opt = ((app.Cin - app.SimResults.C_EtOAc) / app.Cin) * 100;
            
            [ss_time_NaOH_opt, ss_conv_NaOH_opt] = app.detectSteadyStateHelper(app.SimResults.t, conv_NaOH_opt);
            [ss_time_EtOAc_opt, ss_conv_EtOAc_opt] = app.detectSteadyStateHelper(app.SimResults.t, conv_EtOAc_opt);
            
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 10px;"><strong style="color: #FF9500;">🔶 Optimized Results:</strong></p>'];
            if ~isnan(ss_time_NaOH_opt)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ NaOH: t=%.1f min, Conv=%.1f%%</p>', ss_time_NaOH_opt, ss_conv_NaOH_opt)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ NaOH: No steady-state detected</p>'];
            end
            
            if ~isnan(ss_time_EtOAc_opt)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ EtOAc: t=%.1f min, Conv=%.1f%%</p>', ss_time_EtOAc_opt, ss_conv_EtOAc_opt)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ EtOAc: No steady-state detected</p>'];
            end
        end
        
        % Check manual reactant steady states
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            conv_NaOH_manual = ((app.Cin - app.ManualSimResults.C_NaOH) / app.Cin) * 100;
            conv_EtOAc_manual = ((app.Cin - app.ManualSimResults.C_EtOAc) / app.Cin) * 100;
            
            [ss_time_NaOH_manual, ss_conv_NaOH_manual] = app.detectSteadyStateHelper(app.ManualSimResults.t, conv_NaOH_manual);
            [ss_time_EtOAc_manual, ss_conv_EtOAc_manual] = app.detectSteadyStateHelper(app.ManualSimResults.t, conv_EtOAc_manual);
            
            infoHTML = [infoHTML '<p style="margin: 4px 0 2px 0; font-size: 10px;"><strong style="color: #9366D6;">🔷 Manual Results:</strong></p>'];
            if ~isnan(ss_time_NaOH_manual)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ NaOH: t=%.1f min, Conv=%.1f%%</p>', ss_time_NaOH_manual, ss_conv_NaOH_manual)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ NaOH: No steady-state detected</p>'];
            end
            
            if ~isnan(ss_time_EtOAc_manual)
                infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #34C759;">✓ EtOAc: t=%.1f min, Conv=%.1f%%</p>', ss_time_EtOAc_manual, ss_conv_EtOAc_manual)];
            else
                infoHTML = [infoHTML '<p style="margin: 1px 0; font-size: 9px; color: #FF3B30;">✗ EtOAc: No steady-state detected</p>'];
            end
        end
        
        infoHTML = [infoHTML '</div>'];
        infoHTML = [infoHTML '</div>'];
        
        % Add summary note
        if isempty(app.SimResults) && isempty(app.ManualSimResults)
            infoHTML = [infoHTML '<p style="text-align: center; margin: 10px 0 0 0; font-size: 10px; color: #8E8E93; font-style: italic;">Run simulations to see steady-state analysis</p>'];
        else
            infoHTML = [infoHTML '<p style="text-align: center; margin: 8px 0 0 0; font-size: 9px; color: #8E8E93;">Steady-state threshold: <0.1% relative change</p>'];
        end
        
        infoHTML = [infoHTML '</div>'];
        app.ConversionInfo.HTMLSource = infoHTML;
        
    catch
        app.ConversionInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No conversion data available.</div>';
    end
end

% Helper function for steady-state detection (can be reused)
function [ss_time, ss_value] = detectSteadyStateHelper(app, time_data, conversion_data)
    try
        ss_time = NaN;
        ss_value = NaN;
        
        if length(time_data) < 10
            return;
        end
        
        % Steady-state threshold: relative change < 0.1% for conversion
        ss_threshold = 0.0005;; % 0.1% relative change threshold
        
        % Check from 10th point onwards
        for i = 10:length(conversion_data)
            if i > 1
                rel_change = abs(conversion_data(i) - conversion_data(i-1)) / ...
                             max(eps, abs(conversion_data(i-1))) * 100;
                
                if rel_change < ss_threshold
                    ss_time = time_data(i);
                    ss_value = conversion_data(i);
                    break;
                end
            end
        end
        
    catch
        ss_time = NaN;
        ss_value = NaN;
    end
end
        
        function setupDesignSystem(app)
            % Apply consistent design theme
            app.UIFigure.Color = app.BackgroundColor;
            
            % CORRECTED: Calculate overall heat transfer coefficient with CONSISTENT UNITS
            U = 1 / (1/app.h_reactor + app.delta_wall/app.k_wall + 1/app.h_coolant); % [W/(m²·K)]
            U = U * 60;                % Convert to [J/(min·m²·K)] - consistent with minute time scale
            app.UA = U * app.a;        % [J/(min·K)] - ALL in minutes now
            
            % Configure panels
            panels = [app.LeftPanel, app.RightPanel, app.DataPanel, app.ParamsPanel, ...
                     app.ActionsPanel, app.DisplayPanel, app.ExportPanel, app.StatusPanel];
            for i = 1:length(panels)
                panels(i).BackgroundColor = app.SurfaceColor;
                panels(i).BorderType = 'none';
            end
            
            % Configure buttons
            app.UploadButton.BackgroundColor = app.AccentColor;
            app.UploadButton.FontColor = app.SurfaceColor;
            
            app.ManualKButton.BackgroundColor = [147 112 219]/255; % Medium Purple
            app.ManualKButton.FontColor = app.SurfaceColor;
            
            app.OptimizeButton.BackgroundColor = app.WarningColor;
            app.OptimizeButton.FontColor = app.SurfaceColor;
            
            app.OptimizeUltraButton.BackgroundColor = [255 20 147]/255; % Deep Pink for ultra precision
            app.OptimizeUltraButton.FontColor = app.SurfaceColor;
            
            app.ClearAllButton.BackgroundColor = app.ErrorColor; % Red color for destructive action
            app.ClearAllButton.FontColor = app.SurfaceColor;
            
            app.SaveButton.BackgroundColor = app.AccentColor;
            app.SaveButton.FontColor = app.SurfaceColor;
            
            app.ExportDataButton.BackgroundColor = app.AccentColor;
            app.ExportDataButton.FontColor = app.SurfaceColor;
            
            app.ExportPlotsButton.BackgroundColor = [88 86 214]/255; % Indigo
            app.ExportPlotsButton.FontColor = app.SurfaceColor;
            
            % Configure Quick Sim button
            if isfield(app, 'RunQuickSimButton') && ~isempty(app.RunQuickSimButton)
                app.RunQuickSimButton.BackgroundColor = [34 139 34]/255; % Forest Green
                app.RunQuickSimButton.FontColor = app.SurfaceColor;
            end
            
            % Configure switches
            app.JacketSwitch.FontColor = app.TextColor;
            app.ProductsSwitch.FontColor = app.TextColor;

if isfield(app, 'ValidateDataButton') && ~isempty(app.ValidateDataButton)
    app.ValidateDataButton.BackgroundColor = [34 139 34]/255; % Forest Green
    app.ValidateDataButton.FontColor = app.SurfaceColor;
end     

        end
        
        function animateLogo(app)
            % Animated reactor logo
            html = [...
                '<div style="display: flex; justify-content: center; align-items: center; height: 100%; background: white; border-radius: 12px;">' ...
                '<div class="reactor-container" style="position: relative; width: 80px; height: 80px;">' ...
                '<svg width="80" height="80" viewBox="0 0 80 80">' ...
                '<defs>' ...
                '<linearGradient id="reactorGradient" x1="0%" y1="0%" x2="100%" y2="100%">' ...
                '<stop offset="0%" style="stop-color:#007AFF;stop-opacity:1" />' ...
                '<stop offset="100%" style="stop-color:#5856D6;stop-opacity:1" />' ...
                '</linearGradient>' ...
                '</defs>' ...
                '<circle cx="40" cy="40" r="30" fill="none" stroke="url(#reactorGradient)" stroke-width="3" opacity="0.8">' ...
                '<animate attributeName="r" values="30;35;30" dur="2s" repeatCount="indefinite"/>' ...
                '</circle>' ...
                '<path d="M40,10 L40,70 M10,40 L70,40" stroke="url(#reactorGradient)" stroke-width="2" opacity="0.6"/>' ...
                '<circle cx="40" cy="40" r="10" fill="url(#reactorGradient)">' ...
                '<animate attributeName="r" values="10;12;10" dur="2s" repeatCount="indefinite"/>' ...
                '</circle>' ...
                '<circle cx="20" cy="20" r="3" fill="#007AFF" opacity="0.6">' ...
                '<animateTransform attributeName="transform" type="rotate" from="0 40 40" to="360 40 40" dur="8s" repeatCount="indefinite"/>' ...
                '</circle>' ...
                '<circle cx="60" cy="60" r="3" fill="#5856D6" opacity="0.6">' ...
                '<animateTransform attributeName="transform" type="rotate" from="0 40 40" to="-360 40 40" dur="6s" repeatCount="indefinite"/>' ...
                '</circle>' ...
                '</svg>' ...
                '<div style="text-align: center; margin-top: 5px; font-family: -apple-system, BlinkMacSystemFont, ''Segoe UI'', Roboto, Helvetica, Arial; font-size: 10px; color: #8E8E93;">CSTR Simulator</div>' ...
                '</div>' ...
                '</div>'];
            
            app.LogoHTML.HTMLSource = html;
        end
        
        function showLoadingDialog(app, title, message)
            % Show loading dialog with progress
            try
                if ~isempty(app.LoadingDialog) && isvalid(app.LoadingDialog)
                    close(app.LoadingDialog);
                end
                app.LoadingDialog = uiprogressdlg(app.UIFigure, 'Title', title, ...
                                                 'Message', message, 'Indeterminate', 'off', ...
                                                 'Value', 0);
            catch
                % Fallback to status updates if dialog fails
                app.updateStatus(message, 'info');
            end
        end
        
        function updateLoadingDialog(app, value, message)
            % Update loading dialog progress
            try
                if ~isempty(app.LoadingDialog) && isvalid(app.LoadingDialog)
                    app.LoadingDialog.Value = value;
                    if nargin > 2
                        app.LoadingDialog.Message = message;
                    end
                end
            catch
                % Fallback to status updates
                app.updateProgress(value * 100);
                if nargin > 2
                    app.updateStatus(message, 'info');
                end
            end
        end
        
        function closeLoadingDialog(app)
            % Close loading dialog
            try
                if ~isempty(app.LoadingDialog) && isvalid(app.LoadingDialog)
                    close(app.LoadingDialog);
                    app.LoadingDialog = [];
                end
            catch
                % Silent fail
            end
        end
        
        function updateStatus(app, message, type)
            % Update status with appropriate styling
            try
                switch type
                    case 'success'
                        color = app.SuccessColor;
                        icon = '✓';
                    case 'error'
                        color = app.ErrorColor;
                        icon = '✗';
                    case 'warning'
                        color = app.WarningColor;
                        icon = '⚠';
                    case 'info'
                        color = app.AccentColor;
                        icon = 'ℹ';
                    otherwise
                        color = app.TextColor;
                        icon = '•';
                end
                
                app.StatusLabel.Text = sprintf('%s %s', icon, message);
                app.StatusLabel.FontColor = color;
                
                % Auto-hide after 5 seconds
                pause(0.1);
                if ~isempty(app.LoadingTimer) && isvalid(app.LoadingTimer)
                    stop(app.LoadingTimer);
                    delete(app.LoadingTimer);
                end
                app.LoadingTimer = timer('TimerFcn', @(~,~) set(app.StatusLabel, 'Text', ''), ...
                                        'StartDelay', 5, 'ExecutionMode', 'singleShot');
                start(app.LoadingTimer);
            catch
                % Silent fail
            end
        end
        
        function updateProgress(app, value)
            % Update progress label instead of gauge
            try
                if value >= 100
                    app.ProgressLabel.Text = sprintf('Progress: %d%% Complete', round(value));
                    app.ProgressLabel.FontColor = app.SuccessColor;
                else
                    app.ProgressLabel.Text = sprintf('Progress: %d%%', round(value));
                    app.ProgressLabel.FontColor = app.AccentColor;
                end
            catch
                % Silent fail
            end
        end
        
        function flashSuccess(app, component)
            % Flash component with success color
            try
                originalColor = component.BackgroundColor;
                for i = 1:2
                    component.BackgroundColor = app.SuccessColor;
                    pause(0.1);
                    component.BackgroundColor = originalColor;
                    pause(0.1);
                end
            catch
                % Silent fail
            end
        end
        
        function enableControls(app, enable)
            % Enable/disable all input controls
            try
                controls = [app.UploadButton, app.TempDropDown, app.VolumeField, ...
                           app.FlowrateField, app.CinField, app.SimTimeField, ...
                           app.ManualKButton, app.OptimizeButton, app.OptimizeUltraButton, ...
                           app.ClearAllButton, ...
                           app.SaveButton, app.ExportDataButton, app.ExportPlotsButton];
                
                for i = 1:length(controls)
                    if isvalid(controls(i))
                        controls(i).Enable = enable;
                    end
                end
            catch
                % Silent fail
            end
        end
        
        % Clear all data and reset the application
        function clearAllData(app)
            try
                app.showLoadingDialog('Clearing Data', 'Resetting application...');
                
                % Clear all data tables
                app.Data = table();
                app.ComparisonTableData = table();
                app.ComparisonTableDataManual = table();
                
                % Clear all simulation results
                app.SimResults = struct();
                app.ManualSimResults = struct();
                app.QuickSimResults = struct();
                
                % Clear all convergence data
                app.ConvergenceData = struct();
                app.ManualConvergenceData = struct();
                app.QuickSimConvergenceData = struct();
                
                % Reset all convergence status flags
                app.IsConverged = false;
                app.IsManualConverged = false;
                app.IsQuickSimConverged = false;
                
                % Reset separate convergence tracking
                app.ManualReactantsConverged = false;
                app.ManualProductsConverged = false;
                app.ManualTemperatureConverged = false;
                app.OptimizedReactantsConverged = false;
                app.OptimizedProductsConverged = false;
                app.OptimizedTemperatureConverged = false;
                app.QuickSimReactantsConverged = false;
                app.QuickSimProductsConverged = false;
                app.QuickSimTemperatureConverged = false;
                
                % Reset rate constants and R-squared values
                app.k_fit = NaN;
                app.k_manual = NaN;
                app.RSquared = NaN;
                app.RSquaredManual = NaN;
                
                % Clear optimization results
                app.OptimizationResults = struct();
                
                % Reset file label
                app.FileLabel.Text = 'No file selected';
                
                app.updateLoadingDialog(0.3, 'Clearing plots...');
                
                % Clear all plots
                cla(app.UIAxes1);
                cla(app.UIAxes2);
                cla(app.UIAxes3);
                cla(app.UIAxesQuickConc);
                cla(app.UIAxesQuickTemp);
                cla(app.UIAxesQuickProducts);
                cla(app.UIAxesConversionProduct);
cla(app.UIAxesConversionReactant);
                
                % Reset plot titles
                title(app.UIAxes1, 'Concentration Profiles with Differential Convergence');
                xlabel(app.UIAxes1, 'Time (min)');
                ylabel(app.UIAxes1, 'Concentration (mol/L)');
                
                title(app.UIAxes2, 'Temperature Profile with Differential Convergence');
                xlabel(app.UIAxes2, 'Time (min)');
                ylabel(app.UIAxes2, 'Temperature (°C)');
                
                title(app.UIAxes3, 'Product Formation with Differential Convergence');
                xlabel(app.UIAxes3, 'Time (min)');
                ylabel(app.UIAxes3, 'Concentration (mol/L)');
                
                title(app.UIAxesQuickConc, 'Concentration Profiles with Convergence');
                xlabel(app.UIAxesQuickConc, 'Time (min)');
                ylabel(app.UIAxesQuickConc, 'Concentration (mol/L)');
                
                title(app.UIAxesQuickTemp, 'Temperature Profiles');
                xlabel(app.UIAxesQuickTemp, 'Time (min)');
                ylabel(app.UIAxesQuickTemp, 'Temperature (°C)');
                
                title(app.UIAxesQuickProducts, 'Product Formation with Convergence');
                xlabel(app.UIAxesQuickProducts, 'Time (min)');
                ylabel(app.UIAxesQuickProducts, 'Concentration (mol/L)');
                
                app.updateLoadingDialog(0.6, 'Clearing tables...');
                
                % Clear all tables
                app.ComparisonTable.Data = [];
                app.ComparisonTable.ColumnName = {};
                app.ConcentrationTable.Data = [];
                app.ConcentrationTable.ColumnName = {'Time', 'NaOH', 'EtOAc', 'NaOAc', 'EtOH'};
                app.TemperatureTable.Data = [];
                app.TemperatureTable.ColumnName = {'Time', 'Reactor', 'Jacket'};
                
                % Clear analysis text
                app.ConvergenceText.Value = 'No analysis data available. Run a simulation first.';
                
                % Clear info panels
                app.ConcentrationInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
                app.TemperatureInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
                app.ProductsInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
                app.QuickSimInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
                
                % Clear Quick Sim error labels
                app.MaxErrorLabel.Text = 'Max Error: N/A';
                app.MinErrorLabel.Text = 'Min Error: N/A';
                app.MeanErrorLabel.Text = 'Mean Error: N/A';
                
                % Reset parameters to defaults (optional)
                app.SelectedTemp = 30;
                app.TempDropDown.Value = '30';
                app.VolumeField.Value = 2.5;
                app.FlowrateField.Value = 0.18;
                app.CinField.Value = 0.05;
                app.SimTimeField.Value = 30;
                
                % Reset Quick Sim fields
                app.CinQuickField.Value = 0.05;
                app.FlowrateQuickField.Value = 0.18;
                app.TempQuickDropDown.Value = '30';
                app.SimTimeQuickField.Value = 30;
                app.RateConstQuickField.Value = 0.5;
                
                app.updateLoadingDialog(1.0, 'Clear complete!');
                app.closeLoadingDialog();
                
                app.updateStatus('🗑️ All data and plots cleared successfully! Ready for new simulation.', 'success');
                app.flashSuccess(app.ClearAllButton);
                
            catch ME
                app.closeLoadingDialog();
                app.updateStatus(['Clear failed: ' ME.message], 'error');
            end
        end
        
        % CORRECTED: Enhanced ODE function with CONSISTENT UNITS
        function dydt = odefun_individual(app, t, y, k, V0_EtOAc, V1_NaOH, V, Cin, To)
    try
        % Extract variables
        CA = y(1);      % NaOH concentration [mol/L]
        CB = y(2);      % EtOAc concentration [mol/L]
        CC = y(3);      % NaOAc concentration [mol/L]
        CD = y(4);      % EtOH concentration [mol/L]
        T_reactor = y(5);  % Reactor temperature [°C]
        Tj = y(6);      % Jacket temperature [°C]
        
        % Total flowrate [L/min]
        V0_total = V0_EtOAc + V1_NaOH;
        
        % Calculate mixture properties
        [p_avg, Cp_avg] = app.calculate_mixture_properties(CA, CB, CC, CD);
        
        % Reaction rate (second-order) [mol/(L·min)]
        r = k * CA * CB;
        
        % Temperature-dependent heat of reaction [J/mol]
        lambda_T = app.heat_of_reaction(T_reactor + 273.15);
        
        % Heat generation from reaction [J/min]
        Q_rxn = -lambda_T * r * V;
        
        % TEMPERATURE-RESPONSIVE ADAPTIVE COOLING SYSTEM
        if app.IsothermalMode
            % ISOTHERMAL MODE: Temperature-responsive adaptive cooling
            
            % 1. Calculate temperature-dependent parameters
            temp_error = T_reactor - To;
            temp_above_setpoint = max(0, T_reactor - To); % Only positive deviations matter for cooling
            
            % 2. BASE COOLING FLOWRATE - Temperature dependent
            % Establish base flowrate that increases with temperature
            temp_setpoint_celsius = To;
            if temp_setpoint_celsius <= 30
                base_flowrate = 0.8;  % Low base for 30°C
            elseif temp_setpoint_celsius <= 40
                base_flowrate = 1.2;  % Medium base for 30-40°C  
            elseif temp_setpoint_celsius <= 50
                base_flowrate = 1.8;  % Higher base for 40-50°C
            elseif temp_setpoint_celsius <= 60
                base_flowrate = 2.5;  % High base for 50-60°C
            else
                base_flowrate = 3.5;  % Very high base for >60°C
            end
            
            % 3. TEMPERATURE DIFFERENCE SCALING
            temp_difference_reactor_coolant = max(T_reactor - app.Tj0, 0.5);
            
            % Progressive scaling based on temperature difference
            if temp_difference_reactor_coolant > 15
                temp_diff_scaling = 3.0;  % Very high scaling for large differences
            elseif temp_difference_reactor_coolant > 10
                temp_diff_scaling = 2.5;  % High scaling
            elseif temp_difference_reactor_coolant > 7
                temp_diff_scaling = 2.0;  % Medium-high scaling
            elseif temp_difference_reactor_coolant > 5
                temp_diff_scaling = 1.5;  % Medium scaling
            elseif temp_difference_reactor_coolant > 3
                temp_diff_scaling = 1.2;  % Slight scaling
            else
                temp_diff_scaling = 1.0;  % No scaling for small differences
            end
            
            % 4. REACTION HEAT SCALING
            Q_rxn_magnitude = abs(Q_rxn);
            if Q_rxn_magnitude > 5000
                heat_scaling = 2.0;     % High heat generation
            elseif Q_rxn_magnitude > 2000
                heat_scaling = 1.5;     % Medium heat generation
            elseif Q_rxn_magnitude > 500
                heat_scaling = 1.2;     % Low heat generation
            else
                heat_scaling = 1.0;     % Minimal heat generation
            end
            
            % 5. ENHANCED PI CONTROLLER
            % Adjust PI gains based on temperature level
            if temp_setpoint_celsius > 50
                Kp = 2.5;  % Higher gain for high temperatures
                Ki = 0.3;
            elseif temp_setpoint_celsius > 40
                Kp = 2.0;  % Medium gain
                Ki = 0.2;
            else
                Kp = 1.5;  % Lower gain for lower temperatures
                Ki = 0.1;
            end
            
            % Integral windup prevention with temperature-dependent limits
            if abs(temp_error) < app.isothermal_tolerance / 3
                app.integral_error = app.integral_error * 0.8; % Reduce windup
            else
                app.integral_error = app.integral_error + temp_error;
            end
            
            % Limit integral term based on temperature
            max_integral = 5 + temp_setpoint_celsius / 10; % Higher limit for higher temps
            app.integral_error = max(-max_integral, min(max_integral, app.integral_error));
            
            % PI controller output
            control_output = Kp * temp_error + Ki * app.integral_error;
            
            % 6. CALCULATE REQUIRED COOLING FLOWRATE
            % Start with base flowrate for the temperature
            Fj_required = base_flowrate;
            
            % Apply temperature difference scaling
            Fj_required = Fj_required * temp_diff_scaling;
            
            % Apply heat generation scaling
            Fj_required = Fj_required * heat_scaling;
            
            % Apply PI controller correction
            if control_output > 0  % Need more cooling
                control_scaling = 1 + (control_output / 10); % Convert control output to scaling
                Fj_required = Fj_required * min(control_scaling, 2.0); % Limit max scaling
            end
            
            % 7. EMERGENCY COOLING for very high temperatures
            if T_reactor > (To + 5)  % Emergency cooling threshold
                emergency_scaling = 1 + (T_reactor - To - 5) / 5; % Escalating emergency response
                Fj_required = Fj_required * min(emergency_scaling, 3.0);
            end
            
            % 8. APPLY FLOWRATE LIMITS
            if Fj_required > app.Fj_max
                app.Fj_adaptive = app.Fj_max;
                cooling_at_limit = true;
            elseif Fj_required < app.Fj_min
                app.Fj_adaptive = app.Fj_min;
                cooling_at_limit = false;
            else
                app.Fj_adaptive = Fj_required;
                cooling_at_limit = false;
            end
            
            % 9. CALCULATE COOLING PERFORMANCE METRICS
            Q_total_required = abs(Q_rxn) + abs(control_output * 100); % Approximate total cooling need
            app.required_cooling_capacity = Q_total_required;
            app.max_cooling_capacity = app.Fj_max * app.pj/1000 * app.Cpj * temp_difference_reactor_coolant;
            
            if app.max_cooling_capacity > 0
                app.cooling_adequacy_factor = app.max_cooling_capacity / app.required_cooling_capacity;
            else
                app.cooling_adequacy_factor = 0;
            end
            
            % Heat transfer efficiency
            if cooling_at_limit && temp_error > app.isothermal_tolerance
                heat_transfer_efficiency = 0.85; % Reduced efficiency at limits
            else
                heat_transfer_efficiency = 1.0;
            end
            
        else
            % NON-ISOTHERMAL MODE: NO COOLING SYSTEM
            app.Fj_adaptive = 0;
            heat_transfer_efficiency = 0;
            temp_difference_reactor_coolant = T_reactor - app.Tj0;
            Q_total_required = abs(Q_rxn);
            
            app.required_cooling_capacity = abs(Q_rxn);
            app.max_cooling_capacity = 0;
            app.cooling_adequacy_factor = 0;
        end
        
        % MASS BALANCE EQUATIONS (unchanged)
        dCA_dt = (V1_NaOH/V) * (Cin - CA) - r;
        dCB_dt = (V0_EtOAc/V) * (Cin - CB) - r;
        dCC_dt = (V0_total/V) * (0 - CC) + r;
        dCD_dt = (V0_total/V) * (0 - CD) + r;
        
        % HEAT TRANSFER CALCULATION
        if app.IsothermalMode
            temp_difference_jacket = max(T_reactor - Tj, 0.1);
            Q_jacket_theoretical = app.UA * temp_difference_jacket;
            Q_jacket = Q_jacket_theoretical * heat_transfer_efficiency;
            
            % Limit heat transfer to physical maximum
            max_possible_cooling = app.Fj_adaptive * app.pj/1000 * app.Cpj * temp_difference_jacket;
            Q_jacket = min(Q_jacket, max_possible_cooling);
        else
            Q_jacket = 0;
        end
        
        % Convert mixture density
        p_avg_kg_L = p_avg / 1000;
        
        % REACTOR ENERGY BALANCE
        if app.IsothermalMode
            % Temperature-dependent damping
            if temp_setpoint_celsius > 50
                temperature_damping = 12.0; % Strong damping for high temps
            elseif temp_setpoint_celsius > 40
                temperature_damping = 8.0;  % Medium damping
            else
                temperature_damping = 5.0;  % Light damping for low temps
            end
            
            dT_reactor_dt = (V0_total/V) * (To - T_reactor) + ...
                           (Q_rxn) / (p_avg_kg_L * Cp_avg * V) - ...
                           (Q_jacket) / (p_avg_kg_L * Cp_avg * V) - ...
                           temperature_damping * (T_reactor - To);
        else
            dT_reactor_dt = (V0_total/V) * (To - T_reactor) + ...
                           (Q_rxn) / (p_avg_kg_L * Cp_avg * V);
        end
        
        % JACKET ENERGY BALANCE
        if app.IsothermalMode
            pj_kg_L = app.pj / 1000;
            dTj_dt = (app.Fj_adaptive / app.Vj) * (app.Tj0 - Tj) + ...
                    (Q_jacket) / (pj_kg_L * app.Cpj * app.Vj);
        else
            dTj_dt = 0;
        end
        
        % Store enhanced cooling history
        temp_deviation = abs(T_reactor - To);
        app.storeCoolingHistoryEnhanced(t, app.Fj_adaptive, app.required_cooling_capacity, ...
                                       app.cooling_adequacy_factor, temp_deviation, ...
                                       temp_difference_reactor_coolant, Q_total_required);
        
        % Return derivatives
        dydt = [dCA_dt; dCB_dt; dCC_dt; dCD_dt; dT_reactor_dt; dTj_dt];
        
    catch
        dydt = zeros(6,1);
    end
end
        
        % Temperature-dependent heat of reaction function
        function lambda_T = heat_of_reaction(app, T)
            % T in Kelvin, returns J/mol
            lambda_298 = -3.692e7;     % [J/mol] at 298.15 K
            Tref = 298.15;             % [K]
            % Heat capacity difference [J/(mol·K)]
            dCp = (229e3 + 158.8e3) - (59.52e3 + 168.94e3);  % Products - Reactants
            lambda_T = lambda_298 + dCp * (T - Tref);
        end
        
        % Enhanced mixture properties calculation
        function [p_avg, Cp_avg] = calculate_mixture_properties(app, CA, CB, CC, CD)
            % Component properties
            % Densities [kg/m³]
            rho_NaOH = 2100; rho_EtOAc = 900; rho_NaOAc = 1500; rho_EtOH = 790;
            
            % Molar masses [kg/mol]
            M_NaOH = 0.039997; M_EtOAc = 0.08811; M_NaOAc = 0.08203; M_EtOH = 0.04607;
            
            % Molar heat capacities [J/(mol·K)]
            Cp_NaOH = 59.52e3; Cp_EtOAc = 168.94e3; Cp_NaOAc = 229e3; Cp_EtOH = 158.8e3;
            
            % Calculate total moles and mole fractions
            total_moles = max(CA + CB + CC + CD, 1e-12);
            yA = CA / total_moles; yB = CB / total_moles; 
            yC = CC / total_moles; yD = CD / total_moles;
            
            % Convert to mass fractions
            avg_MW = yA*M_NaOH + yB*M_EtOAc + yC*M_NaOAc + yD*M_EtOH;
            xA = (yA * M_NaOH) / avg_MW; xB = (yB * M_EtOAc) / avg_MW;
            xC = (yC * M_NaOAc) / avg_MW; xD = (yD * M_EtOH) / avg_MW;
            
            % Mixture density (Kay's rule) [kg/m³]
            p_avg = 1 / (xA/rho_NaOH + xB/rho_EtOAc + xC/rho_NaOAc + xD/rho_EtOH);
            
            % Mixture heat capacity (mole fraction basis) [J/(mol·K)]
            Cp_avg_molar = yA*Cp_NaOH + yB*Cp_EtOAc + yC*Cp_NaOAc + yD*Cp_EtOH;
            
            % Convert to mass basis [J/(kg·K)]
            Cp_avg = Cp_avg_molar / avg_MW;
        end
        
        % Cost function for optimization
        function cost = calculateCost(app, k, t_exp, C_exp, y0, V0_EtOAc, V1_NaOH, V, Cin, To)
            try
                if k <= 0 || ~isfinite(k)
                    cost = 1e6;
                    return;
                end
                
                % Solve ODE with current k
                options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
                [~, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,k,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                  t_exp, y0, options);
                
                C_sim = y_sim(:,1);  % NaOH concentration
                
                % Calculate cost (sum of squared errors)
                cost = sum((C_exp - C_sim).^2);
                
            catch
                cost = 1e6;
            end
        end
        
        % Enhanced cost function for ultra-precise optimization
        function max_error = calculateMaxErrorWithK(app, params)
            try
                % Extract parameters
                V0_total = params(1);
                k = params(2);
                
                % Equal split
                V0_EtOAc = V0_total / 2;
                V1_NaOH = V0_total / 2;
                
                % Quick calculation for optimization
                t_exp = app.Data.Time;
                C_exp = app.Data.Concentration;
                T0 = app.SelectedTemp;
                To = T0;
                V = app.VolumeField.Value;
                Cin = app.CinField.Value;
                
                % Initial conditions [CA, CB, CC, CD, T_reactor, Tj]
                y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
                
                if k <= 0 || ~isfinite(k)
                    max_error = 1e6;
                    return;
                end
                
                % Use tighter tolerances for ultra-precise simulation
                options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
                [~, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,k,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                  t_exp, y0, options);
                
                C_sim = y_sim(:,1);
                
                % Calculate relative errors
                rel_errors = abs(C_exp - C_sim) ./ C_exp * 100;
                rel_errors(C_exp == 0) = 0;
                
                max_error = max(rel_errors);
                
            catch
                max_error = 1e6;
            end
        end
        
        function max_error = calculateMaxErrorForOptimization(app, params, t_exp, C_exp, y0, V, Cin, To)
            try
                % Extract parameters
                V0_total = params(1);
                k = params(2);
                
                % Equal split for reactants
                V0_EtOAc = V0_total / 2;
                V1_NaOH = V0_total / 2;
                
                if k <= 0 || ~isfinite(k) || V0_total <= 0 || ~isfinite(V0_total)
                    max_error = 1e6;
                    return;
                end
                
                % Solve ODE with current parameters
                options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
                [~, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,k,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                  t_exp, y0, options);
                
                C_sim = y_sim(:,1);  % NaOH concentration
                
                % Calculate relative errors
                rel_errors = abs(C_exp - C_sim) ./ C_exp * 100;
                rel_errors(C_exp == 0) = 0;
                
                max_error = max(rel_errors);
                
                % Add penalty if error is too high to guide optimization toward 1% target
                if max_error > 10
                    max_error = max_error * 2; % Penalize high errors more
                end
                
            catch
                max_error = 1e6;
            end
        end
        
        % ENHANCED: Improved convergence detection with different criteria for manual vs optimized vs quick sim   
function detectConvergence(app, method)
    try
        if strcmp(method, 'manual')
            if isempty(app.ManualSimResults) || ~isfield(app.ManualSimResults, 't')
                return;
            end
            simResults = app.ManualSimResults;
            convergenceData = struct();
            
            % Use manual-specific thresholds (now equal to optimized)
            reactant_threshold = 0.001;  % 5% relative change for reactants
            product_threshold = 0.000005;   % 5% relative change for products
            temp_threshold = 0.0005;      % 2% relative change for temperature
            
        elseif strcmp(method, 'quicksim')
            if isempty(app.QuickSimResults) || ~isfield(app.QuickSimResults, 't')
                return;
            end
            simResults = app.QuickSimResults;
            convergenceData = struct();
            
            % Use quick sim-specific thresholds (medium)
            reactant_threshold = 0.001;  % 7% relative change for reactants
            product_threshold = 0.000005;   % 7% relative change for products
            temp_threshold = 0.0005;      % 3% relative change for temperature
            
        else % optimized
            if isempty(app.SimResults) || ~isfield(app.SimResults, 't')
                return;
            end
            simResults = app.SimResults;
            convergenceData = struct();
            
            % Use optimized-specific thresholds (same as manual now)
            reactant_threshold = 0.001;  % 5% relative change for reactants
            product_threshold  = 0.000005;   % 5% relative change for products
            temp_threshold = 0.0005;      % 2% relative change for temperature
        end
        
        t = simResults.t;
        num_steps = length(t);
        
        % Check if we have enough data points
        if num_steps < 10
            return; % Not enough data points for convergence analysis
        end
        
        % Initialize convergence indices (set to num_steps means not converged)
        naoh_convergence_index = num_steps;
        etoac_convergence_index = num_steps;
        naoac_convergence_index = num_steps;
        etoh_convergence_index = num_steps;
        temp_convergence_index = num_steps;
        
        % Track convergence status separately for reactants and products
        reactantsConverged = true;
        productsConverged = true;
        temperatureConverged = true;
        
        % Step-by-step convergence detection starting from step 2
        for i = 2:num_steps
            
            % Check NaOH convergence (reactant)
            if naoh_convergence_index == num_steps
                rel_change_NaOH = abs(simResults.C_NaOH(i) - simResults.C_NaOH(i-1)) / ...
                                  max(eps, abs(simResults.C_NaOH(i-1))) * 100;
                if rel_change_NaOH < reactant_threshold
                    naoh_convergence_index = i;
                    convergenceData.NaOH_time = t(i);
                    convergenceData.NaOH_conc = simResults.C_NaOH(i);
                    convergenceData.NaOH_converged = true;
                    convergenceData.NaOH_final = simResults.C_NaOH(end);
                end
            end
            
            % Check EtOAc convergence (reactant)
            if etoac_convergence_index == num_steps
                rel_change_EtOAc = abs(simResults.C_EtOAc(i) - simResults.C_EtOAc(i-1)) / ...
                                   max(eps, abs(simResults.C_EtOAc(i-1))) * 100;
                if rel_change_EtOAc < reactant_threshold
                    etoac_convergence_index = i;
                    convergenceData.EtOAc_time = t(i);
                    convergenceData.EtOAc_conc = simResults.C_EtOAc(i);
                    convergenceData.EtOAc_converged = true;
                    convergenceData.EtOAc_final = simResults.C_EtOAc(end);
                end
            end
            
            % Check NaOAc convergence (product)
            if naoac_convergence_index == num_steps
                rel_change_NaOAc = abs(simResults.C_NaOAc(i) - simResults.C_NaOAc(i-1)) / ...
                                   max(eps, abs(simResults.C_NaOAc(i-1))) * 100;
                if rel_change_NaOAc < product_threshold
                    naoac_convergence_index = i;
                    convergenceData.NaOAc_time = t(i);
                    convergenceData.NaOAc_conc = simResults.C_NaOAc(i);
                    convergenceData.NaOAc_converged = true;
                    convergenceData.NaOAc_final = simResults.C_NaOAc(end);
                end
            end
            
            % Check EtOH convergence (product)
            if etoh_convergence_index == num_steps
                rel_change_EtOH = abs(simResults.C_EtOH(i) - simResults.C_EtOH(i-1)) / ...
                                  max(eps, abs(simResults.C_EtOH(i-1))) * 100;
                if rel_change_EtOH < product_threshold
                    etoh_convergence_index = i;
                    convergenceData.EtOH_time = t(i);
                    convergenceData.EtOH_conc = simResults.C_EtOH(i);
                    convergenceData.EtOH_converged = true;
                    convergenceData.EtOH_final = simResults.C_EtOH(end);
                end
            end
            
            % Check Temperature convergence
            if temp_convergence_index == num_steps
                rel_change_T = abs(simResults.T_reactor(i) - simResults.T_reactor(i-1)) / ...
                              max(eps, abs(simResults.T_reactor(i-1))) * 100;
                if rel_change_T < temp_threshold
                    temp_convergence_index = i;
                    convergenceData.T_time = t(i);
                    convergenceData.T_value = simResults.T_reactor(i);
                    convergenceData.T_converged = true;
                    convergenceData.T_final = simResults.T_reactor(end);
                end
            end
        end
        
        % Set convergence flags based on whether indices were found
        if naoh_convergence_index == num_steps
            convergenceData.NaOH_converged = false;
            reactantsConverged = false;
        end
        
        if etoac_convergence_index == num_steps
            convergenceData.EtOAc_converged = false;
            reactantsConverged = false;
        end
        
        if naoac_convergence_index == num_steps
            convergenceData.NaOAc_converged = false;
            productsConverged = false;
        end
        
        if etoh_convergence_index == num_steps
            convergenceData.EtOH_converged = false;
            productsConverged = false;
        end
        
        if temp_convergence_index == num_steps
            convergenceData.T_converged = false;
            temperatureConverged = false;
        end
        
        % Overall convergence requires both reactants and products to converge
        isConverged = reactantsConverged && productsConverged;
        
        % Store convergence data and separate status tracking
        if strcmp(method, 'manual')
            app.ManualConvergenceData = convergenceData;
            app.IsManualConverged = isConverged;
            
            % Store separate convergence status for manual simulation
            app.ManualReactantsConverged = reactantsConverged;
            app.ManualProductsConverged = productsConverged;
            app.ManualTemperatureConverged = temperatureConverged;
            
        elseif strcmp(method, 'quicksim')
            app.QuickSimConvergenceData = convergenceData;
            app.IsQuickSimConverged = isConverged;
            
            % Store separate convergence status for quick sim
            app.QuickSimReactantsConverged = reactantsConverged;
            app.QuickSimProductsConverged = productsConverged;
            app.QuickSimTemperatureConverged = temperatureConverged;
            
            % To detectConvergence function after the 'quicksim' case:
% ENHANCED: Complete the ultra-precision case in detectConvergence function
% Replace the incomplete 'ultra' section in the detectConvergence function with this:

elseif strcmp(method, 'ultra')
    if isempty(app.UltraPrecisionResults) || ~isfield(app.UltraPrecisionResults, 't')
        return;
    end
    simResults = app.UltraPrecisionResults;
    convergenceData = struct();
    
    % Use ultra-specific thresholds (strictest)
    reactant_threshold = 0.001;  % 0.01% relative change for reactants
    product_threshold = 0.00005;   % 0.001% relative change for products  
    temp_threshold = 0.0001;      % 0.01% relative change for temperature
    
    % Store the convergence data and status
    app.UltraConvergenceData = convergenceData;
    
    % Continue with the same convergence detection logic as other methods
    t = simResults.t;
    num_steps = length(t);
    
    % Check if we have enough data points
    if num_steps < 10
        return;
    end
    
    % Initialize convergence indices
    naoh_convergence_index = num_steps;
    etoac_convergence_index = num_steps;
    naoac_convergence_index = num_steps;
    etoh_convergence_index = num_steps;
    temp_convergence_index = num_steps;
    
    % Track convergence status separately for reactants and products
    reactantsConverged = true;
    productsConverged = true;
    temperatureConverged = true;
    
    % Step-by-step convergence detection starting from step 2
    for i = 2:num_steps
        
        % Check NaOH convergence (reactant)
        if naoh_convergence_index == num_steps
            rel_change_NaOH = abs(simResults.C_NaOH(i) - simResults.C_NaOH(i-1)) / ...
                              max(eps, abs(simResults.C_NaOH(i-1))) * 100;
            if rel_change_NaOH < reactant_threshold
                naoh_convergence_index = i;
                convergenceData.NaOH_time = t(i);
                convergenceData.NaOH_conc = simResults.C_NaOH(i);
                convergenceData.NaOH_converged = true;
                convergenceData.NaOH_final = simResults.C_NaOH(end);
            end
        end
        
        % Check EtOAc convergence (reactant)
        if etoac_convergence_index == num_steps
            rel_change_EtOAc = abs(simResults.C_EtOAc(i) - simResults.C_EtOAc(i-1)) / ...
                               max(eps, abs(simResults.C_EtOAc(i-1))) * 100;
            if rel_change_EtOAc < reactant_threshold
                etoac_convergence_index = i;
                convergenceData.EtOAc_time = t(i);
                convergenceData.EtOAc_conc = simResults.C_EtOAc(i);
                convergenceData.EtOAc_converged = true;
                convergenceData.EtOAc_final = simResults.C_EtOAc(end);
            end
        end
        
        % Check NaOAc convergence (product)
        if naoac_convergence_index == num_steps
            rel_change_NaOAc = abs(simResults.C_NaOAc(i) - simResults.C_NaOAc(i-1)) / ...
                               max(eps, abs(simResults.C_NaOAc(i-1))) * 100;
            if rel_change_NaOAc < product_threshold
                naoac_convergence_index = i;
                convergenceData.NaOAc_time = t(i);
                convergenceData.NaOAc_conc = simResults.C_NaOAc(i);
                convergenceData.NaOAc_converged = true;
                convergenceData.NaOAc_final = simResults.C_NaOAc(end);
            end
        end
        
        % Check EtOH convergence (product)
        if etoh_convergence_index == num_steps
            rel_change_EtOH = abs(simResults.C_EtOH(i) - simResults.C_EtOH(i-1)) / ...
                              max(eps, abs(simResults.C_EtOH(i-1))) * 100;
            if rel_change_EtOH < product_threshold
                etoh_convergence_index = i;
                convergenceData.EtOH_time = t(i);
                convergenceData.EtOH_conc = simResults.C_EtOH(i);
                convergenceData.EtOH_converged = true;
                convergenceData.EtOH_final = simResults.C_EtOH(end);
            end
        end
        
        % Check Temperature convergence
        if temp_convergence_index == num_steps
            rel_change_T = abs(simResults.T_reactor(i) - simResults.T_reactor(i-1)) / ...
                          max(eps, abs(simResults.T_reactor(i-1))) * 100;
            if rel_change_T < temp_threshold
                temp_convergence_index = i;
                convergenceData.T_time = t(i);
                convergenceData.T_value = simResults.T_reactor(i);
                convergenceData.T_converged = true;
                convergenceData.T_final = simResults.T_reactor(end);
            end
        end
    end
    
    % Set convergence flags based on whether indices were found
    if naoh_convergence_index == num_steps
        convergenceData.NaOH_converged = false;
        reactantsConverged = false;
    end
    
    if etoac_convergence_index == num_steps
        convergenceData.EtOAc_converged = false;
        reactantsConverged = false;
    end
    
    if naoac_convergence_index == num_steps
        convergenceData.NaOAc_converged = false;
        productsConverged = false;
    end
    
    if etoh_convergence_index == num_steps
        convergenceData.EtOH_converged = false;
        productsConverged = false;
    end
    
    if temp_convergence_index == num_steps
        convergenceData.T_converged = false;
        temperatureConverged = false;
    end
    
    % Overall convergence requires both reactants and products to converge
    isConverged = reactantsConverged && productsConverged;
    
    % Store convergence data and separate status tracking for ultra-precision
    app.UltraConvergenceData = convergenceData;
    app.IsUltraConverged = isConverged;
    
    % Store separate convergence status for ultra-precision simulation
    app.UltraReactantsConverged = reactantsConverged;
    app.UltraProductsConverged = productsConverged;
    app.UltraTemperatureConverged = temperatureConverged;
    
    % Enhanced logging of convergence status
    fprintf('=== ULTRA-PRECISION CONVERGENCE ANALYSIS ===\n');
    fprintf('Reactant Threshold: %.4f%%, Product Threshold: %.4f%%\n', reactant_threshold, product_threshold);
    fprintf('Temperature Threshold: %.4f%%\n', temp_threshold);
    
    if reactantsConverged
        fprintf('Ultra Reactants Converged: NaOH at %.1f min, EtOAc at %.1f min\n', ...
               convergenceData.NaOH_time, convergenceData.EtOAc_time);
    else
        fprintf('Ultra Reactants: Not fully converged\n');
    end
    
    if productsConverged
        fprintf('Ultra Products Converged: NaOAc at %.1f min, EtOH at %.1f min\n', ...
               convergenceData.NaOAc_time, convergenceData.EtOH_time);
    else
        fprintf('Ultra Products: Not fully converged\n');
    end
    
    fprintf('Ultra Overall Converged: %s\n', mat2str(isConverged));
    fprintf('==============================\n');
            
        else % optimized
            app.ConvergenceData = convergenceData;
            app.IsConverged = isConverged;
            
            % Store separate convergence status for optimized simulation
            app.OptimizedReactantsConverged = reactantsConverged;
            app.OptimizedProductsConverged = productsConverged;
            app.OptimizedTemperatureConverged = temperatureConverged;
        end
        
        % Enhanced logging of convergence status with new method
        fprintf('=== RELATIVE CHANGE CONVERGENCE ANALYSIS (%s) ===\n', upper(method));
        fprintf('Reactant Threshold: %.1f%%, Product Threshold: %.1f%%\n', reactant_threshold, product_threshold);
        fprintf('Temperature Threshold: %.1f%%\n', temp_threshold);
        
        if reactantsConverged
            fprintf('Reactants Converged: NaOH at %.1f min, EtOAc at %.1f min\n', ...
                   convergenceData.NaOH_time, convergenceData.EtOAc_time);
        else
            fprintf('Reactants: Not fully converged\n');
        end
        
        if productsConverged
            fprintf('Products Converged: NaOAc at %.1f min, EtOH at %.1f min\n', ...
                   convergenceData.NaOAc_time, convergenceData.EtOH_time);
        else
            fprintf('Products: Not fully converged\n');
        end
        
        fprintf('Overall Converged: %s\n', mat2str(isConverged));
        fprintf('==============================\n');
        
    catch ME
        fprintf('Convergence detection error: %s\n', ME.message);
        if strcmp(method, 'manual')
            app.ManualConvergenceData = struct();
            app.IsManualConverged = false;
            app.ManualReactantsConverged = false;
            app.ManualProductsConverged = false;
            app.ManualTemperatureConverged = false;
        elseif strcmp(method, 'quicksim')
            app.QuickSimConvergenceData = struct();
            app.IsQuickSimConverged = false;
            app.QuickSimReactantsConverged = false;
            app.QuickSimProductsConverged = false;
            app.QuickSimTemperatureConverged = false;
        else
            app.ConvergenceData = struct();
            app.IsConverged = false;
            app.OptimizedReactantsConverged = false;
            app.OptimizedProductsConverged = false;
            app.OptimizedTemperatureConverged = false;
        end
    end
    if strcmp(method, 'manual')
    fprintf('Manual convergence data stored:\n');
    if isfield(app.ManualConvergenceData, 'NaOH_converged')
        fprintf('  NaOH converged: %s\n', mat2str(app.ManualConvergenceData.NaOH_converged));
        if app.ManualConvergenceData.NaOH_converged
            fprintf('  NaOH time: %.2f min\n', app.ManualConvergenceData.NaOH_time);
            fprintf('  NaOH conc: %.6f mol/L\n', app.ManualConvergenceData.NaOH_conc);
        end
    end
elseif strcmp(method, 'optimized')
    fprintf('Optimized convergence data stored:\n');
    if isfield(app.ConvergenceData, 'NaOH_converged')
        fprintf('  NaOH converged: %s\n', mat2str(app.ConvergenceData.NaOH_converged));
        if app.ConvergenceData.NaOH_converged
            fprintf('  NaOH time: %.2f min\n', app.ConvergenceData.NaOH_time);
            fprintf('  NaOH conc: %.6f mol/L\n', app.ConvergenceData.NaOH_conc);
        end
    end
    end
end
   % 1. UPDATE addConvergenceMarkers function - REPLACE EXISTING FUNCTION
function addConvergenceMarkers(app, axes, simResults, convergenceData, color, method)
    try
        if isempty(convergenceData)
            return;
        end
        
        % Check if NaOH convergence data exists and is valid
        if isfield(convergenceData, 'NaOH_time') && isfield(convergenceData, 'NaOH_conc') && ...
           isfield(convergenceData, 'NaOH_converged') && convergenceData.NaOH_converged && ...
           ~isnan(convergenceData.NaOH_time) && ~isnan(convergenceData.NaOH_conc)
            
            conv_time = convergenceData.NaOH_time;
            conv_conc = convergenceData.NaOH_conc;
            final_conc = convergenceData.NaOH_final;
            
            % Define marker styles for different methods
            if strcmp(method, 'Manual')
                linestyle = '-.';  % Dash-dot for manual
                alpha = 0.6;
                marker = 's';      % Square marker
                size = 80;
                lineWidth = 2;
            elseif strcmp(method, 'QuickSim')
                linestyle = ':';   % Dotted for quick sim
                alpha = 0.7;
                marker = 'd';      % Diamond marker
                size = 90;
                lineWidth = 2;
            elseif strcmp(method, 'Ultra-Precision')
                linestyle = '-';   % Solid for ultra-precision
                alpha = 0.95;      % Higher alpha for prominence
                marker = 'p';      % Pentagon marker
                size = 160;        % Larger size for ultra-precision
                lineWidth = 3.5;   % Much thicker lines
            else % Optimized
                linestyle = '--';  % Dashed for optimized
                alpha = 0.8;
                marker = 'o';      % Circle marker
                size = 100;
                lineWidth = 2;
            end
            
            % Add vertical line for convergence time
            xline(axes, conv_time, linestyle, 'Color', color, 'Alpha', alpha, ...
                  'LineWidth', lineWidth, 'DisplayName', sprintf('%s NaOH Convergence Time', method));
            
            % Add horizontal line for steady-state concentration
            yline(axes, final_conc, ':', 'Color', color, 'Alpha', alpha, ...
                  'LineWidth', lineWidth-0.5, 'DisplayName', sprintf('%s NaOH Steady State', method));
            
            % Add convergence point marker
            scatter(axes, conv_time, conv_conc, size, ...
                   'filled', 'Marker', marker, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 2.5, 'DisplayName', sprintf('%s NaOH Convergence', method));
            
            % Add text annotation
            xlims = xlim(axes);
            ylims = ylim(axes);
            text_x = conv_time + 0.02 * (xlims(2) - xlims(1));
            text_y = conv_conc + 0.02 * (ylims(2) - ylims(1));
            
            % Offset annotations to avoid overlap
            if strcmp(method, 'Manual')
                text_y = text_y + 0.05 * (ylims(2) - ylims(1));
                fontSize = 8;
            elseif strcmp(method, 'QuickSim')
                text_x = text_x - 0.08 * (xlims(2) - xlims(1));
                fontSize = 8;
            elseif strcmp(method, 'Ultra-Precision')
                text_x = text_x + 0.05 * (xlims(2) - xlims(1));
                text_y = text_y - 0.03 * (ylims(2) - ylims(1));
                fontSize = 10;  % Larger font for ultra-precision
            else
                fontSize = 8;
            end
            
            % Create annotation text with emoji for ultra-precision
            if strcmp(method, 'Ultra-Precision')
                annotation_text = sprintf('🏆 %s Conv\nt=%.1fmin\nC=%.6f mol/L', method, conv_time, conv_conc);
            else
                annotation_text = sprintf('%s Conv\nt=%.1fmin\nC=%.4f mol/L', method, conv_time, conv_conc);
            end
            
            text(axes, text_x, text_y, annotation_text, ...
                 'FontSize', fontSize, 'Color', color, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', 'EdgeColor', color, 'Margin', 3);
        end
        
        % Also check for EtOAc convergence
        if isfield(convergenceData, 'EtOAc_time') && isfield(convergenceData, 'EtOAc_conc') && ...
           isfield(convergenceData, 'EtOAc_converged') && convergenceData.EtOAc_converged && ...
           ~isnan(convergenceData.EtOAc_time) && ~isnan(convergenceData.EtOAc_conc)
            
            conv_time = convergenceData.EtOAc_time;
            conv_conc = convergenceData.EtOAc_conc;
            
            % Add markers for EtOAc convergence (slightly different style)
            if strcmp(method, 'Ultra-Precision')
                xline(axes, conv_time, ':', 'Color', color*0.8, 'Alpha', 0.9, ...
                      'LineWidth', 2.5, 'DisplayName', sprintf('%s EtOAc Convergence Time', method));
                
                scatter(axes, conv_time, conv_conc, 120, ...
                       'filled', 'Marker', '^', 'MarkerFaceColor', color*0.8, 'MarkerEdgeColor', 'black', ...
                       'LineWidth', 2, 'DisplayName', sprintf('%s EtOAc Convergence', method));
            else
                xline(axes, conv_time, ':', 'Color', color*0.8, 'Alpha', alpha*0.7, ...
                      'LineWidth', 1.5, 'DisplayName', sprintf('%s EtOAc Convergence Time', method));
                
                scatter(axes, conv_time, conv_conc, size*0.7, ...
                       'filled', 'Marker', '^', 'MarkerFaceColor', color*0.8, 'MarkerEdgeColor', 'black', ...
                       'LineWidth', 1, 'DisplayName', sprintf('%s EtOAc Convergence', method));
            end
        end
        
    catch ME
        fprintf('Error adding convergence markers: %s\n', ME.message);
    end
end
     


function addProductConvergenceMarkers(app, axes, simResults, convergenceData, method)
    try
        if isempty(convergenceData)
            return;
        end
        
        % Define colors and styles for different methods
        if strcmp(method, 'Manual')
            naoacColor = [255 204 0]/255 * 0.7; % More muted yellow for manual
            etohColor = [88 86 214]/255 * 0.7;  % More muted indigo for manual
            linestyle = '-.';  % Dash-dot for manual
            alpha = 0.6;
            marker = 's';      % Square marker
            size = 80;
            lineWidth = 2;
            fontSize = 7;
        elseif strcmp(method, 'QuickSim')
            naoacColor = [255 165 0]/255;       % Orange for quick sim
            etohColor = [30 144 255]/255;       % Dodger blue for quick sim
            linestyle = ':';   % Dotted for quick sim
            alpha = 0.7;
            marker = 'd';      % Diamond marker
            size = 90;
            lineWidth = 2;
            fontSize = 8;
        elseif strcmp(method, 'Ultra-Precision')
            naoacColor = [255 20 147]/255;       % Deep pink for ultra-precision
            etohColor = [139 0 139]/255;         % Dark magenta for ultra-precision
            linestyle = '-';   % Solid for ultra-precision
            alpha = 0.95;      % Higher alpha for prominence
            marker = 'p';      % Pentagon marker
            size = 180;        % Much larger size for ultra-precision
            lineWidth = 3.5;   % Much thicker lines
            fontSize = 10;     % Larger font
        else
            naoacColor = [255 204 0]/255;       % Bright yellow for optimized
            etohColor = [88 86 214]/255;        % Bright indigo for optimized
            linestyle = '--';  % Dashed for optimized
            alpha = 0.8;
            marker = 'o';      % Circle marker
            size = 100;
            lineWidth = 2;
            fontSize = 8;
        end
        
        % Add convergence markers for NaOAc (product)
        if isfield(convergenceData, 'NaOAc_time') && isfield(convergenceData, 'NaOAc_conc') && convergenceData.NaOAc_converged
            conv_time = convergenceData.NaOAc_time;
            conv_conc = convergenceData.NaOAc_conc;
            final_conc = convergenceData.NaOAc_final;
            
            % Vertical line for convergence time
            xline(axes, conv_time, linestyle, 'Color', naoacColor, 'Alpha', alpha, ...
                  'LineWidth', lineWidth, 'DisplayName', sprintf('%s NaOAc Convergence', method));
            
            % Horizontal line for steady-state
            yline(axes, final_conc, ':', 'Color', naoacColor, 'Alpha', alpha, ...
                  'LineWidth', lineWidth, 'DisplayName', sprintf('%s NaOAc SS', method));
            
            % Point marker
            scatter(axes, conv_time, conv_conc, size, ...
                   'filled', 'Marker', marker, 'MarkerFaceColor', naoacColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 2.5, 'DisplayName', sprintf('%s NaOAc Point', method));
            
            % Text annotation
            xlims = xlim(axes);
            ylims = ylim(axes);
            text_x = conv_time + 0.02 * (xlims(2) - xlims(1));
            text_y = final_conc + 0.02 * (ylims(2) - ylims(1));
            
            % Offset annotations to avoid overlap
            if strcmp(method, 'Manual')
                text_x = text_x - 0.05 * (xlims(2) - xlims(1));
                text_y = text_y + 0.03 * (ylims(2) - ylims(1));
            elseif strcmp(method, 'QuickSim')
                text_x = text_x - 0.08 * (xlims(2) - xlims(1));
                text_y = text_y - 0.02 * (ylims(2) - ylims(1));
            elseif strcmp(method, 'Ultra-Precision')
                text_x = text_x + 0.03 * (xlims(2) - xlims(1));
                text_y = text_y + 0.05 * (ylims(2) - ylims(1));
            end
            
            % Enhanced annotation for ultra-precision
            if strcmp(method, 'Ultra-Precision')
                annotation_text = sprintf('🏆 NaOAc %s\nt=%.1f\nC=%.6f', method, conv_time, final_conc);
            else
                annotation_text = sprintf('NaOAc %s\nt=%.1f\nC=%.4f', method, conv_time, final_conc);
            end
            
            text(axes, text_x, text_y, annotation_text, ...
                 'FontSize', fontSize, 'Color', naoacColor, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', 'EdgeColor', naoacColor, 'Margin', 2);
        end
        
        % Add convergence markers for EtOH (product)
        if isfield(convergenceData, 'EtOH_time') && isfield(convergenceData, 'EtOH_conc') && convergenceData.EtOH_converged
            conv_time = convergenceData.EtOH_time;
            conv_conc = convergenceData.EtOH_conc;
            final_conc = convergenceData.EtOH_final;
            
            % Similar markers for EtOH
            xline(axes, conv_time, linestyle, 'Color', etohColor, 'Alpha', alpha, ...
                  'LineWidth', lineWidth, 'DisplayName', sprintf('%s EtOH Convergence', method));
            
            yline(axes, final_conc, ':', 'Color', etohColor, 'Alpha', alpha, ...
                  'LineWidth', lineWidth, 'DisplayName', sprintf('%s EtOH SS', method));
            
            scatter(axes, conv_time, conv_conc, size, ...
                   'filled', 'Marker', marker, 'MarkerFaceColor', etohColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 2.5, 'DisplayName', sprintf('%s EtOH Point', method));
            
            % Text annotation for EtOH
            xlims = xlim(axes);
            ylims = ylim(axes);
            text_x = conv_time + 0.02 * (xlims(2) - xlims(1));
            text_y = final_conc - 0.02 * (ylims(2) - ylims(1));
            
            % Offset annotations to avoid overlap
            if strcmp(method, 'Manual')
                text_x = text_x - 0.05 * (xlims(2) - xlims(1));
                text_y = text_y - 0.03 * (ylims(2) - ylims(1));
            elseif strcmp(method, 'QuickSim')
                text_x = text_x - 0.08 * (xlims(2) - xlims(1));
                text_y = text_y + 0.02 * (ylims(2) - ylims(1));
            elseif strcmp(method, 'Ultra-Precision')
                text_x = text_x + 0.03 * (xlims(2) - xlims(1));
                text_y = text_y - 0.05 * (ylims(2) - ylims(1));
            end
            
            % Enhanced annotation for ultra-precision
            if strcmp(method, 'Ultra-Precision')
                annotation_text = sprintf('🏆 EtOH %s\nt=%.1f\nC=%.6f', method, conv_time, final_conc);
            else
                annotation_text = sprintf('EtOH %s\nt=%.1f\nC=%.4f', method, conv_time, final_conc);
            end
            
            text(axes, text_x, text_y, annotation_text, ...
                 'FontSize', fontSize, 'Color', etohColor, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', 'EdgeColor', etohColor, 'Margin', 2);
        end
        
    catch ME
        fprintf('Error adding product convergence markers: %s\n', ME.message);
    end
end

        
        % NEW METHOD: Quick Simulation
        function runQuickSim(app, ~)
            % Clear cooling history for new simulation
app.CoolingTimeHistory = [];
app.CoolingFlowrateHistory = [];
app.CoolingCapacityHistory = [];
app.CoolingAdequacyHistory = [];
app.TemperatureDeviationHistory = [];
            try
                % Read inputs
                Cin = app.CinQuickField.Value;
                V0 = app.FlowrateQuickField.Value;
                T0 = str2double(app.TempQuickDropDown.Value);
                tEnd = app.SimTimeQuickField.Value;
                k = app.RateConstQuickField.Value;
                
                % Validate inputs
                if k <= 0 || ~isfinite(k)
                    app.updateStatus('Invalid rate constant!', 'error');
                    return;
                end
                
                app.updateStatus('Running Quick Simulation with convergence analysis...', 'info');
                
                % Set flowrates
                V0_EtOAc = V0/2;
                V1_NaOH = V0/2;
                V = app.VolumeField.Value;
                
                % Initial conditions
                y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
                
                % Simulate
                tspan = linspace(0, tEnd, 500);
                options = odeset('RelTol',1e-8,'AbsTol',1e-10);
                [t,y] = ode45(@(t,y) app.odefun_individual(t,y,k,V0_EtOAc,V1_NaOH,V,Cin,T0), tspan, y0, options);
                
                % Store results for convergence analysis
                app.QuickSimResults = struct();
                app.QuickSimResults.t = t;
                app.QuickSimResults.C_NaOH = y(:,1);
                app.QuickSimResults.C_EtOAc = y(:,2);
                app.QuickSimResults.C_NaOAc = y(:,3);
                app.QuickSimResults.C_EtOH = y(:,4);
                app.QuickSimResults.T_reactor = y(:,5);
                app.QuickSimResults.T_jacket = y(:,6);
                
                % Detect convergence for quick simulation
                app.detectConvergence('quicksim');
                
                % Plot concentrations with convergence markers
                cla(app.UIAxesQuickConc);
                hold(app.UIAxesQuickConc, 'on');
                
                % Plot concentration profiles
                plot(app.UIAxesQuickConc, t, y(:,1), '-', 'LineWidth', 2, 'Color', [255 59 48]/255, 'DisplayName', 'NaOH');
                plot(app.UIAxesQuickConc, t, y(:,2), '--', 'LineWidth', 2, 'Color', [255 149 0]/255, 'DisplayName', 'EtOAc');
                plot(app.UIAxesQuickConc, t, y(:,3), '-.', 'LineWidth', 2, 'Color', [255 204 0]/255, 'DisplayName', 'NaOAc');
                plot(app.UIAxesQuickConc, t, y(:,4), ':', 'LineWidth', 2, 'Color', [88 86 214]/255, 'DisplayName', 'EtOH');
                
                % Add convergence markers
                if ~isempty(app.QuickSimConvergenceData)
                    app.addConvergenceMarkers(app.UIAxesQuickConc, app.QuickSimResults, app.QuickSimConvergenceData, [34 139 34]/255, 'QuickSim');
                end
                
                title(app.UIAxesQuickConc, 'Concentration Profiles with Quick Sim Convergence Analysis', 'FontWeight', 'bold');
                xlabel(app.UIAxesQuickConc, 'Time (min)');
                ylabel(app.UIAxesQuickConc, 'Concentration (mol/L)');
                legend(app.UIAxesQuickConc, 'Location', 'best');
                grid(app.UIAxesQuickConc, 'off');
                hold(app.UIAxesQuickConc, 'off');
                
                % Plot temperatures
                cla(app.UIAxesQuickTemp);
                hold(app.UIAxesQuickTemp, 'on');
                plot(app.UIAxesQuickTemp, t, y(:,5), '-', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Reactor T');
                plot(app.UIAxesQuickTemp, t, y(:,6), '--', 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Jacket T');
                title(app.UIAxesQuickTemp, 'Temperature Profiles', 'FontWeight', 'bold');
                xlabel(app.UIAxesQuickTemp, 'Time (min)');
                ylabel(app.UIAxesQuickTemp, 'Temperature (°C)');
                legend(app.UIAxesQuickTemp, 'Location', 'best');
                grid(app.UIAxesQuickTemp, 'off');
                hold(app.UIAxesQuickTemp, 'off');
                
                % Plot products with convergence markers
                cla(app.UIAxesQuickProducts);
                hold(app.UIAxesQuickProducts, 'on');
                plot(app.UIAxesQuickProducts, t, y(:,3), '-', 'LineWidth', 2, 'Color', [255 204 0]/255, 'DisplayName', 'NaOAc');
                plot(app.UIAxesQuickProducts, t, y(:,4), '--', 'LineWidth', 2, 'Color', [88 86 214]/255, 'DisplayName', 'EtOH');
                
                % Add product convergence markers
                if ~isempty(app.QuickSimConvergenceData)
                    app.addProductConvergenceMarkers(app.UIAxesQuickProducts, app.QuickSimResults, app.QuickSimConvergenceData, 'QuickSim');
                end
                
                title(app.UIAxesQuickProducts, 'Product Formation with Quick Sim Convergence Analysis', 'FontWeight', 'bold');
                xlabel(app.UIAxesQuickProducts, 'Time (min)');
                ylabel(app.UIAxesQuickProducts, 'Concentration (mol/L)');
                legend(app.UIAxesQuickProducts, 'Location', 'best');
                grid(app.UIAxesQuickProducts, 'off');
                hold(app.UIAxesQuickProducts, 'off');
                
                % Calculate manual-method errors if experimental data exists
                if ~isempty(app.Data) && height(app.Data) > 0
                    C_exp = app.Data.Concentration;
                    t_exp = app.Data.Time;
                    C_sim = interp1(t, y(:,1), t_exp, 'linear', 'extrap');
                    errors = abs(C_exp - C_sim) ./ C_exp * 100;
                    
                    maxE = max(errors);
                    minE = min(errors);
                    meanE = mean(errors);
                    
                    app.MaxErrorLabel.Text = sprintf('Max Error: %.2f%%', maxE);
                    app.MinErrorLabel.Text = sprintf('Min Error: %.2f%%', minE);
                    app.MeanErrorLabel.Text = sprintf('Mean Error: %.2f%%', meanE);
                else
                    app.MaxErrorLabel.Text = 'Max Error: N/A';
                    app.MinErrorLabel.Text = 'Min Error: N/A';
                    app.MeanErrorLabel.Text = 'Mean Error: N/A';
                end
                
                % Update info panel
                app.updateQuickSimInfo();
                
                % Update conversion plots
app.plotConversions();
                
                % Update status with convergence information
                if app.IsQuickSimConverged
                    app.updateStatus(sprintf('Quick simulation completed! Converged at medium thresholds (React: %.1e, Prod: %.1e)', ...
                        app.QuickSimReactantThreshold, app.QuickSimProductThreshold), 'success');
                else
                    app.updateStatus('Quick simulation completed! Not fully converged with medium thresholds', 'warning');
                end
                
            catch ME
                app.updateStatus(['Quick simulation failed: ' ME.message], 'error');
            end
        end
        
        % NEW METHOD: Update Quick Sim Info Panel
        function updateQuickSimInfo(app)
            try
                infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
                infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">Quick Sim Convergence Analysis</h3>'];
                
                % Quick sim convergence thresholds
                infoHTML = [infoHTML '<strong>Quick Sim Thresholds (Medium Criteria):</strong><br>'];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">Reactants: %.1e</p>', app.QuickSimReactantThreshold)];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">Products: %.1e</p>', app.QuickSimProductThreshold)];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">Temperature: %.3f°C/min</p>', app.QuickSimTemperatureThreshold)];
                
                % Quick sim convergence status
                infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                infoHTML = [infoHTML '<strong>Quick Sim Convergence Status:</strong><br>'];
                
                if app.QuickSimReactantsConverged
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Reactants: Converged</p>'];
                else
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Reactants: Not Converged</p>'];
                end
                
                if app.QuickSimProductsConverged
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Products: Converged</p>'];
                else
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Products: Not Converged</p>'];
                end
                
                if app.QuickSimTemperatureConverged
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Temperature: Converged</p>'];
                else
                    infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Temperature: Not Converged</p>'];
                end
                
                % Overall status with color coding
                if app.IsQuickSimConverged
                    infoHTML = [infoHTML '<p style="margin: 5px 0; color: #34C759; font-weight: bold;">✅ Overall: CONVERGED</p>'];
                else
                    infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500; font-weight: bold;">⚠ Overall: NOT CONVERGED</p>'];
                end
                
                % Comparison with other methods
                if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
                    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                    infoHTML = [infoHTML '<strong>Convergence Comparison:</strong><br>'];
                    
                    if ~isempty(app.ManualSimResults)
                        if app.IsManualConverged
                            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">Manual (relaxed): ✓ Converged</p>'];
                        else
                            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">Manual (relaxed): ⚠ Not Converged</p>'];
                        end
                    end
                    
                    if app.IsQuickSimConverged
                        infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">Quick Sim (medium): ✓ Converged</p>'];
                    else
                        infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">Quick Sim (medium): ⚠ Not Converged</p>'];
                    end
                    
                    if ~isempty(app.SimResults)
                        if app.IsConverged
                            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">Optimized (strict): ✓ Converged</p>'];
                        else
                            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">Optimized (strict): ⚠ Not Converged</p>'];
                        end
                    end
                end
                
                infoHTML = [infoHTML '</div>'];
                app.QuickSimInfo.HTMLSource = infoHTML;
                
            catch
                % Silent fail
            end
        end
        
        % Core functionality methods
        function runSimulation(app)
            % Clear cooling history for new simulation
app.CoolingTimeHistory = [];
app.CoolingFlowrateHistory = [];
app.CoolingCapacityHistory = [];
app.CoolingAdequacyHistory = [];
app.TemperatureDeviationHistory = [];
            try
                app.enableControls('off');
                app.showLoadingDialog('Running Optimization', 'Initializing simulation...');
                
                % Validate inputs
                if isempty(app.Data) || isempty(app.SelectedTemp)
                    app.closeLoadingDialog();
                    app.enableControls('on');
                    app.updateStatus('Please upload data and select temperature', 'error');
                    uialert(app.UIFigure, 'Please upload data and select a temperature.', 'Input Required', 'Icon', 'error');
                    return;
                end
                
                        
        % Initialize optimization results structure  % <-- ADD THIS SECTION
        app.OptimizationResults = struct();
        app.OptimizationResults.InitialFlowrate = app.FlowrateField.Value;
        app.OptimizationResults.OptimalFlowrate = NaN;
        app.OptimizationResults.OptimalK = NaN;
        app.OptimizationResults.MaxError = NaN;
        app.OptimizationResults.RSquared = NaN;
                
                % Get parameters
                t_exp = app.Data.Time;
                C_exp = app.Data.Concentration;
                T0 = app.SelectedTemp;
                To = T0;
                V = app.VolumeField.Value;
                initial_V0 = app.FlowrateField.Value;
                Cin = app.CinField.Value;
                extTime = app.SimTimeField.Value;
                
                app.updateLoadingDialog(0.1, 'Updating parameters...');
                
                % Initial conditions [CA, CB, CC, CD, T_reactor, Tj]
                y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
                
                % IMPROVED OPTIMIZATION: Optimize both flowrate and rate constant
                app.updateLoadingDialog(0.2, 'Optimizing flowrate and rate constant for 1% error target...');
                
                % Define bounds for optimization
                V0_min = 0.01;   % [L/min]
                V0_max = 1.0;    % [L/min] 
                k_min = 0.001;   % [L/mol·min]
                k_max = 10;      % [L/mol·min]
                
                % Enhanced objective function that returns maximum relative error
                objectiveFunc = @(params) app.calculateMaxErrorForOptimization(params, t_exp, C_exp, y0, V, Cin, To);
                
                % Initial guess: [V0_total, k]
                x0 = [initial_V0, 0.5];
                lb = [V0_min, k_min];
                ub = [V0_max, k_max];
                
                % Use constrained optimization with enhanced settings
                options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 300, ...
                                      'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                                      'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', 1500);
                
                [x_optimal, max_error] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
                
                V0_optimal = x_optimal(1);
                k_fit = x_optimal(2);
                
                if k_fit <= 0 || ~isfinite(k_fit)
                    error('Invalid rate constant fitted: k = %.4f', k_fit);
                end
                
                % Update parameters with optimal values
                app.FlowrateField.Value = V0_optimal;
                V0_EtOAc = V0_optimal / 2;
                V1_NaOH = V0_optimal / 2;
                
                app.k_fit = k_fit;
                app.updateLoadingDialog(0.4, sprintf('Optimal parameters found. Max error: %.3f%%', max_error));
                % UPDATE OPTIMIZATION RESULTS STRUCTURE
app.OptimizationResults.OptimalFlowrate = V0_optimal;
app.OptimizationResults.OptimalK = k_fit;
app.OptimizationResults.MaxError = max_error;
                
                % Update properties
                app.V = V;
                app.V0 = V0_optimal;
                app.V0_EtOAc = V0_EtOAc;
                app.V1_NaOH = V1_NaOH;
                app.Cin = Cin;
                app.SimTime = extTime;
                
                app.updateLoadingDialog(0.5, 'Running final simulation...');
                
                % Check if this is an extension and adjust initial conditions
if isfield(app, 'ExtensionMode') && app.ExtensionMode
    % Modify the experimental data for optimization to start from extension point
    originalIdx = find(app.Data.Time <= app.OriginalSimTime, 1, 'last');
    if ~isempty(originalIdx)
        % Use data from extension point forward
        t_exp = app.Data.Time(originalIdx:end) - app.OriginalSimTime; % Reset time to start from 0
        C_exp = app.Data.Concentration(originalIdx:end);
        
        % Update initial conditions to start from extension point
        y0 = [app.OriginalEndConcentration, app.OriginalEndConcentration, 0, 0, T0, app.Tj0];
        
        app.updateStatus('Optimizing for extended simulation from continuation point', 'info');
    end
    app.ExtensionMode = false; % Reset flag
end
                
                % Simulate with optimal parameters
                t_span = linspace(0, extTime, 2000);
                odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-11); % Tighter tolerances
                [t_sim, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,k_fit,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                      t_span, y0, odeOptions);
                
                app.updateLoadingDialog(0.6, 'Storing results...');
                
                % Store results
                app.SimResults = struct();
                app.SimResults.t = t_sim;
                app.SimResults.C_NaOH = y_sim(:,1);
                app.SimResults.C_EtOAc = y_sim(:,2);
                app.SimResults.C_NaOAc = y_sim(:,3);
                app.SimResults.C_EtOH = y_sim(:,4);
                app.SimResults.T_reactor = y_sim(:,5);
                app.SimResults.T_jacket = y_sim(:,6);
                % DEBUG: Check simulation results storage
fprintf('=== SIMULATION RESULTS STORAGE DEBUG ===\n');
fprintf('Time points stored: %d\n', length(t_sim));
fprintf('y_sim dimensions: %dx%d\n', size(y_sim, 1), size(y_sim, 2));
fprintf('NaOH final: %.8f, initial: %.8f\n', y_sim(end,1), y_sim(1,1));
fprintf('EtOAc final: %.8f, initial: %.8f\n', y_sim(end,2), y_sim(1,2));
fprintf('NaOAc final: %.8f, initial: %.8f\n', y_sim(end,3), y_sim(1,3));
fprintf('EtOH final: %.8f, initial: %.8f\n', y_sim(end,4), y_sim(1,4));
fprintf('=======================================\n');
                
                app.updateLoadingDialog(0.7, 'Analyzing convergence...');
                
                % Detect convergence using optimized thresholds
                app.detectConvergence('optimized');
                
                % Check convergence and suggest extension if needed
                app.checkAndSuggestTimeExtension('optimized');
                
                % Create comparison table for optimized results
                C_sim_interp = interp1(t_sim, y_sim(:,1), t_exp, 'linear', 'extrap');
                app.createComparisonTable(t_exp, C_exp, t_sim, y_sim(:,1), 'optimized');
                
                app.updateLoadingDialog(0.8, 'Calculating metrics...');
                
                % Calculate R²
                app.RSquared = app.calculateRSquared(C_exp, C_sim_interp);
                % UPDATE R-SQUARED IN OPTIMIZATION RESULTS
app.OptimizationResults.RSquared = app.RSquared;
                
                app.updateLoadingDialog(0.9, 'Updating visualizations...');
                
                % Update all plots
                app.plotConcentrations();
                app.plotTemperatures();
                app.plotProducts();
                app.plotConversions();
                app.updateAnalysisDisplay();
                
                app.updateLoadingDialog(1.0, 'Complete!');              
                
% ADD THIS LINE TO CLOSE THE DIALOG
app.closeLoadingDialog();
                
                % Calculate final error metrics
              % Calculate final error metrics
error_percent = abs(C_exp - C_sim_interp) ./ C_exp * 100;
mean_error = mean(error_percent);
max_error_final = max(error_percent);
% FINAL UPDATE OF OPTIMIZATION RESULTS WITH CALCULATED METRICS
app.OptimizationResults.MaxError = max_error_final;  % Update with final calculated error
app.OptimizationResults.MeanError = mean_error;     % Add mean error
app.OptimizationResults.Method = 'Constrained Optimization (fmincon)';

% Calculate temperature performance
if app.IsothermalMode
    temp_deviation = max(abs(app.SimResults.T_reactor - app.SelectedTemp));
    if app.cooling_adequacy_factor > 1.2
        cooling_status = '✅ ADEQUATE';
    elseif app.cooling_adequacy_factor > 0.8
        cooling_status = '⚠️ MARGINAL';
    else
        cooling_status = '❌ INADEQUATE';
    end
    
    app.updateStatus(sprintf('🌡️ Isothermal simulation complete! Max T deviation: %.3f°C | Cooling: %s', ...
        temp_deviation, cooling_status), 'success');
    
    % Show cooling system performance dialog
    app.showCoolingSystemReport(temp_deviation, cooling_status);
else
    % Original success message
    app.updateStatus(sprintf('🎯 Optimization complete: Mean: %.3f%%, Max: %.3f%%, R²: %.4f', ...
        mean_error, max_error_final, app.RSquared), 'success');
end
                
                app.enableControls('on');
                
                % Switch to concentration tab
                app.TabGroup.SelectedTab = app.ConcentrationTab;
                
            catch ME
                app.closeLoadingDialog();
                app.enableControls('on');
                app.updateProgress(0);
                app.updateStatus(['Simulation failed: ' ME.message], 'error');
                uialert(app.UIFigure, ['Simulation failed: ', ME.message], 'Error', 'Icon', 'error');
            end
            % Check if this is an extension and adjust initial conditions
if isfield(app, 'ExtensionMode') && app.ExtensionMode
    % Modify the experimental data for optimization to start from extension point
    originalIdx = find(app.Data.Time <= app.OriginalSimTime, 1, 'last');
    if ~isempty(originalIdx)
        % Use data from extension point forward
        t_exp = app.Data.Time(originalIdx:end) - app.OriginalSimTime; % Reset time to start from 0
        C_exp = app.Data.Concentration(originalIdx:end);
        
        % Update initial conditions to start from extension point
        y0 = [app.OriginalEndConcentration, app.OriginalEndConcentration, 0, 0, T0, app.Tj0];
        
        app.updateStatus('Optimizing for extended simulation from continuation point', 'info');
    end
    app.ExtensionMode = false; % Reset flag
end
            app.updateProfileTables();
        end
        
    function calculateManualRateConstant(app, startTime, startConc)
    try
        if nargin < 3
            % Original calculation from beginning
            if isempty(app.Data) || height(app.Data) == 0
                error('No experimental data available');
            end
            CA_out = app.Data.Concentration(end);
            CA0 = app.CinField.Value;
            tau = app.VolumeField.Value / app.FlowrateField.Value;
        else
            % Extended calculation from specific time point
            CA_out = app.Data.Concentration(end);
            CA0 = startConc; % Use concentration at extension point as new initial
            extendedTime = app.SimTimeField.Value - startTime;
            tau = extendedTime; % Use extended time period
        end
        
        % Manual calculation for second-order reaction in CSTR
        if CA_out <= 0 || CA_out >= CA0
            error('Invalid steady-state concentration');
        end
        
        % k [L/mol·min] = (CA0 - CA_out) / (tau * CA_out^2)
        k_manual = (CA0 - CA_out) / (tau * CA_out^2);
        
        if k_manual <= 0 || ~isfinite(k_manual)
            error('Invalid manual rate constant: k = %.4f', k_manual);
        end
        
        app.k_manual = k_manual;
        
        if nargin >= 3
            app.updateStatus(sprintf('Extended Manual k = %.4f L/mol·min (from t=%.1f)', k_manual, startTime), 'success');
        else
            app.updateStatus(sprintf('Manual k = %.4f L/mol·min', k_manual), 'success');
        end
        
    catch ME
        app.updateStatus(['Manual calculation failed: ' ME.message], 'error');
        app.k_manual = NaN;
    end
end
        
        function runManualSimulation(app)
            
            % Clear cooling history for new simulation
app.CoolingTimeHistory = [];
app.CoolingFlowrateHistory = [];
app.CoolingCapacityHistory = [];
app.CoolingAdequacyHistory = [];
app.TemperatureDeviationHistory = [];
            try
                app.showLoadingDialog('Manual Simulation', 'Calculating manual rate constant...');
                
                if isnan(app.k_manual) || app.k_manual <= 0
                    app.updateLoadingDialog(0.2, 'Calculating rate constant...');
                    app.calculateManualRateConstant();
                    if isnan(app.k_manual)
                        app.closeLoadingDialog();
                        return;
                    end
                end
                
                % Simulation parameters
                T0 = app.SelectedTemp;
                To = T0;
                V = app.VolumeField.Value;
                V0 = app.FlowrateField.Value;
                V0_EtOAc = V0 / 2;
                V1_NaOH = V0 / 2;
                Cin = app.CinField.Value;
                extTime = app.SimTimeField.Value;
                
                % Initial conditions
                y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
                
                app.updateLoadingDialog(0.4, 'Running manual simulation...');
                
                % Simulate
                t_span = linspace(0, extTime, 2000);
                options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
                [t_sim, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,app.k_manual,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                      t_span, y0, options);
                
                app.updateLoadingDialog(0.6, 'Storing manual results...');
                
                % Store results
                app.ManualSimResults = struct();
                app.ManualSimResults.t = t_sim;
                app.ManualSimResults.C_NaOH = y_sim(:,1);
                app.ManualSimResults.C_EtOAc = y_sim(:,2);
                app.ManualSimResults.C_NaOAc = y_sim(:,3);
                app.ManualSimResults.C_EtOH = y_sim(:,4);
                app.ManualSimResults.T_reactor = y_sim(:,5);
                app.ManualSimResults.T_jacket = y_sim(:,6);
                
                app.updateLoadingDialog(0.7, 'Analyzing manual convergence...');
                
                % Detect convergence for manual simulation using manual thresholds
                app.detectConvergence('manual');
                
                % Check convergence and suggest extension if needed
                app.checkAndSuggestTimeExtension('manual');
                
                app.updateLoadingDialog(0.8, 'Calculating manual metrics...');
                
                % Calculate R² for manual method and create comparison table
                if ~isempty(app.Data)
                    t_exp = app.Data.Time;
                    C_exp = app.Data.Concentration;
                    C_manual_interp = interp1(t_sim, y_sim(:,1), t_exp, 'linear', 'extrap');
                    app.RSquaredManual = app.calculateRSquared(C_exp, C_manual_interp);
                    
                    % Create comparison table for manual results
                    app.createComparisonTable(t_exp, C_exp, t_sim, y_sim(:,1), 'manual');
                end
                
                app.updateLoadingDialog(1.0, 'Manual simulation complete!');
                app.closeLoadingDialog();
                
                app.updateStatus(sprintf('Manual simulation complete! k=%.4f, R²=%.4f', ...
                                       app.k_manual, app.RSquaredManual), 'success');
                
            catch ME
                app.closeLoadingDialog();
                app.updateStatus(['Manual simulation failed: ' ME.message], 'error');
            end
        end
        
        function optimizeAllParameters(app)
            try
                if isempty(app.Data) || isempty(app.SelectedTemp)
                    app.updateStatus('Upload data first!', 'error');
                    return;
                end
                
                app.enableControls('off');
                app.showLoadingDialog('Optimizing Parameters', 'Starting optimization...');
                
                % Initial values
                initial_V0 = app.V0;
                
                % Bounds
                V0_min = 0.001;  % [L/min]
                V0_max = 1;      % [L/min]
                k_min = 0.00001; % [L/mol·min]
                k_max = 10;      % [L/mol·min]
                
                app.updateLoadingDialog(0.2, 'Setting up optimization...');
                
                % Objective function
                objectiveFunc = @(params) app.calculateMaxErrorWithK(params);
                
                % Initial guess
                x0 = [initial_V0, 0.5];
                lb = [V0_min, k_min];
                ub = [V0_max, k_max];
                
                app.updateLoadingDialog(0.5, 'Running optimization...');
                
                % Optimization with enhanced settings for better precision
                options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 500, ...
                                      'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                                      'FunctionTolerance', 1e-10, 'MaxFunctionEvaluations', 2000);
                [x_optimal, max_error] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
                
                V0_optimal = x_optimal(1);
                k_optimal = x_optimal(2);
                
                app.updateLoadingDialog(0.8, 'Applying optimal parameters...');
                
                % Update fields
                app.FlowrateField.Value = V0_optimal;
                app.k_fit = k_optimal;
                
                % Run final simulation
                app.closeLoadingDialog();
                app.runSimulation();
                
                % Store results
                app.OptimizationResults = struct();
                app.OptimizationResults.InitialFlowrate = initial_V0;
                app.OptimizationResults.OptimalFlowrate = V0_optimal;
                app.OptimizationResults.OptimalK = k_optimal;
                app.OptimizationResults.MaxError = max_error;
                app.OptimizationResults.RSquared = app.RSquared;
                
                app.enableControls('on');
                
                % Show results
                if max_error < 1.0
                    emoji = '🎉';
                    status = 'EXCELLENT';
                else
                    emoji = '💡';
                    status = 'GOOD';
                end
                
                msgText = sprintf(['%s Optimization %s!\n\n' ...
                                  'Optimal flowrate: %.4f L/min\n' ...
                                  'Optimal k: %.4f L/mol·min\n' ...
                                  'Maximum error: %.3f%%\n' ...
                                  'R²: %.4f'], ...
                                  emoji, status, V0_optimal, k_optimal, max_error, app.RSquared);
                
                uialert(app.UIFigure, msgText, 'Optimization Results', 'Icon', 'success');
                
            catch ME
                app.closeLoadingDialog();
                app.enableControls('on');
                app.updateStatus('Optimization failed!', 'error');
            end
        end
        
        function optimizeUltraPrecision(app)
    try
        if isempty(app.Data) || isempty(app.SelectedTemp)
            app.updateStatus('Upload data first!', 'error');
            uialert(app.UIFigure, 'Please upload data and select a temperature first.', 'Missing Input');
            return;
        end
        
        app.enableControls('off');
        
        % Create enhanced progress dialog
        d = uiprogressdlg(app.UIFigure, 'Title', '🎯 Ultra-Precision Optimization', ...
                          'Message', 'Initializing ultra-precision optimization for <0.5% error...', ...
                          'Indeterminate', 'off', ...
                          'ShowPercentage', 'on');
        
        % Store initial values
        initial_V0 = app.V0;
        
        % Define enhanced search ranges for ultra precision
        V0_min = 0.000001;  % [L/min] - Even finer minimum
        V0_max = 5;       % [L/min] - Broader maximum range
        k_min = 0.0000001; % [L/mol·min] - Finer minimum for k
        k_max = 100;       % [L/mol·min] - Broader maximum for k
        
        d.Value = 0.1;
        d.Message = 'Setting up ultra-precision optimization problem...';
        
        % Enhanced objective function for ultra precision
        objectiveFunc = @(params) app.calculateMaxErrorWithK(params);
        
        % Multiple starting points for global optimization
        nStartPoints = 5;
        best_result = inf;
        best_x = [];
        
        d.Value = 0.2;
        d.Message = 'Testing multiple starting points for global optimum...';
        
        for i = 1:nStartPoints
            % Generate different starting points
            if i == 1
                x0 = [initial_V0, 0.5]; % Original guess
            else
                x0 = [V0_min + (V0_max - V0_min) * rand(), k_min + (k_max - k_min) * rand()];
            end
            
            % Bounds
            lb = [V0_min, k_min];
            ub = [V0_max, k_max];
            
            d.Value = 0.2 + 0.5 * i / nStartPoints;
            d.Message = sprintf('Running optimization from starting point %d/%d...', i, nStartPoints);
            
            % Ultra-precise optimization options
            options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 1000, ...
                                  'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, ...
                                  'FunctionTolerance', 1e-15, 'MaxFunctionEvaluations', 5000, ...
                                  'Algorithm', 'interior-point');
            
            try
                [x_temp, fval_temp] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
                
                if fval_temp < best_result && isfinite(fval_temp)
                    best_result = fval_temp;
                    best_x = x_temp;
                end
            catch
                % Continue with next starting point if this one fails
                continue;
            end
        end
        
        if isempty(best_x)
            error('All optimization attempts failed');
        end
        
        V0_optimal = best_x(1);
        k_optimal = best_x(2);
        max_error = best_result;
        
        d.Value = 0.8;
        d.Message = 'Running ultra-precision simulation...';
        
        % Store ultra-precision rate constant
        app.k_ultra = k_optimal;
        
        % Run ultra-precision simulation
        T0 = app.SelectedTemp;
        To = T0;
        V = app.VolumeField.Value;
        V0_EtOAc = V0_optimal / 2;
        V1_NaOH = V0_optimal / 2;
        Cin = app.CinField.Value;
        extTime = app.SimTimeField.Value;
        
        % Initial conditions
        y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
        
        % Simulate with ultra-precise parameters
        t_span = linspace(0, extTime, 2000);
        odeOptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-12); % Even tighter tolerances
        [t_sim, y_sim] = ode45(@(t,y) app.odefun_individual(t,y,k_optimal,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                              t_span, y0, odeOptions);
        
        % Store ultra-precision results
        app.UltraPrecisionResults = struct();
        app.UltraPrecisionResults.t = t_sim;
        app.UltraPrecisionResults.C_NaOH = y_sim(:,1);
        app.UltraPrecisionResults.C_EtOAc = y_sim(:,2);
        app.UltraPrecisionResults.C_NaOAc = y_sim(:,3);
        app.UltraPrecisionResults.C_EtOH = y_sim(:,4);
        app.UltraPrecisionResults.T_reactor = y_sim(:,5);
        app.UltraPrecisionResults.T_jacket = y_sim(:,6);
        % DEBUG: Check ultra-precision results storage
fprintf('=== ULTRA-PRECISION RESULTS STORAGE DEBUG ===\n');
fprintf('Time points stored: %d\n', length(t_sim));
fprintf('y_sim dimensions: %dx%d\n', size(y_sim, 1), size(y_sim, 2));
fprintf('NaOH final: %.8f, initial: %.8f\n', y_sim(end,1), y_sim(1,1));
fprintf('EtOAc final: %.8f, initial: %.8f\n', y_sim(end,2), y_sim(1,2));
fprintf('NaOAc final: %.8f, initial: %.8f\n', y_sim(end,3), y_sim(1,3));
fprintf('EtOH final: %.8f, initial: %.8f\n', y_sim(end,4), y_sim(1,4));
fprintf('==============================================\n');
        
        % Detect convergence for ultra-precision
        app.detectConvergence('ultra');
        
        % Calculate R² for ultra-precision and create comparison table
        if ~isempty(app.Data)
            t_exp = app.Data.Time;
            C_exp = app.Data.Concentration;
            C_ultra_interp = interp1(t_sim, y_sim(:,1), t_exp, 'linear', 'extrap');
            app.RSquaredUltra = app.calculateRSquared(C_exp, C_ultra_interp);
            
            % Create comparison table for ultra-precision results
            app.createComparisonTable(t_exp, C_exp, t_sim, y_sim(:,1), 'ultra');
        end
        
        % Update flowrate field with optimal value
        app.FlowrateField.Value = V0_optimal;
        
        % Close progress dialog
        close(d);
        
        % Update all plots to include ultra-precision results
        app.plotConcentrations();
        app.plotTemperatures();
        app.plotProducts();
        app.plotConversions();
        app.updateAnalysisDisplay();
        
        % Store enhanced optimization results
        app.OptimizationResults = struct();
        app.OptimizationResults.InitialFlowrate = initial_V0;
        app.OptimizationResults.OptimalFlowrate = V0_optimal;
        app.OptimizationResults.OptimalK = k_optimal;
        app.OptimizationResults.MaxError = max_error;
        app.OptimizationResults.AllErrorsBelow05Percent = max_error < 0.5;
        app.OptimizationResults.AllErrorsBelow1Percent = max_error < 1.0;
        app.OptimizationResults.RSquared = app.RSquaredUltra;
        app.OptimizationResults.Method = 'Ultra-Precision Multi-Start';
        
        app.enableControls('on');
        
        % Enhanced results display
        if max_error < 0.5
            emoji = '🏆';
            status = 'ULTRA-PRECISE SUCCESS';
            color = 'success';
            msgText = sprintf(['%s %s!\n\n' ...
                              '✨ TARGET ACHIEVED: All errors below 0.5%%\n\n' ...
                              '🎯 Ultra flowrate: %.5f L/min\n' ...
                              '🎯 Each reactant: %.5f L/min\n' ...
                              '🎯 Ultra k: %.6f L/mol·min\n' ...
                              '🎯 Maximum error: %.4f%%\n' ...
                              '🎯 R² correlation: %.6f\n' ...
                              '🎯 Ultra-precision achieved!'], ...
                              emoji, status, V0_optimal, V0_optimal/2, k_optimal, max_error, app.RSquaredUltra);
            app.updateStatus(sprintf('🏆 Ultra-precision success! Max error: %.4f%%', max_error), 'success');
        elseif max_error < 1.0
            emoji = '🎉';
            status = 'HIGH PRECISION SUCCESS';
            color = 'success';
            msgText = sprintf(['%s %s!\n\n' ...
                              '✓ High precision achieved: All errors below 1.0%%\n\n' ...
                              '• Ultra flowrate: %.5f L/min\n' ...
                              '• Each reactant: %.5f L/min\n' ...
                              '• Ultra k: %.6f L/mol·min\n' ...
                              '• Maximum error: %.4f%%\n' ...
                              '• R² correlation: %.6f'], ...
                              emoji, status, V0_optimal, V0_optimal/2, k_optimal, max_error, app.RSquaredUltra);
            app.updateStatus(sprintf('🎉 High precision achieved! Max error: %.4f%%', max_error), 'success');
        else
            emoji = '💡';
            status = 'OPTIMIZATION COMPLETE';
            color = 'warning';
            msgText = sprintf(['%s %s\n\n' ...
                              'Best achievable precision with current data:\n\n' ...
                              '• Best flowrate: %.5f L/min\n' ...
                              '• Each reactant: %.5f L/min\n' ...
                              '• Best k: %.6f L/mol·min\n' ...
                              '• Minimum error: %.4f%%\n' ...
                              '• R² correlation: %.6f\n\n' ...
                              'Consider data quality or model refinements'], ...
                              emoji, status, V0_optimal, V0_optimal/2, k_optimal, max_error, app.RSquaredUltra);
            app.updateStatus(sprintf('💡 Best achievable: %.4f%% error', max_error), 'warning');
        end
        
        % Flash appropriate color
        if max_error < 0.5
            app.flashSuccess(app.OptimizeUltraButton);
        end
        
        uialert(app.UIFigure, msgText, 'Ultra-Precision Results', 'Icon', color);
        
    catch ME
        if exist('d', 'var') && isvalid(d)
            close(d);
        end
        app.enableControls('on');
        app.updateStatus('Ultra-precision optimization failed!', 'error');
        uialert(app.UIFigure, ['Ultra-precision optimization failed: ', ME.message], 'Error');
    end
end
        
        function r_squared = calculateRSquared(~, experimental, simulated)
            try
                valid_idx = ~isnan(experimental) & ~isnan(simulated);
                exp_clean = experimental(valid_idx);
                sim_clean = simulated(valid_idx);
                
                if length(exp_clean) < 2
                    r_squared = NaN;
                    return;
                end
                
                exp_mean = mean(exp_clean);
                ss_tot = sum((exp_clean - exp_mean).^2);
                ss_res = sum((exp_clean - sim_clean).^2);
                
                if ss_tot == 0
                    r_squared = 1;
                else
                    r_squared = 1 - (ss_res / ss_tot);
                end
                
                r_squared = max(0, min(1, r_squared));
                
            catch
                r_squared = NaN;
            end
        end
        
 function checkAndSuggestTimeExtension(app, method)
    try
        if strcmp(method, 'manual')
            isConverged = app.IsManualConverged;
            simType = 'Manual';
            simResults = app.ManualSimResults;
        elseif strcmp(method, 'quicksim')
            isConverged = app.IsQuickSimConverged;
            simType = 'Quick Sim';
            simResults = app.QuickSimResults;
        else
            isConverged = app.IsConverged;
            simType = 'Optimized';
            simResults = app.SimResults;
        end
        
        if ~isConverged
            % Suggest time extension
            currentTime = app.SimTimeField.Value;
            suggestedTime = currentTime * 2;
            
            msg = sprintf(['%s simulation has not fully converged within %.1f minutes.\n\n' ...
                'Some species may still be changing significantly.\n\n' ...
                          'Would you like to extend the simulation time to %.1f minutes?\n' ...
                          'Note: Rate constant will be recalculated from the extension point.'], ...
                          simType, currentTime, suggestedTime);
            
            selection = uiconfirm(app.UIFigure, msg, 'Convergence Warning', ...
                                 'Options', {'Extend Time', 'Keep Current', 'Custom Time'}, ...
                                 'DefaultOption', 1, 'Icon', 'warning');
            
            switch selection
                case 'Extend Time'
                    % Get concentration at current end time for rate constant recalculation
                    if ~isempty(simResults) && isfield(simResults, 't')
                        endConcentration = simResults.C_NaOH(end);
                        app.recalculateForExtension(method, currentTime, endConcentration, suggestedTime);
                    end
                    app.SimTimeField.Value = suggestedTime;
                    app.updateStatus(sprintf('Extended simulation time to %.1f minutes with new rate constant', suggestedTime), 'info');
                    
                case 'Custom Time'
                    answer = inputdlg('Enter new simulation time (minutes):', 'Custom Time', 1, {num2str(suggestedTime)});
                    if ~isempty(answer)
                        try
                            newTime = str2double(answer{1});
                            if isfinite(newTime) && newTime > currentTime
                                % Get concentration at current end time
                                if ~isempty(simResults) && isfield(simResults, 't')
                                    endConcentration = simResults.C_NaOH(end);
                                    app.recalculateForExtension(method, currentTime, endConcentration, newTime);
                                end
                                app.SimTimeField.Value = newTime;
                                app.updateStatus(sprintf('Set simulation time to %.1f minutes with new rate constant', newTime), 'info');
                            end
                        catch
                            app.updateStatus('Invalid time entered', 'error');
                        end
                    end
                case 'Keep Current'
                    app.updateStatus('Keeping current simulation time', 'info');
            end
        end
    catch
        % Silent fail
    end
end
function recalculateForExtension(app, method, originalTime, endConcentration, newTime)
    try
        if strcmp(method, 'manual')
            % Recalculate manual rate constant from extension point
            app.calculateManualRateConstant(originalTime, endConcentration);
            app.updateStatus(sprintf('Manual rate constant recalculated for extension from t=%.1f min', originalTime), 'info');
            
        elseif strcmp(method, 'optimized')
            % Store original values for reference
            app.OriginalSimTime = originalTime;
            app.OriginalEndConcentration = endConcentration;
            app.ExtensionMode = true;
            app.updateStatus(sprintf('Optimization will recalculate for extension from t=%.1f min', originalTime), 'info');
            
        elseif strcmp(method, 'quicksim')
            % For quick sim, we might want to suggest a different rate constant
            % but since it's user-defined, we'll just inform them
            app.updateStatus(sprintf('Quick sim extended. Consider adjusting rate constant for time after %.1f min', originalTime), 'warning');
        end
        
    catch ME
        app.updateStatus(['Extension recalculation failed: ' ME.message], 'error');
    end
end
        
       function createComparisonTable(app, t_exp, C_exp, t_sim, C_sim, method)
    try
        if isempty(t_exp) || isempty(C_exp) || isempty(t_sim) || isempty(C_sim)
            if strcmp(method, 'manual')
                app.ComparisonTableDataManual = table();
            elseif strcmp(method, 'ultra')
                app.ComparisonTableDataUltra = table();
            else
                app.ComparisonTableData = table();
            end
            return;
        end
        
        C_sim_interp = interp1(t_sim, C_sim, t_exp, 'linear', 'extrap');
        C_sim_interp(~isfinite(C_sim_interp)) = 0;
        rel_error = abs(C_exp - C_sim_interp) ./ C_exp * 100;
        rel_error(C_exp == 0) = 0;
        
        if strcmp(method, 'manual')
            app.ComparisonTableDataManual = table(t_exp, C_exp, C_sim_interp, rel_error, ...
                'VariableNames', {'Time_min', 'Experimental_mol_L', 'Manual_Simulated_mol_L', 'Manual_Error_percent'});
        elseif strcmp(method, 'ultra')
            app.ComparisonTableDataUltra = table(t_exp, C_exp, C_sim_interp, rel_error, ...
                'VariableNames', {'Time_min', 'Experimental_mol_L', 'Ultra_Simulated_mol_L', 'Ultra_Error_percent'});
        else
            app.ComparisonTableData = table(t_exp, C_exp, C_sim_interp, rel_error, ...
                'VariableNames', {'Time_min', 'Experimental_mol_L', 'Optimized_Simulated_mol_L', 'Optimized_Error_percent'});
        end
        
        % Update table display with combined data
        app.updateComparisonTableDisplay();
        
    catch
        if strcmp(method, 'manual')
            app.ComparisonTableDataManual = table();
        elseif strcmp(method, 'ultra')
            app.ComparisonTableDataUltra = table();
        else
            app.ComparisonTableData = table();
        end
    end
       end
        
        function updateComparisonTableDisplay(app)
    try
        % Combine all comparison data including ultra-precision
        if ~isempty(app.ComparisonTableDataManual) && ~isempty(app.ComparisonTableData) && ~isempty(app.ComparisonTableDataUltra)
            % All three methods exist
            combinedData = [app.ComparisonTableDataManual.Time_min, ...
                           app.ComparisonTableDataManual.Experimental_mol_L, ...
                           app.ComparisonTableDataManual.Manual_Simulated_mol_L, ...
                           app.ComparisonTableDataManual.Manual_Error_percent, ...
                           app.ComparisonTableData.Optimized_Simulated_mol_L, ...
                           app.ComparisonTableData.Optimized_Error_percent, ...
                           app.ComparisonTableDataUltra.Ultra_Simulated_mol_L, ...
                           app.ComparisonTableDataUltra.Ultra_Error_percent];
            app.ComparisonTable.Data = combinedData;
            app.ComparisonTable.ColumnName = {'Time', 'Exp', 'Manual Sim', 'Manual Error%', 'Opt Sim', 'Opt Error%', 'Ultra Sim', 'Ultra Error%'};
        elseif ~isempty(app.ComparisonTableDataManual) && ~isempty(app.ComparisonTableData)
            % Manual and optimized exist
            combinedData = [app.ComparisonTableDataManual.Time_min, ...
                           app.ComparisonTableDataManual.Experimental_mol_L, ...
                           app.ComparisonTableDataManual.Manual_Simulated_mol_L, ...
                           app.ComparisonTableDataManual.Manual_Error_percent, ...
                           app.ComparisonTableData.Optimized_Simulated_mol_L, ...
                           app.ComparisonTableData.Optimized_Error_percent];
            app.ComparisonTable.Data = combinedData;
            app.ComparisonTable.ColumnName = {'Time', 'Exp', 'Manual Sim', 'Manual Error%', 'Opt Sim', 'Opt Error%'};
        elseif ~isempty(app.ComparisonTableDataUltra)
            % Only ultra-precision data exists
            app.ComparisonTable.Data = app.ComparisonTableDataUltra;
            app.ComparisonTable.ColumnName = {'Time', 'Exp', 'Ultra Sim', 'Ultra Error%'};
        elseif ~isempty(app.ComparisonTableDataManual)
            % Only manual data exists
            app.ComparisonTable.Data = app.ComparisonTableDataManual;
            app.ComparisonTable.ColumnName = {'Time', 'Exp', 'Manual Sim', 'Manual Error%'};
        elseif ~isempty(app.ComparisonTableData)
            % Only optimized data exists
            app.ComparisonTable.Data = app.ComparisonTableData;
            app.ComparisonTable.ColumnName = {'Time', 'Exp', 'Opt Sim', 'Opt Error%'};
        else
            % No data
            app.ComparisonTable.Data = [];
            app.ComparisonTable.ColumnName = {};
        end
        
    catch
        % Silent fail
    end
end
        
        % Plotting methods with enhanced convergence markers
        function plotConcentrations(app)
            try
                cla(app.UIAxes1);
                hold(app.UIAxes1, 'on');
                
                legendEntries = {};
                
                % Plot experimental data first (if available)
                if ~isempty(app.Data)
                    scatter(app.UIAxes1, app.Data.Time, app.Data.Concentration, 80, ...
                           'MarkerEdgeColor', app.ErrorColor, ...
                           'MarkerFaceColor', app.BackgroundColor, ...
                           'LineWidth', 2.5, 'DisplayName', 'Experimental Data');
                    legendEntries{end+1} = 'Experimental Data';
                end
                
                % Plot manual simulation results (if available)
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                    plot(app.UIAxes1, app.ManualSimResults.t, app.ManualSimResults.C_NaOH, ':', ...
                         'LineWidth', 3, 'Color', [147 112 219]/255, 'DisplayName', 'Manual k Simulation');
                    legendEntries{end+1} = 'Manual k Simulation';
                    
                    % Add convergence markers for manual simulation
                    if ~isempty(app.ManualConvergenceData)
                        app.addConvergenceMarkers(app.UIAxes1, app.ManualSimResults, app.ManualConvergenceData, [147 112 219]/255, 'Manual');
                    end
                end
                
                % Plot optimized simulation results (if available)
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    plot(app.UIAxes1, app.SimResults.t, app.SimResults.C_NaOH, '-', ...
                         'LineWidth', 3, 'Color', app.AccentColor, 'DisplayName', 'Optimized Simulation');
                    legendEntries{end+1} = 'Optimized Simulation';
                    
                    % Add convergence markers for optimized simulation
                    if ~isempty(app.ConvergenceData)
                        app.addConvergenceMarkers(app.UIAxes1, app.SimResults, app.ConvergenceData, app.AccentColor, 'Optimized');
                    end
                end
                 
% Plot ultra-precision simulation results (if available)
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
    plot(app.UIAxes1, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, '-', ...
         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision Simulation');
    legendEntries{end+1} = 'Ultra-Precision Simulation';
    
    % ADD CONVERGENCE MARKERS FOR ULTRA-PRECISION SIMULATION
    if ~isempty(app.UltraConvergenceData)
        app.addConvergenceMarkers(app.UIAxes1, app.UltraPrecisionResults, app.UltraConvergenceData, [255 20 147]/255, 'Ultra-Precision');
    end
end


if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't') && app.ShowProducts
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra NaOAc');
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
         'LineWidth', 4, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra EtOH');
    
    % Add convergence markers for ultra-precision products
    if ~isempty(app.UltraConvergenceData)
        app.addProductConvergenceMarkers(app.UIAxes3, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra');
    end
end
                
                % Enhanced styling with error information
                if ~isnan(app.k_fit) && app.k_fit > 0
                    if ~isempty(app.ComparisonTableData)
                        max_error = max(app.ComparisonTableData.Optimized_Error_percent);
                        titleStr = sprintf('CSTR Concentration Profiles (k_{opt} = %.4f L/mol·min, Max Error: %.3f%%)', app.k_fit, max_error);
                    else
                        titleStr = sprintf('CSTR Concentration Profiles (k_{opt} = %.4f L/mol·min)', app.k_fit);
                    end
                elseif ~isnan(app.k_manual) && app.k_manual > 0
                    titleStr = sprintf('CSTR Concentration Profiles (k_{manual} = %.4f L/mol·min)', app.k_manual);
                else
                    titleStr = 'CSTR Concentration Profiles';
                end
                
                title(app.UIAxes1, titleStr, 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(app.UIAxes1, 'Time (min)', 'FontWeight', 'bold');
                ylabel(app.UIAxes1, 'NaOH Concentration (mol/L)', 'FontWeight', 'bold');
                
                % Only show legend if there are entries
                if ~isempty(legendEntries)
                    legend(app.UIAxes1, 'Location', 'best', 'Box', 'on', 'FontSize', 11);
                end
                
                grid(app.UIAxes1, 'off');
                app.UIAxes1.GridAlpha = 0.3;
                app.UIAxes1.Box = 'on';
                app.UIAxes1.FontSize = 11;
                
                hold(app.UIAxes1, 'off');
                
                % Update info panel
                app.updateConcentrationInfo();
                
            catch ME
                app.updateStatus(['Plot error: ' ME.message], 'error');
            end
            
            
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't') && app.ShowProducts
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra NaOAc');
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
         'LineWidth', 4, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra EtOH');
    
    % Add convergence markers for ultra-precision products
    if ~isempty(app.UltraConvergenceData)
        app.addProductConvergenceMarkers(app.UIAxes3, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra');
    end
end
        end
        
   % REPLACE the updateConcentrationInfo function with this enhanced version:

function updateConcentrationInfo(app)
    try
        infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
        infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">Performance Metrics & Enhanced Convergence</h3>'];
        
        % Performance metrics section
        if ~isnan(app.RSquared) && app.RSquared > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Optimized R²:</strong> %.6f</p>', app.RSquared)];
        end
        
        if ~isnan(app.RSquaredManual) && app.RSquaredManual > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Manual R²:</strong> %.6f</p>', app.RSquaredManual)];
        end
        
        % NEW: Ultra-precision R² if available
        if ~isnan(app.RSquaredUltra) && app.RSquaredUltra > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #FF1493;"><strong>🏆 Ultra-Precision R²:</strong> %.8f</p>', app.RSquaredUltra)];
        end
        
        if ~isnan(app.k_fit) && app.k_fit > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Optimized k:</strong> %.6f L/mol·min</p>', app.k_fit)];
        end
        
        if ~isnan(app.k_manual) && app.k_manual > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Manual k:</strong> %.6f L/mol·min</p>', app.k_manual)];
        end
        
        % NEW: Ultra-precision k if available
        if ~isnan(app.k_ultra) && app.k_ultra > 0
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #FF1493;"><strong>🏆 Ultra-Precision k:</strong> %.8f L/mol·min</p>', app.k_ultra)];
        end
        
        if ~isempty(app.Data)
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Data Points:</strong> %d</p>', height(app.Data))];
        end
        
        % Enhanced error information
        if ~isempty(app.ComparisonTableData)
            max_error = max(app.ComparisonTableData.Optimized_Error_percent);
            mean_error = mean(app.ComparisonTableData.Optimized_Error_percent);
            if max_error < 0.5
                errorColor = '#34C759'; % Green
                errorIcon = '🏆';
            elseif max_error < 1.0
                errorColor = '#FF9500'; % Orange
                errorIcon = '🎯';
            else
                errorColor = '#FF3B30'; % Red
                errorIcon = '⚠';
            end
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: %s;"><strong>%s Optimized Max Error:</strong> %.4f%%</p>', errorColor, errorIcon, max_error)];
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: %s;"><strong>Optimized Mean Error:</strong> %.4f%%</p>', errorColor, mean_error)];
        end
        
        % NEW: Ultra-precision error information
        if ~isempty(app.ComparisonTableDataUltra)
            max_error_ultra = max(app.ComparisonTableDataUltra.Ultra_Error_percent);
            mean_error_ultra = mean(app.ComparisonTableDataUltra.Ultra_Error_percent);
            if max_error_ultra < 0.1
                errorColorUltra = '#32CD32'; % Lime Green
                errorIconUltra = '🏆🏆';
            elseif max_error_ultra < 0.5
                errorColorUltra = '#34C759'; % Green
                errorIconUltra = '🏆';
            else
                errorColorUltra = '#FF9500'; % Orange
                errorIconUltra = '🎯';
            end
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: %s;"><strong>%s Ultra-Precision Max Error:</strong> %.6f%%</p>', errorColorUltra, errorIconUltra, max_error_ultra)];
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: %s;"><strong>Ultra-Precision Mean Error:</strong> %.6f%%</p>', errorColorUltra, mean_error_ultra)];
        end
        
        % Enhanced convergence status with separate tracking including Ultra-Precision
        infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
        infoHTML = [infoHTML '<strong>Multi-Level Convergence Status:</strong><br>'];
        
        % Manual convergence status
        if app.ManualReactantsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Manual Reactants: Converged (Relaxed)</p>'];
        elseif ~isempty(app.ManualSimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Manual Reactants: Not Converged</p>'];
        end
        
        if app.ManualProductsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Manual Products: Converged (Relaxed)</p>'];
        elseif ~isempty(app.ManualSimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Manual Products: Not Converged</p>'];
        end
        
        % Optimized convergence status
        if app.OptimizedReactantsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Optimized Reactants: Converged (Standard)</p>'];
        elseif ~isempty(app.SimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Optimized Reactants: Not Converged</p>'];
        end
        
        if app.OptimizedProductsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Optimized Products: Converged (Standard)</p>'];
        elseif ~isempty(app.SimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Optimized Products: Not Converged</p>'];
        end
        
        % NEW: Ultra-precision convergence status
        if app.UltraReactantsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #32CD32;"><strong>🏆 Ultra-Precision Reactants: Converged (Strict)</strong></p>'];
        elseif ~isempty(app.UltraPrecisionResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF1493;">⚠ Ultra-Precision Reactants: Not Converged</p>'];
        end
        
        if app.UltraProductsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #32CD32;"><strong>🏆 Ultra-Precision Products: Converged (Strict)</strong></p>'];
        elseif ~isempty(app.UltraPrecisionResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF1493;">⚠ Ultra-Precision Products: Not Converged</p>'];
        end
        
        if app.UltraTemperatureConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #32CD32;"><strong>🏆 Ultra-Precision Temperature: Converged (Strict)</strong></p>'];
        elseif ~isempty(app.UltraPrecisionResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF1493;">⚠ Ultra-Precision Temperature: Not Converged</p>'];
        end
        
        % Quick sim convergence status
        if app.QuickSimReactantsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Quick Sim Reactants: Converged (Medium)</p>'];
        elseif ~isempty(app.QuickSimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Quick Sim Reactants: Not Converged</p>'];
        end
        
        if app.QuickSimProductsConverged
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #34C759;">✓ Quick Sim Products: Converged (Medium)</p>'];
        elseif ~isempty(app.QuickSimResults)
            infoHTML = [infoHTML '<p style="margin: 2px 0; color: #FF9500;">⚠ Quick Sim Products: Not Converged</p>'];
        end
        
        % NEW: Add convergence criteria summary
        infoHTML = [infoHTML '<hr style="margin: 8px 0; border: none; border-top: 1px solid #ddd;">'];
        infoHTML = [infoHTML '<strong>Convergence Criteria Summary:</strong><br>'];
        infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px;">Manual (Relaxed): React %.3f%%, Prod %.3f%%</p>', app.ManualReactantThreshold*100, app.ManualProductThreshold*100)];
        infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px;">Optimized (Standard): React %.3f%%, Prod %.3f%%</p>', app.OptimizedReactantThreshold*100, app.OptimizedProductThreshold*100)];
        if ~isempty(app.UltraPrecisionResults)
            infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px; color: #FF1493;"><strong>Ultra-Precision (Strict): React %.4f%%, Prod %.5f%%</strong></p>', app.UltraReactantThreshold*100, app.UltraProductThreshold*100)];
        end
        infoHTML = [infoHTML sprintf('<p style="margin: 1px 0; font-size: 9px;">Quick Sim (Medium): React %.3f%%, Prod %.3f%%</p>', app.QuickSimReactantThreshold*100, app.QuickSimProductThreshold*100)];
        
        infoHTML = [infoHTML '</div>'];
        app.ConcentrationInfo.HTMLSource = infoHTML;
        
    catch
        % Silent fail
    end
end
        
     function plotTemperatures(app)
    try
        cla(app.UIAxes2);
        hold(app.UIAxes2, 'on');
        
        % Plot temperature profiles with enhanced visualization
        legendEntries = {};
        
        % Plot manual simulation temperatures (if available)
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            plot(app.UIAxes2, app.ManualSimResults.t, app.ManualSimResults.T_reactor, ':', ...
                 'LineWidth', 3, 'Color', [255 102 102]/255, 'DisplayName', 'Manual Reactor T');
            legendEntries{end+1} = 'Manual Reactor T';
            
            if app.ShowJacket
                plot(app.UIAxes2, app.ManualSimResults.t, app.ManualSimResults.T_jacket, ':', ...
                     'LineWidth', 2.5, 'Color', [102 178 255]/255, 'DisplayName', 'Manual Jacket T');
                legendEntries{end+1} = 'Manual Jacket T';
            end
        end
        
        % Re-simulate using current kinetics for temperature profile visualization
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            % Use current simulation parameters
            T0 = app.SelectedTemp;
            To = T0;
            V = app.VolumeField.Value;
            V0 = app.FlowrateField.Value;
            V0_EtOAc = V0 / 2;
            V1_NaOH = V0 / 2;
            Cin = app.CinField.Value;
            extTime = app.SimTimeField.Value;
            k_current = app.k_fit;
            
            % Initial conditions
            y0 = [Cin, Cin, 0, 0, T0, app.Tj0];
            
            % Simulate with current kinetics for temperature profile
            t_span = linspace(0, extTime, 2000);
            options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
            [t_current, y_current] = ode45(@(t,y) app.odefun_individual(t,y,k_current,V0_EtOAc,V1_NaOH,V,Cin,To), ...
                                          t_span, y0, options);
            
            % Temperature data with current kinetics
            x = t_current;
            y_reactor = y_current(:,5);
            y_jacket = y_current(:,6);
            
            % Enhanced visualization based on isothermal mode
            if app.IsothermalMode
                % Isothermal mode: Show controlled temperature with gradients
                fill_color = app.SuccessColor;  % Green for controlled
                line_style = '-';
                line_width = 3;
                temp_label = 'Reactor T (Isothermal)';
            else
                % Non-isothermal mode: Show temperature rise with warning colors
                fill_color = app.WarningColor;  % Orange/red for rising temperature
                line_style = '-';
                line_width = 4;  % Thicker line to emphasize temperature rise
                temp_label = 'Reactor T (Rising!)';
            end
            
            % Gradient fill for reactor temperature
            fill_x = [x; flipud(x)];
            fill_y = [y_reactor; ones(size(y_reactor))*min(y_reactor)];
            fill(app.UIAxes2, fill_x, fill_y, fill_color, ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Temperature Range');
            
            % Main temperature lines with enhanced styling
            plot(app.UIAxes2, x, y_reactor, line_style, 'LineWidth', line_width, 'Color', app.ErrorColor, ...
                 'DisplayName', temp_label);
            legendEntries{end+1} = temp_label;
            
       % Jacket temperature if enabled
if app.ShowJacket
    if app.IsothermalMode
        jacket_color = app.AccentColor;
    else
        jacket_color = [255 140 0]/255; % Orange if not controlled
    end
    plot(app.UIAxes2, x, y_jacket, '-', ...
        'LineWidth', 2.5, 'Color', jacket_color, 'DisplayName', 'Jacket T');
    legendEntries{end+1} = 'Jacket T';
end
            
            % Temperature rise markers and annotations
            max_temp = max(y_reactor);
            temp_rise = max_temp - y_reactor(1);
            [max_temp_val, max_idx] = max(y_reactor);
            
            % Peak temperature marker
            scatter(app.UIAxes2, x(max_idx), max_temp_val, 100, 'filled', ...
                   'MarkerFaceColor', app.ErrorColor, 'MarkerEdgeColor', 'none', ...
                   'DisplayName', sprintf('Peak: %.1f°C', max_temp_val));
            legendEntries{end+1} = sprintf('Peak: %.1f°C', max_temp_val);
            
            % Add temperature rise annotation
            if temp_rise > 1  % Only show if significant rise
                annotation_text = sprintf('+%.1f°C', temp_rise);
                text(app.UIAxes2, x(max_idx), max_temp_val + 0.5, annotation_text, ...
                     'FontSize', 12, 'FontWeight', 'bold', 'Color', app.ErrorColor, ...
                     'HorizontalAlignment', 'center');
            end
            
            % Reference lines
            yline(app.UIAxes2, app.SelectedTemp, ':', 'Color', app.SecondaryTextColor, ...
                  'LineWidth', 1.5, 'Label', 'Setpoint', 'LabelHorizontalAlignment', 'left');
            
            if ~app.IsothermalMode && temp_rise > 5
                % Add warning line for significant temperature rise
                yline(app.UIAxes2, app.SelectedTemp + 10, '--', 'Color', app.ErrorColor, ...
                      'LineWidth', 1, 'Label', 'Caution Zone', 'LabelHorizontalAlignment', 'right');
            end
        end
        
        % Enhanced title with mode information
        if app.IsothermalMode
            mode_status = '🌡️ ISOTHERMAL CONTROL';
        else
            mode_status = '🔥 NON-ISOTHERMAL (Temperature Rising)';
        end
        
        titleStr = sprintf('%s | Enhanced Temperature Control', mode_status);
        title(app.UIAxes2, titleStr, 'FontWeight', 'bold', 'FontSize', 14);
        xlabel(app.UIAxes2, 'Time (min)', 'FontWeight', 'bold');
        ylabel(app.UIAxes2, 'Temperature (°C)', 'FontWeight', 'bold');
        
        if ~isempty(legendEntries)
            legend(app.UIAxes2, 'Location', 'best', 'Box', 'on', 'FontSize', 11);
        end
        
        grid(app.UIAxes2, 'off');
        app.UIAxes2.GridAlpha = 0.3;
        app.UIAxes2.Box = 'on';
        app.UIAxes2.FontSize = 11;
        
        hold(app.UIAxes2, 'off');
        
        % Update enhanced info panel
        app.updateTemperatureInfoEnhanced();
        
    catch ME
        app.updateStatus(['Temperature plot error: ' ME.message], 'error');
    end
end
        
        function updateTemperatureInfoArrhenius(app)
    try
        infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
        infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">Enhanced Temperature Analysis with Arrhenius Kinetics</h3>'];
        
        % Arrhenius parameters section
        infoHTML = [infoHTML '<div style="background: #E3F2FD; padding: 8px; border-radius: 6px; margin-bottom: 10px;">'];
        infoHTML = [infoHTML '<h4 style="margin: 0 0 5px 0; color: #1976D2;">⚗️ Arrhenius Parameters</h4>'];
        infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;"><strong>Activation Energy:</strong> %.1f kJ/mol</p>', app.Ea/1000)];
        infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;"><strong>Pre-exponential Factor:</strong> %.2e L/mol·min</p>', app.A_preexp)];
        
        % Calculate rate constant at current temperature
        k_current = app.calculateArrheniusRateConstant(app.SelectedTemp);
        infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;"><strong>k at %.0f°C:</strong> %.4f L/mol·min</p>', app.SelectedTemp, k_current)];
        infoHTML = [infoHTML '</div>'];
        
        % Temperature analysis for Arrhenius simulation
        if ~isempty(app.SimResults) && isfield(app.SimResults, 'T_reactor')
            y_reactor = app.SimResults.T_reactor;
            y_jacket = app.SimResults.T_jacket;
            
            max_temp = max(y_reactor);
            temp_rise = max_temp - y_reactor(1);
            final_temp = y_reactor(end);
            
            % Calculate rate constants at key temperatures
            k_initial = app.calculateArrheniusRateConstant(y_reactor(1));
            k_max = app.calculateArrheniusRateConstant(max_temp);
            k_final = app.calculateArrheniusRateConstant(final_temp);
            
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Peak Temperature:</strong> %.3f°C (k=%.4f)</p>', max_temp, k_max)];
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Temperature Rise:</strong> +%.3f°C</p>', temp_rise)];
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Final Temperature:</strong> %.3f°C (k=%.4f)</p>', final_temp, k_final)];
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Rate Constant Ratio (max/initial):</strong> %.2f</p>', k_max/k_initial)];
            
            % Add jacket temperature info if shown
            if app.ShowJacket
                temp_diff = max(y_reactor) - max(y_jacket);
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Max T Difference:</strong> %.3f°C</p>', temp_diff)];
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Final Jacket T:</strong> %.3f°C</p>', y_jacket(end))];
            end
        end
        
        % Add note about Arrhenius kinetics
        infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
        infoHTML = [infoHTML '<p style="margin: 5px 0; font-style: italic; color: #666; font-size: 10px;">'];
        infoHTML = [infoHTML '📝 Temperature profile calculated using Arrhenius equation: k = A·exp(-Ea/RT)</p>'];
        
        infoHTML = [infoHTML '</div>'];
        app.TemperatureInfo.HTMLSource = infoHTML;
        
    catch
        % Silent fail
    end
end

        
        function updateTemperatureInfo(app)
            try
                infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
                infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">Enhanced Temperature Analysis</h3>'];
                
                % Optimized simulation temperature info
                if ~isempty(app.SimResults) && isfield(app.SimResults, 'T_reactor')
                    y_reactor = app.SimResults.T_reactor;
                    y_jacket = app.SimResults.T_jacket;
                    
                    max_temp = max(y_reactor);
                    temp_rise = max_temp - y_reactor(1);
                    final_temp = y_reactor(end);
                    
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Peak Temperature:</strong> %.3f°C</p>', max_temp)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Temperature Rise:</strong> +%.3f°C</p>', temp_rise)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Final Temperature:</strong> %.3f°C</p>', final_temp)];
                    
                    % Add jacket temperature info if shown
                    if app.ShowJacket
                        temp_diff = max(y_reactor) - max(y_jacket);
                        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Max T Difference:</strong> %.3f°C</p>', temp_diff)];
                        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Final Jacket T:</strong> %.3f°C</p>', y_jacket(end))];
                    end
                    
                    % Add convergence info with different criteria
                    if app.OptimizedTemperatureConverged
                        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #34C759;"><strong>Opt T Convergence:</strong> %.1f min</p>', app.ConvergenceData.T_time)];
                    else
                        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500;"><strong>Opt T Status:</strong> Not Converged</p>'];
                    end
                end
                
                % Manual simulation temperature comparison
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 'T_reactor')
                    manual_max_temp = max(app.ManualSimResults.T_reactor);
                    manual_temp_rise = manual_max_temp - app.ManualSimResults.T_reactor(1);
                    
                    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Manual Peak T:</strong> %.3f°C</p>', manual_max_temp)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Manual T Rise:</strong> +%.3f°C</p>', manual_temp_rise)];
                    
                    % Manual temperature convergence
                    if app.ManualTemperatureConverged
                        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #34C759;"><strong>Manual T Convergence:</strong> %.1f min</p>', app.ManualConvergenceData.T_time)];
                    else
                        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500;"><strong>Manual T Status:</strong> Not Converged</p>'];
                    end
                end
                
                % Quick sim temperature comparison
                if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 'T_reactor')
                    quicksim_max_temp = max(app.QuickSimResults.T_reactor);
                    quicksim_temp_rise = quicksim_max_temp - app.QuickSimResults.T_reactor(1);
                    
                    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Quick Sim Peak T:</strong> %.3f°C</p>', quicksim_max_temp)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Quick Sim T Rise:</strong> +%.3f°C</p>', quicksim_temp_rise)];
                    
                    % Quick sim temperature convergence
                    if app.QuickSimTemperatureConverged
                        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #34C759;"><strong>Quick Sim T Convergence:</strong> %.1f min</p>', app.QuickSimConvergenceData.T_time)];
                    else
                        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500;"><strong>Quick Sim T Status:</strong> Not Converged</p>'];
                    end
                end
                % Add this section after the Quick sim temperature comparison:

% Ultra-precision simulation temperature comparison
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 'T_reactor')
    ultra_max_temp = max(app.UltraPrecisionResults.T_reactor);
    ultra_temp_rise = ultra_max_temp - app.UltraPrecisionResults.T_reactor(1);
    
    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Ultra-Precision Peak T:</strong> %.3f°C</p>', ultra_max_temp)];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Ultra-Precision T Rise:</strong> +%.3f°C</p>', ultra_temp_rise)];
    
    % Ultra-precision temperature convergence
    if app.UltraTemperatureConverged
        infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #34C759;"><strong>Ultra-Precision T Convergence:</strong> %.1f min</p>', app.UltraConvergenceData.T_time)];
    else
        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500;"><strong>Ultra-Precision T Status:</strong> Not Converged</p>'];
    end
end
if app.IsothermalMode && ~isempty(app.SimResults)
    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>Cooling System:</strong></p>')];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Adaptive Flow: %.3f L/min</p>', app.Fj_adaptive)];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Adequacy: %.2fx</p>', app.cooling_adequacy_factor)];
    
    if app.cooling_adequacy_factor > 1.2
        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #34C759;">Status: ✅ Adequate</p>'];
    elseif app.cooling_adequacy_factor > 0.8
        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF9500;">Status: ⚠️ Marginal</p>'];
    else
        infoHTML = [infoHTML '<p style="margin: 5px 0; color: #FF3B30;">Status: ❌ Inadequate</p>'];
    end
end
                
                infoHTML = [infoHTML '</div>'];
                app.TemperatureInfo.HTMLSource = infoHTML;
                
            catch
                % Silent fail
            end
        end
        
        % ENHANCED TEMPERATURE INFO WITH ADAPTIVE FLOWRATE DISPLAY
function updateTemperatureInfoEnhanced(app)
    try
        infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
        infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">🌡️ Temperature-Responsive Cooling Analysis</h3>'];
        
        % Display current control mode with flowrate info
        if app.IsothermalMode
            infoHTML = [infoHTML '<div style="background: #E3F2FD; padding: 8px; border-radius: 6px; margin-bottom: 10px;">'];
            infoHTML = [infoHTML '<h4 style="margin: 0 0 5px 0; color: #1976D2;">✅ ADAPTIVE ISOTHERMAL CONTROL</h4>'];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 12px; font-weight: bold; color: #1976D2;">🌊 Current Adaptive Flowrate: %.3f L/min</p>', app.Fj_adaptive)];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">📊 Flowrate Utilization: %.1f%% of maximum</p>', (app.Fj_adaptive/app.Fj_max)*100)];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">🎯 Target Temperature: %.1f°C</p>', app.SelectedTemp)];
            infoHTML = [infoHTML '</div>'];
        else
            infoHTML = [infoHTML '<div style="background: #FFEBEE; padding: 8px; border-radius: 6px; margin-bottom: 10px;">'];
            infoHTML = [infoHTML '<h4 style="margin: 0 0 5px 0; color: #D32F2F;">❌ NO COOLING SYSTEM</h4>'];
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 11px; color: #D32F2F;"><strong>🌊 Adaptive Flowrate: 0.000 L/min (DISCONNECTED)</strong></p>'];
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 10px; color: #FF5722;">⚠️ All reaction heat accumulates - risk of thermal runaway!</p>'];
            infoHTML = [infoHTML '</div>'];
        end
        
        % Temperature analysis for current simulation
        if ~isempty(app.SimResults) && isfield(app.SimResults, 'T_reactor')
            y_reactor = app.SimResults.T_reactor;
            y_jacket = app.SimResults.T_jacket;
            
            max_temp = max(y_reactor);
            temp_rise = max_temp - y_reactor(1);
            final_temp = y_reactor(end);
            
            % Display flowrate progression if available
            if ~isempty(app.CoolingFlowrateHistory) && app.IsothermalMode
                avg_flowrate = mean(app.CoolingFlowrateHistory);
                max_flowrate_used = max(app.CoolingFlowrateHistory);
                min_flowrate_used = min(app.CoolingFlowrateHistory);
                
                infoHTML = [infoHTML '<div style="background: #E8F5E8; padding: 8px; border-radius: 6px; margin-bottom: 10px;">'];
                infoHTML = [infoHTML '<h4 style="margin: 0 0 5px 0; color: #2E7D32;">📈 Adaptive Flowrate Performance</h4>'];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;"><strong>🌊 Average Flowrate Used: %.3f L/min</strong></p>', avg_flowrate)];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">📊 Peak Flowrate: %.3f L/min (%.1f%% of max)</p>', max_flowrate_used, (max_flowrate_used/app.Fj_max)*100)];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">📉 Minimum Flowrate: %.3f L/min</p>', min_flowrate_used)];
                infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">🎚️ Flowrate Range: %.3f L/min</p>', max_flowrate_used - min_flowrate_used)];
                infoHTML = [infoHTML '</div>'];
            end
            
            % Temperature performance with flowrate context
            if app.IsothermalMode
                temp_control_quality = 100 * sum(abs(y_reactor - app.SelectedTemp) <= app.isothermal_tolerance) / length(y_reactor);
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>🎯 Temperature Control Quality:</strong> %.1f%%</p>', temp_control_quality)];
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>📊 Max Deviation:</strong> %.3f°C</p>', max(abs(y_reactor - app.SelectedTemp)))];
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>❄️ Cooling Adequacy:</strong> %.2fx</p>', app.cooling_adequacy_factor)];
                
                % Flowrate efficiency analysis
                if ~isempty(app.CoolingFlowrateHistory)
                    flowrate_efficiency = temp_control_quality / ((avg_flowrate/app.Fj_max)*100) * 100;
                    if flowrate_efficiency > 150
                        efficiency_color = '#2E7D32'; % Green - very efficient
                        efficiency_text = 'EXCELLENT';
                    elseif flowrate_efficiency > 100
                        efficiency_color = '#388E3C'; % Light green - good
                        efficiency_text = 'GOOD';
                    elseif flowrate_efficiency > 50
                        efficiency_color = '#F57C00'; % Orange - fair
                        efficiency_text = 'FAIR';
                    else
                        efficiency_color = '#D32F2F'; % Red - poor
                        efficiency_text = 'POOR';
                    end
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: %s;"><strong>⚡ Cooling Efficiency:</strong> %s (%.0f%%)</p>', efficiency_color, efficiency_text, flowrate_efficiency)];
                end
            else
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #D32F2F;"><strong>🔥 Temperature Rise:</strong> +%.2f°C</p>', temp_rise)];
                infoHTML = [infoHTML sprintf('<p style="margin: 5px 0; color: #FF5722;"><strong>🌡️ Peak Temperature:</strong> %.2f°C</p>', max_temp)];
                
                % Safety warnings for high temperatures
                if temp_rise > 15
                    infoHTML = [infoHTML '<p style="margin: 8px 0; padding: 6px; background: #FFCDD2; border-left: 4px solid #F44336; color: #D32F2F; font-weight: bold; font-size: 10px;">'];
                    infoHTML = [infoHTML '🚨 DANGER: Extreme temperature rise! Immediate cooling required!</p>'];
                elseif temp_rise > 10
                    infoHTML = [infoHTML '<p style="margin: 8px 0; padding: 6px; background: #FFF3E0; border-left: 4px solid #FF9800; color: #F57C00; font-weight: bold; font-size: 10px;">'];
                    infoHTML = [infoHTML '⚠️ WARNING: Significant temperature rise detected!</p>'];
                end
            end
            
            infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;"><strong>🏁 Final Temperature:</strong> %.2f°C</p>', final_temp)];
        end
        
        % Cooling system capacity summary with temperature-responsive info
        infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
        if app.IsothermalMode
            infoHTML = [infoHTML '<strong>🔧 Temperature-Responsive Cooling System:</strong><br>'];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">• Flowrate range: %.2f - %.1f L/min</p>', app.Fj_min, app.Fj_max)];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px; color: #1976D2;"><strong>• Current adaptive rate: %.3f L/min</strong></p>', app.Fj_adaptive)];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">• Coolant inlet: %.1f°C</p>', app.Tj0)];
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">• Control tolerance: ±%.1f°C</p>', app.isothermal_tolerance)];
            
            % Temperature-specific base flowrates
            temp_setpoint = app.SelectedTemp;
            if temp_setpoint <= 30
                base_info = '0.8 L/min (Low temp base)';
            elseif temp_setpoint <= 40
                base_info = '1.2 L/min (Medium temp base)';
            elseif temp_setpoint <= 50
                base_info = '1.8 L/min (High temp base)';
            elseif temp_setpoint <= 60
                base_info = '2.5 L/min (Very high temp base)';
            else
                base_info = '3.5 L/min (Extreme temp base)';
            end
            infoHTML = [infoHTML sprintf('<p style="margin: 2px 0; font-size: 11px;">• Base flowrate for %.0f°C: %s</p>', temp_setpoint, base_info)];
        else
            infoHTML = [infoHTML '<strong style="color: #D32F2F;">🚫 No Cooling System:</strong><br>'];
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 11px; color: #D32F2F;">• Jacket: COMPLETELY DISCONNECTED</p>'];
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 11px; color: #D32F2F;">• Adaptive flowrate: 0.000 L/min</p>'];
            infoHTML = [infoHTML '<p style="margin: 2px 0; font-size: 11px; color: #FF5722;">• Temperature will rise significantly!</p>'];
        end
        
        infoHTML = [infoHTML '</div>'];
        app.TemperatureInfo.HTMLSource = infoHTML;
        
    catch
        app.TemperatureInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">Temperature info unavailable.</div>';
    end
end
        
        % Products plot with enhanced convergence markers
        function plotProducts(app)
            try
                cla(app.UIAxes3);
                hold(app.UIAxes3, 'on');
                
                % Plot products from both manual and optimized simulations
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't') && app.ShowProducts
                    plot(app.UIAxes3, app.ManualSimResults.t, app.ManualSimResults.C_NaOAc, ':', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Manual NaOAc');
                    plot(app.UIAxes3, app.ManualSimResults.t, app.ManualSimResults.C_EtOH, ':', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Manual EtOH');
                     % Add convergence markers for manual products using product thresholds
                    if ~isempty(app.ManualConvergenceData)
                        app.addProductConvergenceMarkers(app.UIAxes3, app.ManualSimResults, app.ManualConvergenceData, 'Manual');
                    end
                end
                
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't') && app.ShowProducts
                    plot(app.UIAxes3, app.SimResults.t, app.SimResults.C_NaOAc, '-', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Optimized NaOAc');
                    plot(app.UIAxes3, app.SimResults.t, app.SimResults.C_EtOH, '-', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Optimized EtOH');
                    
                    % Add convergence markers for optimized products using product thresholds
                    if ~isempty(app.ConvergenceData)
                        app.addProductConvergenceMarkers(app.UIAxes3, app.SimResults, app.ConvergenceData, 'Optimized');
                    end
                end
                
                % Plot ultra-precision simulation products (if available)
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't') && app.ShowProducts
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision NaOAc');
    plot(app.UIAxes3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
         'LineWidth', 4, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra-Precision EtOH');
    
    % Add convergence markers for ultra-precision products using product thresholds
    if ~isempty(app.UltraConvergenceData)
        app.addProductConvergenceMarkers(app.UIAxes3, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra-Precision');
    end
end
                
                % Styling
                title(app.UIAxes3, 'Product Formation with Differential Convergence Analysis', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(app.UIAxes3, 'Time (min)', 'FontWeight', 'bold');
                ylabel(app.UIAxes3, 'Concentration (mol/L)', 'FontWeight', 'bold');
                legend(app.UIAxes3, 'Location', 'best', 'Box', 'on', 'FontSize', 11);
                grid(app.UIAxes3, 'off');
                app.UIAxes3.GridAlpha = 0.3;
                app.UIAxes3.Box = 'on';
                app.UIAxes3.FontSize = 11;
                
                hold(app.UIAxes3, 'off');
                
                % Update info panel with convergence information
                app.updateProductsInfo();
                
            catch
                % Silent fail
            end
        end
        
        % Enhanced products info with separate convergence details
        function updateProductsInfo(app)
            
            try
                infoHTML = '<div style="font-family: -apple-system; padding: 10px; background: #F2F2F7; border-radius: 8px;">';
                infoHTML = [infoHTML '<h3 style="margin: 0 0 10px 0; color: #000;">Product Yields & Differential Convergence</h3>'];
                
                % Manual products analysis
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 'C_NaOAc')
                    final_NaOAc_manual = app.ManualSimResults.C_NaOAc(end);
                    final_EtOH_manual = app.ManualSimResults.C_EtOH(end);
                    conversion_manual = (app.Cin - app.ManualSimResults.C_NaOH(end)) / app.Cin * 100;
                    infoHTML = [infoHTML '<strong>Manual Products:</strong><br>'];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">NaOAc: %.6f mol/L</p>', final_NaOAc_manual)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">EtOH: %.6f mol/L</p>', final_EtOH_manual)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Conversion: %.3f%%</p>', conversion_manual)];
                end
                
                % Optimized products analysis
                if ~isempty(app.SimResults) && isfield(app.SimResults, 'C_NaOAc')
                    final_NaOAc = app.SimResults.C_NaOAc(end);
                    final_EtOH = app.SimResults.C_EtOH(end);
                    conversion = (app.Cin - app.SimResults.C_NaOH(end)) / app.Cin * 100;
                    
                    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                    infoHTML = [infoHTML '<strong>Optimized Products:</strong><br>'];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">NaOAc: %.6f mol/L</p>', final_NaOAc)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">EtOH: %.6f mol/L</p>', final_EtOH)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Conversion: %.3f%%</p>', conversion)];
                end
                
                % Quick sim products analysis
                if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 'C_NaOAc')
                    final_NaOAc_quick = app.QuickSimResults.C_NaOAc(end);
                    final_EtOH_quick = app.QuickSimResults.C_EtOH(end);
                    conversion_quick = (app.Cin - app.QuickSimResults.C_NaOH(end)) / app.Cin * 100;
                    
                    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
                    infoHTML = [infoHTML '<strong>Quick Sim Products:</strong><br>'];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">NaOAc: %.6f mol/L</p>', final_NaOAc_quick)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">EtOH: %.6f mol/L</p>', final_EtOH_quick)];
                    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Conversion: %.3f%%</p>', conversion_quick)];
                end
                % Add this section after the Quick sim products analysis:

% Ultra-precision products analysis
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 'C_NaOAc')
    final_NaOAc_ultra = app.UltraPrecisionResults.C_NaOAc(end);
    final_EtOH_ultra = app.UltraPrecisionResults.C_EtOH(end);
    conversion_ultra = (app.Cin - app.UltraPrecisionResults.C_NaOH(end)) / app.Cin * 100;
    
    infoHTML = [infoHTML '<hr style="margin: 10px 0; border: none; border-top: 1px solid #ccc;">'];
    infoHTML = [infoHTML '<strong>Ultra-Precision Products:</strong><br>'];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">NaOAc: %.6f mol/L</p>', final_NaOAc_ultra)];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">EtOH: %.6f mol/L</p>', final_EtOH_ultra)];
    infoHTML = [infoHTML sprintf('<p style="margin: 5px 0;">Conversion: %.3f%%</p>', conversion_ultra)];
end
                
                infoHTML = [infoHTML '</div>'];
                app.ProductsInfo.HTMLSource = infoHTML;
                
            catch
                % Silent fail
            end
        end
        
        
        % Plot conversion analysis
    % Enhanced Plot conversion analysis with steady-state markers
function plotConversions(app)
    try
        % Clear both conversion plots
        cla(app.UIAxesConversionProduct);
        cla(app.UIAxesConversionReactant);
        
        % Plot conversion of reactant to product (top graph)
        hold(app.UIAxesConversionProduct, 'on');
        
        % Manual simulation conversion (if available)
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            % Calculate product conversion (NaOAc and EtOH formation)
            conv_NaOAc_manual = (app.ManualSimResults.C_NaOAc / app.Cin) * 100;
            conv_EtOH_manual = (app.ManualSimResults.C_EtOH / app.Cin) * 100;
            
            plot(app.UIAxesConversionProduct, app.ManualSimResults.t, conv_NaOAc_manual, ':', ...
                 'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Manual NaOAc Conversion');
            plot(app.UIAxesConversionProduct, app.ManualSimResults.t, conv_EtOH_manual, ':', ...
                 'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Manual EtOH Conversion');
            
            % Add steady-state markers for manual products
            [ss_time_NaOAc, ss_conv_NaOAc] = detectSteadyState(app.ManualSimResults.t, conv_NaOAc_manual);
            [ss_time_EtOH, ss_conv_EtOH] = detectSteadyState(app.ManualSimResults.t, conv_EtOH_manual);
            
            if ~isnan(ss_time_NaOAc)
                scatter(app.UIAxesConversionProduct, ss_time_NaOAc, ss_conv_NaOAc, 120, ...
                       'filled', 'Marker', 's', 'MarkerFaceColor', [255 204 0]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Manual NaOAc SS (%.1f min, %.1f%%)', ss_time_NaOAc, ss_conv_NaOAc));
                
                % Add annotation
                text(app.UIAxesConversionProduct, ss_time_NaOAc + 1, ss_conv_NaOAc + 2, ...
                     sprintf('Manual NaOAc SS\n(%.1f, %.1f%%)', ss_time_NaOAc, ss_conv_NaOAc), ...
                     'FontSize', 8, 'Color', [255 204 0]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [255 204 0]/255, 'Margin', 2);
            end
            
            if ~isnan(ss_time_EtOH)
                scatter(app.UIAxesConversionProduct, ss_time_EtOH, ss_conv_EtOH, 120, ...
                       'filled', 'Marker', 's', 'MarkerFaceColor', [88 86 214]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Manual EtOH SS (%.1f min, %.1f%%)', ss_time_EtOH, ss_conv_EtOH));
                
                % Add annotation
                text(app.UIAxesConversionProduct, ss_time_EtOH + 1, ss_conv_EtOH - 3, ...
                     sprintf('Manual EtOH SS\n(%.1f, %.1f%%)', ss_time_EtOH, ss_conv_EtOH), ...
                     'FontSize', 8, 'Color', [88 86 214]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [88 86 214]/255, 'Margin', 2);
            end
        end
        
        % Optimized simulation conversion (if available)
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            % Calculate product conversion (NaOAc and EtOH formation)
            conv_NaOAc_opt = (app.SimResults.C_NaOAc / app.Cin) * 100;
            conv_EtOH_opt = (app.SimResults.C_EtOH / app.Cin) * 100;
            
            plot(app.UIAxesConversionProduct, app.SimResults.t, conv_NaOAc_opt, '-', ...
                 'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Optimized NaOAc Conversion');
            plot(app.UIAxesConversionProduct, app.SimResults.t, conv_EtOH_opt, '-', ...
                 'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Optimized EtOH Conversion');
             % Add steady-state markers for optimized products - ADD THIS COMPLETE SECTION
    [ss_time_NaOAc_opt, ss_conv_NaOAc_opt] = detectSteadyState(app.SimResults.t, conv_NaOAc_opt);
    [ss_time_EtOH_opt, ss_conv_EtOH_opt] = detectSteadyState(app.SimResults.t, conv_EtOH_opt);
    
    if ~isnan(ss_time_NaOAc_opt)
        scatter(app.UIAxesConversionProduct, ss_time_NaOAc_opt, ss_conv_NaOAc_opt, 120, ...
               'filled', 'Marker', 'o', 'MarkerFaceColor', [255 204 0]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
               'DisplayName', sprintf('Opt NaOAc SS (%.1f min, %.1f%%)', ss_time_NaOAc_opt, ss_conv_NaOAc_opt));
        
        % Add vertical and horizontal lines
        xline(app.UIAxesConversionProduct, ss_time_NaOAc_opt, '--', 'Color', [255 204 0]/255, ...
              'Alpha', 0.6, 'LineWidth', 1.5);
        yline(app.UIAxesConversionProduct, ss_conv_NaOAc_opt, '--', 'Color', [255 204 0]/255, ...
              'Alpha', 0.6, 'LineWidth', 1.5);
        
        % Add annotation
        text(app.UIAxesConversionProduct, ss_time_NaOAc_opt + 1, ss_conv_NaOAc_opt + 5, ...
             sprintf('Opt NaOAc SS\n(%.1f, %.1f%%)', ss_time_NaOAc_opt, ss_conv_NaOAc_opt), ...
             'FontSize', 8, 'Color', [255 204 0]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [255 204 0]/255, 'Margin', 2);
    end
    
    if ~isnan(ss_time_EtOH_opt)
        scatter(app.UIAxesConversionProduct, ss_time_EtOH_opt, ss_conv_EtOH_opt, 120, ...
               'filled', 'Marker', 'o', 'MarkerFaceColor', [88 86 214]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
               'DisplayName', sprintf('Opt EtOH SS (%.1f min, %.1f%%)', ss_time_EtOH_opt, ss_conv_EtOH_opt));
        
        % Add vertical and horizontal lines
        xline(app.UIAxesConversionProduct, ss_time_EtOH_opt, '--', 'Color', [88 86 214]/255, ...
              'Alpha', 0.6, 'LineWidth', 1.5);
        yline(app.UIAxesConversionProduct, ss_conv_EtOH_opt, '--', 'Color', [88 86 214]/255, ...
              'Alpha', 0.6, 'LineWidth', 1.5);
        
        % Add annotation
        text(app.UIAxesConversionProduct, ss_time_EtOH_opt + 1, ss_conv_EtOH_opt - 3, ...
             sprintf('Opt EtOH SS\n(%.1f, %.1f%%)', ss_time_EtOH_opt, ss_conv_EtOH_opt), ...
             'FontSize', 8, 'Color', [88 86 214]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [88 86 214]/255, 'Margin', 2);
    end
end
        
        % Quick simulation conversion (if available)
        if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
            % Calculate product conversion (NaOAc and EtOH formation)
            conv_NaOAc_quick = (app.QuickSimResults.C_NaOAc / app.Cin) * 100;
            conv_EtOH_quick = (app.QuickSimResults.C_EtOH / app.Cin) * 100;
            
            plot(app.UIAxesConversionProduct, app.QuickSimResults.t, conv_NaOAc_quick, '--', ...
                 'LineWidth', 2.5, 'Color', [255 165 0]/255, 'DisplayName', 'Quick Sim NaOAc Conversion');
            plot(app.UIAxesConversionProduct, app.QuickSimResults.t, conv_EtOH_quick, '--', ...
                 'LineWidth', 2.5, 'Color', [30 144 255]/255, 'DisplayName', 'Quick Sim EtOH Conversion');
            
            % Add steady-state markers for quick sim products
            [ss_time_NaOAc_quick, ss_conv_NaOAc_quick] = detectSteadyState(app.QuickSimResults.t, conv_NaOAc_quick);
            [ss_time_EtOH_quick, ss_conv_EtOH_quick] = detectSteadyState(app.QuickSimResults.t, conv_EtOH_quick);
            
            if ~isnan(ss_time_NaOAc_quick)
                scatter(app.UIAxesConversionProduct, ss_time_NaOAc_quick, ss_conv_NaOAc_quick, 100, ...
                       'filled', 'Marker', 'd', 'MarkerFaceColor', [255 165 0]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Quick NaOAc SS (%.1f min, %.1f%%)', ss_time_NaOAc_quick, ss_conv_NaOAc_quick));
            end
            
            if ~isnan(ss_time_EtOH_quick)
                scatter(app.UIAxesConversionProduct, ss_time_EtOH_quick, ss_conv_EtOH_quick, 100, ...
                       'filled', 'Marker', 'd', 'MarkerFaceColor', [30 144 255]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Quick EtOH SS (%.1f min, %.1f%%)', ss_time_EtOH_quick, ss_conv_EtOH_quick));
            end
        end
        
        % Ultra-precision simulation conversion (if available) - ADD THIS SECTION
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
    % Calculate product conversion (NaOAc and EtOH formation)
    conv_NaOAc_ultra = (app.UltraPrecisionResults.C_NaOAc / app.Cin) * 100;
    conv_EtOH_ultra = (app.UltraPrecisionResults.C_EtOH / app.Cin) * 100;
    
    plot(app.UIAxesConversionProduct, app.UltraPrecisionResults.t, conv_NaOAc_ultra, '-', ...
         'LineWidth', 5, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision NaOAc Conversion');
    plot(app.UIAxesConversionProduct, app.UltraPrecisionResults.t, conv_EtOH_ultra, '-', ...
         'LineWidth', 5, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra-Precision EtOH Conversion');
    
    % Add steady-state markers for ultra-precision products
    [ss_time_NaOAc_ultra, ss_conv_NaOAc_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_NaOAc_ultra);
    [ss_time_EtOH_ultra, ss_conv_EtOH_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_EtOH_ultra);
    
    if ~isnan(ss_time_NaOAc_ultra)
        scatter(app.UIAxesConversionProduct, ss_time_NaOAc_ultra, ss_conv_NaOAc_ultra, 180, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [255 20 147]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 3, ...
               'DisplayName', sprintf('Ultra NaOAc SS (%.1f min, %.1f%%)', ss_time_NaOAc_ultra, ss_conv_NaOAc_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionProduct, ss_time_NaOAc_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        yline(app.UIAxesConversionProduct, ss_conv_NaOAc_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        
        % Add enhanced annotation for ultra-precision
        text(app.UIAxesConversionProduct, ss_time_NaOAc_ultra + 1, ss_conv_NaOAc_ultra + 8, ...
             sprintf('🏆 Ultra NaOAc SS\n(%.1f, %.1f%%)', ss_time_NaOAc_ultra, ss_conv_NaOAc_ultra), ...
             'FontSize', 10, 'Color', [255 20 147]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [255 20 147]/255, 'Margin', 4);
    end
    
    if ~isnan(ss_time_EtOH_ultra)
        scatter(app.UIAxesConversionProduct, ss_time_EtOH_ultra, ss_conv_EtOH_ultra, 180, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [139 0 139]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 3, ...
               'DisplayName', sprintf('Ultra EtOH SS (%.1f min, %.1f%%)', ss_time_EtOH_ultra, ss_conv_EtOH_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionProduct, ss_time_EtOH_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        yline(app.UIAxesConversionProduct, ss_conv_EtOH_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        
        % Add enhanced annotation for ultra-precision
        text(app.UIAxesConversionProduct, ss_time_EtOH_ultra + 1, ss_conv_EtOH_ultra - 6, ...
             sprintf('🏆 Ultra EtOH SS\n(%.1f, %.1f%%)', ss_time_EtOH_ultra, ss_conv_EtOH_ultra), ...
             'FontSize', 10, 'Color', [139 0 139]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [139 0 139]/255, 'Margin', 4);
    end
end
        
        title(app.UIAxesConversionProduct, 'Conversion: Reactant to Product Formation with Steady-State Markers', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel(app.UIAxesConversionProduct, 'Time (min)', 'FontWeight', 'bold');
        ylabel(app.UIAxesConversionProduct, 'Conversion (%)', 'FontWeight', 'bold');
        legend(app.UIAxesConversionProduct, 'Location', 'best', 'Box', 'on', 'FontSize', 10);
        grid(app.UIAxesConversionProduct, 'off');
        app.UIAxesConversionProduct.GridAlpha = 0.3;
        app.UIAxesConversionProduct.Box = 'on';
        app.UIAxesConversionProduct.FontSize = 11;
        ylim(app.UIAxesConversionProduct, [0 100]);
        hold(app.UIAxesConversionProduct, 'off');
        
        % Plot unreacted reactant conversion (bottom graph)
        hold(app.UIAxesConversionReactant, 'on');
        
        % Manual simulation unreacted conversion (if available)
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            % Calculate unreacted conversion for both reactants
            conv_NaOH_manual = ((app.Cin - app.ManualSimResults.C_NaOH) / app.Cin) * 100;
            conv_EtOAc_manual = ((app.Cin - app.ManualSimResults.C_EtOAc) / app.Cin) * 100;
            
            plot(app.UIAxesConversionReactant, app.ManualSimResults.t, conv_NaOH_manual, ':', ...
                 'LineWidth', 3, 'Color', [255 59 48]/255, 'DisplayName', 'Manual NaOH Conversion');
            plot(app.UIAxesConversionReactant, app.ManualSimResults.t, conv_EtOAc_manual, ':', ...
                 'LineWidth', 3, 'Color', [255 149 0]/255, 'DisplayName', 'Manual EtOAc Conversion');
            
            % Add steady-state markers for manual reactants
            [ss_time_NaOH_manual, ss_conv_NaOH_manual] = detectSteadyState(app.ManualSimResults.t, conv_NaOH_manual);
            [ss_time_EtOAc_manual, ss_conv_EtOAc_manual] = detectSteadyState(app.ManualSimResults.t, conv_EtOAc_manual);
            
            if ~isnan(ss_time_NaOH_manual)
                scatter(app.UIAxesConversionReactant, ss_time_NaOH_manual, ss_conv_NaOH_manual, 120, ...
                       'filled', 'Marker', 's', 'MarkerFaceColor', [255 59 48]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Manual NaOH SS (%.1f min, %.1f%%)', ss_time_NaOH_manual, ss_conv_NaOH_manual));
                
                % Add annotation
                text(app.UIAxesConversionReactant, ss_time_NaOH_manual + 1, ss_conv_NaOH_manual + 2, ...
                     sprintf('Manual NaOH SS\n(%.1f, %.1f%%)', ss_time_NaOH_manual, ss_conv_NaOH_manual), ...
                     'FontSize', 8, 'Color', [255 59 48]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [255 59 48]/255, 'Margin', 2);
            end
            
            if ~isnan(ss_time_EtOAc_manual)
                scatter(app.UIAxesConversionReactant, ss_time_EtOAc_manual, ss_conv_EtOAc_manual, 120, ...
                       'filled', 'Marker', 's', 'MarkerFaceColor', [255 149 0]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Manual EtOAc SS (%.1f min, %.1f%%)', ss_time_EtOAc_manual, ss_conv_EtOAc_manual));
                
                % Add annotation
                text(app.UIAxesConversionReactant, ss_time_EtOAc_manual + 1, ss_conv_EtOAc_manual - 3, ...
                     sprintf('Manual EtOAc SS\n(%.1f, %.1f%%)', ss_time_EtOAc_manual, ss_conv_EtOAc_manual), ...
                     'FontSize', 8, 'Color', [255 149 0]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [255 149 0]/255, 'Margin', 2);
            end
        end
        
        % Optimized simulation unreacted conversion (if available)
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            % Calculate unreacted conversion for both reactants
            conv_NaOH_opt = ((app.Cin - app.SimResults.C_NaOH) / app.Cin) * 100;
            conv_EtOAc_opt = ((app.Cin - app.SimResults.C_EtOAc) / app.Cin) * 100;
            
            plot(app.UIAxesConversionReactant, app.SimResults.t, conv_NaOH_opt, '-', ...
                 'LineWidth', 3, 'Color', [255 59 48]/255, 'DisplayName', 'Optimized NaOH Conversion');
            plot(app.UIAxesConversionReactant, app.SimResults.t, conv_EtOAc_opt, '-', ...
                 'LineWidth', 3, 'Color', [255 149 0]/255, 'DisplayName', 'Optimized EtOAc Conversion');
            
            % Add steady-state markers for optimized reactants
            [ss_time_NaOH_opt, ss_conv_NaOH_opt] = detectSteadyState(app.SimResults.t, conv_NaOH_opt);
            [ss_time_EtOAc_opt, ss_conv_EtOAc_opt] = detectSteadyState(app.SimResults.t, conv_EtOAc_opt);
            
            if ~isnan(ss_time_NaOH_opt)
                scatter(app.UIAxesConversionReactant, ss_time_NaOH_opt, ss_conv_NaOH_opt, 120, ...
                       'filled', 'Marker', 'o', 'MarkerFaceColor', [255 59 48]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Opt NaOH SS (%.1f min, %.1f%%)', ss_time_NaOH_opt, ss_conv_NaOH_opt));
                
                % Add vertical and horizontal lines
                xline(app.UIAxesConversionReactant, ss_time_NaOH_opt, '--', 'Color', [255 59 48]/255, ...
                      'Alpha', 0.6, 'LineWidth', 1.5);
                yline(app.UIAxesConversionReactant, ss_conv_NaOH_opt, '--', 'Color', [255 59 48]/255, ...
                      'Alpha', 0.6, 'LineWidth', 1.5);
                
                % Add annotation
                text(app.UIAxesConversionReactant, ss_time_NaOH_opt + 1, ss_conv_NaOH_opt + 5, ...
                     sprintf('Opt NaOH SS\n(%.1f, %.1f%%)', ss_time_NaOH_opt, ss_conv_NaOH_opt), ...
                     'FontSize', 8, 'Color', [255 59 48]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [255 59 48]/255, 'Margin', 2);
            end
            
            if ~isnan(ss_time_EtOAc_opt)
                scatter(app.UIAxesConversionReactant, ss_time_EtOAc_opt, ss_conv_EtOAc_opt, 120, ...
                       'filled', 'Marker', 'o', 'MarkerFaceColor', [255 149 0]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Opt EtOAc SS (%.1f min, %.1f%%)', ss_time_EtOAc_opt, ss_conv_EtOAc_opt));
                
                % Add vertical and horizontal lines
                xline(app.UIAxesConversionReactant, ss_time_EtOAc_opt, '--', 'Color', [255 149 0]/255, ...
                      'Alpha', 0.6, 'LineWidth', 1.5);
                yline(app.UIAxesConversionReactant, ss_conv_EtOAc_opt, '--', 'Color', [255 149 0]/255, ...
                      'Alpha', 0.6, 'LineWidth', 1.5);
                
                % Add annotation
                text(app.UIAxesConversionReactant, ss_time_EtOAc_opt + 1, ss_conv_EtOAc_opt - 3, ...
                     sprintf('Opt EtOAc SS\n(%.1f, %.1f%%)', ss_time_EtOAc_opt, ss_conv_EtOAc_opt), ...
                     'FontSize', 8, 'Color', [255 149 0]/255, 'FontWeight', 'bold', ...
                     'BackgroundColor', 'white', 'EdgeColor', [255 149 0]/255, 'Margin', 2);
            end
        end
        
        % Quick simulation unreacted conversion (if available)
        if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
            % Calculate unreacted conversion for both reactants
            conv_NaOH_quick = ((app.Cin - app.QuickSimResults.C_NaOH) / app.Cin) * 100;
            conv_EtOAc_quick = ((app.Cin - app.QuickSimResults.C_EtOAc) / app.Cin) * 100;
            
            plot(app.UIAxesConversionReactant, app.QuickSimResults.t, conv_NaOH_quick, '--', ...
                 'LineWidth', 2.5, 'Color', [255 80 80]/255, 'DisplayName', 'Quick Sim NaOH Conversion');
            plot(app.UIAxesConversionReactant, app.QuickSimResults.t, conv_EtOAc_quick, '--', ...
                 'LineWidth', 2.5, 'Color', [255 180 50]/255, 'DisplayName', 'Quick Sim EtOAc Conversion');
            
            % Add steady-state markers for quick sim reactants
            [ss_time_NaOH_quick, ss_conv_NaOH_quick] = detectSteadyState(app.QuickSimResults.t, conv_NaOH_quick);
            [ss_time_EtOAc_quick, ss_conv_EtOAc_quick] = detectSteadyState(app.QuickSimResults.t, conv_EtOAc_quick);
            
            if ~isnan(ss_time_NaOH_quick)
                scatter(app.UIAxesConversionReactant, ss_time_NaOH_quick, ss_conv_NaOH_quick, 100, ...
                       'filled', 'Marker', 'd', 'MarkerFaceColor', [255 80 80]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Quick NaOH SS (%.1f min, %.1f%%)', ss_time_NaOH_quick, ss_conv_NaOH_quick));
            end
            
            if ~isnan(ss_time_EtOAc_quick)
                scatter(app.UIAxesConversionReactant, ss_time_EtOAc_quick, ss_conv_EtOAc_quick, 100, ...
                       'filled', 'Marker', 'd', 'MarkerFaceColor', [255 180 50]/255, ...
                       'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
                       'DisplayName', sprintf('Quick EtOAc SS (%.1f min, %.1f%%)', ss_time_EtOAc_quick, ss_conv_EtOAc_quick));
            end
        end
% Ultra-precision simulation unreacted conversion (if available) - ADD THIS COMPLETE SECTION
if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
    % Calculate unreacted conversion for both reactants
    conv_NaOH_ultra = ((app.Cin - app.UltraPrecisionResults.C_NaOH) / app.Cin) * 100;
    conv_EtOAc_ultra = ((app.Cin - app.UltraPrecisionResults.C_EtOAc) / app.Cin) * 100;
    
    plot(app.UIAxesConversionReactant, app.UltraPrecisionResults.t, conv_NaOH_ultra, '-', ...
         'LineWidth', 5, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision NaOH Conversion');
    plot(app.UIAxesConversionReactant, app.UltraPrecisionResults.t, conv_EtOAc_ultra, '-', ...
         'LineWidth', 5, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra-Precision EtOAc Conversion');
    
    % Add steady-state markers for ultra-precision reactants
    [ss_time_NaOH_ultra, ss_conv_NaOH_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_NaOH_ultra);
    [ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_EtOAc_ultra);
    
    if ~isnan(ss_time_NaOH_ultra)
        scatter(app.UIAxesConversionReactant, ss_time_NaOH_ultra, ss_conv_NaOH_ultra, 180, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [255 20 147]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 3, ...
               'DisplayName', sprintf('Ultra NaOH SS (%.1f min, %.1f%%)', ss_time_NaOH_ultra, ss_conv_NaOH_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionReactant, ss_time_NaOH_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        yline(app.UIAxesConversionReactant, ss_conv_NaOH_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        
        % Add enhanced annotation for ultra-precision
        text(app.UIAxesConversionReactant, ss_time_NaOH_ultra + 1, ss_conv_NaOH_ultra + 8, ...
             sprintf('🏆 Ultra NaOH SS\n(%.1f, %.1f%%)', ss_time_NaOH_ultra, ss_conv_NaOH_ultra), ...
             'FontSize', 10, 'Color', [255 20 147]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [255 20 147]/255, 'Margin', 4);
    end
    
    if ~isnan(ss_time_EtOAc_ultra)
        scatter(app.UIAxesConversionReactant, ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra, 180, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [139 0 139]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 3, ...
               'DisplayName', sprintf('Ultra EtOAc SS (%.1f min, %.1f%%)', ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionReactant, ss_time_EtOAc_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        yline(app.UIAxesConversionReactant, ss_conv_EtOAc_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.95, 'LineWidth', 3.5);
        
        % Add enhanced annotation for ultra-precision
        text(app.UIAxesConversionReactant, ss_time_EtOAc_ultra + 1, ss_conv_EtOAc_ultra - 6, ...
             sprintf('🏆 Ultra EtOAc SS\n(%.1f, %.1f%%)', ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra), ...
             'FontSize', 10, 'Color', [139 0 139]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [139 0 139]/255, 'Margin', 4);
    end
     
    % Add steady-state markers for ultra-precision reactants
    [ss_time_NaOH_ultra, ss_conv_NaOH_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_NaOH_ultra);
    [ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra] = detectSteadyState(app.UltraPrecisionResults.t, conv_EtOAc_ultra);
    
    if ~isnan(ss_time_NaOH_ultra)
        scatter(app.UIAxesConversionReactant, ss_time_NaOH_ultra, ss_conv_NaOH_ultra, 140, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [255 20 147]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
               'DisplayName', sprintf('Ultra NaOH SS (%.1f min, %.1f%%)', ss_time_NaOH_ultra, ss_conv_NaOH_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionReactant, ss_time_NaOH_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.8, 'LineWidth', 2.5);
        yline(app.UIAxesConversionReactant, ss_conv_NaOH_ultra, '-', 'Color', [255 20 147]/255, ...
              'Alpha', 0.8, 'LineWidth', 2.5);
        
        % Add annotation
        text(app.UIAxesConversionReactant, ss_time_NaOH_ultra + 1, ss_conv_NaOH_ultra + 6, ...
             sprintf('Ultra NaOH SS\n(%.1f, %.1f%%)', ss_time_NaOH_ultra, ss_conv_NaOH_ultra), ...
             'FontSize', 8, 'Color', [255 20 147]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [255 20 147]/255, 'Margin', 2);
    end
    
    if ~isnan(ss_time_EtOAc_ultra)
        scatter(app.UIAxesConversionReactant, ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra, 140, ...
               'filled', 'Marker', 'p', 'MarkerFaceColor', [139 0 139]/255, ...
               'MarkerEdgeColor', 'black', 'LineWidth', 2, ...
               'DisplayName', sprintf('Ultra EtOAc SS (%.1f min, %.1f%%)', ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra));
        
        % Add enhanced markers for ultra-precision
        xline(app.UIAxesConversionReactant, ss_time_EtOAc_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.8, 'LineWidth', 2.5);
        yline(app.UIAxesConversionReactant, ss_conv_EtOAc_ultra, '-', 'Color', [139 0 139]/255, ...
              'Alpha', 0.8, 'LineWidth', 2.5);
        
        % Add annotation
        text(app.UIAxesConversionReactant, ss_time_EtOAc_ultra + 1, ss_conv_EtOAc_ultra - 4, ...
             sprintf('Ultra EtOAc SS\n(%.1f, %.1f%%)', ss_time_EtOAc_ultra, ss_conv_EtOAc_ultra), ...
             'FontSize', 8, 'Color', [139 0 139]/255, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', [139 0 139]/255, 'Margin', 2);
    end
end
        title(app.UIAxesConversionReactant, 'Conversion: Unreacted Reactants with Steady-State Markers', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel(app.UIAxesConversionReactant, 'Time (min)', 'FontWeight', 'bold');
        ylabel(app.UIAxesConversionReactant, 'Conversion (%)', 'FontWeight', 'bold');
        legend(app.UIAxesConversionReactant, 'Location', 'best', 'Box', 'on', 'FontSize', 10);
        grid(app.UIAxesConversionReactant, 'off');
        app.UIAxesConversionReactant.GridAlpha = 0.3;
        app.UIAxesConversionReactant.Box = 'on';
        app.UIAxesConversionReactant.FontSize = 11;
        ylim(app.UIAxesConversionReactant, [0 100]);
        hold(app.UIAxesConversionReactant, 'off');
        
    catch ME
        app.updateStatus(['Conversion plot error: ' ME.message], 'error');
    end
    
    % Nested function for steady-state detection
    function [ss_time, ss_value] = detectSteadyState(time_data, conversion_data)
        try
            ss_time = NaN;
            ss_value = NaN;
            
            if length(time_data) < 10
                return;
            end
            
            % Steady-state threshold: relative change < 0.1% for conversion
            ss_threshold = 0.00005; % 0.1% relative change threshold
            
            % Check from 10th point onwards
            for i = 10:length(conversion_data)
                if i > 1
                    rel_change = abs(conversion_data(i) - conversion_data(i-1)) / ...
                                 max(eps, abs(conversion_data(i-1))) * 100;
                    
                    if rel_change < ss_threshold
                        ss_time = time_data(i);
                        ss_value = conversion_data(i);
                        break;
                    end
                end
            end
            
        catch
            ss_time = NaN;
            ss_value = NaN;
        end
    end
    % Add this line at the very end of the plotConversions function, just before the final 'end':
% Update conversion info panel
app.updateConversionInfo();
end
        function updateAnalysisDisplay(app)
    try
        % Enhanced convergence analysis text with relative change criteria
       info = sprintf('ENHANCED CONVERGENCE ANALYSIS - MULTI-LEVEL CRITERIA\n\n');
        
        % Display convergence thresholds first
        info = sprintf('%sCONVERGENCE THRESHOLDS (Relative Change %%):\n', info);
        info = sprintf('%sManual & Optimized Reactants: %.1f%% (EQUAL)\n', info, app.ManualReactantThreshold);
        info = sprintf('%sManual & Optimized Products: %.1f%% (EQUAL)\n', info, app.ManualProductThreshold);
        info = sprintf('%sManual & Optimized Temperature: %.1f%% (EQUAL)\n', info, app.ManualTemperatureThreshold);
        info = sprintf('%sQuick Sim Reactants: %.1f%% (relaxed)\n', info, app.QuickSimReactantThreshold);
        info = sprintf('%sQuick Sim Products: %.1f%% (relaxed)\n', info, app.QuickSimProductThreshold);
        info = sprintf('%sQuick Sim Temperature: %.1f%% (relaxed)\n', info, app.QuickSimTemperatureThreshold);
        
        % NEW: Ultra-precision thresholds
if ~isempty(app.UltraPrecisionResults)
    info = sprintf('%s🏆 Ultra-Precision Reactants: %.4f%% (STRICT)\n', info, app.UltraReactantThreshold*100);
    info = sprintf('%s🏆 Ultra-Precision Products: %.5f%% (STRICT)\n', info, app.UltraProductThreshold*100);
    info = sprintf('%s🏆 Ultra-Precision Temperature: %.4f%% (STRICT)\n', info, app.UltraTemperatureThreshold*100);
end

info = sprintf('%sQuick Sim Reactants: %.3f%% (Medium)\n', info, app.QuickSimReactantThreshold*100);
info = sprintf('%sQuick Sim Products: %.3f%% (Medium)\n', info, app.QuickSimProductThreshold*100);

info = sprintf('%s\n=== CONVERGENCE RESULTS ===\n', info);

% Manual convergence results
if ~isempty(app.ManualSimResults)
    if app.IsManualConverged
        info = sprintf('%sMANUAL (Relaxed): ✓ CONVERGED\n', info);
    else
        info = sprintf('%sMANUAL (Relaxed): ⚠ NOT CONVERGED\n', info);
    end
    
    if ~isempty(app.ManualConvergenceData)
        if isfield(app.ManualConvergenceData, 'NaOH_converged') && app.ManualConvergenceData.NaOH_converged
            info = sprintf('%s  - NaOH: %.2f min (%.6f mol/L)\n', info, app.ManualConvergenceData.NaOH_time, app.ManualConvergenceData.NaOH_conc);
        end
        if isfield(app.ManualConvergenceData, 'EtOAc_converged') && app.ManualConvergenceData.EtOAc_converged
            info = sprintf('%s  - EtOAc: %.2f min (%.6f mol/L)\n', info, app.ManualConvergenceData.EtOAc_time, app.ManualConvergenceData.EtOAc_conc);
        end
    end
end

% Optimized convergence results
if ~isempty(app.SimResults)
    if app.IsConverged
        info = sprintf('%sOPTIMIZED (Standard): ✓ CONVERGED\n', info);
    else
        info = sprintf('%sOPTIMIZED (Standard): ⚠ NOT CONVERGED\n', info);
    end
    
    if ~isempty(app.ConvergenceData)
        if isfield(app.ConvergenceData, 'NaOH_converged') && app.ConvergenceData.NaOH_converged
            info = sprintf('%s  - NaOH: %.2f min (%.6f mol/L)\n', info, app.ConvergenceData.NaOH_time, app.ConvergenceData.NaOH_conc);
        end
        if isfield(app.ConvergenceData, 'EtOAc_converged') && app.ConvergenceData.EtOAc_converged
            info = sprintf('%s  - EtOAc: %.2f min (%.6f mol/L)\n', info, app.ConvergenceData.EtOAc_time, app.ConvergenceData.EtOAc_conc);
        end
    end
end

% NEW: Ultra-precision convergence results
if ~isempty(app.UltraPrecisionResults)
    if app.IsUltraConverged
        info = sprintf('%s🏆 ULTRA-PRECISION (Strict): ✓ CONVERGED\n', info);
    else
        info = sprintf('%s🏆 ULTRA-PRECISION (Strict): ⚠ NOT CONVERGED\n', info);
    end
    
    if ~isempty(app.UltraConvergenceData)
        if isfield(app.UltraConvergenceData, 'NaOH_converged') && app.UltraConvergenceData.NaOH_converged
            info = sprintf('%s  - NaOH: %.2f min (%.8f mol/L)\n', info, app.UltraConvergenceData.NaOH_time, app.UltraConvergenceData.NaOH_conc);
        end
        if isfield(app.UltraConvergenceData, 'EtOAc_converged') && app.UltraConvergenceData.EtOAc_converged
            info = sprintf('%s  - EtOAc: %.2f min (%.8f mol/L)\n', info, app.UltraConvergenceData.EtOAc_time, app.UltraConvergenceData.EtOAc_conc);
        end
        if isfield(app.UltraConvergenceData, 'NaOAc_converged') && app.UltraConvergenceData.NaOAc_converged
            info = sprintf('%s  - NaOAc: %.2f min (%.8f mol/L)\n', info, app.UltraConvergenceData.NaOAc_time, app.UltraConvergenceData.NaOAc_conc);
        end
        if isfield(app.UltraConvergenceData, 'EtOH_converged') && app.UltraConvergenceData.EtOH_converged
            info = sprintf('%s  - EtOH: %.2f min (%.8f mol/L)\n', info, app.UltraConvergenceData.EtOH_time, app.UltraConvergenceData.EtOH_conc);
        end
        if isfield(app.UltraConvergenceData, 'T_converged') && app.UltraConvergenceData.T_converged
            info = sprintf('%s  - Temperature: %.2f min (%.4f °C)\n', info, app.UltraConvergenceData.T_time, app.UltraConvergenceData.T_value);
        end
    end
end

% Quick sim convergence results
if ~isempty(app.QuickSimResults)
    if app.IsQuickSimConverged
        info = sprintf('%sQUICK SIM (Medium): ✓ CONVERGED\n', info);
    else
        info = sprintf('%sQUICK SIM (Medium): ⚠ NOT CONVERGED\n', info);
    end
    
    if ~isempty(app.QuickSimConvergenceData)
        if isfield(app.QuickSimConvergenceData, 'NaOH_converged') && app.QuickSimConvergenceData.NaOH_converged
            info = sprintf('%s  - NaOH: %.2f min (%.6f mol/L)\n', info, app.QuickSimConvergenceData.NaOH_time, app.QuickSimConvergenceData.NaOH_conc);
        end
        if isfield(app.QuickSimConvergenceData, 'EtOAc_converged') && app.QuickSimConvergenceData.EtOAc_converged
            info = sprintf('%s  - EtOAc: %.2f min (%.6f mol/L)\n', info, app.QuickSimConvergenceData.EtOAc_time, app.QuickSimConvergenceData.EtOAc_conc);
        end
    end
end

% Performance comparison if multiple methods available
if (~isempty(app.SimResults) || ~isempty(app.ManualSimResults)) && ~isempty(app.UltraPrecisionResults)
    info = sprintf('%s\n=== PERFORMANCE COMPARISON ===\n', info);
    
    if ~isnan(app.RSquaredManual)
        info = sprintf('%sManual R²: %.6f\n', info, app.RSquaredManual);
    end
    if ~isnan(app.RSquared)
        info = sprintf('%sOptimized R²: %.6f\n', info, app.RSquared);
    end
    if ~isnan(app.RSquaredUltra)
        info = sprintf('%s🏆 Ultra-Precision R²: %.8f\n', info, app.RSquaredUltra);
    end
    
    % Error comparison
    if ~isempty(app.ComparisonTableData)
        opt_max_error = max(app.ComparisonTableData.Optimized_Error_percent);
        info = sprintf('%sOptimized Max Error: %.4f%%\n', info, opt_max_error);
    end
    
    if ~isempty(app.ComparisonTableDataUltra)
        ultra_max_error = max(app.ComparisonTableDataUltra.Ultra_Error_percent);
        info = sprintf('%s🏆 Ultra-Precision Max Error: %.6f%%\n', info, ultra_max_error);
    end
end

info = sprintf('%s\n=== RECOMMENDATIONS ===\n', info);
if ~isempty(app.UltraPrecisionResults) && app.IsUltraConverged
    info = sprintf('%s🏆 Ultra-precision achieved strict convergence!\n', info);
    info = sprintf('%s✨ Highest accuracy simulation available.\n', info);
elseif ~isempty(app.SimResults) && app.IsConverged
    info = sprintf('%s✓ Standard optimization achieved convergence.\n', info);
    info = sprintf('%s💡 Try Ultra-Precision for even better accuracy.\n', info);
else
    info = sprintf('%s⚠ Consider extending simulation time or adjusting parameters.\n', info);
    info = sprintf('%s🔄 Try different convergence criteria or methods.\n', info);
end

        
        app.ConvergenceText.Value = info;
        app.updateComparisonTableDisplay();
        
    catch
        app.ConvergenceText.Value = 'Error displaying enhanced analysis';
    end
    app.updateProfileTables();
end
        
        function updateProfileTables(app)
            try
                % Update concentration profile table
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    % Create concentration profile data
                    concData = [app.SimResults.t, app.SimResults.C_NaOH, app.SimResults.C_EtOAc, ...
                               app.SimResults.C_NaOAc, app.SimResults.C_EtOH];
                    app.ConcentrationTable.Data = concData;
                    app.ConcentrationTable.ColumnName = {'Time (min)', 'NaOH', 'EtOAc', 'NaOAc', 'EtOH'};
                end
                
                % Update temperature profile table
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    % Create temperature profile data
                    tempData = [app.SimResults.t, app.SimResults.T_reactor, app.SimResults.T_jacket];
                    app.TemperatureTable.Data = tempData;
                    app.TemperatureTable.ColumnName = {'Time (min)', 'Reactor (°C)', 'Jacket (°C)'};
                end
                
            catch
                % Silent fail
            end
            
        end
        
        % Export methods
       function exportData(app)
    try
        if isempty(app.SimResults) && isempty(app.ManualSimResults) && isempty(app.QuickSimResults)
            app.updateStatus('No data to export!', 'error');
            return;
        end
        
        app.showLoadingDialog('Exporting Data', 'Preparing comprehensive data export...');
        
        [file, path] = uiputfile({'*.xlsx', 'Excel Files (*.xlsx)'}, ...
                                'Export Comprehensive Data', 'CSTR_Comprehensive_Results.xlsx');
        if isequal(file, 0)
            app.closeLoadingDialog();
            return;
        end
        
        fullPath = fullfile(path, file);
        
        % Delete existing file if it exists to avoid conflicts
        if exist(fullPath, 'file')
            delete(fullPath);
        end
        
        %% 1. REACTOR PARAMETERS & SUMMARY
        app.updateLoadingDialog(0.1, 'Exporting reactor parameters...');
        
        % Calculate residence times
        tau_total = app.V / app.V0;  % Total residence time
        tau_individual = app.V / (app.V0/2);  % Individual reactant residence time
        
        % Create parameters summary table
        paramNames = {'Reactor Volume (L)'; 'Total Flowrate (L/min)'; 'EtOAc Flowrate (L/min)'; 
                     'NaOH Flowrate (L/min)'; 'Inlet Concentration (mol/L)'; 'Temperature (°C)';
                     'Simulation Time (min)'; 'Total Residence Time (min)'; 'Individual Residence Time (min)';
                     'Manual Rate Constant (L/mol·min)'; 'Optimized Rate Constant (L/mol·min)';
                     'Manual R-squared'; 'Optimized R-squared'};
        
        paramValues = {app.V; app.V0; app.V0_EtOAc; app.V1_NaOH; app.Cin; app.SelectedTemp;
                      app.SimTime; tau_total; tau_individual;
                      app.k_manual; app.k_fit; app.RSquaredManual; app.RSquared};
        
        parametersTable = table(paramNames, paramValues, 'VariableNames', {'Parameter', 'Value'});
        
        try
            writetable(parametersTable, fullPath, 'Sheet', '1_Parameters', 'WriteMode', 'overwritesheet');
        catch ME
            fprintf('Warning: Failed to write parameters: %s\n', ME.message);
        end
        
        %% 2. EXPERIMENTAL DATA
        app.updateLoadingDialog(0.15, 'Exporting experimental data...');
        
        if ~isempty(app.Data) && height(app.Data) > 0
            try
                writetable(app.Data, fullPath, 'Sheet', '2_Experimental', 'WriteMode', 'overwritesheet');
            catch ME
                fprintf('Warning: Failed to write experimental data: %s\n', ME.message);
            end
        end
        
        %% 3. MANUAL SIMULATION RESULTS
        app.updateLoadingDialog(0.25, 'Exporting manual simulation results...');
        
        if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
            try
                manualData = table(app.ManualSimResults.t, app.ManualSimResults.C_NaOH, app.ManualSimResults.C_EtOAc, ...
                                  app.ManualSimResults.C_NaOAc, app.ManualSimResults.C_EtOH, app.ManualSimResults.T_reactor, ...
                                  app.ManualSimResults.T_jacket, ...
                                  'VariableNames', {'Time_min', 'C_NaOH_molL', 'C_EtOAc_molL', ...
                                  'C_NaOAc_molL', 'C_EtOH_molL', 'T_reactor_C', 'T_jacket_C'});
                writetable(manualData, fullPath, 'Sheet', '3_Manual', 'WriteMode', 'overwritesheet');
            catch ME
                fprintf('Warning: Failed to write manual simulation data: %s\n', ME.message);
            end
        end
        
        %% 4. OPTIMIZED SIMULATION RESULTS
        app.updateLoadingDialog(0.35, 'Exporting optimized simulation results...');
        
        if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
            try
                optimizedData = table(app.SimResults.t, app.SimResults.C_NaOH, app.SimResults.C_EtOAc, ...
                                     app.SimResults.C_NaOAc, app.SimResults.C_EtOH, app.SimResults.T_reactor, ...
                                     app.SimResults.T_jacket, ...
                                     'VariableNames', {'Time_min', 'C_NaOH_molL', 'C_EtOAc_molL', ...
                                     'C_NaOAc_molL', 'C_EtOH_molL', 'T_reactor_C', 'T_jacket_C'});
                writetable(optimizedData, fullPath, 'Sheet', '4_Optimized', 'WriteMode', 'overwritesheet');
            catch ME
                fprintf('Warning: Failed to write optimized simulation data: %s\n', ME.message);
            end
        end
        
        %% 5. QUICK SIMULATION RESULTS
        app.updateLoadingDialog(0.45, 'Exporting Quick Sim results...');
        
        if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
            try
                quickSimData = table(app.QuickSimResults.t, app.QuickSimResults.C_NaOH, app.QuickSimResults.C_EtOAc, ...
                                    app.QuickSimResults.C_NaOAc, app.QuickSimResults.C_EtOH, app.QuickSimResults.T_reactor, ...
                                    app.QuickSimResults.T_jacket, ...
                                    'VariableNames', {'Time_min', 'C_NaOH_molL', 'C_EtOAc_molL', ...
                                    'C_NaOAc_molL', 'C_EtOH_molL', 'T_reactor_C', 'T_jacket_C'});
                writetable(quickSimData, fullPath, 'Sheet', '5_QuickSim', 'WriteMode', 'overwritesheet');
            catch ME
                fprintf('Warning: Failed to write quick sim data: %s\n', ME.message);
            end
        end
        
        %% 6. ULTRA-PRECISION SIMULATION RESULTS
        app.updateLoadingDialog(0.47, 'Exporting ultra-precision simulation results...');
        
        if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
            try
                ultraPrecisionData = table(app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, app.UltraPrecisionResults.C_EtOAc, ...
                                          app.UltraPrecisionResults.C_NaOAc, app.UltraPrecisionResults.C_EtOH, app.UltraPrecisionResults.T_reactor, ...
                                          app.UltraPrecisionResults.T_jacket, ...
                                          'VariableNames', {'Time_min', 'C_NaOH_molL', 'C_EtOAc_molL', ...
                                          'C_NaOAc_molL', 'C_EtOH_molL', 'T_reactor_C', 'T_jacket_C'});
                writetable(ultraPrecisionData, fullPath, 'Sheet', '6_UltraPrecision', 'WriteMode', 'overwritesheet');
            catch ME
                fprintf('Warning: Failed to write ultra-precision simulation data: %s\n', ME.message);
            end
        end
        
        %% 7. ERROR ANALYSIS & COMPARISON
        app.updateLoadingDialog(0.55, 'Exporting error analysis...');
        
        if ~isempty(app.Data) && height(app.Data) > 0
            try
                % Create comprehensive comparison table
                t_exp = app.Data.Time;
                C_exp = app.Data.Concentration;
                
                % Initialize comparison data
                comparisonData = table(t_exp, C_exp, 'VariableNames', {'Time_min', 'Experimental_molL'});
                
                % Add manual simulation data if available
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                    C_manual_interp = interp1(app.ManualSimResults.t, app.ManualSimResults.C_NaOH, t_exp, 'linear', 'extrap');
                    manual_errors = abs(C_exp - C_manual_interp) ./ C_exp * 100;
                    comparisonData.Manual_Simulated_molL = C_manual_interp;
                    comparisonData.Manual_Error_percent = manual_errors;
                end
                
                % Add optimized simulation data if available
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    C_opt_interp = interp1(app.SimResults.t, app.SimResults.C_NaOH, t_exp, 'linear', 'extrap');
                    opt_errors = abs(C_exp - C_opt_interp) ./ C_exp * 100;
                    comparisonData.Optimized_Simulated_molL = C_opt_interp;
                    comparisonData.Optimized_Error_percent = opt_errors;
                end
                
                % Add ultra-precision simulation data if available
                if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                    C_ultra_interp = interp1(app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, t_exp, 'linear', 'extrap');
                    ultra_errors = abs(C_exp - C_ultra_interp) ./ C_exp * 100;
                    comparisonData.Ultra_Simulated_molL = C_ultra_interp;
                    comparisonData.Ultra_Error_percent = ultra_errors;
                end
                
                % Add quick sim data if available
                if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                    C_quick_interp = interp1(app.QuickSimResults.t, app.QuickSimResults.C_NaOH, t_exp, 'linear', 'extrap');
                    quick_errors = abs(C_exp - C_quick_interp) ./ C_exp * 100;
                    comparisonData.QuickSim_Simulated_molL = C_quick_interp;
                    comparisonData.QuickSim_Error_percent = quick_errors;
                end
                
                writetable(comparisonData, fullPath, 'Sheet', '7_Error_Analysis', 'WriteMode', 'overwritesheet');
                
                % Create error statistics summary
                errorStats = {};
                errorValues = {};
                
                if exist('manual_errors', 'var')
                    errorStats = [errorStats; 'Manual Max Error (%)'; 'Manual Min Error (%)'; 'Manual Mean Error (%)'];
                    errorValues = [errorValues; max(manual_errors); min(manual_errors); mean(manual_errors)];
                end
                
                if exist('opt_errors', 'var')
                    errorStats = [errorStats; 'Optimized Max Error (%)'; 'Optimized Min Error (%)'; 'Optimized Mean Error (%)'];
                    errorValues = [errorValues; max(opt_errors); min(opt_errors); mean(opt_errors)];
                end
                
                if exist('ultra_errors', 'var')
                    errorStats = [errorStats; 'Ultra Max Error (%)'; 'Ultra Min Error (%)'; 'Ultra Mean Error (%)'];
                    errorValues = [errorValues; max(ultra_errors); min(ultra_errors); mean(ultra_errors)];
                end
                
                if exist('quick_errors', 'var')
                    errorStats = [errorStats; 'Quick Sim Max Error (%)'; 'Quick Sim Min Error (%)'; 'Quick Sim Mean Error (%)'];
                    errorValues = [errorValues; max(quick_errors); min(quick_errors); mean(quick_errors)];
                end
                
                if ~isempty(errorStats)
                    errorSummary = table(errorStats, errorValues, 'VariableNames', {'Error_Metric', 'Value'});
                    writetable(errorSummary, fullPath, 'Sheet', '8_Error_Stats', 'WriteMode', 'overwritesheet');
                end
                
            catch ME
                fprintf('Warning: Failed to write error analysis: %s\n', ME.message);
            end
        end
        
        %% 8. CONVERGENCE ANALYSIS DATA
        app.updateLoadingDialog(0.55, 'Exporting convergence analysis...');
        
        try
            convergenceInfo = {};
            convergenceValues = {};
            
            % Manual convergence data
            if ~isempty(app.ManualConvergenceData)
                convergenceInfo = [convergenceInfo; 'Manual Overall Converged'; 'Manual Reactants Converged'; 
                                  'Manual Products Converged'; 'Manual Temperature Converged'];
                convergenceValues = [convergenceValues; app.IsManualConverged; app.ManualReactantsConverged; 
                                    app.ManualProductsConverged; app.ManualTemperatureConverged];
                
                if isfield(app.ManualConvergenceData, 'NaOH_time') && app.ManualConvergenceData.NaOH_converged
                    convergenceInfo = [convergenceInfo; 'Manual NaOH Convergence Time (min)'; 'Manual NaOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ManualConvergenceData.NaOH_time; app.ManualConvergenceData.NaOH_conc];
                end
                
                if isfield(app.ManualConvergenceData, 'EtOAc_time') && app.ManualConvergenceData.EtOAc_converged
                    convergenceInfo = [convergenceInfo; 'Manual EtOAc Convergence Time (min)'; 'Manual EtOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ManualConvergenceData.EtOAc_time; app.ManualConvergenceData.EtOAc_conc];
                end
                
                if isfield(app.ManualConvergenceData, 'NaOAc_time') && app.ManualConvergenceData.NaOAc_converged
                    convergenceInfo = [convergenceInfo; 'Manual NaOAc Convergence Time (min)'; 'Manual NaOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ManualConvergenceData.NaOAc_time; app.ManualConvergenceData.NaOAc_conc];
                end
                
                if isfield(app.ManualConvergenceData, 'EtOH_time') && app.ManualConvergenceData.EtOH_converged
                    convergenceInfo = [convergenceInfo; 'Manual EtOH Convergence Time (min)'; 'Manual EtOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ManualConvergenceData.EtOH_time; app.ManualConvergenceData.EtOH_conc];
                end
                
                if isfield(app.ManualConvergenceData, 'T_time') && app.ManualConvergenceData.T_converged
                    convergenceInfo = [convergenceInfo; 'Manual Temp Convergence Time (min)'; 'Manual Temp Convergence Value (°C)'];
                    convergenceValues = [convergenceValues; app.ManualConvergenceData.T_time; app.ManualConvergenceData.T_value];
                end
            end
            
            % Optimized convergence data
            if ~isempty(app.ConvergenceData)
                convergenceInfo = [convergenceInfo; 'Optimized Overall Converged'; 'Optimized Reactants Converged'; 
                                  'Optimized Products Converged'; 'Optimized Temperature Converged'];
                convergenceValues = [convergenceValues; app.IsConverged; app.OptimizedReactantsConverged; 
                                    app.OptimizedProductsConverged; app.OptimizedTemperatureConverged];
                
                if isfield(app.ConvergenceData, 'NaOH_time') && app.ConvergenceData.NaOH_converged
                    convergenceInfo = [convergenceInfo; 'Optimized NaOH Convergence Time (min)'; 'Optimized NaOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ConvergenceData.NaOH_time; app.ConvergenceData.NaOH_conc];
                end
                
                if isfield(app.ConvergenceData, 'EtOAc_time') && app.ConvergenceData.EtOAc_converged
                    convergenceInfo = [convergenceInfo; 'Optimized EtOAc Convergence Time (min)'; 'Optimized EtOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ConvergenceData.EtOAc_time; app.ConvergenceData.EtOAc_conc];
                end
                
                if isfield(app.ConvergenceData, 'NaOAc_time') && app.ConvergenceData.NaOAc_converged
                    convergenceInfo = [convergenceInfo; 'Optimized NaOAc Convergence Time (min)'; 'Optimized NaOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ConvergenceData.NaOAc_time; app.ConvergenceData.NaOAc_conc];
                end
                
                if isfield(app.ConvergenceData, 'EtOH_time') && app.ConvergenceData.EtOH_converged
                    convergenceInfo = [convergenceInfo; 'Optimized EtOH Convergence Time (min)'; 'Optimized EtOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.ConvergenceData.EtOH_time; app.ConvergenceData.EtOH_conc];
                end
                
                if isfield(app.ConvergenceData, 'T_time') && app.ConvergenceData.T_converged
                    convergenceInfo = [convergenceInfo; 'Optimized Temp Convergence Time (min)'; 'Optimized Temp Convergence Value (°C)'];
                    convergenceValues = [convergenceValues; app.ConvergenceData.T_time; app.ConvergenceData.T_value];
                end
            end
            
            % Ultra-precision convergence data
            if ~isempty(app.UltraConvergenceData)
                convergenceInfo = [convergenceInfo; 'Ultra Overall Converged'; 'Ultra Reactants Converged'; 
                                  'Ultra Products Converged'; 'Ultra Temperature Converged'];
                convergenceValues = [convergenceValues; app.IsUltraConverged; app.UltraReactantsConverged; 
                                    app.UltraProductsConverged; app.UltraTemperatureConverged];
                
                if isfield(app.UltraConvergenceData, 'NaOH_time') && app.UltraConvergenceData.NaOH_converged
                    convergenceInfo = [convergenceInfo; 'Ultra NaOH Convergence Time (min)'; 'Ultra NaOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.UltraConvergenceData.NaOH_time; app.UltraConvergenceData.NaOH_conc];
                end
                
                if isfield(app.UltraConvergenceData, 'EtOAc_time') && app.UltraConvergenceData.EtOAc_converged
                    convergenceInfo = [convergenceInfo; 'Ultra EtOAc Convergence Time (min)'; 'Ultra EtOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.UltraConvergenceData.EtOAc_time; app.UltraConvergenceData.EtOAc_conc];
                end
                
                if isfield(app.UltraConvergenceData, 'NaOAc_time') && app.UltraConvergenceData.NaOAc_converged
                    convergenceInfo = [convergenceInfo; 'Ultra NaOAc Convergence Time (min)'; 'Ultra NaOAc Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.UltraConvergenceData.NaOAc_time; app.UltraConvergenceData.NaOAc_conc];
                end
                
                if isfield(app.UltraConvergenceData, 'EtOH_time') && app.UltraConvergenceData.EtOH_converged
                    convergenceInfo = [convergenceInfo; 'Ultra EtOH Convergence Time (min)'; 'Ultra EtOH Convergence Value (mol/L)'];
                    convergenceValues = [convergenceValues; app.UltraConvergenceData.EtOH_time; app.UltraConvergenceData.EtOH_conc];
                end
                
                if isfield(app.UltraConvergenceData, 'T_time') && app.UltraConvergenceData.T_converged
                    convergenceInfo = [convergenceInfo; 'Ultra Temp Convergence Time (min)'; 'Ultra Temp Convergence Value (°C)'];
                    convergenceValues = [convergenceValues; app.UltraConvergenceData.T_time; app.UltraConvergenceData.T_value];
                end
            end
            
            % Add convergence thresholds for all methods with field existence checks
            if isprop(app, 'ManualReactantThreshold')
                convergenceInfo = [convergenceInfo; 'Manual Reactant Threshold (%)'];
                convergenceValues = [convergenceValues; app.ManualReactantThreshold*100];
            end
            if isprop(app, 'ManualProductThreshold')
                convergenceInfo = [convergenceInfo; 'Manual Product Threshold (%)'];
                convergenceValues = [convergenceValues; app.ManualProductThreshold*100];
            end
            if isprop(app, 'ManualTemperatureThreshold')
                convergenceInfo = [convergenceInfo; 'Manual Temperature Threshold (%)'];
                convergenceValues = [convergenceValues; app.ManualTemperatureThreshold*100];
            end
            if isprop(app, 'OptimizedReactantThreshold')
                convergenceInfo = [convergenceInfo; 'Optimized Reactant Threshold (%)'];
                convergenceValues = [convergenceValues; app.OptimizedReactantThreshold*100];
            end
            if isprop(app, 'OptimizedProductThreshold')
                convergenceInfo = [convergenceInfo; 'Optimized Product Threshold (%)'];
                convergenceValues = [convergenceValues; app.OptimizedProductThreshold*100];
            end
            if isprop(app, 'OptimizedTemperatureThreshold')
                convergenceInfo = [convergenceInfo; 'Optimized Temperature Threshold (%)'];
                convergenceValues = [convergenceValues; app.OptimizedTemperatureThreshold*100];
            end
            
            if ~isempty(convergenceInfo)
                convergenceTable = table(convergenceInfo, convergenceValues, 'VariableNames', {'Convergence_Metric', 'Value'});
                writetable(convergenceTable, fullPath, 'Sheet', '8_Convergence', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write convergence analysis: %s\n', ME.message);
        end
        
        %% 9. YIELD SUMMARY
        app.updateLoadingDialog(0.65, 'Exporting yield analysis...');
        
        try
            yieldSummary = {};
            yieldValues = {};
            
            % Manual yields
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                final_NaOAc_manual = app.ManualSimResults.C_NaOAc(end);
                final_EtOH_manual = app.ManualSimResults.C_EtOH(end);
                final_conversion_manual = (app.Cin - app.ManualSimResults.C_NaOH(end)) / app.Cin * 100;
                final_NaOH_remaining_manual = app.ManualSimResults.C_NaOH(end);
                final_EtOAc_remaining_manual = app.ManualSimResults.C_EtOAc(end);
                
                yieldSummary = [yieldSummary; 'Manual Final NaOAc (mol/L)'; 'Manual Final EtOH (mol/L)'; 'Manual Final Conversion (%)';
                               'Manual Final NaOH Remaining (mol/L)'; 'Manual Final EtOAc Remaining (mol/L)'];
                yieldValues = [yieldValues; final_NaOAc_manual; final_EtOH_manual; final_conversion_manual;
                              final_NaOH_remaining_manual; final_EtOAc_remaining_manual];
            end
            
            % Optimized yields
            if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                final_NaOAc_opt = app.SimResults.C_NaOAc(end);
                final_EtOH_opt = app.SimResults.C_EtOH(end);
                final_conversion_opt = (app.Cin - app.SimResults.C_NaOH(end)) / app.Cin * 100;
                final_NaOH_remaining_opt = app.SimResults.C_NaOH(end);
                final_EtOAc_remaining_opt = app.SimResults.C_EtOAc(end);
                
                yieldSummary = [yieldSummary; 'Optimized Final NaOAc (mol/L)'; 'Optimized Final EtOH (mol/L)'; 'Optimized Final Conversion (%)';
                               'Optimized Final NaOH Remaining (mol/L)'; 'Optimized Final EtOAc Remaining (mol/L)'];
                yieldValues = [yieldValues; final_NaOAc_opt; final_EtOH_opt; final_conversion_opt;
                              final_NaOH_remaining_opt; final_EtOAc_remaining_opt];
            end
            
            % Ultra-precision yields
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                final_NaOAc_ultra = app.UltraPrecisionResults.C_NaOAc(end);
                final_EtOH_ultra = app.UltraPrecisionResults.C_EtOH(end);
                final_conversion_ultra = (app.Cin - app.UltraPrecisionResults.C_NaOH(end)) / app.Cin * 100;
                final_NaOH_remaining_ultra = app.UltraPrecisionResults.C_NaOH(end);
                final_EtOAc_remaining_ultra = app.UltraPrecisionResults.C_EtOAc(end);
                
                yieldSummary = [yieldSummary; 'Ultra Final NaOAc (mol/L)'; 'Ultra Final EtOH (mol/L)'; 'Ultra Final Conversion (%)';
                               'Ultra Final NaOH Remaining (mol/L)'; 'Ultra Final EtOAc Remaining (mol/L)'];
                yieldValues = [yieldValues; final_NaOAc_ultra; final_EtOH_ultra; final_conversion_ultra;
                              final_NaOH_remaining_ultra; final_EtOAc_remaining_ultra];
            end
            
            % Quick sim yields
            if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                final_NaOAc_quick = app.QuickSimResults.C_NaOAc(end);
                final_EtOH_quick = app.QuickSimResults.C_EtOH(end);
                final_conversion_quick = (app.Cin - app.QuickSimResults.C_NaOH(end)) / app.Cin * 100;
                final_NaOH_remaining_quick = app.QuickSimResults.C_NaOH(end);
                final_EtOAc_remaining_quick = app.QuickSimResults.C_EtOAc(end);
                
                yieldSummary = [yieldSummary; 'Quick Sim Final NaOAc (mol/L)'; 'Quick Sim Final EtOH (mol/L)'; 'Quick Sim Final Conversion (%)';
                               'Quick Sim Final NaOH Remaining (mol/L)'; 'Quick Sim Final EtOAc Remaining (mol/L)'];
                yieldValues = [yieldValues; final_NaOAc_quick; final_EtOH_quick; final_conversion_quick;
                              final_NaOH_remaining_quick; final_EtOAc_remaining_quick];
            end
            
            if ~isempty(yieldSummary)
                yieldTable = table(yieldSummary, yieldValues, 'VariableNames', {'Yield_Metric', 'Value'});
                writetable(yieldTable, fullPath, 'Sheet', '9_Yield_Summary', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write yield analysis: %s\n', ME.message);
        end
        
        %% 10. CONVERSION TO PRODUCT DATA
        app.updateLoadingDialog(0.75, 'Exporting conversion to product data...');
        
        try
            conversionProductData = table();
            
            % Manual conversion data
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                conv_NaOAc_manual = (app.ManualSimResults.C_NaOAc / app.Cin) * 100;
                conv_EtOH_manual = (app.ManualSimResults.C_EtOH / app.Cin) * 100;
                
                conversionProductData.Time_min = app.ManualSimResults.t;
                conversionProductData.Manual_NaOAc_Conv_percent = conv_NaOAc_manual;
                conversionProductData.Manual_EtOH_Conv_percent = conv_EtOH_manual;
            end
            
            % Optimized conversion data
            if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                conv_NaOAc_opt = (app.SimResults.C_NaOAc / app.Cin) * 100;
                conv_EtOH_opt = (app.SimResults.C_EtOH / app.Cin) * 100;
                % DEBUG: Check optimized conversion values
    fprintf('=== OPTIMIZED CONVERSION DEBUG ===\n');
    fprintf('Time points: %d\n', length(app.SimResults.t));
    fprintf('C_NaOAc range: [%.8f, %.8f]\n', min(app.SimResults.C_NaOAc), max(app.SimResults.C_NaOAc));
    fprintf('C_EtOH range: [%.8f, %.8f]\n', min(app.SimResults.C_EtOH), max(app.SimResults.C_EtOH));
    fprintf('app.Cin value: %.8f\n', app.Cin);
    fprintf('conv_NaOAc_opt range: [%.2f, %.2f]%%\n', min(conv_NaOAc_opt), max(conv_NaOAc_opt));
    fprintf('conv_EtOH_opt range: [%.2f, %.2f]%%\n', min(conv_EtOH_opt), max(conv_EtOH_opt));
    fprintf('Any NaN in conv_NaOAc_opt: %d\n', sum(isnan(conv_NaOAc_opt)));
    fprintf('Any NaN in conv_EtOH_opt: %d\n', sum(isnan(conv_EtOH_opt)));
    fprintf('================================\n');
                
                if isempty(conversionProductData)
                    conversionProductData.Time_min = app.SimResults.t;
                end
                conversionProductData.Optimized_NaOAc_Conv_percent = conv_NaOAc_opt;
                conversionProductData.Optimized_EtOH_Conv_percent = conv_EtOH_opt;
            end
            
            % Ultra-precision conversion data
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                conv_NaOAc_ultra = (app.UltraPrecisionResults.C_NaOAc / app.Cin) * 100;
                conv_EtOH_ultra = (app.UltraPrecisionResults.C_EtOH / app.Cin) * 100;
                % DEBUG: Check ultra-precision conversion values
    fprintf('=== ULTRA-PRECISION CONVERSION DEBUG ===\n');
    fprintf('Time points: %d\n', length(app.UltraPrecisionResults.t));
    fprintf('C_NaOAc range: [%.8f, %.8f]\n', min(app.UltraPrecisionResults.C_NaOAc), max(app.UltraPrecisionResults.C_NaOAc));
    fprintf('C_EtOH range: [%.8f, %.8f]\n', min(app.UltraPrecisionResults.C_EtOH), max(app.UltraPrecisionResults.C_EtOH));
    fprintf('app.Cin value: %.8f\n', app.Cin);
    fprintf('conv_NaOAc_ultra range: [%.2f, %.2f]%%\n', min(conv_NaOAc_ultra), max(conv_NaOAc_ultra));
    fprintf('conv_EtOH_ultra range: [%.2f, %.2f]%%\n', min(conv_EtOH_ultra), max(conv_EtOH_ultra));
    fprintf('Any NaN in conv_NaOAc_ultra: %d\n', sum(isnan(conv_NaOAc_ultra)));
    fprintf('Any NaN in conv_EtOH_ultra: %d\n', sum(isnan(conv_EtOH_ultra)));
    fprintf('========================================\n');
                
                if isempty(conversionProductData)
                    conversionProductData.Time_min = app.UltraPrecisionResults.t;
                end
                conversionProductData.Ultra_NaOAc_Conv_percent = conv_NaOAc_ultra;
                conversionProductData.Ultra_EtOH_Conv_percent = conv_EtOH_ultra;
            end
            
            % Quick sim conversion data
            if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                conv_NaOAc_quick = (app.QuickSimResults.C_NaOAc / app.Cin) * 100;
                conv_EtOH_quick = (app.QuickSimResults.C_EtOH / app.Cin) * 100;
                
                if isempty(conversionProductData)
                    conversionProductData.Time_min = app.QuickSimResults.t;
                end
                conversionProductData.QuickSim_NaOAc_Conv_percent = conv_NaOAc_quick;
                conversionProductData.QuickSim_EtOH_Conv_percent = conv_EtOH_quick;
            end
            
            if ~isempty(conversionProductData)
                writetable(conversionProductData, fullPath, 'Sheet', '10_Conv_Product', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write conversion to product data: %s\n', ME.message);
        end
        
        %% 11. CONVERSION UNREACTED DATA
        app.updateLoadingDialog(0.85, 'Exporting unreacted conversion data...');
        
        try
            conversionUnreactedData = table();
            
            % Manual unreacted conversion data
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                conv_NaOH_manual = ((app.Cin - app.ManualSimResults.C_NaOH) / app.Cin) * 100;
                conv_EtOAc_manual = ((app.Cin - app.ManualSimResults.C_EtOAc) / app.Cin) * 100;
                
                conversionUnreactedData.Time_min = app.ManualSimResults.t;
                conversionUnreactedData.Manual_NaOH_Conv_percent = conv_NaOH_manual;
                conversionUnreactedData.Manual_EtOAc_Conv_percent = conv_EtOAc_manual;
            end
            
            % Optimized unreacted conversion data
            if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                conv_NaOH_opt = ((app.Cin - app.SimResults.C_NaOH) / app.Cin) * 100;
                conv_EtOAc_opt = ((app.Cin - app.SimResults.C_EtOAc) / app.Cin) * 100;
                
                if isempty(conversionUnreactedData)
                    conversionUnreactedData.Time_min = app.SimResults.t;
                end
                conversionUnreactedData.Optimized_NaOH_Conv_percent = conv_NaOH_opt;
                conversionUnreactedData.Optimized_EtOAc_Conv_percent = conv_EtOAc_opt;
            end
            
            % Ultra-precision unreacted conversion data
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                conv_NaOH_ultra = ((app.Cin - app.UltraPrecisionResults.C_NaOH) / app.Cin) * 100;
                conv_EtOAc_ultra = ((app.Cin - app.UltraPrecisionResults.C_EtOAc) / app.Cin) * 100;
                
                if isempty(conversionUnreactedData)
                    conversionUnreactedData.Time_min = app.UltraPrecisionResults.t;
                end
                conversionUnreactedData.Ultra_NaOH_Conv_percent = conv_NaOH_ultra;
                conversionUnreactedData.Ultra_EtOAc_Conv_percent = conv_EtOAc_ultra;
            end
            
            % Quick sim unreacted conversion data
            if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                conv_NaOH_quick = ((app.Cin - app.QuickSimResults.C_NaOH) / app.Cin) * 100;
                conv_EtOAc_quick = ((app.Cin - app.QuickSimResults.C_EtOAc) / app.Cin) * 100;
                
                if isempty(conversionUnreactedData)
                    conversionUnreactedData.Time_min = app.QuickSimResults.t;
                end
                conversionUnreactedData.QuickSim_NaOH_Conv_percent = conv_NaOH_quick;
                conversionUnreactedData.QuickSim_EtOAc_Conv_percent = conv_EtOAc_quick;
            end
            
            if ~isempty(conversionUnreactedData)
                writetable(conversionUnreactedData, fullPath, 'Sheet', '11_Conv_Unreacted', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write unreacted conversion data: %s\n', ME.message);
        end
        
        %% NEW: 12. CONVERSION TO PRODUCT WITH CONVERGENCE DATA
app.updateLoadingDialog(0.875, 'Exporting conversion to product with convergence...');

try
    conversionProductConvergenceData = table();
    
    % Manual conversion to product data with convergence info
    if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
        conv_NaOAc_manual = (app.ManualSimResults.C_NaOAc / app.Cin) * 100;
        conv_EtOH_manual = (app.ManualSimResults.C_EtOH / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOAc_converged_manual = zeros(size(conv_NaOAc_manual));
        conv_EtOH_converged_manual = zeros(size(conv_EtOH_manual));
        
        if ~isempty(app.ManualConvergenceData)
            if isfield(app.ManualConvergenceData, 'NaOAc_converged') && app.ManualConvergenceData.NaOAc_converged
                conv_idx = find(app.ManualSimResults.t >= app.ManualConvergenceData.NaOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOAc_converged_manual(conv_idx:end) = 1;
                end
            end
            if isfield(app.ManualConvergenceData, 'EtOH_converged') && app.ManualConvergenceData.EtOH_converged
                conv_idx = find(app.ManualSimResults.t >= app.ManualConvergenceData.EtOH_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOH_converged_manual(conv_idx:end) = 1;
                end
            end
        end
        
        conversionProductConvergenceData.Time_min = app.ManualSimResults.t;
        conversionProductConvergenceData.Manual_NaOAc_Conv_percent = conv_NaOAc_manual;
        conversionProductConvergenceData.Manual_EtOH_Conv_percent = conv_EtOH_manual;
        conversionProductConvergenceData.Manual_NaOAc_Converged = conv_NaOAc_converged_manual;
        conversionProductConvergenceData.Manual_EtOH_Converged = conv_EtOH_converged_manual;
    end
    
    % Optimized conversion to product data with convergence info
    if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
        conv_NaOAc_opt = (app.SimResults.C_NaOAc / app.Cin) * 100;
        conv_EtOH_opt = (app.SimResults.C_EtOH / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOAc_converged_opt = zeros(size(conv_NaOAc_opt));
        conv_EtOH_converged_opt = zeros(size(conv_EtOH_opt));
        
        if ~isempty(app.ConvergenceData)
            if isfield(app.ConvergenceData, 'NaOAc_converged') && app.ConvergenceData.NaOAc_converged
                conv_idx = find(app.SimResults.t >= app.ConvergenceData.NaOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOAc_converged_opt(conv_idx:end) = 1;
                end
            end
            if isfield(app.ConvergenceData, 'EtOH_converged') && app.ConvergenceData.EtOH_converged
                conv_idx = find(app.SimResults.t >= app.ConvergenceData.EtOH_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOH_converged_opt(conv_idx:end) = 1;
                end
            end
        end
        
        if isempty(conversionProductConvergenceData)
            conversionProductConvergenceData.Time_min = app.SimResults.t;
        end
        conversionProductConvergenceData.Optimized_NaOAc_Conv_percent = conv_NaOAc_opt;
        conversionProductConvergenceData.Optimized_EtOH_Conv_percent = conv_EtOH_opt;
        conversionProductConvergenceData.Optimized_NaOAc_Converged = conv_NaOAc_converged_opt;
        conversionProductConvergenceData.Optimized_EtOH_Converged = conv_EtOH_converged_opt;
    end
    
    % Ultra-precision conversion to product data with convergence info
    if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
        conv_NaOAc_ultra = (app.UltraPrecisionResults.C_NaOAc / app.Cin) * 100;
        conv_EtOH_ultra = (app.UltraPrecisionResults.C_EtOH / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOAc_converged_ultra = zeros(size(conv_NaOAc_ultra));
        conv_EtOH_converged_ultra = zeros(size(conv_EtOH_ultra));
        
        if ~isempty(app.UltraConvergenceData)
            if isfield(app.UltraConvergenceData, 'NaOAc_converged') && app.UltraConvergenceData.NaOAc_converged
                conv_idx = find(app.UltraPrecisionResults.t >= app.UltraConvergenceData.NaOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOAc_converged_ultra(conv_idx:end) = 1;
                end
            end
            if isfield(app.UltraConvergenceData, 'EtOH_converged') && app.UltraConvergenceData.EtOH_converged
                conv_idx = find(app.UltraPrecisionResults.t >= app.UltraConvergenceData.EtOH_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOH_converged_ultra(conv_idx:end) = 1;
                end
            end
        end
        
        if isempty(conversionProductConvergenceData)
            conversionProductConvergenceData.Time_min = app.UltraPrecisionResults.t;
        end
        conversionProductConvergenceData.Ultra_NaOAc_Conv_percent = conv_NaOAc_ultra;
        conversionProductConvergenceData.Ultra_EtOH_Conv_percent = conv_EtOH_ultra;
        conversionProductConvergenceData.Ultra_NaOAc_Converged = conv_NaOAc_converged_ultra;
        conversionProductConvergenceData.Ultra_EtOH_Converged = conv_EtOH_converged_ultra;
    end
    
    if ~isempty(conversionProductConvergenceData)
        writetable(conversionProductConvergenceData, fullPath, 'Sheet', '12_Conv_Product_Convergence', 'WriteMode', 'overwritesheet');
    end
    
catch ME
    fprintf('Warning: Failed to write conversion to product with convergence data: %s\n', ME.message);
end

%% NEW: 13. CONVERSION UNREACTED WITH CONVERGENCE DATA
app.updateLoadingDialog(0.88, 'Exporting conversion unreacted with convergence...');

try
    conversionUnreactedConvergenceData = table();
    
    % Manual unreacted conversion data with convergence info
    if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
        conv_NaOH_manual = ((app.Cin - app.ManualSimResults.C_NaOH) / app.Cin) * 100;
        conv_EtOAc_manual = ((app.Cin - app.ManualSimResults.C_EtOAc) / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOH_converged_manual = zeros(size(conv_NaOH_manual));
        conv_EtOAc_converged_manual = zeros(size(conv_EtOAc_manual));
        
        if ~isempty(app.ManualConvergenceData)
            if isfield(app.ManualConvergenceData, 'NaOH_converged') && app.ManualConvergenceData.NaOH_converged
                conv_idx = find(app.ManualSimResults.t >= app.ManualConvergenceData.NaOH_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOH_converged_manual(conv_idx:end) = 1;
                end
            end
            if isfield(app.ManualConvergenceData, 'EtOAc_converged') && app.ManualConvergenceData.EtOAc_converged
                conv_idx = find(app.ManualSimResults.t >= app.ManualConvergenceData.EtOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOAc_converged_manual(conv_idx:end) = 1;
                end
            end
        end
        
        conversionUnreactedConvergenceData.Time_min = app.ManualSimResults.t;
        conversionUnreactedConvergenceData.Manual_NaOH_Conv_percent = conv_NaOH_manual;
        conversionUnreactedConvergenceData.Manual_EtOAc_Conv_percent = conv_EtOAc_manual;
        conversionUnreactedConvergenceData.Manual_NaOH_Converged = conv_NaOH_converged_manual;
        conversionUnreactedConvergenceData.Manual_EtOAc_Converged = conv_EtOAc_converged_manual;
    end
    
    % Optimized unreacted conversion data with convergence info
    if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
        conv_NaOH_opt = ((app.Cin - app.SimResults.C_NaOH) / app.Cin) * 100;
        conv_EtOAc_opt = ((app.Cin - app.SimResults.C_EtOAc) / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOH_converged_opt = zeros(size(conv_NaOH_opt));
        conv_EtOAc_converged_opt = zeros(size(conv_EtOAc_opt));
        
        if ~isempty(app.ConvergenceData)
            if isfield(app.ConvergenceData, 'NaOH_converged') && app.ConvergenceData.NaOH_converged
                conv_idx = find(app.SimResults.t >= app.ConvergenceData.NaOH_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOH_converged_opt(conv_idx:end) = 1;
                end
            end
            if isfield(app.ConvergenceData, 'EtOAc_converged') && app.ConvergenceData.EtOAc_converged
                conv_idx = find(app.SimResults.t >= app.ConvergenceData.EtOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOAc_converged_opt(conv_idx:end) = 1;
                end
            end
        end
        
        if isempty(conversionUnreactedConvergenceData)
            conversionUnreactedConvergenceData.Time_min = app.SimResults.t;
        end
        conversionUnreactedConvergenceData.Optimized_NaOH_Conv_percent = conv_NaOH_opt;
        conversionUnreactedConvergenceData.Optimized_EtOAc_Conv_percent = conv_EtOAc_opt;
        conversionUnreactedConvergenceData.Optimized_NaOH_Converged = conv_NaOH_converged_opt;
        conversionUnreactedConvergenceData.Optimized_EtOAc_Converged = conv_EtOAc_converged_opt;
    end
    
    % Ultra-precision unreacted conversion data with convergence info
    if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
        conv_NaOH_ultra = ((app.Cin - app.UltraPrecisionResults.C_NaOH) / app.Cin) * 100;
        conv_EtOAc_ultra = ((app.Cin - app.UltraPrecisionResults.C_EtOAc) / app.Cin) * 100;
        
        % Add convergence markers
        conv_NaOH_converged_ultra = zeros(size(conv_NaOH_ultra));
        conv_EtOAc_converged_ultra = zeros(size(conv_EtOAc_ultra));
        
        if ~isempty(app.UltraConvergenceData)
            if isfield(app.UltraConvergenceData, 'NaOH_converged') && app.UltraConvergenceData.NaOH_converged
                conv_idx = find(app.UltraPrecisionResults.t >= app.UltraConvergenceData.NaOH_time, 1);
                if ~isempty(conv_idx)
                    conv_NaOH_converged_ultra(conv_idx:end) = 1;
                end
            end
            if isfield(app.UltraConvergenceData, 'EtOAc_converged') && app.UltraConvergenceData.EtOAc_converged
                conv_idx = find(app.UltraPrecisionResults.t >= app.UltraConvergenceData.EtOAc_time, 1);
                if ~isempty(conv_idx)
                    conv_EtOAc_converged_ultra(conv_idx:end) = 1;
                end
            end
        end
        
        if isempty(conversionUnreactedConvergenceData)
            conversionUnreactedConvergenceData.Time_min = app.UltraPrecisionResults.t;
        end
        conversionUnreactedConvergenceData.Ultra_NaOH_Conv_percent = conv_NaOH_ultra;
        conversionUnreactedConvergenceData.Ultra_EtOAc_Conv_percent = conv_EtOAc_ultra;
        conversionUnreactedConvergenceData.Ultra_NaOH_Converged = conv_NaOH_converged_ultra;
        conversionUnreactedConvergenceData.Ultra_EtOAc_Converged = conv_EtOAc_converged_ultra;
    end
    
    if ~isempty(conversionUnreactedConvergenceData)
        writetable(conversionUnreactedConvergenceData, fullPath, 'Sheet', '13_Conv_Unreacted_Convergence', 'WriteMode', 'overwritesheet');
    end
    
catch ME
    fprintf('Warning: Failed to write conversion unreacted with convergence data: %s\n', ME.message);
end

%% NEW: 14. CONVERGENCE SUMMARY TABLE
app.updateLoadingDialog(0.885, 'Exporting convergence summary...');

try
    convergenceSummaryName = {};
    convergenceSummaryValue = {};
    
    % Manual convergence summary
    if ~isempty(app.ManualConvergenceData)
        convergenceSummaryName = [convergenceSummaryName; {'Manual NaOH Converged'; 'Manual NaOH Convergence Time (min)'; 'Manual NaOH Convergence Concentration (mol/L)'}];
        if isfield(app.ManualConvergenceData, 'NaOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.NaOH_converged}];
            if app.ManualConvergenceData.NaOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.NaOH_time; app.ManualConvergenceData.NaOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Manual EtOAc Converged'; 'Manual EtOAc Convergence Time (min)'; 'Manual EtOAc Convergence Concentration (mol/L)'}];
        if isfield(app.ManualConvergenceData, 'EtOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.EtOAc_converged}];
            if app.ManualConvergenceData.EtOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.EtOAc_time; app.ManualConvergenceData.EtOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Manual NaOAc Converged'; 'Manual NaOAc Convergence Time (min)'; 'Manual NaOAc Convergence Concentration (mol/L)'}];
        if isfield(app.ManualConvergenceData, 'NaOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.NaOAc_converged}];
            if app.ManualConvergenceData.NaOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.NaOAc_time; app.ManualConvergenceData.NaOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Manual EtOH Converged'; 'Manual EtOH Convergence Time (min)'; 'Manual EtOH Convergence Concentration (mol/L)'}];
        if isfield(app.ManualConvergenceData, 'EtOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.EtOH_converged}];
            if app.ManualConvergenceData.EtOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ManualConvergenceData.EtOH_time; app.ManualConvergenceData.EtOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
    end
    
    % Optimized convergence summary
    if ~isempty(app.ConvergenceData)
        convergenceSummaryName = [convergenceSummaryName; {'Optimized NaOH Converged'; 'Optimized NaOH Convergence Time (min)'; 'Optimized NaOH Convergence Concentration (mol/L)'}];
        if isfield(app.ConvergenceData, 'NaOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.NaOH_converged}];
            if app.ConvergenceData.NaOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.NaOH_time; app.ConvergenceData.NaOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Optimized EtOAc Converged'; 'Optimized EtOAc Convergence Time (min)'; 'Optimized EtOAc Convergence Concentration (mol/L)'}];
        if isfield(app.ConvergenceData, 'EtOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.EtOAc_converged}];
            if app.ConvergenceData.EtOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.EtOAc_time; app.ConvergenceData.EtOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Optimized NaOAc Converged'; 'Optimized NaOAc Convergence Time (min)'; 'Optimized NaOAc Convergence Concentration (mol/L)'}];
        if isfield(app.ConvergenceData, 'NaOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.NaOAc_converged}];
            if app.ConvergenceData.NaOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.NaOAc_time; app.ConvergenceData.NaOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Optimized EtOH Converged'; 'Optimized EtOH Convergence Time (min)'; 'Optimized EtOH Convergence Concentration (mol/L)'}];
        if isfield(app.ConvergenceData, 'EtOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.EtOH_converged}];
            if app.ConvergenceData.EtOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.ConvergenceData.EtOH_time; app.ConvergenceData.EtOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
    end
    
    % Ultra-precision convergence summary
    if ~isempty(app.UltraConvergenceData)
        convergenceSummaryName = [convergenceSummaryName; {'Ultra NaOH Converged'; 'Ultra NaOH Convergence Time (min)'; 'Ultra NaOH Convergence Concentration (mol/L)'}];
        if isfield(app.UltraConvergenceData, 'NaOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.NaOH_converged}];
            if app.UltraConvergenceData.NaOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.NaOH_time; app.UltraConvergenceData.NaOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Ultra EtOAc Converged'; 'Ultra EtOAc Convergence Time (min)'; 'Ultra EtOAc Convergence Concentration (mol/L)'}];
        if isfield(app.UltraConvergenceData, 'EtOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.EtOAc_converged}];
            if app.UltraConvergenceData.EtOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.EtOAc_time; app.UltraConvergenceData.EtOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Ultra NaOAc Converged'; 'Ultra NaOAc Convergence Time (min)'; 'Ultra NaOAc Convergence Concentration (mol/L)'}];
        if isfield(app.UltraConvergenceData, 'NaOAc_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.NaOAc_converged}];
            if app.UltraConvergenceData.NaOAc_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.NaOAc_time; app.UltraConvergenceData.NaOAc_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
        
        convergenceSummaryName = [convergenceSummaryName; {'Ultra EtOH Converged'; 'Ultra EtOH Convergence Time (min)'; 'Ultra EtOH Convergence Concentration (mol/L)'}];
        if isfield(app.UltraConvergenceData, 'EtOH_converged')
            convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.EtOH_converged}];
            if app.UltraConvergenceData.EtOH_converged
                convergenceSummaryValue = [convergenceSummaryValue; {app.UltraConvergenceData.EtOH_time; app.UltraConvergenceData.EtOH_conc}];
            else
                convergenceSummaryValue = [convergenceSummaryValue; {'N/A'; 'N/A'}];
            end
        else
            convergenceSummaryValue = [convergenceSummaryValue; {false; 'N/A'; 'N/A'}];
        end
    end
    
    if ~isempty(convergenceSummaryName)
        convergenceSummaryTable = table(convergenceSummaryName, convergenceSummaryValue, 'VariableNames', {'Convergence_Parameter', 'Value'});
        writetable(convergenceSummaryTable, fullPath, 'Sheet', '14_Convergence_Summary', 'WriteMode', 'overwritesheet');
    end
    
catch ME
    fprintf('Warning: Failed to write convergence summary: %s\n', ME.message);
end

        
        %% NEW: 15. TEMPERATURE PROFILES DATA
        app.updateLoadingDialog(0.87, 'Exporting temperature profiles...');
        
        try
            temperatureProfileData = table();
            
            % Manual simulation temperature data
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                temperatureProfileData.Time_min = app.ManualSimResults.t;
                temperatureProfileData.Manual_Reactor_Temp_C = app.ManualSimResults.T_reactor;
                temperatureProfileData.Manual_Jacket_Temp_C = app.ManualSimResults.T_jacket;
            end
            
            % Optimized simulation temperature data
            if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                if isempty(temperatureProfileData)
                    temperatureProfileData.Time_min = app.SimResults.t;
                end
                temperatureProfileData.Optimized_Reactor_Temp_C = app.SimResults.T_reactor;
                temperatureProfileData.Optimized_Jacket_Temp_C = app.SimResults.T_jacket;
            end
            
            % Ultra-precision simulation temperature data
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                if isempty(temperatureProfileData)
                    temperatureProfileData.Time_min = app.UltraPrecisionResults.t;
                end
                temperatureProfileData.Ultra_Reactor_Temp_C = app.UltraPrecisionResults.T_reactor;
                temperatureProfileData.Ultra_Jacket_Temp_C = app.UltraPrecisionResults.T_jacket;
            end
            
            % Quick simulation temperature data
            if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                if isempty(temperatureProfileData)
                    temperatureProfileData.Time_min = app.QuickSimResults.t;
                end
                temperatureProfileData.QuickSim_Reactor_Temp_C = app.QuickSimResults.T_reactor;
                temperatureProfileData.QuickSim_Jacket_Temp_C = app.QuickSimResults.T_jacket;
            end
            
            if ~isempty(temperatureProfileData)
                writetable(temperatureProfileData, fullPath, 'Sheet', '12_Temperature_Profiles', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write temperature profiles: %s\n', ME.message);
        end
        
        %% NEW: 16. TEMPERATURE STATISTICS
        app.updateLoadingDialog(0.89, 'Exporting temperature statistics...');
        
        try
            tempStatsName = {};
            tempStatsValue = {};
            
            % Manual simulation temperature statistics
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 'T_reactor')
                tempStatsName = [tempStatsName; 'Manual Reactor Max Temp (°C)'; 'Manual Reactor Min Temp (°C)'; 'Manual Reactor Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.ManualSimResults.T_reactor); min(app.ManualSimResults.T_reactor); mean(app.ManualSimResults.T_reactor)];
                
                tempStatsName = [tempStatsName; 'Manual Jacket Max Temp (°C)'; 'Manual Jacket Min Temp (°C)'; 'Manual Jacket Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.ManualSimResults.T_jacket); min(app.ManualSimResults.T_jacket); mean(app.ManualSimResults.T_jacket)];
                
                % Additional manual statistics
                manual_temp_rise = max(app.ManualSimResults.T_reactor) - app.ManualSimResults.T_reactor(1);
                manual_temp_diff = max(app.ManualSimResults.T_reactor) - max(app.ManualSimResults.T_jacket);
                tempStatsName = [tempStatsName; 'Manual Max Temperature Rise (°C)'; 'Manual Max Reactor-Jacket Diff (°C)'];
                tempStatsValue = [tempStatsValue; manual_temp_rise; manual_temp_diff];
            end
            
            % Optimized simulation temperature statistics
            if ~isempty(app.SimResults) && isfield(app.SimResults, 'T_reactor')
                tempStatsName = [tempStatsName; 'Optimized Reactor Max Temp (°C)'; 'Optimized Reactor Min Temp (°C)'; 'Optimized Reactor Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.SimResults.T_reactor); min(app.SimResults.T_reactor); mean(app.SimResults.T_reactor)];
                
                tempStatsName = [tempStatsName; 'Optimized Jacket Max Temp (°C)'; 'Optimized Jacket Min Temp (°C)'; 'Optimized Jacket Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.SimResults.T_jacket); min(app.SimResults.T_jacket); mean(app.SimResults.T_jacket)];
                
                % Additional optimized statistics
                opt_temp_rise = max(app.SimResults.T_reactor) - app.SimResults.T_reactor(1);
                opt_temp_diff = max(app.SimResults.T_reactor) - max(app.SimResults.T_jacket);
                tempStatsName = [tempStatsName; 'Optimized Max Temperature Rise (°C)'; 'Optimized Max Reactor-Jacket Diff (°C)'];
                tempStatsValue = [tempStatsValue; opt_temp_rise; opt_temp_diff];
                
                % Temperature control metrics (if isothermal mode was used)
                if app.IsothermalMode
                    temp_deviation = max(abs(app.SimResults.T_reactor - app.SelectedTemp));
                    tempStatsName = [tempStatsName; 'Optimized Max Temp Deviation from Setpoint (°C)'; 'Optimized Isothermal Control Quality'];
                    control_quality = temp_deviation <= app.isothermal_tolerance;
                    tempStatsValue = [tempStatsValue; temp_deviation; control_quality];
                end
            end
            
            % Ultra-precision simulation temperature statistics
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 'T_reactor')
                tempStatsName = [tempStatsName; 'Ultra Reactor Max Temp (°C)'; 'Ultra Reactor Min Temp (°C)'; 'Ultra Reactor Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.UltraPrecisionResults.T_reactor); min(app.UltraPrecisionResults.T_reactor); mean(app.UltraPrecisionResults.T_reactor)];
                
                tempStatsName = [tempStatsName; 'Ultra Jacket Max Temp (°C)'; 'Ultra Jacket Min Temp (°C)'; 'Ultra Jacket Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.UltraPrecisionResults.T_jacket); min(app.UltraPrecisionResults.T_jacket); mean(app.UltraPrecisionResults.T_jacket)];
                
                % Additional ultra-precision statistics
                ultra_temp_rise = max(app.UltraPrecisionResults.T_reactor) - app.UltraPrecisionResults.T_reactor(1);
                ultra_temp_diff = max(app.UltraPrecisionResults.T_reactor) - max(app.UltraPrecisionResults.T_jacket);
                tempStatsName = [tempStatsName; 'Ultra Max Temperature Rise (°C)'; 'Ultra Max Reactor-Jacket Diff (°C)'];
                tempStatsValue = [tempStatsValue; ultra_temp_rise; ultra_temp_diff];
            end
            
            % Quick simulation temperature statistics
            if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 'T_reactor')
                tempStatsName = [tempStatsName; 'QuickSim Reactor Max Temp (°C)'; 'QuickSim Reactor Min Temp (°C)'; 'QuickSim Reactor Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.QuickSimResults.T_reactor); min(app.QuickSimResults.T_reactor); mean(app.QuickSimResults.T_reactor)];
                
                tempStatsName = [tempStatsName; 'QuickSim Jacket Max Temp (°C)'; 'QuickSim Jacket Min Temp (°C)'; 'QuickSim Jacket Avg Temp (°C)'];
                tempStatsValue = [tempStatsValue; max(app.QuickSimResults.T_jacket); min(app.QuickSimResults.T_jacket); mean(app.QuickSimResults.T_jacket)];
                
                % Additional quick sim statistics
                quick_temp_rise = max(app.QuickSimResults.T_reactor) - app.QuickSimResults.T_reactor(1);
                quick_temp_diff = max(app.QuickSimResults.T_reactor) - max(app.QuickSimResults.T_jacket);
                tempStatsName = [tempStatsName; 'QuickSim Max Temperature Rise (°C)'; 'QuickSim Max Reactor-Jacket Diff (°C)'];
                tempStatsValue = [tempStatsValue; quick_temp_rise; quick_temp_diff];
            end
            
            % Cooling system statistics (if adaptive cooling was used)
            if app.IsothermalMode && ~isempty(app.SimResults)
                tempStatsName = [tempStatsName; 'Adaptive Coolant Flowrate (L/min)'; 'Cooling Adequacy Factor'; 'Required Cooling Capacity (J/min)'; 'Max Cooling Capacity (J/min)'];
                tempStatsValue = [tempStatsValue; app.Fj_adaptive; app.cooling_adequacy_factor; app.required_cooling_capacity; app.max_cooling_capacity];
            end
            
            if ~isempty(tempStatsName)
                temperatureStatsTable = table(tempStatsName, tempStatsValue, 'VariableNames', {'Temperature_Statistic', 'Value'});
                writetable(temperatureStatsTable, fullPath, 'Sheet', '13_Temperature_Stats', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write temperature statistics: %s\n', ME.message);
        end
        
        %% 17. OPTIMIZED PARAMETERS (Updated sheet number)
        app.updateLoadingDialog(0.95, 'Exporting optimized parameters...');
        
        try
            optimizedParams = {};
            optimizedValues = {};
            
            % Optimization results if available
            if ~isempty(app.OptimizationResults)
                if isfield(app.OptimizationResults, 'InitialFlowrate')
                    optimizedParams = [optimizedParams; 'Initial Flowrate (L/min)'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.InitialFlowrate];
                end
                if isfield(app.OptimizationResults, 'OptimalFlowrate')
                    optimizedParams = [optimizedParams; 'Optimal Flowrate (L/min)'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.OptimalFlowrate];
                end
                if isfield(app.OptimizationResults, 'OptimalK')
                    optimizedParams = [optimizedParams; 'Optimal Rate Constant (L/mol·min)'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.OptimalK];
                end
                if isfield(app.OptimizationResults, 'MaxError')
                    optimizedParams = [optimizedParams; 'Maximum Error (%)'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.MaxError];
                end
                if isfield(app.OptimizationResults, 'RSquared')
                    optimizedParams = [optimizedParams; 'R-Squared'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.RSquared];
                end
                if isfield(app.OptimizationResults, 'Method')
                    optimizedParams = [optimizedParams; 'Optimization Method'];
                    optimizedValues = [optimizedValues; {app.OptimizationResults.Method}];
                end
                if isfield(app.OptimizationResults, 'AllErrorsBelow1Percent')
                    optimizedParams = [optimizedParams; 'All Errors Below 1%'; 'All Errors Below 0.5%'];
                    optimizedValues = [optimizedValues; app.OptimizationResults.AllErrorsBelow1Percent;
                                      app.OptimizationResults.AllErrorsBelow05Percent];
                end
            end
            
            % Add current optimized values even if optimization results not stored
            if ~isnan(app.k_fit)
                optimizedParams = [optimizedParams; 'Current Optimized k (L/mol·min)'; 'Current Optimized R²'];
                optimizedValues = [optimizedValues; app.k_fit; app.RSquared];
            end
            
            if ~isnan(app.k_ultra)
                optimizedParams = [optimizedParams; 'Ultra-Precision k (L/mol·min)'; 'Ultra-Precision R²'];
                optimizedValues = [optimizedValues; app.k_ultra; app.RSquaredUltra];
            end
            
            if ~isnan(app.k_manual)
                optimizedParams = [optimizedParams; 'Manual k (L/mol·min)'; 'Manual R²'];
                optimizedValues = [optimizedValues; app.k_manual; app.RSquaredManual];
            end
            
            % Add reactor operating conditions
            optimizedParams = [optimizedParams; 'Operating Temperature (°C)'; 'Reactor Volume (L)'; 
                              'Total Flowrate (L/min)'; 'Residence Time (min)'; 'Inlet Concentration (mol/L)'];
            optimizedValues = [optimizedValues; app.SelectedTemp; app.V; app.V0; app.V/app.V0; app.Cin];
            
            if ~isempty(optimizedParams)
                optimizedTable = table(optimizedParams, optimizedValues, 'VariableNames', {'Parameter', 'Value'});
                writetable(optimizedTable, fullPath, 'Sheet', '14_Optimized_Params', 'WriteMode', 'overwritesheet');
            end
            
        catch ME
            fprintf('Warning: Failed to write optimized parameters: %s\n', ME.message);
        end
        
        app.updateLoadingDialog(1.0, 'Export complete!');
        app.closeLoadingDialog();
        
   app.updateStatus('📊 Comprehensive data exported successfully with cooling system analysis!', 'success');
app.flashSuccess(app.ExportDataButton);

msgText = sprintf(['✅ Comprehensive Export Complete!\n\n' ...
                  'Exported 17 sheets with:\n' ...
                  '• Reactor parameters & residence times\n' ...
                  '• All simulation results (Manual, Optimized, Ultra, Quick)\n' ...
                  '• Complete error analysis & statistics\n' ...
                  '• Detailed convergence analysis with thresholds\n' ...
                  '• Yield summary for all methods\n' ...
                  '• Conversion to product data\n' ...
                  '• Conversion unreacted data\n' ...
                  '• 🌡️ Temperature profiles\n' ...
                  '• 📊 Temperature statistics\n' ...
                  '• ❄️ Cooling system analysis\n' ...
                  '• 📈 Cooling flowrate time series\n' ...
                  '• ⚡ Energy & Mass Balance Analysis (NEW)\n' ...
                  '• Optimized parameters summary\n\n' ...
                  'Energy & Mass Balance Analysis Includes:\n' ...
                  '• Reactor energy balance with heat generation/removal\n' ...
                  '• Component mass balance verification\n' ...
                  '• Stoichiometric consistency checks\n' ...
                  '• Reactor performance indicators (conversion, yield, selectivity)\n' ...
                  '• Space-time yield calculations\n\n' ...
                  'File saved: %s'], file);
        
        uialert(app.UIFigure, msgText, 'Export Complete', 'Icon', 'success');
        
    catch ME
        app.closeLoadingDialog();
        app.updateStatus(['Export failed: ' ME.message], 'error');
        uialert(app.UIFigure, ['Export failed: ', ME.message], 'Export Error', 'Icon', 'error');
    end
    
    %% NEW: 18. COOLING SYSTEM ANALYSIS
app.updateLoadingDialog(0.97, 'Exporting cooling system analysis...');

try
    coolingAnalysisName = {};
    coolingAnalysisValue = {};
    
    % Cooling system summary statistics
    if app.IsothermalMode && ~isempty(app.CoolingFlowrateHistory)
        coolingAnalysisName = [coolingAnalysisName; 'Isothermal Mode Active'; 'Average Cooling Flowrate (L/min)'; 
                              'Max Cooling Flowrate (L/min)'; 'Min Cooling Flowrate (L/min)'];
        coolingAnalysisValue = [coolingAnalysisValue; true; mean(app.CoolingFlowrateHistory); 
                               max(app.CoolingFlowrateHistory); min(app.CoolingFlowrateHistory)];
        
        coolingAnalysisName = [coolingAnalysisName; 'Average Required Capacity (J/min)'; 'Max Required Capacity (J/min)';
                              'Average Adequacy Factor'; 'Min Adequacy Factor'];
        coolingAnalysisValue = [coolingAnalysisValue; mean(app.CoolingCapacityHistory); max(app.CoolingCapacityHistory);
                               mean(app.CoolingAdequacyHistory); min(app.CoolingAdequacyHistory)];
        
        coolingAnalysisName = [coolingAnalysisName; 'Average Temp Deviation (°C)'; 'Max Temp Deviation (°C)';
                              'Temperature Control Quality (%)'; 'Cooling System Utilization (%)'];
        avg_temp_dev = mean(app.TemperatureDeviationHistory);
        max_temp_dev = max(app.TemperatureDeviationHistory);
        temp_control_quality = 100 * sum(app.TemperatureDeviationHistory <= app.isothermal_tolerance) / length(app.TemperatureDeviationHistory);
        cooling_utilization = 100 * mean(app.CoolingFlowrateHistory) / app.Fj_max;
        coolingAnalysisValue = [coolingAnalysisValue; avg_temp_dev; max_temp_dev; temp_control_quality; cooling_utilization];
    else
        coolingAnalysisName = [coolingAnalysisName; 'Isothermal Mode Active'; 'Fixed Cooling Flowrate (L/min)'];
        coolingAnalysisValue = [coolingAnalysisValue; false; app.Fj];
    end
    
    % Add Arrhenius parameters
    coolingAnalysisName = [coolingAnalysisName; 'Activation Energy (kJ/mol)'; 'Pre-exponential Factor (L/mol·min)'];
    coolingAnalysisValue = [coolingAnalysisValue; app.Ea/1000; app.A_preexp];
    
    % Add cooling system design parameters
    coolingAnalysisName = [coolingAnalysisName; 'Min Coolant Flow (L/min)'; 'Max Coolant Flow (L/min)';
                          'Coolant Inlet Temp (°C)'; 'Control Tolerance (°C)'; 'PI Controller Gain'];
    coolingAnalysisValue = [coolingAnalysisValue; app.Fj_min; app.Fj_max; app.Tj0; 
                           app.isothermal_tolerance; app.cooling_control_gain];
    
    if ~isempty(coolingAnalysisName)
        coolingAnalysisTable = table(coolingAnalysisName, coolingAnalysisValue, 'VariableNames', {'Cooling_Parameter', 'Value'});
        writetable(coolingAnalysisTable, fullPath, 'Sheet', '15_Cooling_Analysis', 'WriteMode', 'overwritesheet');
    end
    
catch ME
    fprintf('Warning: Failed to write cooling analysis: %s\n', ME.message);
end

%% NEW: 19. COOLING FLOWRATE TIME SERIES DATA
app.updateLoadingDialog(0.98, 'Exporting cooling flowrate time series...');

try
    if ~isempty(app.CoolingTimeHistory) && ~isempty(app.CoolingFlowrateHistory)
        % Create time series data for cooling system
        coolingTimeSeriesData = table(app.CoolingTimeHistory(:), app.CoolingFlowrateHistory(:), ...
                                     app.CoolingCapacityHistory(:), app.CoolingAdequacyHistory(:), ...
                                     app.TemperatureDeviationHistory(:), ...
                                     'VariableNames', {'Time_min', 'Cooling_Flowrate_Lmin', ...
                                     'Required_Capacity_Jmin', 'Adequacy_Factor', 'Temp_Deviation_C'});
        
        % Add calculated rate constants using Arrhenius equation if temperature data exists
        if ~isempty(app.SimResults) && isfield(app.SimResults, 'T_reactor')
            % Interpolate temperature data to cooling time points
            T_reactor_interp = interp1(app.SimResults.t, app.SimResults.T_reactor, ...
                                      app.CoolingTimeHistory, 'linear', 'extrap');
            
            % Calculate Arrhenius rate constants
            k_arrhenius = zeros(size(T_reactor_interp));
            for i = 1:length(T_reactor_interp)
                k_arrhenius(i) = app.calculateArrheniusRateConstant(T_reactor_interp(i));
            end
            
            coolingTimeSeriesData.Reactor_Temp_C = T_reactor_interp(:);
            coolingTimeSeriesData.Arrhenius_k_Lmolmin = k_arrhenius(:);
        end
        
        writetable(coolingTimeSeriesData, fullPath, 'Sheet', '16_Cooling_TimeSeries', 'WriteMode', 'overwritesheet');
    end
    
catch ME
    fprintf('Warning: Failed to write cooling time series: %s\n', ME.message);
end

app.updateLoadingDialog(1.0, 'Export complete!');
        app.closeLoadingDialog();
    
       end

    end
    
    % Component initialization
    methods (Access = private)
        
        function storeCoolingHistory(app, time, flowrate, required_capacity, adequacy_factor, temp_deviation)
    try
        % Initialize arrays if this is the first call
        if isempty(app.CoolingTimeHistory)
            app.CoolingTimeHistory = time;
            app.CoolingFlowrateHistory = flowrate;
            app.CoolingCapacityHistory = required_capacity;
            app.CoolingAdequacyHistory = adequacy_factor;
            app.TemperatureDeviationHistory = temp_deviation;
        else
            % Append new data (limit storage to avoid memory issues)
            if length(app.CoolingTimeHistory) < 10000  % Limit to 10k points
                app.CoolingTimeHistory(end+1) = time;
                app.CoolingFlowrateHistory(end+1) = flowrate;
                app.CoolingCapacityHistory(end+1) = required_capacity;
                app.CoolingAdequacyHistory(end+1) = adequacy_factor;
                app.TemperatureDeviationHistory(end+1) = temp_deviation;
            end
        end
    catch
        % Silent fail to avoid disrupting simulation
    end
        end
        
        function storeCoolingHistoryEnhanced(app, time, flowrate, required_capacity, adequacy_factor, ...
                                    temp_deviation, temp_difference, total_cooling_required)
    try
        % Initialize arrays if this is the first call
        if isempty(app.CoolingTimeHistory)
            app.CoolingTimeHistory = time;
            app.CoolingFlowrateHistory = flowrate;
            app.CoolingCapacityHistory = required_capacity;
            app.CoolingAdequacyHistory = adequacy_factor;
            app.TemperatureDeviationHistory = temp_deviation;
            % Add new tracking variables
            app.TempDifferenceHistory = temp_difference;
            app.TotalCoolingRequiredHistory = total_cooling_required;
        else
            % Append new data (limit storage to avoid memory issues)
            if length(app.CoolingTimeHistory) < 10000  % Limit to 10k points
                app.CoolingTimeHistory(end+1) = time;
                app.CoolingFlowrateHistory(end+1) = flowrate;
                app.CoolingCapacityHistory(end+1) = required_capacity;
                app.CoolingAdequacyHistory(end+1) = adequacy_factor;
                app.TemperatureDeviationHistory(end+1) = temp_deviation;
                app.TempDifferenceHistory(end+1) = temp_difference;
                app.TotalCoolingRequiredHistory(end+1) = total_cooling_required;
            end
        end
    catch
        % Silent fail to avoid disrupting simulation
    end
end

        
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1600 1000];  % Increased width for new tabs
            app.UIFigure.Name = 'Modeling and Simulation of Saponification of Ethyl Acetate and Sodium Hydroxide in Continuos Stirred Tank Reactor for Final Year Project (FYP2) : MUHAMMAD LUQMAN HAKIMI BIN MOHD ZAMRI_2021817262_CEEH2208C';
            
            % Create main grid layout
            app.MainGridLayout = uigridlayout(app.UIFigure);
            app.MainGridLayout.ColumnWidth = {'0.5x', '2.5x'};
            app.MainGridLayout.RowHeight = {'1x'};
            app.MainGridLayout.ColumnSpacing = 10;
            app.MainGridLayout.Padding = [10 10 10 10];
            
            % Create left panel
            app.LeftPanel = uipanel(app.MainGridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Title = '';
            
            app.LeftGridLayout = uigridlayout(app.LeftPanel);
            app.LeftGridLayout.RowHeight = {45, 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.LeftGridLayout.ColumnWidth = {'1x'};
            app.LeftGridLayout.RowSpacing = 15;
            app.LeftGridLayout.Padding = [15 15 15 15];
            
            % Create logo section
            app.LogoHTML = uihtml(app.LeftGridLayout);
            app.LogoHTML.Layout.Row = 1;
            app.LogoHTML.Layout.Column = 1;
            app.LogoHTML.Position = [0 0 300 30];
            
            % Create data upload section
            app.DataPanel = uipanel(app.LeftGridLayout);
            app.DataPanel.Layout.Row = 2;
            app.DataPanel.Layout.Column = 1;
            app.DataPanel.Title = 'Data Upload';
            app.DataPanel.FontWeight = 'bold';
            
            app.DataGrid = uigridlayout(app.DataPanel);
           app.DataGrid.RowHeight = {'fit', 'fit', 'fit'};
            app.DataGrid.ColumnWidth = {'1x'};
            app.DataGrid.Padding = [10 10 10 10];
            
            app.UploadButton = uibutton(app.DataGrid, 'push');
            app.UploadButton.Layout.Row = 1;
            app.UploadButton.Layout.Column = 1;
            app.UploadButton.Text = '📁 Upload Experimental Data';
            app.UploadButton.FontWeight = 'bold';
            app.UploadButton.ButtonPushedFcn = createCallbackFcn(app, @UploadButtonPushed, true);
            
            app.FileLabel = uilabel(app.DataGrid);
            app.FileLabel.Layout.Row = 2;
            app.FileLabel.Layout.Column = 1;
            app.FileLabel.Text = 'No file selected';
            app.FileLabel.FontAngle = 'italic';
            app.FileLabel.HorizontalAlignment = 'center';
            
            % ADD AFTER app.FileLabel creation:
app.ValidateDataButton = uibutton(app.DataGrid, 'push');
app.ValidateDataButton.Layout.Row = 3;
app.ValidateDataButton.Layout.Column = 1;
app.ValidateDataButton.Text = '🔍 Validate Data Quality';
app.ValidateDataButton.FontWeight = 'bold';
app.ValidateDataButton.ButtonPushedFcn = createCallbackFcn(app, @ValidateDataButtonPushed, true);
            
            % Create parameters section
            app.ParamsPanel = uipanel(app.LeftGridLayout);
            app.ParamsPanel.Layout.Row = 3;
            app.ParamsPanel.Layout.Column = 1;
            app.ParamsPanel.Title = 'Reactor Parameters';
            app.ParamsPanel.FontWeight = 'bold';
            
            app.ParamsGrid = uigridlayout(app.ParamsPanel);
            app.ParamsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
            app.ParamsGrid.ColumnWidth = {'1x', '1x'};
            app.ParamsGrid.Padding = [10 10 10 10];
            
            % Temperature dropdown
            app.TempLabel = uilabel(app.ParamsGrid);
            app.TempLabel.Layout.Row = 1;
            app.TempLabel.Layout.Column = 1;
            app.TempLabel.Text = 'Temperature (°C)';
            app.TempLabel.FontWeight = 'bold';
            
            app.TempDropDown = uidropdown(app.ParamsGrid);
            app.TempDropDown.Layout.Row = 1;
            app.TempDropDown.Layout.Column = 2;
            app.TempDropDown.Items = {'30', '35', '40', '45', '50', '55', '60', '65', '70'};
            app.TempDropDown.Value = '30';
            app.TempDropDown.ValueChangedFcn = createCallbackFcn(app, @TempDropDownValueChanged, true);
            
            % Volume field
            app.VolumeLabel = uilabel(app.ParamsGrid);
            app.VolumeLabel.Layout.Row = 2;
            app.VolumeLabel.Layout.Column = 1;
            app.VolumeLabel.Text = 'Volume (L)';
            app.VolumeLabel.FontWeight = 'bold';
            
            app.VolumeField = uieditfield(app.ParamsGrid, 'numeric');
            app.VolumeField.Layout.Row = 2;
            app.VolumeField.Layout.Column = 2;
            app.VolumeField.Value = 2.5;
            app.VolumeField.Limits = [0.1 10];
            
            % Flowrate field
            app.FlowrateLabel = uilabel(app.ParamsGrid);
            app.FlowrateLabel.Layout.Row = 3;
            app.FlowrateLabel.Layout.Column = 1;
            app.FlowrateLabel.Text = 'Flowrate (L/min)';
            app.FlowrateLabel.FontWeight = 'bold';
            
            app.FlowrateField = uieditfield(app.ParamsGrid, 'numeric');
            app.FlowrateField.Layout.Row = 3;
            app.FlowrateField.Layout.Column = 2;
            app.FlowrateField.Value = 0.18;
            app.FlowrateField.Limits = [0.01 2];
            
            % Inlet concentration field
            app.CinLabel = uilabel(app.ParamsGrid);
            app.CinLabel.Layout.Row = 4;
            app.CinLabel.Layout.Column = 1;
            app.CinLabel.Text = 'C_in (mol/L)';
            app.CinLabel.FontWeight = 'bold';
            
            app.CinField = uieditfield(app.ParamsGrid, 'numeric');
            app.CinField.Layout.Row = 4;
            app.CinField.Layout.Column = 2;
            app.CinField.Value = 0.05;
            app.CinField.Limits = [0.001 1];
            
            % Simulation time field
            app.SimTimeLabel = uilabel(app.ParamsGrid);
            app.SimTimeLabel.Layout.Row = 5;
            app.SimTimeLabel.Layout.Column = 1;
            app.SimTimeLabel.Text = 'Sim Time (min)';
            app.SimTimeLabel.FontWeight = 'bold';
            
            app.SimTimeField = uieditfield(app.ParamsGrid, 'numeric');
            app.SimTimeField.Layout.Row = 5;
            app.SimTimeField.Layout.Column = 2;
            app.SimTimeField.Value = 30;
            app.SimTimeField.Limits = [5 200];
            
            % Create actions section with 4 rows now
            app.ActionsPanel = uipanel(app.LeftGridLayout);
            app.ActionsPanel.Layout.Row = 4;
            app.ActionsPanel.Layout.Column = 1;
            app.ActionsPanel.Title = 'Simulation Actions';
            app.ActionsPanel.FontWeight = 'bold';
            
            app.ActionsGrid = uigridlayout(app.ActionsPanel);
         app.ActionsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
            app.ActionsGrid.ColumnWidth = {'1x'};
            app.ActionsGrid.Padding = [10 10 10 10];
            
            app.ManualKButton = uibutton(app.ActionsGrid, 'push');
            app.ManualKButton.Layout.Row = 1;
            app.ManualKButton.Layout.Column = 1;
            app.ManualKButton.Text = '🧮 Manual k (Relaxed Conv.)';
            app.ManualKButton.FontWeight = 'bold';
            app.ManualKButton.ButtonPushedFcn = createCallbackFcn(app, @ManualKButtonPushed, true);
            
            app.OptimizeButton = uibutton(app.ActionsGrid, 'push');
            app.OptimizeButton.Layout.Row = 2;
            app.OptimizeButton.Layout.Column = 1;
            app.OptimizeButton.Text = '🎯 Optimize (Standard Conv.)';
            app.OptimizeButton.FontWeight = 'bold';
            app.OptimizeButton.ButtonPushedFcn = createCallbackFcn(app, @OptimizeButtonPushed, true);
            
            app.OptimizeUltraButton = uibutton(app.ActionsGrid, 'push');
            app.OptimizeUltraButton.Layout.Row = 3;
            app.OptimizeUltraButton.Layout.Column = 1;
            app.OptimizeUltraButton.Text = '🏆 Ultra-Precision (Strict Conv.)';
            app.OptimizeUltraButton.FontWeight = 'bold';
            app.OptimizeUltraButton.ButtonPushedFcn = createCallbackFcn(app, @OptimizeUltraButtonPushed, true);
            
            % NEW CLEAR ALL BUTTON
            app.ClearAllButton = uibutton(app.ActionsGrid, 'push');
            app.ClearAllButton.Layout.Row = 4;
            app.ClearAllButton.Layout.Column = 1;
            app.ClearAllButton.Text = '🗑️ Clear All Data';
            app.ClearAllButton.FontWeight = 'bold';
            app.ClearAllButton.ButtonPushedFcn = createCallbackFcn(app, @ClearAllButtonPushed, true);
            
            
            % Create display options section
            app.DisplayPanel = uipanel(app.LeftGridLayout);
            app.DisplayPanel.Layout.Row = 5;
            app.DisplayPanel.Layout.Column = 1;
            app.DisplayPanel.Title = 'Display Options';
            app.DisplayPanel.FontWeight = 'bold';
            
            app.DisplayGrid = uigridlayout(app.DisplayPanel);
            app.DisplayGrid.RowHeight = {'fit', 'fit','fit'};
            app.DisplayGrid.ColumnWidth = {'1x', 'fit'};
            app.DisplayGrid.Padding = [10 10 10 10];
            
            app.JacketLabel = uilabel(app.DisplayGrid);
            app.JacketLabel.Layout.Row = 1;
            app.JacketLabel.Layout.Column = 1;
            app.JacketLabel.Text = 'Show Jacket Temperature';
            
            app.JacketSwitch = uiswitch(app.DisplayGrid, 'slider');
            app.JacketSwitch.Layout.Row = 1;
            app.JacketSwitch.Layout.Column = 2;
            app.JacketSwitch.ValueChangedFcn = createCallbackFcn(app, @JacketSwitchValueChanged, true);
            
            app.ProductsLabel = uilabel(app.DisplayGrid);
            app.ProductsLabel.Layout.Row = 2;
            app.ProductsLabel.Layout.Column = 1;
            app.ProductsLabel.Text = 'Show Products';
            
            app.ProductsSwitch = uiswitch(app.DisplayGrid, 'slider');
            app.ProductsSwitch.Layout.Row = 2;
            app.ProductsSwitch.Layout.Column = 2;
            app.ProductsSwitch.Value = 'On';
            app.ProductsSwitch.ValueChangedFcn = createCallbackFcn(app, @ProductsSwitchValueChanged, true);
            
            % Isothermal Control Toggle
app.IsothermalLabel = uilabel(app.DisplayGrid);
app.IsothermalLabel.Layout.Row = 3;
app.IsothermalLabel.Layout.Column = 1;
app.IsothermalLabel.Text = 'Isothermal Control';

app.IsothermalSwitch = uiswitch(app.DisplayGrid, 'slider');
app.IsothermalSwitch.Layout.Row = 3;
app.IsothermalSwitch.Layout.Column = 2;
app.IsothermalSwitch.Value = 'On';
app.IsothermalSwitch.ValueChangedFcn = createCallbackFcn(app, @IsothermalSwitchValueChanged, true);
            
            
            
            
            
            
            % Create export section
            app.ExportPanel = uipanel(app.LeftGridLayout);
            app.ExportPanel.Layout.Row = 6;
            app.ExportPanel.Layout.Column = 1;
            app.ExportPanel.Title = 'Export & Save';
            app.ExportPanel.FontWeight = 'bold';
            
            app.ExportGrid = uigridlayout(app.ExportPanel);
           app.ExportGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
            app.ExportGrid.ColumnWidth = {'1x'};
            app.ExportGrid.Padding = [10 10 10 10];
            
            app.SaveButton = uibutton(app.ExportGrid, 'push');
            app.SaveButton.Layout.Row = 1;
            app.SaveButton.Layout.Column = 1;
            app.SaveButton.Text = '💾 Save Project';
            app.SaveButton.FontWeight = 'bold';
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            
            app.ExportDataButton = uibutton(app.ExportGrid, 'push');
            app.ExportDataButton.Layout.Row = 2;
            app.ExportDataButton.Layout.Column = 1;
            app.ExportDataButton.Text = '📊 Export Data';
            app.ExportDataButton.FontWeight = 'bold';
            app.ExportDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportDataButtonPushed, true);
            
            app.ExportPlotsButton = uibutton(app.ExportGrid, 'push');
            app.ExportPlotsButton.Layout.Row = 3;
            app.ExportPlotsButton.Layout.Column = 1;
            app.ExportPlotsButton.Text = '🖼️ Export Plots';
            app.ExportPlotsButton.FontWeight = 'bold';
            app.ExportPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @ExportPlotsButtonPushed, true);
            
            % ADD AFTER app.ExportPlotsButton creation:
app.ExportEnhancedButton = uibutton(app.ExportGrid, 'push');
app.ExportEnhancedButton.Layout.Row = 4;
app.ExportEnhancedButton.Layout.Column = 1;
app.ExportEnhancedButton.Text = '🎯 Export Enhanced Analysis';
app.ExportEnhancedButton.FontWeight = 'bold';
app.ExportEnhancedButton.ButtonPushedFcn = createCallbackFcn(app, @ExportEnhancedButtonPushed, true);
            
            % Create status section
            app.StatusPanel = uipanel(app.LeftGridLayout);
            app.StatusPanel.Layout.Row = 7;
            app.StatusPanel.Layout.Column = 1;
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.FontWeight = 'bold';
            
            app.StatusGrid = uigridlayout(app.StatusPanel);
            app.StatusGrid.RowHeight = {'fit', 'fit'};
            app.StatusGrid.ColumnWidth = {'1x'};
            app.StatusGrid.Padding = [10 10 10 10];
            
            app.StatusLabel = uilabel(app.StatusGrid);
            app.StatusLabel.Layout.Row = 1;
            app.StatusLabel.Layout.Column = 1;
            app.StatusLabel.Text = 'Ready for Quick Sim and convergence analysis';
            app.StatusLabel.FontWeight = 'bold';
            
            app.ProgressLabel = uilabel(app.StatusGrid);
            app.ProgressLabel.Layout.Row = 2;
            app.ProgressLabel.Layout.Column = 1;
            app.ProgressLabel.Text = 'Progress: 0%';
            
            % Create right panel for visualization
            app.RightPanel = uipanel(app.MainGridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Title = '';
            
            app.RightGrid = uigridlayout(app.RightPanel);
            app.RightGrid.RowHeight = {'1x'};
            app.RightGrid.ColumnWidth = {'1x'};
            app.RightGrid.Padding = [10 10 10 10];
            
            % Create tab group
            app.TabGroup = uitabgroup(app.RightGrid);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;
            
            % Create concentration tab
            app.ConcentrationTab = uitab(app.TabGroup);
            app.ConcentrationTab.Title = 'Concentrations';
            
            app.ConcentrationGrid = uigridlayout(app.ConcentrationTab);
            app.ConcentrationGrid.RowHeight = {'1x', 'fit'};
            app.ConcentrationGrid.ColumnWidth = {'1x'};
            
            app.UIAxes1 = uiaxes(app.ConcentrationGrid);
            app.UIAxes1.Layout.Row = 1;
            app.UIAxes1.Layout.Column = 1;
            title(app.UIAxes1, 'Concentration Profiles');
            xlabel(app.UIAxes1, 'Time (min)');
            ylabel(app.UIAxes1, 'Concentration (mol/L)');
            
            app.ConcentrationInfo = uihtml(app.ConcentrationGrid);
            app.ConcentrationInfo.Layout.Row = 2;
            app.ConcentrationInfo.Layout.Column = 1;
            app.ConcentrationInfo.Position = [0 0 400 150];
            
            % Create temperature tab
            app.TemperatureTab = uitab(app.TabGroup);
            app.TemperatureTab.Title = 'Temperature';
            
            app.TemperatureGrid = uigridlayout(app.TemperatureTab);
            app.TemperatureGrid.RowHeight = {'1x', 'fit'};
            app.TemperatureGrid.ColumnWidth = {'1x'};
            
            app.UIAxes2 = uiaxes(app.TemperatureGrid);
            app.UIAxes2.Layout.Row = 1;
            app.UIAxes2.Layout.Column = 1;
            title(app.UIAxes2, 'Temperature Profile');
            xlabel(app.UIAxes2, 'Time (min)');
            ylabel(app.UIAxes2, 'Temperature (°C)');
            
            app.TemperatureInfo = uihtml(app.TemperatureGrid);
            app.TemperatureInfo.Layout.Row = 2;
            app.TemperatureInfo.Layout.Column = 1;
            app.TemperatureInfo.Position = [0 0 400 150];
            
            % Create products tab
            app.ProductsTab = uitab(app.TabGroup);
            app.ProductsTab.Title = 'Products';
            
            app.ProductsGrid = uigridlayout(app.ProductsTab);
            app.ProductsGrid.RowHeight = {'1x', 'fit'};
            app.ProductsGrid.ColumnWidth = {'1x'};
            
            app.UIAxes3 = uiaxes(app.ProductsGrid);
            app.UIAxes3.Layout.Row = 1;
            app.UIAxes3.Layout.Column = 1;
            title(app.UIAxes3, 'Product Formation');
            xlabel(app.UIAxes3, 'Time (min)');
            ylabel(app.UIAxes3, 'Concentration (mol/L)');
            
            app.ProductsInfo = uihtml(app.ProductsGrid);
            app.ProductsInfo.Layout.Row = 2;
            app.ProductsInfo.Layout.Column = 1;
            app.ProductsInfo.Position = [0 0 400 150];
            
            % Create analysis tab
            app.AnalysisTab = uitab(app.TabGroup);
            app.AnalysisTab.Title = 'Analysis';

            app.AnalysisGrid = uigridlayout(app.AnalysisTab);
            app.AnalysisGrid.RowHeight = {'1x', '1x'};
            app.AnalysisGrid.ColumnWidth = {'1x', '1x'};
            app.AnalysisGrid.RowSpacing = 10;
            app.AnalysisGrid.ColumnSpacing = 10;

            % Convergence panel with proper grid layout
            app.ConvergencePanel = uipanel(app.AnalysisGrid);
            app.ConvergencePanel.Layout.Row = 1;
            app.ConvergencePanel.Layout.Column = 1;
            app.ConvergencePanel.Title = 'Differential Convergence Analysis';
            app.ConvergencePanel.FontWeight = 'bold';

            % Create grid inside convergence panel
            convergenceGrid = uigridlayout(app.ConvergencePanel);
            convergenceGrid.RowHeight = {'1x'};
            convergenceGrid.ColumnWidth = {'1x'};
            convergenceGrid.Padding = [5 5 5 5];

            app.ConvergenceText = uitextarea(convergenceGrid);
            app.ConvergenceText.Layout.Row = 1;
            app.ConvergenceText.Layout.Column = 1;
            app.ConvergenceText.Editable = 'off';
            app.ConvergenceText.FontName = 'Courier';
            app.ConvergenceText.FontSize = 9;

           % Comparison table with proper titled panel
app.ComparisonTablePanel = uipanel(app.AnalysisGrid);
app.ComparisonTablePanel.Layout.Row = 1;
app.ComparisonTablePanel.Layout.Column = 2;
app.ComparisonTablePanel.Title = 'Experimental data vs Simulated data';
app.ComparisonTablePanel.FontWeight = 'bold';

% Create grid inside comparison table panel
comparisonGrid = uigridlayout(app.ComparisonTablePanel);
comparisonGrid.RowHeight = {'1x'};
comparisonGrid.ColumnWidth = {'1x'};
comparisonGrid.Padding = [5 5 5 5];

            app.ComparisonTable = uitable(comparisonGrid);
            app.ComparisonTable.Layout.Row = 1;
            app.ComparisonTable.Layout.Column = 1;
            app.ComparisonTable.ColumnName = {'Time', 'Experimental', 'Simulated', 'Error %'};
            app.ComparisonTable.RowName = {};

            % Concentration table panel with proper grid layout
            app.ConcentrationTablePanel = uipanel(app.AnalysisGrid);
            app.ConcentrationTablePanel.Layout.Row = 2;
            app.ConcentrationTablePanel.Layout.Column = 1;
            app.ConcentrationTablePanel.Title = 'Concentration Profile';
            app.ConcentrationTablePanel.FontWeight = 'bold';

            concentrationGrid = uigridlayout(app.ConcentrationTablePanel);
            concentrationGrid.RowHeight = {'1x'};
            concentrationGrid.ColumnWidth = {'1x'};
            concentrationGrid.Padding = [5 5 5 5];

            app.ConcentrationTable = uitable(concentrationGrid);
            app.ConcentrationTable.Layout.Row = 1;
            app.ConcentrationTable.Layout.Column = 1;
            app.ConcentrationTable.ColumnName = {'Time', 'NaOH', 'EtOAc', 'NaOAc', 'EtOH'};
            app.ConcentrationTable.RowName = {};

            % Temperature table panel with proper grid layout
            app.TemperatureTablePanel = uipanel(app.AnalysisGrid);
            app.TemperatureTablePanel.Layout.Row = 2;
            app.TemperatureTablePanel.Layout.Column = 2;
            app.TemperatureTablePanel.Title = 'Temperature Profile';
            app.TemperatureTablePanel.FontWeight = 'bold';

            temperatureGrid = uigridlayout(app.TemperatureTablePanel);
            temperatureGrid.RowHeight = {'1x'};
            temperatureGrid.ColumnWidth = {'1x'};
            temperatureGrid.Padding = [5 5 5 5];

            app.TemperatureTable = uitable(temperatureGrid);
            app.TemperatureTable.Layout.Row = 1;
            app.TemperatureTable.Layout.Column = 1;
            app.TemperatureTable.ColumnName = {'Time', 'Reactor', 'Jacket'};
            app.TemperatureTable.RowName = {};
            
            % NEW TAB 5: Quick Sim
            app.QuickSimTab = uitab(app.TabGroup);
            app.QuickSimTab.Title = 'Quick Sim';
            
            app.QuickSimGrid = uigridlayout(app.QuickSimTab);
            app.QuickSimGrid.RowHeight = {'1x'};
            app.QuickSimGrid.ColumnWidth = {'1x', '2x'};
            app.QuickSimGrid.ColumnSpacing = 10;
            app.QuickSimGrid.Padding = [10 10 10 10];
            
            % Left column: Input controls
            inputPanel = uipanel(app.QuickSimGrid);
            inputPanel.Layout.Row = 1;
            inputPanel.Layout.Column = 1;
            inputPanel.Title = 'Quick Simulation Parameters';
            inputPanel.FontWeight = 'bold';
            
            inputGrid = uigridlayout(inputPanel);
            inputGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            inputGrid.ColumnWidth = {'1x', '1x'};
            inputGrid.RowSpacing = 10;
            inputGrid.Padding = [10 10 10 10];
            
            % Cin field
            cinLabel = uilabel(inputGrid);
            cinLabel.Layout.Row = 1;
            cinLabel.Layout.Column = 1;
            cinLabel.Text = 'C_in (mol/L):';
            cinLabel.FontWeight = 'bold';
            
            app.CinQuickField = uieditfield(inputGrid, 'numeric');
            app.CinQuickField.Layout.Row = 1;
            app.CinQuickField.Layout.Column = 2;
            app.CinQuickField.Value = 0.05;
            app.CinQuickField.Limits = [0.001 1];
            
            % Flowrate field
            flowrateQuickLabel = uilabel(inputGrid);
            flowrateQuickLabel.Layout.Row = 2;
            flowrateQuickLabel.Layout.Column = 1;
            flowrateQuickLabel.Text = 'Flowrate (L/min):';
            flowrateQuickLabel.FontWeight = 'bold';
            
            app.FlowrateQuickField = uieditfield(inputGrid, 'numeric');
            app.FlowrateQuickField.Layout.Row = 2;
            app.FlowrateQuickField.Layout.Column = 2;
            app.FlowrateQuickField.Value = 0.18;
            app.FlowrateQuickField.Limits = [0.01 2];
            
            % Temperature dropdown
            tempQuickLabel = uilabel(inputGrid);
            tempQuickLabel.Layout.Row = 3;
            tempQuickLabel.Layout.Column = 1;
            tempQuickLabel.Text = 'Temperature (°C):';
            tempQuickLabel.FontWeight = 'bold';
            
            app.TempQuickDropDown = uidropdown(inputGrid);
            app.TempQuickDropDown.Layout.Row = 3;
            app.TempQuickDropDown.Layout.Column = 2;
            app.TempQuickDropDown.Items = {'30', '35', '40', '45', '50', '55', '60', '65', '70'};
            app.TempQuickDropDown.Value = '30';
            
            % Sim time field
            simTimeQuickLabel = uilabel(inputGrid);
            simTimeQuickLabel.Layout.Row = 4;
            simTimeQuickLabel.Layout.Column = 1;
            simTimeQuickLabel.Text = 'Sim Time (min):';
            simTimeQuickLabel.FontWeight = 'bold';
            
            app.SimTimeQuickField = uieditfield(inputGrid, 'numeric');
            app.SimTimeQuickField.Layout.Row = 4;
            app.SimTimeQuickField.Layout.Column = 2;
            app.SimTimeQuickField.Value = 30;
            app.SimTimeQuickField.Limits = [5 200];
            
            % Rate constant field
            rateConstLabel = uilabel(inputGrid);
            rateConstLabel.Layout.Row = 5;
            rateConstLabel.Layout.Column = 1;
            rateConstLabel.Text = 'Rate Const (L/mol·min):';
            rateConstLabel.FontWeight = 'bold';
            
            app.RateConstQuickField = uieditfield(inputGrid, 'numeric');
            app.RateConstQuickField.Layout.Row = 5;
            app.RateConstQuickField.Layout.Column = 2;
            app.RateConstQuickField.Value = 0.5;
            app.RateConstQuickField.Limits = [0.001 10];
            
            % Run button
            app.RunQuickSimButton = uibutton(inputGrid, 'push');
            app.RunQuickSimButton.Layout.Row = 6;
            app.RunQuickSimButton.Layout.Column = [1 2];
            app.RunQuickSimButton.Text = '🚀 Run Quick Sim';
            app.RunQuickSimButton.FontWeight = 'bold';
            app.RunQuickSimButton.ButtonPushedFcn = createCallbackFcn(app, @runQuickSim, true);
            
            % Error metrics
            errorPanel = uipanel(inputGrid);
            errorPanel.Layout.Row = 7;
            errorPanel.Layout.Column = [1 2];
            errorPanel.Title = 'Error Metrics vs Experimental';
            errorPanel.FontWeight = 'bold';
            
            errorGrid = uigridlayout(errorPanel);
            errorGrid.RowHeight = {'fit', 'fit', 'fit'};
            errorGrid.ColumnWidth = {'1x'};
            errorGrid.Padding = [5 5 5 5];
            
            app.MaxErrorLabel = uilabel(errorGrid);
            app.MaxErrorLabel.Layout.Row = 1;
            app.MaxErrorLabel.Layout.Column = 1;
            app.MaxErrorLabel.Text = 'Max Error: N/A';
            app.MaxErrorLabel.FontWeight = 'bold';
            
            app.MinErrorLabel = uilabel(errorGrid);
            app.MinErrorLabel.Layout.Row = 2;
            app.MinErrorLabel.Layout.Column = 1;
            app.MinErrorLabel.Text = 'Min Error: N/A';
            app.MinErrorLabel.FontWeight = 'bold';
            
            app.MeanErrorLabel = uilabel(errorGrid);
            app.MeanErrorLabel.Layout.Row = 3;
            app.MeanErrorLabel.Layout.Column = 1;
            app.MeanErrorLabel.Text = 'Mean Error: N/A';
            app.MeanErrorLabel.FontWeight = 'bold';
            
            % Info panel
            app.QuickSimInfo = uihtml(inputGrid);
            app.QuickSimInfo.Layout.Row = 8;
            app.QuickSimInfo.Layout.Column = [1 2];
            app.QuickSimInfo.Position = [0 0 400 120];
            % Right column: Plots
            plotPanel = uipanel(app.QuickSimGrid);
            plotPanel.Layout.Row = 1;
            plotPanel.Layout.Column = 2;
            plotPanel.Title = 'Quick Simulation Results with Convergence Analysis';
            plotPanel.FontWeight = 'bold';
            
            plotGrid = uigridlayout(plotPanel);
            plotGrid.RowHeight = {'1x', '1x', '1x'};
            plotGrid.ColumnWidth = {'1x'};
            plotGrid.RowSpacing = 10;
            plotGrid.Padding = [10 10 10 10];
            
            % Concentration plot
            app.UIAxesQuickConc = uiaxes(plotGrid);
            app.UIAxesQuickConc.Layout.Row = 1;
            app.UIAxesQuickConc.Layout.Column = 1;
            title(app.UIAxesQuickConc, 'Concentration Profiles with Convergence');
            xlabel(app.UIAxesQuickConc, 'Time (min)');
            ylabel(app.UIAxesQuickConc, 'Concentration (mol/L)');
            
            % Temperature plot
            app.UIAxesQuickTemp = uiaxes(plotGrid);
            app.UIAxesQuickTemp.Layout.Row = 2;
            app.UIAxesQuickTemp.Layout.Column = 1;
            title(app.UIAxesQuickTemp, 'Temperature Profiles');
            xlabel(app.UIAxesQuickTemp, 'Time (min)');
            ylabel(app.UIAxesQuickTemp, 'Temperature (°C)');
            
            % Products plot
            app.UIAxesQuickProducts = uiaxes(plotGrid);
            app.UIAxesQuickProducts.Layout.Row = 3;
            app.UIAxesQuickProducts.Layout.Column = 1;
            title(app.UIAxesQuickProducts, 'Product Formation with Convergence');
            xlabel(app.UIAxesQuickProducts, 'Time (min)');
            ylabel(app.UIAxesQuickProducts, 'Concentration (mol/L)');
            
            % Create context menus
            app.PlotContextMenu = uicontextmenu(app.UIFigure);
            
            app.CopyDataMenuItem = uimenu(app.PlotContextMenu);
            app.CopyDataMenuItem.Text = 'Copy Data to Clipboard';
            app.CopyDataMenuItem.MenuSelectedFcn = createCallbackFcn(app, @copyDataMenuSelected, true);
            
            app.SaveImageMenuItem = uimenu(app.PlotContextMenu);
            app.SaveImageMenuItem.Text = 'Save Image';
            app.SaveImageMenuItem.MenuSelectedFcn = createCallbackFcn(app, @saveImageMenuSelected, true);
            
            % Assign context menus to axes
            app.UIAxes1.ContextMenu = app.PlotContextMenu;
            app.UIAxes2.ContextMenu = app.PlotContextMenu;
            app.UIAxes3.ContextMenu = app.PlotContextMenu;
            app.UIAxesQuickConc.ContextMenu = app.PlotContextMenu;
            app.UIAxesQuickTemp.ContextMenu = app.PlotContextMenu;
            app.UIAxesQuickProducts.ContextMenu = app.PlotContextMenu;
    
            
            % NEW TAB 6: Conversion Analysis
app.ConversionTab = uitab(app.TabGroup);
app.ConversionTab.Title = 'Conversion';
app.ConversionGrid = uigridlayout(app.ConversionTab);
app.ConversionGrid.RowHeight = {'1x', '1x', 'fit'};
app.ConversionGrid.ColumnWidth = {'1x'};
app.ConversionGrid.RowSpacing = 10;
app.ConversionGrid.Padding = [10 10 10 10];

% Top graph: Conversion of reactant to product vs time
app.UIAxesConversionProduct = uiaxes(app.ConversionGrid);
app.UIAxesConversionProduct.Layout.Row = 1;
app.UIAxesConversionProduct.Layout.Column = 1;
title(app.UIAxesConversionProduct, 'Conversion: Reactant to Product Formation vs Time');
xlabel(app.UIAxesConversionProduct, 'Time (min)');
ylabel(app.UIAxesConversionProduct, 'Conversion (%)');

% Bottom graph: Conversion of unreacted reactant vs time  
app.UIAxesConversionReactant = uiaxes(app.ConversionGrid);
app.UIAxesConversionReactant.Layout.Row = 2;
app.UIAxesConversionReactant.Layout.Column = 1;
title(app.UIAxesConversionReactant, 'Conversion: Unreacted Reactant vs Time');
xlabel(app.UIAxesConversionReactant, 'Time (min)');
ylabel(app.UIAxesConversionReactant, 'Conversion (%)');
            
% Assign context menus to conversion axes
app.UIAxesConversionProduct.ContextMenu = app.PlotContextMenu;
app.UIAxesConversionReactant.ContextMenu = app.PlotContextMenu;
% ADD THIS MISSING SECTION: Conversion convergence info panel
app.ConversionInfo = uihtml(app.ConversionGrid);
app.ConversionInfo.Layout.Row = 3;
app.ConversionInfo.Layout.Column = 1;
app.ConversionInfo.Position = [0 0 400 120];
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
           % 1. PARAMETER UNCERTAINTY QUANTIFICATION
    function calculateParameterUncertainty(app)
        try
            if isempty(app.Data) || isnan(app.k_fit)
                app.updateStatus('No optimized parameters to analyze', 'warning');
                return;
            end
            
            app.showCancellableDialog('Uncertainty Analysis', 'Calculating parameter confidence intervals...');
            
            % Bootstrap analysis for parameter uncertainty
            nBootstrap = 500;
            t_exp = app.Data.Time;
            C_exp = app.Data.Concentration;
            n_data = length(C_exp);
            
            % Store bootstrap results
            k_bootstrap = zeros(nBootstrap, 1);
            V0_bootstrap = zeros(nBootstrap, 1);
            
            % Original parameters
            V0_original = app.FlowrateField.Value;
            k_original = app.k_fit;
            
            app.updateCancellableDialog(10, 'Running bootstrap analysis...');
            
            for i = 1:nBootstrap
                if app.OperationCancelled
                    return;
                end
                
                % Bootstrap resample
                boot_indices = randi(n_data, n_data, 1);
                t_boot = t_exp(boot_indices);
                C_boot = C_exp(boot_indices);
                
                % Sort by time
                [t_boot, sort_idx] = sort(t_boot);
                C_boot = C_boot(sort_idx);
                
                try
                    % Objective function for bootstrap sample
                    objectiveFunc = @(params) app.calculateMaxErrorForOptimization(params, t_boot, C_boot, ...
                        [app.Cin, app.Cin, 0, 0, app.SelectedTemp, app.Tj0], app.V, app.Cin, app.SelectedTemp);
                    
                    % Optimization bounds
                    lb = [0.01, 0.001];
                    ub = [1.0, 10];
                    
                    % Use original parameters as starting point with small perturbation
                    x0 = [V0_original + 0.01*randn(), k_original + 0.1*randn()];
                    x0 = max(min(x0, ub), lb);
                    
                    options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 100);
                    [x_opt, ~] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
                    
                    V0_bootstrap(i) = x_opt(1);
                    k_bootstrap(i) = x_opt(2);
                    
                catch
                    V0_bootstrap(i) = V0_original;
                    k_bootstrap(i) = k_original;
                end
                
                if mod(i, 50) == 0
                    app.updateCancellableDialog(10 + 70 * i / nBootstrap, ...
                        sprintf('Bootstrap progress: %d/%d', i, nBootstrap));
                end
            end
            
            app.updateCancellableDialog(80, 'Calculating confidence intervals...');
            
            % Calculate confidence intervals (95%)
            alpha = 0.05;
            k_ci = [prctile(k_bootstrap, 100*alpha/2), prctile(k_bootstrap, 100*(1-alpha/2))];
            V0_ci = [prctile(V0_bootstrap, 100*alpha/2), prctile(V0_bootstrap, 100*(1-alpha/2))];
            
            % Calculate standard errors
            k_std = std(k_bootstrap);
            V0_std = std(V0_bootstrap);
            
            % Store uncertainty results
            app.UncertaintyResults = struct();
            app.UncertaintyResults.k_bootstrap = k_bootstrap;
            app.UncertaintyResults.V0_bootstrap = V0_bootstrap;
            app.UncertaintyResults.k_ci = k_ci;
            app.UncertaintyResults.V0_ci = V0_ci;
            app.UncertaintyResults.k_std = k_std;
            app.UncertaintyResults.V0_std = V0_std;
            app.UncertaintyResults.k_original = k_original;
            app.UncertaintyResults.V0_original = V0_original;
            
            app.updateCancellableDialog(100, 'Uncertainty analysis complete!');
            app.closeCancellableDialog();
            
            % Display results
            msgText = sprintf(['📊 Parameter Uncertainty Analysis (95%% CI):\n\n' ...
                              '🎯 Rate Constant (k):\n' ...
                              '  • Value: %.6f ± %.6f L/mol·min\n' ...
                              '  • 95%% CI: [%.6f, %.6f]\n' ...
                              '  • Relative Error: ±%.2f%%\n\n' ...
                              '🌊 Flowrate (V₀):\n' ...
                              '  • Value: %.6f ± %.6f L/min\n' ...
                              '  • 95%% CI: [%.6f, %.6f]\n' ...
                              '  • Relative Error: ±%.2f%%'], ...
                              k_original, k_std, k_ci(1), k_ci(2), 100*k_std/k_original, ...
                              V0_original, V0_std, V0_ci(1), V0_ci(2), 100*V0_std/V0_original);
            
            uialert(app.UIFigure, msgText, 'Parameter Uncertainty Results', 'Icon', 'info');
            app.updateStatus('📊 Parameter uncertainty analysis completed', 'success');
            
        catch ME
            app.closeCancellableDialog();
            app.updateStatus(['Uncertainty analysis failed: ' ME.message], 'error');
        end
    end
    
    % 2. RESIDUAL ANALYSIS AND DIAGNOSTICS
    function performResidualsAnalysis(app)
        try
            if isempty(app.SimResults) || isempty(app.Data)
                app.updateStatus('Need simulation results and experimental data', 'warning');
                return;
            end
            
            % Calculate residuals
            t_exp = app.Data.Time;
            C_exp = app.Data.Concentration;
            C_sim = interp1(app.SimResults.t, app.SimResults.C_NaOH, t_exp, 'linear', 'extrap');
            
            residuals = C_exp - C_sim;
            standardized_residuals = residuals / std(residuals);
            
            % Create residuals analysis figure
            fig = figure('Name', 'Residuals Analysis', 'Position', [100 100 1200 800]);
            
            % 1. Residuals vs Fitted
            subplot(2, 3, 1);
            scatter(C_sim, residuals, 60, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
            xlabel('Fitted Values (mol/L)');
            ylabel('Residuals (mol/L)');
            title('Residuals vs Fitted');
            grid on;
            hold on;
            yline(0, '--k', 'LineWidth', 1.5);
            
            % Add trend line
            p = polyfit(C_sim, residuals, 1);
            x_trend = linspace(min(C_sim), max(C_sim), 100);
            y_trend = polyval(p, x_trend);
            plot(x_trend, y_trend, 'r--', 'LineWidth', 2);
            
            % 2. Q-Q Plot for normality
            subplot(2, 3, 2);
            qqplot(standardized_residuals);
            title('Q-Q Plot (Normality Test)');
            grid on;
            
            % 3. Residuals vs Time
            subplot(2, 3, 3);
            plot(t_exp, residuals, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
            xlabel('Time (min)');
            ylabel('Residuals (mol/L)');
            title('Residuals vs Time');
            grid on;
            hold on;
            yline(0, '--k', 'LineWidth', 1.5);
            
            % 4. Histogram of Residuals
            subplot(2, 3, 4);
            histogram(standardized_residuals, 'Normalization', 'pdf', 'FaceColor', [0.8500 0.3250 0.0980]);
            hold on;
            x_norm = linspace(-3, 3, 100);
            y_norm = normpdf(x_norm, 0, 1);
            plot(x_norm, y_norm, 'k-', 'LineWidth', 2);
            xlabel('Standardized Residuals');
            ylabel('Density');
            title('Residuals Distribution');
            legend('Residuals', 'Normal', 'Location', 'best');
            grid on;
            
            % 5. Leverage vs Residuals
            subplot(2, 3, 5);
            n = length(t_exp);
            leverage = 1/n + (t_exp - mean(t_exp)).^2 / sum((t_exp - mean(t_exp)).^2);
            scatter(leverage, abs(standardized_residuals), 60, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
            xlabel('Leverage');
            ylabel('|Standardized Residuals|');
            title('Leverage vs |Standardized Residuals|');
            grid on;
            
            % Add Cook's distance contours
            hold on;
            cooks_distance = 0.5;
            x_cook = linspace(0, max(leverage), 100);
            y_cook = sqrt(cooks_distance * 2 * (1 - x_cook) ./ x_cook);
            plot(x_cook, y_cook, 'r--', 'LineWidth', 1.5);
            text(max(leverage)*0.7, max(y_cook)*0.8, "Cook's d = 0.5", 'Color', 'red');
            
            % 6. Autocorrelation of Residuals
            subplot(2, 3, 6);
            [autocorr_vals, lags] = xcorr(residuals, min(10, floor(length(residuals)/4)), 'normalized');
            stem(lags, autocorr_vals, 'LineWidth', 1.5);
            xlabel('Lag');
            ylabel('Autocorrelation');
            title('Residuals Autocorrelation');
            grid on;
            hold on;
            % 95% confidence bounds
            conf_bound = 1.96/sqrt(length(residuals));
            yline(conf_bound, '--r', 'LineWidth', 1);
            yline(-conf_bound, '--r', 'LineWidth', 1);
            
            % Statistical tests
            [h_jb, p_jb] = jbtest(residuals);
            [h_bp, p_bp] = app.breuschPaganTest(C_sim, residuals);
            dw_stat = app.durbinWatsonTest(residuals);
            
            % Store diagnostics
            app.DiagnosticsResults = struct();
            app.DiagnosticsResults.residuals = residuals;
            app.DiagnosticsResults.standardized_residuals = standardized_residuals;
            app.DiagnosticsResults.leverage = leverage;
            app.DiagnosticsResults.jarque_bera_p = p_jb;
            app.DiagnosticsResults.breusch_pagan_p = p_bp;
            app.DiagnosticsResults.durbin_watson = dw_stat;
            app.DiagnosticsResults.rmse = sqrt(mean(residuals.^2));
            app.DiagnosticsResults.mae = mean(abs(residuals));
            
            % Add summary text
            sgtitle(sprintf(['Model Diagnostics Summary\n' ...
                           'RMSE: %.6f | MAE: %.6f | Jarque-Bera p: %.4f | Breusch-Pagan p: %.4f | Durbin-Watson: %.3f'], ...
                           app.DiagnosticsResults.rmse, app.DiagnosticsResults.mae, p_jb, p_bp, dw_stat), ...
                           'FontSize', 14, 'FontWeight', 'bold');
            
            app.updateStatus('📈 Residuals analysis completed - check new figure', 'success');
            
        catch ME
            app.updateStatus(['Residuals analysis failed: ' ME.message], 'error');
        end
    end
    
    % Helper function: Breusch-Pagan Test
    function [h, p] = breuschPaganTest(app, x, residuals)
        try
            n = length(residuals);
            squared_residuals = residuals.^2;
            
            % Regression of squared residuals on fitted values
            X = [ones(n,1), x(:)];
            beta = X \ squared_residuals;
            fitted_sq_res = X * beta;
            
            % Test statistic
            SSR = sum((fitted_sq_res - mean(squared_residuals)).^2);
            SST = sum((squared_residuals - mean(squared_residuals)).^2);
            R2 = SSR / SST;
            
            LM_stat = n * R2;
            p = 1 - chi2cdf(LM_stat, 1);
            h = p < 0.05;
        catch
            h = false;
            p = 1;
        end
    end
    
    % Helper function: Durbin-Watson Test
    function dw = durbinWatsonTest(app, residuals)
        try
            n = length(residuals);
            if n < 2
                dw = NaN;
                return;
            end
            
            diff_residuals = diff(residuals);
            dw = sum(diff_residuals.^2) / sum(residuals.^2);
        catch
            dw = NaN;
        end
    end
    
    % 3. ENHANCED DATA VALIDATION AND OUTLIER DETECTION
    function enhancedDataValidation(app)
        try
            if isempty(app.Data)
                return;
            end
            
            originalData = app.Data;
            t = app.Data.Time;
            C = app.Data.Concentration;
            n = length(t);
            
            app.showCancellableDialog('Data Validation', 'Analyzing data quality...');
            
            % 1. Basic validation
            validation_results = struct();
            validation_results.n_points = n;
            validation_results.time_span = max(t) - min(t);
            validation_results.concentration_range = [min(C), max(C)];
            
            app.updateCancellableDialog(20, 'Detecting outliers...');
            
            % 2. Outlier detection using multiple methods
            outliers_iqr = app.detectOutliersIQR(C);
            outliers_zscore = abs(zscore(C)) > 3;
            outliers_modified_zscore = app.detectOutliersModifiedZScore(C);
            
            % Combine outlier detection methods
            outliers_combined = outliers_iqr | outliers_zscore | outliers_modified_zscore;
            n_outliers = sum(outliers_combined);
            
            validation_results.outliers_indices = find(outliers_combined);
            validation_results.n_outliers = n_outliers;
            validation_results.outlier_percentage = 100 * n_outliers / n;
            
            app.updateCancellableDialog(40, 'Analyzing data smoothness...');
            
            % 3. Data smoothness analysis
            if n > 2
                % Calculate derivatives (rate of change)
                dt = diff(t);
                dC = diff(C);
                dC_dt = dC ./ dt;
                
                % Detect sudden changes
                dC_dt_threshold = 3 * std(dC_dt);
                sudden_changes = abs(dC_dt) > dC_dt_threshold;
                
                validation_results.sudden_changes_indices = find(sudden_changes) + 1;
                validation_results.n_sudden_changes = sum(sudden_changes);
                validation_results.max_rate_change = max(abs(dC_dt));
            end
            
            app.updateCancellableDialog(60, 'Checking data gaps...');
            
            % 4. Data gaps analysis
            median_dt = median(diff(t));
            large_gaps = diff(t) > 3 * median_dt;
            validation_results.large_gaps_indices = find(large_gaps) + 1;
            validation_results.n_large_gaps = sum(large_gaps);
            
            % 5. Signal-to-noise ratio estimation
            if n > 10
                try
                    window_size = min(5, floor(n/3));
                    if window_size >= 3
                        C_smooth = movmean(C, window_size);
                        noise = C - C_smooth;
                        signal_power = var(C_smooth);
                        noise_power = var(noise);
                        snr_db = 10 * log10(signal_power / noise_power);
                        validation_results.snr_db = snr_db;
                    end
                catch
                    validation_results.snr_db = NaN;
                end
            end
            
            app.updateCancellableDialog(80, 'Generating validation report...');
            
            % 6. Generate validation report
            if n_outliers > 0 || validation_results.n_sudden_changes > 0 || validation_results.n_large_gaps > 0
                % Create validation figure
                fig = figure('Name', 'Data Quality Analysis', 'Position', [100 100 1200 600]);
                
                subplot(1, 2, 1);
                plot(t, C, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
                hold on;
                
                if n_outliers > 0
                    scatter(t(outliers_combined), C(outliers_combined), 100, 'red', 'filled', 'MarkerEdgeColor', 'black');
                end
                
                xlabel('Time (min)');
                ylabel('Concentration (mol/L)');
                title('Data Quality Overview');
                grid on;
                legend('Data', 'Outliers', 'Location', 'best');
                
                subplot(1, 2, 2);
                if n > 2
                    plot(t(2:end), dC_dt, 'g-', 'LineWidth', 1.5);
                    hold on;
                    if validation_results.n_sudden_changes > 0
                        scatter(t(validation_results.sudden_changes_indices), ...
                               dC_dt(validation_results.sudden_changes_indices-1), ...
                               100, 'red', 'filled');
                    end
                    xlabel('Time (min)');
                    ylabel('Rate of Change (mol/L/min)');
                    title('Rate of Change Analysis');
                    grid on;
                end
            end
            
            % Store validation results
            app.ValidationResults = validation_results;
            
            app.updateCancellableDialog(100, 'Data validation complete!');
            app.closeCancellableDialog();
            
            % Display validation summary
            msgText = sprintf(['📋 Data Quality Report:\n\n' ...
                              '📊 Dataset: %d points over %.1f minutes\n' ...
                              '🎯 Concentration range: [%.6f, %.6f] mol/L\n\n' ...
                              '⚠️ Quality Issues:\n' ...
                              '  • Outliers: %d (%.1f%%)\n' ...
                              '  • Sudden changes: %d\n' ...
                              '  • Large time gaps: %d\n\n' ...
                              '📡 Signal quality:\n' ...
                              '  • SNR: %.1f dB'], ...
                              n, validation_results.time_span, ...
                              validation_results.concentration_range(1), validation_results.concentration_range(2), ...
                              n_outliers, validation_results.outlier_percentage, ...
                              validation_results.n_sudden_changes, validation_results.n_large_gaps, ...
                              validation_results.snr_db);
            
            if n_outliers > 0
                selection = uiconfirm(app.UIFigure, [msgText '\n\nWould you like to remove detected outliers?'], ...
                                     'Data Quality Report', ...
                                     'Options', {'Remove Outliers', 'Keep All Data'}, ...
                                     'DefaultOption', 2, 'Icon', 'warning');
                
                if strcmp(selection, 'Remove Outliers')
                    % Remove outliers
                    clean_indices = ~outliers_combined;
                    app.Data = app.Data(clean_indices, :);
                    app.updateStatus(sprintf('🧹 Removed %d outliers from dataset', n_outliers), 'success');
                    app.FileLabel.Text = sprintf('%s (%d points, %d removed)', ...
                                                app.FileLabel.Text, sum(clean_indices), n_outliers);
                end
            else
                uialert(app.UIFigure, msgText, 'Data Quality Report', 'Icon', 'info');
            end
            
        catch ME
            app.closeCancellableDialog();
            app.updateStatus(['Data validation failed: ' ME.message], 'error');
        end
    end
    
    % Helper: IQR outlier detection
    function outliers = detectOutliersIQR(app, data)
        Q1 = prctile(data, 25);
        Q3 = prctile(data, 75);
        IQR = Q3 - Q1;
        lower_bound = Q1 - 1.5 * IQR;
        upper_bound = Q3 + 1.5 * IQR;
        outliers = data < lower_bound | data > upper_bound;
    end
    
    % Helper: Modified Z-score outlier detection
    function outliers = detectOutliersModifiedZScore(app, data)
        median_data = median(data);
        mad_data = median(abs(data - median_data));
        modified_z_scores = 0.6745 * (data - median_data) / mad_data;
        outliers = abs(modified_z_scores) > 3.5;
    end
    
    % 4. ENHANCED PROGRESS TRACKING WITH CANCELLATION
    function showCancellableDialog(app, title, message)
        try
            % Create a custom dialog with cancellation capability
            app.CancellableDialog = uifigure('Name', title, 'Position', [500 400 400 200], ...
                                           'WindowStyle', 'modal', 'Resize', 'off');
            
            grid = uigridlayout(app.CancellableDialog, [4 1]);
            grid.Padding = [20 20 20 20];
            grid.RowHeight = {'fit', '1x', 'fit', 'fit'};
            
            % Title
            titleLabel = uilabel(grid);
            titleLabel.Layout.Row = 1;
            titleLabel.Text = title;
            titleLabel.FontSize = 16;
            titleLabel.FontWeight = 'bold';
            titleLabel.HorizontalAlignment = 'center';
            
            % Message
            app.CancellableMessage = uilabel(grid);
            app.CancellableMessage.Layout.Row = 2;
            app.CancellableMessage.Text = message;
            app.CancellableMessage.HorizontalAlignment = 'center';
            app.CancellableMessage.WordWrap = 'on';
            
            % Progress bar
            app.CancellableProgress = uislider(grid);
            app.CancellableProgress.Layout.Row = 3;
            app.CancellableProgress.Limits = [0 100];
            app.CancellableProgress.Value = 0;
            app.CancellableProgress.Enable = 'off';
            
            % Cancel button
            app.CancelButton = uibutton(grid, 'push');
            app.CancelButton.Layout.Row = 4;
            app.CancelButton.Text = 'Cancel';
            app.CancelButton.ButtonPushedFcn = @(~,~) app.cancelOperation();
            
            % Initialize cancellation flag
            app.OperationCancelled = false;
            
        catch ME
            app.updateStatus(['Failed to create cancellable dialog: ' ME.message], 'error');
        end
    end
    
    function updateCancellableDialog(app, progress, message)
        try
            if isfield(app, 'CancellableDialog') && isvalid(app.CancellableDialog) && ~app.OperationCancelled
                app.CancellableProgress.Value = progress;
                if nargin > 2
                    app.CancellableMessage.Text = message;
                end
                drawnow;
            end
        catch
            % Silent fail
        end
    end
    
    function cancelOperation(app)
        app.OperationCancelled = true;
        app.closeCancellableDialog();
        app.updateStatus('🛑 Operation cancelled by user', 'warning');
    end
    
    function closeCancellableDialog(app)
        try
            if isfield(app, 'CancellableDialog') && isvalid(app.CancellableDialog)
                close(app.CancellableDialog);
                app.CancellableDialog = [];
            end
        catch
            % Silent fail
        end
    end
    
    % 5. ENHANCED EXPORT WITH UNCERTAINTY AND DIAGNOSTICS
    function exportEnhancedResults(app)
        try
            if isempty(app.SimResults) && isempty(app.ManualSimResults)
                app.updateStatus('No simulation results to export!', 'error');
                return;
            end
            
            [file, path] = uiputfile({'*.xlsx', 'Excel Files (*.xlsx)'}, ...
                                    'Export Enhanced Results', 'CSTR_Enhanced_Analysis.xlsx');
            if isequal(file, 0)
                return;
            end
            
            fullPath = fullfile(path, file);
            app.showCancellableDialog('Enhanced Export', 'Preparing comprehensive enhanced export...');
            
            % Delete existing file
            if exist(fullPath, 'file')
                delete(fullPath);
            end
            
            % Export with enhanced data
            app.exportDataWithUncertainty(fullPath);
            
            app.closeCancellableDialog();
            app.updateStatus('📊 Enhanced results exported with uncertainty analysis!', 'success');
            
        catch ME
            app.closeCancellableDialog();
            app.updateStatus(['Enhanced export failed: ' ME.message], 'error');
        end
    end
    
    function exportDataWithUncertainty(app, fullPath)
        try
            % 1. Enhanced parameters
            app.updateCancellableDialog(10, 'Exporting enhanced parameters...');
            
            paramNames = {'Reactor Volume (L)'; 'Total Flowrate (L/min)'; 'EtOAc Flowrate (L/min)'; 
                         'NaOH Flowrate (L/min)'; 'Inlet Concentration (mol/L)'; 'Temperature (°C)';
                         'Simulation Time (min)'; 'Manual Rate Constant (L/mol·min)'; 'Optimized Rate Constant (L/mol·min)';
                         'Manual R-squared'; 'Optimized R-squared'};
            
            paramValues = {app.V; app.V0; app.V0_EtOAc; app.V1_NaOH; app.Cin; app.SelectedTemp;
                          app.SimTime; app.k_manual; app.k_fit; app.RSquaredManual; app.RSquared};
            
            % Add uncertainty information if available
            if isfield(app, 'UncertaintyResults') && ~isempty(app.UncertaintyResults)
                paramNames = [paramNames; {'Rate Constant Std Error (L/mol·min)'; 'Rate Constant 95% CI Lower'; 
                             'Rate Constant 95% CI Upper'; 'Flowrate Std Error (L/min)'; 'Flowrate 95% CI Lower'; 
                             'Flowrate 95% CI Upper'}];
                paramValues = [paramValues; {app.UncertaintyResults.k_std; app.UncertaintyResults.k_ci(1); 
                              app.UncertaintyResults.k_ci(2); app.UncertaintyResults.V0_std; 
                              app.UncertaintyResults.V0_ci(1); app.UncertaintyResults.V0_ci(2)}];
            end
            
            parametersTable = table(paramNames, paramValues, 'VariableNames', {'Parameter', 'Value'});
            writetable(parametersTable, fullPath, 'Sheet', '1_Enhanced_Parameters', 'WriteMode', 'overwritesheet');
            
            % 2. Diagnostics results
            app.updateCancellableDialog(30, 'Exporting diagnostics...');
            
            if isfield(app, 'DiagnosticsResults') && ~isempty(app.DiagnosticsResults)
                diagNames = {'RMSE'; 'MAE'; 'Jarque-Bera p-value'; 'Breusch-Pagan p-value'; 'Durbin-Watson Statistic'};
                diagValues = {app.DiagnosticsResults.rmse; app.DiagnosticsResults.mae; 
                             app.DiagnosticsResults.jarque_bera_p; app.DiagnosticsResults.breusch_pagan_p; 
                             app.DiagnosticsResults.durbin_watson};
                
                diagnosticsTable = table(diagNames, diagValues, 'VariableNames', {'Diagnostic', 'Value'});
                writetable(diagnosticsTable, fullPath, 'Sheet', '2_Model_Diagnostics', 'WriteMode', 'overwritesheet');
                
                % Detailed residuals
                if ~isempty(app.Data)
                    t_exp = app.Data.Time;
                    C_exp = app.Data.Concentration;
                    residualsTable = table(t_exp, C_exp, app.DiagnosticsResults.residuals, ...
                                          app.DiagnosticsResults.standardized_residuals, app.DiagnosticsResults.leverage, ...
                                          'VariableNames', {'Time_min', 'Experimental_molL', 'Residuals_molL', ...
                                          'Standardized_Residuals', 'Leverage'});
                    writetable(residualsTable, fullPath, 'Sheet', '3_Residuals_Analysis', 'WriteMode', 'overwritesheet');
                end
            end
            
            % 3. Data validation results
            app.updateCancellableDialog(50, 'Exporting data validation...');
            
            if isfield(app, 'ValidationResults') && ~isempty(app.ValidationResults)
                valNames = {'Number of Points'; 'Time Span (min)'; 'Concentration Range Min'; 'Concentration Range Max';
                           'Number of Outliers'; 'Outlier Percentage'; 'Number of Sudden Changes'; 'Number of Large Gaps';
                           'SNR (dB)'};
                valValues = {app.ValidationResults.n_points; app.ValidationResults.time_span; 
                            app.ValidationResults.concentration_range(1); app.ValidationResults.concentration_range(2);
                            app.ValidationResults.n_outliers; app.ValidationResults.outlier_percentage;
                            app.ValidationResults.n_sudden_changes; app.ValidationResults.n_large_gaps;
                            app.ValidationResults.snr_db};
                
                validationTable = table(valNames, valValues, 'VariableNames', {'Metric', 'Value'});
                writetable(validationTable, fullPath, 'Sheet', '4_Data_Quality', 'WriteMode', 'overwritesheet');
            end
            
            % 4. Bootstrap results (if available)
            app.updateCancellableDialog(70, 'Exporting uncertainty analysis...');
            
            if isfield(app, 'UncertaintyResults') && ~isempty(app.UncertaintyResults)
                bootstrapTable = table(app.UncertaintyResults.k_bootstrap, app.UncertaintyResults.V0_bootstrap, ...
                                      'VariableNames', {'k_bootstrap_LmolMin', 'V0_bootstrap_Lmin'});
                writetable(bootstrapTable, fullPath, 'Sheet', '5_Bootstrap_Results', 'WriteMode', 'overwritesheet');
            end
            
            % 5. Continue with regular export sheets...
            app.updateCancellableDialog(90, 'Exporting simulation data...');
            
            % Export simulation results
            if ~isempty(app.SimResults)
                simData = table(app.SimResults.t, app.SimResults.C_NaOH, app.SimResults.C_EtOAc, ...
                               app.SimResults.C_NaOAc, app.SimResults.C_EtOH, app.SimResults.T_reactor, ...
                               'VariableNames', {'Time_min', 'C_NaOH_molL', 'C_EtOAc_molL', ...
                               'C_NaOAc_molL', 'C_EtOH_molL', 'T_reactor_C'});
                writetable(simData, fullPath, 'Sheet', '6_Simulation_Results', 'WriteMode', 'overwritesheet');
            end
            
            app.updateCancellableDialog(100, 'Export complete!');
            
        catch ME
            error('Enhanced export failed: %s', ME.message);
        end
    end
    end
    
    % Component initialization callbacks
    methods (Access = private)
        
  function exportPlots(app)
    try
        if isempty(app.SimResults) && isempty(app.ManualSimResults) && isempty(app.QuickSimResults)
            app.updateStatus('No plots to export!', 'error');
            return;
        end
        
        app.showLoadingDialog('Exporting Plots', 'Preparing comprehensive plot export...');
        
        folder = uigetdir(pwd, 'Select Folder for Plot Images');
        if isequal(folder, 0)
            app.closeLoadingDialog();
            return;
        end
        
        % Create timestamp for unique filenames
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        exportedCount = 0;
        
        app.updateLoadingDialog(0.1, 'Exporting concentration plot...');
        
        % 1. Export main concentration plot
        try
            fig1 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
            ax1 = axes(fig1, 'Position', [0.1, 0.15, 0.8, 0.75]);
            
            % Recreate the concentration plot
            hold(ax1, 'on');
            
            % Plot experimental data if available
            if ~isempty(app.Data)
                scatter(ax1, app.Data.Time, app.Data.Concentration, 80, ...
                       'MarkerEdgeColor', app.ErrorColor, ...
                       'MarkerFaceColor', app.BackgroundColor, ...
                       'LineWidth', 2.5, 'DisplayName', 'Experimental Data');
            end
            
            % Plot manual simulation if available
            if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                plot(ax1, app.ManualSimResults.t, app.ManualSimResults.C_NaOH, ':', ...
                     'LineWidth', 3, 'Color', [147 112 219]/255, 'DisplayName', 'Manual k Simulation');
            end
            
            % Plot optimized simulation if available
            if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                plot(ax1, app.SimResults.t, app.SimResults.C_NaOH, '-', ...
                     'LineWidth', 3, 'Color', app.AccentColor, 'DisplayName', 'Optimized Simulation');
            end
            
        % Plot ultra-precision simulation if available
            if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                plot(ax1, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, '-', ...
                     'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision Simulation');
                
                % *** ADD ULTRA-PRECISION CONVERGENCE MARKERS TO EXPORT ***
                if ~isempty(app.UltraConvergenceData)
                    app.addConvergenceMarkers(ax1, app.UltraPrecisionResults, app.UltraConvergenceData, [255 20 147]/255, 'Ultra-Precision');
                end
            end

            
            % Add title with error information
            if ~isnan(app.k_fit) && app.k_fit > 0
                if ~isempty(app.ComparisonTableData)
                    max_error = max(app.ComparisonTableData.Optimized_Error_percent);
                    titleStr = sprintf('CSTR Concentration Profiles (k_{opt} = %.4f L/mol·min, Max Error: %.3f%%)', app.k_fit, max_error);
                else
                    titleStr = sprintf('CSTR Concentration Profiles (k_{opt} = %.4f L/mol·min)', app.k_fit);
                end
            else
                titleStr = 'CSTR Concentration Profiles with Differential Convergence';
            end
            
            title(ax1, titleStr, 'FontWeight', 'bold', 'FontSize', 14);
            xlabel(ax1, 'Time (min)', 'FontWeight', 'bold');
            ylabel(ax1, 'NaOH Concentration (mol/L)', 'FontWeight', 'bold');
            legend(ax1, 'Location', 'best', 'Box', 'on');
            grid(ax1, 'off'); % No grid
            
            % Save the plot
            filename = sprintf('Concentration_Profiles_%s.png', timestamp);
            if exist('exportgraphics', 'file')
                exportgraphics(ax1, fullfile(folder, filename), 'Resolution', 300);
            else
                print(fig1, fullfile(folder, filename), '-dpng', '-r300');
            end
            close(fig1);
            exportedCount = exportedCount + 1;
        catch ME
            fprintf('Concentration plot export failed: %s\n', ME.message);
        end
        
        app.updateLoadingDialog(0.2, 'Exporting temperature plot...');
        
        % 2. Export temperature plot
        if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
            try
                fig2 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax2 = axes(fig2, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax2, 'on');
                
                % Plot manual simulation temperatures if available
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 'T_reactor')
                    plot(ax2, app.ManualSimResults.t, app.ManualSimResults.T_reactor, ':', ...
                         'LineWidth', 3, 'Color', [255 102 102]/255, 'DisplayName', 'Manual Reactor T');
                    if app.ShowJacket
                        plot(ax2, app.ManualSimResults.t, app.ManualSimResults.T_jacket, ':', ...
                             'LineWidth', 2.5, 'Color', [102 178 255]/255, 'DisplayName', 'Manual Jacket T');
                    end
                end
                
                % Plot optimized temperature data
                if ~isempty(app.SimResults)
                    plot(ax2, app.SimResults.t, app.SimResults.T_reactor, '-', ...
                         'LineWidth', 3, 'Color', app.ErrorColor, 'DisplayName', 'Reactor Temperature');
                    
                    if app.ShowJacket
                        plot(ax2, app.SimResults.t, app.SimResults.T_jacket, '-', ...
                             'LineWidth', 2.5, 'Color', app.AccentColor, 'DisplayName', 'Jacket Temperature');
                    end
                end
                
                % Plot ultra-precision temperatures if available
                if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 'T_reactor')
                    plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.T_reactor, '-', ...
                         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision Reactor T');
                    if app.ShowJacket
                        plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.T_jacket, '-', ...
                             'LineWidth', 3.5, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra-Precision Jacket T');
                    end
                end
                
                title(ax2, 'Enhanced Temperature Profile with Heat Transfer', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax2, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax2, 'Temperature (°C)', 'FontWeight', 'bold');
                legend(ax2, 'Location', 'best', 'Box', 'on');
                grid(ax2, 'off'); % No grid
                
                filename = sprintf('Temperature_Profiles_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax2, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig2, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig2);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Temperature plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.3, 'Exporting products plot...');
        
        % 3. Export products plot
        if app.ShowProducts && (~isempty(app.SimResults) || ~isempty(app.ManualSimResults))
            try
                fig3 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax3 = axes(fig3, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax3, 'on');
                
                % Plot manual products
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                    plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_NaOAc, ':', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Manual NaOAc');
                    plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_EtOH, ':', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Manual EtOH');
                end
                
                % Plot optimized products
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    plot(ax3, app.SimResults.t, app.SimResults.C_NaOAc, '-', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Optimized NaOAc');
                    plot(ax3, app.SimResults.t, app.SimResults.C_EtOH, '-', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Optimized EtOH');
                end
                
                % Plot ultra-precision products
                if ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't')
                    plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
                         'LineWidth', 4, 'Color', [255 20 147]/255, 'DisplayName', 'Ultra-Precision NaOAc');
                    plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
                         'LineWidth', 4, 'Color', [139 0 139]/255, 'DisplayName', 'Ultra-Precision EtOH');
                end
                
                title(ax3, 'Product Formation with Differential Convergence Analysis', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax3, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax3, 'Concentration (mol/L)', 'FontWeight', 'bold');
                legend(ax3, 'Location', 'best', 'Box', 'on');
                grid(ax3, 'off'); % No grid
                
                filename = sprintf('Product_Formation_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax3, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig3, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig3);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Products plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.4, 'Exporting Quick Sim concentration plot...');
        
        % 4. Export Quick Sim concentration plot
        if ~isempty(app.QuickSimResults)
            try
                fig4 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax4 = axes(fig4, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax4, 'on');
                
                plot(ax4, app.QuickSimResults.t, app.QuickSimResults.C_NaOH, '-', ...
                     'LineWidth', 2, 'Color', [255 59 48]/255, 'DisplayName', 'NaOH');
                plot(ax4, app.QuickSimResults.t, app.QuickSimResults.C_EtOAc, '--', ...
                     'LineWidth', 2, 'Color', [255 149 0]/255, 'DisplayName', 'EtOAc');
                plot(ax4, app.QuickSimResults.t, app.QuickSimResults.C_NaOAc, '-.', ...
                     'LineWidth', 2, 'Color', [255 204 0]/255, 'DisplayName', 'NaOAc');
                plot(ax4, app.QuickSimResults.t, app.QuickSimResults.C_EtOH, ':', ...
                     'LineWidth', 2, 'Color', [88 86 214]/255, 'DisplayName', 'EtOH');
                
                title(ax4, 'Quick Sim Concentration Profiles with Convergence Analysis', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax4, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax4, 'Concentration (mol/L)', 'FontWeight', 'bold');
                legend(ax4, 'Location', 'best', 'Box', 'on');
                grid(ax4, 'off'); % No grid
                
                filename = sprintf('QuickSim_Concentrations_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax4, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig4, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig4);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Quick Sim concentration plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.5, 'Exporting Quick Sim temperature plot...');
        
        % 5. Export Quick Sim temperature plot
        if ~isempty(app.QuickSimResults)
            try
                fig5 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax5 = axes(fig5, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax5, 'on');
                plot(ax5, app.QuickSimResults.t, app.QuickSimResults.T_reactor, '-', ...
                     'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Reactor T');
                plot(ax5, app.QuickSimResults.t, app.QuickSimResults.T_jacket, '--', ...
                     'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Jacket T');
                
                title(ax5, 'Quick Sim Temperature Profiles', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax5, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax5, 'Temperature (°C)', 'FontWeight', 'bold');
                legend(ax5, 'Location', 'best', 'Box', 'on');
                grid(ax5, 'off'); % No grid
                
                filename = sprintf('QuickSim_Temperature_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax5, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig5, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig5);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Quick Sim temperature plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.6, 'Exporting Quick Sim products plot...');
        
        % 6. Export Quick Sim products plot
        if ~isempty(app.QuickSimResults)
            try
                fig6 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax6 = axes(fig6, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax6, 'on');
                plot(ax6, app.QuickSimResults.t, app.QuickSimResults.C_NaOAc, '-', ...
                     'LineWidth', 2, 'Color', [255 204 0]/255, 'DisplayName', 'NaOAc');
                plot(ax6, app.QuickSimResults.t, app.QuickSimResults.C_EtOH, '--', ...
                     'LineWidth', 2, 'Color', [88 86 214]/255, 'DisplayName', 'EtOH');
                
                title(ax6, 'Quick Sim Product Formation with Convergence Analysis', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax6, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax6, 'Concentration (mol/L)', 'FontWeight', 'bold');
                legend(ax6, 'Location', 'best', 'Box', 'on');
                grid(ax6, 'off'); % No grid
                
                filename = sprintf('QuickSim_Products_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax6, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig6, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig6);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Quick Sim products plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.7, 'Exporting conversion plots...');
        
        % 7. Export conversion product formation plot
        if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults) || ~isempty(app.QuickSimResults)
            try
                fig7 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax7 = axes(fig7, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax7, 'on');
                
                % Manual simulation conversion
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                    conv_NaOAc_manual = (app.ManualSimResults.C_NaOAc / app.Cin) * 100;
                    conv_EtOH_manual = (app.ManualSimResults.C_EtOH / app.Cin) * 100;
                    
                    plot(ax7, app.ManualSimResults.t, conv_NaOAc_manual, ':', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Manual NaOAc Conversion');
                    plot(ax7, app.ManualSimResults.t, conv_EtOH_manual, ':', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Manual EtOH Conversion');
                end
                
                % Optimized simulation conversion
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    conv_NaOAc_opt = (app.SimResults.C_NaOAc / app.Cin) * 100;
                    conv_EtOH_opt = (app.SimResults.C_EtOH / app.Cin) * 100;
                    
                    plot(ax7, app.SimResults.t, conv_NaOAc_opt, '-', ...
                         'LineWidth', 3, 'Color', [255 204 0]/255, 'DisplayName', 'Optimized NaOAc Conversion');
                    plot(ax7, app.SimResults.t, conv_EtOH_opt, '-', ...
                         'LineWidth', 3, 'Color', [88 86 214]/255, 'DisplayName', 'Optimized EtOH Conversion');
                end
                
                % Quick simulation conversion
                if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                    conv_NaOAc_quick = (app.QuickSimResults.C_NaOAc / app.Cin) * 100;
                    conv_EtOH_quick = (app.QuickSimResults.C_EtOH / app.Cin) * 100;
                    
                    plot(ax7, app.QuickSimResults.t, conv_NaOAc_quick, '--', ...
                         'LineWidth', 2.5, 'Color', [255 165 0]/255, 'DisplayName', 'Quick Sim NaOAc Conversion');
                    plot(ax7, app.QuickSimResults.t, conv_EtOH_quick, '--', ...
                         'LineWidth', 2.5, 'Color', [30 144 255]/255, 'DisplayName', 'Quick Sim EtOH Conversion');
                end
                
                title(ax7, 'Conversion: Reactant to Product Formation vs Time', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax7, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax7, 'Conversion (%)', 'FontWeight', 'bold');
                legend(ax7, 'Location', 'best', 'Box', 'on');
                grid(ax7, 'off'); % No grid
                ylim(ax7, [0 100]);
                
                filename = sprintf('Conversion_Product_Formation_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax7, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig7, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig7);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Conversion product plot export failed: %s\n', ME.message);
            end
        end
        
        app.updateLoadingDialog(0.9, 'Exporting unreacted reactants plot...');
        
        % 8. Export conversion unreacted reactants plot
        if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults) || ~isempty(app.QuickSimResults)
            try
                fig8 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
                ax8 = axes(fig8, 'Position', [0.1, 0.15, 0.8, 0.75]);
                
                hold(ax8, 'on');
                
                % Manual simulation unreacted conversion
                if ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't')
                    conv_NaOH_manual = ((app.Cin - app.ManualSimResults.C_NaOH) / app.Cin) * 100;
                    conv_EtOAc_manual = ((app.Cin - app.ManualSimResults.C_EtOAc) / app.Cin) * 100;
                    
                    plot(ax8, app.ManualSimResults.t, conv_NaOH_manual, ':', ...
                         'LineWidth', 3, 'Color', [255 59 48]/255, 'DisplayName', 'Manual NaOH Conversion');
                    plot(ax8, app.ManualSimResults.t, conv_EtOAc_manual, ':', ...
                         'LineWidth', 3, 'Color', [255 149 0]/255, 'DisplayName', 'Manual EtOAc Conversion');
                end
                
                % Optimized simulation unreacted conversion
                if ~isempty(app.SimResults) && isfield(app.SimResults, 't')
                    conv_NaOH_opt = ((app.Cin - app.SimResults.C_NaOH) / app.Cin) * 100;
                    conv_EtOAc_opt = ((app.Cin - app.SimResults.C_EtOAc) / app.Cin) * 100;
                    
                    plot(ax8, app.SimResults.t, conv_NaOH_opt, '-', ...
                         'LineWidth', 3, 'Color', [255 59 48]/255, 'DisplayName', 'Optimized NaOH Conversion');
                    plot(ax8, app.SimResults.t, conv_EtOAc_opt, '-', ...
                         'LineWidth', 3, 'Color', [255 149 0]/255, 'DisplayName', 'Optimized EtOAc Conversion');
                end
                
                % Quick simulation unreacted conversion
                if ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't')
                    conv_NaOH_quick = ((app.Cin - app.QuickSimResults.C_NaOH) / app.Cin) * 100;
                    conv_EtOAc_quick = ((app.Cin - app.QuickSimResults.C_EtOAc) / app.Cin) * 100;
                    
                    plot(ax8, app.QuickSimResults.t, conv_NaOH_quick, '--', ...
                         'LineWidth', 2.5, 'Color', [255 80 80]/255, 'DisplayName', 'Quick Sim NaOH Conversion');
                    plot(ax8, app.QuickSimResults.t, conv_EtOAc_quick, '--', ...
                         'LineWidth', 2.5, 'Color', [255 180 50]/255, 'DisplayName', 'Quick Sim EtOAc Conversion');
                end
                
                title(ax8, 'Conversion: Unreacted Reactants vs Time', 'FontWeight', 'bold', 'FontSize', 14);
                xlabel(ax8, 'Time (min)', 'FontWeight', 'bold');
                ylabel(ax8, 'Conversion (%)', 'FontWeight', 'bold');
                legend(ax8, 'Location', 'best', 'Box', 'on');
                grid(ax8, 'off'); % No grid
                ylim(ax8, [0 100]);
                
                filename = sprintf('Conversion_Unreacted_Reactants_%s.png', timestamp);
                if exist('exportgraphics', 'file')
                    exportgraphics(ax8, fullfile(folder, filename), 'Resolution', 300);
                else
                    print(fig8, fullfile(folder, filename), '-dpng', '-r300');
                end
                close(fig8);
                exportedCount = exportedCount + 1;
            catch ME
                fprintf('Conversion reactants plot export failed: %s\n', ME.message);
            end
        end
          app.updateLoadingDialog(0.95, 'Exporting combined concentration plots...');
        
        % 9. Export combined concentration plots (NEW)
        app.exportCombinedConcentrationPlots(folder, timestamp);
        exportedCount = exportedCount + 2; % Add 2 more plots to the count
        app.updateLoadingDialog(1.0, 'Export complete!');
        app.closeLoadingDialog();
        app.updateLoadingDialog(0.96, 'Exporting concentration profiles with convergence markers...');
        
        % 10. Export concentration profiles with detailed convergence analysis (NEW)
        try
            app.exportConcentrationWithConvergenceMarkers(folder, timestamp);
            exportedCount = exportedCount + 1; % Add 1 more plot to the count
        catch ME
            fprintf('Concentration convergence analysis plot export failed: %s\n', ME.message);
        end
        
        % Then continue with the existing final message code...
        app.updateLoadingDialog(1.0, 'Export complete!');
        
        if exportedCount > 0
            app.updateStatus(sprintf('📷 Successfully exported %d enhanced plots including convergence analysis!', exportedCount), 'success');
            app.flashSuccess(app.ExportPlotsButton);
            
            % Show detailed success message (MODIFY existing message to include convergence plot)
            msgText = sprintf(['🎉 Comprehensive Plot Export Complete!\n\n' ...
                              '📁 Exported to: %s\n\n' ...
                              '📊 Successfully exported %d plots:\n' ...
                              '  • Concentration Profiles\n' ...
                              '  • Temperature Profiles\n' ...
                              '  • Product Formation\n' ...
                              '  • Quick Sim Concentrations\n' ...
                              '  • Quick Sim Temperature\n' ...
                              '  • Quick Sim Products\n' ...
                              '  • Conversion Product Formation\n' ...
                              '  • Conversion Unreacted Reactants\n' ...
                              '  • Combined Concentration (15 min)\n' ...
                              '  • Combined Concentration (Full Time)\n' ...
                              '  • Concentration with Convergence Analysis\n\n' ...
                              '✨ Features:\n' ...
                              '  • High resolution (300 DPI)\n' ...
                              '  • Convergence markers for all species\n' ...
                              '  • Multiple simulation methods\n' ...
                              '  • Differential convergence criteria\n' ...
                              '  • Professional formatting\n' ...
                              '  • Reactants & products analysis'], ...
                              folder, exportedCount);
            
            uialert(app.UIFigure, msgText, 'Plot Export Complete', 'Icon', 'success');
        else
            app.updateStatus('⚠️ No plots were successfully exported', 'warning');
            uialert(app.UIFigure, 'No plots were successfully exported. Please check that simulations have been run and try again.', 'Export Warning', 'Icon', 'warning');
        end
        
    catch ME
        app.closeLoadingDialog();
        app.updateStatus(['Plot export failed: ' ME.message], 'error');
        fprintf('Plot export error: %s\n', ME.message);
        uialert(app.UIFigure, ['Plot export failed: ', ME.message], 'Export Error', 'Icon', 'error');
    end
end

function exportConcentrationWithConvergenceMarkers(app, folder, timestamp)
    % Export concentration profiles with detailed convergence markers for all species
    % Shows reactants and products with their respective convergence points and lines
    
    try
        exportedCount = 0;
        
        % Check if we have simulation data and convergence data to export
        hasOptimized = ~isempty(app.SimResults) && isfield(app.SimResults, 't');
        hasManual = ~isempty(app.ManualSimResults) && isfield(app.ManualSimResults, 't');
        hasUltra = ~isempty(app.UltraPrecisionResults) && isfield(app.UltraPrecisionResults, 't');
        hasQuickSim = ~isempty(app.QuickSimResults) && isfield(app.QuickSimResults, 't');
        
        if ~(hasOptimized || hasManual || hasUltra || hasQuickSim)
            return;
        end
        
        % 1. COMPREHENSIVE CONVERGENCE ANALYSIS PLOT
        app.updateLoadingDialog(0.92, 'Creating convergence analysis plot...');
        
        try
            fig = figure('Visible', 'off', 'Position', [50, 50, 1600, 1000]);
            
            % Create subplots: 2x2 layout
            % Top left: Reactants with convergence
            ax1 = subplot(2, 2, 1);
            hold(ax1, 'on');
            
            % Top right: Products with convergence  
            ax2 = subplot(2, 2, 2);
            hold(ax2, 'on');
            
            % Bottom: Combined view with all convergence markers
            ax3 = subplot(2, 1, 2);
            hold(ax3, 'on');
            
            %% REACTANTS PLOT (Top Left)
            % Plot experimental data if available
            if ~isempty(app.Data)
                scatter(ax1, app.Data.Time, app.Data.Concentration, 100, ...
                       'MarkerEdgeColor', [0.8 0.2 0.2], 'MarkerFaceColor', [1 0.9 0.9], ...
                       'LineWidth', 2.5, 'DisplayName', 'Experimental NaOH', 'Marker', 'o');
            end
            
            % Manual Simulation Reactants
            if hasManual
                plot(ax1, app.ManualSimResults.t, app.ManualSimResults.C_NaOH, ':', ...
                     'LineWidth', 3, 'Color', [0.9 0.3 0.3], 'DisplayName', 'Manual NaOH');
                plot(ax1, app.ManualSimResults.t, app.ManualSimResults.C_EtOAc, ':', ...
                     'LineWidth', 3, 'Color', [0.9 0.6 0.2], 'DisplayName', 'Manual EtOAc');
                
                % Add manual convergence markers
                if ~isempty(app.ManualConvergenceData)
                    app.addReactantConvergenceMarkersDetailed(ax1, app.ManualSimResults, app.ManualConvergenceData, 'Manual');
                end
            end
            
            % Optimized Simulation Reactants
            if hasOptimized
                plot(ax1, app.SimResults.t, app.SimResults.C_NaOH, '-', ...
                     'LineWidth', 3, 'Color', [0.2 0.4 0.8], 'DisplayName', 'Optimized NaOH');
                plot(ax1, app.SimResults.t, app.SimResults.C_EtOAc, '-', ...
                     'LineWidth', 3, 'Color', [0.3 0.7 0.3], 'DisplayName', 'Optimized EtOAc');
                
                % Add optimized convergence markers
                if ~isempty(app.ConvergenceData)
                    app.addReactantConvergenceMarkersDetailed(ax1, app.SimResults, app.ConvergenceData, 'Optimized');
                end
            end
            
            % Ultra-Precision Simulation Reactants
            if hasUltra
                plot(ax1, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, '-', ...
                     'LineWidth', 4, 'Color', [1 0.08 0.58], 'DisplayName', 'Ultra-Precision NaOH');
                plot(ax1, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOAc, '-', ...
                     'LineWidth', 4, 'Color', [0.54 0 0.54], 'DisplayName', 'Ultra-Precision EtOAc');
                
                % Add ultra-precision convergence markers
                if ~isempty(app.UltraConvergenceData)
                    app.addReactantConvergenceMarkersDetailed(ax1, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra-Precision');
                end
                
            end
            
            % Quick Simulation Reactants
            if hasQuickSim
                plot(ax1, app.QuickSimResults.t, app.QuickSimResults.C_NaOH, '-', ...
                     'LineWidth', 2.5, 'Color', [0.0 0.5 0.0], 'DisplayName', 'Quick Sim NaOH');
                plot(ax1, app.QuickSimResults.t, app.QuickSimResults.C_EtOAc, '-', ...
                     'LineWidth', 2.5, 'Color', [1.0 0.65 0], 'DisplayName', 'Quick Sim EtOAc');
                
                % Add quick sim convergence markers
                if ~isempty(app.QuickSimConvergenceData)
                    app.addReactantConvergenceMarkersDetailed(ax1, app.QuickSimResults, app.QuickSimConvergenceData, 'QuickSim');
                end
            end
            
            title(ax1, 'Reactants Consumption with Convergence Markers', 'FontWeight', 'bold', 'FontSize', 14);
            xlabel(ax1, 'Time (min)', 'FontWeight', 'bold');
            ylabel(ax1, 'Concentration (mol/L)', 'FontWeight', 'bold');
            legend(ax1, 'Location', 'best', 'FontSize', 10);
            grid(ax1, 'off'); ax1.GridAlpha = 0.3;
            
            %% PRODUCTS PLOT (Top Right)
            % Manual Simulation Products
            if hasManual
                plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_NaOAc, ':', ...
                     'LineWidth', 3, 'Color', [0.8 0.7 0.2], 'DisplayName', 'Manual NaOAc');
                plot(ax2, app.ManualSimResults.t, app.ManualSimResults.C_EtOH, ':', ...
                     'LineWidth', 3, 'Color', [0.6 0.4 0.8], 'DisplayName', 'Manual EtOH');
                
                % Add manual product convergence markers
                if ~isempty(app.ManualConvergenceData)
                    app.addProductConvergenceMarkersDetailed(ax2, app.ManualSimResults, app.ManualConvergenceData, 'Manual');
                end
            end
            
            % Optimized Simulation Products
            if hasOptimized
                plot(ax2, app.SimResults.t, app.SimResults.C_NaOAc, '-', ...
                     'LineWidth', 3, 'Color', [0.9 0.6 0.0], 'DisplayName', 'Optimized NaOAc');
                plot(ax2, app.SimResults.t, app.SimResults.C_EtOH, '-', ...
                     'LineWidth', 3, 'Color', [0.5 0.3 0.7], 'DisplayName', 'Optimized EtOH');
                
                % Add optimized product convergence markers
                if ~isempty(app.ConvergenceData)
                    app.addProductConvergenceMarkersDetailed(ax2, app.SimResults, app.ConvergenceData, 'Optimized');
                end
            end
            
            % Ultra-Precision Simulation Products
            if hasUltra
                plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '-', ...
                     'LineWidth', 4, 'Color', [1 0.27 0], 'DisplayName', 'Ultra-Precision NaOAc');
                plot(ax2, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '-', ...
                     'LineWidth', 4, 'Color', [0.13 0.55 0.13], 'DisplayName', 'Ultra-Precision EtOH');
                
                % Add ultra-precision product convergence markers
                if ~isempty(app.UltraConvergenceData)
                    app.addProductConvergenceMarkersDetailed(ax2, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra-Precision');
                end
                 if ~isempty(app.UltraConvergenceData)
                    app.addProductConvergenceMarkers(ax3, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra-Precision');
                end
            end
            
            % Quick Simulation Products
            if hasQuickSim
                plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_NaOAc, '-', ...
                     'LineWidth', 2.5, 'Color', [0.18 0.31 0.31], 'DisplayName', 'Quick Sim NaOAc');
                plot(ax2, app.QuickSimResults.t, app.QuickSimResults.C_EtOH, '-', ...
                     'LineWidth', 2.5, 'Color', [0.5 0.0 0.5], 'DisplayName', 'Quick Sim EtOH');
                
                % Add quick sim product convergence markers
                if ~isempty(app.QuickSimConvergenceData)
                    app.addProductConvergenceMarkersDetailed(ax2, app.QuickSimResults, app.QuickSimConvergenceData, 'QuickSim');
                end
            end
            
            title(ax2, 'Products Formation with Convergence Markers', 'FontWeight', 'bold', 'FontSize', 14);
            xlabel(ax2, 'Time (min)', 'FontWeight', 'bold');
            ylabel(ax2, 'Concentration (mol/L)', 'FontWeight', 'bold');
            legend(ax2, 'Location', 'best', 'FontSize', 10);
            grid(ax2, 'off'); ax2.GridAlpha = 0.3;
            
            %% COMBINED VIEW (Bottom)
            % Plot experimental data
            if ~isempty(app.Data)
                scatter(ax3, app.Data.Time, app.Data.Concentration, 80, ...
                       'MarkerEdgeColor', [0.8 0.2 0.2], 'MarkerFaceColor', [1 0.9 0.9], ...
                       'LineWidth', 2, 'DisplayName', 'Experimental NaOH', 'Marker', 'o');
            end
            
            % Plot all species from best available method (prefer ultra, then optimized, then manual)
            if hasUltra
                % Ultra-precision all species
                plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOH, '-', ...
                     'LineWidth', 3, 'Color', [1 0.08 0.58], 'DisplayName', 'Ultra NaOH');
                plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOAc, '-', ...
                     'LineWidth', 3, 'Color', [0.54 0 0.54], 'DisplayName', 'Ultra EtOAc');
                plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_NaOAc, '--', ...
                     'LineWidth', 3, 'Color', [1 0.27 0], 'DisplayName', 'Ultra NaOAc');
                plot(ax3, app.UltraPrecisionResults.t, app.UltraPrecisionResults.C_EtOH, '--', ...
                     'LineWidth', 3, 'Color', [0.13 0.55 0.13], 'DisplayName', 'Ultra EtOH');
                
                % Add all convergence markers
                if ~isempty(app.UltraConvergenceData)
                    app.addAllConvergenceMarkersDetailed(ax3, app.UltraPrecisionResults, app.UltraConvergenceData, 'Ultra-Precision');
                end
                method_used = 'Ultra-Precision';
                
            elseif hasOptimized
                % Optimized all species
                plot(ax3, app.SimResults.t, app.SimResults.C_NaOH, '-', ...
                     'LineWidth', 3, 'Color', [0.2 0.4 0.8], 'DisplayName', 'Opt NaOH');
                plot(ax3, app.SimResults.t, app.SimResults.C_EtOAc, '-', ...
                     'LineWidth', 3, 'Color', [0.3 0.7 0.3], 'DisplayName', 'Opt EtOAc');
                plot(ax3, app.SimResults.t, app.SimResults.C_NaOAc, '--', ...
                     'LineWidth', 3, 'Color', [0.9 0.6 0.0], 'DisplayName', 'Opt NaOAc');
                plot(ax3, app.SimResults.t, app.SimResults.C_EtOH, '--', ...
                     'LineWidth', 3, 'Color', [0.5 0.3 0.7], 'DisplayName', 'Opt EtOH');
                
                % Add all convergence markers
                if ~isempty(app.ConvergenceData)
                    app.addAllConvergenceMarkersDetailed(ax3, app.SimResults, app.ConvergenceData, 'Optimized');
                end
                method_used = 'Optimized';
                
            elseif hasManual
                % Manual all species
                plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_NaOH, '-', ...
                     'LineWidth', 3, 'Color', [0.9 0.3 0.3], 'DisplayName', 'Manual NaOH');
                plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_EtOAc, '-', ...
                     'LineWidth', 3, 'Color', [0.9 0.6 0.2], 'DisplayName', 'Manual EtOAc');
                plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_NaOAc, '--', ...
                     'LineWidth', 3, 'Color', [0.8 0.7 0.2], 'DisplayName', 'Manual NaOAc');
                plot(ax3, app.ManualSimResults.t, app.ManualSimResults.C_EtOH, '--', ...
                     'LineWidth', 3, 'Color', [0.6 0.4 0.8], 'DisplayName', 'Manual EtOH');
                
                % Add all convergence markers
                if ~isempty(app.ManualConvergenceData)
                    app.addAllConvergenceMarkersDetailed(ax3, app.ManualSimResults, app.ManualConvergenceData, 'Manual');
                end
                method_used = 'Manual';
            end
            
            % Create comprehensive title
            if ~isnan(app.k_fit)
                titleStr = sprintf('Complete CSTR Analysis', ...
                                  method_used, app.k_fit);
            elseif ~isnan(app.k_manual)
                titleStr = sprintf('Complete CSTR Analysis', ...
                                  method_used, app.k_manual);
            else
                titleStr = sprintf('Complete CSTR Analysis', method_used);
            end
            
            title(ax3, titleStr, 'FontWeight', 'bold', 'FontSize', 15);
            xlabel(ax3, 'Time (min)', 'FontWeight', 'bold', 'FontSize', 12);
            ylabel(ax3, 'Concentration (mol/L)', 'FontWeight', 'bold', 'FontSize', 12);
            legend(ax3, 'Location', 'eastoutside', 'FontSize', 10);
            grid(ax3, 'off'); ax3.GridAlpha = 0.3;
            
            % Add overall figure title
            sgtitle(sprintf('CSTR Saponification: Concentration Profiles Analysis\nReactor: %.1fL, T: %.0f°C, Flowrate: %.3f L/min', ...
                           app.V, app.SelectedTemp, app.FlowrateField.Value), ...
                           'FontSize', 16, 'FontWeight', 'bold');
            
            % Adjust subplot spacing
            set(fig, 'Units', 'normalized');
            subplot(2, 2, 1); set(gca, 'Position', [0.08, 0.55, 0.38, 0.35]);
            subplot(2, 2, 2); set(gca, 'Position', [0.55, 0.55, 0.38, 0.35]);
            subplot(2, 1, 2); set(gca, 'Position', [0.08, 0.08, 0.75, 0.35]);
            
            % Export the comprehensive convergence plot
            filename = sprintf('Concentration_Convergence_Analysis_%s.png', timestamp);
            if exist('exportgraphics', 'file')
                exportgraphics(fig, fullfile(folder, filename), 'Resolution', 300);
            else
                print(fig, fullfile(folder, filename), '-dpng', '-r300');
            end
            close(fig);
            exportedCount = exportedCount + 1;
            
            fprintf('Successfully exported convergence analysis plot: %s\n', filename);
            
        catch ME
            fprintf('Convergence analysis plot export failed: %s\n', ME.message);
        end
        
    catch ME
        fprintf('Export concentration with convergence markers failed: %s\n', ME.message);
    end
end

% Helper function: Add detailed reactant convergence markers
function addReactantConvergenceMarkersDetailed(app, axes, simResults, convergenceData, method)
    try
        if isempty(convergenceData)
            return;
        end
        
        % Define method-specific styling
        switch method
            case 'Manual'
                lineAlpha = 0.6; markerSize = 70; lineStyle = '-.';
                naohColor = [0.9 0.3 0.3]; etoacColor = [0.9 0.6 0.2];
            case 'Optimized'
                lineAlpha = 0.8; markerSize = 80; lineStyle = '--';
                naohColor = [0.2 0.4 0.8]; etoacColor = [0.3 0.7 0.3];
            case 'Ultra-Precision'
                lineAlpha = 0.9; markerSize = 100; lineStyle = '-';
                naohColor = [1 0.08 0.58]; etoacColor = [0.54 0 0.54];
            case 'QuickSim'
                lineAlpha = 0.7; markerSize = 75; lineStyle = ':';
                naohColor = [0.0 0.5 0.0]; etoacColor = [1.0 0.65 0];
            otherwise
                lineAlpha = 0.8; markerSize = 80; lineStyle = '--';
                naohColor = [0.5 0.5 0.5]; etoacColor = [0.7 0.7 0.7];
        end
        
        % NaOH convergence
        if isfield(convergenceData, 'NaOH_converged') && convergenceData.NaOH_converged
            conv_time = convergenceData.NaOH_time;
            conv_conc = convergenceData.NaOH_conc;
            
            xline(axes, conv_time, lineStyle, 'Color', naohColor, 'Alpha', lineAlpha, ...
                  'LineWidth', 2, 'DisplayName', sprintf('%s NaOH Conv.', method));
            scatter(axes, conv_time, conv_conc, markerSize, 'filled', ...
                   'MarkerFaceColor', naohColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 1.5, 'DisplayName', sprintf('%s NaOH Point', method));
        end
        
        % EtOAc convergence
        if isfield(convergenceData, 'EtOAc_converged') && convergenceData.EtOAc_converged
            conv_time = convergenceData.EtOAc_time;
            conv_conc = convergenceData.EtOAc_conc;
            
            xline(axes, conv_time, lineStyle, 'Color', etoacColor, 'Alpha', lineAlpha, ...
                  'LineWidth', 2, 'DisplayName', sprintf('%s EtOAc Conv.', method));
            scatter(axes, conv_time, conv_conc, markerSize, 'filled', ...
                   'MarkerFaceColor', etoacColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 1.5, 'DisplayName', sprintf('%s EtOAc Point', method));
        end
        
    catch ME
        fprintf('Error adding reactant convergence markers: %s\n', ME.message);
    end
end

% Helper function: Add detailed product convergence markers
function addProductConvergenceMarkersDetailed(app, axes, simResults, convergenceData, method)
    try
        if isempty(convergenceData)
            return;
        end
        
        % Define method-specific styling
        switch method
            case 'Manual'
                lineAlpha = 0.6; markerSize = 70; lineStyle = '-.';
                naoacColor = [0.8 0.7 0.2]; etohColor = [0.6 0.4 0.8];
            case 'Optimized'
                lineAlpha = 0.8; markerSize = 80; lineStyle = '--';
                naoacColor = [0.9 0.6 0.0]; etohColor = [0.5 0.3 0.7];
            case 'Ultra-Precision'
                lineAlpha = 0.9; markerSize = 100; lineStyle = '-';
                naoacColor = [1 0.27 0]; etohColor = [0.13 0.55 0.13];
            case 'QuickSim'
                lineAlpha = 0.7; markerSize = 75; lineStyle = ':';
                naoacColor = [0.18 0.31 0.31]; etohColor = [0.5 0.0 0.5];
            otherwise
                lineAlpha = 0.8; markerSize = 80; lineStyle = '--';
                naoacColor = [0.5 0.5 0.5]; etohColor = [0.7 0.7 0.7];
        end
        
        % NaOAc convergence
        if isfield(convergenceData, 'NaOAc_converged') && convergenceData.NaOAc_converged
            conv_time = convergenceData.NaOAc_time;
            conv_conc = convergenceData.NaOAc_conc;
            
            xline(axes, conv_time, lineStyle, 'Color', naoacColor, 'Alpha', lineAlpha, ...
                  'LineWidth', 2, 'DisplayName', sprintf('%s NaOAc Conv.', method));
            scatter(axes, conv_time, conv_conc, markerSize, 'filled', ...
                   'MarkerFaceColor', naoacColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 1.5, 'DisplayName', sprintf('%s NaOAc Point', method));
        end
        
        % EtOH convergence
        if isfield(convergenceData, 'EtOH_converged') && convergenceData.EtOH_converged
            conv_time = convergenceData.EtOH_time;
            conv_conc = convergenceData.EtOH_conc;
            
            xline(axes, conv_time, lineStyle, 'Color', etohColor, 'Alpha', lineAlpha, ...
                  'LineWidth', 2, 'DisplayName', sprintf('%s EtOH Conv.', method));
            scatter(axes, conv_time, conv_conc, markerSize, 'filled', ...
                   'MarkerFaceColor', etohColor, 'MarkerEdgeColor', 'black', ...
                   'LineWidth', 1.5, 'DisplayName', sprintf('%s EtOH Point', method));
        end
        
    catch ME
        fprintf('Error adding product convergence markers: %s\n', ME.message);
    end
end

% Helper function: Add all convergence markers for combined view
function addAllConvergenceMarkersDetailed(app, axes, simResults, convergenceData, method)
    try
        if isempty(convergenceData)
            return;
        end
        
        % Define colors for all species based on method
        switch method
            case 'Manual'
                colors = struct('NaOH', [0.9 0.3 0.3], 'EtOAc', [0.9 0.6 0.2], ...
                               'NaOAc', [0.8 0.7 0.2], 'EtOH', [0.6 0.4 0.8]);
                lineAlpha = 0.6; markerSize = 60;
            case 'Optimized'
                colors = struct('NaOH', [0.2 0.4 0.8], 'EtOAc', [0.3 0.7 0.3], ...
                               'NaOAc', [0.9 0.6 0.0], 'EtOH', [0.5 0.3 0.7]);
                lineAlpha = 0.8; markerSize = 70;
            case 'Ultra-Precision'
                colors = struct('NaOH', [1 0.08 0.58], 'EtOAc', [0.54 0 0.54], ...
                               'NaOAc', [1 0.27 0], 'EtOH', [0.13 0.55 0.13]);
                lineAlpha = 0.9; markerSize = 90;
            otherwise
                colors = struct('NaOH', [0.5 0.5 0.5], 'EtOAc', [0.6 0.6 0.6], ...
                               'NaOAc', [0.7 0.7 0.7], 'EtOH', [0.8 0.8 0.8]);
                lineAlpha = 0.7; markerSize = 65;
        end
        
        species = {'NaOH', 'EtOAc', 'NaOAc', 'EtOH'};
        
        for i = 1:length(species)
            speciesName = species{i};
            convergenceField = sprintf('%s_converged', speciesName);
            timeField = sprintf('%s_time', speciesName);
            concField = sprintf('%s_conc', speciesName);
            
            if isfield(convergenceData, convergenceField) && convergenceData.(convergenceField)
                conv_time = convergenceData.(timeField);
                conv_conc = convergenceData.(concField);
                
                % Add vertical line
                xline(axes, conv_time, '--', 'Color', colors.(speciesName), 'Alpha', lineAlpha, ...
                      'LineWidth', 1.5, 'DisplayName', sprintf('%s %s Conv.', method, speciesName));
                
                % Add convergence point
                scatter(axes, conv_time, conv_conc, markerSize, 'filled', ...
                       'MarkerFaceColor', colors.(speciesName), 'MarkerEdgeColor', 'black', ...
                       'LineWidth', 1, 'DisplayName', sprintf('%s %s Point', method, speciesName));
                
                % Add small text label
                text(axes, conv_time + 0.5, conv_conc, sprintf('%s', speciesName), ...
                     'FontSize', 8, 'Color', colors.(speciesName), 'FontWeight', 'bold');
            end
        end
        
    catch ME
        fprintf('Error adding all convergence markers: %s\n', ME.message);
    end
end
        
function startupFcn(app)
    
    % Initialize cooling history tracking
app.CoolingTimeHistory = [];
app.CoolingFlowrateHistory = [];
app.CoolingCapacityHistory = [];
app.CoolingAdequacyHistory = [];
app.TemperatureDeviationHistory = [];
    % Initialize extension handling
app.OriginalSimTime = NaN;
app.OriginalEndConcentration = NaN;
app.ExtensionMode = false;
    try
        % Initialize design system
        app.setupDesignSystem();
        
        % Setup animated logo
        app.animateLogo();
        
        % Initialize parameters
        app.SelectedTemp = 30;
        app.V = 2.5;
        app.V0 = 0.18;
        app.V0_EtOAc = 0.09;
        app.V1_NaOH = 0.09;
        app.Cin = 0.05;
        app.SimTime = 30;
        
        % Initialize convergence status
        app.IsConverged = false;
        app.IsManualConverged = false;
        app.IsQuickSimConverged = false;
        
        % Initialize separate convergence tracking
        app.ManualReactantsConverged = false;
        app.ManualProductsConverged = false;
        app.ManualTemperatureConverged = false;
        app.OptimizedReactantsConverged = false;
        app.OptimizedProductsConverged = false;
        app.OptimizedTemperatureConverged = false;
        app.QuickSimReactantsConverged = false;
        app.QuickSimProductsConverged = false;
        app.QuickSimTemperatureConverged = false;
        
        % Initialize tables with proper empty structure
        app.Data = table();
        app.ComparisonTableData = table();
        app.ComparisonTableDataManual = table();
        
        % Initialize simulation results
        app.SimResults = struct();
        app.ManualSimResults = struct();
        app.QuickSimResults = struct();
        app.ConvergenceData = struct();
        app.ManualConvergenceData = struct();
        app.QuickSimConvergenceData = struct();
        app.OptimizationResults = struct();
        
        % Initialize rate constants
        app.k_fit = NaN;
        app.k_manual = NaN;
        app.RSquared = NaN;
        app.RSquaredManual = NaN;
        
        % Set default display options
        app.ShowJacket = false;
        app.ShowProducts = true;
        
        % Initialize info panels
        app.ConcentrationInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
        app.TemperatureInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
        app.ProductsInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
        app.QuickSimInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No simulation data available.</div>';
        % Add this line in the startupFcn after other info panel initializations:
        app.ConversionInfo.HTMLSource = '<div style="font-family: -apple-system; padding: 10px;">No conversion data available.</div>';
        % Set initial convergence text
        app.ConvergenceText.Value = 'No analysis data available. Run a simulation first.';
        
        % Initialize table displays
        app.ComparisonTable.Data = [];
        app.ComparisonTable.ColumnName = {};
        app.ConcentrationTable.Data = [];
        app.ConcentrationTable.ColumnName = {'Time', 'NaOH', 'EtOAc', 'NaOAc', 'EtOH'};
        app.TemperatureTable.Data = [];
        app.TemperatureTable.ColumnName = {'Time', 'Reactor', 'Jacket'};
        
        % Initialize Quick Sim error labels
        app.MaxErrorLabel.Text = 'Max Error: N/A';
        app.MinErrorLabel.Text = 'Min Error: N/A';
        app.MeanErrorLabel.Text = 'Mean Error: N/A';
        
        % NEW: Welcome message with look-ahead convergence info
        app.updateStatus('🚀 Enhanced CSTR Simulator ready! NEW: Look-ahead convergence (next 5 pts <5% error)', 'success');
        % ADD these properties to the startupFcn function after other initializations:

% Initialize ultra-precision convergence tracking
app.IsUltraConverged = false;
app.UltraReactantsConverged = false;
app.UltraProductsConverged = false;
app.UltraTemperatureConverged = false;
app.UltraConvergenceData = struct();

% Initialize ultra-precision results
app.UltraPrecisionResults = struct();
app.ComparisonTableDataUltra = table();
app.RSquaredUltra = NaN;
app.k_ultra = NaN;

% Add ultra-precision thresholds (strictest criteria)
app.UltraReactantThreshold = 0.0001;   % 0.01% relative change for reactants
app.UltraProductThreshold = 0.00001;   % 0.001% relative change for products
app.UltraTemperatureThreshold = 0.0001; % 0.01% relative change for temperature
        
        % ADD THIS LINE before the final catch:
        app.initializeAdaptiveCooling();

    catch ME
        app.updateStatus(['Startup failed: ' ME.message], 'error');
    end
end

        
        % Value changed function: TempDropDown
        function TempDropDownValueChanged(app, ~)
            app.SelectedTemp = str2double(app.TempDropDown.Value);
            app.updateStatus(sprintf('Temperature set to %d°C', app.SelectedTemp), 'info');
        end
        
        % Value changed function: JacketSwitch
        function JacketSwitchValueChanged(app, ~)
            app.ShowJacket = strcmp(app.JacketSwitch.Value, 'On');
            if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
                app.plotTemperatures();
                app.updateStatus(sprintf('Jacket display %s', lower(app.JacketSwitch.Value)), 'info');
            end
        end
        
        % Value changed function: ProductsSwitch
        function ProductsSwitchValueChanged(app, ~)
            app.ShowProducts = strcmp(app.ProductsSwitch.Value, 'On');
            if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
                app.plotProducts();
                app.updateStatus(sprintf('Products display %s', lower(app.ProductsSwitch.Value)), 'info');
            end
        end
        
function UploadButtonPushed(app, ~)
    try
        % Clear any existing data first to prevent conflicts
        app.Data = table();
        app.ComparisonTableData = table();
        app.ComparisonTableDataManual = table();
        
        [file, path] = uigetfile({'*.xlsx;*.xls;*.csv;*.txt', 'Data Files (*.xlsx,*.xls,*.csv,*.txt)'}, ...
                                'Select Experimental Data File');
        if isequal(file, 0)
            return;
        end
        
        app.showLoadingDialog('Loading Data', 'Reading experimental data...');
        
        fullPath = fullfile(path, file);
        [~, ~, ext] = fileparts(file);
        
        % Initialize data variable
        data = table();
        
        switch lower(ext)
            case {'.xlsx', '.xls'}
                data = readtable(fullPath);
            case '.csv'
                data = readtable(fullPath);
            case '.txt'
                data = readtable(fullPath, 'Delimiter', '\t');
            otherwise
                error('Unsupported file format');
        end
        
        app.updateLoadingDialog(0.5, 'Validating data...');
        
        % Validate data format
        if width(data) < 2
            error('Data must have at least 2 columns (Time, Concentration)');
        end
        
        % Ensure we have numeric data
        if ~isnumeric(data{:,1}) || ~isnumeric(data{:,2})
            error('First two columns must contain numeric data');
        end
        
        % Assume first column is time, second is concentration
        data.Properties.VariableNames{1} = 'Time';
        data.Properties.VariableNames{2} = 'Concentration';
        
        % Remove invalid rows
        validRows = ~isnan(data.Time) & ~isnan(data.Concentration) & ...
                   isfinite(data.Time) & isfinite(data.Concentration) & ...
                   data.Time >= 0 & data.Concentration >= 0;
        data = data(validRows, :);
        
        if height(data) < 3
            error('Need at least 3 valid data points after filtering');
        end
        
        % Sort by time
        data = sortrows(data, 'Time');
        
        % Store the data
        app.Data = data;
        app.FileLabel.Text = sprintf('%s (%d points)', file, height(data));
        
        % Clear any existing simulation results to prevent conflicts
        app.SimResults = struct();
        app.ManualSimResults = struct();
        app.QuickSimResults = struct();
        app.ConvergenceData = struct();
        app.ManualConvergenceData = struct();
        app.QuickSimConvergenceData = struct();
        
        % Reset convergence flags
        app.IsConverged = false;
        app.IsManualConverged = false;
        app.IsQuickSimConverged = false;
        
        app.updateLoadingDialog(1.0, 'Data loaded successfully!');
        app.closeLoadingDialog();
        
        app.updateStatus(sprintf('📈 Data loaded: %d points from %s', height(data), file), 'success');
        app.flashSuccess(app.UploadButton);
        
    catch ME
        app.closeLoadingDialog();
        % Clear data on error to ensure clean state
        app.Data = table();
        app.FileLabel.Text = 'No file selected';
        app.updateStatus(['Upload failed: ' ME.message], 'error');
        uialert(app.UIFigure, ['Failed to load data: ', ME.message], 'Data Load Error', 'Icon', 'error');
    end
end
        
        % Button pushed function: ManualKButton
        function ManualKButtonPushed(app, ~)
            if isempty(app.Data)
                app.updateStatus('Upload experimental data first!', 'error');
                uialert(app.UIFigure, 'Please upload experimental data first.', 'Missing Data', 'Icon', 'warning');
                return;
            end
            app.runManualSimulation();
            app.plotConcentrations();
            app.plotTemperatures();
            app.plotProducts();
             app.plotConversions(); 
            app.updateAnalysisDisplay();
        end
        
        % Button pushed function: OptimizeButton
        function OptimizeButtonPushed(app, ~)
            app.runSimulation();
        end
        
        % Button pushed function: OptimizeUltraButton
        function OptimizeUltraButtonPushed(app, ~)
            app.optimizeUltraPrecision();
        end
        
        % Button pushed function: ClearAllButton
        function ClearAllButtonPushed(app, ~)
            % Confirmation dialog
            selection = uiconfirm(app.UIFigure, ...
                'This will clear all data, simulations, and plots. Are you sure?', ...
                'Clear All Data', ...
                'Options', {'Yes, Clear All', 'Cancel'}, ...
                'DefaultOption', 2, ...
                'Icon', 'warning');
            
            if strcmp(selection, 'Yes, Clear All')
                app.clearAllData();
            end
        end
        
        % Button pushed function: SaveButton
        function SaveButtonPushed(app, ~)
            app.saveProject();
        end
        
        % Button pushed function: ExportDataButton
        function ExportDataButtonPushed(app, ~)
            app.exportData();
        end
        
        % Button pushed function: ExportPlotsButton
        function ExportPlotsButtonPushed(app, ~)
            app.exportPlots();
        end

function ValidateDataButtonPushed(app, ~)
    if isempty(app.Data)
        app.updateStatus('No data to validate', 'warning');
        uialert(app.UIFigure, 'Please upload experimental data first.', 'No Data');
        return;
    end
    app.enhancedDataValidation();
end

function ExportEnhancedButtonPushed(app, ~)
    app.exportEnhancedResults();
end

% ENHANCED COOLING SYSTEM INITIALIZATION WITH TEMPERATURE-RESPONSIVE PARAMETERS
function initializeAdaptiveCooling(app)
    try
        % Temperature-responsive cooling parameters
        app.Fj_min = 0.05;               % Minimum flowrate [L/min]
        app.Fj_max = 20;                 % Increased maximum flowrate [L/min] 
        app.Fj_adaptive = 1.0;           % Start with reasonable flowrate
        app.isothermal_tolerance = 0.1;   % ±0.1°C control tolerance
        app.cooling_control_gain = 1.5;   % Base PI gain
        app.integral_error = 0;
        app.IsothermalMode = true;
        
        % Initialize tracking arrays
        app.TempDifferenceHistory = [];
        app.TotalCoolingRequiredHistory = [];
        
        % Calculate maximum cooling capacity
        T_max_expected = max(str2double(app.TempDropDown.Items));
        app.max_cooling_capacity = app.Fj_max * app.pj/1000 * app.Cpj * ...
                                  (T_max_expected - app.Tj0);
        
        app.updateStatus('🌡️ Temperature-responsive adaptive cooling system initialized', 'success');
        
    catch ME
        app.updateStatus(['Temperature-responsive cooling initialization failed: ' ME.message], 'error');
    end
end
    
   function IsothermalSwitchValueChanged(app, ~)
    app.IsothermalMode = strcmp(app.IsothermalSwitch.Value, 'On');
    
    % Reset cooling system parameters when switching modes
    app.integral_error = 0;  % Reset PI controller integral term
    
    if app.IsothermalMode
        % Enable isothermal control
        app.Fj_adaptive = app.Fj;  % Start with nominal flowrate
        app.updateStatus('🌡️ Isothermal control ENABLED - Adaptive cooling active', 'success');
        
        % Initialize cooling system for isothermal operation
        app.required_cooling_capacity = 0;
        app.cooling_adequacy_factor = 1.0;
        
    else
        % Disable isothermal control - NO COOLING AT ALL
        app.Fj_adaptive = 0;  % CHANGED: No cooling flowrate (was app.Fj_min)
        app.updateStatus('🔥 Isothermal control DISABLED - NO COOLING! Temperature will rise!', 'warning');
        
        % Set cooling system status for no cooling
        app.required_cooling_capacity = 0;
        app.max_cooling_capacity = 0;
        app.cooling_adequacy_factor = 0;
        
        % Enhanced warning about temperature rise
        if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
            msgText = ['⚠️ CRITICAL: Non-Isothermal Mode Activated\n\n' ...
                      'Cooling jacket COMPLETELY DISCONNECTED:\n' ...
                      '• Zero coolant flowrate (0 L/min)\n' ...
                      '• No heat removal capability\n' ...
                      '• No active temperature control\n' ...
                      '• Risk of thermal runaway!\n\n' ...
                      'Reactor temperature will rise significantly!\n' ...
                      'Monitor carefully for safety.'];
            uialert(app.UIFigure, msgText, 'No Cooling Warning', 'Icon', 'warning');
        end
    end
    
    % Clear cooling history for new mode
    app.CoolingTimeHistory = [];
    app.CoolingFlowrateHistory = [];
    app.CoolingCapacityHistory = [];
    app.CoolingAdequacyHistory = [];
    app.TemperatureDeviationHistory = [];
    app.TempDifferenceHistory = [];
app.TotalCoolingRequiredHistory = [];
    
    % Re-run simulations if data exists to show effect of mode change
    if ~isempty(app.SimResults) || ~isempty(app.ManualSimResults)
        app.updateStatus('♻️ Re-running simulations with new temperature control mode...', 'info');
        
        % Re-run optimized simulation if it exists
        if ~isempty(app.SimResults) && ~isnan(app.k_fit)
            app.runSimulation();
        end
        
        % Re-run manual simulation if it exists  
        if ~isempty(app.ManualSimResults) && ~isnan(app.k_manual)
            app.runManualSimulation();
            app.plotConcentrations();
            app.plotTemperatures();
            app.plotProducts();
            app.plotConversions();
            app.updateAnalysisDisplay();
        end
        
        % Update all plots to reflect new temperature behavior
        if ~isempty(app.SimResults)
            app.plotConcentrations();
            app.plotTemperatures();
            app.plotProducts();
            app.plotConversions();
            app.updateAnalysisDisplay();
        end
    else
        % Just update temperature plot to show mode change
        app.plotTemperatures();
    end
end
    
   function showCoolingSystemReport(app, temp_deviation, cooling_status)
    try
        % Calculate flowrate statistics if available
        if ~isempty(app.CoolingFlowrateHistory)
            avg_flowrate = mean(app.CoolingFlowrateHistory);
            max_flowrate_used = max(app.CoolingFlowrateHistory);
            min_flowrate_used = min(app.CoolingFlowrateHistory);
            flowrate_range = max_flowrate_used - min_flowrate_used;
        else
            avg_flowrate = app.Fj_adaptive;
            max_flowrate_used = app.Fj_adaptive;
            min_flowrate_used = app.Fj_adaptive;
            flowrate_range = 0;
        end
        
        msgText = sprintf(['🌡️ Temperature-Responsive Cooling Performance\n\n' ...
                          '📊 Temperature Control:\n' ...
                          '  • Setpoint: %.1f°C\n' ...
                          '  • Max deviation: %.3f°C\n' ...
                          '  • Control tolerance: ±%.1f°C\n\n' ...
                          '🌊 Adaptive Flowrate Performance:\n' ...
                          '  • Current flowrate: %.3f L/min\n' ...
                          '  • Average flowrate: %.3f L/min\n' ...
                          '  • Peak flowrate: %.3f L/min (%.1f%% of max)\n' ...
                          '  • Flowrate range used: %.3f L/min\n' ...
                          '  • System utilization: %.1f%%\n\n' ...
                          '❄️ Cooling System Status: %s\n' ...
                          '  • Required capacity: %.1f J/min\n' ...
                          '  • Available capacity: %.1f J/min\n' ...
                          '  • Adequacy factor: %.2fx\n\n' ...
                          '🎚️ System Configuration:\n' ...
                          '  • Flowrate range: %.2f - %.1f L/min\n' ...
                          '  • Coolant inlet: %.1f°C'], ...
                          app.SelectedTemp, temp_deviation, app.isothermal_tolerance, ...
                          app.Fj_adaptive, avg_flowrate, max_flowrate_used, ...
                          (max_flowrate_used/app.Fj_max)*100, flowrate_range, ...
                          (avg_flowrate/app.Fj_max)*100, cooling_status, ...
                          app.required_cooling_capacity, app.max_cooling_capacity, ...
                          app.cooling_adequacy_factor, ...
                          app.Fj_min, app.Fj_max, app.Tj0);
        
        % Add temperature-specific recommendations
        if app.SelectedTemp > 50 && avg_flowrate < 2.0
            msgText = [msgText '\n\n💡 RECOMMENDATION: High temperature operation detected. Current flowrate usage is optimal.'];
        elseif app.SelectedTemp <= 30 && avg_flowrate > 1.5
            msgText = [msgText '\n\n💡 INFO: Low temperature operation. Higher flowrate indicates good safety margin.'];
        end
        
        if app.cooling_adequacy_factor < 0.8
            msgText = [msgText '\n\n⚠️ WARNING: Cooling system approaching limits! Consider lower temperature or higher coolant capacity.'];
        end
        
        if flowrate_range > 3.0
            msgText = [msgText '\n\n📈 INFO: Wide flowrate range indicates good adaptive response to temperature changes.'];
        end
        
        uialert(app.UIFigure, msgText, 'Temperature-Responsive Cooling Report', 'Icon', 'info');
        
    catch
        % Silent fail
    end
end
    
    function max_error = calculateMaxErrorIsothermal(app, params, t_exp, C_exp, y0, V, Cin, T_setpoint)
        try
            V0_total = params(1);
            k = params(2);
            
            if k <= 0 || ~isfinite(k) || V0_total <= 0
                max_error = 1e6;
                return;
            end
            
            V0_EtOAc = V0_total / 2;
            V1_NaOH = V0_total / 2;
            
            % Reset integral error for each optimization iteration
            app.integral_error = 0;
            
            % Simulate with isothermal control
            options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
            [~, y_sim] = ode45(@(t,y) app.odefun_individual_isothermal(t,y,k,V0_EtOAc,V1_NaOH,V,Cin,T_setpoint), ...
                              t_exp, y0, options);
            
            C_sim = y_sim(:,1);
            
            % Calculate relative errors
            rel_errors = abs(C_exp - C_sim) ./ C_exp * 100;
            rel_errors(C_exp == 0) = 0;
            
            max_error = max(rel_errors);
            
            % Add penalty for poor temperature control (isothermal requirement)
            temp_deviations = abs(y_sim(:,5) - T_setpoint);
            max_temp_deviation = max(temp_deviations);
            
            if max_temp_deviation > app.isothermal_tolerance
                % Penalize poor temperature control
                temp_penalty = (max_temp_deviation / app.isothermal_tolerance) * 5;
                max_error = max_error + temp_penalty;
            end
            
        catch
            max_error = 1e6;
        end
    end
    function displayEnhancedCoolingReport(app)
    try
        if isempty(app.CoolingTimeHistory)
            return;
        end
        
        % Calculate cooling performance statistics
        avg_flowrate = mean(app.CoolingFlowrateHistory);
        max_flowrate = max(app.CoolingFlowrateHistory);
        avg_temp_diff = mean(app.TempDifferenceHistory);
        max_temp_diff = max(app.TempDifferenceHistory);
        avg_cooling_required = mean(app.TotalCoolingRequiredHistory);
        max_cooling_required = max(app.TotalCoolingRequiredHistory);
        
        % Calculate efficiency metrics
        flowrate_utilization = avg_flowrate / app.Fj_max * 100;
        peak_utilization = max_flowrate / app.Fj_max * 100;
        
        % Temperature control quality
        temp_control_quality = 100 * sum(app.TemperatureDeviationHistory <= app.isothermal_tolerance) / ...
                              length(app.TemperatureDeviationHistory);
        
        msgText = sprintf(['🌡️ Enhanced Cooling System Performance Report\n\n' ...
                          '📊 Flowrate Performance:\n' ...
                          '  • Average flowrate: %.3f L/min (%.1f%% utilization)\n' ...
                          '  • Peak flowrate: %.3f L/min (%.1f%% utilization)\n' ...
                          '  • Flowrate range: %.3f - %.3f L/min\n\n' ...
                          '🌡️ Temperature Control:\n' ...
                          '  • Average temp difference: %.2f°C\n' ...
                          '  • Maximum temp difference: %.2f°C\n' ...
                          '  • Control quality: %.1f%% within tolerance\n\n' ...
                          '❄️ Cooling Demand:\n' ...
                          '  • Average cooling required: %.1f J/min\n' ...
                          '  • Peak cooling required: %.1f J/min\n' ...
                          '  • System adequacy: %.2fx average'], ...
                          avg_flowrate, flowrate_utilization, ...
                          max_flowrate, peak_utilization, ...
                          min(app.CoolingFlowrateHistory), max(app.CoolingFlowrateHistory), ...
                          avg_temp_diff, max_temp_diff, temp_control_quality, ...
                          avg_cooling_required, max_cooling_required, ...
                          mean(app.CoolingAdequacyHistory));
        
        % Add recommendations
        if peak_utilization > 90
            msgText = [msgText '\n\n⚠️ RECOMMENDATION: Peak cooling utilization is high. Consider increasing maximum flowrate.'];
        end
        
        if temp_control_quality < 80
            msgText = [msgText '\n\n⚠️ RECOMMENDATION: Temperature control quality is below 80%. Consider adjusting control parameters.'];
        end
        
        if max_temp_diff > 10
            msgText = [msgText '\n\n⚠️ RECOMMENDATION: Large temperature differences detected. Monitor for thermal runaway risk.'];
        end
        
        uialert(app.UIFigure, msgText, 'Enhanced Cooling Performance', 'Icon', 'info');
        
    catch ME
        app.updateStatus(['Cooling report failed: ' ME.message], 'error');
    end
end


    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CSTR_Saponification_Modern
            
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
            
            % Stop any running timers
            try
                if ~isempty(app.LoadingTimer) && isvalid(app.LoadingTimer)
                    stop(app.LoadingTimer);
                    delete(app.LoadingTimer);
                end
            catch
                %diam ja
            end
            
            % Close any open dialogs
            try
                if ~isempty(app.LoadingDialog) && isvalid(app.LoadingDialog)
                    close(app.LoadingDialog);
                end
            catch
                % Sdiam
            end

try
    if ~isempty(app.CancellableDialog) && isvalid(app.CancellableDialog)
        close(app.CancellableDialog);
    end
catch
    % ignore fail ja
end

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
