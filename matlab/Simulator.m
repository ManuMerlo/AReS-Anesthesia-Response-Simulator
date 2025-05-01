% matlab version of the simulator
classdef Simulator < handle
    properties ( Access = protected )
        patients = Patient.empty;       % List of simulated patients
        
        % List of arrays to store patients' propofol, remifentanil,
        % norepinehrpne, and rocuronoium infusion rates.
        u_prop_all = [];    % [mg * sec^-1]
        u_remi_all = [];    % [µg * sec^-1]
        u_nore_all = [];    % [µg * sec^-1]
        u_rocu_all = [];    % [mg * sec^-1]

        % List of lists to store patients' propofol, remifentanil, 
        % norepinephrine, and rocuronium plasma concentrations and
        % effect site concentrations of propofol, remifentanil, and rocuronium
        
        cp_prop_all = [];   % [µg * mL^-1]
        ce_prop_all = [];   % [µg * mL^-1]
        cp_remi_all = [];   % [ng * mL^-1]
        ce_remi_all = [];   % [ng * mL^-1]
        c_nore_all = [];    % [ng * mL^-1]
        cp_rocu_all = [];   % [µg * mL^-1]
        ce_rocu_all = [];   % [µg * mL^-1]
        

        ce_del_all = [];     % The variable that represents the effect of delay in pd model on propofol effect-site concentraion of all patients
        ce_wav_all = [];     % The variable that represents the effect of WAV filter on propofol effect-site concentraion of all patients
        ce_bis_all = [];     % The variable that represents the effect of BIS filter on propofol effect-site concentraion of all patients
        
        % List of lists to store the WAV and BIS values of patients
        wav_all = [];
        bis_all = [];
        
        % List of lists to store the hemodynamics variables of patients
        co_all = [];        % Cardiac Output [L/min]
        map_all = [];       % Mean Arterial Pressure [mmHg]
        hr_all = [];        % Heart Rate [beats min^-1]
        sv_all = [];        % Stroke Volume [mL]
        tpr_all = []        % Total Peripheral Resistance [mmHg mL^-1 min]
        
        % List of lists to store the probabilities of neuromuscular blockade categories of patients
        nmb_m0_all = [];
        nmb_m1_all = [];
        nmb_m2_all = [];
        nmb_m3_all = [];

        t_sim = 0                % The simulation time in seconds for the current patient
        t_s = 1                  % The sampling time in seconds for the simulation of the current patient
        steps = 0                % Stores the simulation steps. It is updated once a simulation step is completed.               
        current_time = 0;

        nore_molweight = 169.18; % Molecular weight of norepinephrine in [g/mol]
    end

    properties(Access = public)
        current_patient = Patient.empty; % The current patient being simulated
    end
    
    methods
        
        function obj = Simulator()
            % Constructor for the Simulator class.
        end
        
        function reset(obj)
        % Reset all the values of the simulator.
            obj.patients = [];
            obj.u_prop_all = [];
            obj.u_remi_all = [];
            obj.u_nore_all = [];
            obj.u_rocu_all = [];
            obj.cp_prop_all = [];
            obj.cp_remi_all = [];
            obj.c_nore_all = [];
            obj.cp_rocu_all = [];
            obj.ce_prop_all = [];
            obj.ce_remi_all = [];
            obj.ce_rocu_all = [];
            obj.ce_del_all = [];
            obj.ce_wav_all = [];
            obj.ce_bis_all = [];
            obj.wav_all = [];
            obj.bis_all = [];
            obj.co_all = [];
            obj.map_all = [];
            obj.hr_all = [];
            obj.sv_all = [];
            obj.tpr_all = [];
            obj.nmb_m0_all = [];
            obj.nmb_m1_all = [];
            obj.nmb_m2_all = [];
            obj.nmb_m3_all = [];
            obj.current_patient = [];
            obj.t_sim = 0;
            obj.t_s = 0;
            obj.steps = 0;
            obj.current_time = 0;
        end

        function patient_data = get_patient_demographics(obj)
        %  Return the patient data for the current patient: [age, height, weight,gender, bmi, lbm]
            patient_data = obj.current_patient.get_patient_demographics();
        end

        function patient_state = get_patient_state(obj)
        % Return a dictionary with the patient state for the current patient:
        % [WAV, BIS, MAP, CO, HR, SV,NMB_m0, NMB_m1, NMB_m2, NMB_m3, cp_prop, 
        % cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, ce_rocu, ce_del, ce_wav, ce_bis]
            state = obj.current_patient.get_patient_state();
            c_nore = state('c_nore');
            c_nore = c_nore * obj.nore_molweight / 1000; % Norepinephrine blood concentration in [ng/ml]
            state('c_nore') = c_nore;
            patient_state = state;
        end

        function patient_state_history = get_patient_state_history(obj)
        % Return a dictionary with the patient state history from the start
        % of the simulation to the current time for the current patient:
        % [WAV, BIS, MAP, CO, HR, SV,NMB_m0, NMB_m1, NMB_m2, NMB_m3, cp_prop,
        % cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, ce_rocu, ce_del, ce_wav, ce_bis]
            state = obj.current_patient.get_patient_state_history();
            c_nore = state('c_nore');
            c_nore = c_nore * obj.nore_molweight / 1000; % Norepinephrine blood concentration in [ng/ml]
            state('c_nore') = c_nore;
            patient_state_history = state;
        end

        function patient_input_history = get_patient_input_history(obj)
        % Return a dictionary with the patient input history from the start
        % of the simulation to the current time for the current patient:
            input = obj.current_patient.get_patient_input_history();
            input_nore = input('u_nore');
            input_nore = input_nore * obj.nore_molweight / 1000; % Norepinephrine infusion rate in [µg/s]
            input('u_nore') = input_nore;
            patient_input_history = input;
        end

        function patient_phase = get_patient_phase(obj)
        % Return the current phase of the patient:
            patient_phase = obj.current_patient.get_patient_phase();
        end

        function disturbances = get_patient_disturbances(obj)
        % Return the disturbances of the current patient:
            disturbances = obj.current_patient.get_patient_disturbances();
        end

        function privateProps = getPrivateProperties(obj)
            % Get metadata for the class
            mc = metaclass(obj);
            privateProps = {};
            
            % Loop through properties and find private ones
            for prop = mc.PropertyList'
                if strcmp(prop.GetAccess, 'protected')
                    privateProps{end+1} = prop.Name;
                end
            end
        end

        function patient_results = get_patient_results(obj, num_simulation)
        %Get patient records for the specified simulation.
        % If no simulation index is given, return all the patient records.
        % Parameters:
        %   num_simulation (optional): The index of the simulation to retrieve data from. If not provided, return all data.
        % Returns:
        %   patient_results: A containers.Map with the requested variables or all patient data.

        if nargin < 2 || isempty(num_simulation)
            num_simulation = [];
        end

        if ~isempty(num_simulation)
            if ~(num_simulation >= 1 && num_simulation <= length(obj.patients)) 
                error('IndexError: Simulation index out of range');
            end
        end

        % Get all the attributes ending with '_all' and that are lists and return the data for the given simulation
        patient_results = containers.Map();
        variables = obj.getPrivateProperties();

        for i = 1:length(variables)
            var = variables{i};
            if endsWith(var, '_all')
                if isempty(num_simulation)
                    patient_results(var) = obj.(var);
                else
                    patient_results(var) = obj.(var)(num_simulation, :);
                end
            end
        end
        end


        function initialize_simulation(obj, patient_data, t_sim, t_s, pk_models, pd_models, varargin)
        % Initialize the simulation with the given patient data and simulation parameters.
        % Parameters:
        %   patient_data: patient data [height, weight, age, gender, E0, k_d , delay, C50P, gammaP, rms_nonlin][array]
        %   t_sim: simulation time [double]
        %   t_s: simulation step size [double]
        %   pk_models: pharmacokinetic models for propofol, remifentanil, norepinephrine and rocuronium [map]
        %   pd_models: pharmacodynamic models for propofol and remifentanil [map]
        %   varargin: optional arguments 
        %       opiates: boolean to include opiates [boolean]
        %       blood_sampling: type of blood sampling [string]
        %       internal_states : The internal states of the patient
        %       output_init: The initial values for map, hr, sv
        %       interaction: type of interaction [Interaction]
        %       doh_measure: type of depth of hypnosis measure [DoHMeasure]
        %       stimuli: stimuli for the simulation [map]
        %       volume_status: volume status of the patient [map]

            % Create input parser
            p = inputParser;

            % Required inputs
            addRequired(p, 'patient_data');
            addRequired(p, 't_sim');
            addRequired(p, 't_s'); 
            addRequired(p, 'pk_models');  % This are already checked
            addRequired(p, 'pd_models');  % This are already checked
 
            % Optional inputs with default values
            addOptional(p, 'opiates', []);
            addOptional(p, 'blood_sampling', []);
            addOptional(p, 'internal_states', []);
            addOptional(p, 'interaction', []);
            addOptional(p, 'doh_measure', []);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed_variability', []);
            addOptional(p, 'seed_disturbance', []);
            addOptional(p, 'worst_case', []);
            addOptional(p, 'output_init', []);


            % Parse inputs (handle optional arguments via varargin)
            parse(p, patient_data, t_sim, t_s, pk_models, pd_models, varargin{:});

            patient_data = obj.refactor_patient_data(patient_data);
            
            % Initialize the simulation with given patient data and simulation parameters.
            if ~isempty(p.Results.stimuli)
                disturbance_model = Disturbance(t_sim, p.Results.stimuli,p.Results.seed_disturbance, p.Results.worst_case);
            else
                disturbance_model = [];
            end

            obj.current_patient = Patient(patient_data, ...
                'opiates', p.Results.opiates, ...
                'blood_sampling', p.Results.blood_sampling, ...
                'internal_states', p.Results.internal_states, ...
                'pk_models', p.Results.pk_models, ...
                'pd_models', p.Results.pd_models, ...
                'interaction', p.Results.interaction, ...
                'dohMeasure', p.Results.doh_measure, ...
                'disturbance_model', disturbance_model, ...
                'volume_status', p.Results.volume_status, ...
                'seed_variability', p.Results.seed_variability,...
                'output_init', p.Results.output_init);
            
            obj.patients(end+1) = obj.current_patient;
            
            obj.t_sim = t_sim;
            obj.t_s = t_s;
            obj.steps = int32(obj.t_sim / obj.t_s);
            obj.current_time = 0;

        end

        function init_simulation_from_file(obj, id_patient, t_sim, t_s, varargin)
        % Initialize the simulation using patient data from a file.
        % Parameters:
        %   id_patient: The index of the patient to simulate [double]
        %   t_sim: simulation time [double]
        %   t_s: simulation step size [double]
        %   varargin: optional arguments
        %      opiates: boolean to include opiates [boolean]
        %      blood_sampling: type of blood sampling [string]
        %      internal_states : The internal states of the patient
        %      output_init: The initial values for map, hr, sv
        %      pk_models: pharmacokinetic models for propofol, remifentanil, norepinephrine and rocuronium [map]
        %      pd_models: pharmacodynamic models for propofol and remifentanil [map]
        %      interaction: type of interaction [Interaction]
        %      doh_measure: type of depth of hypnosis measure [DoHMeasure]
        %      stimuli: stimuli for the simulation [map]
        %      volume_status: volume status of the patient [map]
        
            % Create input parser
            p = inputParser;
            
            % Required inputs
            addRequired(p, 'id_patient');
            addRequired(p, 't_sim');
            addRequired(p, 't_s');

            % Optional inputs with default values
            addOptional(p, 'opiates', []);
            addOptional(p, 'blood_sampling', []);
            addOptional(p, 'internal_states', []);
            addOptional(p, 'pk_models',[]);
            addOptional(p, 'pd_models',[]);
            addOptional(p, 'interaction', []);
            addOptional(p, 'doh_measure', []);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed_variability', []);
            addOptional(p, 'seed_disturbance', []);
            addOptional(p, 'worst_case', []);
            addOptional(p, 'output_init', []);

            % Parse inputs, including varargin
            parse(p, id_patient, t_sim, t_s, varargin{:});

            if id_patient < 1 || id_patient > 44
                error('Patient index out of range. Please select a patient between 1 and 10.');
            end

            [pk_models, pd_models] = obj.parse_pk_pd_models(p.Results.pk_models, p.Results.pd_models);

            % Initialize the simulation using patient data from a file.
            if pk_models('prop') == Model.Schnider
                path = 'pd_data_schnider.mat';
                data = load(path);
                patient_data = data.pd_data_schnider(id_patient, :);
            elseif pk_models('prop') == Model.Eleveld
                path = 'pd_data_eleveld.mat';
                data = load(path);
                patient_data = data.pd_data_eleveld(id_patient, :);
            else
                error('Model not supported');
            end

            var_param = {'opiates', p.Results.opiates, 'blood_sampling',...
                p.Results.blood_sampling, 'internal_states', p.Results.internal_states, 'interaction', p.Results.interaction, ...
                'doh_measure', p.Results.doh_measure, 'stimuli', p.Results.stimuli,...
                'volume_status', p.Results.volume_status, 'seed_variability', p.Results.seed_variability, ...
                'seed_disturbance', p.Results.seed_disturbance, 'worst_case', p.Results.worst_case,...
                'output_init', p.Results.output_init};

            obj.initialize_simulation(patient_data, t_sim, t_s, pk_models, pd_models, var_param{:});
        end

        function init_simulation_from_data(obj, data, t_sim, t_s, varargin)
        % Initialize the simulation using patient data from a matrix.
        % Parameters:
        %   data: patient data [height, weight, age, gender, E0, k_d , delay, CE50P, gammaP, rms_nonlin][array]
        %   t_sim: simulation time [double]
        %   t_s: simulation step size [double]
        %   varargin: optional arguments
        %      opiates: boolean to include opiates [boolean]
        %      blood_sampling: type of blood sampling [string]
        %      internal_states : The internal states of the patient
        %      output_init: The initial values for map, hr, sv
        %      pk_models: pharmacokinetic models for propofol, remifentanil, norepinephrine and rocuronium [map]
        %      pd_models: pharmacodynamic models for propofol and remifentanil [map]
        %      interaction: type of interaction [Interaction]
        %      doh_measure: type of depth of hypnosis measure [DoHMeasure]
        %      stimuli: stimuli for the simulation [map]
        %      volume_status: volume status of the patient [map]

            % Create input parser
            p = inputParser;
            
            % Required inputs
            addRequired(p, 'data');
            addRequired(p, 't_sim');
            addRequired(p, 't_s');

            % Optional inputs with default values
            addOptional(p, 'opiates', []);
            addOptional(p, 'blood_sampling', []);
            addOptional(p, 'internal_states', []);
            addOptional(p, 'pk_models', []);
            addOptional(p, 'pd_models', []);
            addOptional(p, 'interaction',[]);
            addOptional(p, 'doh_measure',[]);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed_variability', []);
            addOptional(p, 'seed_disturbance', []);
            addOptional(p, 'worst_case', []);
            addOptional(p, 'output_init', []);


            % Parse inputs
            parse(p, data, t_sim, t_s, varargin{:});

            pk_models, pd_models = obj.parse_pk_pd_models(p.Results.pk_models, p.Results.pd_models);

            var_param = {'opiates', p.Results.opiates, 'blood_sampling',...
                p.Results.blood_sampling, 'internal_states', p.Results.internal_states,...
                'interaction', p.Results.interaction,...
                'doh_measure', p.Results.doh_measure, 'stimuli', p.Results.stimuli, ...
                'volume_status', 'seed_variability', p.Results.seed_variability, ...
                'seed_disturbance', p.Results.seed_disturbance, 'worst_case', p.Results.worst_case,...
                'output_init', p.Results.output_init};
            obj.initialize_simulation(patient_data, t_sim, t_s, pk_models, pd_models, var_param{:});

        end

        function run_complete_simulation(obj, u_prop, u_remi, u_nore, u_rocu)
        % Simulates the patient for the given inputs from the start to the end of the simulation.
        % Parameters: 
        %   u_prop: propofol input [list]       [mg * sec^-1]
        %   u_remi: remifentanil input [list]   [µg * sec^-1]
        %   u_nore: norepinephrine input [list] [µg * sec^-1]
        %   u_rocu: rocuronium input [list]     [mg * sec^-1]

            obj.check_run_inputs(u_prop, u_remi, u_nore, u_rocu);
            
            for time = 1:obj.steps
                obj.one_step_simulation(u_prop(time * obj.t_s), u_remi(time * obj.t_s), u_nore(time * obj.t_s), u_rocu(time * obj.t_s));
            end

        end

        function one_step_simulation(obj, u_prop, u_remi, u_nore, u_rocu)
        % Simulates the patient for one step with the given inputs.
        % Parameters:
        %   u_prop: propofol input [double]       [mg * sec^-1]
        %   u_remi: remifentanil input [double]   [µg * sec^-1]
        %   u_nore: norepinephrine input [double] [µg * sec^-1]
        %   u_rocu: rocuronium input [double]     [mg * sec^-1]

            if isempty(obj.current_patient)
                error('No patient data found. Please add the data for a simulation first.');
            end
            
            u_nore = u_nore * 1000 / obj.nore_molweight; % Norepinephrine infusion rate in [nmol/s]

            obj.current_patient.step(u_prop, u_remi, u_nore, u_rocu, obj.t_s);
            obj.current_time = obj.current_time + obj.t_s;
        end

        function simulation_finished = is_simulation_finished(obj)
        % Check if the simulation is finished.
            simulation_finished = obj.current_time >= obj.t_sim;
        end

        function save_simulation(obj)
        % Save the simulation results to the list of all patients' data.
            if isempty(obj.current_patient)
                error('No patient data found. Please add a simulation first.');
            end

            inputs = obj.get_patient_input_history();
            state = obj.get_patient_state_history();
            obj.add_result_to_all(inputs, state);
        end

        function add_result_to_all(obj, inputs, state)
        % Add the current patient data to the list of all patients' data.

            % Add the input data
            obj.u_prop_all = [obj.u_prop_all, inputs('u_prop')];
            obj.u_remi_all = [obj.u_remi_all, inputs('u_remi')];
            obj.u_nore_all = [obj.u_nore_all, inputs('u_nore')];
            obj.u_rocu_all = [obj.u_rocu_all, inputs('u_rocu')];

            % Add the doh
            obj.wav_all = [obj.wav_all, state('wav')];
            obj.bis_all = [obj.bis_all, state('bis')];

            % Add the Hemodynamic variables
            obj.map_all = [obj.map_all, state('map')];
            obj.co_all = [obj.co_all, state('co')];
            obj.hr_all = [obj.hr_all, state('hr')];
            obj.sv_all = [obj.sv_all, state('sv')];
            obj.tpr_all = [obj.tpr_all, state('tpr')];

            % Add the neuro-muscular blockade 
            obj.nmb_m0_all = [obj.nmb_m0_all, state('nmb_m0')];
            obj.nmb_m1_all = [obj.nmb_m1_all, state('nmb_m1')];
            obj.nmb_m2_all = [obj.nmb_m2_all, state('nmb_m2')];
            obj.nmb_m3_all = [obj.nmb_m3_all, state('nmb_m3')];

            % Add the plasma and effect site concentrations
            obj.cp_prop_all = [obj.cp_prop_all, state('cp_prop')];
            obj.cp_remi_all = [obj.cp_remi_all, state('cp_remi')];

            obj.c_nore_all = [obj.c_nore_all, state('c_nore')];
            obj.cp_rocu_all = [obj.cp_rocu_all, state('cp_rocu')];
            obj.ce_prop_all = [obj.ce_prop_all, state('ce_prop')];
            obj.ce_remi_all = [obj.ce_remi_all, state('ce_remi')];
            obj.ce_rocu_all = [obj.ce_rocu_all, state('ce_rocu')];

            obj.ce_del_all = [obj.ce_del_all, state('ce_del')];
            obj.ce_wav_all = [obj.ce_wav_all, state('ce_wav')];
            obj.ce_bis_all = [obj.ce_bis_all, state('ce_bis')];
        end

        function patient_data = refactor_patient_data(obj, patient_data)
            % Refactor the patient data to include the BMI and the lean body mass.
            height = patient_data(1);
            weight = patient_data(2);
            gender = patient_data(4);  % 0 for male and 1 for female

            bmi = weight / (0.01 * height) ^ 2;
            % fat free mass
            % Janmahasatian S, Duffull SB, Ash S, Ward LC, Byrne NM, Green B.
            % Quantification of lean bodyweight. Clin Pharmacokinet. 2005;44(10):1051-65.
            % doi: 10.2165/00003088-200544100-00004. PMID: 16176118.
            num = 9.27 * 1000 * weight;
            if gender == 0 % male
                den = 6.68 * 1000 + 216 * bmi;
            else % female
                den = 8.78 * 1000 + 244 * bmi;
            end
            lbm = num / den;
            patient_data = [patient_data(1:4), bmi, lbm, patient_data(5:end)];
        end

        function check_run_inputs(obj, u_prop, u_remi, u_nore, u_rocu)
            if isempty(obj.current_patient)
                error('No patient data found. Please add the data for a simulation first.');
            end

            % Check if the input lists have the same length as the simulation time
            % If not, adjust the simulation time and number of steps accordingly
            if length(u_prop) ~= obj.t_sim || length(u_remi) ~= obj.t_sim || length(u_nore) ~= obj.t_sim || length(u_rocu) ~= obj.t_sim
                obj.t_sim = min([length(u_prop), length(u_remi), length(u_nore), length(u_rocu), obj.t_sim]);
                obj.steps = int32(obj.t_sim / obj.t_s);
            end
        end

        function save_to_csv(obj, path, filename)
            % Save the simulation results to a CSV file with columns mapped to new variable names

            % Create the directory if it does not exist
            if ~exist(path, 'dir')
                mkdir(path);
            end

            % Define the mapping for variable names
            columns_mapping = struct(...
                'u_prop_all', 'u_prop', 'u_remi_all', 'u_remi', 'u_nore_all', 'u_nore', ...
                'u_rocu_all', 'u_rocu', 'cp_prop_all', 'cp_prop', 'cp_remi_all', 'cp_remi', ...
                'c_nore_all', 'c_nore', 'cp_rocu_all', 'cp_rocu', 'ce_prop_all', 'ce_prop', ...
                'ce_remi_all', 'ce_remi', 'ce_rocu_all', 'ce_rocu', 'ce_del_all', 'ce_del',...
                'ce_wav_all', 'ce_wav', 'ce_bis_all', 'ce_bis', 'wav_all', 'wav', 'bis_all', 'bis', ...
                'map_all', 'map', 'co_all', 'co', 'hr_all', 'hr', 'sv_all', 'sv', ...
                'nmb_m0_all', 'nmb_m0', 'nmb_m1_all', 'nmb_m1', 'nmb_m2_all', 'nmb_m2', 'nmb_m3_all', 'nmb_m3');

            % Initialize data storage
            data = struct();
            fieldNames = fieldnames(columns_mapping);

            % Extract data and rename columns based on mapping
            for i = 1:numel(fieldNames)
                field = fieldNames{i};
                if isprop(obj, field)  % Check if the property exists
                    variableName = columns_mapping.(field);
                    data.(variableName) = obj.(field);  % Directly store the property as a column
                end
            end

            % Convert the structure to a table
            T = struct2table(data);

            % Save the table to a CSV file
            fullPath = fullfile(path, filename);
            writetable(T, fullPath);

            % Notify the user
            fprintf('Data saved to %s\n', fullPath);
        end

        function plot_simulation(obj,varargin)
        % Plot the simulation results.

            t_input = linspace(0,length(obj.u_prop_all),length(obj.u_prop_all));
            t_input = t_input/60;   % minutes
            t_final = obj.t_sim/60;         % minutes

            %% DoH
            figure
            subplot(4,2,1)
            plot(t_input, obj.u_prop_all)
            ylabel('$u_{prop} [mg/s]$','Interpreter','latex')
            % ylim([0 0.5])
            xlim([0, t_final])
            subplot(4,2,2)
            plot(t_input,  obj.u_remi_all)
            % ylim([0 0.5])
            ylabel('$u_{remi} [\mu g/s]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(4,2,3)
            plot(t_input,  obj.cp_prop_all)
            ylabel('$C_{p, prop} [\mu g/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(4,2,4)
            plot(t_input,  obj.cp_remi_all)
            ylabel('$C_{p,remi} [ng/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(4,2,5)
            plot(t_input,  obj.ce_prop_all)
            ylabel('$C_{e,prop} [\mu g/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(4,2,6)
            plot(t_input,  obj.ce_remi_all)
            ylabel('$C_{e,remi} [ng/ml]$','Interpreter','latex')
            xlim([0, t_final])
            xlabel({'Time (min)'},'Interpreter','latex')
            subplot(4,2,7)
            plot(t_input,  obj.wav_all)
            ylabel('wav$_{cns}$ ','Interpreter','latex')
            xlim([0, t_final])
            subplot(4,2,8)
            plot(t_input,  obj.bis_all)
            ylabel('bis ','Interpreter','latex')
            xlim([0, t_final])
            xlabel({'Time (min)'},'Interpreter','latex')



            %% hemodynamic plots
            figure
            subplot(5,2,1)
            plot(t_input, obj.u_prop_all)
            ylabel('$u_{prop} [mg/s]$','Interpreter','latex')
            % ylim([0 0.5])
            xlim([0, t_final])
            subplot(5,2,2)
            plot(t_input, obj.u_remi_all)
            % ylim([0 0.5])
            ylabel('$u_{remi} [\mu g/s]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,3)
            plot(t_input, obj.cp_prop_all)
            ylabel('$C_{p, prop} [\mu g/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,4)
            plot(t_input, obj.cp_remi_all)
            ylabel('$C_{p,remi} [ng/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,5)
            plot(t_input, obj.u_nore_all)
            ylim([0 1])
            ylabel('$u_{nore} [\mu g/s]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,6)
            plot(t_input, obj.c_nore_all)
            ylabel('$C_{blood,nore} [ng/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(5, 2, 7)
            plot(t_input, obj.map_all)
            ylabel('map [mmHg]','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,8)
            plot(t_input, obj.co_all)
            ylabel('co [L/min]','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,9)
            plot(t_input, obj.hr_all)
            ylabel('hr [beats/min]','Interpreter','latex')
            xlim([0, t_final])
            subplot(5,2,10)
            plot(t_input, obj.sv_all)
            ylabel('sv [ml]','Interpreter','latex')
            xlim([0, t_final])
            xlabel({'Time (min)'},'Interpreter','latex')

            %% nmb
            figure
            subplot(6,1,1)
            plot(t_input, obj.u_rocu_all)
            ylabel('$u_{rocu} [mg/s]$','Interpreter','latex')
            % ylim([0 0.1])
            xlim([0, t_final])
            subplot(6,1,2)
            plot(t_input, obj.cp_rocu_all)
            ylabel('$C_{pRo} [\mu g/ml]$','Interpreter','latex')
            xlim([0, t_final])
            subplot(6, 1, 3)
            plot(t_input, obj.nmb_m0_all)
            ylabel('Probability (nmb)','Interpreter','latex')
            legend('m = 0')
            subplot(6, 1, 4)
            plot(t_input, obj.nmb_m1_all)
            ylabel('Probability (nmb)','Interpreter','latex')
            legend('m = 1')
            subplot(6, 1, 5)
            plot(t_input, obj.nmb_m2_all)
            ylabel('Probability (nmb)','Interpreter','latex')
            legend('m = 2')
            subplot(6, 1, 6)
            plot(t_input, obj.nmb_m3_all)
            ylabel('Probability (nmb)','Interpreter','latex')
            legend( 'm = 3')
            xlim([0, t_final])
            xlabel({'Time (min)'},'Interpreter','latex')  
    end

    function [pk_models, pd_models] = parse_pk_pd_models(~,pk,pd)
    % Parse the pharmacokinetic and pharmacodynamic models from the input arguments.
    % Parameters:
    %   varargin: optional arguments

    pk_models = containers.Map({'prop', 'remi', 'nore', 'rocu'},{Model.Eleveld, Model.Minto, Model.Joachim, Model.Dahe});
    pd_models = containers.Map({'prop', 'remi'},{Model.Eleveld, Model.Minto});

    if ~isempty(pk)
        for i = keys(pk)
            field_name = i{1}; 
            if isKey(pk_models, field_name) % check for wrong field names
                pk_models(field_name) = pk(field_name);
            end
        end
    end

    if ~isempty(pd)
        for i = keys(pd)
            field_name = i{1}; 
            if isKey(pd_models, field_name) % check for wrong field names
                pd_models(field_name) = pd(field_name);
            end
        end
    end
    end

    end 
end
