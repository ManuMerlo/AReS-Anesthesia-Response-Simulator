classdef SimulatorConcentration < Simulator & handle
    %{
    In this class, we use TCI objects to compute infusion rates of each drug based on the target concentration.
    %}
    properties
        tci_prop = [];
        tci_remi = [];
        tci_nore = [];
        tci_rocu = [];
    end

    methods
        function obj = SimulatorConcentration()
        % Call the parent class constructor
            obj = obj@Simulator();
        end
        
        function reset(obj)
        % Reset the simulator to its initial state.
            reset@Simulator(obj);
            obj.tci_prop = [];
            obj.tci_remi = [];
            obj.tci_nore = [];
            obj.tci_rocu = [];
        end

    
        function initialize_simulation(obj, patient_data, t_sim, t_s, varargin)
        % Initialize the simulation with the given parameters.
        % Parameters:
        %   patient_data: patient data [height, weight, age, gender, E0, k_d , delay, C50P, gammaP] [array]
        %   t_sim: simulation time [double]
        %   t_s: sampling time [double]
        %   varargin: optional arguments
        %      opiates: opiates [bool]
        %      blood_sampling: blood sampling [string]
        %      internal_states : The internal states of the patient
        %      output_init: The initial values for map, hr, sv
        %      interaction: interaction [Interaction]
        %      doh_measure: depth of hypnosis measure [DoHMeasure]
        %      stimuli: stimuli [array]
        %      volume_status: volume status [VolumeStatus]
        %      pk_models: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium [Model]
        %      pd_models: pharmacodynamic model for propofol, remifentanil [Model]
        %      pk_models_TCI: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium in TCI [Model]
        %      pd_models_TCI: pharmacodynamic model for propofol, remifentanil in TCI [Model]
        %      patient_data_TCI: patient data for TCI [array]

            
            % Parse input parameters (same as before)
            p = inputParser;
            addRequired(p, 'patient_data');
            addRequired(p, 't_sim');
            addRequired(p, 't_s');
            addOptional(p, 'opiates', []);
            addOptional(p, 'blood_sampling', []);
            addOptional(p, 'internal_states', []);
            addOptional(p, 'interaction', []);
            addOptional(p, 'doh_measure', []);
            addOptional(p, 'pk_models', []);
            addOptional(p, 'pd_models', []);
            addOptional(p, 'pk_models_TCI', []);
            addOptional(p, 'pd_models_TCI', []);
            addOptional(p, 'modes_TCI', []);
            addOptional(p, 'patient_data_TCI', []);
            addOptional(p, 'limits_TCI', []);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed', []);
            addOptional(p, 'output_init', []);
            parse(p, patient_data, t_sim, t_s, varargin{:});

            if isempty(p.Results.patient_data_TCI)
                error('No patient data for TCI found.');
            end
            patient_data_TCI = p.Results.patient_data_TCI;
            patient_data_TCI = refactor_patient_data(obj, patient_data_TCI);

            % Reconstruct varargin for TCI settings
            tci_varargin = {
                'pk_models', p.Results.pk_models_TCI, ...
                'pd_models', p.Results.pd_models_TCI, ...
                'patient_data', patient_data_TCI, ...
                'opiates', p.Results.opiates, ...
                'blood_sampling', p.Results.blood_sampling
            };

            % Default values
            limits = containers.Map({'cp_limit_prop','infusion_limit_prop', ...
                                    'cp_limit_remi','infusion_limit_remi', ...
                                    'cp_limit_nore','infusion_limit_nore', ...
                                    'cp_limit_rocu','infusion_limit_rocu'}, ...
                                    {10, 2, 10, 0.5, 15, 0.5, 10, 0.1});
            

            modes_TCI = containers.Map({'prop','remi','nore','rocu'}, ...
                                    {TciMode.EffectSite, TciMode.EffectSite, TciMode.Plasma, TciMode.EffectSite});

            
            % Check if p.Results.limits_TCI is not empty and update the limits
            if ~isempty(p.Results.limits_TCI)
                for i = keys(p.Results.limits_TCI)
                    field_name = i{1}; 
                    if isKey(limits, field_name) % check for wrong field names
                        limits(field_name) = p.Results.limits_TCI(field_name);
                    end
                end
            end

            if ~isempty(p.Results.modes_TCI)
                for i = keys(p.Results.modes_TCI)
                    field_name = i{1}; 
                    if isKey(modes_TCI, field_name) % check for wrong field names
                        modes_TCI(field_name) = p.Results.modes_TCI(field_name);
                    end
                end
            end

            % Call TCI settings for various drugs
            obj.tci_prop = TCI();
            obj.tci_prop.tci_setting(Drug.Propofol,modes_TCI('prop'), limits('cp_limit_prop'), limits('infusion_limit_prop'), t_s, tci_varargin{:});
            obj.tci_prop.reset_state();

            obj.tci_remi = TCI();
            obj.tci_remi.tci_setting(Drug.Remifentanil, modes_TCI('remi'), limits('cp_limit_remi'), limits('infusion_limit_remi'), t_s, tci_varargin{:});
            obj.tci_remi.reset_state();

            obj.tci_nore = TCI();
            obj.tci_nore.tci_setting(Drug.Norepinephrine, modes_TCI('nore'), limits('cp_limit_nore'), limits('infusion_limit_nore'), t_s, tci_varargin{:});
            obj.tci_nore.reset_state();

            obj.tci_rocu = TCI();
            obj.tci_rocu.tci_setting(Drug.Rocuronium, modes_TCI('rocu'), limits('cp_limit_rocu'), limits('infusion_limit_rocu'), t_s, tci_varargin{:});
            obj.tci_rocu.reset_state();

            % Prepare varargin for parent method
            sim_varargin = {
                'internal_states', p.Results.internal_states, ...
                'interaction', p.Results.interaction, ...
                'doh_measure', p.Results.doh_measure, ...
                'opiates', p.Results.opiates, ...
                'blood_sampling', p.Results.blood_sampling, ...
                'stimuli', p.Results.stimuli, ...
                'volume_status', p.Results.volume_status, ...
                'seed', p.Results.seed,...
                'output_init', p.Results.output_init
            };

            % Call parent method
            initialize_simulation@Simulator(obj, patient_data, t_sim, t_s,p.Results.pk_models,p.Results.pd_models,sim_varargin{:});
        end
        
        function init_simulation_from_file(obj, id_patient, t_sim, t_s, varargin)
        % Initialize the simulation with the given patient ID.
        % Parameters:
        %   id_patient: Id of the patient to simulate [double]
        %   t_sim: simulation time [double]
        %   t_s: sampling time [double]
        %   varargin: optional arguments
        %      opiates: opiates [bool]
        %      blood_sampling: blood sampling [string]
        %      interaction: interaction [Interaction]
        %      doh_measure: depth of hypnosis measure [DoHMeasure]
        %      stimuli: stimuli [array]
        %      volume_status: volume status [VolumeStatus]
        %      pk_models: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium [Model]
        %      pd_models: pharmacodynamic model for propofol, remifentanil [Model]
        %      pk_models_TCI: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium in TCI [Model]
        %      pd_models_TCI: pharmacodynamic model for propofol, remifentanil in TCI [Model]

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
            addOptional(p, 'interaction', []);
            addOptional(p, 'doh_measure', []);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'pk_models', []);
            addOptional(p, 'pd_models', []);
            addOptional(p, 'pk_models_TCI', []);
            addOptional(p, 'pd_models_TCI', []);
            addOptional(p, 'modes_TCI', []);
            addOptional(p, 'patient_data_TCI', []);
            addOptional(p, 'limits_TCI', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed', []);
            addOptional(p, 'output_init', []);

            % Parse inputs, including varargin
            parse(p, id_patient, t_sim, t_s, varargin{:});

            [pk_models, pd_models] = obj.parse_pk_pd_models(p.Results.pk_models,p.Results.pd_models);
            [pk_models_TCI, pd_models_TCI] = obj.parse_pk_pd_models(p.Results.pk_models_TCI, p.Results.pd_models_TCI);
            
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

            if pk_models_TCI('prop') == Model.Schnider
                path = 'pd_data_schnider.mat';
                data = load(path);
                patient_data_TCI = data.pd_data_schnider(id_patient, :);
            elseif pk_models_TCI('prop') == Model.Eleveld
                path = 'pd_data_eleveld.mat';
                data = load(path);
                patient_data_TCI = data.pd_data_eleveld(id_patient, :);
            else
                error('Model not supported');
            end

            var_TCI = {'pk_models_TCI', pk_models_TCI,'pd_models_TCI',pd_models_TCI,...
                'modes_TCI', p.Results.modes_TCI, 'patient_data_TCI', patient_data_TCI, ...
                'limits_TCI', p.Results.limits_TCI};
            var_patient = {'opiates', p.Results.opiates, 'blood_sampling', ...
                p.Results.blood_sampling, 'internal_states', p.Results.internal_states, ...
                'pk_models', pk_models,'pd_models',...
                pd_models,  'interaction', p.Results.interaction, 'doh_measure', ...
                p.Results.doh_measure, 'stimuli', p.Results.stimuli,  ...
                'seed', p.Results.seed, ...
                'output_init', p.Results.output_init};
            obj.initialize_simulation(patient_data, t_sim, t_s, var_patient{:}, var_TCI{:});
        end

        function init_simulation_from_data(obj, data, t_sim, t_s, varargin)
        % Initialize the simulation with the given patient data.
        % Parameters:
        %   data: patient data [height, weight, age, gender, E0, k_d , delay, C50P, gammaP] [array]
        %   t_sim: simulation time [double]
        %   t_s: sampling time [double]
        %   varargin: optional arguments 
        %       opiates: opiates [bool]
        %       blood_sampling: blood sampling [string]
        %       internal_states: internal states [containers.Map]
        %       interaction: interaction [Interaction]
        %       doh_measure: depth of hypnosis measure [DoHMeasure]
        %       stimuli: stimuli [array]
        %       volume_status: volume status [VolumeStatus]
        %       pk_models: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium [Model]
        %       pd_models: pharmacodynamic model for propofol, remifentanil [Model]
        %       pk_models_TCI: pharmacokinetic model for propofol, remifentanil, norepinephrine and rocuronium in TCI [Model]
        %       pd_models_TCI: pharmacodynamic model for propofol, remifentanil in TCI [Model]
        %       patient_data_TCI: patient data for TCI [array]

           
            % Create input parser
            p = inputParser;
            
            % Required inputs
            addRequired(p, 'data');
            addRequired(p, 't_sim');
            addRequired(p, 't_s');

            % Optional inputs with default values
            addOptional(p, 'opiates',[]);
            addOptional(p, 'blood_sampling', []);
            addOptional(p, 'internal_states', []);
            addOptional(p, 'interaction', []);
            addOptional(p, 'doh_measure', []);
            addOptional(p, 'stimuli', []);
            addOptional(p, 'pk_models', []);
            addOptional(p, 'pd_models', []);
            addOptional(p, 'pk_models_TCI', []);
            addOptional(p, 'pd_models_TCI', []);
            addOptional(p, 'modes_TCI', TciMode.EffectSite);
            addOptional(p, 'patient_data_TCI', []);
            addOptional(p, 'limits_TCI', []);
            addOptional(p, 'volume_status',[]);
            addOptional(p, 'seed', []);
            addOptional(p, 'output_init', []);

            % Parse inputs
            parse(p, data, t_sim, t_s, varargin{:});

            [pk_models, pd_models] = obj.parse_pk_pd_models(p.Results.pk_models, p.Results.pd_models);
            [pk_models_TCI, pd_models_TCI] = obj.parse_pk_pd_models(p.Results.pk_models_TCI, p.Results.pd_models_TCI);

            var_TCI = {'pk_models_TCI', pk_models_TCI,'pd_models_TCI', pd_models_TCI, ...
                'modes_TCI', p.Results.modes_TCI, 'patient_data_TCI', patient_data_TCI, ...
                'limits_TCI', p.Results.limits_TCI};
            var_patient = {'opiates', p.Results.opiates, 'blood_sampling',...
                p.Results.blood_sampling, 'pk_models', pk_models,'pd_models',...
                pd_models,  'interaction', p.Results.interaction, 'doh_measure',...
                p.Results.doh_measure, 'stimuli', p.Results.stimuli, ...
                'seed', p.Results.seed, ...
                'output_init', p.Results.output_init};

            obj.initialize_simulation( patient_data, t_sim, t_s, var_patient{:}, var_TCI{:});

        end

        function run_complete_simulation(obj, t_prop, t_remi, t_nore, t_rocu)
        % Run the complete simulation with the given parameters.
        % Parameters:
        %   t_prop: propofol concentrations targets for the TCI in [µg/mL] [array]
        %   t_remi: remifentanil concentrations targets for the TCI in [ng/mL] [array]
        %   t_nore: norepinephrine concentrations targets for the TCI in [ng/mL] [array]
        %   t_rocu: rocuronium concentrations targets for the TCI in [µg/mL] [array]

            obj.check_run_inputs(t_prop, t_remi, t_nore, t_rocu);

            for time = 1:obj.steps
                obj.one_step_simulation(t_prop(time * obj.t_s), t_remi(time * obj.t_s), t_nore(time * obj.t_s), t_rocu(time * obj.t_s));
            end
         end
    
        function one_step_simulation(obj, t_prop, t_remi, t_nore, t_rocu)
        % Run one step of the simulation with the given parameters.
        % Parameters:
        %   t_prop: propofol concentrations targets for the TCI in [µg/mL] [float]
        %   t_remi: remifentanil concentrations targets for the TCI in [ng/mL] [float]
        %   t_nore: norepinephrine concentrations targets for the TCI in [ng/mL] [float]
        %   t_rocu: rocuronium concentrations targets for the TCI in [µg/mL] [float]

            if isempty(obj.current_patient)
                error('No patient data found. Please add the data for a simulation first.');
            end

            t_nore = t_nore * 1000 / obj.nore_molweight; % Norepinephrine plasma concentration in [nmol/L]
             
            obj.tci_prop.compute_infusion(1, t_prop);
            obj.tci_remi.compute_infusion(1, t_remi);
            obj.tci_nore.compute_infusion(1, t_nore);
            obj.tci_rocu.compute_infusion(1, t_rocu);

            u_prop = obj.tci_prop.u;
            u_remi = obj.tci_remi.u;
            u_nore = obj.tci_nore.u;
            u_rocu = obj.tci_rocu.u;
            
            obj.current_patient.step(u_prop, u_remi, u_nore, u_rocu, obj.t_s);
            
            obj.current_time = obj.current_time + obj.t_s;

        end
    end
end
