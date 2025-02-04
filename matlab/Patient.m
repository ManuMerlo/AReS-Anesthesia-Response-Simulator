classdef Patient < handle
    properties (Access = private)
        height
        weight
        age
        gender
        bmi
        lbm
        
        interaction
        dohMeasure
        time = 0

        pk_prop
        pk_remi
        pk_nore
        pk_rocu

        pd_doh
        pd_hemo
        pd_nmb

        disturbance_model
        doh_dis = zeros(1, 1)
        co_dis = zeros(1, 1)

        patient_phase = PatientPhase.Induction
        maintenance_time = 0
        steady_state = false
        steady_state_values = []

        volume_status = []
        volume_status_coeff

        x_0_prop = zeros(3,1) 
        ce_0_prop = 0
        x_wav_filtered0 = zeros(2,1)
        x_bis_lti_0 = zeros(4,1)
        x_bis_delay_0 = 0
        x_0_delay = 0

        x_0_remi = zeros(3,1)
        ce_0_remi = 0
        
        x_0_nore = zeros(3,1)
        x0_nore_delayed = 0
        x_0_rocu = zeros(3,1)
        ce_0_rocu = 0
        
        hemo_init

        cp_prop_arr = []
        ce_prop_arr = []
        ce_del_arr = []
        ce_wav_arr = []
        ce_bis_arr = []

        cp_remi_arr = []
        ce_remi_arr = []

        c_nore_arr = []
        cp_rocu_arr = []

        u_prop_arr = []
        u_remi_arr = []
        u_nore_arr = []
        u_rocu_arr = []

        wav = []
        bis = []

        co = []
        map = []
        hr = []
        sv = []

        nmb_m0 = []
        nmb_m1 = []
        nmb_m2 = []
        nmb_m3 = []

    end

    methods (Static)
        function val = get_last_value_or_default(arr, default)
        % Helper function to get the last value from an array or a default value if the array is empty
            if ~isempty(arr)
                val = arr(end);
            else
                val = default;
            end
        end
    end

    methods
        function obj = Patient(data, varargin)
        % Constructor for the Patient class

            % Create an input parser
            p = inputParser;

            % Required arguments
            addRequired(p, 'data');

            % Optional parameters 
            addOptional(p, 'internal_states', []);
            addOptional(p, 'output_init', []);
            addOptional(p, 'pk_models',[]);
            addOptional(p, 'pd_models',[]);
            addOptional(p, 'opiates', true);
            addOptional(p, 'blood_sampling', 'arterial');
            addOptional(p, 'interaction', Interaction.Surface);
            addOptional(p, 'dohMeasure', DoHMeasure.Both);
            addOptional(p, 'disturbance_model', []);
            addOptional(p, 'volume_status', []);
            addOptional(p, 'seed', []);

            % Parse the inputs
            parse(p, data, varargin{:});

            % Patient demographic parameters
            obj.height = data(1);
            obj.weight = data(2);
            obj.age = data(3);
            obj.gender = data(4); % 1 for female and 0 for male
            obj.bmi = data(5);
            obj.lbm = data(6);

            opiates = p.Results.opiates; % If opiates are used or not during the surgery
            blood_sampling = p.Results.blood_sampling; % The blood sampling site: arterial or venous
            
            % Default internal states
            internal_states = containers.Map({'x_0_prop','ce_0_prop','x_0_remi','ce_0_remi','x_0_nore','x_0_nore_delayed','x_0_rocu','ce_0_rocu'}, {zeros(3,1), 0, zeros(3,1), 0, zeros(3,1), 0, zeros(3,1), 0});
            
            % # Initialize the internal states of the patient
            if ~isempty(p.Results.internal_states)
                for key = keys(internal_states)
                    if isKey(p.Results.internal_states, key)
                        internal_states(key) = p.Results.internal_states(key);
                    end
                end
            end

            obj.x_0_prop = internal_states('x_0_prop');
            obj.ce_0_prop = internal_states('ce_0_prop');
            obj.x_0_remi = internal_states('x_0_remi');
            obj.ce_0_remi = internal_states('ce_0_remi');
            obj.x_0_nore = internal_states('x_0_nore');
            obj.x0_nore_delayed = internal_states('x_0_nore_delayed');
            obj.x_0_rocu = internal_states('x_0_rocu');
            obj.ce_0_rocu = internal_states('ce_0_rocu');

            % Pk-PD models names for propofol, remifentanil, norepinephrine, and rocuronium
            % Note: this values are checked in the Simulators class
            pk_models = p.Results.pk_models;
            pd_models = p.Results.pd_models;

            var_pk = {obj.age, obj.weight, obj.height, obj.gender, obj.bmi, obj.lbm, opiates, blood_sampling};
            var_pd = {obj.age, obj.weight, blood_sampling};

            if pd_models('prop') == Model.Wav
                % Parameters to define the Propofol PD model
                e0 = data(7);
                ke0_prop = data(8);
                delay = data(9);
                ec50_prop = 0.01 * data(7) * data(10);
                gammap = data(11);
                pd_data = [e0, ke0_prop, delay, ec50_prop, gammap];
                var_pd = [var_pd, pd_data];
            end

            % Pharmacokinetic models
            obj.pk_prop = PharmacokineticModel;
            obj.pk_prop = obj.pk_prop.pk_model(Drug.Propofol, pk_models('prop'), var_pk{:});

            obj.pk_remi = PharmacokineticModel;
            obj.pk_remi = obj.pk_remi.pk_model(Drug.Remifentanil, pk_models('remi'), var_pk{:});

            obj.pk_nore = PharmacokineticModel;
            obj.pk_nore = obj.pk_nore.pk_model(Drug.Norepinephrine, pk_models('nore'), var_pk{:});

            obj.pk_rocu = PharmacokineticModel;
            obj.pk_rocu = obj.pk_rocu.pk_model(Drug.Rocuronium, pk_models('rocu'), var_pk{:});

            % Pharmacodynamic models
            obj.pd_doh = PharmacodynamicDoH;
            obj.pd_doh = obj.pd_doh.pd_model_prop(pd_models('prop'), var_pd{:});
            obj.pd_doh = obj.pd_doh.pd_model_remi(pd_models('remi'), var_pd{:});
            if ~isempty(p.Results.output_init) && isKey(p.Results.output_init, 'doh')
                obj.pd_doh.e0 = p.Results.output_init('doh');
            end

            % Initialize the hemodynamic model
            obj.pd_hemo = PharmacodynamicHemodynamics(obj.age);

            % Initialize the hemodyanmic variables.
            obj.hemo_init(4) = 0;
            obj.hemo_init(5) = 0;
            
            if ~isempty(p.Results.output_init)
                output_init = p.Results.output_init;
                
                if isKey(output_init, 'hr')
                    hr = output_init('hr');
                else
                    hr = obj.pd_hemo.base_hr;
                end
                base_hr = hr - obj.pd_hemo.base_hr * obj.pd_hemo.ltde_hr;

                if isKey(output_init, 'sv')
                    sv = output_init('sv');
                else
                    sv = obj.pd_hemo.base_sv;
                end

                base_sv = sv / (1 - obj.pd_hemo.hr_sv * log(hr / (obj.pd_hemo.base_hr * (1 + obj.pd_hemo.ltde_hr)))) - obj.pd_hemo.base_sv * obj.pd_hemo.ltde_sv;

                if isKey(output_init, 'co')
                    co = output_init('co');
                    if ~(co == base_hr * base_sv / 1000)
                        warning('The initial values of CO, HR, and SV are not consistent.')
                    end
                else
                    co = obj.pd_hemo.base_sv * obj.pd_hemo.base_hr / 1000;
                end

                if isKey(output_init, 'map')
                    map = output_init('map');
                else
                    map = obj.pd_hemo.base_sv * obj.pd_hemo.base_hr * obj.pd_hemo.base_tpr;
                end

                base_tpr = map / (co * 1000);

                obj.hemo_init(1) = base_tpr;
                obj.hemo_init(2) = base_sv;
                obj.hemo_init(3) = base_hr;

            else
                %  # Add interpatient variability to the hemodynamic variables if seed is not None
                if ~isempty(p.Results.seed)
                    obj.pd_hemo = obj.pd_hemo.interpatient_variability(p.Results.seed);
                end
                obj.hemo_init(1) = obj.pd_hemo.base_tpr;
                obj.hemo_init(2) = obj.pd_hemo.base_sv;
                obj.hemo_init(3) = obj.pd_hemo.base_hr;
            end

            % Initialize the neuromuscular blockade model
            obj.pd_nmb = PharmacodynamicNMB;
            obj.pd_nmb = obj.pd_nmb.pd_model(Drug.Rocuronium);

            %  Patient simulation settings
            obj.interaction = p.Results.interaction;
            obj.dohMeasure = p.Results.dohMeasure;

            % Disturbance model
            obj.disturbance_model = p.Results.disturbance_model;

            % Volume status
            obj.volume_status = p.Results.volume_status;
            obj.volume_status_coeff = containers.Map({'co', 'hr', 'sv', 'map'}, {1, 1, 1, 1});

            % Note: the other variables are initialized in the defition of the properties

        end

        function patient_data = get_patient_demographics(obj)
        % Get the patient demographics as a map (age, height, weight, gender, bmi, lbm)
        % Returns:
        %   patient_data: A map containing the patient demographics  [containers.Map]
    
            keys = {'age', 'height', 'weight', 'gender', 'bmi', 'lbm'};
            values = {obj.age, obj.height, obj.weight, obj.gender, obj.bmi, obj.lbm};
            patient_data = containers.Map(keys, values);
        end

        function patient_state = get_patient_state(obj)
        % Get the current patient state
        % Returns:
        %   patient_state: A map containing the patient current state [containers.Map]  
        
            % Get the base values
            base_wav = obj.pd_doh.e0;
            base_bis = obj.pd_doh.e0;
            base_co = obj.pd_hemo.base_sv * obj.pd_hemo.base_hr / 1000;
            base_map = obj.pd_hemo.base_sv * obj.pd_hemo.base_hr * obj.pd_hemo.base_tpr;
            base_hr = obj.pd_hemo.base_hr;
            base_sv = obj.pd_hemo.base_sv;

            % Use computed values if available, otherwise use base values
            wav_ = obj.get_last_value_or_default(obj.wav, base_wav);
            bis_ = obj.get_last_value_or_default(obj.bis, base_bis);
            map_ = obj.get_last_value_or_default(obj.map, base_map);
            co_ = obj.get_last_value_or_default(obj.co, base_co);
            hr_ = obj.get_last_value_or_default(obj.hr, base_hr);
            sv_ = obj.get_last_value_or_default(obj.sv, base_sv);
            nmb_m0_ = obj.get_last_value_or_default(obj.nmb_m0, 1);
            nmb_m1_ = obj.get_last_value_or_default(obj.nmb_m1, 0);
            nmb_m2_ = obj.get_last_value_or_default(obj.nmb_m2, 0);
            nmb_m3_ = obj.get_last_value_or_default(obj.nmb_m3, 0);
            cp_prop = obj.get_last_value_or_default(obj.cp_prop_arr, 0);
            cp_remi = obj.get_last_value_or_default(obj.cp_remi_arr, 0);
            c_nore = obj.get_last_value_or_default(obj.c_nore_arr, 0);
            cp_rocu = obj.get_last_value_or_default(obj.cp_rocu_arr, 0);
            ce_prop = obj.get_last_value_or_default(obj.ce_prop_arr, 0);
            ce_remi = obj.get_last_value_or_default(obj.ce_remi_arr, 0);
            ce_del = obj.get_last_value_or_default(obj.ce_del_arr, 0);
            ce_wav = obj.get_last_value_or_default(obj.ce_wav_arr, 0);
            ce_bis = obj.get_last_value_or_default(obj.ce_bis_arr, 0);

            keys = {'wav', 'bis', 'map', 'co', 'hr', 'sv','nmb_m0', 'nmb_m1', 'nmb_m2', 'nmb_m3','cp_prop', 'cp_remi', 'c_nore', 'cp_rocu', 'ce_prop', 'ce_remi', 'ce_del', 'ce_wav', 'ce_bis'};
            values = {wav_, bis_, map_, co_, hr_, sv_, nmb_m0_, nmb_m1_, nmb_m2_, nmb_m3_, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, ce_del, ce_wav, ce_bis};
            patient_state = containers.Map(keys, values);
        end

        function patient_state_history = get_patient_state_history(obj)
        % Get the patient state history
        % Returns:
        %   patient_state_history: A map containing the patient state history from the  beginning of the simulation [containers.Map]
            keys = {'wav', 'bis', 'map', 'co', 'hr', 'sv','nmb_m0', 'nmb_m1', 'nmb_m2', 'nmb_m3', 'cp_prop', 'cp_remi', 'c_nore', 'cp_rocu', 'ce_prop', 'ce_remi', 'ce_del', 'ce_wav', 'ce_bis'};
            values = { obj.wav, obj.bis, obj.map, obj.co, obj.hr, obj.sv, obj.nmb_m0, obj.nmb_m1, obj.nmb_m2, obj.nmb_m3, obj.cp_prop_arr, obj.cp_remi_arr, obj.c_nore_arr, obj.cp_rocu_arr, obj.ce_prop_arr, obj.ce_remi_arr, obj.ce_del_arr, obj.ce_wav_arr, obj.ce_bis_arr};
            patient_state_history = containers.Map(keys, values);
        end

        function patient_input_history = get_patient_input_history(obj)
        % Get the patient input history
        % Returns:
        %   patient_input_history: A map containing the patient input history from the beginning of the simulation [containers.Map]

            keys = {'u_prop', 'u_remi', 'u_nore', 'u_rocu'};
            values = {obj.u_prop_arr, obj.u_remi_arr, obj.u_nore_arr, obj.u_rocu_arr};
            patient_input_history = containers.Map(keys, values);
        end

        function disturbances = get_patient_disturbances(obj)
        % Get the patient disturbances
        % Returns:  
        %   disturbances: A map containing the patient disturbances [containers.Map]
            disturbances = containers.Map({'doh_dis', 'co_dis'}, {obj.doh_dis, obj.co_dis});
        end

        function [phase, steady_state] = get_patient_phase(obj)
        % Get the current patient phase and the steady state status
        % Returns:
        %   phase: The current patient phase [PatientPhase]
            phase = obj.patient_phase;
            steady_state = obj.steady_state;
        end

        function interal_state = get_patient_internal_states(obj)
        % Get the current patient internal states
        % Returns:
        %   internal_state: A map containing the patient internal states [containers.Map]
            keys = {'x_0_prop', 'ce_0_prop', 'x_0_remi', 'ce_0_remi', 'x_0_nore', 'x_0_nore_delayed', 'x_0_rocu', 'ce_0_rocu'};
            values = {obj.x_0_prop, obj.ce_0_prop, obj.x_0_remi, obj.ce_0_remi, obj.x_0_nore, obj.x0_nore_delayed, obj.x_0_rocu, obj.ce_0_rocu};
            interal_state = containers.Map(keys, values);
        end

        function set_time(obj, time)
        % Set the current time of the simulation
            obj.time = time;
        end

        function set_disturbance_model(obj, disturbance_model)
        % Set the disturbance model
        % Parameters:
        %   disturbance_model: The disturbance model [DisturbanceModel]
            obj.disturbance_model = disturbance_model;
        end

        function set_volume_status(obj, volume_status)
        % Set the volume status
        % Parameters:
        %   volume_status: The volume status [VolumeStatus]
            obj.volume_status = volume_status;
        end

        function clear_data(obj)
        % Clear the patient data except for the last values. Needed to start the simulation in maintenance phase.
            obj.cp_prop_arr = obj.cp_prop_arr(end);
            obj.ce_prop_arr = obj.ce_prop_arr(end);
            obj.ce_del_arr = obj.ce_del_arr(end);
            obj.ce_wav_arr = obj.ce_wav_arr(end);
            obj.ce_bis_arr = obj.ce_bis_arr(end);

            obj.cp_remi_arr = obj.cp_remi_arr(end);
            obj.ce_remi_arr = obj.ce_remi_arr(end);

            obj.cp_rocu_arr = obj.cp_rocu_arr(end);
            obj.c_nore_arr = obj.c_nore_arr(end);

            obj.u_prop_arr = obj.u_prop_arr(end);
            obj.u_remi_arr = obj.u_remi_arr(end);
            obj.u_nore_arr = obj.u_nore_arr(end);

            obj.wav = obj.wav(end);
            obj.bis = obj.bis(end);

            obj.co = obj.co(end);
            obj.map = obj.map(end);
            obj.hr = obj.hr(end);
            obj.sv = obj.sv(end);

            obj.nmb_m0 = obj.nmb_m0(end);
            obj.nmb_m1 = obj.nmb_m1(end);
            obj.nmb_m2 = obj.nmb_m2(end);
            obj.nmb_m3 = obj.nmb_m3(end);

        end

        function [cp_prop_sim, cp_remi_sim, c_nore_sim, cp_rocu_sim] = compute_plasma_concentration(obj, u_prop, u_remi, u_nore, u_rocu, t)
        % Compute the plasma concentration of the drugs for a given time step given the infusion rates
        % Parameters:
        %   u_prop: Propofol infusion rate [double]
        %   u_remi: Remifentanil infusion rate [double]
        %   u_nore: Norepinephrine infusion rate [double]
        %   u_rocu: Rocuronium infusion rate [double]
        %   t: Time step [array]

            [cp_prop_sim, ~, x_prop] = lsim(obj.pk_prop.pk,u_prop,t,obj.x_0_prop);   
            [cp_remi_sim, ~, x_remi] = lsim(obj.pk_remi.pk,u_remi,t,obj.x_0_remi);
            
            [c_nore_sim, ~, x_nore] = lsim(obj.pk_nore.pk, u_nore, t, obj.x_0_nore);
            [cp_rocu_sim, ~, x_rocu] = lsim(obj.pk_rocu.pk, u_rocu, t, obj.x_0_rocu);

            cp_prop_sim(cp_prop_sim<0) = 0;
            cp_remi_sim(cp_remi_sim<0) = 0;
            c_nore_sim(c_nore_sim<0) = 0;
            cp_rocu_sim(cp_rocu_sim<0) = 0;

            % Update initial states
            obj.x_0_prop = x_prop(end,:);
            obj.x_0_remi = x_remi(end,:);
            obj.x_0_nore = x_nore(end,:);
            obj.x_0_rocu = x_rocu(end,:);

            obj.cp_prop_arr = [obj.cp_prop_arr; cp_prop_sim(1:end)];
            obj.cp_remi_arr = [obj.cp_remi_arr; cp_remi_sim(1:end)];
            obj.c_nore_arr = [obj.c_nore_arr; c_nore_sim(1:end)];
            obj.cp_rocu_arr = [obj.cp_rocu_arr; cp_rocu_sim(1:end)];

        end

        function wav_interv = compute_WAV(obj, ce_delayed, ce_remi_sim, t)
        % Compute the WAV value for a given time step
        % Parameters:
        %   ce_delayed: Delayed effect site concentration [double]
        %   ce_remi_sim: Remifentanil effect site concentration [double]
        %   t: Time step [array]

            [ce_wav_filtered,~, x_wav_filtered] = lsim(obj.pd_doh.neurow_sensor_dynamics(),ce_delayed,t,obj.x_wav_filtered0); 

            ce_wav_filtered(ce_wav_filtered<0) = 0;
            obj.x_wav_filtered0 = x_wav_filtered(end, :);
            obj.ce_wav_arr = [obj.ce_wav_arr; ce_wav_filtered(1:end)];
        
            if obj.interaction == Interaction.Surface
                wav_interv = obj.pd_doh.responseSurfaceModel(ce_wav_filtered, ce_remi_sim);
            elseif obj.interaction == Interaction.No_interaction
                wav_interv = obj.pd_doh.hillfun(ce_wav_filtered);
            else
                error('Interaction mode not available: %s',inter_mode)
            end

        end

        function bis_interv = compute_BIS(obj, ce_delayed, ce_remi_sim, t)
        % Compute the BIS value for a given time step
        % Parameters:
        %   ce_delayed: Delayed effect site concentration [array]
        %   ce_remi_sim: Remifentanil effect site concentration [array]
        %   t: Time step [array]

            [bis_delay, bis_lti] = obj.pd_doh.bis_sensor_dynamics();

            [ce_bis_filtered, ~, x_bis_lti] = lsim(bis_lti,ce_delayed,t,obj.x_bis_lti_0);
            [ce_bis_filtered,~, x_bis_delay] = lsim(bis_delay, ce_bis_filtered,t,obj.x_bis_delay_0);     

            ce_bis_filtered(ce_bis_filtered<0) = 0;
            obj.ce_bis_arr = [obj.ce_bis_arr; ce_bis_filtered(1:end)];

            obj.x_bis_lti_0 = x_bis_lti(end, :);
            obj.x_bis_delay_0 = x_bis_delay(end, :);
            
            if obj.interaction == Interaction.Surface
                bis_interv = obj.pd_doh.responseSurfaceModel(ce_bis_filtered, ce_remi_sim);
            elseif obj.interaction == Interaction.No_interaction
                bis_interv = obj.pd_doh.hillfun(ce_bis_filtered);
            else
                error('Interaction mode not available: %s',inter_mode)
            end
        end

        function [co_interv, map_interv, hr_interv, sv_interv] = compute_hemodynamic_variables(obj, cp_prop, cp_remi, c_nore, t_sim)
        % Compute the hemodynamic variables for a given time step
        % Parameters:
        %   cp_prop: Propofol plasma concentration [array]
        %   cp_remi: Remifentanil plasma concentration [array]
        %   c_nore: Norepinephrine plasma concentration [array]
        %   t_sim: Time step [array]

            % Norepinephrine's effect on hemodynamic variables
            [c_nore_delayed, ~, x_nore_delayed] = lsim(obj.pd_hemo.delay_ss(), c_nore, t_sim, obj.x0_nore_delayed);
            
            c_nore_delayed(c_nore_delayed<0) = 0;
            obj.x0_nore_delayed = x_nore_delayed(end,:);

            map_nore = obj.pd_hemo.hillfun(c_nore_delayed, 'map');
            co_nore = obj.pd_hemo.hillfun(c_nore_delayed, 'co');

            % General pharmacodynamic interaction (GPDI) model between Propofol and Remifentanil
            % Initialize the hemodynamic parameters
            tpr_interv = zeros(length(t_sim),1);
            sv_star_interv = zeros(length(t_sim),1);
            hr_star_interv = zeros(length(t_sim),1);
            tde_hr_interv = zeros(length(t_sim),1);
            tde_sv_interv = zeros(length(t_sim),1);

            tpr_interv(1) = obj.hemo_init(1);
            sv_star_interv(1) = obj.hemo_init(2);
            hr_star_interv(1) = obj.hemo_init(3);
            tde_hr_interv(1) = obj.hemo_init(4) + obj.pd_hemo.base_hr*obj.pd_hemo.ltde_hr;
            tde_sv_interv(1) = obj.hemo_init(5) + obj.pd_hemo.base_sv*obj.pd_hemo.ltde_sv;

            % Solving the nonlinear state space model
            for i=1:length(t_sim)-1
                tspan = [t_sim(i) t_sim(i+1)];
                ic = obj.hemo_init;
                opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
                [~, y_ode] = ode45(@(t, y) obj.pd_hemo.ode_prop_hemodynamic(t, y, cp_prop(i), cp_remi(i)),...
                    tspan, ic, opts);
                tde_hr = y_ode(end, 4) + obj.pd_hemo.base_hr*obj.pd_hemo.ltde_hr;
                tde_sv = y_ode(end, 5) + obj.pd_hemo.base_sv*obj.pd_hemo.ltde_sv;

                tpr_interv(i+1) = y_ode(end, 1);
                sv_star_interv(i+1) = y_ode(end, 2);
                hr_star_interv(i+1) = y_ode(end, 3);
                tde_hr_interv(i+1) = tde_hr;
                tde_sv_interv(i+1) = tde_sv;
                obj.hemo_init = y_ode(end,:);
            end

            hr_interv = hr_star_interv + tde_hr_interv;
            sv_interv = (sv_star_interv+tde_sv_interv).*(1 - obj.pd_hemo.hr_sv*log(hr_interv/(obj.pd_hemo.base_hr*(1+obj.pd_hemo.ltde_hr))));
            co_interv = hr_interv.*sv_interv/1000;

            % Adding the effect of disturbances
            if ~isempty(obj.disturbance_model)
                interval = min(t_sim, length(obj.co_dis) - obj.time);
                co_interv = co_interv + obj.co_dis(obj.time + 1 : obj.time + 1 + interval);
            end
            
            map_interv = tpr_interv.*co_interv*1000;
            hr_interv = co_interv./ sv_interv * 1000;
            sv_interv = co_interv*1000./hr_interv;
            
            % Linear interaction with norepinephrine
            co_interv = co_interv + co_nore;
            map_interv = map_interv + map_nore;

            if ~isempty(obj.volume_status) && isKey(obj.volume_status, obj.time)
                status = obj.volume_status(obj.time);
                if status == VolumeStatus.Hypovolemia
                    obj.volume_status_coeff = containers.Map({'sv', 'co', 'map', 'hr'}, {0.74, 0.868, 1.04,1.16});
                elseif status == VolumeStatus.Normovolemia
                    obj.volume_status_coeff = containers.Map({'sv', 'co', 'map', 'hr'}, {1, 1, 1, 1});
                end
            end

            co_interv = co_interv * obj.volume_status_coeff('co');
            map_interv = map_interv * obj.volume_status_coeff('map');
            hr_interv = hr_interv * obj.volume_status_coeff('hr');
            sv_interv = sv_interv * obj.volume_status_coeff('sv');
        end

        function nmb_interv = compute_nmb(obj, cp_rocu, t)
        % Compute the neuromuscular blockade for a given time step
        % Parameters:
        %   cp_rocu: Rocuronium plasma concentration [array]
        %   t: Time step [array]
        
            [ce_rocu, ~, x_ce_rocu] = lsim(obj.pd_nmb.pd_ce, cp_rocu, t, obj.ce_0_rocu);  

            ce_rocu(ce_rocu<0) = 0;
            obj.ce_0_rocu = x_ce_rocu(end);
            
            % Compute the probability for each category of NMB
            nmb_interv = obj.pd_nmb.hillfun(ce_rocu);
            
        end

        function step(obj, u_prop, u_remi, u_nore, u_rocu, t_s)
        % Simulate the patient for a given time step
        % Parameters:
        %   u_prop: Propofol infusion rate [double]
        %   u_remi: Remifentanil infusion rate [double]
        %   u_nore: Norepinephrine infusion rate [double]
        %   u_rocu: Rocuronium infusion rate [double]
        %   t_s: Time step [double]

            t_sim = 0:1:t_s-1;
            % Generate constant input arrays for the simulation that lasts t_s seconds
            % One step of the simulation is t_s seconds

            u_sim_prop = u_prop*ones(length(t_sim),1);
            u_sim_remi = u_remi*ones(length(t_sim),1);
            u_sim_nore = u_nore*ones(length(t_sim),1);
            u_sim_rocu = u_rocu*ones(length(t_sim),1);

            obj.u_prop_arr = [obj.u_prop_arr; u_sim_prop(1:end)];
            obj.u_remi_arr = [obj.u_remi_arr; u_sim_remi(1:end)];
            obj.u_nore_arr = [obj.u_nore_arr; u_sim_nore(1:end)];
            obj.u_rocu_arr = [obj.u_rocu_arr; u_sim_rocu(1:end)];

            % Compute plasma concentrations

            [cp_prop_sim, cp_remi_sim, c_nore_sim, cp_rocu_sim] = obj.compute_plasma_concentration(u_sim_prop, u_sim_remi, u_sim_nore, u_sim_rocu, t_sim);
            
            % Compute the disturbances caused by stimulations which also affected by
            % the Propofol and Remifentnail plasma concentrations

            if ~isempty(obj.disturbance_model)
                [obj.doh_dis, obj.co_dis] = obj.disturbance_model.get_disturbances(obj.time, cp_prop_sim(end), cp_remi_sim(end));
            end

            % Compute effect-site concentrations of Propofol and  Remifentanil

            [ce_prop_sim, ~, x_ce_sim] = lsim(obj.pd_doh.pd_prop_ce, cp_prop_sim, t_sim, obj.ce_0_prop);              
            [ce_remi_sim, ~, x_ce_remi] = lsim(obj.pd_doh.pd_remi_ce, cp_remi_sim, t_sim, obj.ce_0_remi);  

            ce_prop_sim(ce_prop_sim<0) = 0;
            ce_remi_sim(ce_remi_sim<0) = 0;
            
            % Compute depth of hypnosis
            if obj.pd_doh.delay == 0
                ce_delayed = ce_prop_sim;
                x_delay = x_ce_sim;
            else
                [ce_delayed,~,x_delay] = lsim(obj.pd_doh.delay_ss(), ce_prop_sim,t_sim, obj.x_0_delay);  
            end
            ce_delayed(ce_delayed<0) = 0;

            obj.ce_0_prop = x_ce_sim(end);
            obj.ce_0_remi = x_ce_remi(end);
            obj.x_0_delay = x_delay(end,:);

            obj.ce_prop_arr = [obj.ce_prop_arr; ce_prop_sim(1:end)];
            obj.ce_remi_arr = [obj.ce_remi_arr; ce_remi_sim(1:end)];
            obj.ce_del_arr = [obj.ce_del_arr; ce_delayed(1:end)];


            if obj.dohMeasure == DoHMeasure.Wav || obj.dohMeasure == DoHMeasure.Both
                [wav_interv] = obj.compute_WAV(ce_delayed, ce_remi_sim, t_sim);
                if ~isempty(obj.disturbance_model)
                    interval = min(t_s, length(obj.doh_dis) - obj.time);
                    wav_interv = wav_interv + obj.doh_dis(obj.time + 1: obj.time + 1 + interval);
                    wav_interv = max(0, min(100, wav_interv));
                end
                elems = wav_interv;
            else 
                wav_interv = zeros(length(t_sim),1);
            end

            obj.wav = [obj.wav; wav_interv(1:end)];

            if obj.dohMeasure == DoHMeasure.Bis || obj.dohMeasure == DoHMeasure.Both
                [bis_interv] = obj.compute_BIS(ce_delayed, ce_remi_sim, t_sim);
                if ~isempty(obj.disturbance_model)
                    [dbis, ~ ]= obj.pd_doh.bis_sensor_dynamics();
                    t = 0:1:size(obj.doh_dis)-1;
                    [delayed_disturb,~,~ ] = lsim(dbis, obj.doh_dis, t);
                    interval = min(t_s, length(delayed_disturb) - obj.time);
                    bis_interv = bis_interv + delayed_disturb(obj.time + 1: obj.time + 1 +interval);
                    bis_interv = max(0, min(100, bis_interv));
                end
                elems = bis_interv;
            else
                bis_interv = zeros(length(t_sim),1);
            end
            obj.bis = [obj.bis; bis_interv(1:end)];
            
            if obj.patient_phase == PatientPhase.Induction
                for elem = elems'
                    if elem < 60
                        obj.maintenance_time = obj.maintenance_time + 1;
                    else 
                        obj.maintenance_time = 0;
                    end
                    if obj.maintenance_time > 180
                        obj.patient_phase = PatientPhase.Maintenance;
                        disp('Patient has entered the maintenance phase');
                        break;
                    end
                end
            end
                
            if obj.patient_phase == PatientPhase.Maintenance && obj.time > 1000 && ~obj.steady_state
                lower_bound = 50 - 50 * 0.05;
                upper_bound = 50 + 50 * 0.05;
                for elem = elems'
                    if isempty(obj.steady_state_values)
                        obj.steady_state_values = elem;
                    else
                        if abs(elem - obj.steady_state_values(end)) < 0.2 && elem > lower_bound && elem < upper_bound
                            obj.steady_state_values = [obj.steady_state_values, elem]; 
                        else
                            obj.steady_state_values = []; % Clear the array
                        end
                    end
            
                    % Check if we have at least 180 steady state values
                    if length(obj.steady_state_values) >= 180
                        if abs(obj.steady_state_values(end - 179) - obj.steady_state_values(end)) < 0.5
                            obj.steady_state = true;
                        end
                    end
                end
            end


            % Compute the hemodynamic variables
            [co_interv, map_interv, hr_interv, sv_interv] = obj.compute_hemodynamic_variables(cp_prop_sim, cp_remi_sim, c_nore_sim, t_sim);

            % Compute the probabilities of the NMB categories
            nmb_interv = obj.compute_nmb(cp_rocu_sim, t_sim);
            
            % Save the stats in the recording arrays
            obj.co = [obj.co; co_interv(1:end)];
            obj.map = [obj.map; map_interv(1:end)];
            obj.hr = [obj.hr; hr_interv(1:end)];
            obj.sv = [obj.sv; sv_interv(1:end)];
            
            obj.nmb_m0 = [obj.nmb_m0; nmb_interv(1:end,1)];
            obj.nmb_m1 = [obj.nmb_m1; nmb_interv(1:end,2)];
            obj.nmb_m2 = [obj.nmb_m2; nmb_interv(1:end,3)];
            obj.nmb_m3 = [obj.nmb_m3; nmb_interv(1:end,4)];
 
            obj.time = obj.time + t_s;

        end

    end
    
end 

        