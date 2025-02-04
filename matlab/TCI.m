classdef TCI < handle
    properties
        drug;
        mode = TciMode.EffectSite;
        time_step = 1;
        cp_limit;
        infusion_limit;

        pkpd_model;                 % A discrete state space model (four states)
        pk_model;
        pd_model;

        u;
        target;
        x_internal;                 % internal state of the PK-PD model
    end
    properties (Access = private)
        zero_input_plasma;
        step_plasma;
        zero_input_effect;
        impluse_effect;
    end

    methods

        function tci_setting(obj, drug, mode, cp_limit, infusion_limit, time_step, varargin)
        % Set the TCI settings for the drug.
        % Parameters:
        %   drug: drug [Drug]
        %   mode: target concentration mode [TciMode]
        %   cp_limit: maximum plasma concentration [double]
        %   infusion_limit: maximum infusion rate [double]
        %   time_step: time step [double]
        %   varargin: optional arguments
        %       model_prop_pk_TCI: pharmacokinetic model for propofol [Model]
        %       model_prop_pd_TCI: pharmacodynamic model for propofol [Model]
        %       patient_data: patient data [array]
        %       opiates: opiates [bool]
        %       blood_sampling: blood sampling [string]


            p = inputParser;
            addRequired(p,'drug');
            addRequired(p,'mode');
            addRequired(p,'cp_limit');
            addRequired(p,'infusion_limit');
            addRequired(p,'time_step');
            addOptional(p,'pk_models', []);
            addOptional(p,'pd_models', []);
            addOptional(p,'patient_data', []);
            addOptional(p,'opiates', true);
            addOptional(p,'blood_sampling', 'arterial');

            p.parse(drug, mode, cp_limit, infusion_limit, time_step, varargin{:});
            
            obj.drug = drug;
            obj.mode = mode;
            obj.time_step = time_step;
            obj.cp_limit = cp_limit;
            obj.infusion_limit = infusion_limit;

            patient_data = p.Results.patient_data;

            % Patient demographic parameters
            height = patient_data(1);
            weight = patient_data(2);
            age = patient_data(3);
            gender = patient_data(4);
            bmi = patient_data(5);
            lbm = patient_data(6);

            opiates = p.Results.opiates;
            blood_sampling = p.Results.blood_sampling;

            var_pk = {age, weight, height, gender, bmi, lbm, opiates, blood_sampling};
            var_pd = {age, weight, blood_sampling};

            if p.Results.pd_models('prop') == Model.Wav
                % Parameters to define the Propofol PD model
                e0 = patient_data(7);
                ke0_prop = patient_data(8);
                delay = patient_data(9);
                ec50_prop = 0.01 * patient_data(7) * patient_data(10);
                gammap = patient_data(11);
                pd_data = [e0, ke0_prop, delay, ec50_prop, gammap];
                var_pd = [var_pd, pd_data];
            end

            switch obj.drug
                case Drug.Propofol
                    obj.tci_setting_propofol(p.Results.pk_models('prop'), p.Results.pd_models('prop'), var_pk, var_pd);
                case Drug.Remifentanil
                    obj.tci_setting_remifentanil(p.Results.pk_models('remi'), p.Results.pd_models('remi'), var_pk, var_pd);
                case  Drug.Norepinephrine
                    obj.tci_setting_norepinephrine(p.Results.pk_models('nore'), var_pk);
                    obj.mode = TciMode.Plasma;
                case Drug.Rocuronium
                    obj.tci_setting_rocuronium(p.Results.pk_models('rocu'), var_pk);
            end

            pred_horizon_plasma = 7;
            pred_horizon_effect = 50;

            if obj.mode == TciMode.Plasma
                obj.pkpd_model = obj.pk_model;
                obj.pkpd_model = c2d(obj.pkpd_model, obj.time_step);
                obj.zero_input_plasma = obj.zeroInputResponse(obj.pkpd_model.C(1,:), pred_horizon_plasma);
                obj.step_plasma = obj.stepResponse(obj.pkpd_model.C(1,:), pred_horizon_plasma);
            
            elseif obj.mode == TciMode.EffectSite
                A = [obj.pk_model.A, zeros(3,1); obj.pd_model.B, 0, 0, obj.pd_model.A];
                B = [obj.pk_model.B; 0];
                C = eye(4);
                D = 0;
                obj.pkpd_model = ss(A, B, C, D);
                obj.pkpd_model = c2d(obj.pkpd_model, obj.time_step);
                obj.zero_input_plasma = obj.zeroInputResponse(obj.pkpd_model.C(1,:), pred_horizon_plasma);
                obj.step_plasma = obj.stepResponse(obj.pkpd_model.C(1,:), pred_horizon_plasma);
                obj.zero_input_effect = obj.zeroInputResponse(obj.pkpd_model.C(4,:), pred_horizon_effect);
                obj.impluse_effect = obj.impulseResponse_ce(obj.pkpd_model.C(4,:), pred_horizon_effect);

            else
                error('The target concentration %s is not supported', obj.mode);
            end

            obj.x_internal = zeros(size(obj.pkpd_model.A,2),1);
            
        end

        function tci_setting_propofol(obj, model_prop_pk, model_prop_pd, var_pk, var_pd)
        % Set the TCI settings for propofol.
        % Parameters:
        %   model_prop_pk: pharmacokinetic model for propofol [Model]
        %   model_prop_pd: pharmacodynamic model for propofol [Model]
        %   patient_data: patient data [array]
        %   opiates: opiates [bool]
        %   blood_sampling: blood sampling [string

            pk_object = PharmacokineticModel;
            pd_object = PharmacodynamicDoH;

            pk_object = pk_object.pk_model(Drug.Propofol, model_prop_pk, var_pk{:});
            pd_object = pd_object.pd_model_prop(model_prop_pd, var_pd{:});

            obj.pk_model = pk_object.pk;
            obj.pd_model = pd_object.pd_prop_ce;

        end

        function tci_setting_remifentanil(obj,model_remi_pk, model_remi_pd, var_pk, var_pd)
        % Set the TCI settings for remifentanil.
        % Parameters:
        %   model: pharmacokinetic model for remifentanil [Model]

            pk_object = PharmacokineticModel;
            pk_object = pk_object.pk_model(Drug.Remifentanil, model_remi_pk, var_pk{:});

            pd_object = PharmacodynamicDoH;
            pd_object = pd_object.pd_model_remi(model_remi_pd, var_pd{:});

            obj.pk_model = pk_object.pk;
            obj.pd_model = pd_object.pd_remi_ce;
        end

        function tci_setting_rocuronium(obj, model, var_pk)
        % Set the TCI settings for rocuronium.
            pk_object = PharmacokineticModel;
            pk_object = pk_object.pk_model(Drug.Rocuronium, model, var_pk{:});
            pd_object = PharmacodynamicNMB;
            pd_object = pd_object.pd_model(Drug.Rocuronium);

            obj.pk_model = pk_object.pk;
            obj.pd_model = pd_object.pd_ce;
        end
        function tci_setting_norepinephrine(obj, model, var_pk)
            pk_object = PharmacokineticModel;
            pk_object = pk_object.pk_model(Drug.Norepinephrine, Model.Joachim, var_pk{:});

            obj.pk_model = pk_object.pk;
            obj.mode = TciMode.Plasma;
        end
        
        function reset_state(obj)
            obj.x_internal = zeros(size(obj.pkpd_model.A,2),1);
        end

        function compute_infusion(obj, n_step, target)
        % Compute the infusion rate for the target concentration.
        % Parameters:
        %   n_step: number of steps [double]
        %   target: target concentration [double]

            obj.target = target;
            obj.u = [];
            for i = 1:n_step
                u_interv = tci_interv(obj,obj.zero_input_plasma,obj.step_plasma, ...
                    obj.zero_input_effect,obj.impluse_effect,obj.x_internal);
                if u_interv > obj.infusion_limit
                    u_interv = obj.infusion_limit;
                end
                obj.x_internal= obj.pkpd_model.A*obj.x_internal + obj.pkpd_model.B*u_interv;
                obj.u = [obj.u; u_interv];
            end

        end

        function u_interv = tci_interv(obj,zero_input_plasma,step_input_plasma,zero_input_effect,...
                impulse_effect,x0)
        % Compute the infusion rate for the target concentration.
        % Parameters:
        %   zero_input_plasma: zero input response for plasma [array]
        %   step_input_plasma: step response for plasma [array]
        %   zero_input_effect: zero input response for effect site [array]
        %   impulse_effect: impulse response for effect site [array]
        %   x0: initial state [array]
        % Returns the infusion rate [double]
        
            if obj.mode == TciMode.Plasma || (obj.mode == TciMode.EffectSite && abs(x0(4)-obj.target) < obj.target*0.05) 
                 % If error < 5 percent,  go to plasma control
                    yfree = zero_input_plasma*x0;
                    target_arr = obj.target*ones(length(yfree), 1);
                    u_interv = 1/(step_input_plasma'*step_input_plasma)*step_input_plasma'*(target_arr-yfree); % Least squares solution. (H_step_Plasma'*H_step_Plasma) is a scalar because the control horizon is 1. This is a matrix multiplication, simple algebra.
            elseif obj.mode == TciMode.EffectSite % Effect site control
                % Otherwise go to effect site control
                yfree = zero_input_effect*x0;
                peakfreetmp = find(yfree == max(yfree), 1);
                peakfree = peakfreetmp(1);
                peak = yfree(peakfree);
                if peak > obj.target
                    u_interv = 0;
                else
                    % simultaneously solve for jpeak and infusion (u) as proposed in
                    % schafer and gregg, 1992d
                    tpeak = peakfree; tpeakold =0;
                    infusion = 0; 
                    while tpeak ~= tpeakold
                        tpeakold = tpeak;
                        infusion = (obj.target-yfree(tpeak))/impulse_effect(tpeak);
                        tpeaktmp = find( (yfree+impulse_effect*infusion) == max(yfree+impulse_effect*infusion), 1);
                        tpeak = tpeaktmp(1);
                    end
                    u_interv = infusion(1);
                    x_new = obj.pkpd_model.A*x0 + obj.pkpd_model.B*u_interv;
                    c_p0 = x0(1);
                    c_p_max = x_new(1);
                    if c_p_max > obj.cp_limit
                        u_interv = u_interv *((obj.cp_limit - c_p0)/(c_p_max - c_p0));
                    end
                end
            else
                error('The target concentration %s is not supported', obj.mode);
            end

            u_interv = max(u_interv, 0); % If you're over, set to zero

        end
        function impulse_resp = impulseResponse_ce(obj, c, predictionhorizon)
            % The impluse response to an initial input is calculated and
            % report the evolution of the concentration corresoinding to
            % matrix "c".
            a = obj.pkpd_model.A; 
            b = obj.pkpd_model.B;
            vDtemp = b;
            % The following code assumes single input single output.
            impulse_resp = zeros(predictionhorizon, 1);
            impulse_resp(1,1) = c*b;
            for i = 2:predictionhorizon
                vDtemp =  a*vDtemp;
                impulse_resp(i, 1) = c*vDtemp;
            end
            impulse_resp = impulse_resp(1:end-1); % starting at the response at the next step
        end
        function stepResp = stepResponse(obj, c, predictionhorizon)
            % The step response to a step input is calculated and
            % report the evolution of the concentration corresoinding to
            % matrix "c".
            a = obj.pkpd_model.A; 
            b = obj.pkpd_model.B;
            stepResp = zeros(predictionhorizon, 1);
            x = zeros(size(a, 2), 1);
            for i=1:predictionhorizon
                x = a*x + b*1;
                stepResp(i, 1) = c*x;
            end
            stepResp = stepResp(1:end-1); % starting at the response at the next step
        end
        function zero_resp  = zeroInputResponse(obj, c, predictionhorizon)
            % The response to zero input and an initial non-zero
            % concentration corresponding to matrix "c" is calculated and
            % report the evolution of the concentrations.
            a = obj.pkpd_model.A;
            zero_resp = zeros(predictionhorizon, size(a, 2));
            zero_resp(1,:) = c;
            for i = 2:predictionhorizon
                zero_resp(i, :) = zero_resp(i-1,:)*a;
            end
            zero_resp  = zero_resp(2:end, :); % starting with the response at the next step.
        end


    end
end