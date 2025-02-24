classdef PharmacodynamicHemodynamics
    %{
    TPR: total peripheral resistance
    SV: stroke volume
    HR: heart rate
    MAP: mean arterial pressure
    CO: cardiac output
    TDE: time-dependent effect

    Please refer to the following two papers for this PD model:
    1) Pharmacodynamic mechanism-based interaction model for the haemodynamic effects of remifentanil and propofol
       in healthy volunteers. Su et al., Br J Anaesthesia, 2023.
    2) Design of a pharmacokinetic/pharmacodynamic model for administration of low dose peripheral norepinephrine during
       general anaesthesia. Joachim et a., Br J Clin Pharmacol, 2024
    %}
    properties (Access = private)
        emax_tpr_prop = -0.778; % maximum effect of propofol on TPR
        ec50_tpr_prop = 3.21;   % c_e that produce half of the maximal propofol effect of TPR [µg ml^-1]
        emax_sv_prop = -0.154;  % maximum effect of propofol on SV
        ec50_sv_prop = 0.44;    % c_e that produce half of the maximal propofol effect of SV [µg ml^-1]
        emax_tpr_remi = -1;     % maximum effect of remifentanil on TPR
        ec50_tpr_remi = 4.59;   % c_e that produce half of the maximal remifentanil effect of TPR [ng ml^-1]
        sl_sv_remi = 0.0581;    % slope of remifentanil effect on SV [ng ml^-1]
        sl_hr_remi = 0.0327;    % slope of remifentanil effect on HR [ng ml^-1]
        
        % maximum magnitude changes of slope of remifentanil on TPR, SV, and HR caused by propofol
        int_tpr = 1;
        int_hr = -0.119;        % [ng ml^-1]
        int_sv = -0.212;        % [ng ml^-1] 
        ec50_int_hr = 0.196;    % Interaction potency of propofol on the effect of remifentanil on HR [µg ml^-1]
        
        % Hill coefficient for TPR sigmoid
        gamma_tpr_prop = 1.83; 
        gamma_tpr_remi = 1;

        fb = 0.661;             % magnitude of the feedback
        kout = 0.072/60;        % first-order dissipation rate constant for SV, HR, and TPR [second^-1]
        k_tde = 0.067/60;       % first-order dissipation rate constant for TDE [second^-1]

        % Inter-patient variability
        omega_base_sv = sqrt(0.0328);         
        omega_base_tpr = sqrt(0.0528);       
        omega_base_hr = sqrt(0.0242);    
        omega_ec50_tpr_prop = sqrt(0.44); 
        omega_emax_tpr_remi = sqrt(0.449); 
        omega_sl_hr_remi = sqrt(0.00382);      
        omega_sl_sv_remi = sqrt(0.00868);
        emax_sv_age = 0.0333;

        % Hill function parameters for norepinephrine effect on hemodynamics variables
        ec50_nore_map = 15.27;  %[nmol l^-1]
        ec50_nore_co = 36;      %[nmol l^-1]
        gamma_nore_map = 1.46;       
        gamma_nore_co = 2.3;
        t_lag = 25              % [second]
    end
    
    properties (Access = public)
        hr_sv = 0.312;          % magnitude of the inverse effect of HR on SV
        ltde_hr = 0.121;        % percentage of increased baseline HR caused by the time-dependent effect
        ltde_sv = 0.0899;       % percentage of increased baseline HR caused by the time-dependent effect
        base_sv = 82.2;         % base value of SV [ml]
        base_hr = 56.1;         % base value of HR [beats min^-1]
        base_tpr = 0.0163;      % base value of TPR [mmHg ml^-1 min]   
    end

    methods

        function obj = PharmacodynamicHemodynamics(age)
            obj.emax_sv_prop = obj.emax_sv_prop*exp(obj.emax_sv_age*(age-35));
        end

        function obj = interpatient_variability(obj, seed)
            % create model's parameters with interpatient_variability
            if ~isempty(seed)
                rng(seed,'twister');
            end
            var_rand = rand(1,1) * 0.5;  
            obj.base_sv = obj.base_sv*exp(obj.omega_base_sv*var_rand);         
            obj.base_hr = obj.base_hr*exp(obj.omega_base_hr*var_rand);        
            obj.base_tpr = obj.base_tpr*exp(obj.omega_base_tpr*var_rand);
            obj.ec50_tpr_prop = obj.ec50_tpr_prop*exp(obj.omega_ec50_tpr_prop*var_rand);
            obj.emax_tpr_remi = obj.emax_tpr_remi + obj.omega_emax_tpr_remi*var_rand;
            obj.sl_hr_remi = obj.sl_hr_remi + obj.omega_sl_hr_remi*var_rand;
            obj.sl_sv_remi = obj.sl_sv_remi + obj.omega_sl_sv_remi*var_rand;
        end

        function e = hillfun(obj, x, output)
        % Hill function for the pharmacodynamic model.
        % Parameters:
        %   x: input concentration [array]
        % Returns the effect of the drug.
            switch output
                case 'map'
                    gamma = obj.gamma_nore_map;
                    ec50 = obj.ec50_nore_map;                                                    %[mmHg]
                    emax = obj.base_hr*(1+obj.ltde_hr)*obj.base_sv*(1+obj.ltde_sv)*obj.base_tpr; %[mmHg]
                case 'co'
                    gamma = obj.gamma_nore_co;
                    ec50 = obj.ec50_nore_co;                                                     %[L/min]
                    emax = obj.base_hr*(1+obj.ltde_hr)*obj.base_sv*(1+obj.ltde_sv)/1000;         %[L/min]
                otherwise
                    error('Norepinephrine is not considered to affect: %s', output);
            end
            x(x<0) = 0;
            % Hill function
            e = x.^gamma./(ec50.^gamma + x.^gamma);
            e = emax*e;
        end

        function [tpr_prop, tpr_remi, sv_prop, sv_remi, hr_remi] = prop_remi_inter(obj, cp_prop, cp_remi)
        % Interaction between propofol and remifentanil.
        % Parameters:
        %   cp_prop: propofol concentration [array]
        %   cp_remi: remifentanil concentration [array]
        % Returns the effect of the drug.
            tpr_prop = (obj.emax_tpr_prop + obj.int_tpr*cp_remi./(obj.ec50_tpr_remi + cp_remi)).*...
                cp_prop.^obj.gamma_tpr_prop./(obj.ec50_tpr_prop^obj.gamma_tpr_prop + cp_prop.^obj.gamma_tpr_prop);
            tpr_prop(tpr_prop < -0.999) = -0.999;

            sv_prop = (obj.emax_sv_prop*cp_prop)./(obj.ec50_sv_prop + cp_prop);
            sv_prop(sv_prop<-0.999) = -0.999;

            tpr_remi = obj.emax_tpr_remi*cp_remi.^obj.gamma_tpr_remi./...
                (obj.ec50_tpr_remi^obj.gamma_tpr_remi + cp_remi.^obj.gamma_tpr_remi);
            tpr_remi(tpr_remi>0.999) = 0.999;

            sv_remi = ( obj.sl_sv_remi + obj.int_sv*cp_prop./(obj.ec50_sv_prop + cp_prop)).*cp_remi;
            sv_remi(sv_remi>0.999) = 0.999;

            hr_remi = (obj.sl_hr_remi + obj.int_hr*cp_prop./(obj.ec50_int_hr + cp_prop)).*cp_remi;
            hr_remi(hr_remi>0.999) = 0.999;
        end

        function dydt = ode_prop_hemodynamic(obj, ~ , y, cp_prop, cp_remi)
        % ODE for the pharmacodynamic model.
        % Parameters:
        %   t: time [double]
        %   y: state variables [array]
        %   cp_prop: propofol concentration [array]
        %   cp_remi: remifentanil concentration [array]
        % Returns the derivative of the state variables.

            dydt = zeros(5,1);
            
            [tpr_prop, tpr_remi, sv_prop, sv_remi, hr_remi] = obj.prop_remi_inter(cp_prop, cp_remi);

            % tpr_prop = u(1); tpr_remi = u(2); sv_prop = u(3); sv_remi = u(4); hr_remi = u(5);

            k_in_hr = obj.kout*(obj.base_hr);
            k_in_sv = obj.kout*(obj.base_sv);
            k_in_tpr = obj.kout*(obj.base_tpr);
          
            tde_hr = y(4) + obj.base_hr*obj.ltde_hr;
            tde_sv = y(5) + obj.base_sv*obj.ltde_sv;

            hr = y(3) + tde_hr;
            sv = (y(2) + tde_sv)*(1 - obj.hr_sv*log(hr/(obj.base_hr*(1+obj.ltde_hr))));

            rmap = y(1)*hr*sv/(obj.base_hr*(1+obj.ltde_hr)*obj.base_sv*(1+obj.ltde_sv)*obj.base_tpr);

            dydt(4) = -obj.k_tde*y(4);
            dydt(5) = -obj.k_tde*y(5);
            dydt(1) = k_in_tpr*rmap^(-obj.fb)*(1 + tpr_prop) - obj.kout*y(1)*(1-tpr_remi);
            dydt(2) = k_in_sv*rmap^(-obj.fb)*(1 + sv_prop) - obj.kout*y(2)*(1 - sv_remi);
            dydt(3) = k_in_hr*rmap^(-obj.fb) - obj.kout*y(3)*(1 - hr_remi);
        end

        function nore_delay = delay_ss(obj)     % delay corresponding to depot compartment
            s = tf('s');
            nore_delay = ss(pade(exp(-obj.t_lag*s),1));
        end
    end
    
end