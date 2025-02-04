classdef PharmacodynamicDoH
    properties (Access = public)
        propofol_model 
        remifentanil_model
        delay = 0;          % default
        e0 = 100;           % default

        pd_prop_ce
        pd_remi_ce
    end

    properties (Access = private)
        % Drug-specific parameters 
        ec50_prop
        ec50_remi
        gam_high_ce;        % default = gam
        gam
        ke0_prop
        ke0_remi
    end

    methods (Static)
        function gns_sim = neurow_sensor_dynamics()
        % Create the neuromuscular sensor dynamics.
        % Returns the continuous-time transfer function of the sensor dynamics.
            %% Sensor dynamics
            % Af =[ -1.76499380516919   0.77880078307141]; %30 sec IIR filter
            % Bf =[  0.00719098459233   0.00661599330988];
            % Gns = tf(Bf, [1 Af], 1, 'variable', 'z^-1');

            % gns_sim = d2c(ss(Gns));
            A = [0.875000000000001	-0.882496902584599;
            1.13314845306682	-1.12499999999999];
            B = [0.0678724069605523;
            -0.0737711496728008];
            C = [0.154464292547276	-0.0448027554524884];
            D = 0.00719098459233000;
            gns_sim= ss(A,B,C,D);
        end

        function [bis_delay_ss, bis_lti] = bis_sensor_dynamics()
        % Create the BIS sensor dynamics.
        % Returns the continuous-time state-space model of the BIS delay and the LTI model of the BIS sensor dynamics.
            s = tf('s');
            bis_delay = 5;      % 5 seconds delay
            bis_delay_ss = ss(pade(exp(-bis_delay*s),1));
            num = [-0.049 0.066 0.976 0.056];
            den = [1 2.27 6.353 1.487 0.056];
            [A, B, C, D]= tf2ss(num, den);
            bis_lti = ss(A, B, C, D);
        end

        function pd_ce = getSS(ke0)
        % Create the pharmacodynamic model for propofol.
        % Returns the continuous-time state-space model of the pharmacodynamic model.
            a = -ke0/60;
            b = ke0/60;
            c = 1;
            d = 0;
            pd_ce = ss(a, b, c, d);
        end

    end

    methods

        function obj = pd_model_prop(obj, propofol_model, varargin)
        % Set the pharmacodynamic model for propofol.
        % Parameters:
        %   propofol_model: propofol model [Model]
        %   varargin: additional parameters for the propofol model

            obj.propofol_model = propofol_model;
            switch obj.propofol_model
                case Model.Wav
                    data = varargin{4};
                    obj.e0 = data(1);
                    obj.ke0_prop = data(2);
                    obj.delay = data(3);
                    obj.ec50_prop = data(4);
                    obj.gam = data(5);
                    obj.gam_high_ce = obj.gam;
     
                case Model.Schnider
                    age = varargin{1};
                    obj.ke0_prop = 0.456; %[min^(-1)]
                    obj.ec50_prop = 2.9 - 0.022*age;
                    obj.delay = 0;
                    obj.gam = 1.43;
                    obj.gam_high_ce = obj.gam;
                    
                case Model.Eleveld
                    age = varargin{1};
                    wgt = varargin{2};
                    blood_sampling = varargin{3};

                    thetaPD_1 = 3.08;                     % Ce50 [microgram/ml]
                    thetaPD_2 = 0.146;                    % ke0 for arterial samples [1/min]
                    thetaPD_3 = 93;                       % Baseline BIS value
                    thetaPD_4 = 1.47;                     % PD sigmoid slope (Ce > Ce50)
                    % thetaPD_5 = 8.03;                   % Residual error [BIS]
                    % thetaPD_6 = 0.0517;                 % Increase in delay with age
                    thetaPD_7 = -0.00635;                 % Decrease in Ce50 with age
                    thetaPD_8 = 1.24;                     % ke0 for veous samples [1/min]
                    thetaPD_9 = 1.89;                     % PD sigmoid slope (Ce < Ce50)

                    % The effect of the blood sampling site on the volume and clearance rates
                    if strcmp(blood_sampling, 'arterial')
                        theta_ke0 = thetaPD_2;
                    elseif strcmp(blood_sampling, 'venous')
                        theta_ke0 = thetaPD_8;
                    end

                    f_agingTheta7 = exp(thetaPD_7*(age - 35));
                    obj.ec50_prop = thetaPD_1*f_agingTheta7;                          % [microgram/ml]
                    obj.ke0_prop = theta_ke0*(wgt/70)^(-0.25);                        % [1/min]
                    obj.e0 = thetaPD_3;
                    obj.delay = 0;
                    obj.gam = thetaPD_9;
                    obj.gam_high_ce = thetaPD_4;                    
                otherwise
                    error('Propofol PD model is not supported: %s', obj.propofol_model);
            end
            
            obj.pd_prop_ce = obj.getSS(obj.ke0_prop);

        end

        function obj = pd_model_remi(obj, remifentanil_model, varargin)
        % Set the pharmacodynamic model for remifentanil.
        % Parameters:
        %   remifentanil_model: remifentanil model [Model]
        %   varargin: additional parameters for the remifentanil model

            age = varargin{1};
            if remifentanil_model == Model.Minto
                obj.ke0_remi = 0.595 - 0.007*(age - 40);
                obj.ec50_remi = 13.1 - 0.148*(age - 40);
            elseif remifentanil_model == Model.Eleveld
                theta1 = -0.0289;
                obj.ec50_remi = 12.7;                        % [ng/ml]
                obj.ke0_remi = 1.09 * exp(theta1 * (age - 35));     % [1/min]
            end

            obj.pd_remi_ce = obj.getSS(obj.ke0_remi);

        end

        function e = hillfun(obj, x)
        % Hill function for the pharmacodynamic model.
        % Parameters:
        %   x: input concentration [array]
        % Returns the effect of the drug.

            x(x<0) = 0;
            e = zeros(length(x), 1);
            for i = 1:length(x)
                if x(i) < obj.ec50_prop
                    gamma = obj.gam;
                else
                    gamma = obj.gam_high_ce;
                end
                e(i) = x(i)^gamma/(obj.ec50_prop^gamma + x(i)^gamma);
            end
            e = obj.e0 - obj.e0.*e;
        end

        function e = responseSurfaceModel(obj, x_prop, x_remi)
        % Response surface model for the pharmacodynamic model.
        % Parameters:
        %   x_prop: propofol concentration [array]
        %   x_remi: remifentanil concentration [array]
        % Returns the effect of the drug.

            beta = 1;
            inter_prop = x_prop./obj.ec50_prop; 
            inter_remi = x_remi./obj.ec50_remi;
            % eps is a small number added to avoid zero in denominator.
            theta = inter_prop./(inter_prop + inter_remi + eps);
            inter = (inter_prop + inter_remi)./(1 - beta*theta+ beta*theta.^2);
            inter(inter<0) = 0;
            e = obj.e0 - obj.e0.*(inter.^obj.gam./(1 + inter.^obj.gam));

        end
       
        function pd_ce_delayed = delay_ss(obj)
        % Create the delayed pharmacodynamic model.
        % Returns the continuous-time state-space model of the delayed pharmacodynamic model.
            s = tf('s');
            pd_ce_delayed = ss(pade(exp(-obj.delay*s),1));  
        end
    end

end