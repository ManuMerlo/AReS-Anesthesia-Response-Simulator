classdef PharmacokineticModel
    properties (Access = private)
        
        % Compartment models parameters
        v1comp3             % Volume of compartment 1  [L]
        v2comp3             % Volume of compartment 2  [L]
        v3comp3             % Volume of compartment 3  [L]
        cl1comp3            % Clearance rate of compartment 1 [L/min]
        cl2comp3            % Clearance rate of compartment 2 [L/min]
        cl3comp3            % Clearance rate of compartment 3 [L/min]
        k10                 % Drug flow rate from compartment 1 to outside of body [1/min]
        k12                 % Drug flow rate from compartment 1 to 2 [1/min]
        k13                 % Drug flow rate from compartment 1 to 3 [1/min]
        k21                 % Drug flow rate from compartment 2 to 1 [1/min]
        k31                 % Drug flow rate from compartment 3 to 1 [1/min]
        ka                  % Drug flow rate from depot compartment to central one [1/min]
    end

    properties (Access = public)
        drug
        model
        pk
    end
    methods(Static)
        % The helper functions defined by Eleveld et al.(2018) for propfol
        % pk model to account for demographic parameters
        function par = f_aging(x, age)
            par = exp(x .* (age - 35));
        end

        function par = f_sigmoid(x, e50, gamma)
            par = (x .^ gamma) ./ (x .^ gamma + e50 .^ gamma);
        end
    end

    methods
        function obj = pk_model(obj, drug, model, varargin)
        % Set the pharmacokinetic model for the drug.
        % Parameters:
        %   drug: drug [Drug]
        %   varargin: additional parameters for the drug model

            obj.drug = drug;
            obj.model = model;

            switch obj.drug
                case Drug.Propofol

                    age = varargin{1};
                    weight = varargin{2};

                    if obj.model == Model.Schnider
                        height = varargin{3};
                        lbm = varargin{6};

                        % Propofol parameters (Schnider model)
                        obj.v1comp3 = 4.27;
                        obj.v2comp3 = 18.9 - 0.391 * (age - 53);
                        obj.v3comp3 = 238;  % [l]
                        obj.cl1comp3 = 1.89 + 0.0456 * (weight - 77) - 0.0681 * (lbm - 59) + 0.0264 * (height - 177);
                        obj.cl2comp3 = 1.29 - 0.024 * (age - 53);
                        obj.cl3comp3 = 0.836;  % [l/min]
                    
                    elseif obj.model == Model.Eleveld
                        gender = varargin{4};
                        bmi = varargin{5};
                        opiates = varargin{7};
                        blood_sampling = varargin{8};

                        % Propofol parameters (Eleveld model)

                        theta_1 = 6.28;                     % V1_ref [litre]
                        theta_2 = 25.5;                     % V2_ref [litre]
                        theta_3 = 273;                      % V3_ref [litre]
                        theta_4 = 1.79;                     % CL_ref (male) [Litre/min]
                        theta_5 = 1.75;                     % Q2_ref [Litre/min]
                        theta_6 = 1.11;                     % Q3_ref [Litre/min]
                        % theta_7 = 0.191;                    % Typical residual error
                        theta_8 = 42.3;                     % CL maturation E50 [weeks]
                        theta_9 = 9.06;                     % CL maturation E50 slope
                        theta_10 = -0.0156;                 % Smaller V2 with age
                        theta_11 = -0.00286;                % Lower CL with age
                        theta_12 = 33.6;                    % Weight for 50% of maximal V1 [kg]
                        theta_13 = -0.0138;                 % Smaller V3 with age
                        theta_14 = 68.3;                    % Maturation of Q3 [weeks]
                        theta_15 = 2.10;                    % CL_ref (female) [Litre/min]
                        theta_16 = 1.3;                     % Higher Q2 for maturation of Q3
                        theta_17 = 1.42;                    % V1 venous samples (children)
                        theta_18 = 0.68;                    % Higher Q2 venous samples

        
                        % Intermediate Coefficients
                        pma = age*365.25/7 + 40;          % Post Menstural Age in weeks
                        f_agingTheta10 = obj.f_aging(theta_10, age);
                        f_centralWGT = weight/(weight + theta_12);
                        f_centralWGTRef = 70/(70 + theta_12);
                        f_CLmaturation = pma^(theta_9)/(pma^(theta_9) + theta_8^(theta_9));
                        f_CLmaturation_ref = 1866^(theta_9)/(1866^(theta_9) + theta_8^(theta_9));
                        f_Q3maturation = (age*365.25/7 + 40)/(age*365.25/7 + 40 + theta_14);
                        f_Q3maturation_ref = (35*365.25/7 + 40)/(35*365.25/7 + 40 + theta_14);

                        % The impact of gender on the intermediate coefficients
                        if gender   % female = 1
                            f_AlSallami = (1.11 + (1-1.11)/(1 + (age/7.1)^(-1.1)))...
                                *(9270*weight/(8780 + 244*bmi));
                            theta_CL = theta_15;
                        else        % male = 0
                            f_AlSallami = (0.88 + (1-0.88)/(1 + (age/13.4)^(-12.7)))...
                                *(9270*weight/(6680 + 216*bmi));
                            theta_CL = theta_4;
                        end
                        f_AlSallami_ref = (0.88 + (1-0.88)/(1 + (35/13.4)^(-12.7)))...
                            *(9270*70/(6680 + 216*24.2));

                        % The impact of opiates' presence on the intermediate coefficients,
                        if opiates
                            f_opiates_theta13 = exp(theta_13*age);
                            f_opiates_theta11 = exp(theta_11*age);
                        else
                            f_opiates_theta13 = 1;
                            f_opiates_theta11 = 1;
                        end

                        % Volume of compartments [Litre]
                        V1_arterial = theta_1*f_centralWGT/f_centralWGTRef;
                        V1_venous = V1_arterial*(1 + theta_17*(1-f_centralWGT));
                        obj.v2comp3 = theta_2*(weight/70)*f_agingTheta10;
                        obj.v3comp3 = theta_3*(f_AlSallami/f_AlSallami_ref)*f_opiates_theta13;

                        % Elimination and inter-compartment clearance rates [litre/min]
                        obj.cl1comp3 = theta_CL*(weight/70)^0.75*f_CLmaturation*f_opiates_theta11/f_CLmaturation_ref;
                        Q2_arterial = theta_5*(obj.v2comp3/theta_2)^0.75*(1+theta_16*(1- f_Q3maturation));
                        Q2_venous = Q2_arterial*theta_18;
                        obj.cl3comp3 = theta_6*(obj.v3comp3/theta_3)^0.75*f_Q3maturation/f_Q3maturation_ref;

                        % The effect of the blood sampling site on the volume and clearance rates
                        if blood_sampling == BloodSampling.ARTERIAL
                            obj.v1comp3 = V1_arterial;
                            obj.cl2comp3 = Q2_arterial;
                        elseif blood_sampling == BloodSampling.VENOUS
                            obj.v1comp3 = V1_venous;
                            obj.cl2comp3 = Q2_venous;
                        end
                    end 

                case Drug.Remifentanil
                    age = varargin{1};

                    if obj.model == Model.Minto
                        lbm = varargin{6};

                        % Remifentanil parameters (Minto model)

                        obj.v1comp3 = 5.1 - 0.0201 * (age - 40) + 0.072 * (lbm - 55);
                        obj.v2comp3 = 9.82 - 0.0811 * (age - 40) + 0.108 * (lbm - 55);
                        obj.v3comp3 = 5.42;
                        obj.cl1comp3 = 2.6 - 0.0162 * (age - 40) + 0.0191 * (lbm - 55);
                        obj.cl2comp3 = 2.05 - 0.0301 * (age - 40);
                        obj.cl3comp3 = 0.076 - 0.00113 * (age - 40);

                    elseif obj.model == Model.Eleveld
                        weight = varargin{2};
                        gender = varargin{4};
                        bmi = varargin{5};

                        % Remifentanil parameters (Eleveld model)

                        v1_ref = 5.81;   % [L]
                        v2_ref = 8.82;   % [L]
                        v3_ref = 5.03;   % [L]
                        cl_ref = 2.58;   % [L/min]
                        q2_ref = 1.72;   % [L/min]
                        q3_ref = 0.124;  % [L/min]
                        theta1 = 2.88;
                        theta2 = -0.00554;
                        theta3 = -0.00327;
                        theta4 = -0.0315;
                        theta5 = 0.47;
                        theta6 = -0.0260;

                        if gender  % female = 1
                            ffm = (1.11 + (1 - 1.11) / (1 + (age / 7.1) ^ (-1.1))) * (9270 * weight / (8780 + 244 * bmi));
                        else
                            ffm = (0.88 + (1 - 0.88) / (1 + (age / 13.4) ^ (-12.7))) * (9270 * weight / (6680 + 216 * bmi));
                        end

                        ffm_ref = (1.11 + (1 - 1.11) / (1 + (35 / 7.1) ^ (-1.1))) * (9270 * 70 / (8780 + 244 * 24.2));
                        size_wg = ffm/ffm_ref;

                        kmat = obj.f_sigmoid(weight, theta1, 2);
                        kmat_ref = obj.f_sigmoid(70, theta1, 2);
                        if gender
                            ksex = 1 + theta5 * obj.f_sigmoid(age, 12, 6) * (1 - obj.f_sigmoid(age, 45, 6));
                        else
                            ksex = 1;
                        end

                        obj.v1comp3 = v1_ref * size_wg * obj.f_aging(theta2,age);
                        obj.v2comp3 = v2_ref * size_wg * obj.f_aging(theta3,age) * ksex;
                        obj.v3comp3 = v3_ref * size_wg * obj.f_aging(theta4,age) * exp(theta6 * weight - 70);
                        obj.cl1comp3 = cl_ref * size_wg ^ 0.75 * (kmat / kmat_ref);
                        obj.cl2comp3 = q2_ref * (obj.v2comp3 / v2_ref) ^ 0.75 * obj.f_aging(theta2,age) * ksex;
                        obj.cl3comp3 = q3_ref * (obj.v3comp3 / v3_ref) ^ 0.75 * obj.f_aging(theta2,age);

                    else
                        error('Model not supported: %s', obj.model);
                    end

                case Drug.Norepinephrine

                    % Norepinephrine parameters (Joachim model)

                    if obj.model == Model.Joachim
                        lbm = varargin{6};
                        obj.k10 = exp(0.64*log(lbm) - 5.52);    % [1/s]
                        obj.ka = 0.02;                          % [1/s]
                        obj.k12 = 0.06;                         % [1/s]
                        obj.k21 = 0.04;                         % [1/s]
                        obj.v1comp3 = 0.49;                     % [L]
                    else 
                        error('Model not supported: %s', obj.model);
                    end
                        
                case Drug.Rocuronium
                    if obj.model == Model.Dahe

                        % Rocuronium parameters (Da Haes Model)

                        weight = varargin{2};
                        obj.v1comp3 = 42*weight/1000;     % [L]
                        obj.v2comp3 = 40*weight/1000;
                        obj.v3comp3 = 69*weight/1000;
                        obj.cl1comp3 = 3.2*weight/1000;     %[l/min]
                        obj.cl2comp3 = 5.2*weight/1000;       
                        obj.cl3comp3 = 0.9*weight/1000;  
                    else
                        error('Model not supported: %s', obj.model);
                    end    

                otherwise
                    error('Drug not supported: %s', obj.drug);
            end

            obj.pk = obj.getSS();
        end

        function pk = getSS(obj)
        % Get the State-Space model for the pharmacokinetic model.
        % Returns the State-Space model.

            % Initialize the state-space matrices based on the drug
            switch obj.drug
                case {Drug.Propofol, Drug.Remifentanil, Drug.Rocuronium} % Three compartment model
                    obj.k10 = obj.cl1comp3 / obj.v1comp3;
                    obj.k12 = obj.cl2comp3 / obj.v1comp3;
                    obj.k13 = obj.cl3comp3 / obj.v1comp3;
                    obj.k21 = obj.cl2comp3 / obj.v2comp3;
                    obj.k31 = obj.cl3comp3 / obj.v3comp3;

                    A = [-(obj.k10 + obj.k12 + obj.k13) obj.k21 obj.k31;...
                        obj.k12 -obj.k21 0; obj.k13 0 -obj.k31]/60;      % simulation is in seconds
                    B = [1/obj.v1comp3;0;0];
                    C = [1 0 0];
                    D = 0;
                    pk = ss(A, B, C, D);


                case Drug.Norepinephrine                                 % Two compartment plus depot PK model
                    A = [-(obj.k10 + obj.k12) obj.k21 obj.ka;...
                         obj.k12 -obj.k21 0; 0 0 -obj.ka];
                    B = [0; 0; 1/obj.v1comp3];
                    C = [1 0 0];
                    D = 0;
                    pk = ss(A, B, C, D);

                otherwise
                    error('Drug not supported: %s', obj.drug);

            end
        end
    end
end
