classdef PharmacodynamicNMB
     %{
    NMB: Nueromuscular Blockade
    m = 0 (NMB is moderate)
    m = 1 (NMB is deep)
    m = 2 (NMB is profound)
    m = 3 (NMB is very profound)

    Please refer to the following paper for this PD model:

    Comparison of two pharmacokinetic–pharmacodynamic models of rocuronium bromide
    during profound neuromuscular block: analysis of estimated and measured 
    post-tetanic count effect, M. Couto et al., BJ Anaesthesia (2022)
    %}
    properties (Access = public)
        drug
        pd_ce       % Linear part of PD model to calculate the rocuronium effect-site concentration
    end
    
    properties (Access = private)
        % Drug-specific parameters 
        gam;        % Steepness of the effect-concentration curve 
        ke0;        % Flow rates from effect-site compartment to the central one for rocuronium [1/min]
        
        ec50_m1;    % Rocuronium effect-site concentration associated with 50% of probability for the m=1 category [µg/mL]
        ec50_m2;    % Rocuronium effect-site concentration associated with 50% of probability for the m=2 category [µg/mL]
        ec50_m3;    % Rocuronium effect-site concentration associated with 50% of probability for the m=3 category [µg/mL]
    end

    methods

        function obj = pd_model(obj, drug)
        % Set the pharmacodynamic model for the drug.
        % Parameters:
        %   drug: drug [Drug]

            obj.drug = drug;
            switch obj.drug
                case Drug.Rocuronium
                    obj.ke0 = 0.134;
                    obj.ec50_m1 = 1.10;
                    obj.ec50_m2 = 1.5;
                    obj.ec50_m3 = 2.2;
                    obj.gam = 4.5;

                otherwise
                    error('Drug not supported: %s', obj.drug);
            end
            
            obj.pd_ce = obj.getSS();
        end

        function e = hillfun(obj, x)
        % Hill function for the pharmacodynamic model.
        % Parameters:
        %   x: input concentration [float]
        % Returns the effect of the drug.

            x(x<0) = 0;
            p1 = x.^obj.gam./(obj.ec50_m1.^obj.gam + x.^obj.gam);
            p2 = x.^obj.gam./(obj.ec50_m2.^obj.gam + x.^obj.gam);
            p3 = x.^obj.gam./(obj.ec50_m3.^obj.gam + x.^obj.gam);
            
            p_m0 = 1 - p1;   % Probability for the m=0 category
            p_m1 = p1 - p2;  % Probability for the m=1 category
            p_m2 = p2 - p3;  % Probability for the m=2 category
            p_m3 = p3;       % Probability for the m=3 category
            e = [p_m0, p_m1, p_m2, p_m3];
        end

        function pd_ce = getSS(obj)
        % Get the State-Space model for the pharmacodynamic model.
        % Returns the State-Space model.
            a = -obj.ke0/60;
            b = obj.ke0/60;
            c = 1;
            d = 0;
            pd_ce = ss(a, b, c, d);
        end

    end

end
