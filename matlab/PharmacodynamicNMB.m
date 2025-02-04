classdef PharmacodynamicNMB
    properties (Access = public)
        drug
        pd_ce
    end
    
    properties (Access = private)
        % Drug-specific parameters 
        gam;
        ke0;        %  [min ^ (-1)]
        
        % Coute model
        ec50_m1;    % [µg ml ^ (-1)]
        ec50_m2;    % [µg ml ^ (-1)]
        ec50_m3;    % [µg ml ^ (-1)]
    end

    methods

        function obj = pd_model(obj, drug)
        % Set the pharmacodynamic model for the drug.
        % Parameters:
        %   drug: drug [Drug]

            obj.drug = drug;
            switch obj.drug
                case Drug.Rocuronium
                    % Rocuronium parameters (Coute Model)
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
            
            p_m0 = 1 - p1;   % Probability between [0 1] if m=0
            p_m1 = p1 - p2;  % Probability between [0 1] if m=1
            p_m2 = p2 - p3;  % Probability between [0 1] if m=2
            p_m3 = p3;       % Probability between [0 1] if m=3
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
