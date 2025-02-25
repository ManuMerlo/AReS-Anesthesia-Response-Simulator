classdef Factory
    % Factory class to create the simulator using different input setting
    % modes
    methods (Static)
        function obj = createSimulator(mode)
            if mode == SimulatorMode.Infusion
                obj = Simulator();
            elseif mode == SimulatorMode.Concentration
                obj =  SimulatorConcentration();
            end
        end
    end
end
