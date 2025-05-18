% matlab version of the main function
function main()

    function complete_simulation(stimuli,volume_status)
        % Complete simulation of a patient with a complete sequence of propofol and remifentanil infusion rates.

        simulator = Factory.createSimulator(SimulatorMode.Infusion);

        for k = 1:44
            disp('Patient id:' + string(k));
            interaction = Interaction.Surface;
            dohMeasure = DoHMeasure.Both;
            vs=volume_status;
            t_sim = 30 * 60;
            t_s = 5;

            pk_models = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Eleveld});
            pd_models = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Eleveld});

            var = {'stimuli',stimuli,'volume_status', vs,'opiates', true,...
                'blood_sampling', BloodSampling.ARTERIAL,'interaction', interaction,...
                'doh_measure', dohMeasure,'pk_models', pk_models, ...
                'pd_models', pd_models,'seed_disturbance',42};

            simulator.init_simulation_from_file(k,t_sim,t_s,var{:});

            u_prop = ones(t_sim,1) * 0.1;
            u_remi = ones(t_sim,1) * 0.1;
            u_nore = ones(t_sim,1) * 0.01;
            u_rocu = ones(t_sim,1) * 0.1;

            simulator.run_complete_simulation(u_prop, u_remi, u_nore, u_rocu);
            simulator.save_simulation()

        end

        results = simulator.get_patient_results();

        % initial letter of the pk_mode:
        if interaction == Interaction.Surface
            inter = 'Yes';
        elseif interaction == Interaction.No_interaction
            inter = 'No';
        else
            error('Invalid interaction type');
        end

        % first letter of the pk_mode:
        model_prop_pk_str = char(pk_models('prop'));
        pd_model_remi_pk_str = char(pk_models('remi'));
        model_prop_pd_str = char(pd_models('prop'));
        pd_model_remi_pd_str = char(pd_models('remi'));
        
        path = './csv/infusion/';
        filename = strcat('pk_', model_prop_pk_str(1),'_', pd_model_remi_pk_str(1), '_pd_', model_prop_pd_str(1), '_', pd_model_remi_pd_str(1), '_', inter, '.csv');
        simulator.save_to_csv(path,filename);
        simulator.plot_simulation();
    end

    function step_by_step_simulation(stimuli,volume_status)
        % Step by step simulation of a patient with a single value of propofol and remifentanil infusion rates at each step.

        simulator = Factory.createSimulator(SimulatorMode.Infusion);

        patient_id = 1;
        disp('Patient id:' + string(patient_id));
        
        pk_models = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Minto});
        pd_models = containers.Map({'prop', 'remi'}, {Model.PATIENT_SPECIFIC, Model.Minto});
        interaction = Interaction.Surface;
        dohMeasure = DoHMeasure.Both;
        vs=volume_status;

        t_sim = 30 * 60;
        t_s = 5;
        var = {'opiates', true, 'blood_sampling', BloodSampling.ARTERIAL,...
            'interaction', interaction, 'doh_measure', dohMeasure,...
            'stimuli', stimuli, 'volume_status', vs, 'pk_models', ...
            pk_models, 'pd_models', pd_models};
        simulator.init_simulation_from_file(patient_id,t_sim,t_s,var{:});
        
        u_prop = ones(t_sim,1) * 0.2;
        u_remi = ones(t_sim,1) * 0.1;
        u_nore = ones(t_sim,1) * 0.0;
        u_nore(1:10*60,1) = 0.0;
        u_rocu = ones(t_sim,1) * 0.02;
        
        steps = t_sim / t_s;
        
        for time = 1:steps
            simulator.one_step_simulation(u_prop(time * t_s), u_remi(time * t_s), u_nore(time * t_s), u_rocu(time * t_s));
        end

        simulator.save_simulation()

        model_prop_pk_str = char(pk_models('prop'));
        model_prop_pd_str = char(pd_models('prop'));  

        % path = './csv/infusion/';
        % filename = strcat('pk_', model_prop_pk_str, '_pd_', model_prop_pd_str, '_', inter,'.csv');
        % simulator.save_to_csv(path, filename);
        simulator.plot_simulation();
    end

    function complete_simulation_TCI(stimuli,volume_status)
        % Complete simulation of a patient with a complete sequence of propofol and remifentanil infusion rates.

        simulator = Factory.createSimulator(SimulatorMode.Concentration);
       
        for k = 1:44
            disp('Patient id:' + string(k));
            interaction = Interaction.Surface;
            dohMeasure = DoHMeasure.Both;
            vs = volume_status;
            modes_TCI =  containers.Map({'prop', 'remi'}, {TciMode.EffectSite, TciMode.EffectSite});

            limits_TCI = containers.Map({'cp_limit_prop', 'cp_limit_remi', ...
                'infusion_limit_prop', 'infusion_limit_remi'}, {10, 10, 2 , 0.5});

            pk_models = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Eleveld});
            pd_models = containers.Map({'prop', 'remi'}, {Model.PATIENT_SPECIFIC, Model.Eleveld});

            t_sim = 30 * 60;
            t_s = 5;

            var = {'opiates', true, 'blood_sampling', BloodSampling.ARTERIAL,...
                'interaction',interaction, 'doh_measure', dohMeasure,...
                'stimuli', stimuli, 'modes_TCI',modes_TCI,'volume_status', vs,...
                'pk_models', pk_models, 'pd_models', pd_models,...
                'pk_models_TCI', pk_models, 'pd_models_TCI', pd_models,...
                'limits_TCI',limits_TCI};

            
            simulator.init_simulation_from_file(k,t_sim,t_s,var{:});

            target_prop = 2.0*ones(t_sim, 1);
            target_remi = 3.6*ones(t_sim, 1);
            target_nore = 0.5*ones(t_sim, 1);
            target_rocu = 2*ones(t_sim, 1);
        
            simulator.run_complete_simulation(target_prop, target_remi, target_nore, target_rocu);
            simulator.save_simulation()

        end
        
        % initial letter of the pk_mode:
        if interaction == Interaction.Surface
            inter = 'Yes';
        elseif interaction == Interaction.No_interaction
            inter = 'No';
        else
            error('Invalid interaction type');
        end


        % first letter of the pk_mode:
        model_prop_pk_str = char(pk_models('prop'));
        pd_model_remi_pk_str = char(pk_models('remi'));
        model_prop_pd_str = char(pd_models('prop'));
        pd_model_remi_pd_str = char(pd_models('remi'));

        path = './csv/concentration/';

        filename = strcat('pk_', model_prop_pk_str(1),'_', pd_model_remi_pk_str(1), '_pd_', model_prop_pd_str(1), '_', pd_model_remi_pd_str(1), '_', inter, '.csv');
        simulator.save_to_csv(path, filename);
        
        simulator.plot_simulation();

    end

    function step_by_step_simulation_TCI(stimuli,volume_status)
        % Step by step simulation of a patient with a single value of propofol and remifentanil infusion rates at each step.

        simulator = Factory.createSimulator(SimulatorMode.Concentration);
        patient_id = 16;
        
        disp('Patient id:' + string(patient_id));

        interaction = Interaction.Surface;   
        dohMeasure = DoHMeasure.Both;
        vs = volume_status;
        modes_TCI =  containers.Map({'prop', 'remi'}, {TciMode.EffectSite, TciMode.EffectSite});

        limits_TCI = containers.Map({'cp_limit_prop', 'cp_limit_remi',...
            'cp_limit_nore', 'cp_limit_rocu','infusion_limit_prop', ...
            'infusion_limit_remi', 'infusion_limit_nore', 'infusion_limit_rocu'},...
            {10, 10, 15, 10, 2, 0.5, 0.5, 0.1});

        pk_models = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Minto});
        pd_models = containers.Map({'prop', 'remi'}, {Model.PATIENT_SPECIFIC, Model.Minto});
        pk_models_TCI = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Minto});
        pd_models_TCI = containers.Map({'prop', 'remi'}, {Model.Eleveld, Model.Minto});
        
        % Manueally initialize the outputs
        output_init = containers.Map({'doh','map', 'co','hr','sv'},{90,150,5.6,80,70});

        t_sim = 30 * 60;
        t_s = 5;
        var = {'opiates', true, 'blood_sampling', BloodSampling.ARTERIAL, 'interaction',...
            interaction, 'doh_measure', dohMeasure, 'stimuli', stimuli, 'modes_TCI',...
            modes_TCI, 'volume_status', vs, 'pk_models', pk_models, 'pd_models',...
            pd_models, 'pk_models_TCI', pk_models_TCI, 'pd_models_TCI', ...
            pd_models_TCI, 'limits_TCI',limits_TCI, 'output_init',output_init};

        simulator.init_simulation_from_file(patient_id,t_sim,t_s,var{:});

        target_prop = 2*ones(t_sim, 1);
        target_remi = 3.6*ones(t_sim, 1);
        target_nore = 0.5*ones(t_sim, 1) * 0;
        target_rocu = 1.5*ones(t_sim, 1);
        
        steps = t_sim / t_s;
        
        for time = 1:steps
            simulator.one_step_simulation(target_prop(time * t_s), target_remi(time * t_s), target_nore(time * t_s), target_rocu(time * t_s));
        end

        simulator.save_simulation()
        
        if interaction == Interaction.Surface
            inter = 'Yes';
        elseif interaction == Interaction.No_interaction
            inter = 'No';
        else
            error('Invalid interaction type');
        end

        % first letter of the pk_mode:
        model_prop_pk_str = char(pk_models('prop'));
        pd_model_remi_pk_str = char(pk_models('remi'));
        model_prop_pd_str = char(pd_models('prop'));
        pd_model_remi_pd_str = char(pd_models('remi'));

        path = './csv/concentration/';
        filename = strcat('pk_', model_prop_pk_str(1),'_', pd_model_remi_pk_str(1), '_pd_', model_prop_pd_str(1), '_', pd_model_remi_pd_str(1), '_', inter, '.csv');
        simulator.save_to_csv(path, filename);

        simulator.plot_simulation();
    end

    function compare_files(file1, file2)
        % Compare two files line by line and print the first difference found.

        fid1 = fopen(file1, 'r');
        fid2 = fopen(file2, 'r');

        if fid1 == -1 || fid2 == -1
            error('Error opening one or both files.');
        end

        line_num = 1;

        while true
            line1 = fgetl(fid1);
            line2 = fgetl(fid2);

            % Check if both files have reached the end
            if ~ischar(line1) && ~ischar(line2)
                disp('Files are identical.');
                break;
            end

            if ~strcmp(line1, line2)
                fprintf('Difference found at line %d:\nFile 1: %s\nFile 2: %s\n', line_num, line1, line2);
                break;
            end

            line_num = line_num + 1;
        end

        % Check if one file has more lines than the other
        if ischar(line1) || ischar(line2)
            disp('Files have different lengths.');
        end

        fclose(fid1);
        fclose(fid2);
    end


% Add all the necessary paths
addpath('enums');
addpath('parameters');

% Initialize the containers.Map
stimuli = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Add stimuli entries as structure
stimuli(5*60)  = struct('disturbanceType', DisturbanceType.INCISION,          'duration', 2 * 60,  'deltas', [10, 10, 20]);
stimuli(1*60)  = struct('disturbanceType', DisturbanceType.INTUBATION,        'duration', 2 * 60,  'deltas', [10, 10, 20]);
stimuli(10*60) = struct('disturbanceType', DisturbanceType.SKIN_MANIPULATION, 'duration', 10 *60, 'deltas', [5, 5, 10]);
stimuli(21*60) = struct('disturbanceType', DisturbanceType.SUTURE,            'duration', 5 * 60,  'deltas', [2, 2, 4]);

volume_status = containers.Map('KeyType', 'double', 'ValueType', 'any');
volume_status(200) = VolumeStatus.Hypovolemia;
volume_status(800) = VolumeStatus.Normovolemia;
volume_status(1200) = VolumeStatus.Hypovolemia;

% stimuli = [];
% volume_status = [];

% step_by_step_simulation(stimuli, volume_status);
complete_simulation(stimuli, volume_status);
% complete_simulation_TCI(stimuli, volume_status);
% step_by_step_simulation_TCI(stimuli, volume_status);

%% Use the following lines to compare results of two simulations.
% file1 = 'results_of_simulation_1.csv';
% file2 = 'results_of_simulation_2.csv';
% compare_files(file1, file2);

end





