% matlab version of the disturbance class
classdef Disturbance < handle
    properties (Access = private)
        t_sim
        doh_values
        hr_values
        map_values

        disturbances
        min_dis_doh = 0.0
        min_dis_hr = 0.0
        min_dis_map = 0.0
        starts

        C50_prop_intubation = 8.69
        C50_rem_intubation = 4.95
        gamma_intubation = 3.22
        epsilon_intubation = 0.11

        worst_case = false
    end
    
    methods
        function obj = Disturbance(t_sim, disturbances, seed, worst_case)
        % Initialize the Disturbance class
        % Parameters:
        %   t_sim: the total simulation time
        %   disturbances: a cell array containing the disturbances
        %   seed: the seed for the random number generator
        %   worst_case: a boolean indicating if the worst case scenario is used
        % Returns:
        %   obj: the Disturbance object

            if ~isempty(seed)
                % disp('Using seed');
                % disp(seed);
                rand('twister', seed);
            end

            obj.t_sim = t_sim;
            obj.doh_values = zeros(t_sim, 1);
            obj.hr_values = zeros(t_sim, 1);
            obj.map_values = zeros(t_sim, 1);

            obj.disturbances = disturbances;
            obj.starts = containers.Map('KeyType', 'double', 'ValueType', 'any');

            % % Distinguish the validation of the sequence of disturbances between induction and maintenance
            % obj.in_maintenance = in_maintenance;

            % Validate the sequence of disturbances
            obj.validate_sequence_disturbances();

            obj.worst_case = worst_case;

        end
        
        function [doh_values, hr_values, map_values] = get_disturbances(obj, start, cp_prop, cp_remi)
        % Get the disturbances, if the time point 
        % Parameters:
        %   time: the time point to get the disturbances
        %   cp_prop: plasma concentration of propofol
        %   cp_remi: plasma concentration of remifentanil
        % Returns:
        %   doh_values: the disturbances for depth of hypnosis
        %   hr_values: the disturbances for heart rate
        %   map_values: the disturbances for mean arterial pressure
        
            if isKey(obj.disturbances, start)
                type = obj.disturbances(start).disturbanceType;
                duration = obj.disturbances(start).duration;
                deltas = obj.disturbances(start).deltas;
                obj.min_dis_doh = 0.2 * deltas(1);
                obj.min_dis_hr = 0.2 * deltas(2);
                obj.min_dis_map = 0.2* deltas(3);
                if obj.worst_case
                    w = 1;
                else
                    w = obj.generate_coeff(cp_prop, cp_remi);
                end
                delta_doh = obj.min_dis_doh + (deltas(1) - obj.min_dis_doh) * w;
                delta_hr = obj.min_dis_hr + (deltas(2) - obj.min_dis_hr) * w;
                delta_map = obj.min_dis_map + (deltas(3) - obj.min_dis_map) * w;
                deltas = [delta_doh, delta_hr, delta_map];

                start = start + 1;

                if type == DisturbanceType.INTUBATION
                    obj.compute_disturbance_intubation(start, duration, deltas);
                elseif type == DisturbanceType.INCISION
                    obj.compute_disturbance_incision(start, duration, deltas);
                elseif type == DisturbanceType.SKIN_MANIPULATION
                    obj.compute_disturbance_skin_manipulation(start, duration, deltas, w);
                elseif type == DisturbanceType.SUTURE
                    obj.compute_disturbance_suture(start, duration, deltas);
                end
            end
            
            doh_values = obj.doh_values;
            hr_values = obj.hr_values;
            map_values = obj.map_values;
        end

        function prob = compute_responce_probability_I(c_prop, c_remi, C50_prop, C50_rem, gamma, epsilon)
        % Compute the probability of no response to a stimulus.
        % Parameters:
        %   c_prop: plasma concentration of propofol
        %   c_remi: plasma concentration of remifentanil
        %   C50_prop: C50 of propofol
        %   C50_rem: C50 of remifentanil
        %   gamma: parameter that describes the steepness of the dose-response curve
        %   epsilon: interaction parameter

            num = (c_prop / C50_prop + c_remi / C50_rem + epsilon * c_prop / C50_prop * c_remi / C50_rem) ^ gamma;
            den = 1 + (c_prop / C50_prop + c_remi / C50_rem + epsilon * c_prop / C50_prop * c_remi / C50_rem) ^ gamma;
            prob = num / den;
        end

        function prob = compute_responce_probability_II(c_prop, c_remi, C50_prop, gamma, epsilon)
        % Compute the probability of no response to a stimulus.
        % Parameters:
        %   c_prop: plasma concentration of propofol
        %   c_remi: plasma concentration of remifentanil
        %   C50_prop: C50 of propofol
        %   gamma: parameter that describes the steepness of the dose-response curve
        %   epsilon: interaction parameter

            num = (c_prop / C50_prop + epsilon * c_prop / C50_prop * c_remi) ^ gamma;
            den = 1 + (c_prop / C50_prop + epsilon * c_prop / C50_prop * c_remi) ^ gamma;
            prob =  num / den;
        end

        function prob = compute_response_probability_III(~,c_prop, c_remi, gamma, epsilon)
        % Compute the probability of no response to a stimulus.
        % Parameters:
        %   c_prop: plasma concentration of propofol
        %   c_remi: plasma concentration of remifentanil
        %   gamma: parameter that describes the steepness of the dose-response curve

            num = (epsilon * c_prop * c_remi) ^ gamma;
            den = 1 + (epsilon * c_prop * c_remi) ^ gamma;
            prob = num / den;
        end

        function pdf = pdf_custom(x, p)
        % Custom PDF defined as:
        % f(x) = (2 * (1 - p) / p^2) * x,      for 0 <= x <= p
        % f(x) = (2 * p / (1 - p)**2) * (1 - x), for p < x <= 1
        % f(x) = 0, otherwise

            if x >= 0 && x <= p
                pdf = (2 * (1 - p) / p ^ 2) * x;
            elseif x > p && x <= 1
                pdf = (2 * p / (1 - p) ^ 2) * (1 - x);
            end
            
        end

        function icdf = inverse_cdf_custom(~,u, p)
        % Analytical inverse CDF for the custom distribution.
        % Parameters:
        %   u: uniform random number in [0, 1]
        %   p: probability of response

            if u >= 0 && u <= (1-p)
                icdf = p * sqrt(u / (1 - p));
            elseif u >= (1 - p) && u <= 1
                k = p - (p ^ 2) / 2 + ((1 - p) ^ 2) / (2 * p) * (u - (1 - p));
                icdf = 1 - sqrt(1 - 2 * k);
            else
                icdf = 0;
            end
        end

        function coeff = generate_coeff(obj, cp_prop, cp_remi)
        % Generates one sample from the custom distribution using the inverse CDF.
        % Parameters:
        %   cp_prop: plasma concentration of propofol
        %   cp_remi: plasma concentration of remifentanil
            prob_no_response = obj.compute_response_probability_III(cp_prop, cp_remi, obj.gamma_intubation, obj.epsilon_intubation);
            p = 1 - prob_no_response;
            % disp('p');
            % disp(p);

            % [true, false] with probability [p, 1-p]
            % 1 specifies the number of samples to return
            % true if the sample is with replacement [don't care]
            toss = rand() < p; % true with probability p
            if ~toss
                coeff = 0;
                % disp('coeff');
                % disp(coeff);
            else
                u = rand();
                % disp('u');
                % disp(u)
                coeff = obj.inverse_cdf_custom(u, p);
                % disp('coeff');
                % disp(coeff);
            end
        end

        function compute_disturbance_intubation(obj, start, duration, delta_values)
        % Compute the disturbances for intubation event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbance values for depth of hypnosis, heart rate, and mean arterial pressure

            sigma = 0.2 + ((0.1-(-0.1)) * rand() - 0.1);
            mode = duration;
            mu = log(mode) + sigma ^ 2;
        
            upper_range = ceil(exp(mu + 4 * sigma));

            end_idx = min(upper_range, obj.t_sim - start);
            time = linspace(start, start + end_idx, end_idx + 1);
            lognorm = lognpdf(time-start, mu, sigma);

            peak_doh = lognorm / max(lognorm) * delta_values(1);
            peak_hr = lognorm / max(lognorm) * delta_values(2);
            peak_map = lognorm / max(lognorm) * delta_values(3);

            range = start:(start + end_idx);
            obj.doh_values(range) = obj.doh_values(range) + peak_doh' ;
            obj.hr_values(range) = obj.hr_values(range) + peak_hr';
            obj.map_values(range) = obj.map_values(range) + peak_map';

            % disp('range')
            % disp(start)
            % disp(start + end_idx)
        end


        function compute_disturbance_incision(obj, start, duration, delta_values)
        % Compute the disturbances for incision event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbances values for depth of hypnosis, heart rate, and mean arterial pressure
            
            peak_time = duration;
            time = linspace(start, start + duration, duration + 1);
            sigma = 0.1 * duration * 2;
            gaussian = normpdf(time-start, peak_time, sigma);

            range = start:(start + duration);
            obj.doh_values(range) = obj.doh_values(range) + (gaussian / max(gaussian) * delta_values(1))';
            obj.hr_values(range) = obj.hr_values(range) + (gaussian / max(gaussian) * delta_values(2))';
            obj.map_values(range) = obj.map_values(range) + (gaussian / max(gaussian) * delta_values(3))';

            decay_start = start + duration;
            total_time = 3600 * 60;
            inflection_point = 3600 * 48;
            time = linspace(decay_start, decay_start + total_time, total_time + 1);
            steepness = 5;
            s_curve = 0.5 * (1 -tanh((time - (decay_start + inflection_point)) / (inflection_point / steepness)));
            
            obj.doh_values(decay_start:obj.t_sim) = obj.doh_values(decay_start:obj.t_sim) + (s_curve(1:obj.t_sim - decay_start +1) * delta_values(1))';
            obj.hr_values(decay_start:obj.t_sim) = obj.hr_values(decay_start:obj.t_sim) + (s_curve(1:obj.t_sim - decay_start + 1) * delta_values(2))';
            obj.map_values(decay_start:obj.t_sim) = obj.map_values(decay_start:obj.t_sim) + (s_curve(1:obj.t_sim - decay_start + 1) * delta_values(3))';
            obj.doh_values(decay_start) = delta_values(1);
            obj.hr_values(decay_start) = delta_values(2);
            obj.map_values(decay_start) = delta_values(3);

            % disp('range')
            % disp(decay_start)
            % disp(obj.t_sim)

        end

        function compute_disturbance_skin_manipulation(obj, start, duration, delta_values, w)
                % Compute the disturbances for skin manipulation event
                % Parameters:
                %   start: the start time of the event
                %   duration: the duration of the event
                %   delta_values: the range of the noise for the disturbances values for depth of hypnosis, heart rate, and mean arterial pressure
                %   w: the weight of the noise

                disturbance = zeros(duration, 1);

                last_index = -1;
                last_duration = 0;

                % Set baseline values
                baseline_doh = obj.doh_values(start);
                baseline_hr = obj.hr_values(start);
                baseline_map = obj.map_values(start);
 
                % Random start and modifier generation
                random_starts = false(1, duration);
                random_modifiers = false(1, duration);

                for i = 1:duration
                    if rand() < 0.5
                        random_starts(i) = true;
                    end
                end

                for i = 1:duration
                    if rand() < 0.05
                        random_modifiers(i) = true;
                    end
                end

                % disp('Random Starts')
                % disp(random_starts)
                % disp('Random Modifiers')
                % disp(random_modifiers)

                for i = 1:duration
                    j = i - 1; % Adjust for MATLAB's 1-based indexing

                    if disturbance(i) == 0 && random_starts(i)
                        max_possible_duration = min(max(duration - j, 1), 180);

                        if max_possible_duration > 15
                            duration_disturbance = 15 + floor(rand() * (max_possible_duration - 15 + 1));
                        elseif max_possible_duration == 1
                            duration_disturbance = 1;
                        else
                            duration_disturbance = 1 + floor(rand() * (max_possible_duration - 1 + 1));
                        end

                        end_idx = min(j + duration_disturbance , duration);  % Ensure proper end index
                        disturbance(i:end_idx) = 1;
                        last_index = j;
                        last_duration = duration_disturbance;

                        noise_doh = rand() * delta_values(1);
                        noise_doh = obj.min_dis_doh + (noise_doh - obj.min_dis_doh) * w;
                        noise_hr = rand() * delta_values(2);
                        noise_hr = obj.min_dis_hr + (noise_hr - obj.min_dis_hr) * w;
                        noise_map = rand() * delta_values(3);
                        noise_map = obj.min_dis_map + (noise_map - obj.min_dis_map) * w;

                        % disp('Noise doh')
                        % disp(noise_doh)
                        % disp('Noise hr')
                        % disp(noise_hr)
                        % disp('Noise map')
                        % disp(noise_map)

                        peak_time = duration_disturbance / 2;
                        time = linspace(0, duration_disturbance, duration_disturbance + 1);
                        sigma = 0.2 * duration_disturbance;
                        gaussian = normpdf(time, peak_time, sigma);
                        gaussian_normalized = gaussian / max(gaussian);

                        range = (start + j):(start + end_idx);
                        obj.doh_values(range) = obj.doh_values(range) + (gaussian_normalized * noise_doh)';
                        obj.hr_values(range) = obj.hr_values(range) + (gaussian_normalized * noise_hr)';
                        obj.map_values(range) = obj.map_values(range) + (gaussian_normalized * noise_map)';

                        % disp('range - case 1:')
                        % disp(start + j)
                        % disp(start + end_idx)

                    elseif disturbance(i) == 1 && random_modifiers(i)
                        duration_disturbance = max(last_duration - (j - last_index), 1);
                        end_idx = min(j + duration_disturbance , duration);  
                        disturbance(i:end_idx) = 2;

                        range = (start + j):(start + end_idx);

                        max_additional_doh = (baseline_doh + delta_values(1)) - max(obj.doh_values(range));
                        max_decrease_doh = max(obj.doh_values(range) - baseline_doh);
                        noise_doh = (max_additional_doh + max_decrease_doh) * rand() - max_decrease_doh;

                        if noise_doh < 0
                            noise_doh = -obj.min_dis_doh + (noise_doh + obj.min_dis_doh) * w;
                            max_decrease_hr = max(obj.hr_values(range)) - baseline_hr;
                            noise_hr = (0 + max_decrease_hr) * rand() - max_decrease_hr;
                            noise_hr = -obj.min_dis_hr + (noise_hr + obj.min_dis_hr) * w;
                            max_decrease_map = max(obj.map_values(range)) - baseline_map;
                            noise_map = (0 + max_decrease_map) * rand() - max_decrease_map;
                            noise_map = -obj.min_dis_map + (noise_map + obj.min_dis_map) * w;
                        else
                            noise_doh = obj.min_dis_doh + (noise_doh - obj.min_dis_doh) * w;
                            max_additional_hr = (baseline_hr + delta_values(2)) - max(obj.hr_values(range));
                            noise_hr = max_additional_hr * rand();
                            noise_hr = obj.min_dis_hr + (noise_hr - obj.min_dis_hr) * w;
                            max_additional_map = (baseline_map + delta_values(3)) - max(obj.map_values(range));
                            noise_map = max_additional_map * rand();
                            noise_map = obj.min_dis_map + (noise_map - obj.min_dis_map) * w;
                        end

                        peak_time = duration_disturbance / 2;
                        time = linspace(0, duration_disturbance, duration_disturbance + 1);
                        sigma = 0.2 * duration_disturbance;
                        gaussian = normpdf(time, peak_time, sigma);
                        gaussian_normalized = gaussian / max(gaussian);

                        obj.doh_values(range) = obj.doh_values(range) + (gaussian_normalized * noise_doh)';
                        obj.hr_values(range) = obj.hr_values(range) + (gaussian_normalized * noise_hr)';
                        obj.map_values(range) = obj.map_values(range) + (gaussian_normalized * noise_map)';

                        obj.doh_values(range) = max(min(obj.doh_values(range), baseline_doh + delta_values(1)), baseline_doh);
                        obj.hr_values(range) = max(min(obj.hr_values(range), baseline_hr + delta_values(2)), baseline_hr);
                        obj.map_values(range) = max(min(obj.map_values(range), baseline_map + delta_values(3)), baseline_map);
                        

                        % disp('Noise doh')
                        % disp(noise_doh)
                        % disp('Noise hr')
                        % disp(noise_hr)
                        % disp('Noise map')
                        % disp(noise_map)

                        % disp('range - case 2:')
                        % disp(start + j)
                        % disp(start + end_idx)

                    elseif disturbance(i) == 2
                        % Do nothing
                    end
                end
        end



        function compute_disturbance_suture(obj, start, duration, delta_values)
        % Compute the disturbances for suture event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbances values for depth of hypnosis, heart rate, and mean arterial pressure
        %   cp_prop: plasma concentration of propofol
        %   cp_remi: plasma concentration of remifentanil

            obj.compute_disturbance_intubation(start, duration, delta_values);
        end

        function validate_sequence_disturbances(obj)
        % Validate the sequence of disturbances

            starts = cell2mat(obj.disturbances.keys());
            [sorted_starts, ~] = sort(starts); 
    

            if obj.disturbances(sorted_starts(1)).disturbanceType ~= DisturbanceType.INTUBATION
                error('The first disturbance must be an intubation event.');
            end

            first_incision_time = inf;
            last_incision_time = -inf;
            incision_count = 0;
            suture_count = 0;
            intubation_count = 0;

            for i = 1:length(obj.disturbances)
                % Check that any skin manipulation occurs after at least the end of one incision          
                % Check that suture events occur after incision events
                % At every suture must correspond a different incision 
                start = sorted_starts(i);
                
                if i ~= length(obj.disturbances)
                    
                    if start + obj.disturbances(start).duration > sorted_starts(i+1)
                        error('Overlapping disturbances detected.');
                    end 
                end

                if obj.disturbances(start).disturbanceType == DisturbanceType.INCISION
                    incision_count = incision_count + 1;
                    if first_incision_time == inf
                        first_incision_time = start;
                    end
                    last_incision_time = start + obj.disturbances(start).duration;
                elseif obj.disturbances(start).disturbanceType == DisturbanceType.SUTURE
                    if incision_count < suture_count || start < last_incision_time
                        error('Suture must occur after an incision event.');
                    end
                    suture_count = suture_count + 1;
                elseif  obj.disturbances(start).disturbanceType == DisturbanceType.SKIN_MANIPULATION
                    if start < first_incision_time
                        error('Skin manipulation must occur after an incision event.');
                    end
                elseif  obj.disturbances(start).disturbanceType == DisturbanceType.INTUBATION
                    intubation_count = intubation_count + 1;
                    if intubation_count > 1
                        error('There must be exactly one intubation event.');
                    end
                end
            end
            if incision_count ~= suture_count
                error('There must be the same number of incision and suture events.')
            end
        end 
    end 
end 
            