% matlab version of the disturbance class

classdef Disturbance < handle
    properties (Access = private)
        t_sim
        doh_values
        co_values

        disturbances
        min_dis_doh = 0.0
        min_dis_co = 0.0
        starts

        C50_prop_intubation = 8.69
        C50_rem_intubation = 4.95
        gamma_intubation = 3.22
        epsilon_intubation = 0.11

    end
    
    methods
        function obj = Disturbance(t_sim, disturbances)
        % Initialize the Disturbance class
        % Parameters:
        %   t_sim: the total simulation time
        %   disturbances: a cell array containing the disturbances
        % Returns:
        %   obj: the Disturbance object

            obj.t_sim = t_sim;
            obj.doh_values = zeros(t_sim, 1);
            obj.co_values = zeros(t_sim, 1);

            obj.disturbances = disturbances;
            obj.starts = containers.Map('KeyType', 'double', 'ValueType', 'any');

            obj.validate_sequence_disturbances();
        end
        
        function [doh_values, co_values] = get_disturbances(obj, start, cp_prop, cp_remi)
        % Get the disturbances, if the time point 
        % Parameters:
        %   time: the time point to get the disturbances
        %   cp_prop: plasma concentration of propofol
        %   cp_remi: plasma concentration of remifentanil
        % Returns:
        %   doh_values: the disturbances for depth of hypnosis
        %   co_values: the disturbances for cardiac output
        
            if isKey(obj.disturbances, start)
                type = obj.disturbances(start).disturbanceType;
                duration = obj.disturbances(start).duration;
                deltas = obj.disturbances(start).deltas;
                obj.min_dis_doh = 0.1 * deltas(1);
                obj.min_dis_co = 0.1 * deltas(2);
                w = obj.generate_coeff(cp_prop, cp_remi);
                delta_doh = obj.min_dis_doh + (deltas(1) - obj.min_dis_doh) * w;
                delta_co = obj.min_dis_co + (deltas(2) - obj.min_dis_co) * w;
                deltas = [delta_doh, delta_co];
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
            co_values = obj.co_values;
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
            % [true, false] with probability [p, 1-p]
            % 1 specifies the number of samples to return
            % true if the sample is with replacement [don't care]
            if ~randsample([true, false], 1, true, [p, 1 - p])
                coeff = 0;
            else
                u = rand();
                coeff = obj.inverse_cdf_custom(u, p);
            end
        end

        function compute_disturbance_intubation(obj, start, duration, delta_values)
        % Compute the disturbances for intubation event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbance values for depth of hypnosis and cardiac output

            sigma = 0.2 + (0.2 * rand() - 0.1);
            mode = duration;
            mu = log(mode) + sigma ^ 2;
        
            upper_range = ceil(exp(mu + 4 * sigma));

            time = linspace(start, start + upper_range, upper_range);
            end_index = min(upper_range, obj.t_sim - start);
            lognorm = lognpdf(time-start, mu, sigma);

            peak_doh = lognorm / max(lognorm) * delta_values(1);
            peak_co = lognorm / max(lognorm) * delta_values(2);

            peak_doh = peak_doh(1:end_index); 
            peak_co = peak_co(1:end_index);
           
            obj.doh_values(start:start + end_index -1) = obj.doh_values(start:start + end_index-1) + peak_doh' ;
            obj.co_values(start:start + end_index-1) = obj.co_values(start:start + end_index-1) + peak_co';
        end


        function compute_disturbance_incision(obj, start, duration, delta_values)
        % Compute the disturbances for incision event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbances values for depth of hypnosis and cardiac output
            
            peak_time = duration;
            time = linspace(start, start + duration, duration);
            sigma = 0.1 * duration * 2;
            gaussian = normpdf(time-start, peak_time, sigma);

            obj.doh_values(start:start + duration-1) = obj.doh_values(start:start + duration-1) + (gaussian / max(gaussian) * delta_values(1))';
            obj.co_values(start:start + duration-1) = obj.co_values(start:start + duration-1) + (gaussian / max(gaussian) * delta_values(2))';
            
            decay_start = start + duration -1;
            total_time = 3600 * 60;
            inflection_point = 3600 * 48;
            time = linspace(decay_start, decay_start + total_time, total_time);
            steepness = 5;
            s_curve = 0.5 * (1 -tanh((time - (decay_start + inflection_point)) / (inflection_point / steepness)));
            
            obj.doh_values(decay_start:obj.t_sim-1) = obj.doh_values(decay_start:obj.t_sim-1) + (s_curve(1:obj.t_sim - decay_start) * delta_values(1))';
            obj.co_values(decay_start:obj.t_sim-1) = obj.co_values(decay_start:obj.t_sim-1) + (s_curve(1:obj.t_sim - decay_start) * delta_values(2))';
            obj.doh_values(decay_start) = delta_values(1);
            obj.co_values(decay_start) = delta_values(2);

        end

        function compute_disturbance_skin_manipulation(obj, start, duration, delta_values, w)
        % Compute the disturbances for skin manipulation event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the range of the noise for the disturbances values for depth of hypnosis and cardiac output
        %   w: the weight of the noise

            disturbance = zeros(duration, 1);

            last_index = 0;
            last_duration = 0;

            if start == 0
                baseline_doh = obj.doh_values(start);
                baseline_co = obj.co_values(start);
            else
                baseline_doh = obj.doh_values(start - 1);
                baseline_co = obj.co_values(start - 1);
            end

            random_starts  = rand(1,obj.t_sim)< 0.5;
            random_modifiers = rand(1,obj.t_sim) < 0.05;

            for i = 1:duration
                if disturbance(i) == 0 && random_starts(i) == 1 
                    max_duration = min(max(duration - i, 1), 180);

                    if max_duration > 15
                        duration_disturbance = randi([15, max_duration]);
                    elseif max_duration == 1
                        duration_disturbance = 1;
                    else
                        duration_disturbance = randi([1, max_duration]);
                    end

                    end_index = min(i + duration_disturbance - 1, duration - 1);
                    disturbance(i:end_index) = 1;
                    last_index = i;
                    last_duration = duration_disturbance;

                    noise_doh = rand() * delta_values(1);
                    noise_doh = obj.min_dis_doh + (noise_doh - obj.min_dis_doh) * w;
                    noise_co = rand() * delta_values(2);
                    noise_co = obj.min_dis_co + (noise_co - obj.min_dis_co) * w;

                    peak_time = duration_disturbance / 2;
                    time = linspace(0, duration_disturbance, duration_disturbance);
                    sigma = 0.2 * duration_disturbance;
                    gaussian = normpdf(time, peak_time, sigma);
                    gaussian_normalized = gaussian / max(gaussian);

                    end_global = start + end_index;
                    obj.doh_values(start + i:end_global) = obj.doh_values(start + i:end_global) + (gaussian_normalized * noise_doh)';
                    obj.co_values(start + i:end_global) = obj.co_values(start + i:end_global) + (gaussian_normalized * noise_co)';

                elseif disturbance(i) == 1 && random_modifiers(i) == 1
                    duration_disturbance = max(last_duration - (i - last_index), 1);
                    end_idx = min(i + duration_disturbance-1, duration-1);
                    disturbance(i:end_idx) = 2;  % Mark as modified
                   
                    max_additional_doh = (baseline_doh + delta_values(1)) - max(obj.doh_values(start + i:start + end_idx));
                    max_decrease_doh = max(obj.doh_values(start + i:start + end_idx) - baseline_doh);
                    noise_doh = (max_additional_doh - (-max_decrease_doh)) * rand() + (-max_decrease_doh);
                    
                    if noise_doh < 0
                        if w ~= 0
                            noise_doh = -obj.min_dis_doh + (noise_doh + obj.min_dis_doh) * w;
                            max_decrease_co = max(obj.co_values(start + i:start + end_idx)) - baseline_co;
                            noise_co = (-max_decrease_co - 0) * rand() + 0;
                            noise_co = -obj.min_dis_co + (noise_co + obj.min_dis_co) * w;
                        else
                            noise_doh = -obj.min_dis_doh;
                            noise_co = -obj.min_dis_co;
                        end
                    else
                        if w ~= 0
                            noise_doh = obj.min_dis_doh + (noise_doh - obj.min_dis_doh) * w;
                            max_additional_co = (baseline_co + delta_values(2)) - max(obj.co_values(start + i:start + end_idx));
                            noise_co = (max_additional_co - 0) * rand() + 0;
                            noise_co = obj.min_dis_co + (noise_co - obj.min_dis_co) * w;
                        else
                            noise_doh = obj.min_dis_doh;
                            noise_co = obj.min_dis_co;
                        end
                    end

                    peak_time = duration_disturbance / 2;
                    time = linspace(0, duration_disturbance, duration_disturbance);
                    sigma = 0.2 * duration_disturbance;
                    gaussian = normpdf(time, peak_time, sigma);
                    gaussian_normalized = gaussian / max(gaussian);
                    
                    obj.doh_values(start + i:start + end_idx) = obj.doh_values(start + i:start + end_idx) + (gaussian_normalized * noise_doh)';
                    obj.co_values(start + i:start + end_idx) = obj.co_values(start + i:start + end_idx) + (gaussian_normalized * noise_co)';
                    obj.doh_values(start + i:start + end_idx) = max(min(obj.doh_values(start + i:start + end_idx), baseline_doh + delta_values(1)), baseline_doh);
                    obj.co_values(start + i:start + end_idx) = max(min(obj.co_values(start + i:start + end_idx), baseline_co + delta_values(2)), baseline_co);
                end
            end
        end

        function compute_disturbance_suture(obj, start, duration, delta_values)
        % Compute the disturbances for suture event
        % Parameters:
        %   start: the start time of the event
        %   duration: the duration of the event
        %   delta_values: the disturbances values for depth of hypnosis, mean arterial pressure, cardiac output
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
        end 
    end 
end 
            