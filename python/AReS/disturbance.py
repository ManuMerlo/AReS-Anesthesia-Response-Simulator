from scipy import stats

import numpy as np

from .utils.enums import DisturbanceType

class Disturbance:
    def __init__(self, t_sim, disturbances=None, in_maintenance=False, seed=None, worst_case=False):
        """
        Initialize the Disturbance class.
        :param t_sim: simulation time.
        :type t_sim: int
        :param disturbances: dict containing the disturbances.
        :type disturbances: dict
        :param in_maintenance: flag to distinguish the validation of the sequence of disturbances between induction and maintenance.
        :type in_maintenance: bool
        :param seed: seed for reproducibility.
        :type seed: int
        :param worst_case: flag to indicate if the worst case should be used.

        """
        self._t_sim = t_sim
        self.doh_values = np.zeros(self._t_sim)
        self.hr_values = np.zeros(self._t_sim)
        self.map_values = np.zeros(self._t_sim)

        self.disturbances = disturbances
        self._min_dis_doh = 0.0
        self._min_dis_hr = 0.0
        self._min_dis_map = 0.0

        # Distinguish the validation of the sequence of disturbances between induction and maintenance
        self.in_maintenance = in_maintenance

        # Validate the sequence of disturbances
        self._validate_sequence_disturbances()

        self._C50_prop_intubation = 8.69
        self._C50_rem_intubation = 4.95
        self._gamma_intubation = 3.22
        self._epsilon_intubation = 0.11

        self._worst_case = worst_case

        if seed is not None:
            np.random.seed(seed)

    def get_disturbances(self, time: int, cp_prop: float = 0, cp_remi: float = 0):
        """
        Get the disturbances, compute the values of doh and CO if time is the start of a disturbance.
        :param time: time of the simulation.
        :type time: int
        :param cp_prop: plasma concentration of propofol.
        :type cp_prop: float
        :param cp_remi: plasma concentration of remifentanil.
        :type cp_remi: float
        :return: the values of the disturbances for doh and co
            :rtype: tuple of numpy arrays
        """

        if time in self.disturbances.keys():
            dis_type, duration, delta_values = self.disturbances[time]
            self._min_dis_doh = 0.2 * delta_values[0]  # 20% of the delta value
            self._min_dis_hr = 0.2 * delta_values[1]  # 20% of the delta value
            self._min_dis_map = 0.2 * delta_values[2]  # 20% of the delta value
            p, w = self._generate_coeff(cp_prop, cp_remi) if not self._worst_case else (1, 1)
            # If the probability of no response is 1 the disturbance amplitude is the minimum disturbance
            delta_doh = self._min_dis_doh + (delta_values[0] - self._min_dis_doh) * w
            delta_hr = self._min_dis_hr + (delta_values[1] - self._min_dis_hr) * w
            delta_map = self._min_dis_map + (delta_values[2] - self._min_dis_map) * w

            if dis_type == DisturbanceType.INTUBATION:
                self._compute_disturbance_intubation(time, duration, [delta_doh, delta_hr, delta_map])
            elif dis_type == DisturbanceType.INCISION:
                self._compute_disturbance_incision(time, duration, [delta_doh, delta_hr, delta_map])
            elif dis_type == DisturbanceType.SKIN_MANIPULATION:
                self._compute_disturbance_skin_manipulation(time, duration, delta_values, w)
            elif dis_type == DisturbanceType.SUTURE:
                self._compute_disturbance_suture(time, duration, [delta_doh, delta_hr, delta_map])

        return self.doh_values, self.hr_values, self.map_values

    def _compute_response_probability_I(self, c_prop, c_remi, C50_prop, C50_rem, gamma, epsilon):
        """
        Compute the probability of no response to a stimulus.
        Based on the paper: "Propofol Reduces Perioperative Remifentanil Requirements in a Synergistic Manner"
        :param c_prop: plasma concentration of propofol.
        :param c_remi: plasma concentration of remifentanil.
        :param C50_prop: C50 of propofol.
        :param C50_rem: C50 of remifentanil.
        :param gamma: parameter that describes the steepness of the dose-response curve.
        :param epsilon: interaction parameter.
        :return: probability of no response to the stimulus.
            :rtype: float
        """
        num = (c_prop / C50_prop + c_remi / C50_rem + epsilon * c_prop / C50_prop * c_remi / C50_rem) ** gamma
        den = 1 + (c_prop / C50_prop + c_remi / C50_rem + epsilon * c_prop / C50_prop * c_remi / C50_rem) ** gamma
        return num / den

    def _compute_response_probability_II(self, c_prop, c_remi, C50_prop, gamma, epsilon):
        """
        Compute the probability of no response to a stimulus.
        Based on the paper: "Propofol Reduces Perioperative Remifentanil Requirements in a Synergistic Manner"
        :param c_prop: plasma concentration of propofol.
        :type c_prop: float
        :param c_remi: plasma concentration of remifentanil.
        :type c_remi: float
        :param C50_prop: C50 of propofol.
        :type C50_prop: float
        :param gamma: parameter that describes the steepness of the dose-response curve.
        :type gamma: float
        :param epsilon: interaction parameter.
        :type epsilon: float
        :return: probability of no response to the stimulus.
        """
        num = (c_prop / C50_prop + epsilon * c_prop / C50_prop * c_remi) ** gamma
        den = 1 + (c_prop / C50_prop + epsilon * c_prop / C50_prop * c_remi) ** gamma
        return num / den

    def _compute_response_probability_III(self, c_prop, c_remi, gamma, epsilon):
        """
        Compute the probability of no response to a stimulus.
        :param c_prop: plasma concentration of propofol.
        :type c_prop: float
        :param c_remi: plasma concentration of remifentanil.
        :type c_remi: float
        :param gamma: parameter that describes the steepness of the dose-response curve.
        :type gamma: float
        :param epsilon: interaction parameter.
        :type epsilon: float
        :return: probability of no response to the stimulus.
        """
        num = (epsilon * c_prop * c_remi) ** gamma
        den = 1 + (epsilon * c_prop * c_remi) ** gamma
        return num / den

    def pdf_custom(self, x, p):
        """
        Custom PDF defined as:
        f(x) = (2 * (1 - p) / p^2) * x,      for 0 <= x <= p
        f(x) = (2 * p / (1 - p)**2) * (1 - x), for p < x <= 1
        f(x) = 0, otherwise
        """
        pdf = np.zeros_like(x)
        mask1 = (x >= 0) & (x <= p)
        mask2 = (x > p) & (x <= 1)
        pdf[mask1] = (2 * (1 - p) / p ** 2) * x[mask1]
        pdf[mask2] = (2 * p / (1 - p) ** 2) * (1 - x[mask2])
        return pdf

    def inverse_cdf_custom(self, u, p):
        """
        Analytical inverse CDF for the custom distribution.
        :param u: uniform random number in [0, 1]
        :type u: float
        :param p: probability of response
        :type p: float
        :return: sample from the custom distribution
        """
        if 0 <= u <= (1 - p):
            x = p * np.sqrt(u / (1 - p))
        elif (1 - p) < u <= 1:
            k = p - (p ** 2) / 2 + (u - (1 - p)) * ((1 - p) ** 2) / (2 * p)
            x = 1 - np.sqrt(1 - 2 * k)
        else:
            x = 0
        return x

    def _generate_coeff(self, cp_prop, cp_remi):
        """
        Generates one sample from the custom distribution using the inverse CDF.
        :param cp_prop: plasma concentration of propofol.
        :type cp_prop: float
        :param cp_remi: plasma concentration of remifentanil.
        :type cp_remi: float
        :return: sample from the custom distribution

        """
        prob_no_response = self._compute_response_probability_III(cp_prop, cp_remi, self._gamma_intubation,
                                                                  self._epsilon_intubation)
        p = 1 - prob_no_response
        # True with probability p
        # False with probability 1 - p
        if not np.random.choice([True, False], p=[p, 1 - p]):
            return p, 0

        u = np.random.uniform(0, 1)
        coeff = self.inverse_cdf_custom(u, p)
        return p, coeff

    def _compute_disturbance_intubation(self, start, duration, delta_values):
        """
        Compute the disturbance for intubation and suture events.
        :param start: start time of the disturbance.
        :type start: int
        :param duration: duration of the disturbance.
        :type duration: int
        :param delta_values: list of delta values for doh and co.
        :type delta_values: list
        """
        # Mode = e^(mu-sigma^2)
        # mu = log(mode) + sigma^2
        # scale = e^(mu) = e^(log(t_peak) + sigma^2) = t_peak * e^(sigma^2)

        sigma = 0.2 + np.random.uniform(-0.1, 0.1)
        mode = duration
        mu = np.log(mode) + sigma ** 2

        upper_range = int(np.exp(mu + 4 * sigma))

        time = np.linspace(start, start + upper_range, upper_range)
        end = min(upper_range, self._t_sim - start)
        lognorm = stats.lognorm.pdf(time - start, sigma, scale=np.exp(mu))[0:end]

        # max = np.exp(-mu+sigma**2/2)/(sigma*np.sqrt(2*math.pi))
        peak_doh = lognorm / np.max(lognorm) * delta_values[0]
        peak_hr = lognorm / np.max(lognorm) * delta_values[1]
        peak_map = lognorm / np.max(lognorm) * delta_values[2]

        self.doh_values[start:start + end] += peak_doh
        self.hr_values[start:start + end] += peak_hr
        self.map_values[start:start + end] += peak_map

    def _compute_disturbance_incision(self, start, duration, delta_values):
        """
        Compute the disturbance for incision events.
        :param start: start time of the disturbance.
        :type start: int
        :param duration: duration of the disturbance.
        :type duration: int
        :param delta_values: list of delta values for doh and co
        :type delta_values: list
        """
        # The disturbance is modeled as a normal distribution with a peak at the end of the incision
        # The values are computed as the product of the normal distribution and the delta values
        # The disturbance decays as an S-curve after the end of the incision

        # normal distribution with the peak at the end of the incision
        peak_time = duration
        time = np.linspace(start, start + duration, duration)
        sigma = 0.1 * duration * 2
        gaussian = stats.norm.pdf(time - start, peak_time, sigma)

        self.doh_values[start:start + duration] += gaussian / np.max(gaussian) * delta_values[0]
        self.hr_values[start:start + duration] += gaussian / np.max(gaussian) * delta_values[1]
        self.map_values[start:start + duration] += gaussian / np.max(gaussian) * delta_values[2]

        # S-curve decay after the end of incision
        decay_start = start + duration
        total_time = 3600 * 60  # 60 hours in seconds
        inflection_point = 3600 * 48  # 48 hours in seconds
        time = np.linspace(decay_start, decay_start + total_time, total_time)
        # The S-curve is computed as a hyperbolic tangent function that goes from 1 to 0
        # The inflection point is at 48 hours
        steepness = 5
        s_curve = 0.5 * (1 - np.tanh((time - (decay_start + inflection_point)) / (inflection_point / steepness)))

        self.doh_values[decay_start:self._t_sim] += s_curve[:self._t_sim - decay_start] * delta_values[0]
        self.hr_values[decay_start:self._t_sim] += s_curve[:self._t_sim - decay_start] * delta_values[1]
        self.map_values[decay_start:self._t_sim] += s_curve[:self._t_sim - decay_start] * delta_values[2]

    def _compute_disturbance_skin_manipulation(self, start, duration, delta_values, w):
        """
        Compute the disturbance for skin manipulation events.
        :param start: Start time of the disturbance.
        :type start: int
        :param duration: Duration (number of time steps) of the disturbance.
        :type duration: int
        :param delta_values: List containing delta values for DOH and CO disturbances.
        :type delta_values: list of float
        :param w: Weight for the disturbance.
        :type w: float
        """

        disturbance = np.zeros(duration, dtype=int)

        last_index = -1
        last_duration = 0

        baseline_doh = self.doh_values[start - 1] if start > 0 else self.doh_values[start]
        baseline_hr = self.hr_values[start - 1] if start > 0 else self.hr_values[start]
        baseline_map = self.map_values[start - 1] if start > 0 else self.map_values[start]

        random_starts = np.random.rand(duration) < 0.5  # 50% chance to start
        random_modifies = np.random.rand(duration) < 0.05  # 5% chance to modify

        for i in range(duration):
            if disturbance[i] == 0 and random_starts[i]:
                max_possible_duration = min(max(duration - i, 1), 180)

                if max_possible_duration > 15:
                    duration_disturbance = np.random.randint(15, max_possible_duration + 1)
                elif max_possible_duration == 1:
                    duration_disturbance = 1
                else:
                    duration_disturbance = np.random.randint(1, max_possible_duration + 1)

                end_idx = min(i + duration_disturbance, duration)
                disturbance[i:end_idx] = 1  # Mark disturbance
                last_index = i
                last_duration = duration_disturbance

                noise_doh = np.random.uniform(0, delta_values[0])
                noise_doh = self._min_dis_doh + (noise_doh - self._min_dis_doh) * w
                noise_hr = np.random.uniform(0, delta_values[1])
                noise_hr = self._min_dis_hr + (noise_hr - self._min_dis_hr) * w
                noise_map = np.random.uniform(0, delta_values[2])
                noise_map = self._min_dis_map + (noise_map - self._min_dis_map) * w

                peak_time = duration_disturbance / 2
                time = np.linspace(0, duration_disturbance, duration_disturbance)
                sigma = 0.2 * duration_disturbance
                gaussian = stats.norm.pdf(time, loc=peak_time, scale=sigma)
                gaussian_normalized = gaussian / np.max(gaussian)

                end_global = start + end_idx
                self.doh_values[start + i:end_global] += gaussian_normalized * noise_doh
                self.hr_values[start + i:end_global] += gaussian_normalized * noise_hr
                self.map_values[start + i:end_global] += gaussian_normalized * noise_map

            elif disturbance[i] == 1 and random_modifies[i]:
                duration_disturbance = max(last_duration - (i - last_index), 1)
                end_idx = min(i + duration_disturbance, duration)
                disturbance[i:end_idx] = 2  # Mark as modified

                max_additional_doh = (baseline_doh + delta_values[0]) - np.max(
                    self.doh_values[start + i:start + end_idx])
                max_decrease_doh = np.max(self.doh_values[start + i:start + end_idx] - baseline_doh)
                noise_doh = np.random.uniform(-max_decrease_doh, max_additional_doh)

                if noise_doh < 0:
                    noise_doh = -self._min_dis_doh + (
                            noise_doh + self._min_dis_doh) * w if w != 0 else -self._min_dis_doh
                    max_decrease_hr = np.max(self.hr_values[start + i:start + end_idx]) - baseline_hr
                    noise_hr = np.random.uniform(-max_decrease_hr, 0)
                    noise_hr = -self._min_dis_hr + (noise_hr + self._min_dis_hr) * w if w != 0 else -self._min_dis_hr
                    max_decrease_map = np.max(self.map_values[start + i:start + end_idx]) - baseline_map
                    noise_map = np.random.uniform(-max_decrease_map, 0)
                    noise_map = -self._min_dis_map + (
                                noise_map + self._min_dis_map) * w if w != 0 else -self._min_dis_map
                else:
                    noise_doh = self._min_dis_doh + (noise_doh - self._min_dis_doh) * w if w != 0 else self._min_dis_doh
                    max_additional_hr = (baseline_hr + delta_values[1]) - np.max(
                        self.hr_values[start + i:start + end_idx])
                    noise_hr = np.random.uniform(0, max_additional_hr)
                    noise_hr = self._min_dis_hr + (noise_hr - self._min_dis_hr) * w if w != 0 else self._min_dis_hr
                    max_additional_map = (baseline_map + delta_values[2]) - np.max(
                        self.map_values[start + i:start + end_idx])
                    noise_map = np.random.uniform(0, max_additional_map)
                    noise_map = self._min_dis_map + (noise_map - self._min_dis_map) * w if w != 0 else self._min_dis_map

                peak_time = duration_disturbance / 2
                time = np.linspace(0, duration_disturbance, duration_disturbance, endpoint=False)
                sigma = 0.2 * duration_disturbance
                gaussian = stats.norm.pdf(time, loc=peak_time, scale=sigma)
                gaussian_normalized = gaussian / np.max(gaussian) if np.max(gaussian) != 0 else gaussian

                self.doh_values[start + i:start + end_idx] += gaussian_normalized * noise_doh
                self.hr_values[start + i:start + end_idx] += gaussian_normalized * noise_hr
                self.map_values[start + i:start + end_idx] += gaussian_normalized * noise_map

                self.doh_values[start + i:start + end_idx] = np.clip(self.doh_values[start + i:start + end_idx],
                                                                     baseline_doh, baseline_doh + delta_values[0])
                self.hr_values[start + i:start + end_idx] = np.clip(self.hr_values[start + i:start + end_idx],
                                                                    baseline_hr, baseline_hr + delta_values[1])
                self.map_values[start + i:start + end_idx] = np.clip(self.map_values[start + i:start + end_idx],
                                                                     baseline_map, baseline_map + delta_values[2])

            elif disturbance[i] == 2:
                continue

    def _compute_disturbance_suture(self, start, duration, delta_values):
        """
        Compute the disturbance for suture events.
        :param start: start time of the disturbance.
        :type start: int
        :param duration: duration of the disturbance.
        :type duration: int
        :param delta_values: list of delta values for doh and co.
        :type delta_values: list
        """
        self._compute_disturbance_intubation(start, duration, delta_values)

    def _validate_sequence_disturbances(self):
        """
        Validate the sequence of disturbances.
        :raises: ValueError if the sequence of disturbances is not valid.
        """

        starts = sorted(self.disturbances.keys())

        first_event = self.disturbances[starts[0]][0]

        if self.in_maintenance:
            if first_event == DisturbanceType.INTUBATION:
                raise ValueError("The intubation event should occur during the induction phase.")
        elif first_event != DisturbanceType.INTUBATION:
            raise ValueError("The first event must be an intubation event.")

        first_incision_time = float('inf')
        last_incision_time = float('-inf')
        incision_count = 0
        suture_count = 0
        intubation_count = 0

        for i in range(len(starts)):
            # Check that any skin manipulation occurs after at least the end of one incision
            # Check that suture events occur after incision events
            # At every suture must correspond a different incision
            start = starts[i]

            if start < 0 or start > self._t_sim:
                raise ValueError('Disturbances must start at a positive time and before the end of the simulation.')
            if self.disturbances[start][1] < 0:
                raise ValueError('Disturbances must have a positive duration.')
            if start + self.disturbances[start][1] > self._t_sim:
                raise ValueError('Disturbances must end before the end of the simulation.')
            if len(self.disturbances[start][2]) != 3:
                raise ValueError('Disturbances must have 3 delta values.')

            if i != len(starts) - 1:
                # if start + duration > start_next
                if start + self.disturbances[start][1] > starts[i + 1]:
                    raise ValueError('Overlapping disturbances detected.')

            if self.disturbances[start][0] == DisturbanceType.INCISION:
                incision_count += 1
                if first_incision_time == float('inf'):
                    first_incision_time = self.disturbances[start][1]
                last_incision_time = start + self.disturbances[start][1]
            elif self.disturbances[start][0] == DisturbanceType.SUTURE:
                if incision_count < suture_count or start < last_incision_time:
                    raise ValueError('Suture must occur after an incision event.')
                suture_count += 1
            elif self.disturbances[start][0] == DisturbanceType.SKIN_MANIPULATION:
                if start < first_incision_time:
                    raise ValueError('Skin manipulation must occur after an incision event.')
            elif self.disturbances[start][0] == DisturbanceType.INTUBATION:
                intubation_count += 1
                if intubation_count > 1:
                    raise ValueError('There must be exactly one intubation event.')

        if incision_count != suture_count:
            raise ValueError('There must be the same number of incision and suture events.')
