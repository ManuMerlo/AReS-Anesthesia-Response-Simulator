class PID:
    def __init__(self, kp, ki=0, kd=0, setpoint=0.5, sample_time=5, output_limits=(0, 2)):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.setpoint = setpoint
        self.sample_time = sample_time
        self.output_limits = output_limits

        self.integral = 0
        self.last_error = 0
        self.first_computation = True
        self.last_output = 0

        # Compute max integral based on infusion rate limit and ki to prevent overshoot
        if ki != 0:
            self.max_integral_error = (self.output_limits[1] / ki)
        else:
            self.max_integral_error = 0

    def compute(self, process_variable, delta_time):
        if delta_time < self.sample_time:
            return self.last_output  # Return last output

        error = process_variable - self.setpoint

        # Integral term
        self.integral += error * delta_time

        # Clamp the integral term based on the calculated max integral error ( anti-windup )
        self.integral = max(min(self.integral, self.max_integral_error), -self.max_integral_error)

        # Derivative term
        # Avoid derivative kick on the first iteration
        if self.first_computation:
            self.first_computation = False
            self.last_error = error
            derivative = 0
        else:
            derivative = (error - self.last_error) / delta_time

        # Compute PID output = P + I + D
        control_output = self.kp * error + self.ki * self.integral + self.kd * derivative

        # Clamp the output to the defined limits
        control_output = max(self.output_limits[0], min(control_output, self.output_limits[1]))

        # Save the current state for the next iteration
        self.last_error = error
        self.last_output = control_output

        return control_output
