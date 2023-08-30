import numpy as np
from pyswarm import pso
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# Import rcParams from matplotlib
from matplotlib import rcParams
import seaborn as sns
import math

# Constants
S0 = 39.40  # g/L
X0 = 0.10  # g/L
P0 = 0.0  # g/L
mu_max = 1.2  # 1/hr
Ks = 1  # g/L
Yxs = 0.5
Ypx = 15.05  # find the correct value later
Yxo2 = 1.06  # g O2/g X, assumed yield coefficient for oxygen consumption
Yax = 0.92  # g acetic acid/g X, assumed value, find the correct value later
O2_0 = .008
A0 = 0.0  # g/L
V0 = 79.21  # L
F = 201.34  # L/hr, assumed constant feed rate
T_feed = 298.15  # Feed temperature (K), assuming room temperature (25°C)
Sf = 55.73  # g/L, assumed constant substrate concentration in feed
kLa = 258.13  # h^-1, assumed mass transfer coefficient 
O2_sat = 0.0075  # g/L, assumed saturation concentration of oxygen 
Vmax = 60000  # L, maximum reactor volume
aspect_ratio = 3
endtime = 24  # hours
timestep = 0.5  # hours

# Inhibition constants
Ki_acetic_acid = 1.2  # g/L, assumed acetic acid inhibition constant
Ki_glucose = 200  # g/L, assumed glucose inhibition constant
K_o2 = .001 #g/L, assumed half saturation constant of oxygen

# Energy balance constants
heat_of_fermentation_ecoli = -20000  # Heat of fermentation for E. coli (J/g)
heat_transfer_coefficient = 500  # Heat transfer coefficient (W/m^2K)
reactor_area = 400  # Reactor surface area touching heat exchanger (m^2)
ambient_temperature = 300  # Ambient temperature (K)
specific_heat_capacity = 4.186e3  # J/kg·K
rho = 1 #kg/L
T0 = 298.15  # Initial temperature value (K), assuming room temperature (25°C)
initial_internal_energy = V0 * rho * specific_heat_capacity * T0  # J
agitation_power_per_volume = 10  # W/L



# Heat exchanger area
volume_m3 = Vmax / 1000
D = (6 * volume_m3 / (math.pi * aspect_ratio))**(1/3)
initial_H = V0 / (math.pi * (D/2)**2)  # Initial height of the liquid in the bioreactor (m)
H_max = aspect_ratio * D  # Maximum height of the bioreactor (m)
H_exchanger_max = H_max  # Maximum allowed height for the heat exchanger (m)



# Define an objective function that takes the PID parameters as input
def objective_function(pid_params):
    Kp, Ki, Kd = pid_params
    
    # Set up and solve the ODE with the given PID parameters
    # (insert your existing code to simulate the bioreactor with the given PID parameters)

    def heat_exchanger_area(V):
        V_m3 = V / 1000  # Convert the volume from liters to cubic meters
        H = V_m3 / (math.pi * (D/2)**2)  # Calculate the height of the liquid in the bioreactor based on its volume
        return 2 * math.pi * (D/2) * H  # Lateral surface area of the cylindrical bioreactor
  
    def pid_controller(Kp, Ki, Kd, setpoint, current_value, integral, derivative, dt):
        error = setpoint - current_value
        P_term = Kp * error
        I_term = Ki * integral * dt
        D_term = Kd * derivative / dt
        return P_term + I_term + D_term
    
    temperature_errors = []
    
    
    def monod_model(t, y, F, Sf):
        S, X, P, V, O2, A, U, T, int_error = y
    
    
        mu = mu_max * S / (Ks + S) * (1 - A / Ki_acetic_acid) * (1 - S / Ki_glucose) * O2/(K_o2+O2) # Monod equation with inhibition terms ADD OXYGEN 
        q_o2 = Yxo2 * mu  # Oxygen consumption rate
    
        # Mass balance equations
        dSdt = -mu * X / Yxs + F * (Sf - S) / V
        dXdt = mu * X - F * X / V
        dPdt = mu * X * Ypx - F * P / V
        dVdt = F * (1 - V / Vmax)
        dO2dt = -q_o2 * X + kLa * (O2_sat - O2)
        dAdt = mu * X * Yax - F * A / V
    
        # PID controller for the heat exchanger area
        setpoint_temperature = 328.5
        Kp = 10000000  # Proportional gain
        Ki = 0 # Integral gain
        Kd = 0  # Derivative gain
        
        temperature_errors.append(setpoint_temperature - T)
        integral_error = int_error + (setpoint_temperature - T) * timestep
        
        derivative_error = (setpoint_temperature - T) / timestep
        heat_exchanger_area_correction = pid_controller(Kp, Ki, Kd, setpoint_temperature, T, integral_error, derivative_error, timestep)
        current_reactor_area = heat_exchanger_area(V) + heat_exchanger_area_correction
        
        
        
        # Energy balance equation
        dUdt = -X * mu * heat_of_fermentation_ecoli \
           - heat_transfer_coefficient * current_reactor_area * (T - ambient_temperature) \
           + F * rho * specific_heat_capacity * (T_feed - T) \
           + agitation_power_per_volume * V
    
        # Temperature balance equation
        dTdt = dUdt / (V * specific_heat_capacity * rho)
    
        return np.array([dSdt, dXdt, dPdt, dVdt, dO2dt, dAdt, dUdt, dTdt, setpoint_temperature - T])
    
    
    
    # Initial conditions
    y0 = [S0, X0, P0, V0, O2_0, A0, initial_internal_energy, T0, 0]  # Add initial_internal_energy and initial_temperature_error here
     
    
    # Time points
    t = np.arange(0, endtime + timestep, timestep)
    t_array = t  # Store the time points array for use inside the monod_model function
    
    
    
    # Solve the ODE
    from scipy.integrate import solve_ivp
    result = solve_ivp(lambda t, y: monod_model(t, y, F, Sf), (0, endtime), y0, t_eval=t, method='BDF')
    
    # We'll minimize the mean squared error between the desired setpoint and the actual temperature
    setpoint_temperature = 328.5
    T_result = result.y[7]  # Extract temperature values (index 7 corresponds to dTdt)
    mse = np.mean((T_result - setpoint_temperature)**2)

    
    return mse

# Set bounds for the PID parameters
lb = [0, 0, 0]  # Lower bounds for Kp, Ki, and Kd
ub = [1000, 1000, 1000]  # Upper bounds for Kp, Ki, and Kd

# Optimize the PID parameters using the particle swarm optimization algorithm
best_pid_params, _ = pso(objective_function, lb, ub)

# Extract the optimized PID parameters
Kp_optimized, Ki_optimized, Kd_optimized = best_pid_params

print("Optimized Kp:", Kp_optimized)
print("Optimized Ki:", Ki_optimized)
print("Optimized Kd:", Kd_optimized)
