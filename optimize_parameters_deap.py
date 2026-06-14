import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from deap import base, creator, tools, algorithms
import random

# Existing Monod model code (monod_model function and constants) should be placed here

# Constants
S0 = 27.5  # g/L
X0 = 1.8  # g/L
P0 = 0.0  # g/L
mu_max = 1.2  # 1/hr
Ks = 1  # g/L
Yxs = 0.5
# ===========================================================================
# ENGINEERED-STRAIN ASSUMPTIONS — NOT experimentally measured.
# L-asparaginase was never produced at bench scale in this project. The values
# below represent a *hypothetical successfully metabolically-engineered* E. coli
# strain: literature-plausible targets for what such a strain could achieve.
# They are placeholders — anyone running this for real must replace them with
# measured values. They directly set product mass, growth, and therefore profit.
#   Ypx : g enzyme per g biomass  — strong recombinant overexpression (~30% of cell mass)
#   Yps : g enzyme per g glucose  — carbon-to-product efficiency (carbon balance)
#   Yax : g acetate per g biomass — reduced overflow metabolism (e.g. pta/ackA
#         knockout / BL21-type host). This is the key lever: at the original 0.92
#         the culture self-poisons with acetate near Ki_acetic_acid and biomass
#         caps at ~1.7 g/L; at 0.05 it reaches a realistic ~20 g/L (high cell density).
Ypx = 0.30  # ENGINEERED PLACEHOLDER (was 15.05, physically impossible)
Yps = 0.25  # ENGINEERED PLACEHOLDER
# ===========================================================================
Yxo2 = 1.06  # g O2/g X, assumed yield coefficient for oxygen consumption
Yax = 0.05  # ENGINEERED PLACEHOLDER — low-acetate strain (was 0.92; see note above)
O2_0 = .008
A0 = 0.0  # g/L
V0 = 3  # L
F = 170  # L/hr, assumed constant feed rate
Sf = 100  # g/L, assumed constant substrate concentration in feed
kLa = 250  # h^-1, assumed mass transfer coefficient 
O2_sat = 0.0075  # g/L, assumed saturation concentration of oxygen 
Vmax = 60000  # L, maximum reactor volume
endtime = 24  # hours
timestep = 0.5  # hours

# Inhibition constants
Ki_acetic_acid = 1.2  # g/L, assumed acetic acid inhibition constant
Ki_glucose = 200  # g/L, assumed glucose inhibition constant
K_o2 = .001 #g/L, assumed half saturation constant of oxygen

def monod_model(y, t, F, Sf, kLa):
    S, X, P, V, O2, A = y

    mu = mu_max * S / (Ks + S) * (1 - A / Ki_acetic_acid) * (1 - S / Ki_glucose) * O2/(K_o2+O2) # Monod equation with inhibition terms ADD OXYGEN
    mu = max(mu, 0.0)  # growth rate cannot be negative; the linear (Levenspiel) inhibition terms can drive mu < 0, which is unphysical
    q_o2 = Yxo2 * mu  # Oxygen consumption rate
    qp = Ypx * mu     # growth-associated specific product formation rate (g product / g biomass / hr)

    # Fed-batch volume balance: feed adds at rate F until the vessel is full,
    # then feeding stops. Using dV/dt = F (not F*(1 - V/Vmax)) keeps the volume
    # ODE consistent with the F/V dilution terms in the species balances below,
    # so total mass is conserved.
    F_eff = F if V < Vmax else 0.0
    dVdt = F_eff

    # Species mass balances (dilution rate = F_eff / V, consistent with dV/dt = F_eff)
    dSdt = -mu * X / Yxs - qp * X / Yps + F_eff * (Sf - S) / V  # substrate consumed for BOTH biomass and product (carbon balance)
    dXdt = mu * X - F_eff * X / V
    dPdt = qp * X - F_eff * P / V
    dO2dt = -q_o2 * X + kLa * (O2_sat - O2)
    dAdt = mu * X * Yax - F_eff * A / V

    return [dSdt, dXdt, dPdt, dVdt, dO2dt, dAdt]

# Initial conditions
y0 = [S0, X0, P0, V0, O2_0, A0]

# Time points
t = np.arange(0, endtime + timestep, timestep)

# Solve the ODE
result = odeint(monod_model, y0, t, args=(F, Sf, kLa), rtol=1e-6, atol=1e-6)

# Objective function to optimize
def objective_function(individual):
    S0_opt, X0_opt, V0_opt, Sf_opt, kLa_opt, F_opt = individual
    y0 = [S0_opt, X0_opt, P0, V0_opt, O2_0, A0]
    result = odeint(monod_model, y0, t, args=(F_opt, Sf_opt, kLa_opt), rtol=1e-6, atol=1e-6)
    P_result = result[:, 2]
    V_result = result[:, 3]

    # Calculate profit.
    # Concentrations are g/L and volumes are L, so every mass below is in GRAMS.
    # Prices are quoted per kg, so convert grams -> kg (divide by 1000) before
    # multiplying by a price, otherwise all dollar figures are 1000x too large.
    G_PER_KG = 1000.0
    final_product_mass = P_result[-1] * V_result[-1] / G_PER_KG       # kg product (final conc. x final broth volume)
    total_glucose_mass_fed = F_opt * Sf_opt * endtime / G_PER_KG      # kg glucose fed over the batch
    initial_glucose_mass = S0_opt * V0_opt / G_PER_KG                 # kg glucose in the initial charge
    inoculum_mass = X0_opt * V0_opt / G_PER_KG                        # kg seed biomass charged

    # Prices. L-asparaginase is a leukemia therapeutic, not a commodity enzyme:
    # formulated drug (e.g. Oncaspar) runs ~$1,000,000 per GRAM. We model crude,
    # unpurified product and do NOT include downstream purification cost, so we
    # use a deliberately conservative bulk price of $1,000,000 per KG (i.e. $1/mg)
    # — ~1000x below the formulated drug to leave room for purification losses.
    # (Original $127/kg undervalued a therapeutic enzyme by ~7 orders of magnitude.)
    price_product = 1_000_000  # $/kg, conservative bulk therapeutic-grade L-asparaginase
    price_glucose = 0.56       # $/kg, bulk glucose (realistic, unchanged)
    price_inoculum = 50        # $/kg seed biomass — realistic seed-culture cost
                               # (was $78,400/kg, which dominated everything and was not biologically meaningful)

    profit_product = final_product_mass * price_product               # $
    cost_glucose = (initial_glucose_mass + total_glucose_mass_fed) * price_glucose  # $
    cost_biomass = inoculum_mass * price_inoculum                     # $

    profit = profit_product - cost_glucose - cost_biomass
    return profit,

# Set up DEAP components
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("S0", random.uniform, 10, 50)
toolbox.register("X0", random.uniform, 0.1, 10) 
toolbox.register("V0", random.uniform, 0, 15000)
toolbox.register("Sf", random.uniform, 50, 150)
toolbox.register("kLa", random.uniform, 100, 400)
toolbox.register("F", random.uniform, 50, 500)
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.S0, toolbox.X0, toolbox.V0, toolbox.Sf, toolbox.kLa, toolbox.F), n=1)  

toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=10, indpb=0.1)
toolbox.register("select", tools.selBest)
toolbox.register("evaluate", objective_function)

def custom_mutate(individual, mu, sigma, indpb):
    for i in range(len(individual)):
        if random.random() < indpb:
            individual[i] += random.gauss(mu, sigma)
            if i == 0:  # S0
                individual[i] = max(min(individual[i], 50), 10)
            elif i == 1:  # X0
                individual[i] = max(min(individual[i], 10), 0.1)
            elif i == 2:  # V0
                individual[i] = max(min(individual[i], 15000), 0)
            elif i == 3:  # Sf
                individual[i] = max(min(individual[i], 150), 50)
            elif i == 4:  # kLa
                individual[i] = max(min(individual[i], 400), 100)
            elif i == 5:  # F
                individual[i] = max(min(individual[i], 500), 50)
    return individual,

toolbox.register("mutate", custom_mutate, mu=0, sigma=5, indpb=0.1)


def optimize_model(ngen, population_size, cxpb, mutpb):
    pop = toolbox.population(n=population_size)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=cxpb, mutpb=mutpb, ngen=ngen, stats=stats, halloffame=hof, verbose=True)

    return hof[0]

optimized_params = optimize_model(ngen=50, population_size=100, cxpb=0.8, mutpb=0.2)
S0_opt, X0_opt, V0_opt, Sf_opt, kLa_opt, F_opt = optimized_params
total_glucose_mass_fed = F_opt * Sf_opt * endtime
print("Optimized parameters:")
print(f"Initial glucose concentration (S0): {S0_opt:.2f} g/L")
print(f"Initial biomass concentration (X0): {X0_opt:.2f} g/L")  # Print optimized biomass concentration
print(f"Initial volume (V0): {V0_opt:.2f} L")
print(f"Glucose concentration in the feed (Sf): {Sf_opt:.2f} g/L")
print(f"Mass transfer coefficient (kLa): {kLa_opt:.2f} h^-1")
print(f"Optimized feed rate (F): {F_opt:.2f} L/hr")

