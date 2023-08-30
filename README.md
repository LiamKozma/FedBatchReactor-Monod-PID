# Fed-Batch Bioreactor Model

This repository contains a Python implementation of a Fed-batch bioreactor model that takes into account mass and energy balance equations along with PID controllers for optimized process control.

## Table of Contents
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
- [Files in this Repository](#files-in-this-repository)
- [Usage](#usage)
- [Optimization](#optimization)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3.x
- NumPy
- Matplotlib
- SciPy
- Seaborn
- pyswarm (for PID parameter optimization)
- DEAP (for parameter optimization)

You can install the required packages using pip:

```bash
pip install numpy scipy matplotlib seaborn pyswarm deap
```

## Files in this Repository

- `bioreactor_model.py`: The main bioreactor model using Monod kinetics and PID control.
- `optimize_pid_pyswarm.py`: File that uses pyswarm to optimize PID parameters.
- `optimize_parameters_deap.py`: File that uses DEAP to optimize initial concentrations and other parameters.

## Usage

Run the main bioreactor model:

```bash
python bioreactor_model.py
```
Note: Before running this script, you should update the model parameters (e.g., PID gains, initial concentrations) based on the optimized values obtained from running optimize_pid_pyswarm.py and optimize_parameters_deap.py. These values must be manually inserted into the corresponding variables in bioreactor_model.py.

### Optimize PID parameters using pyswarm:

```bash
python optimize_pid_pyswarm.py
```

After running this script, you will obtain optimized PID parameters. Manually update these optimized values in the main bioreactor model (bioreactor_model.py) before running it.

### Optimize initial conditions and parameters using DEAP:

```bash
python optimize_parameters_deap.py
```

After running this script, you will obtain optimized initial conditions and other parameters. Manually update these optimized values in the main bioreactor model (bioreactor_model.py) before running it.

## Optimization

Two types of optimization are performed:

1. PID parameter optimization for controlling the heat exchanger area, done using pyswarm.
2. Optimization of initial conditions and reactor parameters like initial glucose concentration, initial biomass concentration, initial volume, glucose concentration in the feed, mass transfer coefficient, and feed rate using DEAP.

## Contributing

If you'd like to contribute, please fork the repository and use a feature branch. Pull requests are warmly welcome.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

- Thanks to anyone whose code was used as inspiration.
