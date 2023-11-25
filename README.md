# Tuned Mass Damper (TMD) System Simulation
This Python code simulates the behavior of a Tuned Mass Damper (TMD) system using numerical methods. The TMD system consists of two masses connected by springs and dampers. The main purpose of the TMD is to mitigate the effects of vibrations.

## Prerequisites
Make sure you have the following Python libraries installed:

- 'numpy'
- 'sympy'
- 'scipy'
- 'matplotlib'

## Usage
1. Import necessary libraries:
2. Define symbolic parameters and equations of motion
3. Define functions for the system dynamics
4. Set up the initial conditions and simulation parameters
5. Integrate the system equations using 'odeint'
6. Visualize the results with an animation and save it as a GIF

## Parameters
- 'm': Mass of the main structure
- 'w': Natural frequency of the main structure
- 'k': Stiffness of the spring connecting the main structure and the TMD mass
- 'l1': Equilibrium position of the main structure
- 'xi': Damping ratio of the main structure
- 'm_d': Mass of the tuned mass damper (TMD)
- 'w_d': Natural frequency of the TMD
- 'k_d': Stiffness of the spring connecting the TMD and the main structure
- 'l2': Equilibrium position of the TMD
- 'xi_d': Damping ratio of the TMD
- 'f0': Amplitude of the external force
- 'bw': Frequency of the external force

## Results
The simulation provides the displacement of the main structure and the TMD over time. The animation shows the movement of the masses and the springs in the TMD system.
