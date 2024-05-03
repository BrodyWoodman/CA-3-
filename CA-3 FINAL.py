# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 15:48:35 2024

@author: brody
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
u = 1.0  # Elastic constant
a = 0.1  # Quartic coefficient
beta = 1.0  # Inverse temperature
num_steps = 10000  # Number of Monte Carlo steps
num_bins = 50  # Number of bins for histogram

#POTENTIAL FUNCTIONS! --------------------------------

# Function to calculate potential energy for harmonic oscillator
def harmonic_potential(x, u):
    return u**2 * x**2

# Function to calculate potential energy for anharmonic oscillator
def anharmonic_potential(x, u, a):
    return u**2 * x**2 + a * x**4

# Metropols algorithm!
def metropolis(x_old, u, a, beta):
    # Generate a trial state
    x_trial = x_old + np.random.normal(0, 1)
    
    # Calculate potential energy differenc (Markov chain)
    delta_V = (anharmonic_potential(x_trial, u, a) - anharmonic_potential(x_old, u, a))
    
    # Acceptance probability
    if delta_V < 0 or np.random.uniform(0, 1) < np.exp(-beta* delta_V):
        return x_trial
    else:
        return x_old

# MONTE CARLO FOR HARMONIC.OSCILLATOR
def simulate_harmonic_oscillator(u, beta, num_steps, num_bins):
    # Initialize position and energy lists
    positions = []
    energies = []
    
    # Start position
    x = 0
    
    # Perform Monte Carlo steps
    for _ in range(num_steps):
        # Update position using Metropolis algorithm
        x = metropolis(x, u, 0, beta)
        positions.append(x)
        
        # Calculate energy
        energy = harmonic_potential(x, u)
        energies.append(energy)
#ENERGYS! 
    # ground state energy
    ground_state_energy = np.mean(energies)
    
    # first excited state energy
    first_excited_energy = min(energy for energy in energies if energy > ground_state_energy) #this is kind a bullshit function that i dont think works
    
    # Probability distribution from the ground state energis
    plt.hist(positions, bins=num_bins, density=True)
    plt.title("Ground State Probability Distribution (Harmonic Oscillator)")
    plt.xlabel("Position")
    plt.ylabel("Probability Density")
    plt.show()
    
    return ground_state_energy, first_excited_energy

#Monte Carlo for Anharmonic.oscillator!
def simulate_anharmonic_oscillator(u, a, beta, num_steps, num_bins):
    # Initialize position and energy lists
    positions = []
    energies = []
    
    # Initialize position
    x = 0
    
    # Monte carlo steps!
    for z in range(num_steps):
        # Update position using Metropolis algorithm
        x = metropolis(x, u, a, beta)
        positions.append(x)
        
        # Energy stuff again
        energy = anharmonic_potential(x, u, a)
        energies.append(energy)
    
    # ground state 
    ground_state_energy = np.mean(energies)
    
    #first excited state energy
    first_excited_energy = min(energy for energy in energies if energy > ground_state_energy)
    
    # probability distribution for the Anhamronic oscialltor
    plt.hist(positions, bins=num_bins, density=True)
    plt.title("Ground State Probability Distribution (Anharmonic Oscillator)")
    plt.xlabel("Position")
    plt.ylabel("Probability Density")
    plt.show()
    
    return ground_state_energy, first_excited_energy

#SIMULATIONS!

# Simulate harmonic oscillator
ground_state_energy_harmonic, first_excited_energy_harmonic = simulate_harmonic_oscillator(u, beta, num_steps, num_bins)
print("Ground State Energy (Harmonic Oscillator):", ground_state_energy_harmonic, "eV")
print("First Excited State Energy (Harmonic Oscillator):", first_excited_energy_harmonic, "eV")
# Simulate ANharmonic oscillator
ground_state_energy_anharmonic, first_excited_energy_anharmonic = simulate_anharmonic_oscillator(u, a, beta, num_steps, num_bins)
print("Ground State Energy (Anharmonic Oscillator):", ground_state_energy_anharmonic,"eV")
print("First Excited State Energy (Anharmonic Oscillator):", first_excited_energy_anharmonic,"eV")