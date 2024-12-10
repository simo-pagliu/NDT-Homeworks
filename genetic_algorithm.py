import random
import numpy as np
import matplotlib.pyplot as plt
import json
import time
import os
import struct
import pandas as pd  # To store results in a DataFrame for easy comparison
import concurrent.futures
import pandas as pd
import loop
from scipy.interpolate import interp1d

# Set working directory to the location of this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

################################################################################
# FITNESS FUNCTION
# Define the Styblinski-Tang function to minimize
# Set the domain limits to -5 to 5 for both x and y
# The global minimum is at f(-2.903534, -2.903534) = -78.332
def fitness_func(thickness_cladding, plenum_height):
    try:
        T_map, Geometrical_Data, He_percentage, Plenum_Pressure, Coolant_Velocity, Void_Swelling, params, Burnup, coolant_infinity_limit = loop.main(thickness_cladding, plenum_height)

        # Fuel Max Temperature
        Max_Fuel_Temperature = np.max(T_map.T)

        # Cladding Max Temperature
        cladding_midline_radii = [
            (Geometrical_Data.cladding_outer_diameter[i] / 2 +
            (Geometrical_Data.cladding_outer_diameter[i] / 2 - Geometrical_Data.thickness_cladding[i])) / 2
            for i in range(len(Geometrical_Data.h_values))
        ]
        for i, radius in enumerate(cladding_midline_radii):
            # Interpolate temperature at the midline radius
            interp_func = interp1d(T_map.r[i, :], T_map.T[i, :], kind="linear", fill_value="extrapolate")
            temp_at_midline = interp_func(radius)
        Max_Cladding = np.max(temp_at_midline)

        # Plenum Pressure
        Plenum_Pressure

        # Istantaneous cladding plastic strain
        #
        #

        # Maximum cladding volumetric swelling
        Maximum_Cladding_Volumetric_Swelling = np.max(Void_Swelling)

        # Maximum coolant velocity
        Maximum_Coolant_Velocity = np.max(Coolant_Velocity)

        # Gap thickness
        Min_Gap_Thickness = np.min(np.array(Geometrical_Data.cladding_outer_diameter) / 2 - np.array(Geometrical_Data.fuel_outer_diameter) / 2 - np.array(Geometrical_Data.thickness_cladding))

        Fuel_limit = 2600 + 273.15
        Cladding_limit = 650 + 273.15

        Fitness_Function = Max_Fuel_Temperature

        fitness_handicap = 1e3
        fitness_handicap_if_fails = 1e5

        Maximum_Cladding_Volumetric_Swelling_limit = 3
        if Maximum_Cladding_Volumetric_Swelling > Maximum_Cladding_Volumetric_Swelling_limit:
            Fitness_Function += fitness_handicap

        Maximum_Coolant_Velocity_limit = 8
        if Maximum_Coolant_Velocity > Maximum_Coolant_Velocity_limit:
            Fitness_Function += fitness_handicap

        Plenum_Pressure_limit = 5e6
        if Plenum_Pressure > Plenum_Pressure_limit:
            Fitness_Function += fitness_handicap
    except:
        Fitness_Function = fitness_handicap_if_fails
    return Fitness_Function

limit_x_low = 20e-6 # Cladding Thickness lower limit
limit_x = 500e-6 # Cladding Thickness upper limit
limit_y = 1.2 # Plenum Height upper limit
################################################################################

################################################################################
# Probability based on fitness (for minimization)
# The probability is inversely proportional to the fitness value
# The higher the fitness, the lower the probability
def fitness_probability(set):
    fitness_values = [fitness_func(pair[0], pair[1]) for pair in set]
    if max(fitness_values) == min(fitness_values):
        probabilities = [1 / len(set)] * len(set)
    else:
        # max_fitness = max(fitness_values)
        # inverted_fitness = [(max_fitness - fitness) for fitness in fitness_values]
        # total_inverted_fitness = sum(inverted_fitness)
        # probabilities = [fit / total_inverted_fitness for fit in inverted_fitness]
        fitness_values = - np.array(fitness_values)
        # set any negative fitness to 0
        fitness_values = np.where(fitness_values < 0, 0, fitness_values)
        # calculate the sum of the fitness values
        sum_fitness = sum(fitness_values)
        # calculate the probability of each individual
        probabilities = fitness_values / sum_fitness
    return probabilities
################################################################################

################################################################################
# Initialize Population
# Randomly initialize the population within the given domain limits
def initialize_population(population_size, x_min, x_max, y_min, y_max):
    return [[random.uniform(x_min, x_max), random.uniform(y_min, y_max)] for _ in range(population_size)]
################################################################################

################################################################################
# Selection
# Tournament selection is basically doing nothing as it keeps the # the same
# Roulette selection is reducing the population size
def tournament_selection(population, tournament_size):
    """Selects the best individual from a random subset of the population.
    Chooses the one with the highest fitness value."""
    selected_parents = []
    for _ in range(len(population)):
        tournament = random.sample(population, tournament_size)
        selected = min(tournament, key=lambda x: fitness_func(x[0], x[1]))
        selected_parents.append(selected)
    return selected_parents

def roulette_selection(population, tournament_size):
    """Selects the best individual from a random subset of the population.
    Based on a probability related to the fitness value."""
    selected_parents = []
    for _ in range(len(population)):
        tournament = random.sample(population, tournament_size)
        probabilities = fitness_probability(tournament)
        selected = random.choices(tournament, weights=probabilities, k=1)[0]
        selected_parents.append(selected)
    return selected_parents
################################################################################

################################################################################
# Crossover
# Right now I'm replacing the parents with the children
# Is it correct? Should I append the children to the population?
def punnet_crossover(population, probability):
    """ Analogy to the Punnet Square.
    Performs crossover between two parents.
    Combines the x and y values of the parents to create two new children.
    Chosen by probability weighted on fitness."""
    crossed_over = population.copy()
    for _ in range(len(population) // 2):
        if random.random() < probability:
            indexes = random.sample(range(len(population)), 2)
            parent1, parent2 = population[indexes[0]], population[indexes[1]]
            children = [[parent1[0], parent2[1]], [parent2[0], parent1[1]], parent1, parent2]
            probabilities = fitness_probability(children)
            selected = random.choices(children, weights=probabilities, k=2)
            crossed_over[indexes[0]], crossed_over[indexes[1]] = selected
    return crossed_over

def bit_crossover(population, probability):
    """Crosses over a portion of the bits of the binary representation of the two parents."""
    crossed_over = population.copy()

    for _ in range(len(population) // 2):
        if random.random() < probability:
            # Select two random parents
            indexes = random.sample(range(len(population)), 2)
            parent1, parent2 = population[indexes[0]], population[indexes[1]]

            # Initialize new offspring
            offspring1, offspring2 = [], []

            for i in range(2):  # Loop over the two variables of the parents
                # Convert the current variable of the parents to binary
                binary1 = f"{struct.unpack('!I', struct.pack('!f', parent1[i]))[0]:032b}"
                binary2 = f"{struct.unpack('!I', struct.pack('!f', parent2[i]))[0]:032b}"

                # Randomly choose crossover points
                # Only the first bunch of bits are used --> just change the fraction part
                start = random.randint(0, 19)
                end = random.randint(start + 1, 20)

                # Perform crossover
                new_binary1 = binary1[:start] + binary2[start:end] + binary1[end:]
                new_binary2 = binary2[:start] + binary1[start:end] + binary2[end:]

                # Convert back to floats
                new_float1 = struct.unpack('!f', struct.pack('!I', int(new_binary1, 2)))[0]
                new_float2 = struct.unpack('!f', struct.pack('!I', int(new_binary2, 2)))[0]

                # Add to offspring
                offspring1.append(new_float1)
                offspring2.append(new_float2)

            # Update the population
            crossed_over[indexes[0]] = offspring1
            crossed_over[indexes[1]] = offspring2

    return crossed_over
################################################################################

################################################################################
# Mutation
# Right now I'm replacing the parents with the children
# Is it correct? Should I append the children to the population?
def redefine_mutation(population, probability):
    """ Randomly redifines the x and y values, indipendently."""
    mutated = population.copy()
    for i in range(len(population)):
        # Mutate x
        if random.random() < probability:
            mutated[i] = [random.uniform(limit_x_low, limit_x), population[i][1]]
        # Mutate y
        if random.random() < probability:
            mutated[i] = [population[i][0], random.uniform(0, limit_y)]
    return mutated

def bit_mutation(population, probability):
    """
    Randomly flips a bit in the binary representation of the decimal part
    of the gene. Ensures the mutation only affects decimals and numbers
    remain as units.
    """
    mutated = [row.copy() for row in population]  # Create a copy of the population
    for i in range(len(population)):
        if random.random() < probability:
            for j in range(len(population[i])):
                value = population[i][j]

                # Convert float to binary representation as an integer
                packed_value = struct.unpack('!I', struct.pack('!f', value))[0]
                
                # The mantissa occupies bits 0-22 (23 bits total)
                # Only mutate bits within the mantissa
                mantissa_bit_index = random.randint(0, 22)

                # Flip the chosen mantissa bit
                flipped_value = packed_value ^ (1 << mantissa_bit_index)

                # Convert back to float
                mutated_value = struct.unpack('!f', struct.pack('!I', flipped_value))[0]
                mutated[i][j] = mutated_value
    return mutated

################################################################################

################################################################################
# Save and Load State
def save_state(population, avg_fitness_history, filename="ga_state.json"):
    state = {
        "population": population,
        "avg_fitness_history": avg_fitness_history
    }
    with open(filename, "w") as f:
        json.dump(state, f)
    print("State saved to", filename)

def load_state(filename=f"ga_state.json"):
    with open(filename, "r") as f:
        state = json.load(f)
    print("State loaded from", filename)
    return state["population"], state["avg_fitness_history"]
################################################################################

################################################################################
# Genetic Algorithm
def GA_loop(max_iterations, initial_population_size, probability_of_crossover,
              probability_of_mutation, tournament_fraction,
              selection, crossover, mutation,
              plot_enabled, load_previous, save_enabled, print_result, log):
    import time
    # Initialization
    if load_previous:
        population, avg_fitness_history = load_state()
    else:
        population = initialize_population(initial_population_size, limit_x_low, limit_x, 0, limit_y)
        avg_fitness_history = []

    # Initialize plot if enabled
    if plot_enabled:
        plt.ion()
        fig, ax = plt.subplots(figsize=(10, 6))
        line, = ax.plot([], [], label="Average Fitness", color="orange")
        ax.set_xlim(0, max_iterations)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Average Fitness")
        ax.set_title("Convergence of Average Fitness Over Iterations")
        ax.legend()
            
    # Main loop
    converged = False
    time_start = time.time()
    for iteration in range(max_iterations):
        tournament_size = max(2, int(tournament_fraction * len(population)))
        population = selection(population, tournament_size)
        population = crossover(population, probability_of_crossover)
        population = mutation(population, probability_of_mutation)

        # Check for any NaN values and remove them
        population = [pair for pair in population if not any(np.isnan(gene) for gene in pair)]

        # Track average fitness
        fitness_values = [fitness_func(pair[0], pair[1]) for pair in population]
        avg_fitness = sum(fitness_values) / len(fitness_values)
        avg_fitness_history.append(avg_fitness)
        
        # Update plot if enabled
        if plot_enabled:
            line.set_xdata(range(len(avg_fitness_history)))
            line.set_ydata(avg_fitness_history)
            ax.set_xlim(0, max(len(avg_fitness_history), max_iterations))
            fig.canvas.draw()  # Redraw the figure
            fig.canvas.flush_events()  # Process the events
            plt.pause(0.01)

        # Convergence check
        if len(avg_fitness_history) > 1 and abs(avg_fitness_history[-1] - avg_fitness_history[-2]) < 0.001 or len(population) == 1:
            print(f"\033[92mConverged after {iteration + 1} iterations.\033[0m")
            convergence_string = f"CONVERGED @ {iteration + 1}"
            converged = True
            break
    
    if not converged:
        residuals = []
        for i in range(len(avg_fitness_history)):
            residuals.append(abs(avg_fitness_history[i] - avg_fitness_history[i-1]))
        residuals = np.array(residuals)
        residual_threshold = 10
        mask =  residuals < residual_threshold
        residual = sum(residuals[mask]) / (len(mask))
        print(f"\033[93mStopped after reaching the maximum number of iterations without reaching convergence, avg residual is {residual}.\033[0m")
        mask = residuals > residual_threshold
        to_converge = sum(mask) / len(mask)
        convergence_string = f"MAXED @ {residual} - Fraction of Iterations to converge: {to_converge}"

    # Get final result
    best_solution = min(population, key=lambda x: fitness_func(x[0], x[1]))
    best_fitness = fitness_func(best_solution[0], best_solution[1])
    
    exe_time = time.time() - time_start
    # Closing plot if enabled
    if plot_enabled:
        plt.ioff()
        plt.show()

    # Save state if enabled
    if save_enabled:
        save_state(population, avg_fitness_history)
    
    if print_result:
        print(f"Best solution: {best_solution}")
        print(f"Best fitness: {best_fitness}")

    # Logging best result if enabled
    if log:
        log_file = "genetic_algotithm.log"
        with open(log_file, "a") as f:
            selections = "Tournament" if selection == tournament_selection else "Roulette"
            crossovers = "Punnet" if crossover == punnet_crossover else "Bit"
            mutations = "Redefine" if mutation == redefine_mutation else "Bit"
            method = f"{selections} - {crossovers} - {mutations}"
            stats = f"I: {max_iterations}, Pop: {initial_population_size}, P_CS: {probability_of_crossover}, P_M: {probability_of_mutation}, P_S: {tournament_fraction}" 
            results = f"Result: {best_fitness} @ x: {best_solution[0]} y: {best_solution[1]}"
            time = f"Time: {exe_time}"
            f.write(f"{method}, {stats}, {results}, {convergence_string}, {time}\n")
        print(f"Results logged in {log_file}")

    # Return performance metrics
    return {
        "iterations": iteration + 1,
        "final_fitness": best_fitness,
        "best_solution": best_solution,
    }

# Define a wrapper function for GA_loop
def run_ga(params):
    population, max_iter, crossover_prob, mutation_prob, tournament_frac, selection, crossover, mutation = params

    # Start execution time tracking
    start_time = time.time()

    # Run the Genetic Algorithm
    results = GA_loop(
        max_iter, len(population), crossover_prob, mutation_prob, tournament_frac,
        selection, crossover, mutation,
        plot_enabled=False, load_previous=False, save_enabled=False, print_result=False, log=False
    )
    
    # Calculate execution time
    execution_time = time.time() - start_time

    # Prepare detailed results similar to the log format
    return {
        "method": f"{selection.__name__} - {crossover.__name__} - {mutation.__name__}",
        "population_size": len(population),
        "max_iterations": max_iter,
        "crossover_probability": crossover_prob,
        "mutation_probability": mutation_prob,
        "tournament_fraction": tournament_frac,
        "final_iterations": results["iterations"],
        "final_fitness": results["final_fitness"],
        "execution_time": execution_time,
    }

################################################################################

################################################################################
# Main Function
if __name__ == "__main__": 
    loop_iterations = 1
    initial_population_size = 150  # Initialize a population
    max_iterations = 1
    crossover_prob = 0.6
    mutation_prob = 0.6
    tournament_frac = 0.2
    selection = tournament_selection  # Replace with your selection function # Options: tournament_selection, roulette_selection
    crossover = bit_crossover  # Replace with your crossover function # Options: punnet_crossover, bit_crossover
    mutation = redefine_mutation  # Replace with your mutation function # Options: redefine_mutation, bit_mutation

    # Run GA with or without population division
    for i in range(loop_iterations):
        print(f"Iter: {i+1}/{loop_iterations}")
        tic = time.time()
        final_population = GA_loop(
            max_iterations, initial_population_size, crossover_prob, mutation_prob, tournament_frac,
            selection, crossover, mutation,
            plot_enabled=False, load_previous=False, save_enabled=False, print_result=False, log = True
        )  
        toc = time.time()
        print(f"Elapsed time: {toc - tic} seconds")
################################################################################