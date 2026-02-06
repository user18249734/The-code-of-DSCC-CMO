# TC-CMO
The code of the TC-CMO algorithm
# DSCC-CMOP
MATLAB source code for the **Dynamic Scaling Cooperative Co-evolutionary Algorithm for Constrained Multi-objective Problems (DSCC-CMOP)**, a decomposition-based multi-population evolutionary algorithm for constrained multi-objective optimization problems (CMOPs).

## Citation Note
If you use this code in your research, please cite the following paper:
> Yufeng Wang, Kaibin Du. A decomposition-based multi-population evolutionary algorithm for solving constrained multi-objective optimization problems. Expert Systems with Applications.

Your cooperation is highly appreciated!

## Algorithm Overview
This project provides MATLAB implementation of DSCC-CMOP, which is specifically designed for solving complex constrained multi-objective optimization problems. The algorithm integrates three key strategies to enhance optimization performance:
- Opposite search strategy
- Multi-population cooperation strategy
- Angle-based neighborhood search strategy

## File Descriptions
### TC-CMO.m
**Core Function**: Implements the core logic of the TC-CMO algorithm
**Core Features**:
1. Initialize three populations (P₁, P₂ and P₄)
2. Calculate the weight vectors (RW) and the neighborhood relation matrix (B)
3. Execute the evolutionary process in two phases

### ArchiveUpdate.m
**Core Function**: Screens and reserves non-dominated solutions to complete the update and maintenance of the archive set
**Core Features**:
1. Filter out feasible non-dominated solutions from the population
2. Eliminate duplicate solutions and perform normalization on the remaining solutions
3. Maintain population diversity based on crowding distance
4. Automatically eliminate inferior solutions when the archive set scale exceeds the threshold

### CalFitness.m
**Core Function**: Calculate the fitness value of each solution
**Core Features**:
1. Distinguish feasible solutions from infeasible solutions based on the constraint violation degree (CV)
2. Detect the dominance relationships among solutions and calculate the dominance count
3. Adopt the Euclidean distance to maintain the solution diversity and ensure uniform distribution
4. Output fitness values considering both solution quality and diversity

### EnvironmentalSelection.m
**Core Function**: Perform environmental selection to construct the next-generation population of the algorithm
**Core Features**:
1. Call the fitness calculation function to get individual fitness values
2. Rank individuals by fitness values and select the top-N solutions
3. Adopt a truncation mechanism to reduce population size when necessary

### EnvironmentalSelection_LAT.m
**Core Function**: Environmental selection for local auxiliary tasks (combines dynamic constraint boundaries to construct the next-generation population for local auxiliary tasks)
**Core Features**:
1. Calculate the total constraint violation degree of individuals and screen the candidate population according to the dynamic constraint boundary (VAR)
2. Call the fitness calculation function to solve individual fitness values by combining objective values and constraint violation degrees
3. Implement differentiated screening logic based on the relationship between the candidate population size and the target size N
4. Call the truncation function to eliminate individuals with high crowding degree for diversity maintenance when the population scale is excessive
5. Rank the final population by fitness values and screen to output the next-generation population with the specified scale

## Usage Instructions
1. Ensure the MATLAB environment is installed and the algorithm's required operation paths are configured correctly
2. Place all `.m` program files in the **same working directory**
3. Run the algorithm on the **PlatEMO platform**, and define the optimization problem to be solved through the platform's `Problem` class
4. Monitor the optimization results in real time, including:
   - Evolutionary process of the population
   - Final distribution of the Pareto front
