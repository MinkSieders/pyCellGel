# pyCellGel

## Description
Script that can model and solve _in silico_ cell-laden hydrogels based on user input.

## Main Commands
Solver (solver): Computes the required loading concentrations to achieve a specified inter-CFU distribution within the 
hydrogel.

Modeller (modeller): Models hydrogels in different user-defined shapes and concentrations.

Figure Generation (figure): Generates visualization figures based on the model data.

Growth Modeller (model_growth): Models the growths of microbial CFU's randomly dispersed through a cell-laden hydrogel. Assumes the gel restricts mobility or that microbes form spherical microcolonies. 

## Input
Each command is accessed through a command-line argument specifying the mode, followed by additional flags for 
customization. Most commands are fitted with default values, make sure to change according to your needs. 

For in-depth explanation on all the flags within each command, use --help or -h. 

## Usage | Solver

`pyCellGel solver`

Calculates the required bacterial concentration to achieve a given target inter-CFU distance within the hydrogel. This is useful for ensuring uniform distribution in experimental setups.

### Arguments:

`--P_target` (float, required): Population fraction required to satisfy the target distance.

`--r_target` (float, required, μm): Desired inter-CFU distance (converted internally to mm).

`--model_prediction` (bool, optional, default: False): If True, a cilindrical hydrogel is additionally modeled based on the computed parameters.

### Example Usage:

To model a gel in which 95 % of the population satisfies a mean nearest neighbour CFU distance of 50 microns, run the following command:

`python pyCellGel.py solver --P_target 0.95 --r_target 50 --model_prediction True`

## Modeller (modeller)

The modeller command allows users to define hydrogel shapes and dimensions while specifying bacterial concentrations either directly or via an optical density (OD600) measurement.

### Shape Options:

Cylinder

Sphere

CubeSquare

CubeRect

### Arguments:

`--shape` (str, optional, default: Sphere): Shape to be modeled. Can be one of the shape options or All.

`--height` (float, optional, mm, default: 0.5): Height of the hydrogel (for applicable shapes).

`--radius` (float, optional, mm, default: 1.0): Radius of the hydrogel (for spheres and cylinders).

`--width` (float, optional, mm, default: 0.75): Width of the hydrogel (for cuboid shapes).

`--length` (float, optional, mm, default: 1.0): Length of the hydrogel (for cuboid shapes).

`--concentration` (float, optional, CFU/ml): Manually specify bacterial concentration. If provided, OD600 conversion is skipped.

`--optical_density` (float, optional, default: 0.005): Specify concentration via OD600 (requires OD-to-CFU conversion curve; defaults to E. coli if not provided).

`--optical_density_curve` (str, optional): Path to a .tsv file containing OD vs. CFU/ml data. If omitted, the script defaults to an internal E. coli calibration curve.

### Example Usage:

`python pyCellGel.py modeller --shape CubeRect --height 1.0 --width 0.5 --length 2.0 --optical_density 0.01`

## Figure Generation (figure)

The figure command generates visualization figures based on the computed or modeled hydrogel data.

### Arguments:

(No additional arguments required.)

### Example Usage:

`python pyCellGel.py figure`

### Output Handling

For each command, an output directory (output_<command>) is created to store results. If an existing directory is found, it is removed and recreated to ensure clean output generation.

## Growth Modeller (model_growth)

The growth modeller command simulates the growth of microbial CFUs within a cell-laden hydrogel, accounting for their random dispersion. It assumes either that the gel restricts microbial mobility or that microbes aggregate into spherical microcolonies over time. This tool is essential for studying how bacterial colonies expand in a controlled hydrogel environment and helps predict microbial behavior under various conditions.

### Arguments:

`--growth_factor` (float, required): Specifies the rate at which the microbial population grows (fold increase per time step).

`--time_steps` (int, optional, default: 100): The number of time steps to simulate the colony growth.

`--initial_concentration` (float, optional, CFU/ml): The initial microbial concentration at the start of the simulation.

`--hydrogel_shape` (str, optional, default: Sphere): The shape of the hydrogel in which the microbial growth is modeled. Options: Sphere, Cylinder, CubeRect, CubeSquare.

`--density_threshold` (float, optional, default: 0.9): The population density at which growth stops (as a fraction of the hydrogel volume).

`--boundary_conditions` (str, optional, default: Periodic): The boundary conditions of the simulation. Options: Periodic, Reflective.

### Example Usage:

To model the growth of an E. coli colony within a sphere-shaped hydrogel:

`python pyCellGel.py model_growth --growth_factor 1.2 --time_steps 200 --hydrogel_shape Sphere --initial_concentration 1e6 --boundary_conditions Periodic`
