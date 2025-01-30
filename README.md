# pyCellGel

## Description
Script that can model and solve _in silico_ cell-laden hydrogels based on user input.

## Main Commands
**Solver** (_solver_): Computes the required loading concentrations to achieve a specified inter-CFU distribution within the 
hydrogel.

**Modeller** (_modeller_): Models hydrogels in different user-defined shapes and concentrations.

**Growth Modeller** (_model_growth_): Models the growths of microbial CFU's randomly dispersed through a cell-laden hydrogel. Assumes the gel restricts mobility or that microbes form spherical microcolonies. 

## Input
Each command is accessed through a command-line argument specifying the mode, followed by additional flags for 
customization. Most commands are fitted with default values, make sure to change according to your needs. 

For in-depth explanation on all the flags within each command, use --help or -h. 

## Usage | Solver

`python pyCellGel solver`

Calculates the required bacterial concentration to achieve a given target inter-CFU distance within the hydrogel. This is useful for ensuring uniform distribution in experimental setups.

### Arguments:

`--P_target` (float, required): Population fraction required to satisfy the target distance.

`--r_target` (float, required, μm): Desired inter-CFU distance (converted internally to mm).

`--model_prediction` (bool, optional, default: False): If True, a cilindrical hydrogel is additionally modeled based on the computed parameters.

### Example Usage:

To model a gel in which 95 % of the population satisfies a mean nearest neighbour CFU distance of 50 microns, run the following command:

`python pyCellGel.py solver --P_target 0.95 --r_target 50 --model_prediction True`

## Usage | Modeller

`python pyCellGel modeller`

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

To construct a rectangular cube of 1 x 0.5 x 2 mm, with a optical density of 0.01 containing randomly distributed CFUs, use the following arguments:

`python pyCellGel.py modeller --shape CubeRect --height 1.0 --width 0.5 --length 2.0 --optical_density 0.01`

## Usage | Growth Modeller

`python pyCellGel model_growth`

The growth modeller command simulates the growth of microbial CFUs within a cell-laden hydrogel, accounting for their random dispersion. It assumes either that the gel restricts microbial mobility or that microbes aggregate into spherical microcolonies over time. This tool is essential for studying how bacterial colonies expand in a controlled hydrogel environment and helps predict microbial behavior under various conditions.

### Output Handling

For each command, an output directory (output_<command>) is created to store results. If an existing directory is found, it is removed and recreated to ensure clean output generation.

### Example Usage:

To model the growth of an E. coli colony within a sphere-shaped hydrogel:

`python pyCellGel.py model_growth --growth_factor 1.2 --time_steps 200 --hydrogel_shape Sphere --initial_concentration 1e6 --boundary_conditions Periodic`
