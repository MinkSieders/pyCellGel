# pyCellGel
Script that can model and solve _in silico_ cell-laden hydrogels.

1. Solver (solver)

The solver command calculates the required bacterial concentration to achieve a given target inter-CFU distance within the hydrogel. This is useful for ensuring uniform distribution in experimental setups.

Arguments:

--P_target (float, required): Population fraction required to satisfy the target distance.

--r_target (float, required, μm): Desired inter-CFU distance (converted internally to mm).

--model_prediction (bool, optional, default: False): If True, a hydrogel is also modeled based on the computed parameters.

Example Usage:

python pyCellGel.py solver --P_target 0.5 --r_target 50 --model_prediction True

2. Modeller (modeller)

The modeller command allows users to define hydrogel shapes and dimensions while specifying bacterial concentrations either directly or via an optical density (OD600) measurement.

Shape Options:

Cylinder

Sphere

CubeSquare

CubeRect

Arguments:

--shape (str, optional, default: Sphere): Shape to be modeled. Can be one of the shape options or All.

--height (float, optional, mm, default: 0.5): Height of the hydrogel (for applicable shapes).

--radius (float, optional, mm, default: 1.0): Radius of the hydrogel (for spheres and cylinders).

--width (float, optional, mm, default: 0.75): Width of the hydrogel (for cuboid shapes).

--length (float, optional, mm, default: 1.0): Length of the hydrogel (for cuboid shapes).

--concentration (float, optional, CFU/ml): Manually specify bacterial concentration. If provided, OD600 conversion is skipped.

--optical_density (float, optional, default: 0.005): Specify concentration via OD600 (requires OD-to-CFU conversion curve; defaults to E. coli if not provided).

--optical_density_curve (str, optional): Path to a .tsv file containing OD vs. CFU/ml data. If omitted, the script defaults to an internal E. coli calibration curve.

Example Usage:

python pyCellGel.py modeller --shape CubeRect --height 1.0 --width 0.5 --length 2.0 --optical_density 0.01

3. Figure Generation (figure)

The figure command generates visualization figures based on the computed or modeled hydrogel data.

Arguments:

(No additional arguments required.)

Example Usage:

python pyCellGel.py figure

Output Handling

For each command, an output directory (output_<command>) is created to store results. If an existing directory is found, it is removed and recreated to ensure clean output generation.
