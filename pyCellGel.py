import mpl_toolkits.mplot3d.art3d as art3d
from scipy.spatial import distance
from itertools import combinations
from itertools import product
import os
import shutil
from scipy.special import gamma
import scipy.integrate as integrate
from scipy.optimize import root_scalar
import argparse
from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle

"""
pyCellGel
Version: 3
Date: January 30 2025
Created: December 2024
Author: Mink Sieders

Version History:
V0-2: Incorporated Modeller, Solver. 
V3: Incorporated growth-dynamics within a hydrogel, simulating very basic growth on equally dispersed CFU's in a hydrogel.
The above dynamics assume microbes stay in place and are physically restricted to their starting location. 
"""

fig_dDist, ax_dDist = None, None

# Determine constants for theoretical distance determination
k_sphere = (3 / (4 * np.pi)) ** (1 / 3) * gamma((4 / 3))  # We assume this constant also is good for other shapes

def plot_SummaryFig(shapes):

    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import letter

    # Function to add an image (PNG or SVG) to a PDF
    def add_image_to_pdf(pdf_canvas, img_path, x, y, max_width, max_height):
        pdf_canvas.drawImage(img_path, x, y, width=max_width, height=max_height, preserveAspectRatio=True)

    for shape in shapes:
        keywords = [shape, "plot_dcurve", "plot_pdist", "plot_dDist"]
        keywords = [str(item).lower() for item in keywords]
        image_paths = [f for f in os.listdir('.') if
                       any(keyword in f.lower() for keyword in keywords) and
                       (f.endswith('.png') or f.endswith('.svg'))]

        # Separate "whole" image from others
        whole_image = next((img for img in image_paths if "whole" in img.lower()), None)
        other_images = [img for img in image_paths if img != whole_image]

        # Create a PDF canvas
        pdf_path = f"summary_{shape}_model.pdf"
        pdf_canvas = canvas.Canvas(pdf_path, pagesize=letter)
        page_width, page_height = letter

        # Layout for the "whole" image (top half)
        if whole_image:
            top_max_width = page_width - 100
            top_max_height = page_height / 2 - 50
            add_image_to_pdf(pdf_canvas, whole_image, 50, page_height / 2 + 25, top_max_width, top_max_height)

        # Layout for the remaining images (bottom half)
        num_remaining = len(other_images)
        if num_remaining > 0:
            section_height = (page_height / 2 - 50) / num_remaining
            bottom_max_width = page_width - 100
            bottom_max_height = section_height - 20
            y_position = page_height / 2 - 25

            for img_path in other_images:
                add_image_to_pdf(pdf_canvas, img_path, 50, y_position - bottom_max_height, bottom_max_width,
                                 bottom_max_height)
                y_position -= section_height

        pdf_canvas.save()


def performScale(pos, ax):
    x_range = pos[:, 0].max() - pos[:, 0].min()
    y_range = pos[:, 1].max() - pos[:, 1].min()
    z_range = pos[:, 2].max() - pos[:, 2].min()
    x_scale = x_range
    y_scale = y_range
    z_scale = z_range
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj = short_proj

    return ax


def plot_WholeGel(shape, CFU_num, gel_radius, gel_height, gel_width, gel_length, gel_volume, gel_concentration, omit_output=False):

    if CFU_num < 1:
        print(f"Warning: The current simulation for {shape} can not be performed with N (total CFUs) > 0, "
              f"increase gel size or concentration CFUs")
        return

    # Create the 3D plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    if shape == "Cylinder":
        # Randomly generate CFU positions within the cylindrical gel
        CFU_positions = []

        for _ in range(int(CFU_num)):
            # Randomly generate a point inside the cylinder
            z = np.random.uniform(0, gel_height)  # height (z-axis)
            r = np.sqrt(np.random.uniform(0, gel_radius**2))  # radius, squared for uniform distribution
            theta = np.random.uniform(0, 2 * np.pi)  # angle in polar coordinates
            x = r * np.cos(theta)  # convert to Cartesian (x)
            y = r * np.sin(theta)  # convert to Cartesian (y)
            CFU_positions.append((x, y, z))

        # Convert to a NumPy array for convenience
        CFU_positions = np.array(CFU_positions)

        # Plot cylinder surfaces
        resolution = 100
        x = np.linspace(-gel_radius, gel_radius, resolution)
        z = np.linspace(0, gel_height, resolution)
        X, Z = np.meshgrid(x, z)
        Y_top = np.sqrt(gel_radius**2 - X**2)  # Top surface
        Y_bottom = -Y_top  # Bottom surface

        ax.plot_surface(X, Y_top, Z, color='lightpink', alpha=0.5, linewidth=0)
        ax.plot_surface(X, Y_bottom, Z, color='lightpink', alpha=0.5, linewidth=0)

        # Add circular caps
        floor = Circle((0, 0), gel_radius, color='lightpink', alpha=0.5)
        ax.add_patch(floor)
        art3d.pathpatch_2d_to_3d(floor, z=0, zdir="z")

        ceiling = Circle((0, 0), gel_radius, color='lightpink', alpha=0.5)
        ax.add_patch(ceiling)
        art3d.pathpatch_2d_to_3d(ceiling, z=gel_height, zdir="z")

        # Plot the randomly distributed CFUs
        ax.scatter(CFU_positions[:, 0], CFU_positions[:, 1], CFU_positions[:, 2],
                   s=1, alpha=0.8, color='black', label='CFU', edgecolors='none')

        # Set plot limits
        ax.set_xlim([-gel_radius, gel_radius])
        ax.set_ylim([-gel_radius, gel_radius])
        ax.set_zlim([0, gel_height])

        ax = performScale(CFU_positions, ax)

    elif shape == "Sphere":
        # Randomly generate CFU positions within the spherical gel
        CFU_positions = []
        for _ in range(int(CFU_num)):
            # Generate a random point in spherical coordinates
            r = np.cbrt(np.random.uniform(0, gel_radius ** 3))  # radius, cubed root for uniform distribution in sphere
            theta = np.random.uniform(0, 2 * np.pi)  # angle in xy-plane
            phi = np.arccos(np.random.uniform(-1, 1))  # angle from z-axis

            # Convert to Cartesian coordinates
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)

            CFU_positions.append((x, y, z))

        # Convert to a NumPy array for convenience
        CFU_positions = np.array(CFU_positions)

        # Plot sphere surface
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x_sphere = gel_radius * np.outer(np.cos(u), np.sin(v))
        y_sphere = gel_radius * np.outer(np.sin(u), np.sin(v))
        z_sphere = gel_radius * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightpink', alpha=0.5)

        # Plot the randomly distributed CFUs
        ax.scatter(CFU_positions[:, 0], CFU_positions[:, 1], CFU_positions[:, 2],
                   s=1, alpha=0.8, color='black', label='CFU', edgecolors='none')

        # Set plot limits
        ax.set_xlim([-gel_radius, gel_radius])
        ax.set_ylim([-gel_radius, gel_radius])
        ax.set_zlim([-gel_radius, gel_radius])

    elif shape == "CubeSquare":
        # Randomly generate CFU positions within the cubic gel
        CFU_positions = []
        half_side = gel_height / 2  # Half the side length for boundaries
        for _ in range(int(CFU_num)):
            x = np.random.uniform(-half_side, half_side)
            y = np.random.uniform(-half_side, half_side)
            z = np.random.uniform(-half_side, half_side)
            CFU_positions.append((x, y, z))

        # Convert to a NumPy array for convenience
        CFU_positions = np.array(CFU_positions)

        # Plot the cube surface
        r = [-half_side, half_side]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s - e)) == r[1] - r[0]:
                ax.plot3D(*zip(s, e), color="lightpink", alpha=0.5)

        # Plot the randomly distributed CFUs
        ax.scatter(CFU_positions[:, 0], CFU_positions[:, 1], CFU_positions[:, 2],
                   s=1, alpha=0.8, color='black', label='CFU', edgecolors='none')

        # Set plot limits
        ax.set_xlim([-half_side*1.5, half_side*1.5])
        ax.set_ylim([-half_side*1.5, half_side*1.5])
        ax.set_zlim([-half_side, half_side])

    elif shape == "CubeRect":
        # Randomly generate CFU positions within the rectangular cuboid gel
        CFU_positions = []
        for _ in range(int(CFU_num)):
            x = np.random.uniform(-gel_length / 2, gel_length / 2)  # Centered x-dimension
            y = np.random.uniform(-gel_width / 2, gel_width / 2)  # Centered y-dimension
            z = np.random.uniform(-gel_height / 2, gel_height / 2)  # Centered z-dimension
            CFU_positions.append((x, y, z))

        # Convert to a NumPy array for convenience
        CFU_positions = np.array(CFU_positions)

        # Create plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Define centered ranges
        r = [-gel_length / 2, gel_length / 2]
        w = [-gel_width / 2, gel_width / 2]
        h = [-gel_height / 2, gel_height / 2]

        # Plot the cuboid surface
        for s, e in combinations(np.array(list(product(r, w, h))), 2):
            if np.sum(np.abs(s - e)) == gel_length or np.sum(np.abs(s - e)) == gel_width or np.sum(
                    np.abs(s - e)) == gel_height:
                ax.plot3D(*zip(s, e), color="lightpink", alpha=0.5)

        # Plot the randomly distributed CFUs
        ax.scatter(CFU_positions[:, 0], CFU_positions[:, 1], CFU_positions[:, 2],
                   s=1, alpha=0.8, color='black', label='CFU', edgecolors='none')

        # Set plot limits
        ax.set_xlim([(-gel_length / 2)*1.5, (gel_length / 2)*1.5])
        ax.set_ylim([(-gel_width / 2)*1.5, (gel_width / 2)*1.5])
        ax.set_zlim([-gel_height / 2, gel_height / 2])

        ax = performScale(CFU_positions, ax)

    if not omit_output:
        # Add shape-specific parameters as text annotations
        ax.text2D(0.05, 0.95, f"{shape}", transform=ax.transAxes, fontsize=10, color='black',
                verticalalignment='top')
        ax.text2D(0.05, 0.90, f"c: {gel_concentration:.2e} (CFU / ml)", transform=ax.transAxes, fontsize=10, color='black',
                verticalalignment='top')
        ax.text2D(0.05, 0.85, f"N: {CFU_num:.2e} (total CFU)", transform=ax.transAxes, fontsize=10, color='black',
                verticalalignment='top')
        ax.text2D(0.05, 0.80, f"V: {gel_volume:.3f} (mm³ or μl)", transform=ax.transAxes, fontsize=10, color='black',
                verticalalignment='top')

        if shape == "Cylinder":
            ax.text2D(0.05, 0.75, f"r: {gel_radius} mm", transform=ax.transAxes, fontsize=10, color='black',
                    verticalalignment='top')
            ax.text2D(0.05, 0.70, f"h: {gel_height} mm", transform=ax.transAxes, fontsize=10, color='black',
                    verticalalignment='top')

        elif shape == "Sphere":
            ax.text2D(0.05, 0.75, f"r: {gel_radius} mm", transform=ax.transAxes, fontsize=10, color='black',
                      verticalalignment='top')

        elif shape == "CubeSquare":
            ax.text2D(0.05, 0.75, f"l: {gel_height} mm", transform=ax.transAxes, fontsize=10, color='black',
                      verticalalignment='top')

        elif shape == "CubeRect":
            ax.text2D(0.05, 0.75, f"l: {gel_length} mm", transform=ax.transAxes, fontsize=10, color='black',
                      verticalalignment='top')
            ax.text2D(0.05, 0.70, f"w: {gel_width} mm", transform=ax.transAxes, fontsize=10, color='black',
                      verticalalignment='top')
            ax.text2D(0.05, 0.65, f"h: {gel_height} mm", transform=ax.transAxes, fontsize=10, color='black',
                      verticalalignment='top')

        # Finalize Layout
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_zlabel("z (mm)")
        plt.savefig(f'modelled_{shape}_whole.png', dpi=1200)
        plt.close()

    plt.close()
    return CFU_positions


def plot_dDist(shape, positions, concentration, all_shapes, k=k_sphere):

    global fig_dDist, ax_dDist

    theoretical_distance = (k * (1000/concentration) ** (1 / 3)
                            * 1000)     # add latest *1000 to convert to micrometers

    print(f"Theoretical ⟨d⟩ (μm) of perfectly distributed CFUs: {theoretical_distance:.1f}")

    if len(positions) < 50_000:
        distances = distance.cdist(positions, positions, metric='euclidean')
        distances = distances * 1000

        # Find the closest distance for each point
        np.fill_diagonal(distances, np.inf)  # Ignore self-distance (set diagonal to infinity)
        closest_distances = np.min(distances, axis=1)  # Minimum distance for each point

        n_points = len(closest_distances)
        jitter = np.random.normal(1, 0.04, n_points)  # Add small jitter for scatter points

        # Create a boxplot for all distances
        if fig_dDist is None and ax_dDist is None:
            fig_dDist, ax_dDist = plt.subplots(figsize=(6, 8))

        # Add a boxplot for the current shape
        ax_dDist.boxplot(
            [closest_distances],
            positions=[all_shapes.index(shape)],  # Set x location
            widths=0.6,
            vert=True,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(facecolor='none', color='black'),
        )
        ax_dDist.scatter(
            jitter + all_shapes.index(shape) - 1,  # Scatter points offset from the box
            closest_distances,
            color='orange',
            alpha=0.3,
            label='Raw data' if shape == all_shapes[0] else None,  # Add label only once
            s=3,
        )

        ax_dDist.axhline(
            y=theoretical_distance,
            color='red',
            linestyle='-',
            label="Theoretical ⟨d⟩" if shape == all_shapes[0] else None,  # Add label only once
        )

        # Add x-tick for the current shape
        ax_dDist.set_xticks(range(len(all_shapes)))
        ax_dDist.set_xticklabels(all_shapes)

        ax_dDist.set_ylabel("d (μm)")
        ax_dDist.grid(alpha=0.3, axis='y')
        ax_dDist.legend()

        print(
            f"In silico ⟨d⟩ (μm) of {n_points:.2e} CFUs in gel: {np.mean(closest_distances):.2f}")

        # If this is the last shape, save the figure
        if shape == all_shapes[-1]:
            fig_dDist.tight_layout()
            fig_dDist.savefig("plot_dDist.png", dpi=1200)
            plt.close(fig_dDist)
            fig_dDist, ax_dDist = None, None  # Reset for future calls

    else:
        print(f"Warning: too many points {len(positions):.2e} in {shape} simulated hydrogel, in silico ⟨d⟩ can not be determined")


def plot_GelSlice(shape, CFU_positions, available_space_per_microbe, omit_output=False):
    cube_root_space = ((available_space_per_microbe ** (
                1 / 3)) / 2) * 4  # Defines the area where you will find at least 1 CFU (probability based) and do it times 4 in all directions

    # Filter CFU positions within the cube-root bounds
    filtered_CFUs = CFU_positions[
        (CFU_positions[:, 0] >= -cube_root_space) & (CFU_positions[:, 0] <= cube_root_space) &
        (CFU_positions[:, 1] >= -cube_root_space) & (CFU_positions[:, 1] <= cube_root_space) &
        (CFU_positions[:, 2] >= -cube_root_space) & (CFU_positions[:, 2] <= cube_root_space)
        ]

    filtered_CFUs = filtered_CFUs*1000
    cube_root_space = cube_root_space*1000

    if omit_output:
        return filtered_CFUs, cube_root_space

    # Create the 3D plot for the filtered CFUs
    fig_slice = plt.figure(figsize=(8, 6))
    ax_slice = fig_slice.add_subplot(111, projection='3d')

    # Plot the filtered CFUs
    ax_slice.scatter(filtered_CFUs[:, 0], filtered_CFUs[:, 1], filtered_CFUs[:, 2],
                     s=10, alpha=0.8, color='black', label='CFU', edgecolors='none')

    # Iterate through each CFU and find its closest neighbor, then plot a line with distance annotation
    for i in range(len(filtered_CFUs)):
        # Calculate distances from CFU[i] to all other CFUs
        distances = distance.cdist([filtered_CFUs[i]], filtered_CFUs, 'euclidean')
        distances = distances[0]  # Flatten the result

        # Exclude the distance to itself (which is 0)
        distances[i] = np.inf  # Set its own distance to infinity to avoid picking itself

        # Find the index of the closest CFU
        closest_idx = np.argmin(distances)
        closest_CFU = filtered_CFUs[closest_idx]

        # Get the distance to the closest CFU
        min_distance = distances[closest_idx]

        # Plot a line between the CFU and its closest CFU
        ax_slice.plot([filtered_CFUs[i][0], closest_CFU[0]],
                      [filtered_CFUs[i][1], closest_CFU[1]],
                      [filtered_CFUs[i][2], closest_CFU[2]], color='black', alpha=0.5, linestyle="dashed")

        # Annotate the line with the distance
        props = dict(facecolor='white', alpha=0.0, linewidth=0)
        mid_point = (filtered_CFUs[i] + closest_CFU) / 2
        ax_slice.text(mid_point[0], mid_point[1], mid_point[2],
                      f'{min_distance:.0f} μm', color='darkblue', fontsize=4, ha='center', bbox=props)

    # Set plot limits to cube root space
    ax_slice.set_xlim([-cube_root_space, cube_root_space])
    ax_slice.set_ylim([-cube_root_space, cube_root_space])
    ax_slice.set_zlim([-cube_root_space, cube_root_space])

    # Label axes
    ax_slice.set_xlabel("x (μm)")
    ax_slice.set_ylabel("y (μm)")
    ax_slice.set_zlabel("z (μm)")
    ax_slice.legend()

    ax_slice.xaxis._axinfo['grid'].update(color='#FF1493', linewidth=0.5)  # Dark pink grid lines
    ax_slice.yaxis._axinfo['grid'].update(color='#FF1493', linewidth=0.5)
    ax_slice.zaxis._axinfo['grid'].update(color='#FF1493', linewidth=0.5)

    plt.savefig(f'modelled_{shape}_slice.png', dpi=1200)
    plt.close()


def print_Parameters(s, v, tcfu, scfu):
    print(f"Derived Parameters {s} Gel:")
    print(f"Gel volume (mm³ or μl): {v:.3f}")
    print(f"Total CFUs in gel: {tcfu:.2e}")
    print(f"Available space per microbe in gel (mm³): {scfu:.2e}")


def model_Gels(CFU_concentration, height, radius, width, length, shapes, omit_output=False):

    if not isinstance(shapes, list):
        shapes = [shapes]

    for shape in shapes:
        if shape == "Cylinder":
            volume = np.pi * radius ** 2 * height  # mm³

        elif shape == "Sphere":
            volume = (4 / 3) * np.pi * radius ** 3  # mm³

        elif shape == "CubeSquare":
            volume = height ** 3  # mm³

        elif shape == "CubeRect":
            volume = height * width * length  # mm³

        else:
            print("Warning: No shape defined, omitting in silico hydrogel construction")
            return

        total_CFUs = CFU_concentration * volume * 1e-3  # total CFUs
        CFU_space = volume / total_CFUs  # mm³ per CFU
        if not omit_output:
            print_Parameters(shape, volume, total_CFUs, CFU_space)

        CFU_positions = plot_WholeGel(shape, total_CFUs, radius, height, width, length, volume, CFU_concentration, omit_output)

        if CFU_positions is not None:
            if not omit_output:

                plot_dDist(shape, CFU_positions, CFU_concentration, shapes)

                plot_GelSlice(shape, CFU_positions, CFU_space)

                print("\n")

        return CFU_positions, CFU_space


def plot_dCurve(conc, k=k_sphere):
    plt.figure(figsize=(8, 6))
    x = np.linspace(0, 1e8, 10_000)
    x = x[x != 0]
    y = k*(1000/x)**(1/3)*1000
    input = k*(1000/conc)**(1/3)*1000

    plt.plot(x, y, color='k', label='Model')  # or 'Data' if you prefer
    plt.scatter(conc, input, color='darkred', marker='x', s=100, label='concentration model')
    plt.axhline(y=input, color='gray', linestyle='--', linewidth=0.7)
    plt.axvline(x=conc, color='gray', linestyle='--', linewidth=0.7)
    plt.annotate(f'{conc:.2e}', (conc, 0), textcoords="offset points", xytext=(0, 0), ha='center', color='gray')
    plt.annotate(f'{input:.2f}', (0, input), textcoords="offset points", xytext=(0, 0), va='center', color='gray')

    plt.xlabel("c (CFU / mL)")
    plt.ylabel("⟨d⟩ (µm)")
    plt.legend()
    plt.savefig('plot_dCurve.png', dpi=1200)
    plt.close()


def plot_pDist(conc, pt=0.95):
    ft = 1 - pt
    vol = 1.571  # hypothetical volume (mm^3)
    Num = conc * vol * 1e-3  # Total CFUs in the volume
    rho = Num / vol  # Point density

    # Define the probability distribution function P(r)
    def P(r, rho):
        normalization = quad(lambda r: 4 * np.pi * r ** 2 * rho * np.exp(- (4 / 3) * np.pi * r ** 3 * rho), 0, np.inf)[0]
        return (4 * np.pi * r ** 2 * rho * np.exp(- (4 / 3) * np.pi * r ** 3 * rho)) / normalization

    i = 1000
    j = 500_000

    r_values = np.linspace(0, i, j)
    P_values = P(r_values, rho)

    # Cumulative integration to find r_lower
    dx = i / j
    cumulative_area = np.cumsum(P_values) * dx
    r_lower_index = np.argmax(cumulative_area >= ft)
    r_middle_index = np.argmax(cumulative_area >= 0.5)
    r_end_index = np.argmax(cumulative_area >= 0.9999)
    r_lower = r_values[r_lower_index]
    r_middle = r_values[r_middle_index]
    r_end = r_values[r_end_index]

    # Plot the probability distribution
    fig, ax = plt.subplots(figsize=(8, 6))

    r_values_vis = r_values * 1000  # Distance in micrometers
    r_lower_vis = r_lower * 1000  # Distance in micrometers
    r_middle_vis = r_middle * 1000  # Distance in micrometers
    r_end_vis = r_end * 1000  # Distance in micrometers

    plt.plot(r_values_vis, P_values, color="b", label="P(r)")
    plt.fill_between(r_values_vis, 0, P_values, where=(r_values_vis >= r_lower_vis), color='yellow', alpha=0.3)
    plt.axvline(x=r_lower_vis, color="r", linestyle="--", label=f"r_lower = {r_lower_vis:.2f} µm")
    plt.text(r_middle_vis, max(P_values)/2, f"P(r > {r_lower_vis:.2f} µm)\n≈ {pt:.2f}", horizontalalignment='center')

    plt.xlim(0, r_end_vis)
    plt.xlabel("⟨d⟩ (µm)")
    plt.ylabel("P(r)")
    plt.grid(True)
    plt.savefig("plot_pDist.png", dpi=1200)
    plt.close()


def determineConcentration(rt, Pt):
    """
    Units for calculation:
    """
    # Define constants
    F_target = 1 - Pt

    # Define the integrand for the PDF
    def integrand(r_prime, rho):
        return 4 * np.pi * r_prime ** 2 * rho * np.exp(-4 / 3 * np.pi * r_prime ** 3 * rho)

    # Define the CDF function for a given rho
    def cdf(r, rho):
        result, _ = integrate.quad(integrand, 0, r, args=(rho,))
        return result

    # Define the objective function for Brent's method
    def objective(rho):
        return cdf(rt, rho) - F_target

    def brentSolver():
        try:
            # Use Brent's method to find rho
            solution = root_scalar(objective, bracket=[1e-10, 1e10], method='brentq')

            # Extract the solution for rho (in CFU/mm³)
            rho_solution = solution.root

            # Convert rho to concentration (C)
            concentration = rho_solution * 1e3  # Convert CFU/mm³ to CFU/ml
            return concentration

        except ValueError as e:
            # Handle the case where Brent's method fails (e.g., f(a) and f(b) have the same sign)
            print(f"Brent's method failed with error: {e}")
            print("Falling back to alternative method...")

            # Fallback method: Perform a simple search (e.g., binary search)
            return alternativeMethod()

    def alternativeMethod():
        # Perform a simple search for rho by testing different values (binary search approach)
        low, high = 1e-10, 1000
        tolerance = 1e-6
        while high - low > tolerance:
            mid = (low + high) / 2
            if cdf(rt, mid) > F_target:
                high = mid
            else:
                low = mid

        # Estimate the concentration using the midpoint (rho) value
        rho_solution = (low + high) / 2
        concentration = rho_solution * 1e3  # Convert CFU/mm³ to CFU/ml
        return concentration

    # Call brentSolver with fallback
    concentration = brentSolver()

    print(f"Required c (CFU / ml) to satisfy P(r > {rt * 1000:.2f} µm) ≈ {Pt:.2f}: {concentration:.4e}")

    return concentration


def main_ConcentrationSolver(P_target, r_target, make_gel_model):

    concentration_desired = determineConcentration(r_target, P_target)

    plot_pDist(concentration_desired, P_target)

    plot_dCurve(concentration_desired)

    if make_gel_model:
        model_Gels(concentration_desired, height=5.0176, radius=7.124, width=None, length=None, shapes=["Cylinder"])

        plot_SummaryFig(["Cylinder"])


def main_GelModeller(user_defined_shape, input_concentration, input_height, input_radius, input_width, input_length, omit_output = False):
    # Obtain the shapes
    if user_defined_shape in shape_options:
        shapes_to_simulate = [user_defined_shape]
    elif user_defined_shape == "All":
        shapes_to_simulate = shape_options
    else:
        print(f"Warning: No valid shape provided, options are: {shape_options} or 'All' to produce all options at once")
        shapes_to_simulate = None

    if shapes_to_simulate is not None:
        print(f"Simulating gels with shapes {shapes_to_simulate}")

        positions, space = model_Gels(input_concentration, input_height, input_radius, input_width, input_length, shapes_to_simulate, omit_output)

        if not omit_output:
            plot_dCurve(input_concentration)

            plot_pDist(input_concentration)

            plot_SummaryFig(shapes_to_simulate)

    else:
        print(f"Warning: Shape was not provided, exiting modeller")

    return positions, space


def figure_1():
    # Define the ranges of r_target and Pt
    r_targets = np.linspace(0.01, 0.15, 150)  # Radii in mm
    Pt_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]  # Population fraction requirements

    # Prepare storage for results
    results = []

    for Pt in Pt_values:
        concentrations = []
        for rt in r_targets:
            concentration = determineConcentration(rt, Pt)
            concentrations.append(concentration)
        results.append(concentrations)

    # Convert results to a numpy array for easy plotting
    results = np.array(results)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Use viridis colormap to assign colors
    viridis_colormap = plt.cm.viridis(np.linspace(0, 1, len(Pt_values)))

    # Plot each line corresponding to a different Pt, with colors from the viridis colormap
    for idx, Pt in enumerate(Pt_values):
        plt.plot(r_targets * 1000, results[idx], label=f'{Pt:.2f}', color=viridis_colormap[idx])

    # Add labels, title, and legend
    plt.ylim(0, 5e6)
    plt.xlabel('d (μm)', fontsize=12)
    plt.ylabel('Concentration (CFU/ml)', fontsize=12)
    plt.legend(title="P(r > d)", fontsize=10)
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig("d_vs_c_vs_p>r.svg")


def Compute_ODConcentration(concentration, optical_density, file_for_od_curve):
    if concentration != None:
        return(concentration)  # CFUs/mL

    else:
        if file_for_od_curve != None:
            curve_file = file_for_od_curve
        else:
            data = [
                (1, 800000000),
                (0.8, 640000000),
                (0.4, 320000000),
                (0.03, 24000000),
                (0.02, 16000000),
                (0.1, 80000000),
                (0.01, 8000000),
            ]

            with open("curve.tsv", "w") as file:
                for row in data:
                    file.write(f"{row[0]}\t{row[1]}\n")
            curve_file = "curve.tsv"

        curve_data = np.loadtxt(curve_file, delimiter="\t")
        OD_values = curve_data[:, 0]
        concentration_values = curve_data[:, 1]
        interp_func = interp1d(OD_values, concentration_values, kind="linear", fill_value="extrapolate")

        return(interp_func(optical_density))


def plot_InitialPositionMicrobes(positions, length, name):
    # Cube limits
    limit = length
    points = positions

    # Microbe properties
    rod_length = 2
    rod_radius = 0.5

    # Create figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_xlabel('Gel X (μm)')
    ax.set_ylabel('Gel Y (μm)')
    ax.set_zlabel('Gel Z (μm)')
    ax.set_title(f'{name}')

    def plot_sphere(ax, center, radius, resolution=12, color='black'):
        u, v = np.mgrid[0:2 * np.pi:resolution * 1j, 0:np.pi:resolution * 1j]
        x = center[0] + radius * np.cos(u) * np.sin(v)
        y = center[1] + radius * np.sin(u) * np.sin(v)
        z = center[2] + radius * np.cos(v)

        ax.plot_surface(x, y, z, color=color, linewidth=0, alpha=1)

    # Function to create a random unit vector
    def random_unit_vector():
        vec = np.random.randn(3)
        return vec / np.linalg.norm(vec)

    # Function to plot a 3D cylinder between two points
    def plot_3D_cylinder(ax, start, end, radius, resolution=8, color='black'):
        # Cylinder axis direction
        direction = end - start
        height = np.linalg.norm(direction)  # Length of the cylinder

        # Generate cylinder points along the Z-axis (default orientation)
        theta = np.linspace(0, 2 * np.pi, resolution)
        z = np.linspace(0, height, resolution)
        theta, z = np.meshgrid(theta, z)

        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        z = z  # Stays along Z-axis initially

        # Stack arrays for transformation
        cylinder_points = np.vstack((x.flatten(), y.flatten(), z.flatten()))

        # Create rotation matrix to align cylinder along the rod direction
        def rotation_matrix(v1, v2):
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            axis = np.cross(v1, v2)
            angle = np.arccos(np.dot(v1, v2))
            if np.linalg.norm(axis) < 1e-6:
                return np.eye(3)  # No rotation needed
            axis = axis / np.linalg.norm(axis)
            K = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]])
            R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
            return R

        # Rotate cylinder from Z-axis to the actual rod direction
        R = rotation_matrix(np.array([0, 0, 1]), direction)
        rotated_points = R @ cylinder_points

        # Reshape back into meshgrid form
        x_rot = rotated_points[0].reshape(resolution, resolution) + start[0]
        y_rot = rotated_points[1].reshape(resolution, resolution) + start[1]
        z_rot = rotated_points[2].reshape(resolution, resolution) + start[2]

        # Plot the cylinder surface
        ax.plot_surface(x_rot, y_rot, z_rot, color=color, linewidth=0, alpha=1)

        # End caps (circles)
        for pt in [start, end]:
            circle = Circle((0, 0), radius, color=color)
            ax.add_patch(circle)
            art3d.pathpatch_2d_to_3d(circle, z=0, zdir="z")
            R_circle = rotation_matrix(np.array([0, 0, 1]), direction)
            transformed_center = R_circle @ np.array([0, 0, 0])
            circle.center = (pt[0] + transformed_center[0], pt[1] + transformed_center[1])
            art3d.pathpatch_2d_to_3d(circle, z=pt[2] + transformed_center[2], zdir="z")

    # Plot microbes as 3D cylinders
    for point in points:
        direction = random_unit_vector()
        start = point - (direction * rod_length / 2)
        end = point + (direction * rod_length / 2)

        # Draw the rod as a 3D cylinder
        plot_3D_cylinder(ax, start, end, rod_radius, color='black')

        # Draw spheres at both ends
        plot_sphere(ax, start, rod_radius, color='black')
        plot_sphere(ax, end, rod_radius, color='black')

    # Show plot
    plt.savefig(f"Initial_Gel.svg")
    plt.close()


def compute_ColonyRadius(fold_increase, initial_bacteria_count=1, packing=0.64):
    # Microbe dimensions (in µm)
    rod_radius = 0.5  # µm
    rod_length = 2  # µm

    # Compute volume of a single bacterium
    V_cylinder = np.pi * (rod_radius ** 2) * rod_length  # Volume of the cylindrical part
    V_spheres = (4/3) * np.pi * (rod_radius ** 3) * 2  # Two hemispherical caps
    V_microbe = V_cylinder + V_spheres  # Total microbe volume

    # Total number of bacteria after growth
    final_bacteria_count = initial_bacteria_count * fold_increase

    # Total volume occupied by bacteria
    V_total = final_bacteria_count * V_microbe

    # Packing efficiency (random close packing ~64%)
    packing_efficiency = packing
    V_sphere = V_total / packing_efficiency

    # Compute the minimal sphere radius that can house all bacteria
    R_colony = ((3 * V_sphere) / (4 * np.pi)) ** (1/3)

    return R_colony


def plot_GrownMicrobes(positions, colony_radius, name):

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('Gel X (μm)')
    ax.set_ylabel('Gel Y (μm)')
    ax.set_zlabel('Gel Z (μm)')
    ax.set_title(f'{name}')

    # Sphere mesh resolution
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 10)

    # Generate a unit sphere
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot each microcolony as a sphere
    for pos in positions:
        x, y, z = pos
        ax.plot_surface(
            colony_radius * sphere_x + x,
            colony_radius * sphere_y + y,
            colony_radius * sphere_z + z,
            color='lime', alpha=0.5
        )

    # Set limits based on the data
    limits = np.max(np.abs(positions)) + colony_radius
    ax.set_xlim(-limits, limits)
    ax.set_ylim(-limits, limits)
    ax.set_zlim(-limits, limits)

    plt.savefig("Grown_Gel.svg")
    plt.close()


def main_GrowthModeller(c_start, c_end):

    # For ease now we just consider a 1 mm diameter sphere. We want to first obtain the initial positions (t0).
    # We do not require special outputs, we just want the positions, so we omit_output.
    pos_start, space = main_GelModeller("Sphere", c_start, None, 0.5, None, None, omit_output=True)

    # slice we do not see the whole model because we do not need to.
    pos_slice, cube_length = plot_GelSlice(shape="Sphere", CFU_positions=pos_start, available_space_per_microbe=space, omit_output=True) # NOTE: Output is in um space

    plot_InitialPositionMicrobes(pos_slice, cube_length, name=f"Initial Distribution {c_start:.1e} CFU/ml")

    c_fold_increase = round(c_end / c_start)

    print(f"Starting Concentration: {c_start:.1e} CFU/ml")
    print(f"End Concentration: {c_end:.1e} CFU/ml, {c_fold_increase} Fold Increase")

    # Now we determine the minimum radius of the formed microcolonies, based on the fold increase.
    radius_grown_colony = compute_ColonyRadius(c_fold_increase)

    print(f"Minimal grown colony radius: {radius_grown_colony:.2f} µm")

    plot_GrownMicrobes(pos_slice, radius_grown_colony, name=f"Grown Microcolonies {c_end:.1e} CFU/ml")


if __name__ == "__main__":

    shape_options = ["Cylinder", "Sphere", "CubeSquare", "CubeRect"]

    parser = argparse.ArgumentParser(description="A script that supports 'solver' and 'modeller' commands.")

    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # Solver command
    solver_parser = subparsers.add_parser("solver", help="Run the solver command")
    solver_parser.add_argument("--P_target", type=float, required=True, help="Population fraction required "
                                                                            "r_target (μm)")
    solver_parser.add_argument("--r_target", type=float, required=True, help="Specify the test")
    solver_parser.add_argument("--model_prediction", type=bool, required=False, default=False, help="Model "
                                                                                                    "a gel which follows"
                                                                                                    "the predicted "
                                                                                                    "parameters")

    # Modeller command
    modeller_parser = subparsers.add_parser("modeller", help="Run the modeller command")
    modeller_parser.add_argument("--height", type=float, required=False, default=0.5, help="height of gel that is to be modelled (mm). Default = 0.5")
    modeller_parser.add_argument("--radius", type=float, required=False, default=1.0, help="radius of gel that is to be modelled (mm). Default = 1.0")
    modeller_parser.add_argument("--width", type=float, required=False, default=0.75, help="width of gel that is to be modelled (mm). Default = 0.75")
    modeller_parser.add_argument("--length", type=float, required=False, default=1.0, help="length of gel that is to be modelled (mm). Default = 1.0")
    modeller_parser.add_argument("--concentration", type=float, required=False, help="Concentration (CFU/ml), omits calculation via OD if input is given. No default.")
    modeller_parser.add_argument("--optical_density", type=float, required=False, default=0.005, help="Concentration (OD600), default 0.005. Requires a standard OD600 vs. CFU/ml curve, "
                                                                                                      "if not given it defaults to an E. coli curve")
    modeller_parser.add_argument("--optical_density_curve", type=str, required=False, help="Accepts a .tsv file with OD and corresponding CFU/ml values, if not provided the script defaults to using"
                                                                                                            "values obtained from E. coli, the used values are saved as .tsv")
    modeller_parser.add_argument("--shape", type=str, required=False, default="Sphere", help=f"Shape to "
                                                                                             f"be modelled, choose: "
                                                                                             f"{shape_options} or 'All'")

    solver_parser = subparsers.add_parser("figure", help="Run the figure command")

    # model_growth command
    growth_parser = subparsers.add_parser("model_growth", help="Run the growth modeller command")
    growth_parser.add_argument("--concentration1", type=float, required=False,
                               help="Starting concentration (CFU/ml) of the initial hydrogel (with randomly dispersed CFU's), omits calculation via OD if input is given. No default.")
    growth_parser.add_argument("--optical_density1", type=float, required=False, default=0.005,
                               help="Concentration (OD600), default 0.005. Requires a standard OD600 vs. CFU/ml curve, "
                                    "if not given it defaults to an E. coli curve")
    growth_parser.add_argument("--optical_density_curve1", type=str, required=False,
                               help="Accepts a .tsv file with OD and corresponding CFU/ml values, if not provided the script defaults to using"
                                    "values obtained from E. coli, the used values are saved as .tsv")
    growth_parser.add_argument("--concentration2", type=float, required=False, help="Concentration (CFU/ml) after incubation, omits calculation via OD if input is given. No default.")
    growth_parser.add_argument("--optical_density2", type=float, required=False, default=1.0,
                                 help="Concentration (OD600), default 1.0. Requires a standard OD600 vs. CFU/ml curve, "
                                      "if not given it defaults to an E. coli curve")
    growth_parser.add_argument("--optical_density_curve2", type=str, required=False,
                                 help="Accepts a .tsv file with OD and corresponding CFU/ml values. Defaults to curve 1")

    # Parse all arguments
    args = parser.parse_args()

    # Setup environment
    output_dir = f"output_{args.command}"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)
    os.chdir(output_dir)

    # Run the command
    if args.command == "solver":
        in_Pt = args.P_target  # population fraction required to satisfy r_target
        in_rt = args.r_target / 1000   # r_target turned into mm from um
        in_make_model = args.model_prediction
        main_ConcentrationSolver(in_Pt, in_rt, in_make_model)

    elif args.command == "modeller":
        in_h = args.height  # mm
        in_r = args.radius  # mm
        in_w = args.width  # mm
        in_l = args.length  # mm
        in_s = args.shape

        in_c = Compute_ODConcentration(args.concentration, args.optical_density, args.optical_density_curve)

        main_GelModeller(in_s, in_c, in_h, in_r, in_w, in_l)

    elif args.command == "figure":
        figure_1()

    elif args.command == "model_growth":
        # Obtain the concentrations of t_start and t_end
        in_c_start = Compute_ODConcentration(args.concentration1, args.optical_density1, args.optical_density_curve1)

        if args.optical_density_curve2 is not None:
            in_c_end = Compute_ODConcentration(args.concentration2, args.optical_density2, args.optical_density_curve2)
        else:
            in_c_end = Compute_ODConcentration(args.concentration2, args.optical_density2, args.optical_density_curve1)

        main_GrowthModeller(in_c_start, in_c_end)
