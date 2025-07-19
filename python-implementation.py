# -*- coding: utf-8 -*-
"""
Project: Mechanics of Solids Simulation Tool  
Course: ES221 - Mechanics of Solids (Spring 2025)

Authors:
- Kavya Shah  
- Goraksh Bendale  
- Arnav Gogate  
- Apoorv Rane  

Description:
This project was developed as part of coursework for ES221 to simulate and analyze solid mechanics principles.
The original implementation has been refactored using Gemini for improved readability, maintainability, and usability.

Contributions and feedback are welcome.
"""


# Define the coordinates of the nodes for each truss element.
# Each inner list represents an element with [x1, y1, x2, y2]
elementsnew = [[-224.36, 11.42, -210.53, 1.00], [-210.53,1.00, -195.16 , 14.37] , [-195.16 , 14.37 , -180.66 , 0.92 ] , [-180.66 , 0.92 , -165.41 , 17.37] ,[-165.41 , 17.37 , -150.19 , 0.92] , [-150.19 , 0.92 , -135.37 , 20.40] , [-135.37 , 20.40 , -120.67 , 0.91] , [-120.67 , 0.91 , -105.36 , 23.17] , [-105.36 , 23.17 , -90.75 , 0.93 ] , [-90.75 , 0.93 , -75.64 , 20.43] , [-75.64 , 20.43 , -60.86 , 0.95 ] , [-60.86 , 0.95 , -45.86 , 17.42] , [-45.86 , 17.42 , -30.56 , 0.98] , [-30.56 , 0.98 , -15.29 , 14.34] , [-15.29 , 14.34 , 0.59 , 1.00 ] , [0.59 , 1.00 , 13.63 , 11.42] , [13.63 , 11.42 , -15.29 , 14.34] , [-15.29 , 14.34 , -45.86 , 17.42] , [-45.86 , 17.42 , -75.64 , 20.43 ] , [-75.64 , 20.43 , -105.36 , 23.17] , [-105.36 , 23.17 , -135.37 , 20.40] , [-135.37 , 20.40 , -165.41 , 17.37] , [-165.41 , 17.37 , -195.16 , 14.37] , [-195.16 , 14.37 , -224.36 , 11.42] , [-210.53 , 1.00 , -180.66 , 0.92] , [-180.66 , 0.92 , -150.19 , 0.92 ] , [-150.19 , 0.92 , -120.67 , 0.91] , [-120.67 , 0.91 , -90.75 , 0.93] , [-90.75 , 0.93 , -60.86 , 0.95] , [-60.86 , 0.95 , -30.56 , 0.98 ] , [-30.56 , 0.98 , 0.59 , 1.00] ]

import matplotlib.pyplot as plt
import itertools
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.widgets import Button

#==============================================================================
# CONFIGURATION AND SETUP
#==============================================================================

#------------------------------------------------------------------------------
# Global Configuration Parameters
#------------------------------------------------------------------------------

# Material Properties
YOUNGS_MODULUS = 200e3  # Young's Modulus (kN/m²)
CROSS_SECTION_AREA = 100  # Cross-sectional area of truss elements (m²)

# Visualization Parameters
DISPLACEMENT_SCALE = 10  # Scale factor for displaying deformed shape
PLOT_FIGSIZE = (15, 6)  # Default figure size for plots
NODE_MARKER_SIZE = 6  # Size of node markers in plots
ELEMENT_LINE_WIDTH = 2  # Width of element lines in plots
COLORMAP = plt.cm.coolwarm  # Colormap for stress visualization
GRID_LINE_STYLE = '--'  # Style for grid lines
GRID_LINE_WIDTH = 0.5  # Width of grid lines

# Degrees of Freedom
DOF_PER_NODE = 2  # Number of degrees of freedom per node (2 for 2D truss)

# Initialize Storage Variables
stress_nodes = []  # Will store original element coordinates for stress calculation
elements = []  # Will store detailed element information (coordinates, properties, etc.)


# Populate the 'elements' list using the coordinates from 'elementsnew'
# and calculate the length of each element
for i in range(len(elementsnew)):
  x1 = elementsnew[i][0]
  y1 = elementsnew[i][1]
  x2 = elementsnew[i][2]
  y2 = elementsnew[i][3]
  # Calculate the length of the element
  length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

  # Append element data (coordinates, area, Young's modulus, and length)
  elements.append((x1, y1, x2, y2, CROSS_SECTION_AREA, YOUNGS_MODULUS, length))

# Store the original element coordinates for later stress calculations
stress_nodes = elementsnew

#==============================================================================
# PROCESSING - STRUCTURAL ANALYSIS
#==============================================================================

def calculate_angles_and_stiffness(truss_elements):
    """
    Calculates the angle and stiffness matrix for each truss element.

    Args:
        truss_elements: A list of tuples, where each tuple contains
                        (x1, y1, x2, y2, A, E, L) for an element.

    Returns:
        A list of stiffness matrices for each element.
    """
    stiffness_matrices = []
    for (x1, y1, x2, y2, area, youngs_modulus, length) in truss_elements:
        # Calculating angle of the element with respect to global coordinate system
        delta_x = x2 - x1
        delta_y = y2 - y1

        # Calculate cosine and sine of the angle
        # Avoid division by zero if element is vertical
        if length == 0:
            cos_theta = 0
            sin_theta = 0
        else:
            cos_theta = delta_x / length
            sin_theta = delta_y / length

        # Constructing the stiffness matrix for the element in global coordinates
        k = (area * youngs_modulus / length) * np.array([[cos_theta**2, cos_theta*sin_theta, -cos_theta**2, -cos_theta*sin_theta],
                                     [cos_theta*sin_theta, sin_theta**2, -cos_theta*sin_theta, -sin_theta**2],
                                     [-cos_theta**2, -cos_theta*sin_theta, cos_theta**2, cos_theta*sin_theta],
                                     [-cos_theta*sin_theta, -sin_theta**2, cos_theta*sin_theta, sin_theta**2]])

        stiffness_matrices.append(k)

    return stiffness_matrices

def plot_truss(elements_data):
    """
    Plots the truss structure based on element coordinates.

    Args:
        elements_data: A list of tuples, where each tuple contains
                       (x1, y1, x2, y2, A, E, L) for an element.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    # Cycle through colors for better visualization
    colors = itertools.cycle(["r", "g", "b", "c", "m", "y", "k"])
    nodes = set()

    # Plot each element and collect unique node coordinates
    for element in elements_data:
        x1, y1, x2, y2, _, _, _ = element
        color = next(colors)
        ax.plot([x1, x2], [y1, y2], marker='o', linestyle='-', color=color, markersize=8, linewidth=2)
        nodes.add((x1, y1))
        nodes.add((x2, y2))

    # Labeling the nodes
    # Sort nodes for consistent labeling
    sorted_nodes = sorted(list(nodes))
    for i, (x, y) in enumerate(sorted_nodes, start=1):
        ax.text(x, y, f"{i}", fontsize=12, ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.6, edgecolor='black'))

    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_title("Truss Structure Visualization")
    ax.grid(True, linestyle="--", linewidth=0.5)
    ax.set_aspect("equal")
    plt.show()

# Assign the populated elements list to truss_elements for consistency
truss_elements = elements

# Create and store the initial plot
initial_fig = plt.figure(figsize=(8, 6))
ax = initial_fig.add_subplot(111)
# Cycle through colors for better visualization
colors = itertools.cycle(["r", "g", "b", "c", "m", "y", "k"])
nodes = set()

# Plot each element and collect unique node coordinates
for element in truss_elements:
    x1, y1, x2, y2, _, _, _ = element
    color = next(colors)
    ax.plot([x1, x2], [y1, y2], marker='o', linestyle='-', color=color, markersize=8, linewidth=2)
    nodes.add((x1, y1))
    nodes.add((x2, y2))

# Labeling the nodes
sorted_nodes = sorted(list(nodes))
for i, (x, y) in enumerate(sorted_nodes, start=1):
    ax.text(x, y, f"{i}", fontsize=12, ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.6, edgecolor='black'))

ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_title("Initial Truss Structure")
ax.grid(True, linestyle="--", linewidth=0.5)
ax.set_aspect("equal")
plt.draw()  # Draw but don't block
plt.pause(0.1)  # Small pause to ensure plot is shown

# Calculate stiffness matrices for all elements
stiffness_matrices = calculate_angles_and_stiffness(truss_elements)

# Print individual stiffness matrices for verification
for i, (x1, y1, x2, y2, A, E, L) in enumerate(truss_elements, start=1):
    print(f"\nElement {i}:")
    print(f"Coordinates: ({x1}, {y1}) to ({x2}, {y2})")
    print(f"Length (L) = {L:.4f}")
    print(f"Stiffness Matrix for Element {i}:\n{stiffness_matrices[i-1]}")

def assemble_global_stiffness(truss_elements_data, element_stiffness_matrices):
    """
    Assembles the global stiffness matrix for the truss structure.

    Args:
        truss_elements_data: A list of tuples, where each tuple contains
                             (x1, y1, x2, y2, A, E, L) for an element.
        element_stiffness_matrices: A list of individual stiffness matrices
                                    for each element.

    Returns:
        A tuple containing:
            - K_global: The assembled global stiffness matrix.
            - all_nodes: A sorted list of unique node coordinates.
    """
    # Extract all unique node coordinates from the elements
    all_nodes = []
    for (x1, y1, x2, y2, _, _, _) in truss_elements_data:
        if (x1, y1) not in all_nodes:
            all_nodes.append((x1, y1))
        if (x2, y2) not in all_nodes:
            all_nodes.append((x2, y2))

    # Sort the nodes for consistent indexing
    all_nodes = sorted(all_nodes)
    # Create a dictionary mapping node coordinates to their index
    node_index_map = {node: idx for idx, node in enumerate(all_nodes)}

    # Define degrees of freedom per node (2 for 2D truss: x and y displacement)
    dof_per_node = 2
    # Calculate the total degrees of freedom in the system
    total_dofs = len(all_nodes) * dof_per_node
    # Initialize the global stiffness matrix with zeros
    K_global = np.zeros((total_dofs, total_dofs))

    # Assembling the global stiffness matrix by adding contributions from each element's stiffness matrix
    for e_idx, element in enumerate(truss_elements_data):
        x1, y1, x2, y2, _, _, _ = element
        # Get the local stiffness matrix for the current element
        k_local = element_stiffness_matrices[e_idx]
        # Get the global node indices for the current element's start and end nodes
        n1_idx = node_index_map[(x1, y1)]
        n2_idx = node_index_map[(x2, y2)]

        # Define the mapping of local degrees of freedom to global degrees of freedom
        # For a 2D truss element with nodes i and j, the DOFs are (2i, 2i+1, 2j, 2j+1)
        dof_map = [
            n1_idx * dof_per_node, n1_idx * dof_per_node + 1,
            n2_idx * dof_per_node, n2_idx * dof_per_node + 1
        ]

        # Add the elements of the local stiffness matrix to the corresponding positions in the global stiffness matrix
        for i in range(4):
            for j in range(4):
                K_global[dof_map[i], dof_map[j]] += k_local[i, j]

    return K_global, all_nodes

# Assemble the global stiffness matrix and get the list of all nodes
K_global, all_nodes = assemble_global_stiffness(truss_elements, stiffness_matrices)

# Output the global stiffness matrix
print("\nGlobal Stiffness Matrix:")
# Set print options for better readability of the large matrix
np.set_printoptions(precision=3, suppress=True)
print(K_global)

# Initialize a zero vector for generalized displacements Q (this will be populated later)
# This variable Q is not used in this function but is present in the original code.
# It's better to initialize it where it's used. Removing for now for clarity.
# Q = np.zeros(len(K_global))

def get_force_vector(all_nodes_coords):
    """
    Creates a graphical interface for inputting forces applied at each node.

    Args:
        all_nodes_coords: A sorted list of unique node coordinates.

    Returns:
        A numpy array representing the global force vector.
    """
    dof_per_node = 2
    total_dofs = len(all_nodes_coords) * dof_per_node
    F_global = np.zeros(total_dofs)

    # Create the main window
    root = tk.Tk()
    root.title("Input Forces")
    root.geometry("600x400")

    # Create a frame to hold the scrollable content
    main_frame = ttk.Frame(root)
    main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    # Create a canvas and scrollbar
    canvas = tk.Canvas(main_frame)
    scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
    scrollable_frame = ttk.Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    # Add a header label
    header = ttk.Label(scrollable_frame, text=f"Enter forces for {len(all_nodes_coords)} nodes", font=('Arial', 12, 'bold'))
    header.grid(row=0, column=0, columnspan=4, pady=10)

    # Dictionary to store the entry widgets
    entries = {}

    # Create entry fields for each node
    for i, (x, y) in enumerate(all_nodes_coords):
        node_label = ttk.Label(scrollable_frame, text=f"Node {i+1} ({x:.2f}, {y:.2f}):")
        node_label.grid(row=i+1, column=0, padx=5, pady=5)

        fx_label = ttk.Label(scrollable_frame, text="Force X:")
        fx_label.grid(row=i+1, column=1, padx=5)
        fx_entry = ttk.Entry(scrollable_frame, width=10)
        fx_entry.insert(0, "0")
        fx_entry.grid(row=i+1, column=2, padx=5)

        fy_label = ttk.Label(scrollable_frame, text="Force Y:")
        fy_label.grid(row=i+1, column=3, padx=5)
        fy_entry = ttk.Entry(scrollable_frame, width=10)
        fy_entry.insert(0, "0")
        fy_entry.grid(row=i+1, column=4, padx=5)

        entries[i] = (fx_entry, fy_entry)

    # Variable to store whether the window was closed properly
    submit_clicked = False

    def on_submit():
        nonlocal submit_clicked
        try:
            for i, (fx_entry, fy_entry) in entries.items():
                fx = float(fx_entry.get() or 0)
                fy = float(fy_entry.get() or 0)
                F_global[i * dof_per_node] = fx
                F_global[i * dof_per_node + 1] = fy
            submit_clicked = True
            root.destroy()
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numeric values for forces.")

    # Pack the canvas and scrollbar
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    # Add submit button
    submit_btn = ttk.Button(root, text="Submit", command=on_submit)
    submit_btn.pack(pady=10)

    # Start the main loop
    root.mainloop()

    if not submit_clicked:
        raise SystemExit("Force input window was closed without submitting.")

    return F_global

# Take input from user to define the global force vector
F_global = get_force_vector(all_nodes)

# Print the resulting global force vector
print("\nGlobal Force Vector:")
print(F_global)

def apply_boundary_conditions_by_node(K_global_matrix, F_global_vector, all_nodes_coords):
    """
    Interactive plot for fixing DOFs by clicking nodes.
    Left click toggles X DOF, right click toggles Y DOF for a node.
    """
    dof_per_node = 2
    total_dofs = len(all_nodes_coords) * dof_per_node
    fixed_dofs = set()
    coord_to_node_map = {coord: idx + 1 for idx, coord in enumerate(all_nodes_coords)}

    # Helper to plot truss and DOF status
    def plot_interactive(ax, nodes, connectivity, fixed_x, fixed_y):
        ax.clear()
        # Draw truss
        for n1, n2 in connectivity:
            x1, y1 = nodes[n1-1]
            x2, y2 = nodes[n2-1]
            ax.plot([x1, x2], [y1, y2], color='gray', lw=2)
        # Draw nodes
        for i, (x, y) in enumerate(nodes):
            marker = 'o'
            color = 'blue'
            if fixed_x[i] and fixed_y[i]:
                color = 'red'
                marker = 's'
            elif fixed_x[i]:
                color = 'orange'
                marker = '>'
            elif fixed_y[i]:
                color = 'green'
                marker = '^'
            ax.plot(x, y, marker=marker, color=color, markersize=12)
            ax.text(x, y, f"{i+1}", fontsize=10, ha='center', va='center', bbox=dict(facecolor='white', alpha=0.6, edgecolor='black'))
        ax.set_title("Click nodes to fix DOFs: Left=X, Right=Y. Red=Both, Orange=X, Green=Y")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.grid(True)
        ax.set_aspect('equal')

    # Initial DOF status
    fixed_x = [False]*len(all_nodes_coords)
    fixed_y = [False]*len(all_nodes_coords)
    connectivity = []
    for el in stress_nodes:
        node1 = coord_to_node_map[(el[0], el[1])]
        node2 = coord_to_node_map[(el[2], el[3])]
        connectivity.append((node1, node2))

    fig, ax = plt.subplots(figsize=(10, 6))
    plot_interactive(ax, all_nodes_coords, connectivity, fixed_x, fixed_y)

    # Click event handler
    def on_click(event):
        if event.inaxes != ax:
            return
        # Find closest node
        min_dist = float('inf')
        idx = None
        for i, (x, y) in enumerate(all_nodes_coords):
            d = (event.xdata-x)**2 + (event.ydata-y)**2
            if d < min_dist:
                min_dist = d
                idx = i
        if min_dist > 30:  # Only if close to a node
            return
        if event.button == 1:  # Left click toggles X DOF
            fixed_x[idx] = not fixed_x[idx]
        elif event.button == 3:  # Right click toggles Y DOF
            fixed_y[idx] = not fixed_y[idx]
        plot_interactive(ax, all_nodes_coords, connectivity, fixed_x, fixed_y)
        fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', on_click)

    # Add submit button
    submit_ax = plt.axes([0.8, 0.01, 0.15, 0.06])
    submit_btn = Button(submit_ax, 'Submit')
    submit_done = {'clicked': False}
    def submit(event):
        submit_done['clicked'] = True
        plt.close(fig)
    submit_btn.on_clicked(submit)

    plt.show()
    if not submit_done['clicked']:
        raise SystemExit("Boundary conditions window was closed without submitting.")

    # Build fixed_dofs from selections
    for i in range(len(all_nodes_coords)):
        if fixed_x[i]:
            fixed_dofs.add(i*dof_per_node)
        if fixed_y[i]:
            fixed_dofs.add(i*dof_per_node+1)

    fixed_dofs = sorted(list(fixed_dofs))

    # Determine the indices of the free degrees of freedom
    free_dofs = np.setdiff1d(np.arange(total_dofs), fixed_dofs)
    K_reduced = K_global_matrix[np.ix_(free_dofs, free_dofs)]
    F_reduced = F_global_vector[free_dofs]
    Q = np.zeros(total_dofs)
    Q_reduced = np.linalg.solve(K_reduced, F_reduced)
    Q[free_dofs] = Q_reduced

    print("\n📈 Displacement Vector (Q):")
    print(Q)
    for i in range(len(all_nodes_coords)):
        print(f"Node {i+1}: Ux = {Q[i*2]:.6f}, Uy = {Q[i*2+1]:.6f}")

    # Creating connection of nodes for each element
    connectivity = []
    for el in stress_nodes:
        node1 = coord_to_node_map[(el[0], el[1])]
        node2 = coord_to_node_map[(el[2], el[3])]
        connectivity.append((node1, node2))

    print("\nCoordinate to Node Mapping:")
    for coord, node in coord_to_node_map.items():
        print(f"  Node {node}: {coord}")

    print("\nElement Connectivity (Node Pairs):")
    for i, (n1, n2) in enumerate(connectivity, 1):
        print(f"  Element {i}: Node {n1} ↔ Node {n2}")

    return Q, coord_to_node_map, connectivity

# Apply boundary conditions, solve for displacements, and get connectivity
Q, coord_to_node, connectivity = apply_boundary_conditions_by_node(K_global, F_global, all_nodes)

def calculate_element_results(elements_data, displacement_vector, coordinate_to_node_map):
    """
    Calculates the axial stress and force in each truss element.

    Args:
        elements_data: A list of tuples, where each tuple contains
                       (x1, y1, x2, y2, A, E, L) for an element.
        displacement_vector: The global displacement vector (Q).
        coordinate_to_node_map: A dictionary mapping coordinates to node numbers.
    """
    print("\nAxial Stresses and Forces in Each Member:")

    # Iterate through each element to calculate stress and force
    for i, (x1, y1, x2, y2, area, youngs_modulus, _) in enumerate(elements_data):
        # Get node numbers (1-based) from the coordinate map
        node1_num = coordinate_to_node_map[(x1, y1)]
        node2_num = coordinate_to_node_map[(x2, y2)]

        # Get the global degrees of freedom indices for the start and end nodes
        dof1_indices = [(node1_num - 1) * 2, (node1_num - 1) * 2 + 1]
        dof2_indices = [(node2_num - 1) * 2, (node2_num - 1) * 2 + 1]

        # Get the displacement components (Ux, Uy) for the start and end nodes from the global displacement vector
        u1, v1 = displacement_vector[dof1_indices[0]], displacement_vector[dof1_indices[1]]
        u2, v2 = displacement_vector[dof2_indices[0]], displacement_vector[dof2_indices[1]]

        # Recalculate element length (although it's in elements_data, recalculating here for clarity
        # and to ensure consistency with current node positions if needed, though not strictly necessary
        # as we're using original geometry for calculations based on original stiffness)
        length = np.hypot(x2 - x1, y2 - y1)

        # Calculate the direction cosines of the element
        # Avoid division by zero if length is zero
        if length == 0:
            cos_theta = 0
            sin_theta = 0
        else:
            cos_theta = (x2 - x1) / length
            sin_theta = (y2 - y1) / length

        # Calculate the axial strain in the element
        # Strain = (Change in length) / Original length
        # Change in length = (u2-u1)*cos(theta) + (v2-v1)*sin(theta)
        strain = ((u2 - u1) * cos_theta + (v2 - v1) * sin_theta) / length
        # Calculate the axial stress using Hooke's Law (Stress = E * Strain)
        stress = youngs_modulus * strain
        # Calculate the axial force in the element (Force = Stress * Area)
        force = stress * area

        # Print the calculated stress and force for the element
        print(f"  Element {i+1} (Node {node1_num} ↔ Node {node2_num}): Stress = {stress:.6f} Pa, Force = {force:.3f} N")

# Calculate and print the axial stresses and forces in each member
calculate_element_results(elements, Q, coord_to_node)

#==============================================================================
# POSTPROCESSING - VISUALIZATION AND RESULTS
#==============================================================================

# Function to plot truss from node coordinates and connectivity with stress coloring
def plot_truss(nodes_coords, connectivity_list, elements_data=None, displacement_vector=None, linestyle='-', label=None):
    """
    Plots a truss structure with stress-based coloring.

    Args:
        nodes_coords: A list of node coordinates [(x1, y1), (x2, y2), ...].
        connectivity_list: A list of tuples representing element connectivity (node pairs, 1-based).
        elements_data: List of element data for stress calculation.
        displacement_vector: Global displacement vector for stress calculation.
        linestyle: The linestyle for plotting the truss members.
        label: The label for the plot legend.
    """
    # Create figure and axis objects with a single subplot
    fig = plt.gcf()  # Get current figure or create new one
    ax = plt.gca()  # Get current axis or create new one
    
    if elements_data is not None and displacement_vector is not None:
        # Calculate stresses for color mapping
        stresses = []
        for i, (x1, y1, x2, y2, area, youngs_modulus, length) in enumerate(elements_data):
            n1, n2 = connectivity_list[i]
            # Get displacements
            u1, v1 = displacement_vector[(n1-1)*2:(n1-1)*2+2]
            u2, v2 = displacement_vector[(n2-1)*2:(n2-1)*2+2]
            
            # Calculate direction cosines
            cos_theta = (x2 - x1) / length
            sin_theta = (y2 - y1) / length
            
            # Calculate strain and stress
            strain = ((u2 - u1) * cos_theta + (v2 - v1) * sin_theta) / length
            stress = youngs_modulus * strain
            stresses.append(stress)  # Using actual stress for coloring (not absolute)
        
        # Create color map
        norm = plt.Normalize(vmin=min(stresses), vmax=max(stresses))
        cmap = plt.cm.coolwarm
        
        # Plot elements with stress coloring
        for i, (n1, n2) in enumerate(connectivity_list):
            x1, y1 = nodes_coords[n1 - 1]
            x2, y2 = nodes_coords[n2 - 1]
            color = cmap(norm(stresses[i]))
            line = ax.plot([x1, x2], [y1, y2], color=color, linestyle=linestyle, lw=2,
                          label=label if i == 0 else None)
        
        # Add colorbar with proper axes
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Axial Stress (Pa)\nRed = Tension, Blue = Compression')
        
        # Add node markers
        for i, (x, y) in enumerate(nodes_coords):
            ax.plot(x, y, 'ko', markersize=6)
            ax.text(x, y, f'{i+1}', fontsize=8, ha='right', va='bottom')
    else:
        # Default plotting without stress coloring
        for i, (n1, n2) in enumerate(connectivity_list):
            x1, y1 = nodes_coords[n1 - 1]
            x2, y2 = nodes_coords[n2 - 1]
            ax.plot([x1, x2], [y1, y2], color='b', linestyle=linestyle, lw=2,
                   label=label if i == 0 else None)
            
        # Add node markers
        for i, (x, y) in enumerate(nodes_coords):
            ax.plot(x, y, 'ko', markersize=6)
            ax.text(x, y, f'{i+1}', fontsize=8, ha='right', va='bottom')

# Function to compute deflected/scaled positions of nodes
def scaled_displacements(displacement_vector, original_nodes_coords, scale_factor=10):
    """
    Calculates the deflected positions of nodes based on calculated displacements and a scale factor.

    Args:
        displacement_vector: The global displacement vector (Q).
        original_nodes_coords: A list of the original node coordinates.
        scale_factor: The factor by which to scale the displacements for visualization.

    Returns:
        A list of the new (deflected) node coordinates.
    """
    new_positions = []
    # Iterate through each node to calculate its new position
    for i in range(len(original_nodes_coords)):
        # Get original coordinates
        x_orig, y_orig = original_nodes_coords[i]
        # Get displacements (Ux, Uy) for the node from the global displacement vector
        # Adjusting index for 0-based list and 2 DOFs per node
        disp_x = displacement_vector[i * 2]
        disp_y = displacement_vector[i * 2 + 1]

        # Calculate new coordinates by adding scaled displacements to original coordinates
        new_x = x_orig + scale_factor * disp_x
        new_y = y_orig + scale_factor * disp_y
        new_positions.append((new_x, new_y))
    return new_positions

# Calculate the scaled deflected node positions
# The scale_factor was 10 in the original code, keeping it the same.
# The label in the plot title should reflect the actual scale factor used.
scale_factor_used = 10
scaled_nodes = scaled_displacements(Q, all_nodes, scale_factor=scale_factor_used)

# Combine all X and Y values from both original and scaled nodes to dynamically set axis limits
all_x = [x for x, _ in all_nodes + scaled_nodes]
all_y = [y for _, y in all_nodes + scaled_nodes]
# Add a margin around the plotted points for better visualization
x_margin = (max(all_x) - min(all_x)) * 0.1 or 1
y_margin = (max(all_y) - min(all_y)) * 0.1 or 1

# Create new figure with two subplots
deformed_fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
plt.subplots_adjust(wspace=0.3)  # Add some space between subplots

# Plot original truss
plt.sca(ax1)  # Set current axes to first subplot
plot_truss(all_nodes, connectivity, linestyle='-', label='Original')
ax1.set_title("Original Truss Structure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.grid(True)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)
ax1.set_ylim(min(all_y) - y_margin, max(all_y) + y_margin)

# Plot deformed truss with stress coloring
plt.sca(ax2)  # Set current axes to second subplot
plot_truss(scaled_nodes, connectivity, elements_data=elements, displacement_vector=Q, 
           linestyle='-', label=f'Deflected (x{scale_factor_used})')
ax2.set_title(f"Deformed Truss Structure\n(Displacement scaled by {scale_factor_used})")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.grid(True)
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)
ax2.set_ylim(min(all_y) - y_margin, max(all_y) + y_margin)

# Add overall title
deformed_fig.suptitle("Truss Analysis Results", fontsize=14, y=1.02)

# Display the deformed plot and wait for it to close
plt.show()

# After deformed plot is closed, close the initial plot
plt.close(initial_fig)

print(all_nodes)
