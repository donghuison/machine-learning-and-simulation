"""
Visualize the difference between node-centered (Python) and cell-centered (Fortran) grids
used in the Stable Fluids implementations.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def visualize_grid_comparison():
    """Create a visual comparison of the two grid types."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Python: Node-centered grid
    ax1.set_title("Python: Node-Centered Grid\n(41x41 nodes)", fontsize=14, weight='bold')
    ax1.set_xlim(-0.05, 0.3)
    ax1.set_ylim(-0.05, 0.3)
    
    # Draw nodes
    node_x = np.linspace(0.0, 0.25, 11)  # Show first 11 nodes (0.0, 0.025, ..., 0.25)
    node_y = np.linspace(0.0, 0.25, 11)
    
    for i, x in enumerate(node_x):
        for j, y in enumerate(node_y):
            # Draw nodes
            ax1.plot(x, y, 'ro', markersize=6)
            
            # Draw grid lines
            if i < len(node_x) - 1:
                ax1.plot([x, node_x[i+1]], [y, y], 'k-', alpha=0.3, linewidth=0.5)
            if j < len(node_y) - 1:
                ax1.plot([x, x], [y, node_y[j+1]], 'k-', alpha=0.3, linewidth=0.5)
    
    # Highlight boundary
    rect1 = patches.Rectangle((0, 0), 0.25, 0.25, linewidth=2, 
                             edgecolor='blue', facecolor='none')
    ax1.add_patch(rect1)
    
    # Add annotations
    ax1.text(0, -0.03, '(0,0)', ha='center', va='top', fontsize=10)
    ax1.text(0.025, -0.03, '0.025', ha='center', va='top', fontsize=8)
    ax1.text(0.25, -0.03, '0.25', ha='center', va='top', fontsize=10)
    ax1.text(-0.03, 0, '(0,0)', ha='right', va='center', fontsize=10)
    ax1.text(-0.03, 0.25, '0.25', ha='right', va='center', fontsize=10)
    
    ax1.set_xlabel("x", fontsize=12)
    ax1.set_ylabel("y", fontsize=12)
    ax1.grid(False)
    ax1.set_aspect('equal')
    
    # Add text explanations
    ax1.text(0.125, -0.08, "Nodes at boundaries\nDerivatives = 0 at edges", 
             ha='center', va='top', fontsize=9, style='italic')
    
    # Fortran: Cell-centered grid
    ax2.set_title("Fortran: Cell-Centered Grid\n(42x42 cells with ghost cells)", 
                  fontsize=14, weight='bold')
    ax2.set_xlim(-0.05, 0.3)
    ax2.set_ylim(-0.05, 0.3)
    
    # Cell edges (including ghost cell)
    cell_edges_x = np.array([-0.0125, 0.0125, 0.0375, 0.0625, 0.0875, 0.1125, 
                            0.1375, 0.1625, 0.1875, 0.2125, 0.2375, 0.2625])
    cell_edges_y = cell_edges_x.copy()
    
    # Draw cell edges
    for x in cell_edges_x:
        ax2.axvline(x, color='k', alpha=0.3, linewidth=0.5)
    for y in cell_edges_y:
        ax2.axhline(y, color='k', alpha=0.3, linewidth=0.5)
    
    # Draw cell centers
    cell_centers_x = (cell_edges_x[:-1] + cell_edges_x[1:]) / 2
    cell_centers_y = cell_centers_x.copy()
    
    for i, x in enumerate(cell_centers_x):
        for j, y in enumerate(cell_centers_y):
            if i == 0 or j == 0:  # Ghost cells
                ax2.plot(x, y, 'g^', markersize=6, label='Ghost cell' if i == 0 and j == 0 else '')
            else:
                ax2.plot(x, y, 'bs', markersize=6, label='Cell center' if i == 1 and j == 1 else '')
    
    # Highlight physical domain (excluding ghost cells)
    rect2 = patches.Rectangle((0.0125, 0.0125), 0.25, 0.25, linewidth=2, 
                             edgecolor='blue', facecolor='none')
    ax2.add_patch(rect2)
    
    # Highlight ghost cell region
    ghost_rect = patches.Rectangle((-0.0125, -0.0125), 0.275, 0.275, linewidth=1.5, 
                                  edgecolor='green', facecolor='none', linestyle='--')
    ax2.add_patch(ghost_rect)
    
    # Add annotations
    ax2.text(0, -0.03, '0.0', ha='center', va='top', fontsize=10)
    ax2.text(-0.0125, -0.03, '-0.0125', ha='center', va='top', fontsize=8)
    ax2.text(0.2625, -0.03, '0.2625', ha='center', va='top', fontsize=10)
    
    ax2.set_xlabel("x", fontsize=12)
    ax2.set_ylabel("y", fontsize=12)
    ax2.grid(False)
    ax2.set_aspect('equal')
    
    # Add legend
    ax2.legend(loc='upper right', fontsize=8)
    
    # Add text explanations
    ax2.text(0.125, -0.08, "Cell values at centers\nGhost cells for BCs", 
             ha='center', va='top', fontsize=9, style='italic')
    
    plt.suptitle("Grid Structure Comparison: Node-Centered vs Cell-Centered", 
                 fontsize=16, weight='bold')
    plt.tight_layout()
    
    # Save the figure
    plt.savefig('grid_comparison.png', dpi=300, bbox_inches='tight')
    print("Grid comparison visualization saved as 'grid_comparison.png'")
    
    # Create a second figure showing interpolation differences
    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Show interpolation example
    ax3.set_title("Python: Interpolation on Nodes", fontsize=14, weight='bold')
    ax3.set_xlim(0, 0.1)
    ax3.set_ylim(0, 0.1)
    
    # Draw a few nodes
    nodes = [(0, 0), (0.025, 0), (0.05, 0), (0.075, 0), (0.1, 0),
             (0, 0.025), (0.025, 0.025), (0.05, 0.025), (0.075, 0.025), (0.1, 0.025),
             (0, 0.05), (0.025, 0.05), (0.05, 0.05), (0.075, 0.05), (0.1, 0.05)]
    
    for node in nodes:
        ax3.plot(node[0], node[1], 'ro', markersize=8)
    
    # Show an interpolation point
    interp_point = (0.04, 0.03)
    ax3.plot(interp_point[0], interp_point[1], 'k*', markersize=12)
    
    # Highlight interpolation stencil
    stencil_nodes = [(0.025, 0.025), (0.05, 0.025), (0.025, 0.05), (0.05, 0.05)]
    for node in stencil_nodes:
        ax3.plot(node[0], node[1], 'ro', markersize=12, markeredgecolor='blue', 
                markeredgewidth=2)
    
    # Draw interpolation box
    interp_box = patches.Rectangle((0.025, 0.025), 0.025, 0.025, linewidth=2, 
                                  edgecolor='blue', facecolor='blue', alpha=0.1)
    ax3.add_patch(interp_box)
    
    ax3.text(0.04, 0.02, 'Interpolation\npoint', ha='center', va='top', fontsize=8)
    ax3.set_xlabel("x", fontsize=12)
    ax3.set_ylabel("y", fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect('equal')
    
    # Fortran interpolation
    ax4.set_title("Fortran: Interpolation on Cells", fontsize=14, weight='bold')
    ax4.set_xlim(0, 0.1)
    ax4.set_ylim(0, 0.1)
    
    # Draw cell boundaries
    for x in [0, 0.025, 0.05, 0.075, 0.1]:
        ax4.axvline(x, color='k', alpha=0.3)
    for y in [0, 0.025, 0.05, 0.075, 0.1]:
        ax4.axhline(y, color='k', alpha=0.3)
    
    # Draw cell centers
    cell_centers = [(0.0125, 0.0125), (0.0375, 0.0125), (0.0625, 0.0125), (0.0875, 0.0125),
                   (0.0125, 0.0375), (0.0375, 0.0375), (0.0625, 0.0375), (0.0875, 0.0375),
                   (0.0125, 0.0625), (0.0375, 0.0625), (0.0625, 0.0625), (0.0875, 0.0625)]
    
    for center in cell_centers:
        ax4.plot(center[0], center[1], 'bs', markersize=8)
    
    # Show interpolation point
    ax4.plot(interp_point[0], interp_point[1], 'k*', markersize=12)
    
    # Highlight interpolation stencil
    stencil_centers = [(0.0125, 0.0125), (0.0375, 0.0125), (0.0125, 0.0375), (0.0375, 0.0375)]
    for center in stencil_centers:
        ax4.plot(center[0], center[1], 'bs', markersize=12, markeredgecolor='red', 
                markeredgewidth=2)
    
    # Draw interpolation region
    interp_region = patches.Rectangle((0, 0), 0.05, 0.05, linewidth=2, 
                                     edgecolor='red', facecolor='red', alpha=0.1)
    ax4.add_patch(interp_region)
    
    ax4.text(0.04, 0.02, 'Interpolation\npoint', ha='center', va='top', fontsize=8)
    ax4.set_xlabel("x", fontsize=12)
    ax4.set_ylabel("y", fontsize=12)
    ax4.grid(False)
    ax4.set_aspect('equal')
    
    plt.suptitle("Interpolation Scheme Comparison", fontsize=16, weight='bold')
    plt.tight_layout()
    
    plt.savefig('interpolation_comparison.png', dpi=300, bbox_inches='tight')
    print("Interpolation comparison visualization saved as 'interpolation_comparison.png'")
    
    plt.show()

if __name__ == "__main__":
    visualize_grid_comparison()