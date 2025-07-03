"""
Quick comparison of Python stable fluids results with Fortran visualization
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Set dark background
plt.style.use("dark_background")

# Load the Fortran result
fortran_img_path = '/Users/donghuison/workspace/myGit/KHU-STUDY/Project-01/machine-learning-and-simulation/english/simulation_scripts/fortran/HD_stableFluid/python/2D-analysis/fig_41x41_curl/evolution_41x41.png'
fortran_img = Image.open(fortran_img_path)

# Load a Python snapshot if available
python_img_path = '/Users/donghuison/workspace/myGit/KHU-STUDY/Project-01/machine-learning-and-simulation/english/simulation_scripts/stable_fluids_frames/frame_0099.png'

try:
    python_img = Image.open(python_img_path)
    
    # Create side-by-side comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    ax1.imshow(python_img)
    ax1.set_title('Python Implementation (Final Frame)')
    ax1.axis('off')
    
    ax2.imshow(fortran_img)
    ax2.set_title('Fortran Implementation (Evolution)')
    ax2.axis('off')
    
    plt.suptitle('Stable Fluids Comparison: Python vs Fortran', fontsize=16)
    plt.tight_layout()
    plt.savefig('comparison_python_fortran.png', dpi=150, bbox_inches='tight')
    plt.close()
    
except FileNotFoundError:
    print("Python frames not found. Showing only Fortran results.")
    plt.figure(figsize=(10, 8))
    plt.imshow(fortran_img)
    plt.title('Fortran Implementation - Stable Fluids Evolution')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('fortran_results_only.png', dpi=150, bbox_inches='tight')
    plt.close()

print("Comparison complete!")
print("\nKey improvements made to match Python implementation:")
print("1. Grid size: 40x40 → 41x41")
print("2. Grid spacing: Fixed formula to match Python")
print("3. Grid origin: Corrected to start at (0.0, 0.0)")
print("4. Diffusion solver: Properly implemented implicit scheme")
print("\nRemaining differences:")
print("- Linear solver: Jacobi (Fortran) vs CG (Python)")
print("- Interpolation: Custom (Fortran) vs scipy.interpn (Python)")
print("- Boundary conditions: Ghost cells vs implicit")