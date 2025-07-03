"""
Visualization script for Stable Fluids simulation
Matches the visual style of stable_fluids_python_simple.py
Focus on curl/vorticity visualization with dark background
"""

from dac_read import dac_read
import numpy as np
import matplotlib.pyplot as plt
import os
import logging

# Try to import cmasher for redshift colormap, fall back to alternatives
try:
    import cmasher as cmr
    has_cmasher = True
except ImportError:
    has_cmasher = False
    logging.warning("cmasher not installed. Using alternative colormap.")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def compute_curl_2d(vx, vy, x, y):
    """
    Compute 2D curl (vorticity) of velocity field
    curl = ∂vy/∂x - ∂vx/∂y
    """
    ny, nx = vx.shape
    curl = np.zeros_like(vx)
    
    # Compute grid spacing
    dx = x[1] - x[0] if len(x) > 1 else 1.0
    dy = y[1] - y[0] if len(y) > 1 else 1.0
    
    # Central differences for interior points
    curl[1:-1, 1:-1] = (vy[1:-1, 2:] - vy[1:-1, :-2]) / (2 * dx) - \
                       (vx[2:, 1:-1] - vx[:-2, 1:-1]) / (2 * dy)
    
    # Boundary values (one-sided differences)
    # Left and right boundaries
    curl[:, 0] = (vy[:, 1] - vy[:, 0]) / dx - (vx[1:, 0] - vx[:-1, 0]).mean() / dy
    curl[:, -1] = (vy[:, -1] - vy[:, -2]) / dx - (vx[1:, -1] - vx[:-1, -1]).mean() / dy
    
    # Top and bottom boundaries
    curl[0, :] = (vy[0, 1:] - vy[0, :-1]).mean() / dx - (vx[1, :] - vx[0, :]) / dy
    curl[-1, :] = (vy[-1, 1:] - vy[-1, :-1]).mean() / dx - (vx[-1, :] - vx[-2, :]) / dy
    
    return curl


def plot_stable_fluids_snapshot(x, y, vx, vy, time_value, output_path, data_type='40x40'):
    """
    Plot a snapshot matching the style of stable_fluids_python_simple.py
    """
    # Set dark background style
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(5, 5), dpi=160)
    ax = plt.gca()
    
    # Create meshgrid (note: Fortran uses different indexing than Python)
    X, Y = np.meshgrid(x, y)
    
    # Compute curl
    curl = compute_curl_2d(vx, vy, x, y)
    
    # Plot curl using contourf
    if has_cmasher:
        cmap = cmr.redshift
    else:
        # Use a red-based colormap as alternative
        cmap = 'hot'
    
    # contour = ax.contourf(X, Y, curl, cmap=cmap, levels=100)
    contour = ax.contourf(X, Y, curl, cmap=cmap, levels=100)
    
    # Determine skip for quiver based on grid size
    if data_type == '40x40' or data_type == '41x41':
        skip = 1  # Show more vectors for small grid
    elif data_type == '300x300':
        skip = 8
    else:
        skip = 10
    
    # Add velocity vectors
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
              vx[::skip, ::skip], vy[::skip, ::skip],
              color="dimgray", alpha=0.8)
    
    # Set limits and labels
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Stable Fluids - Time: {time_value:.1f}')
    ax.set_aspect('equal')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=160, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Saved: {output_path}")


def create_interactive_plot(x, y, data_vx, data_vy, data_pr, output_dir, data_type='40x40'):
    """
    Create an interactive plot that updates in real-time (similar to the Python example)
    Saves snapshots at regular intervals
    """
    plt.style.use("dark_background")
    plt.ion()  # Turn on interactive mode
    
    fig = plt.figure(figsize=(5, 5), dpi=160)
    ax = plt.gca()
    
    # Create meshgrid
    X, Y = np.meshgrid(x, y)

    # Determine skip for quiver
    if data_type == '40x40' or data_type == '41x41':
        skip = 1
    elif data_type == '300x300':
        skip = 8
    else:
        skip = 10
    
    n_timesteps = data_vx.shape[0]
    # n_timesteps = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    
    for t in range(n_timesteps):
        plt.clf()  # Clear the figure
        
        # Compute curl
        curl = compute_curl_2d(data_vx[t], data_vy[t], x, y)
        
        # Plot
        if has_cmasher:
            plt.contourf(X, Y, curl, cmap=cmr.redshift, levels=100)
        else:
            plt.contourf(X, Y, curl, cmap='hot', levels=100)
        
        plt.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                   data_vx[t, ::skip, ::skip], data_vy[t, ::skip, ::skip],
                   color="dimgray")
        
        plt.xlim(X.min(), X.max())
        plt.ylim(Y.min(), Y.max())
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f'Time: {t * 1.0:.1f}')
        plt.gca().set_aspect('equal')
        
        # Save every 10th frame
        if t % 10 == 0:
            output_path = os.path.join(output_dir, f'interactive_frame_{t:04d}.png')
            plt.savefig(output_path, dpi=160, bbox_inches='tight')
        
        plt.draw()
        plt.pause(0.0001)
    
    plt.ioff()
    plt.show()


def main():
    # Configuration
    base_dir = '/Users/donghuison/workspace/myGit/KHU-STUDY/Project-01/machine-learning-and-simulation/english/simulation_scripts/fortran/HD_stableFluid/'
    # data_type = '41x41'  # Options: '40x40', '41x41', '300x300', '400x400'
    data_type = '40x40'  # Options: '40x40', '41x41', '300x300', '400x400'
    
    # Data paths
    data_path = base_dir + f'md_stablefluid/data_{data_type}/*_{{}}_rank=*.dac'
    
    # Read data
    logging.info(f"Reading data for {data_type}")
    data_rh, x, y = dac_read(data_path.format('ro'), dimension=2)
    data_pr, x, y = dac_read(data_path.format('pr'), dimension=2)
    data_vx, x, y = dac_read(data_path.format('vx'), dimension=2)
    data_vy, x, y = dac_read(data_path.format('vy'), dimension=2)
    
    # Output directory
    output_dir = base_dir + f'python/2D-analysis/fig_{data_type}_curl/'
    os.makedirs(output_dir, exist_ok=True)
    
    # Select time steps to plot
    n_timesteps = data_vx.shape[0]
    logging.info(f"Total timesteps: {n_timesteps}")
    
    # Plot snapshots at regular intervals
    snapshot_interval = max(1, n_timesteps // 10)
    
    for t in range(0, n_timesteps, snapshot_interval):
        output_path = os.path.join(output_dir, f'curl_snapshot_{t:04d}_{data_type}.png')
        plot_stable_fluids_snapshot(x, y, data_vx[t], data_vy[t], 
                                  t * 1.0, output_path, data_type)
    
    # Create a comparison figure showing evolution
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    plt.style.use("dark_background")
    fig.suptitle('Stable Fluids Evolution - Vorticity', fontsize=16)
    
    # Select 6 time points
    time_points = np.linspace(0, n_timesteps-1, 6, dtype=int)
    
    for idx, (ax, t) in enumerate(zip(axes.flat, time_points)):
        X, Y = np.meshgrid(x, y)
        curl = compute_curl_2d(data_vx[t], data_vy[t], x, y)
        
        if has_cmasher:
            contour = ax.contourf(X, Y, curl, cmap=cmr.redshift, levels=100)
        else:
            contour = ax.contourf(X, Y, curl, cmap='hot', levels=50)
        
        # Fewer quiver arrows for clarity
        skip = 3 if data_type == '40x40' or data_type == '40x40_serial' else 15
        ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                  data_vx[t, ::skip, ::skip], data_vy[t, ::skip, ::skip],
                  color="dimgray", alpha=0.7, scale=20)
        
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())
        ax.set_title(f'Time = {t * 1.0:.1f}')
        ax.set_aspect('equal')
        
        if idx >= 3:
            ax.set_xlabel('X')
        if idx % 3 == 0:
            ax.set_ylabel('Y')
    
    plt.tight_layout()
    evolution_path = os.path.join(output_dir, f'evolution_{data_type}.png')
    plt.savefig(evolution_path, dpi=200, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Saved evolution plot: {evolution_path}")
    logging.info("Visualization complete!")
    
    # Optional: Create interactive plot (comment out if not needed)
    # logging.info("Starting interactive visualization...")
    # create_interactive_plot(x, y, data_vx, data_vy, data_pr, output_dir, data_type)


if __name__ == "__main__":
    main()