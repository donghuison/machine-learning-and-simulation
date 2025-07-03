"""
Improved visualization for Stable Fluids simulation
Includes vorticity visualization and animation capabilities
Based on the stable_fluids_python_simple.py visualization approach
"""

from dac_read import dac_read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter, FFMpegWriter
import logging
import os
import glob

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def compute_curl_2d(vx, vy, dx, dy):
    """
    Compute 2D curl (vorticity) of velocity field
    curl = ∂vy/∂x - ∂vx/∂y
    """
    ny, nx = vx.shape
    curl = np.zeros_like(vx)
    
    # Central differences for interior points
    curl[1:-1, 1:-1] = (vy[1:-1, 2:] - vy[1:-1, :-2]) / (2 * dx) - \
                       (vx[2:, 1:-1] - vx[:-2, 1:-1]) / (2 * dy)
    
    return curl


def compute_velocity_magnitude(vx, vy):
    """Compute velocity magnitude from components"""
    return np.sqrt(vx**2 + vy**2)


def plot_stable_fluids_frame(x, y, pr, vx, vy, time_step, output_path=None, 
                           show_curl=True, show_pressure=False, show_streamlines=False,
                           data_type='40x40', colormap='RdBu_r'):
    """
    Plot a single frame of the stable fluids simulation
    Similar to the visualization in stable_fluids_python_simple.py
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
    
    # Compute grid spacing (assuming uniform grid)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    # Create meshgrid for plotting
    X, Y = np.meshgrid(x, y)
    
    if show_curl:
        # Compute and plot vorticity (curl)
        curl = compute_curl_2d(vx, vy, dx, dy)
        
        # Use similar colormap to the Python example (redshift-like)
        levels = 100
        vmin, vmax = np.percentile(curl, [1, 99])  # Use percentiles for better contrast
        
        contour = ax.contourf(X, Y, curl, cmap=colormap, levels=levels,
                             vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.02)
        cbar.set_label('Vorticity (curl)', fontsize=12)
        
    elif show_pressure:
        # Plot pressure field
        levels = 50
        contour = ax.contourf(X, Y, pr, cmap='coolwarm', levels=levels)
        cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.02)
        cbar.set_label('Pressure', fontsize=12)
    
    # Determine skip for quiver plot based on data type
    if data_type == '40x40':
        skip = 2
    elif data_type == '300x300':
        skip = 10
    elif data_type == '400x400':
        skip = 12
    else:
        skip = 4
    
    # Add velocity vectors (quiver plot)
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              vx[::skip, ::skip], vy[::skip, ::skip],
              color='dimgray', alpha=0.6, scale=None, width=0.003)
    
    # Add streamlines if requested
    if show_streamlines:
        speed = compute_velocity_magnitude(vx, vy)
        strm = ax.streamplot(x, y, vx.T, vy.T, density=1.5, 
                           color=speed.T, cmap='plasma', linewidth=1)
    
    # Add forcing region rectangle
    force_region = plt.Rectangle((0.4, 0.1), 0.2, 0.2, 
                               fill=False, edgecolor='yellow', 
                               linewidth=2, linestyle='--', alpha=0.5)
    ax.add_patch(force_region)
    ax.text(0.5, 0.32, 'Forcing Region', ha='center', color='yellow', fontsize=10)
    
    # Set plot properties
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title(f'Stable Fluids Simulation - Time: {time_step:.2f}', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
    else:
        return fig, ax


def create_animation(x, y, data_pr, data_vx, data_vy, time_steps, 
                    output_path, data_type='40x40', fps=10):
    """
    Create an animation of the stable fluids simulation
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
    
    # Compute grid spacing
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    # Create meshgrid
    X, Y = np.meshgrid(x, y)
    
    # Initialize plot elements
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    # Determine skip for quiver
    if data_type == '40x40':
        skip = 2
    elif data_type == '300x300':
        skip = 10
    else:
        skip = 4
    
    # Initial plot
    curl_init = compute_curl_2d(data_vx[0], data_vy[0], dx, dy)
    vmin, vmax = np.percentile(curl_init, [1, 99])
    
    contour = ax.contourf(X, Y, curl_init, cmap='RdBu_r', levels=100)
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('Vorticity', fontsize=12)
    
    quiver = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                      data_vx[0, ::skip, ::skip], data_vy[0, ::skip, ::skip],
                      color='dimgray', alpha=0.6, scale=None, width=0.003)
    
    # Add forcing region
    force_region = plt.Rectangle((0.4, 0.1), 0.2, 0.2, 
                               fill=False, edgecolor='yellow', 
                               linewidth=2, linestyle='--', alpha=0.5)
    ax.add_patch(force_region)
    
    title = ax.set_title(f'Stable Fluids - Time: 0.00', fontsize=14)
    
    def animate(frame):
        # Clear previous contours
        for coll in ax.collections[:-1]:  # Keep the last one (rectangle)
            coll.remove()
        
        # Compute new curl
        curl = compute_curl_2d(data_vx[frame], data_vy[frame], dx, dy)
        
        # Update contour plot
        ax.contourf(X, Y, curl, cmap='RdBu_r', levels=100, vmin=vmin, vmax=vmax)
        
        # Update quiver
        quiver.set_UVC(data_vx[frame, ::skip, ::skip], 
                      data_vy[frame, ::skip, ::skip])
        
        # Update title
        title.set_text(f'Stable Fluids - Time: {frame * 1.0:.2f}')
        
        return ax.collections + [quiver, title]
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(time_steps),
                                 interval=1000/fps, blit=False)
    
    # Save animation
    if output_path.endswith('.gif'):
        writer = PillowWriter(fps=fps)
    else:
        writer = FFMpegWriter(fps=fps, metadata=dict(artist='Stable Fluids'), 
                            bitrate=1800)
    
    anim.save(output_path, writer=writer)
    plt.close()
    
    logging.info(f"Animation saved to {output_path}")


def main():
    # Configuration
    base_dir = '/Users/donghuison/workspace/myGit/KHU-STUDY/Project-01/machine-learning-and-simulation/english/simulation_scripts/fortran/HD_stableFluid/'
    data_type = '40x40'  # Options: '40x40', '300x300', '400x400'
    
    # Data paths
    if data_type == '40x40':
        data_path = base_dir + 'md_stablefluid/data_40x40/*_{}_rank=*.dac'
    elif data_type == '300x300':
        data_path = base_dir + 'md_stablefluid/data_300x300/*_{}_rank=*.dac'
    elif data_type == '400x400':
        data_path = base_dir + 'md_stablefluid/data_400x400/*_{}_rank=*.dac'
    else:
        raise ValueError(f"Unknown data_type: {data_type}")
    
    # Read data
    logging.info(f"Reading data for {data_type}")
    data_rh, x, y = dac_read(data_path.format('ro'), dimension=2)
    data_pr, x, y = dac_read(data_path.format('pr'), dimension=2)
    data_vx, x, y = dac_read(data_path.format('vx'), dimension=2)
    data_vy, x, y = dac_read(data_path.format('vy'), dimension=2)
    
    logging.info(f"Data shape: {data_vx.shape}")
    logging.info(f"Grid: x=[{x.min():.3f}, {x.max():.3f}], y=[{y.min():.3f}, {y.max():.3f}]")
    
    # Output directory
    output_dir = base_dir + f'python/2D-analysis/fig_{data_type}_improved/'
    os.makedirs(output_dir, exist_ok=True)
    
    # Time steps to plot
    n_timesteps = data_vx.shape[0]
    if n_timesteps >= 10:
        time_range = list(range(0, n_timesteps, n_timesteps // 10))
    else:
        time_range = list(range(n_timesteps))
    
    logging.info(f"Plotting time steps: {time_range}")
    
    # Plot individual frames with vorticity
    for t in time_range:
        logging.info(f"Plotting time step {t}")
        
        # Vorticity plot
        output_path = os.path.join(output_dir, f'vorticity_{t:04d}_{data_type}.png')
        plot_stable_fluids_frame(x, y, data_pr[t], data_vx[t], data_vy[t], 
                               t * 1.0, output_path, show_curl=True)
        
        # Pressure plot
        output_path = os.path.join(output_dir, f'pressure_{t:04d}_{data_type}.png')
        plot_stable_fluids_frame(x, y, data_pr[t], data_vx[t], data_vy[t], 
                               t * 1.0, output_path, show_curl=False, show_pressure=True)
    
    # Create animations
    logging.info("Creating animations...")
    
    # Vorticity animation
    anim_path = os.path.join(output_dir, f'stable_fluids_vorticity_{data_type}.gif')
    create_animation(x, y, data_pr, data_vx, data_vy, 
                    list(range(n_timesteps)), anim_path, data_type)
    
    # Also save as MP4 if ffmpeg is available
    try:
        anim_path_mp4 = os.path.join(output_dir, f'stable_fluids_vorticity_{data_type}.mp4')
        create_animation(x, y, data_pr, data_vx, data_vy, 
                        list(range(n_timesteps)), anim_path_mp4, data_type)
    except:
        logging.warning("Could not save MP4 animation. Make sure ffmpeg is installed.")
    
    logging.info("Visualization complete!")


if __name__ == "__main__":
    main()