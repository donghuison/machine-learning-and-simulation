from dac_read import dac_read
import numpy as np
import matplotlib.pyplot as plt
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        
def main():
    
    # Create dictionary of base directories for each scheme
    base_dirs = {
        'LID_Cavity': '/Users/donghuison/workspace/myGit/KHU-STUDY/Project-01/machine-learning-and-simulation/english/simulation_scripts/fortran/HD_LidCavity/'
    }
    
    # Select current scheme (예: WENO5-Z+M 사용)
    current_scheme = 'FDM'

    base_dir = base_dirs['LID_Cavity']
    
    data_type = '40x40'
    # data_type = '300x300'
    # data_type = '400x400'

    if data_type == '40x40':
        data_path = base_dir + 'md_lidcavity/data_40x40/*_{}_rank=*.dac'
        iskip=2
        jskip=2
    elif data_type == '300x300':
        data_path = base_dir + 'md_lidcavity/data_300x300/*_{}_rank=*.dac'
        iskip=10
        jskip=10
    elif data_type == '400x400':
        data_path = base_dir + 'md_lidcavity/data_400x400/*_{}_rank=*.dac'
        iskip=10
        jskip=10
    else:
        raise ValueError(f"Unknown data_type: {data_type}")
    
    # Read data using the selected pattern
    data_rh, x, y = dac_read(data_path.format('ro'), dimension=2)
    data_pr, x, y = dac_read(data_path.format('pr'), dimension=2)
    data_vx, x, y = dac_read(data_path.format('vx'), dimension=2)
    data_vy, x, y = dac_read(data_path.format('vy'), dimension=2)
    

    logging.info(f"Reading data for data_type: {data_type}")
    logging.info("Data reading completed.")

    # time_range = [0, 10, 20, 30, 40, 50]
    time_range = [0, 1, 2, 3, 4, 5]

    # Initialize global min and max for the color scale
    global_min = float('inf')
    global_max = float('-inf')
    
    logging.info(f"Calculating global min and max for data_type: {data_type}")
    for t in time_range:
        pr = data_pr[t,:,:]
        # pr = np.log10(pr)
        global_min = min(global_min, pr.min())
        global_max = max(global_max, pr.max())


    xax2d, yax2d = np.meshgrid(x,y)
    extent = [xax2d.min(), xax2d.max(), yax2d.min(), yax2d.max()]
    
    # Ensure output directory exists
    output_dir = base_dir + f'python/2D-analysis/fig_{data_type}/'
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Starting plotting loop.")
    for t in time_range:
        logging.info(f"Plotting time step: {t}")
        fig, ax = plt.subplots(figsize=(8, 6))
        # c = ax.pcolormesh(x, y, np.log10(data_rh[nt0,:,:]), cmap='jet', shading='auto')
        # c = ax.imshow(data_pr[t,:,:], cmap='jet', origin='lower', extent=extent, vmin=global_min, vmax=global_max)
        c = ax.contourf(x,y, data_pr[t,:,:], cmap='coolwarm', origin='lower', vmin=global_min, vmax=global_max)
        ax.set_xlim(xax2d.min(), xax2d.max())
        ax.set_ylim(yax2d.min(), yax2d.max())
        ax.set_xlabel(r'$X$')
        ax.set_ylabel(r'$Y$')
        # ax.set_title(f'$\mathit{{Lid\\ Cavity\\ Flow\\ at\\ }} t={t*1e-2}$')
        ax.grid(False)

        plt.gca().set_aspect('equal', adjustable='box') 
        plt.colorbar(c, ax=ax, orientation='horizontal', label=r'$P$')

        # 화살표 개수 조절을 위한 샘플링 (데이터 타입에 따라 이미 정의된 iskip, jskip 사용)
        plt.quiver(x[::iskip], y[::jskip], 
                   data_vx[t, ::jskip, ::iskip], 
                   data_vy[t, ::jskip, ::iskip], 
                   color="black")
        
        # Save the plot
        output_filename = os.path.join(output_dir, f'pr_vel_{t}_{data_type}_{current_scheme}.pdf')
        plt.savefig(output_filename, format='pdf', bbox_inches='tight', dpi=600)
        plt.close()
    
        
    logging.info("Plotting loop completed.")

if __name__ == "__main__":
    main()