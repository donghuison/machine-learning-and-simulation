module const

   implicit none


   character(*),parameter :: input_dir = "data_40x40/"
   character(*),parameter :: output_dir = "data_40x40/"
   ! character(*),parameter :: input_dir = "data_40x40_serial/"
   ! character(*),parameter :: output_dir = "data_40x40_serial/"

   real(8),parameter :: pi= acos(-1.0d0)
   real(8),parameter :: pi2=2.0d0*pi

   !   time control parameters
   integer,parameter :: nstop = 100  ! Number of time steps
   real(8),parameter :: tend  = 10.0d0  ! Total simulation time
   real(8),parameter :: dtout = 0.1d0
   real(8),parameter :: safety=0.5d0
   real(8),parameter :: dtmin=1.d-10

   ! CELL & MPI
   integer,parameter :: margin=1

   ! -------------------------------------------------------------------
   ! integer,parameter :: mpisize_x=1, mpisize_y=1 !=> 1 procs : (40x40)
   integer,parameter :: mpisize_x=2, mpisize_y=2 !=> 4 procs : (40x40)
   ! integer,parameter :: mpisize_x=15, mpisize_y=15 !=> 225 procs : (300x300)
   ! -------------------------------------------------------------------

   ! -------------------------------------------------------------------
   ! integer,parameter :: domain_x = 41, domain_y = 41  ! Match Python's 41x41 grid
   integer,parameter :: domain_x = 40, domain_y = 40
   ! integer,parameter :: domain_x = 300, domain_y = 300
   ! -------------------------------------------------------------------

   integer,parameter :: ix = int(domain_x/mpisize_x)+2*margin,jx=int(domain_y/mpisize_y)+2*margin
   integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
   integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
   logical,parameter :: pbcheck(2) = (/.false., .false./)  ! Non-Periodic in both directions

   ! size (uniform grid)
   real(8),parameter :: xmin = 0.0d0, ymin = 0.0d0
   real(8),parameter :: xmax = 1.0d0, ymax = 1.0d0

   real(8),parameter :: dxg0 = (xmax-xmin)/dble(igx-margin*2)
   real(8),parameter :: dyg0 = (ymax-ymin)/dble(jgx-margin*2)

   ! model parameter - Stable Fluids specific
   integer,parameter :: nvar = 5

   real(8),parameter :: g_gamma = 1.0d0
   real(8),parameter :: ro0 = 1.0d0, pr0 = 0.0d0
   real(8),parameter :: vx0 = 0.0d0, vy0 = 0.0d0
   real(8),parameter :: dt0 = 0.1d0 ! Stable fluids uses larger timestep
   real(8),parameter :: KINEMATIC_VISCOSITY = 0.0001d0  ! Lower viscosity for stable fluids
   integer,parameter :: MAX_ITER_CG = 5000  ! Max iterations for conjugate gradient

   ! Forcing parameters
   real(8),parameter :: force_x_min = 0.4d0, force_x_max = 0.6d0
   real(8),parameter :: force_y_min = 0.1d0, force_y_max = 0.3d0
   real(8),parameter :: force_magnitude = 1.0d0
   real(8),parameter :: force_time_decay = 0.5d0

end module const

