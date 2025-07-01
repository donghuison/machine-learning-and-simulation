module const

   implicit none

   character(*),parameter :: input_dir = "data_40x40/"
   character(*),parameter :: output_dir = "data_40x40/"
   ! character(*),parameter :: input_dir = "data_300x300/"
   ! character(*),parameter :: output_dir = "data_300x300/"
   ! character(*),parameter :: input_dir = "data_400x400/"
   ! character(*),parameter :: output_dir = "data_400x400/"

   real(8),parameter :: pi= acos(-1.0d0)
   real(8),parameter :: pi2=2.0d0*pi

   !   time control parameters
   integer,parameter :: nstop = int(2.d7)
   ! integer,parameter :: nstop = 500
   real(8),parameter :: tend  = 0.5d0
   real(8),parameter :: dtout = 0.01d0
   real(8),parameter :: safety=0.5d0
   real(8),parameter :: dtmin=1.d-10

   ! CELL & MPI
   integer,parameter :: margin=1

   ! -------------------------------------------------------------------
   integer,parameter :: mpisize_x=2, mpisize_y=2 !=> 4 procs : (40x40)
   ! integer,parameter :: mpisize_x=15, mpisize_y=15 !=> 225 procs : (300x300)
   ! -------------------------------------------------------------------

   ! -------------------------------------------------------------------
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

   ! model parameter - Lid driven cavity specific
   integer,parameter :: nvar = 5

   real(8),parameter :: g_gamma = 1.0d0
   real(8),parameter :: ro0 = 1.0d0, pr0 = 0.0d0
   real(8),parameter :: vx0 = 0.0d0, vy0 = 0.0d0
   real(8),parameter :: dt0 = 0.001d0 ! 40x40
   ! real(8),parameter :: dt0 = 0.00001d0  ! Reduced for 300x300 grid stability
   ! real(8),parameter :: dt0 = 0.00001d0  ! Reduced for 400x400 grid stability
   real(8),parameter :: KINEMATIC_VISCOSITY = 0.1d0
   integer,parameter :: N_PRESSURE_POISSON_ITERATIONS = 50 ! 40x40
   ! integer,parameter :: N_PRESSURE_POISSON_ITERATIONS = 500 ! 300x300, 400x400

end module const

