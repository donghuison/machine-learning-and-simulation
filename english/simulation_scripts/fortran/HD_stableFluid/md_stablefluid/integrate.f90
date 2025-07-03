module integrate
   implicit none
   private

   public :: integrate__StableFluids

contains

   subroutine integrate__StableFluids(nvar,margin,ix,jx,g_gamma,x,dx,y,dy,dt,time,&
      Vc,Vci,KINEMATIC_VISCOSITY,MAX_ITER_CG,xmin,xmax,ymin,ymax)

      use forcing
      use advection
      use operator
      use bnd
      use mpi_setup

      implicit none

      integer,intent(in) :: nvar, margin, ix, jx, MAX_ITER_CG
      real(8),intent(in) :: g_gamma, dt, time, KINEMATIC_VISCOSITY
      real(8),intent(in) :: xmin, xmax, ymin, ymax
      real(8),dimension(ix),intent(in) :: x, dx
      real(8),dimension(jx),intent(in) :: y, dy
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc, Vci

      ! Step 1: Apply forces
      call apply_forcing(nvar, margin, ix, jx, time, x, y, Vc, dt)

      ! Step 2: Advect velocities (self-advection)
      call advect_velocity(nvar, margin, ix, jx, x, y, dx, dy, dt, Vc, Vci, xmin, xmax, ymin, ymax, g_gamma)

      ! Step 3: Diffuse velocities implicitly
      call diffuse_velocity(nvar, margin, ix, jx, dx, dy, dt, Vc, Vci, &
         KINEMATIC_VISCOSITY, MAX_ITER_CG, g_gamma, x, y)

      ! Step 4.1: Compute pressure correction
      ! Step 4.2: Project velocities to be divergence-free
      call project_velocity(nvar, margin, ix, jx, dx, dy, dt, Vc, Vci, &
         MAX_ITER_CG, g_gamma, x, y)

   end subroutine integrate__StableFluids


   subroutine diffuse_velocity(nvar, margin, ix, jx, dx, dy, dt, Vc, Vci, &
      KINEMATIC_VISCOSITY, MAX_ITER_CG, g_gamma, x, y)
      use operator
      use bnd
      use mpi_setup

      implicit none

      integer, intent(in) :: nvar, margin, ix, jx, MAX_ITER_CG
      real(8), intent(in) :: dt, KINEMATIC_VISCOSITY, g_gamma
      real(8), dimension(ix), intent(in) :: x, dx
      real(8), dimension(jx), intent(in) :: y, dy
      real(8), dimension(ix,jx,nvar), intent(inout) :: Vc, Vci

      integer :: i, j, iter
      real(8) :: alpha, beta, error, error_global
      real(8), dimension(ix,jx) :: u_old, v_old, u_new, v_new
      real(8), dimension(ix,jx) :: lap_u, lap_v
      real(8), parameter :: tolerance = 1.0d-6

      ! Store current velocities as right-hand side
      u_old = Vc(:,:,3)
      v_old = Vc(:,:,4)

      ! Initialize with current velocities
      u_new = u_old
      v_new = v_old

      ! Solve (I - nu*dt*Laplacian)u_new = u_old using Jacobi iteration
      alpha = KINEMATIC_VISCOSITY * dt
      error_global = 100.0d0

      do iter = 1, MAX_ITER_CG
         error = 0.0d0

         ! Update interior points using Jacobi iteration
         ! Solving: u_new - alpha * Laplacian(u_new) = u_old
         ! Jacobi: u_new(i,j) = (u_old(i,j) + alpha*(neighbors)) / (1 + 4*alpha/h^2)
         do j = 2, jx-1
            do i = 2, ix-1
               ! Coefficients
               beta = 1.0d0 / (1.0d0 + 2.0d0*alpha/(dx(i)*dx(i)) + 2.0d0*alpha/(dy(j)*dy(j)))

               ! Jacobi update using neighbors from previous iteration
               u_new(i,j) = beta * (u_old(i,j) + &
                  alpha/(dx(i)*dx(i)) * (Vc(i+1,j,3) + Vc(i-1,j,3)) + &
                  alpha/(dy(j)*dy(j)) * (Vc(i,j+1,3) + Vc(i,j-1,3)))

               v_new(i,j) = beta * (v_old(i,j) + &
                  alpha/(dx(i)*dx(i)) * (Vc(i+1,j,4) + Vc(i-1,j,4)) + &
                  alpha/(dy(j)*dy(j)) * (Vc(i,j+1,4) + Vc(i,j-1,4)))

               error = max(error, abs(Vc(i,j,3) - u_new(i,j)))
               error = max(error, abs(Vc(i,j,4) - v_new(i,j)))
            end do
         end do

         ! Update velocities
         do j=2, jx-1
            do i=2, ix-1
               Vc(i,j,3) = u_new(i,j)
               Vc(i,j,4) = v_new(i,j)
            end do
         end do

         ! Apply boundary conditions
         call bnd__exec_vel(nvar, margin, ix, jx, Vc, g_gamma, x, y, Vci)

         ! Check convergence
         call mpi_allreduce(error, error_global, 1, mdp, mmax, mcomw, merr)
         if (error_global < tolerance) exit
      end do

      if(mpid%mpirank == 0) then
         print *, "Diffuse velocity iterations:", iter
         print *, "Diffuse velocity error:", error_global
      end if

   end subroutine diffuse_velocity


   subroutine project_velocity(nvar, margin, ix, jx, dx, dy, dt, Vc, Vci, &
      MAX_ITER_CG, g_gamma, x, y)
      use operator
      use bnd
      use mpi_setup

      implicit none

      integer, intent(in) :: nvar, margin, ix, jx, MAX_ITER_CG
      real(8), intent(in) :: dt, g_gamma
      real(8), dimension(ix), intent(in) :: x, dx
      real(8), dimension(jx), intent(in) :: y, dy
      real(8), dimension(ix,jx,nvar), intent(inout) :: Vc, Vci

      integer :: i, j, iter
      real(8) :: rhs, coeff, error, error_global
      real(8), dimension(ix,jx) :: div, pr_old, pr_new, lap_p
      real(8), dimension(ix,jx) :: dp_dx, dp_dy
      real(8), parameter :: tolerance = 1.0d-6

      error_global = 100.0d0

      ! Compute divergence of velocity field
      call divergence(nvar, margin, ix, jx, Vc, dx, dy, div)

      ! Initialize pressure
      pr_old = Vc(:,:,2)
      pr_new = pr_old

      ! Solve Poisson equation for pressure: Laplacian(p) = div(u)
      do iter = 1, MAX_ITER_CG
         call laplacian(margin, ix, jx, pr_new, lap_p, dx, dy)

         error = 0.0d0

         ! Jacobi iteration
         do j = 2, jx-1
            do i = 2, ix-1
               coeff = 1.0d0 / (2.0d0/(dx(i)*dx(i)) + 2.0d0/(dy(j)*dy(j)))
               rhs = div(i,j)

               pr_new(i,j) = coeff * ( &
                  (pr_old(i+1,j) + pr_old(i-1,j))/(dx(i)*dx(i)) + &
                  (pr_old(i,j+1) + pr_old(i,j-1))/(dy(j)*dy(j)) - rhs)

               error = max(error, abs(pr_new(i,j) - pr_old(i,j)))
            end do
         end do

         ! Apply pressure boundary conditions
         call bnd__exec_prs(margin, ix, jx, pr_new, g_gamma, x, y)

         pr_old = pr_new

         ! Check convergence
         call mpi_allreduce(error, error_global, 1, mdp, mmax, mcomw, merr)
         if (error_global < tolerance) exit
      end do

      if(mpid%mpirank == 0) then
         print *, "Project velocity iterations:", iter
         print *, "Project velocity error:", error_global
      end if

      ! Update pressure in solution vector
      do j = 2, jx-1
         do i = 2, ix-1
            Vc(i,j,2) = pr_new(i,j)
         end do
      end do

      call bnd__exec_prs(margin, ix, jx, Vc(:,:,2), g_gamma, x, y)

      ! Compute pressure gradient
      call central_difference_x(margin, ix, jx, Vc(:,:,2), dp_dx, dx)
      call central_difference_y(margin, ix, jx, Vc(:,:,2), dp_dy, dy)

      ! Correct velocities to be divergence-free
      do j = 2, jx-1
         do i = 2, ix-1
            Vc(i,j,3) = Vc(i,j,3) - dp_dx(i,j)
            Vc(i,j,4) = Vc(i,j,4) - dp_dy(i,j)
         end do
      end do

      ! Apply boundary conditions
      call bnd__exec_vel(nvar, margin, ix, jx, Vc, g_gamma, x, y, Vci)

   end subroutine project_velocity

end module integrate
