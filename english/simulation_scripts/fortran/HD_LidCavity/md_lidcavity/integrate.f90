module integrate

   implicit none
   private

   public :: integrate__EULER1


contains

   subroutine integrate__EULER1(nvar,margin,ix,jx,g_gamma,x,dx,y,dy,dt,&
      Vc,Vci,&
      KINEMATIC_VISCOSITY,&
      N_PRESSURE_POISSON_ITERATIONS)


      integer,intent(in) :: nvar
      integer,intent(in) :: ix,jx,margin
      integer,intent(in) :: N_PRESSURE_POISSON_ITERATIONS
      real(8),intent(in) :: KINEMATIC_VISCOSITY
      real(8),intent(in) :: g_gamma
      real(8),intent(in) :: dt
      real(8),dimension(ix),intent(in) :: x,dx
      real(8),dimension(jx),intent(in) :: y,dy
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vci
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc

      call calc_momentum_eq(nvar,margin,ix,jx,Vc,dx,dy,dt,Vci,KINEMATIC_VISCOSITY,g_gamma,x,y)
      call calc_pressure_eq(nvar,margin,ix,jx,Vc,dx,dy,dt,N_PRESSURE_POISSON_ITERATIONS,g_gamma,x,y)
      call calc_correct_vel(nvar,margin,ix,jx,Vc,dx,dy,dt,Vci,g_gamma,x,y)


   end subroutine integrate__EULER1

   subroutine calc_momentum_eq(nvar,margin,ix,jx,Vc,dx,dy,dt,Vci,KINEMATIC_VISCOSITY,g_gamma,x,y)
      use flux_calc
      use RightHandSide
      use bnd

      implicit none

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),intent(in) :: dt
      real(8),intent(in) :: KINEMATIC_VISCOSITY
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix),intent(in) :: x,dx
      real(8),dimension(jx),intent(in) :: y,dy
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc,Vci

      integer :: i,j,n
      real(8),dimension(ix,jx,nvar) :: flux
      real(8),dimension(ix,jx,nvar) :: rhs

      call flux_calc__vtent(nvar,margin,ix,jx,Vc,dx,dy,flux,KINEMATIC_VISCOSITY)
      call rhs_calc__vtent(nvar,margin,ix,jx,dt,dx,dy,flux,rhs)

      do n=1,nvar
         do j=2, jx-1
            do i=2, ix-1
               Vc(i,j,n) = Vc(i,j,n) + rhs(i,j,n)
            end do
         end do
      end do

      call bnd__exec_vel(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)

   end subroutine calc_momentum_eq


   subroutine calc_pressure_eq(nvar,margin,ix,jx,Vc,dx,dy,dt,N_PRESSURE_POISSON_ITERATIONS,g_gamma,x,y)
      use bnd
      use operator
      use mpi_setup

      implicit none

      integer,intent(in) :: nvar
      integer,intent(in) :: N_PRESSURE_POISSON_ITERATIONS
      integer,intent(in) :: margin,ix,jx
      real(8),intent(in) :: dt
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix),intent(in) :: x,dx
      real(8),dimension(jx),intent(in) :: y,dy
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc

      integer :: iter,i,j
      real(8),parameter :: tolerance = 1.0d-6
      real(8),dimension(ix,jx) :: div
      real(8),dimension(ix,jx) :: pr_old,pr_new
      real(8) :: coeff,rhs,dl2
      real(8) :: dA,dA_global

      dA_global = 100.0d0
      iter = 0

      ! Compute divergence of tentative velocity
      call divergence(nvar,margin,ix,jx,Vc,dx,dy,div)

      ! Initialize pressure
      do j=1, jx
         do i=1, ix
            pr_old(i,j) = Vc(i,j,2)
         end do
      end do

      ! Jacobi iterations
      ! do iter=1, N_PRESSURE_POISSON_ITERATIONS
      do while(dA_global > tolerance .and. iter <= N_PRESSURE_POISSON_ITERATIONS)

         do j=2, jx-1
            do i=2, ix-1
               dl2 = dx(i)*dy(j)
               coeff = (dl2*Vc(i,j,1))/dt
               rhs = coeff*div(i,j)
               pr_new(i,j) = 0.25d0*(pr_old(i+1,j)+pr_old(i-1,j)+pr_old(i,j+1)+pr_old(i,j-1) - rhs)
            end do
         end do

         ! Apply pressure boundary conditions
         call bnd__exec_prs(margin,ix,jx,pr_new,g_gamma,x,y)

         ! update pressure
         dA=0.0d0
         do j=1, jx
            do i=1, ix
               dA = max(abs(pr_new(i,j) - pr_old(i,j)), dA)
               pr_old(i,j) = pr_new(i,j)
            end do
         end do
         call mpi_allreduce(dA,dA_global,1,mdp,mmax,mcomw,merr)
         iter = iter + 1

      end do

      do j=1, jx
         do i=1, ix
            Vc(i,j,2) = pr_new(i,j)
         end do
      end do

      if(mpid%mpirank == 0) then
         print*, "dA_global: ", dA_global
         print*, "iter: ", iter
      end if

   end subroutine calc_pressure_eq


   subroutine calc_correct_vel(nvar,margin,ix,jx,Vc,dx,dy,dt,Vci,g_gamma,x,y)
      use operator
      use bnd

      implicit none

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),intent(in) :: dt
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix),intent(in) :: x,dx
      real(8),dimension(jx),intent(in) :: y,dy
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc,Vci

      integer :: i,j
      real(8),dimension(ix,jx) :: dp_dx,dp_dy
      real(8) :: inv_ro,coeff

      dp_dx = 0.0d0
      dp_dy = 0.0d0

      ! Compute pressure gradient
      call central_difference_x(margin,ix,jx,Vc(:,:,2), dp_dx, dx)
      call central_difference_y(margin,ix,jx,Vc(:,:,2), dp_dy, dx)

      ! Correct velocities
      do j=2, jx-1
         do i=2, ix-1
            inv_ro = 1.d0 / Vc(i,j,1)
            coeff = dt*inv_ro
            Vc(i,j,3) = Vc(i,j,3) - coeff*dp_dx(i,j)
            Vc(i,j,4) = Vc(i,j,4) - coeff*dp_dy(i,j)
         end do
      end do

      ! Apply boundary conditions
      call bnd__exec_vel(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)

   end subroutine calc_correct_vel

end module integrate
