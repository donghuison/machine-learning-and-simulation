module advection
   implicit none
   private

   public :: advect_velocity

contains

   subroutine advect_velocity(nvar, margin, ix, jx, x, y, dx, dy, dt, Vc, Vci, xmin, xmax, ymin, ymax, g_gamma)
      use bnd
      ! use mpi_setup, only : mpid, mnull

      implicit none

      integer, intent(in) :: nvar, margin, ix, jx
      real(8), intent(in) :: dt
      real(8), intent(in) :: xmin, xmax, ymin, ymax
      real(8), dimension(ix), intent(in) :: x, dx
      real(8), dimension(jx), intent(in) :: y, dy
      real(8), dimension(ix,jx,nvar), intent(inout) :: Vc
      real(8), dimension(ix,jx,nvar), intent(inout) :: Vci
      real(8), intent(in) :: g_gamma

      integer :: i, j
      real(8) :: x_back, y_back
      real(8) :: x_depart, y_depart
      real(8) :: u_interp, v_interp
      real(8), dimension(ix,jx) :: u_new, v_new
      real(8) :: domain_min_x, domain_max_x, domain_min_y, domain_max_y


      ! Define domain boundaries
      domain_min_x = xmin
      domain_max_x = xmax
      domain_min_y = ymin
      domain_max_y = ymax

      ! Initialize with current velocities
      do j=1, jx
         do i=1, ix
            u_new(i,j) = Vc(i,j,3)
            v_new(i,j) = Vc(i,j,4)
         end do
      end do

      ! Semi-Lagrangian advection for all interior points
      do j=2, jx-1
         do i=2, ix-1
            ! Skip boundary points where velocity should remain zero
            ! if (i == 1 .or. i == ix .or. j == 1 .or. j == jx) then
            !    u_new(i,j) = 0.0d0
            !    v_new(i,j) = 0.0d0
            !    cycle
            ! end if

            ! Backtrace position
            x_depart = x(i) - dt*Vc(i,j,3)
            y_depart = y(j) - dt*Vc(i,j,4)

            ! Clip to domain boundaries
            x_back = max(domain_min_x, min(domain_max_x, x_depart))
            y_back = max(domain_min_y, min(domain_max_y, y_depart))

            ! Bilinear interpolation
            call bilinear_interp(ix, jx, x, y, Vc(:,:,3), x_back, y_back, u_interp)
            call bilinear_interp(ix, jx, x, y, Vc(:,:,4), x_back, y_back, v_interp)

            u_new(i,j) = u_interp
            v_new(i,j) = v_interp
         end do
      end do

      ! Update velocities

      do j=2, jx-1
         do i=2, ix-1
            Vc(i,j,3) = u_new(i,j)
            Vc(i,j,4) = v_new(i,j)
         end do
      end do

      call bnd__exec_vel(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)

   end subroutine advect_velocity


   subroutine bilinear_interp(ix, jx, x, y, field, x_back, y_back, result)
      implicit none

      integer, intent(in) :: ix, jx
      real(8), dimension(ix), intent(in) :: x
      real(8), dimension(jx), intent(in) :: y
      real(8), dimension(ix,jx), intent(in) :: field
      real(8), intent(in) :: x_back, y_back
      real(8), intent(out) :: result

      integer :: i, j, i0, j0, i1, j1
      real(8) :: wx, wy, dx_local, dy_local

      ! Find grid cell containing the point
      ! Handle edge cases where xi or yi might equal the last grid point
      i0 = 1
      do i=1, ix-1
         if (x_back >= x(i) .and. x_back <= x(i+1)) then
            i0 = i
            exit
         end if
      end do

      j0 = 1
      do j=1, jx-1
         if (y_back >= y(j) .and. y_back <= y(j+1)) then
            j0 = j
            exit
         end if
      end do

      ! Ensure we're within bounds for interpolation
      i0 = max(1, min(ix-1, i0))
      j0 = max(1, min(jx-1, j0))
      i1 = min(i0+1, ix)
      j1 = min(j0+1, jx)

      ! Handle special cases at boundaries
      if (i0 == ix) then
         i0 = ix - 1
         i1 = ix
      end if
      if (j0 == jx) then
         j0 = jx - 1
         j1 = jx
      end if

      ! Compute local grid spacing
      dx_local = x(i1) - x(i0)
      dy_local = y(j1) - y(j0)

      ! Compute weights with safety checks
      if (dx_local > 1.0d-10) then
         wx = (x_back - x(i0)) / dx_local
      else
         wx = 0.0d0
      end if

      if (dy_local > 1.0d-10) then
         wy = (y_back - y(j0)) / dy_local
      else
         wy = 0.0d0
      end if

      ! Clamp weights to [0,1]
      wx = max(0.0d0, min(1.0d0, wx))
      wy = max(0.0d0, min(1.0d0, wy))

      ! Bilinear interpolation
      result = (1.0d0-wx) * (1.0d0-wy) * field(i0,j0) + &
         wx * (1.0d0-wy) * field(i1,j0) + &
         (1.0d0-wx) * wy * field(i0,j1) + &
         wx * wy * field(i1,j1)

   end subroutine bilinear_interp


end module advection
