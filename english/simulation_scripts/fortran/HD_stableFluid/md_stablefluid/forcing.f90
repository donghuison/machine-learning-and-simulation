module forcing
   use const

   implicit none
   private

   public :: apply_forcing

contains

   subroutine apply_forcing(nvar, margin, ix, jx, time, x, y, Vc, dt)
      implicit none

      integer, intent(in) :: nvar, margin, ix, jx
      real(8), intent(in) :: time, dt
      real(8), dimension(ix), intent(in) :: x
      real(8), dimension(jx), intent(in) :: y
      real(8), dimension(ix,jx,nvar), intent(inout) :: Vc

      integer :: i, j
      real(8) :: time_decay, fx, fy

      ! Calculate time decay factor
      time_decay = max(2.0d0-force_time_decay*time, 0.0d0)

      ! Apply forcing in the specified region
      do j=2, jx-1
         do i=2, ix-1
            ! do j = 1, jx
            !    do i = 1, ix
            if (x(i) > force_x_min .and. x(i) < force_x_max .and. &
               y(j) > force_y_min .and. y(j) < force_y_max) then
               ! Apply upward force
               fx = 0.0d0
               fy = force_magnitude * time_decay

               ! Update velocities with forcing
               Vc(i,j,3) = Vc(i,j,3) + dt*fx  ! vx
               Vc(i,j,4) = Vc(i,j,4) + dt*fy  ! vy
            end if
         end do
      end do

   end subroutine apply_forcing

end module forcing
