module getNewdt

   use mpi_setup
   implicit none
   private

   public :: getNewdt__hd

contains

   subroutine getNewdt__hd(margin,safety,dtmin,ix,jx,g_gamma,Vc,x,dx,y,dy&
      ,dt,min_dx,KINEMATIC_VISCOSITY)

      integer,intent(in) :: ix,jx,margin
      real(8),intent(in) :: safety
      real(8),intent(in) :: dtmin
      real(8),intent(in) :: min_dx
      real(8),intent(in) :: g_gamma
      real(8),intent(in) :: KINEMATIC_VISCOSITY
      real(8),dimension(ix),intent(in) :: x,dx
      real(8),dimension(jx),intent(in) :: y,dy
      real(8),dimension(ix,jx,13),intent(in) :: Vc
      real(8), intent(inout) :: dt

      integer :: i,j
      real(8) :: dtmax

      do j=1, jx
         do i=1, ix
            dtmax = (0.5d0 * dx(i)*dy(j) / KINEMATIC_VISCOSITY)
            if(dt > safety * dtmax) then
               print *, "Stability is not guarenteed"
               stop
            end if
         end do
      end do

      ! dtg = min(1.1d0*beforedt,safety/dtmaxi)
      ! call mpi_allreduce(dtg,dt,1,mdp,mmin,mcomw,merr)

   end subroutine getNewdt__hd


end module getNewdt

