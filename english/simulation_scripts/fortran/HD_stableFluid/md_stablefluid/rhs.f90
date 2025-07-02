module RightHandSide
   implicit none

   public :: rhs_calc__vtent

contains

   subroutine rhs_calc__vtent(nvar,margin,ix,jx,dt,dx,dy,flux,rhs)
      implicit none

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),intent(in) :: dt
      real(8),dimension(ix),intent(in) :: dx
      real(8),dimension(jx),intent(in) :: dy
      real(8),dimension(ix,jx,nvar),intent(in) :: flux
      real(8),dimension(ix,jx,nvar),intent(out) :: rhs

      integer :: i,j,n

      ! Update tentative velocities
      rhs = 0.0d0
      do n=1,nvar
         do j=2, jx-1
            do i=2, ix-1
               rhs(i,j,n) = dt*flux(i,j,n)
            end do
         end do
      end do

   end subroutine rhs_calc__vtent

end module RightHandSide
