module flux_calc

   implicit none
   private

   public :: flux_calc__vtent

contains

   subroutine flux_calc__vtent(nvar,margin,ix,jx,Vc,dx,dy,flux,KINEMATIC_VISCOSITY)
      use operator

      implicit none
      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),intent(in) :: KINEMATIC_VISCOSITY
      real(8),dimension(ix,jx,nvar),intent(in) :: Vc
      real(8),dimension(ix),intent(in) :: dx
      real(8),dimension(jx),intent(in) :: dy
      real(8),dimension(ix,jx,nvar),intent(out) :: flux

      integer :: i,j
      real(8),dimension(ix,jx) :: dvx_dx, dvx_dy
      real(8),dimension(ix,jx) :: dvy_dx, dvy_dy
      real(8),dimension(ix,jx) :: lap_vx, lap_vy

      ! Compute derivatives
      call central_difference_x(margin,ix,jx,Vc(:,:,3),dvx_dx,dx)
      call central_difference_y(margin,ix,jx,Vc(:,:,3),dvx_dy,dy)
      call central_difference_x(margin,ix,jx,Vc(:,:,4),dvy_dx,dx)
      call central_difference_y(margin,ix,jx,Vc(:,:,4),dvy_dy,dy)

      ! Compute laplacian
      call laplacian(margin,ix,jx,Vc(:,:,3),lap_vx,dx,dy)
      call laplacian(margin,ix,jx,Vc(:,:,4),lap_vy,dx,dy)

      flux = 0.0d0
      do j=2, jx-1
         do i=2, ix-1
            flux(i,j,3) = -(Vc(i,j,3)*dvx_dx(i,j) + Vc(i,j,4)*dvx_dy(i,j)) + KINEMATIC_VISCOSITY*lap_vx(i,j)

            flux(i,j,4) = -(Vc(i,j,3)*dvy_dx(i,j) + Vc(i,j,4)*dvy_dy(i,j)) + KINEMATIC_VISCOSITY*lap_vy(i,j)
         end do
      end do

   end subroutine flux_calc__vtent

end module flux_calc
