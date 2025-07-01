module operator

   implicit none
   private
   public :: central_difference_x, central_difference_y, laplacian, divergence

contains

   subroutine central_difference_x(margin,ix,jx,f,df_dx,dx)
      implicit none

      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix,jx),intent(in) :: f
      real(8),dimension(ix),intent(in) :: dx
      real(8),dimension(ix,jx),intent(out) :: df_dx

      integer :: i, j
      real(8) :: inv_dl

      df_dx = 0.0d0
      do j=2, jx-1
         do i=2, ix-1
            inv_dl = 1.d0 / (2.d0*dx(i))
            df_dx(i,j) = (f(i+1,j) - f(i-1,j))*inv_dl
         end do
      end do

   end subroutine central_difference_x

   subroutine central_difference_y(margin,ix,jx,f,df_dy,dy)
      implicit none
      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix,jx),intent(in) :: f
      real(8),dimension(jx),intent(in) :: dy
      real(8),dimension(ix,jx),intent(out) :: df_dy

      integer :: i, j
      real(8) :: inv_dl

      df_dy = 0.0d0
      do j=2, jx-1
         do i=2, ix-1
            inv_dl = 1.d0 / (2.d0*dy(j))
            df_dy(i,j) = (f(i,j+1) - f(i,j-1))*inv_dl
         end do
      end do

   end subroutine central_difference_y

   subroutine laplacian(margin,ix,jx,f,lap_f,dx,dy)
      implicit none

      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix,jx),intent(in) :: f
      real(8),dimension(ix),intent(in) :: dx
      real(8),dimension(jx),intent(in) :: dy
      real(8),dimension(ix,jx),intent(out) :: lap_f

      integer :: i,j
      real(8) :: inv_dx2, inv_dy2
      real(8) :: lap_fx, lap_fy

      lap_f = 0.0d0
      do j=2, jx-1
         do i=2, ix-1
            inv_dx2 = 1.0d0 / (dx(i)*dx(i))
            inv_dy2 = 1.0d0 / (dy(j)*dy(j))

            lap_fx = (f(i+1,j) - 2.d0*f(i,j) + f(i-1,j))*inv_dx2
            lap_fy = (f(i,j+1) - 2.d0*f(i,j) + f(i,j-1))*inv_dy2
            lap_f(i,j) = lap_fx + lap_fy
            !    lap_f(i, j) = (f(i+1,j) + f(i-1,j) + f(i,j+1) + f(i,j-1) - 4.0d0*f(i, j))*inv_dl2
         end do
      end do

   end subroutine laplacian

   subroutine divergence(nvar,margin,ix,jx,Vc,dx,dy,div)
      implicit none

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix),intent(in) :: dx
      real(8),dimension(jx),intent(in) :: dy
      real(8),dimension(ix,jx,nvar),intent(in) :: Vc
      real(8),dimension(ix,jx),intent(out) :: div

      integer :: i,j
      real(8),dimension(ix,jx) :: dvx_dx,dvy_dy

      call central_difference_x(margin,ix,jx,Vc(:,:,3), dvx_dx, dx)
      call central_difference_y(margin,ix,jx,Vc(:,:,4), dvy_dy, dy)

      div = 0.0d0
      do j=2, jx-1
         do i=2, ix-1
            div(i,j) = dvx_dx(i,j) + dvy_dy(i,j)
         end do
      end do

   end subroutine divergence

end module operator
