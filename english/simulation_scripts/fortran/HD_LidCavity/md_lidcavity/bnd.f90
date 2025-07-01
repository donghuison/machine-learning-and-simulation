module bnd

   implicit none
   private

   public :: bnd_exec, bnd__exec_vel, bnd__exec_prs

contains

   subroutine bnd_exec(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)

      use mpi_setup, only : mpid, mnull
      use boundary

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix),intent(in) :: x
      real(8),dimension(jx),intent(in) :: y
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc,Vci

      call boundary__mpi(nvar,margin,ix,jx,Vc)

   end subroutine bnd_exec

   subroutine bnd__exec_vel(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)

      use mpi_setup, only : mpid, mnull
      use boundary

      integer,intent(in) :: nvar
      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix),intent(in) :: x
      real(8),dimension(jx),intent(in) :: y
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix,jx,nvar),intent(inout) :: Vc,Vci

      call boundary__mpi(nvar,margin,ix,jx,Vc)

      if(mpid%l == mnull) then
         call bd_consx(0,margin,0.0d0,Vc(:,:,3),ix,jx)
         call bd_consx(0,margin,0.0d0,Vc(:,:,4),ix,jx)
      end if

      if(mpid%r == mnull) then
         call bd_consx(1,margin,0.0d0,Vc(:,:,3),ix,jx)
         call bd_consx(1,margin,0.0d0,Vc(:,:,4),ix,jx)
      end if

      if(mpid%b == mnull)then
         call bd_consy(0,margin,0.0d0,Vc(:,:,3),ix,jx)
         call bd_consy(0,margin,0.0d0,Vc(:,:,4),ix,jx)
      end if

      if(mpid%f == mnull) then
         call bd_consy(1,margin,1.0d0,Vc(:,:,3),ix,jx)
         call bd_consy(1,margin,0.0d0,Vc(:,:,4),ix,jx)
      end if

   end subroutine bnd__exec_vel

   subroutine bnd__exec_prs(margin,ix,jx,pr,g_gamma,x,y)
      use mpi_setup, only : mpid, mnull
      use boundary

      integer,intent(in) :: margin,ix,jx
      real(8),dimension(ix),intent(in) :: x
      real(8),dimension(jx),intent(in) :: y
      real(8),intent(in) :: g_gamma
      real(8),dimension(ix,jx),intent(inout) :: pr

      call boundary__mpi_qq1(margin,ix,jx,pr)

      if(mpid%l == mnull) then
         call bd_synpx_car(0,margin,pr,ix,jx)
      end if

      if(mpid%r == mnull) then
         call bd_synpx_car(1,margin,pr,ix,jx)
      end if

      if(mpid%b == mnull)then
         call bd_synpy_car(0,margin,pr,ix,jx)
      end if

      if(mpid%f == mnull) then
         call bd_consy(1,margin,0.0d0,pr,ix,jx)
      end if

   end subroutine bnd__exec_prs


end module bnd











