module init
   use const
   implicit none
   private

   public :: initialize

   real(8),public,dimension(ix) :: x,dx
   real(8),public,dimension(jx) :: y,dy

   real(8),public,dimension(ix,jx,nvar) :: Uc, Vc
   real(8),public,dimension(ix,jx,nvar) :: Uci, Vci

   integer,public :: ns,nd
   real(8),public :: min_dx
   real(8),public :: time,timep

   integer,public :: mwflag,mw
   real(8),public :: dt,dtg
   real(8),public :: dtnew

contains

   subroutine initialize
      use model, only : model_setup
      use bnd
      use mpi_setup

      call mpi_setup__init(mpisize_x,mpisize_y,pbcheck)

      call model_setup(Uc,Vc,&
         Uci,Vci,&
         x,dx,y,dy,min_dx)

      call bnd_exec(nvar,margin,ix,jx,Vc,g_gamma,x,y,Vci)
      !   call convert__ptoc(margin,ix,jx,g_gamma,Vc,Uc)

      nd = 0
      time = 0.0d0
      ns = 0
      merr  = 0
      dt = dt0

   end subroutine initialize

end module init

