module model

   use const
   use mpi_setup, only : mpid

   implicit none
   private

   public :: model_setup

contains

   subroutine model_setup(Uc,Vc,&
      Uci,Vci,&
      x,dx,y,dy,min_dx)
      !---Input & Output
      real(8),dimension(ix),intent(out) :: x,dx
      real(8),dimension(jx),intent(out) :: y,dy
      real(8),dimension(ix,jx,nvar),intent(out) :: Uc,Vc
      real(8),dimension(ix,jx,nvar),intent(out) :: Uci,Vci
      real(8),intent(out) :: min_dx

      integer :: i,j
      integer :: ig,jg
      integer :: izero,jzero
      real(8),dimension(0:ix) :: xm
      real(8),dimension(igx)  :: xg,dxg
      real(8),dimension(0:igx):: xmg
      real(8),dimension(0:jx) :: ym
      real(8),dimension(jgx)  :: yg,dyg
      real(8),dimension(0:jgx):: ymg


      Uc  = 0.0d0
      Vc  = 0.0d0
      Uci = 0.0d0
      Vci = 0.0d0

      !---Step 1a.-------------------------------------------------------------|
      ! set global x-grid
      do i=1,igx
         dxg(i)=dxg0
      end do

      ! X origin
      izero=margin-1
      ! izero=igx/2-1
      xmg(izero)= -dxg0*0.5d0
      do i=izero,igx-1
         xmg(i+1) = xmg(i)+dxg(i+1)
      end do
      do i=izero-1,0,-1
         xmg(i) = xmg(i+1)-dxg(i+1)
      end do
      do i=1,igx
         xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
      end do
      xg(1) = 0.5d0*(xmg(1) + (xmg(1)-dxg0))

      !---Step 1b.-------------------------------------------------------------|
      ! set global y-grid
      do j=1,jgx
         dyg(j)=dyg0
      end do

      !  Y origin
      jzero=margin-1
      ! jzero=jgx/2-1
      ymg(jzero)=-dyg0*0.5d0
      do j=jzero,jgx-1
         ymg(j+1) = ymg(j)+dyg(j+1)
      end do
      do j=jzero-1,0,-1
         ymg(j) = ymg(j+1)-dyg(j+1)
      end do

      do j=1,jgx
         yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
      end do
      yg(1) = 0.5d0*(ymg(1) + (ymg(1)-dyg0))

      !---Step 2a.-------------------------------------------------------------|
      ! set individual x-grid
      do i=1,ix
         ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
         x(i) = xg(ig)
         dx(i) = dxg(ig)
      end do
      do i=0,ix
         ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
         xm(i) = xmg(ig)
      end do

      !---Step 2b.-------------------------------------------------------------|
      ! set individual y-grid
      do j=1,jx
         jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
         y(j) = yg(jg)
         dy(j) = dyg(jg)
      enddo
      do j=0,jx
         jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
         ym(j) = ymg(jg)
      end do

      ! calculate min_dx
      min_dx = min(minval(dxg),minval(dyg))

      !----------------------------------------------------------------------|
      ! set initial model - Lid driven cavity

      do j=1,jx
         do i=1,ix
            Vc(i,j,1) = ro0
            Vc(i,j,2) = pr0
            Vc(i,j,3) = vx0
            Vc(i,j,4) = vy0

            Vci(i,j,1) = Vc(i,j,1)
            Vci(i,j,2) = Vc(i,j,2)
            Vci(i,j,3) = Vc(i,j,3)
            Vci(i,j,4) = Vc(i,j,4)
         end do
      end do

   end subroutine model_setup


end module model

