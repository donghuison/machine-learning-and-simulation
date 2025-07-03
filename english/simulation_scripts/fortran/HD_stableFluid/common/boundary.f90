module boundary

   implicit none
   private

   public :: bd_consx, bd_consy, bd_frex, bd_frey, &
      bd_inix, bd_iniy,                     &
      bd_synnx_car, bd_synny_car,           &
      bd_synpx_car, bd_synpy_car,           &
      boundary__mpi, boundary__mpi_qq1


contains


   !=====================================================================
   ! INPUT
   !      mbnd :: margin flag
   !      ix,jx,kx :: array size
   !      margin :: margin
   ! INPUT&OUTPUT
   !      qq
   !=====================================================================

   subroutine bd_consx(mbnd,margin,cons,qq,ix,jx)
      implicit none

      integer,intent(in) :: ix,jx
      integer,intent(in) :: margin
      integer,intent(in) :: mbnd
      real(8),intent(in) :: cons
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: ibnd
      integer :: i,j

      if (mbnd .eq. 0)then
         ibnd = 1+margin


         do j=1,jx
            do i=1,margin
               qq(ibnd-i,j) = cons
            enddo
         enddo

      else
         ibnd = ix-margin


         do j=1,jx
            do i=1,margin
               qq(ibnd+i,j) = cons
            enddo
         enddo

      endif

      return
   end subroutine bd_consx


   subroutine bd_consy(mbnd,margin,cons,qq,ix,jx)
      implicit none

      integer,intent(in) :: mbnd ! boundary flag
      integer,intent(in) :: margin ! margin
      integer,intent(in) :: ix,jx! array size
      real(8),intent(in) :: cons
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: jbnd
      integer :: i,j

      if( mbnd .eq. 0)then
         jbnd = 1+margin


         do j=1,margin
            do i=1,ix
               qq(i,jbnd-j) = cons
            enddo
         enddo

      else
         jbnd = jx-margin


         do j=1,margin
            do i=1,ix
               qq(i,jbnd+j) = cons
            enddo
         enddo

      endif

      return
   end subroutine bd_consy


   subroutine bd_frex(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: mbnd ! boundary flag
      integer,intent(in) :: margin ! margin
      integer,intent(in) :: ix,jx ! array size
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: i,j

      if (mbnd .eq. 0)then

         do j=1,jx
            do i=1,margin
               qq(i,j) = qq(margin+1,j)
            enddo
         enddo

      else

         do j=1,jx
            do i=1,margin
               qq(ix-margin+i,j) = qq(ix-margin,j)
            enddo
         enddo

      endif

      return
   end subroutine bd_frex


   subroutine bd_frey(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: mbnd ! boundary flag
      integer,intent(in) :: margin ! margin
      integer,intent(in) :: ix,jx ! array size
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: jbnd
      integer :: i,j

      if( mbnd .eq. 0)then
         jbnd = 1+margin

         do j=1,margin
            do i=1,ix
               qq(i,jbnd-j) = qq(i,jbnd)
            enddo
         enddo

      else
         jbnd = jx-margin

         do j=1,margin
            do i=1,ix
               qq(i,jbnd+j) = qq(i,jbnd)
            enddo
         enddo

      endif

      return
   end subroutine bd_frey


   subroutine bd_inix(mbnd,margin,qq,qqini,ix,jx)
      implicit none

      integer,intent(in) :: mbnd ! boundary flag
      integer,intent(in) :: margin ! margin
      integer,intent(in) :: ix,jx ! array size
      real(8),dimension(ix,jx),intent(in) :: qqini ! initial quantity
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: ibnd
      integer :: i,j

      if (mbnd .eq. 0)then
         ibnd = 1+margin

         do j=1,jx
            do i=1,margin
               qq(ibnd-i,j) = qqini(ibnd-i,j)
            enddo
         enddo

      else
         ibnd = ix-margin


         do j=1,jx
            do i=1,margin
               qq(ibnd+i,j) = qqini(ibnd+i,j)
            enddo
         enddo

      endif

      return
   end subroutine bd_inix


   subroutine bd_iniy(mbnd,margin,qq,qqini,ix,jx)
      implicit none

      integer,intent(in) :: mbnd ! boundary flag
      integer,intent(in) :: margin ! margin
      integer,intent(in) :: ix,jx ! array size
      real(8),dimension(ix,jx),intent(in) :: qqini ! initial quantity
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: jbnd
      integer :: i,j,k

      if( mbnd .eq. 0)then
         jbnd = 1+margin

         do j=1,margin
            do i=1,ix
               qq(i,jbnd-j) = qqini(i,jbnd-j)
            enddo
         enddo

      else
         jbnd = jx-margin

         do j=1,margin
            do i=1,ix
               qq(i,jbnd+j) = qqini(i,jbnd+j)
            enddo
         enddo

      endif

      return
   end subroutine bd_iniy


   subroutine bd_synnx_car(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: margin,ix,jx
      integer,intent(in) :: mbnd
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: i,j

      if(mbnd .eq. 0)then

         do j=1,jx
            do i=1,margin
               qq(margin-i+1,j) = -(qq(margin+i,j))
            end do
         end do

      else

         do j=1,jx
            do i=1,margin
               qq(ix-margin+i,j) = -qq(ix-margin-i+1,j)
            end do
         end do

      end if

      return
   end subroutine bd_synnx_car


   subroutine bd_synny_car(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: margin,ix,jx
      integer,intent(in) :: mbnd

      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: i,j

      if(mbnd .eq. 0)then

         do j=1,margin
            do i=1,ix
               qq(i,margin-j+1) = -(qq(i,margin+j))
            end do
         end do

      else

         do j=1,margin
            do i=1,ix
               qq(i,jx-margin+j) = -qq(i,jx-margin-j+1)
            end do
         end do

      end if

      return
   end subroutine bd_synny_car


   subroutine bd_synpx_car(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: margin,ix,jx
      integer,intent(in) :: mbnd
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: i,j

      if(mbnd .eq. 0)then

         do j=1,jx
            do i=1,margin
               qq(margin-i+1,j) = (qq(margin+i,j))
            enddo
         end do

      else

         do j=1,jx
            do i=1,margin
               qq(ix-margin+i,j) = qq(ix-margin-i+1,j)
            end do
         end do

      end if

      return
   end subroutine bd_synpx_car


   subroutine bd_synpy_car(mbnd,margin,qq,ix,jx)
      implicit none

      integer,intent(in) :: margin,ix,jx
      integer,intent(in) :: mbnd
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer :: i,j

      if(mbnd .eq. 0)then

         do j=1,margin
            do i=1,ix
               qq(i,margin-j+1) = (qq(i,margin+j))
            end do
         end do

      else

         do j=1,margin
            do i=1,ix
               qq(i,jx-margin+j) = qq(i,jx-margin-j+1)
            end do
         end do

      end if

      return
   end subroutine bd_synpy_car


!    subroutine boundary__mpi(nvar,margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)
   subroutine boundary__mpi(nvar,margin,ix,jx,qq)

      use mpi_setup

      integer,intent(in) :: nvar
      integer,intent(in) :: ix,jx,margin
      real(8),dimension(ix,jx,nvar),intent(inout) :: qq

      !   integer,parameter :: mx = 10
      integer,parameter :: mx = 5
      integer :: i,j,mmx,msend,mrecv
      real(8),dimension(margin,jx,mx) :: bufsnd_x,bufrcv_x
      real(8),dimension(ix,margin,mx) :: bufsnd_y,bufrcv_y

      !=== Step 1.==========================================================
      ! surface exchange
      !=====================================================================
      !--- Step 1-1.---
      ! x-direction
      !----------------
      ! left

      mmx = margin*jx*mx
      msend = mpid%l
      mrecv = mpid%r


      do j=1,jx
         do i=1,margin
            bufsnd_x(i,j,1) = qq(margin+i,j,1)
            bufsnd_x(i,j,2) = qq(margin+i,j,2)
            bufsnd_x(i,j,3) = qq(margin+i,j,3)
            bufsnd_x(i,j,4) = qq(margin+i,j,4)
            bufsnd_x(i,j,5) = qq(margin+i,j,5)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_x,mmx,mdp,msend, 0 &
         ,bufrcv_x,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,jx
            do i=1,margin
               qq(ix-margin+i,j,1) = bufrcv_x(i,j,1)
               qq(ix-margin+i,j,2) = bufrcv_x(i,j,2)
               qq(ix-margin+i,j,3) = bufrcv_x(i,j,3)
               qq(ix-margin+i,j,4) = bufrcv_x(i,j,4)
               qq(ix-margin+i,j,5) = bufrcv_x(i,j,5)
            enddo
         enddo

      endif

      !--- Step 1-2.---
      ! right

      mmx = margin*jx*mx
      msend = mpid%r
      mrecv = mpid%l


      do j=1,jx
         do i=1,margin
            bufsnd_x(i,j,1) = qq(ix-2*margin+i,j,1)
            bufsnd_x(i,j,2) = qq(ix-2*margin+i,j,2)
            bufsnd_x(i,j,3) = qq(ix-2*margin+i,j,3)
            bufsnd_x(i,j,4) = qq(ix-2*margin+i,j,4)
            bufsnd_x(i,j,5) = qq(ix-2*margin+i,j,5)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_x,mmx,mdp,msend, 0 &
         ,bufrcv_x,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,jx
            do i=1,margin
               qq(i,j,1) = bufrcv_x(i,j,1)
               qq(i,j,2) = bufrcv_x(i,j,2)
               qq(i,j,3) = bufrcv_x(i,j,3)
               qq(i,j,4) = bufrcv_x(i,j,4)
               qq(i,j,5) = bufrcv_x(i,j,5)
            enddo
         enddo

      end if

      !--- Step 1-3.---
      ! y-direction
      !----------------
      ! back

      mmx = ix*margin*mx
      msend = mpid%b
      mrecv = mpid%f


      do j=1,margin
         do i=1,ix
            bufsnd_y(i,j,1) = qq(i,margin+j,1)
            bufsnd_y(i,j,2) = qq(i,margin+j,2)
            bufsnd_y(i,j,3) = qq(i,margin+j,3)
            bufsnd_y(i,j,4) = qq(i,margin+j,4)
            bufsnd_y(i,j,5) = qq(i,margin+j,5)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_y,mmx,mdp,msend, 0 &
         ,bufrcv_y,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,margin
            do i=1,ix
               qq(i,jx-margin+j,1) = bufrcv_y(i,j,1)
               qq(i,jx-margin+j,2) = bufrcv_y(i,j,2)
               qq(i,jx-margin+j,3) = bufrcv_y(i,j,3)
               qq(i,jx-margin+j,4) = bufrcv_y(i,j,4)
               qq(i,jx-margin+j,5) = bufrcv_y(i,j,5)
            enddo
         enddo

      end if

      !--- Step 1-4.---
      ! y-direction
      !----------------
      ! forth

      mmx = ix*margin*mx
      msend = mpid%f
      mrecv = mpid%b


      do j=1,margin
         do i=1,ix
            bufsnd_y(i,j,1) = qq(i,jx-2*margin+j,1)
            bufsnd_y(i,j,2) = qq(i,jx-2*margin+j,2)
            bufsnd_y(i,j,3) = qq(i,jx-2*margin+j,3)
            bufsnd_y(i,j,4) = qq(i,jx-2*margin+j,4)
            bufsnd_y(i,j,5) = qq(i,jx-2*margin+j,5)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_y,mmx,mdp,msend, 0 &
         ,bufrcv_y,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,margin
            do i=1,ix
               !    ro(i,j) = bufrcv_y(i,j,1)
               !    pr(i,j) = bufrcv_y(i,j,2)
               !    vx(i,j) = bufrcv_y(i,j,3)
               !    vy(i,j) = bufrcv_y(i,j,4)
               !    vz(i,j) = bufrcv_y(i,j,5)
               !    bx(i,j) = bufrcv_y(i,j,6)
               !    by(i,j) = bufrcv_y(i,j,7)
               !    bz(i,j) = bufrcv_y(i,j,8)
               !    phi(i,j) = bufrcv_y(i,j,9)
               !    eta(i,j) = bufrcv_y(i,j,10)
               qq(i,j,1) = bufrcv_y(i,j,1)
               qq(i,j,2) = bufrcv_y(i,j,2)
               qq(i,j,3) = bufrcv_y(i,j,3)
               qq(i,j,4) = bufrcv_y(i,j,4)
               qq(i,j,5) = bufrcv_y(i,j,5)
            enddo
         enddo

      end if

   end subroutine boundary__mpi


   subroutine boundary__mpi_qq1(margin,ix,jx,qq)

      use mpi_setup

      integer,intent(in) :: ix,jx,margin
      real(8),dimension(ix,jx),intent(inout) :: qq

      integer,parameter :: mx = 1
      integer :: i,j,mmx,msend,mrecv
      real(8),dimension(margin,jx,mx) :: bufsnd_x,bufrcv_x
      real(8),dimension(ix,margin,mx) :: bufsnd_y,bufrcv_y

      !=== Step 1.==========================================================
      ! surface exchange
      !=====================================================================
      !--- Step 1-1.---
      ! x-direction
      !----------------
      ! left

      mmx = margin*jx*mx
      msend = mpid%l
      mrecv = mpid%r


      do j=1,jx
         do i=1,margin
            bufsnd_x(i,j,1) = qq(margin+i,j)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_x,mmx,mdp,msend, 0 &
         ,bufrcv_x,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,jx
            do i=1,margin
               qq(ix-margin+i,j) = bufrcv_x(i,j,1)
            enddo
         enddo

      endif

      !--- Step 1-2.---
      ! right

      mmx = margin*jx*mx
      msend = mpid%r
      mrecv = mpid%l


      do j=1,jx
         do i=1,margin
            bufsnd_x(i,j,1) = qq(ix-2*margin+i,j)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_x,mmx,mdp,msend, 0 &
         ,bufrcv_x,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,jx
            do i=1,margin
               qq(i,j) = bufrcv_x(i,j,1)
            enddo
         enddo

      end if

      !--- Step 1-3.---
      ! y-direction
      !----------------
      ! back

      mmx = ix*margin*mx
      msend = mpid%b
      mrecv = mpid%f


      do j=1,margin
         do i=1,ix
            bufsnd_y(i,j,1) = qq(i,margin+j)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_y,mmx,mdp,msend, 0 &
         ,bufrcv_y,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,margin
            do i=1,ix
               qq(i,jx-margin+j) = bufrcv_y(i,j,1)
            enddo
         enddo

      end if

      !--- Step 1-4.---
      ! y-direction
      !----------------
      ! forth

      mmx = ix*margin*mx
      msend = mpid%f
      mrecv = mpid%b


      do j=1,margin
         do i=1,ix
            bufsnd_y(i,j,1) = qq(i,jx-2*margin+j)
         enddo
      enddo


      call mpi_sendrecv               &
         (bufsnd_y,mmx,mdp,msend, 0 &
         ,bufrcv_y,mmx,mdp,mrecv, 0 &
         ,mcomw,mstat,merr)

      if(mrecv /= mnull)then

         do j=1,margin
            do i=1,ix
               qq(i,j) = bufrcv_y(i,j,1)
            enddo
         enddo

      end if

   end subroutine boundary__mpi_qq1



end module boundary

