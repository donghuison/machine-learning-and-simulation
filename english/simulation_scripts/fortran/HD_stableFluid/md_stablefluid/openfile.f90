module openfile

   use dac_header
   use const, only : input_dir,output_dir

   implicit none
   private

   public :: file_input, file_output, file_output_param, file_output_cal_time, file_output_flux

   character :: cno*4
   character :: cnond*4

contains


   subroutine file_input(nd,mpirank,ro,pr,vx,vy,vz,ix,jx)

      integer,intent(in) :: nd,mpirank,ix,jx
      real(8),intent(inout),dimension(ix,jx) :: ro,pr,vx,vy,vz

      integer :: mfi_ro,mfi_pr
      integer :: mfi_vx,mfi_vy,mfi_vz
      integer :: mt,nx0,ix0,jx0

      mt=6
      write(cnond,'(i4.4)') nd
      write(cno,'(i4.4)') mpirank

      mfi_ro=70
      call dacopnr2s(mfi_ro,input_dir//cnond//'_ro_rank='//cno//'.dac',mt,ix0,jx0,nx0)
      read(mfi_ro) ro
      close(mfi_ro)

      mfi_pr=71
      call dacopnr2s(mfi_pr,input_dir//cnond//'_pr_rank='//cno//'.dac',mt,ix0,jx0,nx0)
      read(mfi_pr) pr
      close(mfi_pr)

      mfi_vx=72
      call dacopnr2s(mfi_vx,input_dir//cnond//'_vx_rank='//cno//'.dac',mt,ix0,jx0,nx0)
      read(mfi_vx) vx
      close(mfi_vx)

      mfi_vy=73
      call dacopnr2s(mfi_vy,input_dir//cnond//'_vy_rank='//cno//'.dac',mt,ix0,jx0,nx0)
      read(mfi_vy) vy
      close(mfi_vy)

      mfi_vz=74
      call dacopnr2s(mfi_vz,input_dir//cnond//'_vz_rank='//cno//'.dac',mt,ix0,jx0,nx0)
      read(mfi_vz) vz
      close(mfi_vz)

   end subroutine file_input


!    subroutine file_output(nvar,nd,mpirank,ro,pr,vx,vy,vz,er,frx,fry,frz,ix,jx)
   subroutine file_output(nvar,nd,mpirank,Vc,ix,jx)

      integer,intent(in) :: nvar
      integer,intent(in) :: nd,mpirank,ix,jx
      real(8),intent(in),dimension(ix,jx,nvar) :: Vc
      !   real(8),intent(in),dimension(ix,jx) :: ro,pr,vx,vy,vz
      !   real(8),intent(in),dimension(ix,jx) :: er,frx,fry,frz

      integer :: mf_ro,mf_pr
      integer :: mf_vx,mf_vy,mf_vz

      write(cnond,'(i4.4)') nd
      write(cno,'(i4.4)') mpirank

      mf_ro=20
      call dacdef2s(mf_ro,output_dir//cnond//'_ro_rank='//cno//'.dac',6,ix,jx)
      write(mf_ro) Vc(:,:,1)
      close(mf_ro)

      mf_pr=21
      call dacdef2s(mf_pr,output_dir//cnond//'_pr_rank='//cno//'.dac',6,ix,jx)
      write(mf_pr) Vc(:,:,2)
      close(mf_pr)

      mf_vx=22
      call dacdef2s(mf_vx,output_dir//cnond//'_vx_rank='//cno//'.dac',6,ix,jx)
      write(mf_vx) Vc(:,:,3)
      close(mf_vx)

      mf_vy=23
      call dacdef2s(mf_vy,output_dir//cnond//'_vy_rank='//cno//'.dac',6,ix,jx)
      write(mf_vy) Vc(:,:,4)
      close(mf_vy)

      mf_vz=24
      call dacdef2s(mf_vz,output_dir//cnond//'_vz_rank='//cno//'.dac',6,ix,jx)
      write(mf_vz) Vc(:,:,5)
      close(mf_vz)

   end subroutine file_output

   subroutine file_output_param(nd,dtout,tend,ix,jx,igx,jgx,margin,mpisize &
      ,mpirank,mpisize_x,mpisize_y,gm,x,y,dx,dy)

      integer,intent(in) :: nd,ix,jx,igx,jgx,margin,mpisize,mpirank
      integer,intent(in) :: mpisize_x,mpisize_y
      real(8),intent(in),dimension(ix) :: x,dx
      real(8),intent(in),dimension(jx) :: y,dy
      real(8),intent(in) :: dtout,tend
      real(8),intent(in) :: gm

      integer :: mf_params
      integer :: mf_x,mf_y

      write(cnond,'(i4.4)') nd
      write(cno,'(i4.4)') mpirank

      mf_params=9
      open (mf_params,file=output_dir//'params_rank='//cno//'.txt',form='formatted')

      call dacputparamc(mf_params,'comment','model_parker')
      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'igx',igx)
      call dacputparami(mf_params,'jgx',jgx)
      call dacputparami(mf_params,'margin',margin)
      call dacputparamd(mf_params,'tend',tend)
      call dacputparami(mf_params,'mpi',1)
      call dacputparamd(mf_params,'dtout',dtout)
      call dacputparami(mf_params,'mpisize',mpisize)
      call dacputparami(mf_params,'mpirank',mpirank)
      call dacputparami(mf_params,'mpix',mpisize_x)
      call dacputparami(mf_params,'mpiy',mpisize_y)
      call dacputparamd(mf_params,'x(1)',x(1))
      call dacputparamd(mf_params,'y(1)',y(1))
      call dacputparamd(mf_params,'dx(1)',dx(1))
      call dacputparamd(mf_params,'dy(1)',dy(1))
      call dacputparamd(mf_params,'gm',gm)
      close(mf_params)

      mf_x=11
      call dacdef1d(mf_x,output_dir//'x_rank='//cno//'.dac',6,ix)
      write(mf_x) x
      close(mf_x)

      mf_y=12
      call dacdef1d(mf_y,output_dir//'y_rank='//cno//'.dac',6,jx)
      write(mf_y) y
      close(mf_y)

   end subroutine file_output_param


   subroutine file_output_cal_time(nd,dtout,tend,ix,jx,igx,jgx,margin,mpisize &
      ,mpirank,mpisize_x,mpisize_y,dt,cal_time,iter,time,x,y,dx,dy)

      integer,intent(in) :: nd,ix,jx,igx,jgx,margin,mpisize,mpirank
      integer,intent(in) :: mpisize_x,mpisize_y,iter
      real(8),intent(in),dimension(ix) :: x,dx
      real(8),intent(in),dimension(jx) :: y,dy
      real(8),intent(in) :: dtout,tend
      real(8),intent(in) :: dt,time
      real(8),intent(in) :: cal_time

      integer :: mf_params
      integer :: mf_x,mf_y

      write(cnond,'(i4.4)') nd
      write(cno,'(i4.4)') mpirank

      mf_params=9
      open (mf_params,file=output_dir//'params_rank='//cno//'.txt',form='formatted')

      call dacputparamc(mf_params,'comment','nonlinear_parker_instability')
      call dacputparami(mf_params,'ix',ix)
      call dacputparami(mf_params,'jx',jx)
      call dacputparami(mf_params,'igx',igx)
      call dacputparami(mf_params,'jgx',jgx)
      call dacputparami(mf_params,'margin',margin)
      call dacputparami(mf_params,'iterations',iter)
      call dacputparamd(mf_params,'tend',tend)
      call dacputparami(mf_params,'mpi',1)
      call dacputparamd(mf_params,'dtout',dtout)
      call dacputparami(mf_params,'mpisize',mpisize)
      call dacputparami(mf_params,'mpirank',mpirank)
      call dacputparami(mf_params,'mpix',mpisize_x)
      call dacputparami(mf_params,'mpiy',mpisize_y)
      call dacputparamd(mf_params,'x(1)',x(1))
      call dacputparamd(mf_params,'y(1)',y(1))
      call dacputparamd(mf_params,'dx(1)',dx(1))
      call dacputparamd(mf_params,'dy(1)',dy(1))
      call dacputparamd(mf_params,'dt',dt)
      call dacputparamd(mf_params,'calculation time',cal_time)
      call dacputparamd(mf_params,'time',time)
      close(mf_params)

      mf_x=11
      call dacdef1d(mf_x,output_dir//'x_rank='//cno//'.dac',6,ix)
      write(mf_x) x
      close(mf_x)

      mf_y=12
      call dacdef1d(mf_y,output_dir//'y_rank='//cno//'.dac',6,jx)
      write(mf_y) y
      close(mf_y)

   end subroutine file_output_cal_time



   subroutine file_output_flux(nd,mpirank,roy,feey,ix,jx)

      integer,intent(in) :: nd,mpirank,ix,jx
      real(8),intent(in),dimension(ix,jx) :: roy,feey

      integer :: mf_roy,mf_feey

      write(cnond,'(i4.4)') nd
      write(cno,'(i4.4)') mpirank


      mf_roy=30
      call dacdef2s(mf_roy,output_dir//cnond//'_roy_rank='//cno//'.dac',6,ix,jx)
      write(mf_roy) roy
      close(mf_roy)

      mf_feey=31
      call dacdef2s(mf_feey,output_dir//cnond//'_feey_rank='//cno//'.dac',6,ix,jx)
      write(mf_feey) feey
      close(mf_feey)

   end subroutine file_output_flux


end module openfile

