!########################################
! Model : 2D Lid-Driven Cavity Flow
! Coordinate system : Cartesian grid
!########################################

program main
   use mpi_setup
   use openfile
   use const
   use init
   use getNewdt
   use integrate, only : integrate__EULER1

   implicit none
   real(8), parameter :: outputInterval = 0.1d0
   real(8) :: nextOutputTime
   ! integer, parameter :: outputInterval = 10
   ! integer :: nextOutputTime

   call initialize
   ! stop

   call file_output_param(nd,dtout,tend,ix,jx,igx,jgx,margin &
      ,mpid%mpisize,mpid%mpirank,mpisize_x,mpisize_y &
      ,g_gamma,x,y,dx,dy)

   call file_output(nvar,nd,mpid%mpirank,Vc,ix,jx)

   nd = nd + 1
   nextOutputTime = outputInterval

   loop: do ns=1,nstop

      !----- time update--------------------------------------------------------|
      call getNewdt__hd(margin,safety,dtmin,ix,jx,g_gamma,Vc,x,dx,y,dy&
         ,dt,min_dx,KINEMATIC_VISCOSITY)

      time = time+dt

      call integrate__EULER1(nvar,margin,ix,jx,g_gamma,x,dx,y,dy,dt,&
         Vc,Vci,&
         KINEMATIC_VISCOSITY,&
         N_PRESSURE_POISSON_ITERATIONS)

      ! if (ns == nextOutputTime) then
      if (time >= nextOutputTime) then
         call file_output(nvar,nd,mpid%mpirank,Vc,ix,jx)
         nextOutputTime = nextOutputTime + outputInterval
         nd = nd + 1
      end if

      if(time > tend) exit loop

      if(mpid%mpirank == 0) then
         print*, 'ns: ', ns
         print*, 'time: ', time
         print*, 'dt: ', dt
      end if

   end do loop

   call mpi_finalize(merr)

end program main

