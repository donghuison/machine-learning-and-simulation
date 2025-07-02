!########################################
! Model : 2D Stable Fluids Simulation
! Coordinate system : Cartesian grid
!########################################

program main
   use mpi_setup
   use openfile
   use const
   use init
   use integrate, only : integrate__StableFluids

   implicit none
   real(8), parameter :: outputInterval = 1.0d0  ! Output every 1.0 time units
   real(8) :: nextOutputTime

   call initialize

   call file_output_param(nd,dtout,tend,ix,jx,igx,jgx,margin &
      ,mpid%mpisize,mpid%mpirank,mpisize_x,mpisize_y &
      ,g_gamma,x,y,dx,dy)

   call file_output(nvar,nd,mpid%mpirank,Vc,ix,jx)

   nd = nd + 1
   nextOutputTime = outputInterval

   loop: do ns=1,nstop

      !----- For stable fluids, we use fixed timestep -----------------------|
      dt = dt0
      time = time + dt

      call integrate__StableFluids(nvar,margin,ix,jx,g_gamma,x,dx,y,dy,dt,time,&
         Vc,Vci,&
         KINEMATIC_VISCOSITY,&
         MAX_ITER_CG,&
         xmin,xmax,ymin,ymax)

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

