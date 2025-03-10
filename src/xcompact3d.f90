!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case
  use variables, only : MUIcommandArgs
  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d
  use MUIcoupledBC, only : pushMUI
  

  implicit none

  real :: tstartMUI,tBCMUI,tPushMUI,timeMUI_push,timeMUI_fetch,timeMUI_fetch_inst,timeMUI_push_inst

  call init_xcompact3d()

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     call simu_stats(2)    ! screen/log output
     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((ifilter.ne.0).and.(ilesmod.ne.0)) then
        if (nrank==0) write(*,*) "Filtering the velocity and scalar fields. ifilter and C_filter =",ifilter, C_filter
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)

#ifdef MUI_COUPLING
      ! Set the order of the coupling operations
      if (trim(MUIcommandArgs).eq.trim('coupled')) then
         
         if (sendReceiveMode==0) then ! No Sending 
         call cpu_time(tstartMUI)
            call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
            call cpu_time(tBCMUI)
            timeMUI_fetch=timeMUI_fetch+(tBCMUI-tstartMUI)
            timeMUI_fetch_inst=(tBCMUI-tstartMUI)
        else if (sendReceiveMode==1) then ! Receive first 

            call cpu_time(tstartMUI)
            call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
            call cpu_time(tBCMUI)
            timeMUI_fetch=timeMUI_fetch+(tBCMUI-tstartMUI)
            timeMUI_fetch_inst=(tBCMUI-tstartMUI)

            call cpu_time(tBCMUI)
            call pushMUI(ux1, uy1, uz1)
            call cpu_time(tPushMUI)
            timeMUI_push=timeMUI_push+(tPushMUI-tBCMUI)
            timeMUI_push_inst=(tPushMUI-tBCMUI)

        else if (sendReceiveMode==2) then  ! push first 
            call cpu_time(tstartMUI)
            call pushMUI(ux1, uy1, uz1)
            call cpu_time(tPushMUI)
            timeMUI_push=timeMUI_push+(tPushMUI-tstartMUI)
            timeMUI_push_inst=(tPushMUI-tstartMUI)

            call cpu_time(tPushMUI)
            call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
            call cpu_time(tBCMUI)
            timeMUI_fetch=timeMUI_fetch+(tBCMUI-tPushMUI)
            timeMUI_fetch_inst=(tBCMUI-tPushMUI)
        else 
            write(*,*) "Wrong selction of sendReceiveMode to ", sendReceiveMode
            stop
        endif
        if (itime .le. ifirst+1) then
         call simu_stats(1)
         timeMUI_fetch = 0.0
         timeMUI_push = 0.0
        endif
        
        
      else
         call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
      endif
#else
      call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)
#endif
        if (imove.eq.1) then ! update epsi for moving objects
          if ((iibm.eq.2).or.(iibm.eq.3)) then
             call genepsi3d(ep1)
          else if (iibm.eq.1) then
             call body(ux1,uy1,uz1,ep1)
          endif
        endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)
#ifdef DEBG
        call check_transients()
#endif
        
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)
     if (nrank==0) write(*,*) trim(domainName)," : MUI inst comunication time is Fetch/push : ", &
      timeMUI_fetch_inst,"/", timeMUI_push_inst, " seconds"
     if (nrank==0) write(*,*) trim(domainName)," : MUI average comunication time is Fetch/push : ", &
      timeMUI_fetch/(itime-ifirst-1)/iadvance_time ,"/", timeMUI_push/(itime-ifirst-1)/iadvance_time, " seconds"

     call simu_stats(3) ! screen/log output

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

  enddo !! End time loop

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces

!   use iso_c_binding
#ifdef MUI_COUPLING
use mpi_f08 , only : MPI_comm
  use iso_c_binding 
  use mui_3d_f
  use mui_general_f
  use variables, only :  MUI_COMM_WORLD, uniface_3d
  
  
#endif

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, compute_cfldiff, &
       init_inflow_outflow

  use param, only : ilesmod, jles,itype
  use param, only : irestart

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col,MUIcommandArgs
  use variables, only : nstat, nvisu, nprobe, ilist

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  use probes, only : init_probes
#ifdef MUI_COUPLING
  use MUIcoupling
#endif
  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, argTemp,FNBase

  !! Initialise MPI
!   
#ifdef MUI_COUPLING
   ! call MPI_INIT(ierr)
   call mui_mpi_split_by_app_f(MUI_COMM_WORLD)
#else
  call MPI_INIT(ierr)
  MUI_COMM_WORLD = MPI_COMM_WORLD
#endif
  call MPI_COMM_RANK(MUI_COMM_WORLD,nrank,ierr) 
  call MPI_COMM_SIZE(MUI_COMM_WORLD,nproc,ierr)
  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
     InputFN='./input.i3d'

   if (nargin == 1) then
     call get_command_argument(1,argTemp,FNLength,status)
     if (trim(argTemp).eq.trim('coupled')) then
      MUIcommandArgs = argTemp
     else
      InputFN = argTemp
      back=.true.
      FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
      DecInd=index(FNBase,'.',back)
      if (DecInd >1) then
         FNBase=FNBase(1:(DecInd-1))
      end if
     endif
   elseif (nargin== 2) then
      Write(*,*) "Error: Need fexing if you are going to use more than one argument in the rum command"
      STOP
  endif
  if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file --> ', trim(InputFN)

#ifdef ADIOS2
  if (nrank .eq. 0) then
     write(*,*) " WARNING === WARNING === WARNING === WARNING === WARNING"
     write(*,*) " WARNING: Running Xcompact3d with ADIOS2"
     write(*,*) "          this is currently experimental"
     write(*,*) "          for safety of results it is recommended"
     write(*,*) "          to run the default build as this feature"
     write(*,*) "          is developed. Thank you for trying it."
     write(*,*) " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)



  call decomp_2d_init(nx,ny,nz,p_row,p_col,MUI_COMM_WORLD)
  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

#ifdef MUI_COUPLING
  if (trim(MUIcommandArgs).eq.trim('coupled')) call MUI_init ()
#endif

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces.eq.1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((ioutflow.eq.1).or.(iin.eq.3)) then
     call init_inflow_outflow()
  end if

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise
  use variables, only : MUIcommandArgs
  use tools, only : simu_stats
  use param, only : itype, jles, ilesmod
  use probes, only : finalize_probes
  use visu, only : visu_finalise
  use les, only: finalise_explicit_les
#ifdef MUI_COUPLING
   use iso_c_binding
   use mui_3d_f
   use mui_general_f
#endif
  implicit none

  integer :: ierr
  real :: start, finish
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif
#ifdef MUI_COUPLING
if (trim(MUIcommandArgs).eq.trim('coupled'))  deallocate(uniface_pointers_3d)
#endif
  call simu_stats(4)
  call finalize_probes()
  call visu_finalise()
  if (ilesmod.ne.0) then
     if (jles.gt.0) call finalise_explicit_les()
  endif
  ! pause for 20 sec before terminating the code.
  call cpu_time(start)
  do
        call cpu_time(finish)
        if (finish - start >= 20.0) exit
  end do
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d

subroutine check_transients()

  use decomp_2d, only : mytype
  use mpi
  use var
  
  implicit none

  real(mytype) :: dep, dep1
  integer :: code
   
  dep=maxval(abs(dux1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MUI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX dux1 ', dep1
 
  dep=maxval(abs(duy1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MUI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duy1 ', dep1
 
  dep=maxval(abs(duz1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MUI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duz1 ', dep1
  
end subroutine check_transients
