!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module MUIcoupling

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: MUI_init

contains

  subroutine MUI_init ()

    use decomp_2d_io
    use param , only : zptwofive
    use MPI

   !  coupling variables 
    use iso_c_binding
    use mui_3d_f
    use mui_general_f

    implicit none
    integer :: ifsIndx
    character(len=1024) :: numberSuffix
    NAMELIST/MUICoupling/domainName,interfaceName,interface_count,interfaceDirection, &
   interfacelocation,MUIBC_ID, groupNumb,groupVort, dataOrgShft,spatialScale,timeScale,tolerance,sendReceiveMode, &
   sptlSmpType,tmpSmpType,rSpatialSamp,sigmaSpatialSamp,rTempSamp,sigmaTempSamp,iSync,nSyncAhead,nForget,tempMeanSampLower,tempMeanSampUpper
    NAMELIST/MUICouplingBC/nclxCPL1,nclxCPLn,nclyCPL1,nclyCPLn,nclzCPL1,nclzCPLn

    if (trim(MUIcommandArgs).eq.trim('coupled')) then
      open(10, file="./input.mui") 
      read(10, nml=MUICoupling); rewind(10)  
      read(10, nml=MUICouplingBC); rewind(10)  
      close(10)
      allocate(character(len_trim(interfaceName)+5) :: interfaces3d(interface_count))
      !For multi-domain function, "uniface_pointers_1d" should be used to collect the array of
      ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
      allocate(uniface_pointers_3d(interface_count),ifsDir(interface_count),ifsLoc(interface_count))
   
      do ifsIndx = 1, interface_count
         !Generate character type of number suffix
         if (ifsIndx < 10) then
            write (numberSuffix, "(I1)") ifsIndx
         else if ((ifsIndx < 100) .and. (ifsIndx > 9)) then
            write (numberSuffix, "(I2)") ifsIndx
         else if ((ifsIndx < 1000) .and. (ifsIndx > 99)) then
            write (numberSuffix, "(I3)") ifsIndx
         else
            write (numberSuffix, "(I4)") ifsIndx
         endif
   
         !Create and collect interface names
         interfaces3d(ifsIndx) = trim(interfaceName) // "_" // trim(numberSuffix)
         ifsDir(ifsIndx)=interfaceDirection(ifsIndx)
         ifsLoc(ifsIndx)=interfaceLocation(ifsIndx)
      end do 
      call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, trim(domainName), interfaces3d, interface_count)
      call MUI_create_sampler()
   endif
    return
  end subroutine MUI_init
  !********************************************************************


subroutine MUI_create_sampler()
#ifdef MUI_COUPLING
   ! Coupling varaibles 
      use iso_c_binding
      use mui_3d_f
      use mui_general_f
#endif
      use decomp_2d_io
      use MPI
      integer :: i
   !! limit the pushing boxs dimensions to the domain size
   do i=1,groupNumb  
      groupVort(i,1) =max(groupVort(i,1),1)
      groupVort(i,3) =max(groupVort(i,3),1)
      groupVort(i,5) =max(groupVort(i,5),1)

      groupVort(i,2) =min(groupVort(i,2),nx)
      groupVort(i,4) =min(groupVort(i,4),ny)
      groupVort(i,6) =min(groupVort(i,6),nz)
   enddo
   !!! Set spatial sampler 
   if (trim(sptlSmpType)==trim('exact')) then
      call mui_create_sampler_exact_3d_f(spatial_sampler, tolerance)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =", trim(sptlSmpType)
         write(*,*)  "     - Tolerance =", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('linear')) then
      call mui_create_sampler_pseudo_n2_linear_3d_f(spatial_sampler,rSpatialSamp)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =   ", trim(sptlSmpType)
         write(*,*)  "     - rSpatialSamp  = ", rSpatialSamp
         write(*,*)  "     - sigmaSpatialSamp  = ", sigmaSpatialSamp
         write(*,*)  "     - Tolerance = ", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('gauss')) then
      call mui_create_sampler_gauss_3d_f(spatial_sampler,rSpatialSamp,sigmaSpatialSamp)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =   ", trim(sptlSmpType)
         write(*,*)  "     - rSpatialSamp  = ", rSpatialSamp
         write(*,*)  "     - sigmaSpatialSamp  = ", sigmaSpatialSamp
         write(*,*)  "     - Tolerance = ", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('RBF')) then
      !! to be added 
      if (nrank==0) then 
         Write(*,*) trim(sptlSmpType), "To be added"
         write(*,*) "The available spatial samplers are exact, RBF"
      endif
      stop
   else
      if (nrank==0) then 
         Write(*,*) trim(sptlSmpType), "is wrong option for MUI spatial Sampler"
         write(*,*) "The available spatial samplers are exact, RBF"
         write(*,*)  "======================================================="
      endif
      stop
   endif

   !!! Set Temporal sampler 
   if (trim(tmpSmpType)==trim('exact')) then
      call mui_create_temporal_sampler_exact_3d_f(temporal_sampler, tolerance)
   elseif (trim(tmpSmpType)==trim('gauss')) then
      call mui_create_temporal_sampler_gauss_3d_f(temporal_sampler, rTempSamp,sigmaTempSamp)
      ! write(*,*) " XXXXXXXXXXXXXXXXXXXXX Warning XXXXXXXXXXXXXXXXXXXXX "
      ! write(*,*) " The Temporal Gauss sampler has not been tested for xcompact3D."
      ! write(*,*) "  Do a proper testing before continue using it"
      ! write(*,*) " XXXXXXXXXXXXXXXXXXXXX Warning XXXXXXXXXXXXXXXXXXXXX "
   elseif (trim(tmpSmpType)==trim('mean')) then
      call mui_create_temporal_sampler_mean_3d_f(temporal_sampler, tempMeanSampLower,tempMeanSampUpper)
      write(*,*) " XXXXXXXXXXXXXXXXXXXXX Warning XXXXXXXXXXXXXXXXXXXXX "
      write(*,*) " The Temporal Mean sampler has not been tested for xcompact3D."
      write(*,*) "  Do a proper testing before continue using it"
      write(*,*) " XXXXXXXXXXXXXXXXXXXXX Warning XXXXXXXXXXXXXXXXXXXXX "
   else
      if (nrank==0) then 
         write(*,*) "MUI Error: The temporal sampler ", tmpSmpType, "is not implemented"
      endif
      stop
   endif

   ! Set the fetch function 
   if (trim(sptlSmpType)==trim('exact') .and. trim(tmpSmpType)==trim('exact') ) then
      mui_fetch => mui_fetch_exact_exact_3d_f
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('exact')) then
      mui_fetch => mui_fetch_gauss_exact_3d_f
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('exact')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_exact_3d_f


   elseif (trim(sptlSmpType)==trim('exact') .and. trim(tmpSmpType)==trim('gauss') ) then
      mui_fetch => mui_fetch_exact_gauss_3d_f
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('gauss')) then
      mui_fetch => mui_fetch_gauss_gauss_3d_f
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('gauss')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_gauss_3d_f

   elseif (trim(sptlSmpType)==trim('exact') .and. trim(tmpSmpType)==trim('mean') ) then
      mui_fetch => mui_fetch_exact_mean_3d_f
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('mean')) then
      mui_fetch => mui_fetch_gauss_mean_3d_f
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('mean')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_mean_3d_f



   else
      if (nrank==0) then 
         write(*,*) "The selected spatial and temporal sampler types have not impelemnted yet"
      endif
      stop
   endif

end subroutine MUI_create_sampler



  subroutine visu_MUIBC_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_MUIBC_init
  

end module MUIcoupling
