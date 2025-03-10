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
    integer :: ifsIndx,code,ierror,i,j,k,iGrp,synchronised=1
    character(len=1024) :: numberSuffix
    real(mytype) :: dataLoc(3),grdSpce(3), veryLarg = 10e34
    real(mytype) :: xMinFetch,yMinFetch,zMinFetch,xMaxFetch,yMaxFetch,zMaxFetch
    real(mytype) :: xMin  ,yMin , zMin, xMax ,yMax, zMax
    NAMELIST/MUICoupling/domainName,interfaceName,interface_count, pushInterfacefsID , fetchInterfacefsID, &
     groupNumb,groupVort, dataOrgShft,spatialScale,timeScale,tolerance,sendReceiveMode, &
   sptlSmpType,tmpSmpType,rSpatialSamp,sigmaSpatialSamp,rTempSamp,sigmaTempSamp,iSync,nSyncAhead,&
   nForget,tempMeanSampLower,tempMeanSampUpper,smartSendReceive
    NAMELIST/MUICouplingBC/nclxCPL1,nclxCPLn,nclyCPL1,nclyCPLn,nclzCPL1,nclzCPLn, &
                           nclxCrr1,nclxCrrn,nclyCrr1,nclyCrrn,nclzCrr1,nclzCrrn

    if (trim(MUIcommandArgs).eq.trim('coupled')) then

      open(10, file="./input.mui") 
      read(10, nml=MUICoupling); rewind(10)  
      read(10, nml=MUICouplingBC); rewind(10)  
      close(10)
      ! sanity checks
      if (pushInterfacefsID .gt. interface_count .or. fetchInterfacefsID .gt. interface_count) then 
         if (nrank == 0) write(*,*) "ERROR: pushInterfacefsID and fetchInterfacefsID can't be larger than interface_count"
         call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
         stop
      endif

      allocate(character(len_trim(interfaceName)+5) :: interfaces3d(interface_count))
      allocate(uniface_pointers_3d(interface_count))
      allocate(groupVortLocal(groupNumb,6))
   
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
      end do 
      call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, trim(domainName), interfaces3d, interface_count)
      call MUI_create_sampler()

      !!! Define the local index of points of each push  group 
      groupVortLocal = 0
      do iGrp = 1, groupNumb
         do j=1,3
            i=j*2-1
            if (groupVort(iGrp,i)>=xstart(j) .and. groupVort(iGrp,i) <=xend(j)) then
               groupVortLocal(iGrp,i)=groupVort(iGrp,i)-xstart(j)+1
            else if (groupVort(iGrp,i)<xstart(j) .and. groupVort(iGrp,i+1) >=xstart(j).and. groupVort(iGrp,i+1)<=xend(j)) then
               groupVortLocal(iGrp,i) = 1
            endif

            if (groupVort(iGrp,i+1)>=xstart(j) .and. groupVort(iGrp,i+1) <=xend(j)) then
               groupVortLocal(iGrp,i+1)=groupVort(iGrp,i+1)-xstart(j)+1
            else if (groupVort(iGrp,i) >= xstart(j).and.groupVort(iGrp,i) <= xend(j).and.groupVort(iGrp,i+1) >=xend(j)) then
               groupVortLocal(iGrp,i+1) = xSize(j)
            endif

            if (groupVort(iGrp,i)<=xstart(j) .and. groupVort(iGrp,i+1) >=xend(j)) then
               groupVortLocal(iGrp,i) =1
               groupVortLocal(iGrp,i+1) = xSize(j)
            endif

         enddo
         ! If a  group does not have a data in this local domain, set it to zero and -1. Use -1 to avoid entering the loop. 
         if (minval(groupVortLocal(iGrp,:)) ==0)  groupVortLocal(iGrp,:) = [0, 0, 0, -1, -1, -1]  
      enddo 


      !!!! Find min and max coordinates of each local group
      if (smartSendReceive .eq. 1) then 
         call MUI_initiateSmartSendRecieve()
      endif
      

   endif
    return
  end subroutine MUI_init
  !********************************************************************

  subroutine MUI_initiateSmartSendRecieve ()

   use decomp_2d_io
   use param , only : zptwofive
   use MPI

  !  coupling variables 
   use iso_c_binding
   use mui_3d_f
   use mui_general_f

   implicit none
   integer :: ifsIndx,code,ierror,i,j,k,iGrp,synchronised=1,sendig_peers_size,ierr
   character(len=1024) :: numberSuffix
   real(mytype) :: dataLoc(3),grdSpce(3), veryLarg = 10e34
   real(mytype) :: xMinFetch,yMinFetch,zMinFetch,xMaxFetch,yMaxFetch,zMaxFetch
   real(mytype) :: xMin  ,yMin , zMin, xMax ,yMax, zMax

   grdSpce(1) = xlx/(nx-1)
   grdSpce(2) = yly/(ny-1)
   grdSpce(3) = zlz/(nz-1)

   xMin = veryLarg  ;yMin = veryLarg ; zMin = veryLarg
   xMax = -veryLarg ;yMax = -veryLarg; zMax = -veryLarg

   do iGrp = 1, groupNumb
      do k=groupVortLocal(iGrp,5),groupVortLocal(iGrp,6)
         do j=groupVortLocal(iGrp,3),groupVortLocal(iGrp,4)
            do i=groupVortLocal(iGrp,1),groupVortLocal(iGrp,2)
               dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1)  
               dataLoc(2)  = yp(j+xstart(2)-1) 
               dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3) 
               dataLoc = dataLoc * spatialScale + dataOrgShft

               if(dataLoc(1) .lt. xMin) xMin = dataLoc(1)
               if(dataLoc(2) .lt. yMin) yMin = dataLoc(2)
               if(dataLoc(3) .lt. zMin) zMin = dataLoc(3)
               if(dataLoc(1) .gt. xMax) xMax = dataLoc(1)
               if(dataLoc(2) .gt. yMax) yMax = dataLoc(2)
               if(dataLoc(3) .gt. zMax) zMax = dataLoc(3)
            enddo
         enddo
      enddo         
   enddo

   xMin = xMin - grdSpce(1)
   yMin = yMin - grdSpce(2)
   zMin = zMin - grdSpce(3)

   xMax = xMax + grdSpce(1)
   yMax = yMax + grdSpce(2)
   zMax = zMax + grdSpce(3)

   !!! if there is no group in this rank, set the dimensions to large value to avoid unnecessary communication
   if (xMin .ge. veryLarg) xMin = veryLarg
   if (yMin .ge. veryLarg) yMin = veryLarg
   if (zMin .ge. veryLarg) zMin = veryLarg

   if (xMax .le. -veryLarg) xMax = veryLarg
   if (yMax .le. -veryLarg) yMax = veryLarg
   if (zMax .le. -veryLarg) zMax = veryLarg

   
   xMinFetch = real((xstart(1)-1),mytype)*grdSpce(1) * spatialScale + dataOrgShft(1)
   yMinFetch = yp(xstart(2))                         * spatialScale + dataOrgShft(2)
   zMinFetch = real((xstart(3)-1),mytype)*grdSpce(3) * spatialScale + dataOrgShft(3)

   xMaxFetch = real((xend(1)-1),mytype)*grdSpce(1) * spatialScale + dataOrgShft(1)
   yMaxFetch = yp(xend(2))                         * spatialScale + dataOrgShft(2)
   zMaxFetch = real((xend(3)-1),mytype)*grdSpce(3) * spatialScale + dataOrgShft(3)     

!! Increase the box size to allow spatial interpolation across ranks
xMinFetch = xMinFetch - grdSpce(1)
yMinFetch = yMinFetch - grdSpce(2)
zMinFetch = zMinFetch - grdSpce(3)

xMaxFetch = xMaxFetch + grdSpce(1)
yMaxFetch = yMaxFetch + grdSpce(2)
zMaxFetch = zMaxFetch + grdSpce(3)
 


   do i=1,interface_count

      

      call mui_announce_recv_span_3d_box_f(uniface_pointers_3d(i)%ptr, &
      xMinFetch,yMinFetch,zMinFetch,xMaxFetch,yMaxFetch,zMaxFetch ,real(ifirst-1,c_double)*dt,real(ilast-ifirst+1, c_double)*dt,synchronised)

      call mui_announce_send_span_3d_box_f(uniface_pointers_3d(i)%ptr, &
      xMin,yMin, zMin,xMax,yMax,zMax,real(ifirst-1,c_double)*dt,real(ilast-ifirst+1, c_double)*dt,synchronised)
      
      call mui_commit_3d_f(uniface_pointers_3d(i)%ptr, real(ifirst-1,c_double)*dt)
   enddo

   
   call mui_sendig_peers_size_3d_f(uniface_pointers_3d(pushInterfacefsID)%ptr, sendig_peers_size)
   if (nrank==0 ) write(*,*) "==================================================================="


   call MPI_Barrier(MUI_COMM_WORLD, ierr)
   write(*,*)"MUI domain ", trim(domainName), " and rank ", nrank," has a sending peer size at pushInterfacefsID " , sendig_peers_size
   call MPI_Barrier(MUI_COMM_WORLD, ierr)


   call mui_sendig_peers_size_3d_f(uniface_pointers_3d(fetchInterfacefsID)%ptr, sendig_peers_size)
   if (nrank==0 ) write(*,*) "==================================================================="
   if (nrank==0 ) write(*,*)"MUI ", trim(domainName), ": Activating smart send and receive"

   call MPI_Barrier(MUI_COMM_WORLD, ierr)
   write(*,*)"MUI domain ", trim(domainName), " and rank ", nrank," has a sending peer size at fetchInterfacefsID " , sendig_peers_size
   call MPI_Barrier(MUI_COMM_WORLD, ierr)
   if (nrank==0 ) write(*,*) "==================================================================="


  
   return
 end subroutine MUI_initiateSmartSendRecieve
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
      integer :: i,j,k,cnt,point_count
      real(mytype) :: grdSpce(3)
      real(c_double), dimension (:), allocatable :: point_x, point_y,point_z

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
      call mui_create_sampler_exact_3d_vector(spatial_sampler, tolerance)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =", trim(sptlSmpType)
         write(*,*)  "     - Tolerance =", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('linear')) then
      call mui_create_sampler_pseudo_n2_linear_3d_vector(spatial_sampler,rSpatialSamp)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =   ", trim(sptlSmpType)
         write(*,*)  "     - rSpatialSamp  = ", rSpatialSamp
         write(*,*)  "     - sigmaSpatialSamp  = ", sigmaSpatialSamp
         write(*,*)  "     - Tolerance = ", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('gauss')) then
      call mui_create_sampler_gauss_3d_vector(spatial_sampler,rSpatialSamp,sigmaSpatialSamp)
      if (nrank==0) then
         write(*,*) "MUI spatial Sampler is created with: "
         write(*,*)  "     - Type      =   ", trim(sptlSmpType)
         write(*,*)  "     - rSpatialSamp  = ", rSpatialSamp
         write(*,*)  "     - sigmaSpatialSamp  = ", sigmaSpatialSamp
         write(*,*)  "     - Tolerance = ", tolerance
         write(*,*)  "======================================================="
      endif
   else if (trim(sptlSmpType)==trim('rbf')) then
      !! to be added 
      point_count = 0
      if (xstart(1)==1.and.nclxCPL1==1) then
         point_count = point_count + xsize(2)*xsize(3)
      endif
      if (xend(1)==nx.and.nclxCPLn==1) then
         point_count = point_count + xsize(2)*xsize(3)
      endif

      if (xstart(2)==1.and.nclyCPL1==1) then
         point_count = point_count + xsize(1)*xsize(3)
      endif
      if (xend(2)==ny.and.nclyCPLn==1) then
         point_count = point_count + xsize(1)*xsize(3)
      endif

      if (xstart(3)==1.and.nclzCPL1==1) then
         point_count = point_count + xsize(2)*xsize(1)
      endif
      if (xend(3)==nz.and.nclzCPLn==1) then
         point_count = point_count + xsize(2)*xsize(1)
      endif


      grdSpce(1) = xlx/(nx-1)
      grdSpce(2) = yly/(ny-1)
      grdSpce(3) = zlz/(nz-1)
      allocate (point_x(point_count),point_y(point_count),point_z(point_count))
      cnt = 0
      ! Left BC 
      if (xstart(1)==1.and.nclxCPL1==1) then
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               cnt = cnt + 1
               point_x(cnt) = 0.0_mytype
               point_y(cnt) = yp(j+xstart(2)-1)
               point_z(cnt) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)
            enddo
         enddo
      endif
      ! Right BC 
      if (xend(1)==nx.and.nclxCPLn==1) then
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               cnt = cnt + 1
               point_x(cnt) = xlx
               point_y(cnt) = yp(j+xstart(2)-1)
               point_z(cnt) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)
            enddo
         enddo
      endif
      ! Bottom BC 
      if (xstart(2)==1.and.nclyCPL1==1) then
         do k = 1, xsize(3)
            do i = 1, xsize(1)
               cnt = cnt + 1
               point_x(cnt) = real((i+xstart(1)-1-1),mytype)*grdSpce(1)
               point_y(cnt) = 0.0_mytype
               point_z(cnt) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)
            enddo
         enddo
      endif
      ! Top BC 
      if (xend(2)==ny.and.nclyCPLn==1) then
         do k = 1, xsize(3)
            do i = 1, xsize(1)
               cnt = cnt + 1
               point_x(cnt) = real((i+xstart(1)-1-1),mytype)*grdSpce(1)
               point_y(cnt) = yly
               point_z(cnt) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)
            enddo
         enddo
      endif
      ! Back BC  z = 0
      if (xstart(3)==1.and.nclzCPL1==1) then
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               cnt = cnt + 1
               point_x(cnt) = real((i+xstart(1)-1-1),mytype)*grdSpce(1)
               point_y(cnt) = yp(j+xstart(2)-1)
               point_z(cnt) = 0.0_mytype
            enddo
         enddo
      endif
      ! Front BC z= zlz 
      if (xend(3)==nz.and.nclzCPLn==1) then
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               cnt = cnt + 1
               point_x(cnt) = real((i+xstart(1)-1-1),mytype)*grdSpce(1)
               point_y(cnt) = yp(j+xstart(2)-1)
               point_z(cnt) = zlz
            enddo
         enddo
      endif

      point_x = point_x * spatialScale + dataOrgShft(1)
      point_y = point_y * spatialScale + dataOrgShft(2)
      point_z = point_z * spatialScale + dataOrgShft(3)


      write(*,*) "*****************************************"
      write(*,*) point_count , cnt ,xsize(1), nclxCPL1,xsize(1),xsize(2),xsize(3)
      write(*,*) "*****************************************"
      


      call mui_create_sampler_gauss_3d_vector(spatial_sampler,rSpatialSamp,sigmaSpatialSamp)
      ! call mui_create_sampler_rbf_3d_f(spatial_sampler,rSpatialSamp,point_x,point_y,point_z &
      !                                   point_count,basisFunc,conservative,smoothFunc, &
      !                                   generateMatrix,TRIM(makedirMString),cutoff,cgSolveTol, &
      !                                   cgMaxIter,pouSize,preconditioner,MUI_COMM_WORLD)
      if (nrank==0) then 
         Write(*,*) trim(sptlSmpType), " to be added"
         write(*,*) "The available spatial samplers are exact, RBF"
      endif
      
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
      mui_fetch => mui_fetch_exact_exact_3d_vector
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('exact')) then
      mui_fetch => mui_fetch_gauss_exact_3d_vector
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('exact')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_exact_3d_vector


   elseif (trim(sptlSmpType)==trim('exact') .and. trim(tmpSmpType)==trim('gauss') ) then
      mui_fetch => mui_fetch_exact_gauss_3d_vector
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('gauss')) then
      mui_fetch => mui_fetch_gauss_gauss_3d_vector
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('gauss')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_gauss_3d_vector

   else if (trim(sptlSmpType)==trim('rbf') .and. trim(tmpSmpType)==trim('gauss')) then
      mui_fetch => mui_fetch_gauss_gauss_3d_vector

   elseif (trim(sptlSmpType)==trim('exact') .and. trim(tmpSmpType)==trim('mean') ) then
      mui_fetch => mui_fetch_exact_mean_3d_vector
   else if (trim(sptlSmpType)==trim('gauss') .and. trim(tmpSmpType)==trim('mean')) then
      mui_fetch => mui_fetch_gauss_mean_3d_vector
   else if (trim(sptlSmpType)==trim('linear') .and. trim(tmpSmpType)==trim('mean')) then   
      mui_fetch => mui_fetch_pseudo_n2_linear_mean_3d_vector



   else
      if (nrank==0) then 
         write(*,*) "The selected spatial and temporal sampler types have not impelemnted yet"
      endif
      stop
   endif

   call mui_create_temporal_sampler_exact_3d_f(temporal_sync, tolerance)
   call mui_create_sampler_exact_3d_f(spatial_sync, tolerance)

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
