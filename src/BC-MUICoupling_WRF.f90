!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module MUIcoupledBC_WRF

  use decomp_2d
  use variables
  use param
  use tbl, only :blasius
  use abl, only : inflow, outflow

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_WRF, boundary_conditions_WRF, postprocess_WRF, &
  visu_WRF, visu_WRF_init

contains

  subroutine init_WRF (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use param , only : zptwofive
    use MPI

   !  coupling variables 
    use iso_c_binding
    use mui_3d_f
    use mui_general_f

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,ii,is,it,code

    integer, dimension (:), allocatable :: seed




    if (iscalar==1) then

       phi1(:,:,:,:) = zptwofive !change as much as you want
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi1(:,xsize(2),:,:) = zptwofive
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero

   !  a blasius profile is created in ecoule and then duplicated for the all domain
    
    call blasius()
   !  if (nclx1==2) then ! Use the orignial TBL inlet conditions 
   !    call blasius()   
   !  elseif (nclx1==3) then  ! Use the orignial get inlet conditions from MUI interface
   !    call recieveMUIBC(bxx1,bxy1,bxz1,0.0_mytype) ! last argument is the x-location
   !  else
   !    print *, "ERROR: nclx1 = ", nclx1, " is wrong BC for this simulation type"
   !    stop
   !  endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_WRF
  !********************************************************************
  subroutine boundary_conditions_WRF (ux,uy,uz,phi)

    use navier, only : tbl_flrt
    use param , only : zero, zptwofive
    

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um,point_x,point_y,point_z,fetchedVal
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx
    real(mytype) :: ut1, ut2, ut3, ut4, ut5, ut6,totalFaceArea
    real(mytype) :: utt1,utt2,utt3,utt4,utt5,utt6,mssFlwDiffFace(6),mssFlwDiff
    real(mytype) :: faceNormalx(2),faceNormaly(2),faceNormalz(2)
    real :: recTime1,recTime2,recTime
    integer :: i, j, k, is,code, reset_log=1
   
    faceNormalx(1) = 1; faceNormalx(2) = -1
    faceNormaly(1) = 1; faceNormaly(2) = -1
    faceNormalz(1) = 1; faceNormalz(2) = -1
    ! write (*,*) "====================================================="
    ! write (*,*) "============== Xcompcat3d is wiating for data ========"
    ! write (*,*) "====================================================="
   !  point_x = -315487.65625!real(i,mytype) !long_u(i,j) 
   !  point_y = -612657.9375  !real(j,mytype) !lat_u(i,j)
   !  point_z = 44.483718872070312 !real(k,mytype)
   

   !  call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, point_y, &
   !  point_z, 40.0_mytype, spatial_sampler, temporal_sampler, fetchedVal)

   !  write (*,*) "====================================================="
   !  write (*,*) "============== Xcompcat3d is received ==============",fetchedVal
   !  write (*,*) "====================================================="

    bxx1 = 0.0_mytype ; bxy1 = 0.0_mytype; bxz1 = 0.0_mytype
    bxxn = 0.0_mytype ; bxyn = 0.0_mytype; bxzn = 0.0_mytype

    byx1 = 0.0_mytype ; byy1 = 0.0_mytype; byz1 = 0.0_mytype
    byxn = 0.0_mytype ; byyn = 0.0_mytype; byzn = 0.0_mytype

    bzx1 = 0.0_mytype ; bzy1 = 0.0_mytype; bzz1 = 0.0_mytype
    bzxn = 0.0_mytype ; bzyn = 0.0_mytype; bzzn = 0.0_mytype

   if (nrank ==0 ) write(*,*) "Fetcheing at ", t * timeScale 
   call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "tempvalue"//c_null_char,0.0_mytype,0.0_mytype,0.0_mytype,real(5.123,mytype))
   call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, t * timeScale )

   call cpu_time(recTime1)
    if (nclxCPL1 .eq. 1) then 
      call recieve_WRF(bxx1,bxy1,bxz1,0.0_mytype,1)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
    else
      call blasius()
    endif

    if (nclxCPLn .eq. 1) then 
      call recieve_WRF(bxxn,bxyn,bxzn,   xlx    ,1)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
    else 
      call outflow (ux,uy,uz,phi)
      
    endif

    do k=1,xsize(3)
      do j=1,xsize(2)
         ux1(1,j,k) = bxx1(j,k) 
         ux1(nx,j,k) = bxxn(j,k)
         enddo
    enddo



    if (xstart(2)==1) then !! CPUs at the bottom
      if (nclyCPL1 .eq. 1)  call recieve_WRF(byx1,byy1,byz1,0.0_mytype,2)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
      do k=1,xsize(3)
         do i=1,xsize(1)
            uy1(i,1,k) = byy1(i,k)
            enddo
      enddo

    endif

    if (xend(2)==ny) then !! CPUs at the top      
      if (nclyCPLn .eq. 1) then
         call recieve_WRF(byxn,byyn,byzn,   yly    ,2)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
      else
        do k = 1, xsize(3)
            do i = 1, xsize(1)
              byxn(i, k) = ux(i, xsize(2) - 1, k)
              byyn(i, k) = uy(i, xsize(2) - 1, k)
              byzn(i, k) = uz(i, xsize(2) - 1, k)
            enddo
        enddo

      endif
      do k=1,xsize(3)
         do i=1,xsize(1)
            uy1(i,xsize(2),k) = byyn(i,k)
         enddo
      enddo

    endif

    if (xstart(3)==1) then !! CPUs at the front
      if (nclzCPL1 .eq. 1) call recieve_WRF(bzx1,bzy1,bzz1,0.0_mytype,3)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
      do j=1,xsize(2)
         do i=1,xsize(1)
            uz1(i,j,1) = bzz1(i,j)
            enddo
      enddo

    endif 

    if (xend(3)==nz) then !! CPUs at the back
      if (nclzCPLn .eq. 1) call recieve_WRF(bzxn,bzyn,bzzn,   zlz    ,3)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
      do j=1,xsize(2)
         do i=1,xsize(1)
            uz1(i,j,xsize(3)) = bzzn(i,j)
         enddo
      enddo
    endif 
    call cpu_time(recTime2)
    if (nrank==0) print *, "MUI: Receiving CPU time ", recTime2-recTime1
    !! Mass flow rate correction
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)

    ut1=zero;utt1=zero
    if (ystart(1)==1) then !! CPUs at the inlet
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut1=ut1+(yp(j+1)-yp(j))*(ux2(1,j+1,k)-half*(ux2(1,j+1,k)-ux2(1,j,k)))
        enddo
      enddo
      ! ut1=ut1/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(ut1,utt1,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    utt1=utt1/real(nz,mytype) * zlz * faceNormalx(1) !! Volume flow rate 
    ! Flow rate at the outlet
    ut2=zero;utt2=zero
    if (yend(1)==nx) then !! CPUs at the outlet
      do k=1,ysize(3)
        do j=1,ysize(2)-1
          ut2=ut2+(yp(j+1)-yp(j))*(ux2(ysize(1),j+1,k)-half*(ux2(ysize(1),j+1,k)-ux2(ysize(1),j,k)))
        enddo
      enddo
      ! ut2=ut2/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(ut2,utt2,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    utt2=utt2/real(nz,mytype)*zlz * faceNormalx(2) !! Volume flow rate

     ! Flow rate at the top and bottom
    ut3=zero
    ut4=zero
    do k=1,ysize(3)
      do i=1,ysize(1)
        ut3=ut3+uy2(i,1,k)
        ut4=ut4+uy2(i,ny,k)
      enddo
    enddo
    call MPI_ALLREDUCE(ut3,utt3,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    call MPI_ALLREDUCE(ut4,utt4,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    utt3=utt3/(real(nx*nz,mytype))*xlx*zlz  * faceNormaly(1) !!! Volume flow rate
    utt4=utt4/(real(nx*nz,mytype))*xlx*zlz  * faceNormaly(2)!!! Volume flow rate

    ut5=zero;utt5=zero
    if (ystart(3)==1) then !! CPUs at the front
      do i=1,ysize( 1)
        do j=1,ysize(2)-1
          ut5=ut5+(yp(j+1)-yp(j))*(uz2(i,j+1,1)-half*(uz2(i,j+1,1)-uz2(i,j,1)))
        enddo
      enddo
      ! ut1=ut1/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(ut5,utt5,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    utt5=utt5/real(nx,mytype) * xlx * faceNormalz(1)  !! Volume flow rate 
    ! Flow rate at the outlet
    ut6=zero;utt6=zero
    if (yend(1)==nx) then !! CPUs at the back
      do i=1,ysize(1)
        do j=1,ysize(2)-1
          ut6=ut6+(yp(j+1)-yp(j))*(uz2(i,j+1,ysize(3))-half*(uz2(i,j+1,ysize(3))-uz2(i,j,ysize(3))))
        enddo
      enddo
      ! ut2=ut2/real(ysize(3),mytype)
    endif
    call MPI_ALLREDUCE(ut6,utt6,1,real_type,MPI_SUM,MUI_COMM_WORLD,code)
    utt6=utt6/real(nx,mytype) * xlx  * faceNormalz(2) !! Volume flow rate 
  
   if ((nrank==0).and.(mod(itime,ilist)==0)) then
      print *, "WRF mass left/right" ,utt1,utt2
      print *, "WRF mass bottm/top" ,utt3,utt4
      print *, "WRF mass front/back" ,utt5,utt6

      print *, "Mass diff ", utt1+utt2+utt3+utt4+utt5+utt6
   endif
   mssFlwDiff = utt1+utt2+utt3+utt4+utt5+utt6
   totalFaceArea = 2*(zlz * yly + xlx * yly)  + zlz * xlx
   mssFlwDiffFace(1) =  mssFlwDiff / totalFaceArea  ! left face
   mssFlwDiffFace(2) =  mssFlwDiff / totalFaceArea  ! right face

   mssFlwDiffFace(3) =  mssFlwDiff / totalFaceArea  ! bottom face
   mssFlwDiffFace(4) =  mssFlwDiff / totalFaceArea  ! top face

   mssFlwDiffFace(5) =  mssFlwDiff / totalFaceArea  ! front face
   mssFlwDiffFace(6) =  mssFlwDiff / totalFaceArea  ! back face 

  !  !  print *, "totalFaceArea = " , totalFaceArea, mssFlwDiffFace
  ! !  bxx1 = bxx1 - mssFlwDiffFace(1) * faceNormalx(1) 
  ! !  bxxn = bxxn - mssFlwDiffFace(2) * faceNormalx(2) 
    
  ! !  ! if (xstart(2)==1) byy1 = byy1 - mssFlwDiffFace(3) * faceNormaly(1) 
  ! !  if (xend(2)==ny)  byyn = byyn - mssFlwDiffFace(4) * faceNormaly(2) 

  ! !  if (xstart(3)==1) bzz1 = bzz1 - mssFlwDiffFace(5) * faceNormalz(1) 
  ! !  if (xend(3)==nz)  bzzn = bzzn - mssFlwDiffFace(6) * faceNormalz(2) 

   ux(1,:,:)=bxx1(:,:)
   ux(nx,:,:)=bxxn(:,:)

   if (xstart(2)==1) uy(:,1,:)=byy1(:,:)
   if (xend(2)==ny) uy(:,xsize(2),:)=byyn(:,:)

   if (xstart(3)==1) uz(:,:,1)=bzz1(:,:)
   if (xend(3)==nz) uz(:,:,xsize(3))=bzzn(:,:)

   

   if ( iTime>=ifirst+nForget) then
      if (nrank==0) write(*,*) "MUI: Forgeting log at time stamp of ", t-dt*nForget
      call mui_forget_upper_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, &
               real(t-dt*nForget,c_double),reset_log)
   endif

   

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet
   !  if (nclx1==2) then ! Use the orignial TBL inlet conditions       
   !    call blasius()
 

   !  endif



   ! !  INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
   !  if (iscalar==1) then
   !     do k=1,xsize(3)
   !        do j=1,xsize(2)
   !           phi(1,:,:,:)=zptwofive
   !           if ((xstart(2)==1)) then
   !              phi(:,1,:,:) = one
   !           endif
   !           if ((xend(2)==ny)) THEN
   !              phi(:,xsize(2),:,:) = zptwofive
   !           endif
   !        enddo
   !     enddo
   !  endif

   !  udx=one/dx
   !  udy=one/dy
   !  udz=one/dz
   !  uddx=half/dx
   !  uddy=half/dy
   !  uddz=half/dz
   !  if (nclxn==2) then !OUTFLOW based on a 1D convection equation
   !    do k=1,xsize(3)
   !       do j=1,xsize(2)

   !          cx=ux(nx,j,k)*gdt(itr)*udx

   !          if (cx<zero) cx=zero
   !          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
   !          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
   !          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
   !          if (iscalar==1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
   !          enddo
   !    enddo


   !  endif



    !! Bottom Boundary
   !  if (ncly1 == 2) then
   !    do k = 1, xsize(3)
   !      do i = 1, xsize(1)
   !        byx1(i, k) = zero
   !        byy1(i, k) = zero
   !        byz1(i, k) = zero
   !      enddo
   !    enddo
   !  endif
   !  !! Top Boundary
   !  if (nclyn == 2) then
   !     do k = 1, xsize(3)
   !        do i = 1, xsize(1)
   !           byxn(i, k) = ux(i, xsize(2) - 1, k)
   !           byyn(i, k) = uy(i, xsize(2) - 1, k)
   !           byzn(i, k) = uz(i, xsize(2) - 1, k)
   !        enddo
   !     enddo

   !  endif

    !SCALAR   
    if (itimescheme/=7) then
    if (iscalar/=0) then
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
    endif
    endif
    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)

    

    return
  end subroutine boundary_conditions_WRF

  !********************************************************************
  !********************************************************************
  subroutine recieve_WRF(bx,by,bz,Loc,dimIndx) ! last argument is the x-location 
#ifdef MUI_COUPLING
  ! Coupling varaibles 
    use iso_c_binding
    use mui_3d_f
    use mui_general_f
#endif
    use decomp_2d_io
    use MPI
    use variables, only : nForget
    

    implicit none

    real(mytype) :: Loc,yLoc,x, y, z
    real(mytype),dimension(:,:), intent(inout):: bx,by,bz

    real(mytype) :: fetch_result_3d,dataLoc(3), grdSpce(3)
    real(mytype) :: cplTime,point_y,point_z
    integer :: i, j, k,dimIndx
    grdSpce(1) = xlx/(nx-1)
    grdSpce(2) = yly/(ny-1)
    grdSpce(3) = zlz/(nz-1)
    
    cplTime = t * timeScale 
    
    ! Obtain the vertical BC, i.e. plane normal to x-direction
    
    if (dimIndx==1) then 
      ! if (nclx1==4 .or. nclxn==4) then 
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               dataLoc(1) = Loc
               dataLoc(2) = yp(j+xstart(2)-1) 
               dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)print *, "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "ux"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bx(j,k))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uy"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, by(j,k))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uz"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bz(j,k))
               ! print*, "Xcompact : ", cplTime, dataLoc(1), dataLoc(3), dataLoc(2), bx(j,k), by(j,k), bz(j,k)

            end do
         end do
      ! endif
   endif

   ! Obtain the Horizontal BC, i.e. plane normal to y-direction
   if (dimIndx==2) then 
      ! if (ncly1==4 .or. nclyn==4) then 
         do i = 1, xsize(1)
            do k = 1, xsize(3)
               dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1) 
               dataLoc(2) = Loc 
               dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)print *, "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z

               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "ux"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bx(i,k))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uy"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, by(i,k))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uz"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bz(i,k))

            end do
         end do
      ! endif
   endif

   ! Obtain the Horizontal BC, i.e. plane normal to y-direction
   if (dimIndx==3) then 
      ! if (nclz1==4 .or. nclzn==4) then 
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1) 
               dataLoc(2) = yp(j+xstart(2)-1) 
               dataLoc(3) = Loc 

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)print *, "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z

               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "ux"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bx(i,j))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uy"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, by(i,j))
               call mui_fetch(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uz"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, bz(i,j))

            end do
         end do
      ! endif
   endif


   

   
   
   !  stop
end subroutine recieve_WRF


  !############################################################################
  subroutine postprocess_WRF(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename

  end subroutine postprocess_WRF

  subroutine visu_WRF_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_WRF_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_MUIBC
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TBL-specific visualization
  !!
  !############################################################################
  subroutine visu_WRF(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", num, flush=.true.) ! Reusing temporary array, force flush

  end subroutine visu_WRF

end module MUIcoupledBC_WRF
