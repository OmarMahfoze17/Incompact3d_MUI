!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module MUIcoupledBC_WRF

  use decomp_2d
  use variables
  use param
  use tbl, only :blasius
  use abl, only : inflow, outflow
  use MUIcoupledBC, only : pushMUI

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

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1,noiseX,noiseY,noiseZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct, grdSpce(3),dataLoc(3),cplTime
    real(mytype) :: cx0,cy0,cz0,hg,lg,temp(3)
    integer :: k,j,i,ierror,ii,is,it,code

    integer, dimension (:), allocatable :: seed

    call random_number(noiseX)
    call random_number(noisey)
    call random_number(noisez)
    um=0.5*(u1+u2)


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
   !    write(*,*) "ERROR: nclx1 = ", nclx1, " is wrong BC for this simulation type"
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

    if (iin.eq.3) then
        grdSpce(1) = xlx/(nx-1)
        grdSpce(2) = yly/(ny-1)
        grdSpce(3) = zlz/(nz-1)

        cplTime = 0.0
      do k = 1, xsize(3)
        write(*,*) "xcomcapt initaite  ", k ,'test', xsize(2),xsize(1)
        do j = 1, xsize(2)
          
          do i=1,xsize(1)
            dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1) 
            dataLoc(2) = yp(j+xstart(2)-1) 
            dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)

            dataLoc = dataLoc * spatialScale + dataOrgShft
            ! if (Loc==xlx)write(*,*) "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z
            call mui_fetch(uniface_pointers_3d(fetchInterfacefsID)%ptr, "velocity"//c_null_char, dataLoc(1), dataLoc(2), &
            dataLoc(3), cplTime, spatial_sampler, temporal_sampler, temp)

            
            ux1(i,j,k)=temp(1)+(two*noiseX(i,j,k)-one)*init_noise*um
            uy1(i,j,k)=temp(2)+(two*noisey(i,j,k)-one)*init_noise*um
            uz1(i,j,k)=temp(3)+(two*noisez(i,j,k)-one)*init_noise*um
            
            if (k==1 .and. j==1 .and. i==1) write(*,*) "Xcompact initiat: ", dataLoc(1), dataLoc(3), dataLoc(2), ux1(i,j,k),uy1(i,j,k),uz1(i,j,k)
          enddo
        end do
          
      end do
      
      !  uz1=-uz1

    endif
    if (nrank  ==  0) write(*,*) 'xcomact3d : init end ok'

    return
  end subroutine init_WRF
  !********************************************************************
  subroutine boundary_conditions_WRF (ux1,uy1,uz1,phi)

    use navier, only : tbl_flrt
    use param , only : zero, zptwofive
    

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: noiseX,noiseY,noiseZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um,point_x,point_y,point_z,fetchedVal
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx
    real(mytype) :: ut1, ut2, ut3, ut4, ut5, ut6,totalFaceArea,massFlWTotal
    real(mytype) :: utt1,utt2,utt3,utt4,utt5,utt6,mssFlwDiffFace(6),mssFlwDiff
    real(mytype) :: faceNormalx(2),faceNormaly(2),faceNormalz(2)
    real :: recTime1,recTime2,recTime
    integer :: i, j, k, is,code, reset_log=1
   
    faceNormalx(1) = 1; faceNormalx(2) = -1
    faceNormaly(1) = 1; faceNormaly(2) = -1
    faceNormalz(1) = 1; faceNormalz(2) = -1

    call random_number(noiseX)
    call random_number(noisey)
    call random_number(noisez)
    um=0.5*(u1+u2)

   !  write (*,*) "====================================================="
   !  write (*,*) "============== Xcompcat3d is received ==============",fetchedVal
   !  write (*,*) "====================================================="

    bxx1 = 0.0_mytype ; bxy1 = 0.0_mytype; bxz1 = 0.0_mytype
    bxxn = 0.0_mytype ; bxyn = 0.0_mytype; bxzn = 0.0_mytype

    byx1 = 0.0_mytype ; byy1 = 0.0_mytype; byz1 = 0.0_mytype
    byxn = 0.0_mytype ; byyn = 0.0_mytype; byzn = 0.0_mytype

    bzx1 = 0.0_mytype ; bzy1 = 0.0_mytype; bzz1 = 0.0_mytype
    bzxn = 0.0_mytype ; bzyn = 0.0_mytype; bzzn = 0.0_mytype

   if (nrank ==0 ) write(*,*) trim(domainName),":Fetcheing at ", t * timeScale 
  !  call mui_push_3d_f(uniface_pointers_3d(fetchInterfacefsID)%ptr, "tempvalue"//c_null_char, &
  !         0.0_mytype,0.0_mytype,0.0_mytype,real(5.123,mytype))
   call mui_commit_3d_f(uniface_pointers_3d(fetchInterfacefsID)%ptr, t * timeScale )

  !  call mui_barrier_3d_f(uniface_pointers_3d(fetchInterfacefsID)%ptr,t * timeScale)
   
   call cpu_time(recTime1)
    if (nclxCPL1 .eq. 1) then 
      call recieve_WRF(bxx1,bxy1,bxz1,0.0_mytype,1)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
    else
      call blasius()
    endif

    if (nclxCPLn .eq. 1) then 
      call recieve_WRF(bxxn,bxyn,bxzn,   xlx    ,1)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
    else 
      call outflow_x (ux1,uy1,uz1,phi)      
    endif


    bxx1(:,:) = bxx1(:,:) +(two*noiseX(1,:,:)-one)*inflow_noise*um
    bxy1(:,:) = bxy1(:,:) +(two*noiseY(1,:,:)-one)*inflow_noise*um
    bxz1(:,:) = bxz1(:,:) +(two*noiseZ(1,:,:)-one)*inflow_noise*um
    bxxn(:,:) = bxxn(:,:) +(two*noiseX(nx,:,:)-one)*inflow_noise*um
    bxyn(:,:) = bxyn(:,:) +(two*noiseY(nx,:,:)-one)*inflow_noise*um
    bxzn(:,:) = bxzn(:,:) +(two*noiseZ(nx,:,:)-one)*inflow_noise*um

    if (xstart(2)==1) then !! CPUs at the bottom
      if (nclyCPL1 .eq. 1)  call recieve_WRF(byx1,byy1,byz1,0.0_mytype,2)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
      ! byx1 = byx1 + (two*noiseX(:,1,:)-one)*inflow_noise*um
      ! byy1 = byy1 + (two*noiseY(:,1,:)-one)*inflow_noise*um
      ! byz1 = byz1 + (two*noiseZ(:,1,:)-one)*inflow_noise*um

      byx1 = byx1
      byy1 = byy1
      byz1 = byz1
    endif

    if (xend(2)==ny) then !! CPUs at the top      
      if (nclyCPLn .eq. 1) then
         call recieve_WRF(byxn,byyn,byzn,   yly    ,2)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
        !  byxn = byxn + (two*noiseX(:,xsize(2),:)-one)*inflow_noise*um
        !  byyn = byyn + (two*noiseY(:,xsize(2),:)-one)*inflow_noise*um
        !  byzn = byzn + (two*noiseZ(:,xsize(2),:)-one)*inflow_noise*um

         byxn = byxn
         byyn = byyn
         byzn = byzn

      else
      byxn(:,:) = ux1(:,xsize(2) - 1,:)
      byyn(:,:) = uy1(:,xsize(2) - 1,:)
      byzn(:,:) = uz1(:,xsize(2) - 1,:)
      endif

    endif
    if (nclzCPL1 .eq. 1) then 
      if (xstart(3)==1) then !! CPUs at the front
        call recieve_WRF(bzx1,bzy1,bzz1,0.0_mytype,3)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
        bzx1(:,:) = bzx1(:,:) + (two*noiseX(:,:,1)-one)*inflow_noise*um
        bzy1(:,:) = bzy1(:,:) + (two*noiseY(:,:,1)-one)*inflow_noise*um
        bzz1(:,:) = bzz1(:,:) + (two*noiseZ(:,:,1)-one)*inflow_noise*um
      endif 
    endif
    
    if (nclzCPLn .eq. 1) then 
      if (xend(3)==nz) then !! CPUs at the back        
        call recieve_WRF(bzxn,bzyn,bzzn,   zlz    ,3)   !!!  dimIndx =1,2,3 to fetch data in x,y,z normal plane
        bzxn(:,:) = bzxn(:,:) + (two*noiseX(:,:,xsize(3))-one)*inflow_noise*um
        bzyn(:,:) = bzyn(:,:) + (two*noisey(:,:,xsize(3))-one)*inflow_noise*um
        bzzn(:,:) = bzzn(:,:) + (two*noiseZ(:,:,xsize(3))-one)*inflow_noise*um
      endif 
    else
      call outflow_z (ux1,uy1,uz1,phi)
    endif

    call cpu_time(recTime2)
    if (nrank==0) write(*,*) "MUI: Receiving CPU time ", recTime2-recTime1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Apply BC to the velocity 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    if (xstart(2)==1 .and. ncly1==2) then !! CPUs at the bottom
      ux1(:,1,:) = byx1(:,:)
      uy1(:,1,:) = byy1(:,:)
      uz1(:,1,:) = byz1(:,:)
    endif
    if (xend(2)==ny .and. nclyn==2) then !! CPUs at the top      
      ux1(:,xsize(2),:) = byxn(:,:)
      uy1(:,xsize(2),:) = byyn(:,:)
      uz1(:,xsize(2),:) = byzn(:,:)
    endif

    if (xstart(3)==1 .and. nclz1==2) then !! CPUs at the front
        ux1(:,:,1) = bzx1(:,:)
        uy1(:,:,1) = bzy1(:,:)
        uz1(:,:,1) = bzz1(:,:)
    endif 
    if (xend(3)==nz .and. nclzn==2) then !! CPUs at the back
        ux1(:,:,xsize(3)) = bzxn(:,:)
        uy1(:,:,xsize(3)) = bzyn(:,:)
        uz1(:,:,xsize(3)) = bzzn(:,:)
    endif 
    if (nclx1==2) then 
      ux1(1,:,:) = bxx1(:,:)
      uy1(1,:,:) = bxy1(:,:)
      uz1(1,:,:) = bxz1(:,:)
    endif
    if (nclxn==2) then 
      ux1(nx,:,:) = bxxn(:,:) 
      uy1(nx,:,:) = bxyn(:,:)
      uz1(nx,:,:) = bxzn(:,:)
    endif 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Mass flow rate correction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    if (yend(3)==nz) then !! CPUs at the back
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
      write(*,*) trim(domainName),' ',"WRF mass left/right" ,utt1,utt2
      write(*,*) trim(domainName),' ',"WRF mass bottm/top" ,utt3,utt4
      write(*,*) trim(domainName),' ',"WRF mass front/back" ,utt5,utt6

      write(*,*) trim(domainName),' ',"Mass diff ", utt1+utt2+utt3+utt4+utt5+utt6
   endif

   mssFlwDiff = utt1+utt2+utt3+utt4+utt5+utt6
   massFlWTotal = abs(utt1*nclxCrr1)+abs(utt2*nclxCrrn)+abs(utt3*nclyCrr1)+abs(utt4*nclyCrrn) &
                  +abs(utt5*nclzCrr1)+abs(utt6*nclzCrrn)
   totalFaceArea = zlz*yly*(nclxCrr1+nclxCrrn) + xlx * yly*(nclyCrr1+nclyCrrn) &
                  + zlz * xlx*(nclzCrr1+nclzCrrn)
   
   mssFlwDiffFace(1) =  mssFlwDiff * nclxCrr1   ! left face
   mssFlwDiffFace(2) =  mssFlwDiff * nclxCrrn   ! right face

   mssFlwDiffFace(3) =  mssFlwDiff * nclyCrr1   ! bottom face
   mssFlwDiffFace(4) =  mssFlwDiff * nclyCrrn  ! top face

   mssFlwDiffFace(5) =  mssFlwDiff * nclzCrr1  ! front face
   mssFlwDiffFace(6) =  mssFlwDiff * nclzCrrn  ! back face 

   bxx1 = bxx1 - abs(utt1)/massFlWTotal * mssFlwDiffFace(1) * faceNormalx(1) /(yly*zlz)
   bxxn = bxxn - abs(utt2)/massFlWTotal * mssFlwDiffFace(2) * faceNormalx(2) /(yly*zlz)
    
   if (xstart(2)==1) byy1 = byy1 - abs(utt3)/massFlWTotal * mssFlwDiffFace(3) * faceNormaly(1) /(xlx*zlz)
   if (xend(2)==ny) byyn = byyn - abs(utt4)/massFlWTotal * mssFlwDiffFace(4) * faceNormaly(2) /(xlx*zlz)

   if (xstart(3)==1) bzz1 = bzz1 - abs(utt5)/massFlWTotal * mssFlwDiffFace(5) * faceNormalz(1) /(yly*xlx)
   if (xend(3)==nz)  bzzn = bzzn - abs(utt6)/massFlWTotal * mssFlwDiffFace(6) * faceNormalz(2) /(yly*xlx)

   !!! Reapply the BC after the mass flow rate correction. 
   ux1(1,:,:)=bxx1(:,:)
   ux1(nx,:,:)=bxxn(:,:)
   if (xstart(2)==1) uy1(:,1,:)=byy1(:,:)
   if (xend(2)==ny) uy1(:,xsize(2),:)=byyn(:,:)

   if (xstart(3)==1) uz1(:,:,1)=bzz1(:,:)
   if (xend(3)==nz) uz1(:,:,xsize(3))=bzzn(:,:)

   if ( iTime>=ifirst+nForget) then
      if (nrank==0) write(*,*) trim(domainName), " Forgeting log at time stamp of ", t-dt*nForget
      call mui_forget_upper_3d_f(uniface_pointers_3d(fetchInterfacefsID)%ptr, &
               real(t-dt*nForget,c_double),reset_log)
   endif   

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

    real(mytype) :: Loc,yLoc,x, y, z, temp(3)
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
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               dataLoc(1) = Loc
               dataLoc(2) = yp(j+xstart(2)-1) 
               dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)write(*,*) "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z
               call mui_fetch(uniface_pointers_3d(fetchInterfacefsID)%ptr, "velocity"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, temp)
               bx(j,k) = temp(1)
               by(j,k) = temp(2)
               bz(j,k) = temp(3)
               if (k==2 .and. j==2 .and. nrank==0 .and. Loc==0.0_mytype) write(*,*) trim(domainName), cplTime, dataLoc(1), dataLoc(2), dataLoc(3), bx(j,k), by(j,k), bz(j,k)

            end do
         end do
      ! endif
   endif

   ! Obtain the Horizontal BC, i.e. plane normal to y-direction
   if (dimIndx==2) then 
         do i = 1, xsize(1)
            do k = 1, xsize(3)
               dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1) 
               dataLoc(2) = Loc 
               dataLoc(3) = real((k+xstart(3)-1-1),mytype)*grdSpce(3)

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)write(*,*) "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z

               call mui_fetch(uniface_pointers_3d(fetchInterfacefsID)%ptr, "velocity"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, temp)
               bx(i,k) = temp(1)
               by(i,k) = temp(2)
               bz(i,k) = temp(3)

            end do
         end do
      ! endif
   endif

   ! Obtain the vertical BC, i.e. plane normal to z-direction
   if (dimIndx==3) then 

         do j = 1, xsize(2)
            do i = 1, xsize(1)
               dataLoc(1) = real((i+xstart(1)-1-1),mytype)*grdSpce(1) 
               dataLoc(2) = yp(j+xstart(2)-1) 
               dataLoc(3) = Loc 

               dataLoc = dataLoc * spatialScale + dataOrgShft
               ! if (Loc==xlx)write(*,*) "MUI domain ",trim(domainName)," is Receivinfg at location", point_x, point_y, point_z

               call mui_fetch(uniface_pointers_3d(fetchInterfacefsID)%ptr, "velocity"//c_null_char, dataLoc(1), dataLoc(2), &
               dataLoc(3), cplTime, spatial_sampler, temporal_sampler, temp)
               bx(i,j) = temp(1) 
               by(i,j) = temp(2)
               bz(i,j) = temp(3)


            end do
         end do
      ! endif
   endif

  !  bz = -bz     

   

   
   
   !  stop
   
end subroutine recieve_WRF

!*******************************************************************************
  !
  subroutine outflow_x (ux,uy,uz,phi)
  !
  !*******************************************************************************

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609.
    uxmin=1609.
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
      enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MUI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MUI_COMM_WORLD,code)

    cx=0.5*(uxmax1+uxmin1)*gdt(itr)*udx
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
        if (iscalar.eq.1) then
          phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        endif
      enddo
    enddo

    return
  end subroutine outflow_x 

  !*******************************************************************************
  !
  subroutine outflow_z (ux,uy,uz,phi)
    !
    !*******************************************************************************
  
      USE param
      USE variables
      USE decomp_2d
      USE MPI
  
      implicit none
  
      integer :: i,j,k,code
      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
      real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
      real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cz,uzmin,uzmax,uzmin1,uzmax1
  
      udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz
  
      uzmax=-1609.
      uzmin=1609.
      if (xend(3)==nz) then
        do j=1,xsize(2)
          do i=1,xsize(1)
            if (uz(i,j,xsize(3)-1).gt.uzmax) uzmax=uz(i,j,xsize(3)-1)
            if (uz(i,j,xsize(3)-1).lt.uzmin) uzmin=uz(i,j,xsize(3)-1)
          enddo
        enddo
      endif
  
      call MPI_ALLREDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,MUI_COMM_WORLD,code)
      call MPI_ALLREDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,MUI_COMM_WORLD,code)
      cz=0.5*(uzmax1+uzmin1)*gdt(itr)*udz
      if (xend(3)==nz) then
        do j=1,xsize(2)
          do i=1,xsize(1)
            bzxn(i,j)=ux(i,j,xsize(3))-cz*(ux(i,j,xsize(3))-ux(i,j,xsize(3)-1))
            bzyn(i,j)=uy(i,j,xsize(3))-cz*(uy(i,j,xsize(3))-uy(i,j,xsize(3)-1))
            bzzn(i,j)=uz(i,j,xsize(3))-cz*(uz(i,j,xsize(3))-uz(i,j,xsize(3)-1))
            if (iscalar.eq.1) then
              phi(i,j,xsize(3),:)=phi(i,j,xsize(3),:)-cz*(phi(i,j,xsize(3),:)-phi(i,j,xsize(3)-1,:))
            endif
          enddo
        enddo
      endif
  
      return
    end subroutine outflow_z 

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
