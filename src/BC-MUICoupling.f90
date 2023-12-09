!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module MUIcoupledBC

  use decomp_2d
  use variables
  use param
  use tbl, only :blasius

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_MUIBC, boundary_conditions_MUIBC, postprocess_MUIBC, &
  visu_MUIBC, visu_MUIBC_init, pushMUI

contains

  subroutine init_MUIBC (ux1,uy1,uz1,ep1,phi1)

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

    !a blasius profile is created in ecoule and then duplicated for the all domain
    if (nclx1==2) then ! Use the orignial TBL inlet conditions 
      call blasius()   
    elseif (nclx1==3) then  ! Use the orignial get inlet conditions from MUI interface
      call recieveMUIBC()
    else
      print *, "ERROR: nclx1 = ", nclx1, " is wrong BC for this simulation type"
      stop
    endif

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
  end subroutine init_MUIBC
  !********************************************************************
  subroutine boundary_conditions_MUIBC (ux,uy,uz,phi)

    use navier, only : tbl_flrt
    use param , only : zero, zptwofive

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: x, y, z, um
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx

    integer :: i, j, k, is
      
    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet
    if (nclx1==2) then ! Use the orignial TBL inlet conditions 
      
      call blasius()
    elseif (nclx1==3) then ! Use the orignial get inlet conditions from MUI interface 
      call recieveMUIBC()
    else
      print *, "ERROR: nclx1 = ", nclx1, " is wrong BC for this simulation type"
      stop
    endif
    
   !  
    
    !INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
    if (iscalar==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(1,:,:,:)=zptwofive
             if ((xstart(2)==1)) then
                phi(:,1,:,:) = one
             endif
             if ((xend(2)==ny)) THEN
                phi(:,xsize(2),:,:) = zptwofive
             endif
          enddo
       enddo
    endif

    !OUTFLOW based on a 1D convection equation

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx<zero) cx=zero
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
          if (iscalar==1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

    !! Bottom Boundary
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
    endif
    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k)
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    endif

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
  end subroutine boundary_conditions_MUIBC

  !********************************************************************
  !********************************************************************
  subroutine recieveMUIBC()
#ifdef MUI_COUPLING
  ! Coupling varaibles 
    use iso_c_binding
    use mui_3d_f
    use mui_general_f
#endif
    use decomp_2d_io
    use MPI
    use param, only : zero, zptwo, zpeight, one, nine
    use dbg_schemes, only: exp_prec, sqrt_prec

    implicit none

    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
    real(mytype) :: delta_int, delta_eta, eps_eta

    real(mytype) :: x, y, z,fetch_result_3d,point_x,point_y,point_z, grdSpce(3)
    integer :: i, j, k, is
    grdSpce(1) = xlx/(nx-1)
    grdSpce(2) = yly/(ny-1)
    grdSpce(3) = zlz/(nz-1)


    print *, "Fetched 3D interface values at time  ",t
    do k = 1, xsize(3)
      do j = 1, xsize(2)
         point_x = 0.0_mytype + dataOrgShft(1)
         point_y = yp(j+xstart(2)-1) + dataOrgShft(2)
         point_z = real((k+xstart(3)-1),mytype)*grdSpce(3) + dataOrgShft(3)
         
         call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "ux"//c_null_char, point_x, point_y, &
         point_z, t, spatial_sampler, temporal_sampler, bxx1(j,k))

         call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uy"//c_null_char, point_x, point_y, &
         point_z, t, spatial_sampler, temporal_sampler, bxy1(j,k))

         call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uz"//c_null_char, point_x, point_y, &
         point_z, t, spatial_sampler, temporal_sampler, bxz1(j,k))

      end do
   end do

    
end subroutine recieveMUIBC

!********************************************************************
!********************************************************************
subroutine pushMUI(ux1, uy1, uz1)
#ifdef MUI_COUPLING
   ! Coupling varaibles 
      use iso_c_binding
      use mui_3d_f
      use mui_general_f
#endif
      use decomp_2d_io
      use MPI
      use param, only : zero, zptwo, zpeight, one, nine
      use dbg_schemes, only: exp_prec, sqrt_prec

      implicit none
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
      real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
      real(mytype) :: delta_int, delta_eta, eps_eta

      real(mytype) :: x, y, z,fetch_result_3d,point_x,point_y,point_z, grdSpce(3)
      integer :: i, j, k, is,iGrp
      integer, allocatable,dimension(:,:) :: groupVortLocal  ! the group is defined with votecies (ix1,ix2,jy1,jy2,kz1,kz2)
      ! integer, allocatable,dimension(:,:,:) :: cornLoc,corn
      groupNumb=1
      allocate(groupVortLocal(groupNumb,6))

      !!! DO NOT dx,dy,dz used in the code as their values depends on the BC not the actual grid spacing.
      grdSpce(1) = xlx/(nx-1)
      grdSpce(2) = yly/(ny-1)
      grdSpce(3) = zlz/(nz-1)

      !!! Define the local index of points of each  group 
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

         enddo
         ! If a  group does not have a data in this local domain, set it to zero and -1. Use -1 to avoid entering the loop. 
         if (minval(groupVortLocal(iGrp,:)) ==0)  groupVortLocal(iGrp,:) = [0, 0, 0, -1, -1, -1]  
      enddo 
      ! write(*,*) nrank, groupVortLocal
      !/ Push data to MUI interface 
      do iGrp = 1, groupNumb
         do k=groupVortLocal(iGrp,5),groupVortLocal(iGrp,6)
            do j=groupVortLocal(iGrp,3),groupVortLocal(iGrp,4)
               do i=groupVortLocal(iGrp,1),groupVortLocal(iGrp,2)
                  point_x = real((i+xstart(1)-1),mytype)*grdSpce(1)  + dataOrgShft(1)
                  point_y = yp(j+xstart(2)-1) + dataOrgShft(2)
                  point_z = real((k+xstart(3)-1),mytype)*grdSpce(3) + dataOrgShft(3)

                  ! print *, "MUI Sending  at location", point_x, point_y, point_z
                  
                  call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, &
                  point_y, point_z, ux1(i,j,k))

                  call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, point_x, &
                  point_y, point_z, uy1(i,j,k))

                  call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, point_x, &
                  point_y, point_z, uz1(i,j,k))

                  ! print *, "At ", point_x, point_y, point_z, "Xcompact 3d Pushed ux", ux1(i,j,k)

               enddo
            enddo
         enddo
         call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, T)
      enddo

      
      ! print *, nrank,"Group Vortext local  ", groupVortLocal(iGrp,:)
      
      ! do k = 1, xsize(3)
      ! do j = 1, xsize(2)
      !    point_x = 0.0_mytype
      !    point_y = yp(j+xstart(2)-1)
      !    point_z = real((k+xstart(3)-1-1),mytype)*dz
         
      !    call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "ux"//c_null_char, point_x, point_y, &
      !    point_z, t, spatial_sampler, temporal_sampler, bxx1(j,k))

      !    call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uy"//c_null_char, point_x, point_y, &
      !    point_z, t, spatial_sampler, temporal_sampler, bxy1(j,k))

      !    call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(MUIBC_ID(1))%ptr, "uz"//c_null_char, point_x, point_y, &
      !    point_z, t, spatial_sampler, temporal_sampler, bxz1(j,k))

      !    end do
      ! end do

      
end subroutine pushMUI
  !############################################################################
  subroutine postprocess_MUIBC(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename

  end subroutine postprocess_MUIBC

  subroutine visu_MUIBC_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_MUIBC_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_MUIBC
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TBL-specific visualization
  !!
  !############################################################################
  subroutine visu_MUIBC(ux1, uy1, uz1, pp3, phi1, ep1, num)

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

  end subroutine visu_MUIBC

end module MUIcoupledBC
