
module local_subroutines
  implicit none

  contains

!********************************************************************
  !********************************************************************
  subroutine blasius(bxx1,bxy1,bxz1)

    implicit none
    integer, parameter :: mytype = KIND(0.0D0)
    integer, parameter ::ny =129, nz=32
    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,f_bl_inf,g_bl_inf
    real(mytype) :: delta_int, delta_eta, eps_eta


    real(mytype) ::  y, dy, yly, re,xnu
    integer :: j, k
    real(mytype), dimension(:), allocatable :: yp
    real(mytype), dimension(ny,nz) :: bxx1,bxy1,bxz1
    allocate(yp(ny))
    yly = 20
    re = 1250.0

    xnu=1.0/re
    call read_yp(yp,ny)
    
    dy=yly/real(ny-1,mytype)
    do k=1,nz
       do j=1,ny

          y=yp(j)

          eta_bl=y*real(4.91,mytype)/9.0

          !OLD POLYNOMIAL FITTING

          delta_eta=0.0
          eps_eta=0.0
          delta_int=0.2

          if (eta_bl>=(real(7.5,mytype)/9.0)) then
             delta_eta=eta_bl-real(7.5,mytype)/9.0
             eta_bl=real(7.5,mytype)/9.0
             eps_eta=real(0.00015,mytype)
          end if

          f_bl=1678.64209592595000_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
               +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
               +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
               +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
               +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
               -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
               -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1

          f_bl=f_bl+(1-dexp(-delta_eta/delta_int))*eps_eta

          if (eta_bl >= (7.15_mytype/9.0)) then
             delta_int=0.8
             delta_eta=eta_bl-7.15_mytype/9.0
             eta_bl   =       7.15_mytype/9.0
             eps_eta  =     0.0005_mytype
          end if

          g_bl=4924.052847797540_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
               +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
               +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
               +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
               +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
               +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
               +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1 ! &

          g_bl=g_bl+(1-dexp(-delta_eta/delta_int))*eps_eta

          x_bl=1.0/(4.91_mytype**2*xnu)

          bxx1(j,k)=f_bl/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity
          bxy1(j,k)=g_bl*sqrt(xnu/x_bl)/1.000546554_mytype
          bxz1(j,k)=0.0
       enddo
    enddo

    ! do j=1,ny

    !   write(*,*) bxx1(j,1)
    ! enddo
    !STORE VALUE F_BL_INF G_BL_INF (ONLY 1.0 MORE TIME)------------------

    y=yly
    eta_bl=y*4.91_mytype/9.0  !The 9 is due to interpolation

    delta_eta=0.0
    eps_eta=0.0
    delta_int=0.2

    if (eta_bl>=(7.5_mytype/9.0)) then
       delta_eta=eta_bl-7.5_mytype/9.0
       eta_bl   =       7.5_mytype/9.0
       eps_eta  =   0.00015_mytype
    end if

    !To assure 1.0 in infinity
    f_bl_inf=1678.6420959259500_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
            +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
            +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
            +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
            +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
            -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
            -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1


    f_bl_inf=f_bl_inf+(1-dexp(-delta_eta/delta_int))*eps_eta
    f_bl_inf=f_bl_inf/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity



    if (eta_bl>= (7.15_mytype/9.0)) then
       delta_int=0.8
       delta_eta=eta_bl-7.15_mytype/9.0
       eta_bl   =       7.15_mytype/9.0
       eps_eta  =     0.0005_mytype
    end if

    g_bl_inf=4924.05284779754_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
            +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
            +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
            +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
            +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
            +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
            +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1


    g_bl_inf=g_bl_inf+(1-dexp(-delta_eta/delta_int))*eps_eta
    g_bl_inf=g_bl_inf/1.000546554_mytype

    return
  end subroutine blasius

!=========================================================================
!=========================================================================
!=========================================================================
subroutine read_yp(yp,ny)
  implicit none

  ! Parameters
  integer, intent(in) :: ny

  ! Input/output variables
  real(8), dimension(ny), intent(out) :: yp

  ! Local variables
  integer :: unit_number =10
  integer :: num_points
  integer :: iostatus
  real(8) :: value

  ! Initialize
  ! allocate(yp(ny))
  yp = 0.0
  num_points = 0

  ! Open the file
  open(unit=unit_number, file='./yp.dat', status='old', action='read', iostat=iostatus)
  if (iostatus /= 0) then
    write(*, *) 'Error opening file: ./yp.dat'
    return
  end if

  ! Read data from the file
  do
    read(unit_number, *, iostat=iostatus) value
    if (iostatus /= 0) exit
    num_points = num_points + 1
    if (num_points > ny) then
      write(*, *) 'Too many data points in the file.'
      return
    end if
    yp(num_points) = value
  end do

  ! Close the file
  close(unit_number)

  ! Resize the array to the actual number of data points
  ! allocate(yp(num_points))

end subroutine read_yp

  end module local_subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  use iso_c_binding
  use, intrinsic:: iso_fortran_env, only: stdout=>output_unit, stdin=>input_unit, stderr=>error_unit
  use mui_3d_f
  use mui_general_f
  use local_subroutines


  implicit none
  integer, parameter :: mytype = KIND(0.0D0)
  integer, parameter ::ny =129, nz=32
  real(mytype), dimension(:), allocatable :: yp,zp
  !Local variables
  character(len=1024) :: domainName = 'dummySendReceive'
  character(len=1024) :: interfaceName = 'interface'
  integer(c_int) :: interface_count
  character(:), allocatable, target :: interfaces3d(:)
  character(:), allocatable :: domain3d
  character(len=1024) :: numberSuffix
  integer(c_int) :: MUI_COMM_WORLD
  real(mytype) :: point_x,point_y,point_z
  integer(kind=8) :: i,j,k, itime,ntime
  real(mytype) :: T,dt,zlz,yly,dz
  real(mytype), dimension(ny,nz) :: bxx1,bxy1,bxz1
  allocate(yp(ny),zp(nz))

  zlz = 5.0_mytype
  yly = 20.0_mytype
  ntime = 5    ! numner of time iteration 
  dt =0.01_mytype
  T= 0.0_mytype

  interface_count =1
  dz=zlz/real(nz,mytype)
call read_yp(yp,ny)
do k=1,nz
  zp(k)=dz*(k-1)
enddo
! write(*,*) dz
  !Call mui_mpi_split_by_app_f() function to init MPI
  call mui_mpi_split_by_app_f(MUI_COMM_WORLD)

  !Allociate memory based on number of interfaces
  allocate(character(len_trim(interfaceName)+5) :: interfaces3d(interface_count))
  !For multi-domain function, "uniface_pointers_3d" should be used to collect the array of
  ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
  allocate(uniface_pointers_3d(interface_count))

  ! Obtain the domain name
  domain3d = trim(domainName)

  ! Create interface names
  do i = 1, interface_count
    !Generate character type of number suffix
    if (i < 10) then
        write (numberSuffix, "(I1)") i
    else if ((i < 100) .and. (i > 9)) then
        write (numberSuffix, "(I2)") i
    else if ((i < 1000) .and. (i > 99)) then
        write (numberSuffix, "(I3)") i
    else
        write (numberSuffix, "(I4)") i
    endif
    interfaces3d(i) = trim(interfaceName) // "_" // trim(numberSuffix)
  end do 

  ! !Create MUI interfaces. MUI interfaces will be collected by the "uniface_pointers_3d" after this subroutine
  call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, trim(domainName), interfaces3d, interface_count)


  call blasius(bxx1,bxy1,bxz1)

  do itime = 0, ntime
    
    do k = 1, nz
      do j=1, ny
        point_x = 0.0_mytype 
        point_y = yp(j)
        point_z = zp(k)
        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, &
        point_y, point_z, bxx1(j,k))

        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, point_x, &
        point_y, point_z, bxy1(j,k))

        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, point_x, &
        point_y, point_z, bxz1(j,k))
      end do
    end do

    call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, T)
    T=T+dt
    
  enddo



end program main


