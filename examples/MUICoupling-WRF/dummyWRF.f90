program read_wrf_nc
    use netcdf
    use iso_c_binding
    use, intrinsic:: iso_fortran_env, only: stdout=>output_unit, stdin=>input_unit, stderr=>error_unit
    use mui_3d_f
    use mui_general_f

    implicit none
    
    ! Define variables
    integer, parameter :: mytype = KIND(0.0D0)
    integer :: ncid, varid_lat, varid_lon, varid_u, varid_v, dimid_lat, dimid_lon,dimid_eta
    integer :: varid_w, varid_ph, varid_phb
    integer :: lon_size, lat_size,eta_size, i, j,k,i1,i2,j1,j2
    real(mytype ), allocatable :: lat(:,:), lon(:,:)
    real(mytype ), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), ph(:,:,:), phb(:,:,:)
    real(mytype ), allocatable :: x(:,:),y(:,:),z(:,:,:)
    real(mytype ) :: ref_lat,ref_lon,T0,dt,TEnd,T
    character(len=256) :: dataTime, dataDir , dataDate,base_filename
    integer :: status,fileIndx,iSync,nSyncAhead
    integer :: dataH_1 , dataH_2 , DataMin_1, DataMin_2,DataMinStep = 1


    character(len=1024) :: interfaceName = 'interface'
    integer(c_int) :: interface_count
    character(:), allocatable, target :: interfaces3d(:)
    character(len=256), allocatable :: fileNames(:)
    type(c_ptr), target :: spatial_sampler=c_null_ptr
    type(c_ptr), target :: temporal_sampler=c_null_ptr
    character(len=1024) :: numberSuffix
    integer(c_int) :: MUI_COMM_WORLD
    real(c_double) :: tolerance=1e-7_c_double
    real(mytype) :: point_x,point_y,point_z,time1,time2,waitTime
    real(mytype) :: TempValue

    interface_count =1

    call mui_mpi_split_by_app_f(MUI_COMM_WORLD)
    allocate(character(len_trim(interfaceName)+5) :: interfaces3d(interface_count))
    allocate(uniface_pointers_3d(interface_count))
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
   
    call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, trim('dummyWRF'), interfaces3d, interface_count)
    call mui_create_sampler_exact_3d_f(spatial_sampler, tolerance)
    call mui_create_temporal_sampler_exact_3d_f(temporal_sampler, tolerance)
    iSync = 1
    nSyncAhead=5
    
    ! ref_lat   =  
    ! ref_lon   =  -7.37267

    ref_lat   = 57.18855 ! 33
    ref_lon   = -7.37267 ! -79
   
    ! Input file name
    dataDir = "/media/omar/MyBook/run/UKTC/AskerveinHill/run_AskerveinHill/grid_2_fineTime/" ! "/home/omar/WORK/codes/WRF/WRF/run_case/"
    ! dataDir = "/media/omar/MyBook_6/ExCalibur/AskerveinHill/WRF_data/Grid_2_fineTime/" ! "/home/omar/WORK/codes/WRF/WRF/run_case/"
    ! dataDir = "/mnt/lustre/a2fs-work3/work/c01/c01/omah/RUN/UKTC/codes/WRF/WRF/run_AskerveinHill/Grid_2_fineTime/" ! "/home/omar/WORK/codes/WRF/WRF/run_case/"

    base_filename = "wrfout_d04_" !'wrfout_d02_'
    dataDate = "1983-10-03" !"2023-11-26"
    dataH_1  = 12 
    dataH_2  = 0
    DataMin_1 = 0
    DataMin_2 = 131
    DataMinStep = 1

    T0=0 ; TEnd= 3600; dT =60;
    ! do i = 1, DataMin_2-DataMin_1+1
    !     !Generate character type of number suffix
        
    !     write (dataTime, "(:I2:)") i

    !     print *, dataTime
        
    !     ! interfaces3d(i) = trim(interfaceName) // "_" // trim(numberSuffix)
    ! end do

    allocate(fileNames(DataMin_2-DataMin_1+1))
    print *, "============ Check the WRF out put files exists =================="
    j=0
    dataH_2=dataH_1
    do i = 1, DataMin_2-DataMin_1+1
        ! Generate the full filename based on the pattern
        
        
        write(dataTime, '(I2.2,A,I2.2,A)') dataH_2,"-", j,"-00"        

        fileNames(i) = trim(dataDir) // trim(base_filename) // trim(dataDate)// "_" // trim(dataTime)
        ! print *, fileNames(i)
        status = nf90_open(trim(fileNames(i)), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, "Error: Unable to open NetCDF file ", fileNames(i)
            stop
        endif
        status = nf90_close(ncid)

        j=j+1;
        if (j==60) then
            dataH_2=dataH_2+1
            j=0
        endif 
    end do
    print *, "============ All the files exits  =================="


    T =T0
    do fileIndx =1 , DataMin_2-DataMin_1+1
            
        ! fileNames(fileIndx)= "/home/omar/WORK/codes/WRF/WRF/run_case/wrfout_d02_2023-11-26_00:00:00"
  
        print*, "reading file : ", fileNames(fileIndx)
        ! Open the NetCDF file
        status = nf90_open(trim(fileNames(fileIndx)), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, "Error: Unable to open NetCDF file!"
            stop
        endif


        ! Read latitude and longitude dimensions
        status = nf90_inq_dimid(ncid, "south_north", dimid_lat)
        status = nf90_inquire_dimension(ncid, dimid_lat, len=lat_size)
        status = nf90_inq_dimid(ncid, "west_east", dimid_lon)
        status = nf90_inquire_dimension(ncid, dimid_lon, len=lon_size)
        status = nf90_inq_dimid(ncid, "bottom_top", dimid_eta)
        status = nf90_inquire_dimension(ncid, dimid_eta, len=eta_size)

        ! print *, "lat_size = ", lat_size
        ! print *, "log_size = ", lon_size
        ! print *, "eta_size = ", eta_size

        

        ! Read velocity variables (assuming 2D fields)
        status = nf90_inq_varid(ncid, "U", varid_u)
        status = nf90_inq_varid(ncid, "V", varid_v)
        status = nf90_inq_varid(ncid, "W", varid_w)
        
        allocate(u(lon_size, lat_size,eta_size),v(lon_size, lat_size,eta_size),w(lon_size, lat_size,eta_size))
        
        
        status = nf90_get_var(ncid, varid_u, u)
        status = nf90_get_var(ncid, varid_v, v)
        status = nf90_get_var(ncid, varid_w, w) 
        
        ! get the cartesian coordinates 
        
        ! Allocate memory for latitude and longitude arrays
        allocate(lat(lon_size, lat_size))
        allocate(lon(lon_size, lat_size))
        allocate(x(lon_size, lat_size),y(lon_size, lat_size),z(lon_size, lat_size,eta_size))
        allocate(ph(lon_size, lat_size,eta_size),phb(lon_size, lat_size,eta_size))

        ! Read latitude and longitude variables
        status = nf90_inq_varid(ncid, "XLAT", varid_lat)
        status = nf90_get_var(ncid, varid_lat, lat)
        status = nf90_inq_varid(ncid, "XLONG", varid_lon)
        status = nf90_get_var(ncid, varid_lon, lon)
        
        ! call ll_to_xy(grid%xlat_u(ips:ipe,jps:jpe),grid%xlong_u(ips:ipe,jps:jpe), &
        !         lat_ref  ,lon_ref  ,ips,ipe,jps,jpe,x_u(ips:ipe,jps:jpe),y_u(ips:ipe,jps:jpe)) 
        call ll_to_xy(lat,lon,ref_lat,ref_lon,1,lon_size,1,lat_size,x,y)

        ! level heights
        
        status = nf90_inq_varid(ncid, "PH", varid_ph)
        status = nf90_inq_varid(ncid, "PHB", varid_phb) 
        status = nf90_get_var(ncid, varid_ph, ph)
        status = nf90_get_var(ncid, varid_phb, phb)
        z= (PH+PHB)/9.81   

        ! print*, z(1,1,eta_size),x(1,1),y(1,1),x(lon_size,lat_size),y(lon_size,lat_size)
        ! Close the NetCDF file
        status = nf90_close(ncid)
        ! ==============================================================================
        !================================= MUI Sending =================================
        ! ==============================================================================
        i1=85
        j1=85
        i2=140
        j2=140
        write(*,*) "Dummy WRF is sending data at time ", T,x(i1,j1),x(i2,j2),y(i1,j1),y(i2,j2)
        ! write(*,*) "Dummy WRF is sending data at time ", x(i2,j2)-x(i1,j1),y(i2,j2)-y(i1,j1)
        ! write(*,*) "Dummy WRF is sending data at time ", T,lat(i1,j1),lon(i1,j1),lat(i2,j2),lon(i2,j2)
        ! write(*,*) "Dummy WRF is sending data at time ", T,x(130,130),x(140,140),y(130,130),y(140,140)
        ! write(*,*) "Dummy WRF Sent at 1 ", u(1,1,1), v(1,1,1), w(1,1,1)
        ! write(*,*) "Dummy WRF Sent at 2 ", u(1,1,2), v(1,1,2), w(1,1,2)
        ! do i=1,lon_size
        !     do j=1, lat_size
        !         do k = 1, eta_size
        do i=80,140
            do j=80, 140
                ! if (i> 90 .and. i<130 .and. j> 90 .and. j<130  ) then 

                ! else
                    call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, point_y, &
                        0.0_mytype, 0.0_mytype)
                    call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, point_x, point_y, &
                        0.0_mytype, 0.0_mytype)
                    call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, point_x, point_y, &
                        0.0_mytype, 0.0_mytype)
                    do k = 1, 20

                        

                        point_x = x(i,j)!real(i,mytype) !long_u(i,j) 
                        point_y = y(i,j)  !real(j,mytype) !lat_u(i,j)
                        point_z = z(i,j,k) !real(k,mytype)                    

                        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, point_y, &
                        point_z, u(i,j,k))
                        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, point_x, point_y, &
                        point_z, v(i,j,k))
                        call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, point_x, point_y, &
                        point_z, w(i,j,k))
                        

                    end do
                ! end if
            end do
        end do

        


        call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, real(T,c_double))  

        ! Deallocate memory
        deallocate(lat, lon,x,y,z, u, v,w,ph,phb)
        T=T+dT
        if (fileIndx > nSyncAhead) then
        if (iSync ==1 ) then 
             
             call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(1)%ptr, "tempvalue"//c_null_char, 0.0_mytype,0.0_mytype, &
            0.0_mytype, T-dt*nSyncAhead , spatial_sampler, temporal_sampler, TempValue)
            print *, "Dummy is reading is syncing at xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ", T-dt*nSyncAhead
        endif
        endif
        ! if (fileIndx > 5) then
        !     call cpu_time(time1)
        !     waitTime = 0.0
        !     do while (waitTime<=20) 
        !         call cpu_time(time2)
        !         waitTime = real(time2-time1)
        !     enddo
        ! endif

    enddo

    print *, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    print *, "xxxxxxxxxxxxxxxxxx Dummy Ended xxxxxxxxxxxxxxxxxxxxxxxx"
    print *, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    

end program read_wrf_nc

SUBROUTINE ll_to_xy(lat,lon,lat_ref,lon_ref,ips,ipe,jps,jpe,x,y)
    !!---------------------------------------------------------------------
    !!              ***  ROUTINE ll_to_xy  ***
    !!
    !! ** Purpose : Convert lat_long to X, Y coordinates
    !!---------------------------------------------------------------------
    integer, parameter :: mytype = KIND(0.0D0)
    INTEGER,              INTENT(IN   ) :: ips,ipe,jps,jpe
    REAL(mytype ), DIMENSION( ips:ipe, jps:jpe ) :: lat, lon
    REAL(mytype ), DIMENSION( ips:ipe, jps:jpe )::  x, y
    INTEGER                             :: i,j
    REAL(mytype )                 :: lon_ref,lat_ref,lat1,lon1,dlon,dlat,lat_ref_rad,lon_ref_rad
    REAL(mytype )                 :: lon2,lat2
    Real(mytype )                 :: pi=3.141592653589793238462643383279502884197


    !! Convert lat and Long to radian 
dvsdsv=3
    lat_ref_rad = lat_ref*pi/180
    lon_ref_rad = lon_ref*pi/180
    
    do i = ips,ipe
       do j = jps,jpe

          lat1 = lat(i,j)*pi/180.0;
          lat2 =  lat_ref_rad

          lon1 = lon_ref_rad;
          lon2 = lon_ref_rad;

          dlon = lon2 - lon1 ;
          dlat = lat2 - lat1;

          Y(i, j) = 2.0 * 6371000.0 * ASIN(SQRT(SIN(dlat / 2.0)**2 + COS(lat1) &
          * COS(lat2) * SIN(dlon / 2.0)**2)) * SIGN(1.0_mytype,-dlat)
         
          lat1 = lat_ref_rad;
          lon1 = lon(i,j)*pi/180.0;
          
          dlon = lon_ref_rad - lon1 ;
          dlat = lat_ref_rad - lat1;

          X(i, j) = 2.0 * 6371000.0 *ASIN(SQRT(SIN(dlat / 2.0)**2 + COS(lat1) * COS(lat_ref_rad) &
          * SIN(dlon / 2.0)**2)) *  SIGN(1.0_mytype,-dlon)
       enddo
    enddo



  
  END SUBROUTINE  ll_to_xy