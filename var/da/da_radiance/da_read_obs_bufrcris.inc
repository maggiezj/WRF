subroutine da_read_obs_bufrcris (obstype,iv,infile)	 
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_cris                  read bufr format cris data

  implicit none


  character(9)      ,  intent (in)  :: obstype
  character(100)    ,  intent (in)  :: infile
  type (iv_type)    ,intent (inout) :: iv
#ifdef BUFR  
  real(kind=8)    ::  obs_time  
  type (datalink_type), pointer  :: head, p, current, prev
  type(info_type)                :: info
  type(model_loc_type)           :: loc 
  type(model_loc_type)           :: loc_fov  

! Subroutine constants
  real(r_kind)                   :: ten      =  10.0_r_kind
  real(r_kind)                   :: r90      =  90.0_r_kind
  real(r_kind)                   :: r360     = 360.0_r_kind

! Number of channels for sensors in BUFR
  integer(i_kind),parameter      :: nchan = 399        !--- 399 subset ch out of 1305 ch for CrIS
  integer(i_kind),parameter      :: n_totchan  = 399
  integer(i_kind),parameter      :: maxinfo    =  33   
  integer(i_kind)                :: inst,platform_id,satellite_id,sensor_id
  integer(i_kind), allocatable   :: ptotal(:), nread(:)
  real(r_kind)                   :: crit, pred,
  integer(i_kind)                :: ifgat, iout, iobs
  logical                        :: outside, outside_all, iuse

! BUFR functions
  integer(i_kind)   :: ireadsb,ireadmg

! Variables for BUFR IO    
  real(r_double),dimension(5)           :: linele
  real(r_double),dimension(14)          :: allspot
  real(r_double),dimension(2,n_totchan) :: allchan
! real(r_double),dimension(3,10)        :: cscale  !cris no scale
! real(r_double),dimension(6)           :: cloud_frac  ! As cloud_frac is missing from BUFR, 'FCPH'  
  real(r_double),dimension(2)           :: cloud_properties 
  integer           :: numbufr,ibufr  
  logical           :: found, head_found 
  real(r_kind)      :: step, start,step_adjust
  character(len=8)  :: subset
  character(len=4)  :: senname
  character(len=80) :: allspotlist
  integer(i_kind)   :: jstart
  integer(i_kind)   :: iret
  character(10)     :: date

! Work variables for time
  integer(i_kind)   :: idate, im, iy, idd, ihh
  integer(i_kind)   :: idate5(6)

! Other work variables
  real(r_kind)                          :: piece
  real(r_kind)                          :: dlon_earth,dlat_earth
  real(r_kind)                          :: sat_look_angle, look_angle_est
  real(r_kind)                          :: sat_zenang
  real(r_kind)                          :: radi
  real(r_kind),dimension(10)            :: sscale
  real(r_kind),dimension(n_totchan)     :: temperature
  real(r_kind),allocatable,dimension(:) :: data_all
  real(r_kind)                          :: Effective_Temperature
  
  logical                               :: clear

! Set standard parameters
  integer(i_kind),parameter :: ichan=-999  ! fov-based surface code is not channel specific for cris
  !integer(i_kind),parameter:: sfc_channel=501 !used in thinning routine if cloud informatino is not available 
  real(r_kind),parameter    :: expansion=one         ! exansion factor for fov-based location code.                                                 ! use one for ir sensors.
  real(r_kind),parameter    :: tbmin  = 50._r_kind
  real(r_kind),parameter    :: tbmax  = 550._r_kind
! real(r_kind),parameter    :: earth_radius = 6371000._r_kind
  real(r_kind),parameter    :: rato   = 0.87997285_r_kind          

   if (trace_use) call da_trace_entry("da_read_obs_bufrcris")
!  0.0  Initialize variables
!-----------------------------------
  platform_id  = 17   ! jpss series, for noaa should be 1
  sensor_id = 27 ! cris, for 
  nreal  = maxinfo
  cris= obstype == 'cris' .or. obstype == 'cris-fsr'
  bad_line=-1
  step  = 3.322
  start = -51.675
!  step_adjust = 0.625_r_kind
  senname = 'CRIS'
  allspotlist= &
   'SIID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAZA BEARAZ SOZA SOLAZI SAID'
  num_bufr(:)=0
  numbufr=0
  allocate(nread(1:rtminit_nsensor))
  allocate(ptotal(0:num_fgat_time))
  nread(1:rtminit_nsensor) = 0
  ptotal(0:num_fgat_time) = 0
  iobs = 0                 ! for thinning, argument is inout  
  num_cris_file    = 0
  num_cris_local   = 0
  num_cris_global  = 0
  num_cris_used    = 0
  num_cris_thinned = 0  
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
            num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
 
   if (numbufr==0) numbufr=1
 
  bufrfile:  do ibufr=1,numbufr
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'
      end if
   end if 
   lnbufr = 96
   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning(__FILE__,__LINE__, &
         (/"Cannot open file "//infile/))
      if (trace_use) call da_trace_exit("da_read_obs_bufrcris")
      return
   end if

! Open BUFR table
   call openbf(lnbufr,'IN',lnbufr) 
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)
   iy=0
   im=0
   idd=0
   ihh=0
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt='(a,4i4,2x,a)') &
      'Bufr file date is ',iy,im,idd,ihh,trim(infile)

   ! Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------
! Allocate arrays to hold data
  nele=nreal+nchan
  allocate(data_all(nele))
   if ( ibufr == 1 ) then
      allocate (head)
      nullify  ( head % next )
      p => head
   end if
  

! Big loop to read data file

 message_loop: do while(ireadmg(lnbufr,subset,idate)>=0)  

     read_loop: do while (ireadsb(lnbufr)==0)
         num_cris_file = num_cris_file + 1	 

!    Read CRIS FOV information
         call ufbint(lnbufr,linele,7,1,iret,'FOVN SLNM QMRKH MJFC HMSL FORN  (CRCHN)')				
         
!    Check that the number of channels in BUFR is what we are expecting
         if (nint(linele(7)) /= n_totchan) then
            write(6,*)'READ_CRIS:  ***ERROR*** CrIS BUFR contains ',&
                nint(linele(7)),' channels, but CRTM expects ',n_totchan
            exit message_loop
         end if

!    Top level QC check
         if ( linele(3) /= zero) then
	         cycle read_loop  
	 end if
           ! problem with profile (QMRKH)
           ! May want to eventually set to 
           ! QMRHK <= 1, as data is, and I
           ! quote, 'slightly suspect'
 
         if ( bad_line == nint(linele(2))) then
!        zenith angle/scan spot mismatch, reject entire line
           cycle read_loop
         else
           bad_line = -1
         end if
         
         ifov = nint(linele(1))               ! field of view
         ifor = nint(linele(6))               ! field of regard
         sat_zenang  = allspot(10)            ! satellite zenith angle
        
!  CRIS field-of-view ranges from 1 to 9, corresponding to the 9 sensors measured
!  per field-of-regard.  The field-of-regard ranges from 1 to 30.  For reference, FOV 
!  pattern within the FOR is :
!                FOV#      7 8 9|7 8 9
!                FOV#      4 5 6|4 5 6
!                FOV#      1 2 3|1 2 3 (spacecraft velocity up the screen)
!                ----------------------
!                FOR#        x    x+1
!  FORs are scanned from left limb (FOR=1) to right limb (FOR=30)
!
!  For now, we will simply choose IFOV=5.  See Fig. 58 of CrIS SDR ATBD (Rev. D) for a picture.

!    Only use central IFOV
          if (ifov /= 5) cycle read_loop

!    Remove data on edges, not contain here like gsi

!    Check field of view (FOVN), field-of-regard (FORN), and satellite zenith angle (SAZA)
       iscn = nint(linele(2))                                  ! scan line
      
        if( ifov < 1 .or. ifov > 9  .or. &                      ! FOVN not betw. 1 & 9
            ifor < 1 .or. ifor > 30 .or. sat_zenang > 90._r_kind )then      ! FORN not betw. 1 & 30
           write(stdout,*) 'READ_CRIS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
              ' STRANGE OBS INFO(FOVN,FORN,SLNM,SAZA):', ifov, ifor, iscn, sat_zenang
           cycle read_loop
        end if
        if ( ifor <= 15 ) sat_zenang = -sat_zenang

          call ufbint(lnbufr,allspot,14,1,iret,allspotlist)
         if(iret /= 1) cycle read_loop

!    Extract satellite id.  If not the one we want, read next record
         ksatid = nint(allspot(14))
         ! SAID 224 is JPSS-NPP (jpss-0)
         ! SAID 225 is NOAA-20 (jpss-1)
         ! SAID 226 is NOAA-21 (jpss-2)
         if ( ( ksatid > 226) .or. ( ksatid < 224) ) then 
            write(unit=message(1),fmt='(A,I6)') 'Unknown platform: ', ksatid
            call da_warning(__FILE__,__LINE__,message(1:1))
         end if
         if ( ksatid == 224 ) then
             platform_id  = 17
             satellite_id = 0
         else if ( ksatid == 225 ) then
             platform_id  = 1
             satellite_id = 20
         else if ( ksatid == 226 ) then
             platform_id  = 1
             satellite_id = 21
         end if

!    Check observing position
         info%lat  =  allspot(8)  ! latitude
         info%lon  =  allspot(9)  ! longitude)
         if( abs(info%lat) > r90  .or. abs(info%lon) > r360 .or. &
         (abs(info%lat) == r90 .and. info%lon /= zero) )then
         write(unit=stdout,fmt=*) &
         'READ_CRIS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
               ' STRANGE OBS POINT (LAT,LON):', info%lat, info%lon
            cycle read_loop
         end if		 
	 		 
         call da_llxy (info, loc, outside, outside_all)	
         if (outside_all) cycle
	     inst = 0	 
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do	 
         if (inst == 0) cycle read_loop		 
		 
!    Check obs time
         idate5(1) = nint(allspot(2)) ! year
         idate5(2) = nint(allspot(3)) ! month
         idate5(3) = nint(allspot(4)) ! day
         idate5(4) = nint(allspot(5)) ! hour
         idate5(5) = nint(allspot(6)) ! minute
         idate5(6) = nint(allspot(7)) ! second		
		
         if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
             idate5(2) < 1    .or. idate5(2) >   12 .or. &
             idate5(3) < 1    .or. idate5(3) >   31 .or. &
             idate5(4) <0     .or. idate5(4) >   24 .or. &
             idate5(5) <0     .or. idate5(5) >   60 )then

            write(6,*)'READ_CRIS:  ### ERROR IN READING ', 'CRIS', ' BUFR DATA:', &
                ' STRANGE OBS TIME (YMDHM):', idate5(1:5)
             cycle read_loop
         end if
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)		
         if ( obs_time < time_slots(0) .or.  &
           obs_time >= time_slots(num_fgat_time) ) cycle read_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat) ) exit
         end do	
         num_cris_global = num_cris_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1  

!    Compare CRIS satellite scan angle and zenith angle
 
        look_angle_est = (start + float(ifor)*step)*deg2rad
        sat_look_angle=asin(rato*sin(sat_zenang*deg2rad))

        if(abs(sat_look_angle)*rad2deg > 90._r_kind) then            ! MAX_SENSOR_ZENITH_ANGLE=90°
          write(6,*)'READ_CRIS WARNING lza error ',sat_look_angle,look_angle_est
          cycle read_loop
        end if

        if (abs(sat_look_angle - look_angle_est)*rad2deg > one) then
           write(6,*)' READ_CRIS WARNING uncertainty in look angle ', &
               look_angle_est*rad2deg,sat_look_angle*rad2deg,sat_zenang,sis,ifor,start,step,allspot(11),allspot(12),allspot(13)
           bad_line = iscn
           cycle read_loop
        endif


!   Clear Amount  (percent clear)


!xxx        call ufbrep(lnbufr,cloud_frac,1,6,iret,'TOCC')
!xxx!    Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
!xxx        pred = cloud_frac(1)
            pred = 100.0_r_kind    ! pred needs to have a value for WCOSS debug execution to work.

! As cloud_frac is missing from BUFR, use proxy of warmest fov over 
! non ice surfaces.  Fixed channels (assuming the 399 set) for now.
! This is moved to below where the radiances are read in.

        ! call ufbrep(lnbufr,cloud_frac,1,6,iret,'FCPH')
        ! if (outside) cycle ! No good for this PE		
        ! num_iasi_local = num_iasi_local + 1
		 
        ! write(unit=info%date_char, &
        ! fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
        ! idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
        ! ':', idate5(5), ':', idate5(6)
        ! info%elv = 0.0  !aquaspot%selv		   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!new algorithm??????? 		!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!          If cloud_properties is missing from BUFR, use proxy of warmest fov. 
!          the surface channel is fixed and set earlier in the code (501).from gsi
!          this is moved to below where the radiances are read in.


 
         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad
            crit = 1. 
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_cris_thinned=num_cris_thinned+1
               cycle
            end if
         end if	
	
! CrIS data unscaled, section removed from read_iasi

!    Read Cris channel number(CHNM) and radiance (SCRA)

         call ufbint(lnbufr,allchan,2,n_totchan,iret,'SCRA CHNM')
         if( iret /= n_totchan)then
              write(unit=stdout,fmt=*) &
			  'READ_CRIS:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
                iret, ' CH DATA IS READ INSTEAD OF ',n_totchan				
           cycle read_loop
         end if

         num_cris_used = num_cris_used + 1		 
         nread(inst) = nread(inst) + 1
								
         iskip = 0
         jstart=1
         do i=1,n_totchan
!     check that channel number is within reason
            if (( allchan(1,i) > zero .and. allchan(1,i) < 99999._r_kind) .and. &  ! radiance bounds
               (allchan(2,i) < 1305._r_kind .and. allchan(2,i) > zero )) then     ! chan # bounds
!         radiance to BT calculation
              radi = allchan(1,i)*1000.0_r_kind    !Concersion from W to mW

              if (rtm_option == rtm_option_crtm) then
#ifdef CRTM
                 call CRTM_Planck_Temperature(inst,i,radi,temperature(i))
#endif
              else if (rtm_option == rtm_option_rttov) then
#ifdef RTTOV
                 Effective_Temperature = LOG( ( coefs(inst)%coef%planck1(i) / radi ) + 1 )
                 Effective_Temperature = coefs(inst)%coef%planck2(i)/Effective_Temperature
                 temperature(i) = ( Effective_Temperature - coefs(inst)%coef%ff_bco(i) ) / &
                                  coefs(inst)%coef%ff_bcs(i)
#endif
              end if
              if(temperature(i) < tbmin .or. temperature(i) > tbmax ) then
                 temperature(i) = min(tbmax,max(zero,temperature(i)))

              end if
            else           ! error with channel number or radiance
!             write(6,*)'READ_IASI:  iasi chan error',i,allchan(1,i), allchan(2,i)
              temperature(i) = min(tbmax,max(zero,temperature(i)))

            end if
         end do
         if(iskip > 0) write(6,*) ' READ_CRIS : iskip > 0 ',iskip  !!!!!! iskip not changed?????????????????????
!         if( iskip >= 10 )cycle read_loop 		


! (Comment copied from above:)
!          If cloud_properties is missing from BUFR, use proxy of warmest fov.
!          the surface channel is fixed and set earlier in the code (501).from gsi
!          this is moved to below where the radiances are read in.
           
!          Cloud / clear tests.
           clear = .false.
           pred = zero

!          Cloud information  may be missing depending on how the VIIRS granules align
!          with the CrIS granules.  
!          Cloud Amount, TOCC is total cloud cover [%], HOCT is cloud height [m]
           call ufbint(lnbufr,cloud_properties,2,1,iret,'TOCC HOCT')
           if ( cloud_properties(1) <= r100 .and. cloud_properties(1) >= zero .and. &
                cloud_properties(2) < 1.0e6_r_kind ) then
!             Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
              if ( cloud_properties(1) < one ) then     !Assume clear
                 clear = .true.
              else                                ! Assume a lapse rate to convert hgt to delta TB.
                 pred = cloud_properties(2) *7.0_r_kind / r1000
              endif
           else
         
!          If cloud_properties is missing from BUFR, use proxy of warmest fov. 
!          the surface channel is fixed and set earlier in the code (501).
              if ( tsavg*0.98_r_kind <= temperature(i)) then   ! 0.98 is a crude estimate of the surface emissivity
                   clear = .true.
                else
                   pred = (tsavg * 0.98_r_kind - temperature(i)) 
                endif
           endif   !clearest FOV check
           


         data_all(5) = sat_zenang             ! satellite zenith angle (deg)
         data_all(6) = allspot(11)            ! satellite azimuth angle (deg)
         data_all(9) = allspot(12)            ! solar zenith angle (deg)
         data_all(10)= allspot(13)            ! solar azimuth angle (deg)
         do l=1,nchan
            data_all(l+nreal) = temperature(l)   ! brightness temerature
         end do
				
!  4.0   assign information to sequential radiance structure
!--------------------------------------------------------------------------
         allocate ( p % tb_inv (1:nchan) )
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = 1
         p%scanpos          = ifovn
         p%satzen           = data_all(5)
         p%satazi           = data_all(6)  ! look angle (deg) ! airsspot%bearaz
         p%solzen           = data_all(9)
         p%solazi           = data_all(10)
         p%tb_inv(1:nchan)  = data_all(nreal+1:nreal+nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat		

         allocate (p%next)   ! add next data
         p => p%next
         nullify (p%next) 		
     end do read_loop
  end do
  call closbf(lnbufr)

  !Deallocate temporary array for next bufrfile do loop
  deallocate(data_all)
end do bufrfile

   if (thinning .and. num_cris_global > 0 ) then

#ifdef DM_PARALLEL 
      
      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do 
      end do
   
      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i)
            end do
         end do 
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i) = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )

#endif

      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_cris_used_tmp = num_cris_used
      do j = 1, num_cris_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_cris_thinned = num_cris_thinned + 1
            num_cris_used = num_cris_used - 1
            nread(n) = nread(n) - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_cris_used
   iv%total_rad_channel = iv%total_rad_channel + num_cris_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_cris_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_cris_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=message(1),fmt='(a)') '   num_cris_file num_cris_global  num_cris_local   num_cris_used num_cris_thinned'
   write(unit=message(2),fmt='(5(6x,i10))') num_cris_file, num_cris_global, num_cris_local, num_cris_used, num_cris_thinned
   call da_message(message(1:2))

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
 

   do i = 1, iv%num_inst 
    
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad  = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad
      call da_allocate_rad_iv (i, nchan, iv)
	  
   end do

   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0
   p => head

   do n = 1, num_cris_used
      i = p%sensor_index 
      nread(i) = nread(i) + 1 
  
      call da_initialize_rad_iv (i, nread(i), iv, p)
  
      current => p
      p => p%next
      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
   deallocate ( p )
   deallocate (nread)
   deallocate (ptotal)

   call closbf(lnbufr)
   close(lnbufr)

!   call da_free_unit(lnbufr)

   if (trace_use) call da_trace_exit("da_read_obs_bufrcris")
#else
   call da_error(__FILE__,__LINE__,(/"Needs to be compiled with a BUFR library"/))
#endif  


end subroutine da_read_obs_bufrcris
