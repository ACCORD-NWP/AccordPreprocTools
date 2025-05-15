PROGRAM SELECT_GPSSOL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Selection package for ground-based GPS ZTD observations
  ! This package:
  ! * reads in all ZTD observations and f.g. departures given by a
  !   file
  ! * sorts and calculates statistics on these
  !     () filters out all the data for which
  !        statid,date,hhmmss,lat,lon,alt,obsvalue,dt is missing
  ! * selects stations which pass a list of criteria
  !    - distribution of f.g. departures must be Gaussian
  !    - most stable center if several centers for the same station
  !    - remove stations which are too close
  ! * writes out a list of stations containing:
  !     () statid,lat,lon,alt,dt,bias,obserr
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! History:
  !  Mar2006  Paul Poli, METEO FRANCE    Original code
  !  Oct2007  Paul Poli, METEO FRANCE    Added possibility to retain
  !                                      only obs closest to central time
  !  Apr2008 Paul Poli, METEO FRANCE     rounding problem in diff method correction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINES CALLED: 5
  !  * read_odb_gpssol
  !  * calc_statistics_gpssol
  !  * print_statistics
  !  * calc_statistics_active_gpssol
  !  * print_statistics_active
  !  * select_gpssol_stations
  !  * write_list_gpssol
  ! EXTERNAL MODULES: parkind1, gpssol_mod
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod

  implicit none

  ! ODB data
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max) :: stat_cen
  TYPE(STAT_CEN_ACTIVE_TYPE), dimension(1:nstat_cen_active_max) :: stat_cen_active
  TYPE(STATISTICS_ACTIVE_TYPE) :: statistics_active

  ! error flags
  integer(KIND=JPIM) :: error_flag, AllocStatus
  character(len=100) :: subrout

  print *,'start'
  nstat_cen=0
  nstat_cen_active=0

  call READ_ODB_GPSSOL(trim(fic_list_files)//trim(ficoutput_ext), &
    stat_cen,stat_cen_active,error_flag)
  subrout='read_odb_gpssol'
  IF (error_flag .ne. no_error) GOTO 900

  call CALC_STATISTICS_GPSSOL(stat_cen,error_flag)
  subrout='calc_statistics'
  IF (error_flag .ne. no_error) GOTO 900

  call PRINT_STATISTICS(stat_cen,'stats.txt'//trim(ficoutput_ext))

  call CALC_STATISTICS_ACTIVE_GPSSOL(stat_cen_active, &
       statistics_active, error_flag)
  subrout='calc_statistics_active'
  IF (error_flag .ne. no_error) GOTO 900

  call PRINT_STATISTICS_ACTIVE(stat_cen_active,statistics_active, &
       'stats_active.txt'//trim(ficoutput_ext))

  call SELECT_GPSSOL_STATIONS(stat_cen,error_flag)
  subrout='select_gpssol_stations'
  IF (error_flag .ne. no_error) GOTO 900

  call WRITE_LIST_GPSSOL(trim(fic_list_gpssol)//trim(ficoutput_ext), &
       stat_cen,error_flag)
  subrout='write_list_gpssol'
  IF (error_flag .ne. no_error) GOTO 900

      print '(A)','select_gpssol: program terminated OK'
      GOTO 901
      ! return
      
  900 print '(A,I10,A,A)','Error ',error_flag,' in ',trim(subrout)
      ! return
  901 CONTINUE
  
END PROGRAM SELECT_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_ODB_GPSSOL(list_files, &
  stat_cen,stat_cen_active,error_flag)
  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod
  implicit none

  ! I/O variables
  character(len=*),   intent(in) :: list_files
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max), &
                      intent(inout) :: stat_cen
  TYPE(STAT_CEN_ACTIVE_TYPE), dimension(1:nstat_cen_active_max), &
                      intent(inout) :: stat_cen_active
  integer(KIND=JPIM), intent(out) :: error_flag

  ! local variables
  integer(KIND=JPIM) :: iostatus, AllocStatus
  integer(KIND=JPIM) :: ific, nfics, nlines
  character(len=100) :: string100
  character(len=300) :: string300
  integer(KIND=JPIM), parameter :: nfics_max = nperiod_max
  character(len=100), dimension(1:nfics_max) :: fics
  logical            :: didnotfindnull
  
  character(len=8)   :: statid1
  character(len=6)   :: str_time1
  integer(KIND=JPIM) :: date1, time1, dt1, tslot1, abs_time_diff1
  real(KIND=JPRB)    :: lat1, lon1, alt1, diff_alt1,alt_mod1, &
                        obsvalue1, fg_depar1, ZHD_calc1, ZWD_calc1, an_depar1
  integer(KIND=JPIM) :: status_active1, status_blacklisted1

  integer(KIND=JPIM) :: istat_cen, idata, jdata, idt, &
                        tperiod1, idata_tperiod, iperiod, i_id(0:15), imins
  logical            :: ll_found, exclude1, do_consider_obs, do_replace_obs

  error_flag = no_error

  ! read list of files
  print '(A)','read_odb_gpssol: Opening file '//trim(list_files)
  OPEN(unit=55,file=trim(list_files),status='OLD',action='READ', &
       iostat=iostatus)
  nfics=0
  IF (iostatus .ne. no_error) THEN
    print '(A)','read_odb_gpssol: Could not open file '//trim(list_files)
    error_flag = error_reading_file
    return
  ENDIF

  readloop: DO WHILE (iostatus == 0 .and. nfics .lt. nfics_max)
    READ(unit=55,fmt='(A100)',iostat=iostatus) string100
    IF (iostatus /= 0) EXIT
    nfics=nfics+1
    fics(nfics)=trim(string100)
  ENDDO readloop

  readif: IF (iostatus > 0) THEN
            print '(A)','read_odb_gpssol: An error has occurred while reading.'
            error_flag = error_reading_file
            CLOSE(55)
            return
          ELSE
            print '(A,I10,A)','read_odb_gpssol: Finished reading. Found',nfics,' records.'
            CLOSE(55)
          ENDIF readif

  ! read file by file
!  open(54,file='BRUS-NKG.dat')
  nperiod=nfics
  nperiod7=nperiod*nslot_max
  fileloop: DO ific=1,nfics

    tperiod1 = ific
    print '(A)','read_odb_gpssol: Opening file '
    print '(A)',trim(fics(ific))
    OPEN(unit=55,file=trim(fics(ific)),status='OLD',action='READ', &
       iostat=iostatus)
    IF (iostatus .ne. no_error) THEN
      print '(A)','read_odb_gpssol: Could not open file'
      error_flag = error_reading_file
      return
    ENDIF

    DO nlines=1,6
      READ(unit=55,fmt='(A300)',iostat=iostatus) string300
      IF (iostatus .ne. no_error) THEN
        error_flag = error_reading_file
        return
      ENDIF
    ENDDO
    nlines=0

    readloop_1: DO WHILE (iostatus == 0)
      READ(unit=55,fmt='(A300)',iostat=iostatus) string300
    IF (iostatus /= 0) EXIT
    !IF (string300(1:14) .ne. ': nb de lignes') THEN
    !    READ(string300,*,iostat=iostatus) &
    !    statid1, date1, time1, dt1, tslot1, &
    !    lat1, lon1, alt1, alt_mod1, &
    !    obsvalue1, fg_depar1, ZHD_calc1, ZWD_calc1, &
    !    status_active1, status_blacklisted1, an_depar1

    IF (string300(1:14) .ne. ': nb de lignes') THEN
        READ(string300,*,iostat=iostatus) &
        statid1, date1, time1, dt1, tslot1,&
        lat1, lon1, alt1, alt_mod1, &
        obsvalue1, fg_depar1, &
        status_active1, status_blacklisted1, an_depar1
       ! tslot1=1

         IF (iostatus .eq. no_error) THEN
               ! check there is no missing value
             IF (obsvalue1 >= 1.8 .and. obsvalue1 <= 3.) THEN
               ! retain line
                 nlines=nlines+1
                 status_blacklisted1=1 !! F. Meier 20240301 
                ! determine if non-active station (status_blacklisted1=1)
                ! selection only perform on non-active data(monitoring)
    ifpassive: IF (status_blacklisted1 .eq. 1) THEN 

                   ! if data are 59 past the hour, separate them from the rest of the
                   ! series
                   IF (flag_separate_59mins) THEN
                      write (str_time1,'(I6.6)') time1
                      IF (str_time1(3:4) .eq. '59') THEN
                      statid1(8:8)='9'
                      dt1=59
                      ENDIF
                   ENDIF

                  ! if data are for certain minutes, exclude them
                  IF (flag_retain_only_select_mins) THEN
                  write (str_time1,'(I6.6)') time1
                  exclude1=.false.
                     DO imins=1,n_retains
                        IF (statid1(5:8) .eq. cen_mins(imins).and. &
                        & str_time1(3:4) .ne. retain_mins(imins)) exclude1=.true.
                     ENDDO
                     IF (exclude1) cycle readloop_1
                  ENDIF
            
                  ! if data are for 59 mins, may exclude them
                  IF (flag_exclude_59mins .and. str_time1(3:4) .eq. '59') THEN
                  cycle readloop_1
                  ENDIF

                  ! determine if new station
                  istat_cen=FIND_STATID(nstat_cen,nstat_cen_max,stat_cen(1:nstat_cen_max),statid1)
                  IF (istat_cen.lt.0) THEN
                     ! new stat_cen: add it to the list
                     IF (nstat_cen .eq. nstat_cen_max) THEN
                     print '(A,I10)','read_odb_gpssol: too many stations ',nstat_cen_max
                     error_flag = error_programming
                     return
                  ENDIF
               ! assign the read values into the correct index of array stat_cen
               nstat_cen=nstat_cen+1
               stat_cen(nstat_cen)%namestr=statid1
               stat_cen(nstat_cen)%lat=lat1
               stat_cen(nstat_cen)%lon=lon1
               stat_cen(nstat_cen)%alt=alt1
               stat_cen(nstat_cen)%alt_mod=alt_mod1 !alt1 !alt_mod1
               stat_cen(nstat_cen)%cst_lat=.TRUE.
               stat_cen(nstat_cen)%cst_lon=.TRUE.
               stat_cen(nstat_cen)%cst_alt=.TRUE.
               stat_cen(nstat_cen)%cst_diff_alt=.TRUE.
               stat_cen(nstat_cen)%ndt=0
               stat_cen(nstat_cen)%ndts(:)=0
               stat_cen(nstat_cen)%tperiod_statistics(:)%ndata=0

               ! each period has it's own statistics
               DO iperiod=1,nperiod
                  stat_cen(nstat_cen)%tperiod_statistics(iperiod)%mask(:)=-1
                  allocate(stat_cen(nstat_cen)%tperiod_statistics(iperiod)% &
                         & data_present7(nslot_max))
                  stat_cen(nstat_cen)%tperiod_statistics(iperiod)% &
                   & data_present7(:)=0
               ENDDO

               stat_cen(nstat_cen)%statistics%ndata=0
               stat_cen(nstat_cen)%statistics%ndata_10days=0
               allocate(stat_cen(nstat_cen)%statistics%data_present(nperiod))
               stat_cen(nstat_cen)%statistics%data_present(:)=0
               allocate(stat_cen(nstat_cen)%statistics%data_present7(nperiod7))
               stat_cen(nstat_cen)%statistics%data_present7(:)=0
               istat_cen=nstat_cen
            ELSE
            ! else for old stat_cen: check station position is constant
               IF (abs(stat_cen(istat_cen)%lat-lat1) .gt. 0.1) &
                  stat_cen(istat_cen)%cst_lat=.FALSE.
               IF (abs(stat_cen(istat_cen)%lon-lon1) .gt. 0.1) &
                  stat_cen(istat_cen)%cst_lon=.FALSE.
               IF (abs(stat_cen(istat_cen)%alt-alt1) .gt. 0.1) &
                  stat_cen(istat_cen)%cst_alt=.FALSE.
               IF (abs(abs(stat_cen(istat_cen)%alt_mod-stat_cen(istat_cen)%alt) &
                  & -abs(alt_mod1-alt1)) .gt. 0.1) &
                  stat_cen(istat_cen)%cst_diff_alt=.FALSE.
            ENDIF

            ! if flag_retain_only_central_time_obs, determine whether new obs 
            ! needs to be considered (either as addition or replacement), or discarded
            ! central time
            call GET_ABS_TIMEDIFF(time1,abs_time_diff1)
            IF (flag_retain_only_central_time_obs .and. &
              & stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata > 0) THEN
               IF (stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata > 1) THEN
                  print *,'ERROR!!!!! flag_retain_only_central_time_obs is active &
                         & but there is more than one observation in tperiod_statistics', tperiod1
                  STOP
               ENDIF
               IF (stat_cen(istat_cen)%tperiod_statistics(tperiod1)%abs_time_diff > &
                 & abs_time_diff1) THEN
                  do_consider_obs=.true.
                  do_replace_obs=.true.
               ELSE
                  do_consider_obs=.false.
                  do_replace_obs=.false.
               ENDIF
            ELSE
               do_consider_obs=.true.
               do_replace_obs=.false.
            ENDIF
            
            ! add data to the data series
            IF (do_consider_obs) THEN ! [
               IF (stat_cen(istat_cen)%statistics%ndata .eq. ndata_max) THEN
                  error_flag=error_programming
                  print '(A)','read_odb_gpssol: too many data'
                  print '(A,A)','station :',statid1
                  return
               ENDIF

               IF (.not.do_replace_obs) stat_cen(istat_cen)%statistics%ndata = &
                                      & stat_cen(istat_cen)%statistics%ndata + 1_JPIM

               IF (dt1 .eq. 0) THEN
                  SELECT CASE (statid1(5:7))
                  CASE ('ACR')
                   dt1=15.
                  CASE ('ASI')
                   dt1=15.
                  CASE ('BKG')
                   dt1=60.
                  CASE ('BME')
                   dt1=60.
                  CASE ('GFZ')
                   dt1=15.
                  CASE ('GF1')
                   dt1=15.
                  CASE ('GOP')
                   dt1=60.
                  CASE ('IEE')
                   dt1=10.
                  CASE ('IGE')
                   dt1=15.
                  CASE ('KNM')
                   dt1=15.
                  CASE ('LPT')
                   dt1=5.
                  CASE ('MET')
                   dt1=15.
                  CASE ('MTG')
                   dt1=15.
                  CASE ('MTR')
                   dt1=15.
                  CASE ('NKG')
                   dt1=5.
                  CASE ('NGA')
                   dt1=15.
                  CASE ('ROB')
                   dt1=15.
                  CASE ('SGN')
                   dt1=15.
                  CASE ('TUA')
                   dt1=60.
                  CASE DEFAULT
                   print '(A,A)','center not known !!!!!!!!!!!!!!! ', &
                       & stat_cen(istat_cen)%namestr
                   stop
                  END SELECT
               ENDIF

               idata=stat_cen(istat_cen)%statistics%ndata
!               stat_cen(istat_cen)%data(idata)%tperiod =tperiod1
!               stat_cen(istat_cen)%data(idata)%tslot   =tslot1
!               stat_cen(istat_cen)%data(idata)%dt      =dt1
               stat_cen(istat_cen)%data(idata)%fg_depar=fg_depar1
!               stat_cen(istat_cen)%data(idata)%obsvalue=obsvalue1
!               stat_cen(istat_cen)%data(idata)%ZHD_calc=ZHD_calc1
!               stat_cen(istat_cen)%data(idata)%ZWD_calc=ZWD_calc1
!               stat_cen(istat_cen)%statistics%mask(idata) = idata

               IF (stat_cen(istat_cen)%ndt.eq.0) THEN
                  ll_found=.FALSE.
               ELSE
                  idt=0
                  ll_found=.FALSE.
                  DO WHILE (idt .lt. stat_cen(istat_cen)%ndt .and. .not.(ll_found))
                     idt=idt+1
                     IF (stat_cen(istat_cen)%dts(idt) .eq. INT(dt1,JPIM)) &
                        ll_found=.TRUE.
                  ENDDO
               ENDIF
               IF (.not.(ll_found)) THEN
                  IF (.not.do_replace_obs) stat_cen(istat_cen)%ndt=stat_cen(istat_cen)%ndt+1
                  idt=stat_cen(istat_cen)%ndt
                  stat_cen(istat_cen)%dts(idt)=INT(dt1,JPIM)
               ENDIF
               IF (.not.do_replace_obs) stat_cen(istat_cen)%ndts(idt)=stat_cen(istat_cen)%ndts(idt)+1

               IF (stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata .eq. &
                 & ndata_period_max) THEN
                  error_flag=error_programming
                  print '(A)','read_odb_gpssol: too many data per time period'
                  print '(A,A)','station :',statid1
                  print '(A,I1)','time period :',tperiod1
                  return
               ENDIF
               IF (.not.do_replace_obs) &
                     stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata = &
                   & stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata +1
               idata_tperiod = stat_cen(istat_cen)%tperiod_statistics(tperiod1)%ndata
               stat_cen(istat_cen)%tperiod_statistics(tperiod1)%mask(idata_tperiod) = &
                & idata

               stat_cen(istat_cen)%statistics%data_present(tperiod1)=1
               stat_cen(istat_cen)%statistics%data_present7 &
                & ((tperiod1-1)*nslot_max+tslot1)=1
               stat_cen(istat_cen)%tperiod_statistics(tperiod1)%data_present7 &
                & (tslot1)=1
                
               stat_cen(istat_cen)%tperiod_statistics(tperiod1)%abs_time_diff=abs_time_diff1
!               IF (statid1 .eq. 'BRUS-NKG') write (54,*) 'BRUS-NKG',date1,time1,tperiod1,tslot1,fg_depar1
            ENDIF ! ] do_consider_obs

            ! ] active data
            ELSEIF (status_blacklisted1 .eq. 0) THEN
            ! [
            
            ! determine if new station
            i_id=FIND_CHAIN(nstat_cen_active,8,nstat_cen_active_max, &
                      stat_cen_active%namestr,statid1)
            istat_cen=i_id(1)
            IF (istat_cen.lt.0) THEN
            ! new stat_cen_active: add it to the list
               IF (nstat_cen_active .eq. nstat_cen_active_max) THEN
                  print '(A,I10)','read_odb_gpssol: too many active stations ',nstat_cen_active_max
                  error_flag = error_programming
                  return
               ENDIF
               nstat_cen_active=nstat_cen_active+1
               stat_cen_active(nstat_cen_active)%namestr=statid1
               print '(A,A,I10)','Added stat_cen_active: ',statid1,nstat_cen_active
               stat_cen_active(nstat_cen_active)%statistics%ndata(:)=0
               stat_cen_active(nstat_cen_active)%statistics%ndata_rejected(:)=0
               istat_cen=nstat_cen_active
            !ELSE
            ! old stat_cen_active: check station position is constant
            ENDIF

            ! add data to the data series
            IF (status_active1 .eq. 1) THEN ! data active
               IF (stat_cen_active(istat_cen)%statistics%ndata(0) .eq. ndata_max) THEN
                  error_flag=error_programming
                  print '(A)','read_odb_gpssol: too many active data'
                  print '(A,A)','station :',statid1
                  return
               ENDIF
               stat_cen_active(istat_cen)%statistics%ndata(0) = &
                & stat_cen_active(istat_cen)%statistics%ndata(0) + 1_JPIM

               idata=stat_cen_active(istat_cen)%statistics%ndata(0)
   !            stat_cen_active(istat_cen)%data(idata)%tperiod =tperiod1
               stat_cen_active(istat_cen)%data(idata)%tslot   =tslot1
               stat_cen_active(istat_cen)%statistics%ndata(tslot1) = &
                & stat_cen_active(istat_cen)%statistics%ndata(tslot1) + 1_JPIM
   !            stat_cen_active(istat_cen)%data(idata)%dt      =dt1
               stat_cen_active(istat_cen)%data(idata)%fg_depar=fg_depar1
   !            stat_cen_active(istat_cen)%data(idata)%obsvalue=obsvalue1
               stat_cen_active(istat_cen)%data(idata)%an_depar=an_depar1
   !            stat_cen_active(istat_cen)%data(idata)%ZHD_calc=ZHD_calc1
   !            stat_cen_active(istat_cen)%data(idata)%ZWD_calc=ZWD_calc1
   !            stat_cen_active(istat_cen)%statistics%mask(idata) = idata

   !            IF (statid1 .eq. 'BRUS-NKG') write (54,*) 'BRUS-NKG active',tperiod1,tslot1,fg_depar1

            ELSE ! data rejected
               stat_cen_active(istat_cen)%statistics%ndata_rejected(0) = &
                & stat_cen_active(istat_cen)%statistics%ndata_rejected(0) + 1_JPIM
               stat_cen_active(istat_cen)%statistics%ndata_rejected(tslot1) = &
                & stat_cen_active(istat_cen)%statistics%ndata_rejected(tslot1) + 1_JPIM
            ENDIF

            ENDIF ifpassive !]

          ENDIF
        ENDIF
        iostatus = 0
      ENDIF
!      fics(nfics)=trim(string100)
    ENDDO readloop_1
  
    readif_1: IF (iostatus > 0) THEN
                print '(A)','read_odb_gpssol: An error has occurred while reading.'
                error_flag = error_reading_file
                CLOSE(55)
                return
              ELSE
                print '(A,I10,A)','read_odb_gpssol: Finished reading. Found',nlines,' records.'
                CLOSE(55)
              ENDIF readif_1
  ENDDO fileloop
!  close(54)

END SUBROUTINE READ_ODB_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION DIDNOTFINDNULL(string) result(flag_ok)
  USE parkind1, only: JPIM
  implicit none

  ! I/O variables
  character(len=*) :: string
  logical          :: flag_ok

  ! local variable
  integer(KIND=JPIM) :: ichar, nchar
  
  flag_ok = .TRUE.
  nchar = len_trim(string)
  IF (nchar .lt. 4) return
  DO ichar=1,len_trim(string)-4
    IF (string(ichar:ichar+4) .eq. 'NULL') THEN
      flag_ok = .FALSE.
      return
    ENDIF
  ENDDO
  
END FUNCTION DIDNOTFINDNULL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CALC_STATISTICS_GPSSOL(stat_cen,error_flag)
  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod
  implicit none

  ! I/O variables
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max), intent(inout) :: stat_cen
  integer(KIND=JPIM), intent(out) :: error_flag

  ! local variables
  real(KIND=JPRB) :: sum1, sum2, mean, stdev, &
                     mean_min, mean_max, stdev_min, stdev_max
  integer(KIND=JPIM) :: istat_cen,ndata,iperiod,idata,iidata, &
                        nperiod_ok,nperiod_10days_ok,nperiod7_10days_ok
  real(KIND=JPRB), dimension(1:ndata_max) :: tab_to_test
  real(KIND=JPRB), dimension(1:nperiod_max) :: tab_mean, tab_stdev
  real(KIND=JPRB), dimension(1:2) :: tab2_mean, tab2_stdev, tab2_10days
  logical :: IS_DISTRIBUTION_GAUSSIAN
  logical, parameter :: flag_ouputdatafiles = .FALSE.
  character(len=12) :: outputdatafiles_filename

  nperiod_10days_ok=minval((/nperiod_10days,nperiod/))
  nperiod7_10days_ok=minval((/nperiod_10days*nslot_max,nperiod7/))

  write(*,*) 'nperiod_10days_ok=',nperiod_10days_ok,'nperiod=',nperiod,'nperiod_10days=',nperiod_10days
  write(*,*) 'nperiod7_10days_ok=',nperiod7_10days_ok,nperiod_10days*nslot_max,'nperiod7=',nperiod7
  write(*,*) 'nstat_cen=',nstat_cen 
  DO istat_cen=1,nstat_cen

     ! statistics for the whole time period
     ndata=stat_cen(istat_cen)%statistics%ndata
     tab_to_test(:)=0.
     IF (ndata .gt. 2) THEN
        sum1=0.
        DO idata=1,ndata
          sum1=sum1+stat_cen(istat_cen)%data(idata)%fg_depar
          tab_to_test(idata)=stat_cen(istat_cen)%data(idata)%fg_depar
        ENDDO
        mean=sum1/REAL(ndata,JPRB)
        sum2=0.
        DO idata=1,ndata
          sum2=sum2+(stat_cen(istat_cen)%data(idata)%fg_depar-mean)**2.
        ENDDO
        stdev=sqrt(sum2/REAL(ndata-1,JPRB))
        if (istat_cen.eq.1) then
         write(*,*) 'haliho=',mean,stdev
        endif
        stat_cen(istat_cen)%statistics%gauss= &
         & IS_DISTRIBUTION_GAUSSIAN(ndata,tab_to_test(1:ndata),mean,stdev)
        IF (flag_ouputdatafiles) THEN
           outputdatafiles_filename=stat_cen(istat_cen)%namestr//'.dat'
           OPEN(UNIT=55,file=outputdatafiles_filename,action='WRITE')
           WRITE(unit=55,fmt='(I8)') ndata
           WRITE(unit=55,fmt='(4G20.12)') tab_to_test(1:ndata)
           CLOSE(55)
           IF (stat_cen(istat_cen)%namestr .eq. 'ALDE-MET') stop
           IF (stat_cen(istat_cen)%namestr .eq. 'ABBS-MET') THEN
                print *,'ABBS-MET stat_cen(istat_cen)%statistics%data_present(:)', &
                stat_cen(istat_cen)%statistics%data_present(:)
                print *,'ABBS-MET SUM ',SUM(stat_cen(istat_cen)%statistics%data_present(:))
                stop
           ENDIF
        ENDIF
     ELSE ! not enough data
        mean=0.
        stdev=0.
        stat_cen(istat_cen)%statistics%gauss=.FALSE.
     ENDIF
     stat_cen(istat_cen)%statistics%mean = mean
     stat_cen(istat_cen)%statistics%stdev = stdev
     stat_cen(istat_cen)%statistics%time_coverage_10days = &
     & REAL(SUM(stat_cen(istat_cen)%statistics%data_present &
     & (1:nperiod_10days_ok)),JPRB)/ &
     & REAL(nperiod_10days_ok)*100.
     stat_cen(istat_cen)%statistics%time_coverage = &
     & REAL(SUM(stat_cen(istat_cen)%statistics%data_present(:)),JPRB)/ &
     & REAL(nperiod)*100.
     stat_cen(istat_cen)%statistics%time_coverage7_10days = &
     & REAL(SUM(stat_cen(istat_cen)%statistics%data_present7 &
     & (1:nperiod7_10days_ok)),JPRB)/ &
     & REAL(nperiod7_10days_ok)*100.
     stat_cen(istat_cen)%statistics%time_coverage7 = &
     & REAL(SUM(stat_cen(istat_cen)%statistics%data_present7(:)),JPRB)/ &
     & REAL(nperiod7)*100.

     ! statistics by time period
     nperiod_ok=0
     DO iperiod=1,nperiod

        ndata=stat_cen(istat_cen)%tperiod_statistics(iperiod)%ndata
        IF (ndata .gt. 0) THEN
           sum1=0.
           DO idata=1,ndata
              iidata=stat_cen(istat_cen)%tperiod_statistics(iperiod)%mask(idata)
              sum1=sum1+stat_cen(istat_cen)%data(iidata)%fg_depar
           ENDDO
           mean=sum1/REAL(ndata,JPRB)
           sum2=0.
           DO idata=1,ndata
              iidata=stat_cen(istat_cen)%tperiod_statistics(iperiod)%mask(idata)
              sum2=sum2+(stat_cen(istat_cen)%data(iidata)%fg_depar-mean)**2.
           ENDDO
           IF (ndata .gt. 2) THEN
              stdev=sqrt(sum2/REAL(ndata-1,JPRB))
           ELSE
              stdev=0
           ENDIF
           nperiod_ok=nperiod_ok+1
           tab_mean(nperiod_ok)= mean
           tab_stdev(nperiod_ok)= stdev
           IF (iperiod .le. nperiod_10days) &
              stat_cen(istat_cen)%statistics%ndata_10days = &
            & stat_cen(istat_cen)%statistics%ndata_10days + ndata
        ELSE
           mean=0.
           stdev=0.
        ENDIF
        stat_cen(istat_cen)%tperiod_statistics(iperiod)%mean = mean
        stat_cen(istat_cen)%tperiod_statistics(iperiod)%stdev = stdev
        stat_cen(istat_cen)%tperiod_statistics(iperiod)%time_coverage7 = &
        & REAL( &
        &SUM(stat_cen(istat_cen)%tperiod_statistics(iperiod)%data_present7(:))&
        & ,JPRB)/ &
        & REAL(nslot_max)*100.
     ENDDO
     IF (nperiod_ok .gt. 2) THEN
        tab2_mean=STATS(nperiod_ok,tab_mean(1:nperiod_ok))
        tab2_stdev=STATS(nperiod_ok,tab_stdev(1:nperiod_ok))
     ELSE
        tab2_mean=(/9999.,9999./)
        tab2_stdev=(/9999.,9999./)
     ENDIF
     mean_min=minval(tab_mean(1:nperiod_ok))
     mean_max=maxval(tab_mean(1:nperiod_ok))
     stdev_min=minval(tab_stdev(1:nperiod_ok))
     stdev_max=maxval(tab_stdev(1:nperiod_ok))
     stat_cen(istat_cen)%nperiod_ok=nperiod_ok
     stat_cen(istat_cen)%mean_min=mean_min
     stat_cen(istat_cen)%mean_max=mean_max
     stat_cen(istat_cen)%mean_mean=tab2_mean(1)
     stat_cen(istat_cen)%mean_stdev=tab2_mean(2)
     stat_cen(istat_cen)%stdev_min=stdev_min
     stat_cen(istat_cen)%stdev_max=stdev_max
     stat_cen(istat_cen)%stdev_mean=tab2_stdev(1)
     stat_cen(istat_cen)%stdev_stdev=tab2_stdev(2)
     
     IF (stat_cen(istat_cen)%statistics%ndata_10days .gt. 2) THEN
        tab_to_test = stat_cen(istat_cen)%data%fg_depar
        tab2_10days=STATS(stat_cen(istat_cen)%statistics%ndata_10days, &
         & tab_to_test(1:stat_cen(istat_cen)%statistics%ndata_10days))
     ELSE
        tab2_10days=(/9999.,9999./)
     ENDIF
     stat_cen(istat_cen)%statistics%mean_10days=tab2_10days(1)
     stat_cen(istat_cen)%statistics%stdev_10days=tab2_10days(2)

  ENDDO

END SUBROUTINE CALC_STATISTICS_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CALC_STATISTICS_ACTIVE_GPSSOL(stat_cen_active, &
  statistics_active, error_flag)
  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod
  implicit none

  ! I/O variables
  TYPE(STAT_CEN_ACTIVE_TYPE), dimension(1:nstat_cen_active_max), &
                      intent(inout) :: stat_cen_active
  TYPE(STATISTICS_ACTIVE_TYPE), intent(inout) :: statistics_active
  integer(KIND=JPIM), intent(out) :: error_flag

  ! local variables
  integer(KIND=JPIM) :: istat_cen,ndata,ndata_rejected,ndata0,idata,idata0,islot,ndata_all
  integer(KIND=JPIM), dimension(0:nslot_max) :: jdata
  real(KIND=JPRB), dimension(0:nslot_max,1:ndata_max) :: tab_hxb_hxa_all, tab_y0_hxb_all, tab_y0_hxa_all
  real(KIND=JPRB), dimension(1:ndata_max) :: tab_hxb_hxa, tab_y0_hxb, tab_y0_hxa
  real(KIND=JPRB), dimension(1:2) :: tab2_hxb_hxa_all, tab2_y0_hxb_all, tab2_y0_hxa_all
  real(KIND=JPRB), dimension(1:2) :: tab2_hxb_hxa, tab2_y0_hxb, tab2_y0_hxa
  real(KIND=JPRB)    :: sum2b,sum2a,sum2o

  statistics_active%ndata(:)=0
  statistics_active%ndata_rejected(:)=0
  jdata(:)=0
  DO istat_cen=1,nstat_cen_active

     ! statistics for all the timeslots
     ndata0=stat_cen_active(istat_cen)%statistics%ndata(0)

     DO islot=0,nslot_max
       ndata=stat_cen_active(istat_cen)%statistics%ndata(islot)
       ndata_rejected=stat_cen_active(istat_cen)%statistics%ndata_rejected(islot)
       IF (ndata .gt. 0) THEN
         idata=0
         statistics_active%ndata(islot)=statistics_active%ndata(islot) + ndata
         statistics_active%ndata_rejected(islot)=statistics_active%ndata_rejected(islot) + ndata_rejected
         DO idata0=1,ndata0
           IF (islot .eq. 0 .or. &
             & stat_cen_active(istat_cen)%data(idata0)%tslot .eq. islot) THEN
           idata=idata+1
           tab_hxb_hxa(idata)=stat_cen_active(istat_cen)%data(idata0)%an_depar - &
                            & stat_cen_active(istat_cen)%data(idata0)%fg_depar
           tab_y0_hxb(idata) =stat_cen_active(istat_cen)%data(idata0)%fg_depar
           tab_y0_hxa(idata) =stat_cen_active(istat_cen)%data(idata0)%an_depar
           jdata(islot)=jdata(islot)+1
           tab_hxb_hxa_all(islot,jdata(islot))=stat_cen_active(istat_cen)%data(idata0)%an_depar - &
                                             & stat_cen_active(istat_cen)%data(idata0)%fg_depar
           tab_y0_hxb_all(islot,jdata(islot)) =stat_cen_active(istat_cen)%data(idata0)%fg_depar
           tab_y0_hxa_all(islot,jdata(islot)) =stat_cen_active(istat_cen)%data(idata0)%an_depar
           ENDIF
         ENDDO
       ENDIF

       IF (ndata .gt. 2) THEN
         tab2_hxb_hxa=STATS(ndata,tab_hxb_hxa(1:ndata))
         tab2_y0_hxb =STATS(ndata,tab_y0_hxb (1:ndata))
         tab2_y0_hxa =STATS(ndata,tab_y0_hxa (1:ndata))
         stat_cen_active(istat_cen)%statistics%mean_hxb_hxa (islot)= &
           tab2_hxb_hxa(1)
         stat_cen_active(istat_cen)%statistics%stdev_hxb_hxa(islot)= &
           tab2_hxb_hxa(2)
         stat_cen_active(istat_cen)%statistics%mean_y0_hxb  (islot)= &
           tab2_y0_hxb (1)
         stat_cen_active(istat_cen)%statistics%stdev_y0_hxb (islot)= &
         tab2_y0_hxb (2)
         stat_cen_active(istat_cen)%statistics%mean_y0_hxa  (islot)= &
           tab2_y0_hxa (1)
         stat_cen_active(istat_cen)%statistics%stdev_y0_hxa (islot)= &
           tab2_y0_hxa (2)
  
         sum2a=0.
         sum2b=0.
         sum2o=0.
         DO idata=1,ndata
            sum2a=sum2a+abs( (tab_hxb_hxa(idata)-tab2_hxb_hxa(1))* &
                      & (tab_y0_hxa (idata)-tab2_y0_hxa (1)) )
            sum2b=sum2b+abs( (tab_hxb_hxa(idata)-tab2_hxb_hxa(1))* &
                      & (tab_y0_hxb (idata)-tab2_y0_hxb (1)) )
            sum2o=sum2o+abs( (tab_y0_hxb (idata)-tab2_y0_hxb (1))* &
                      & (tab_y0_hxa (idata)-tab2_y0_hxa (1)) )
         ENDDO
         stat_cen_active(istat_cen)%statistics%sigb(islot)=sqrt(sum2b/REAL(ndata,JPRB))
         stat_cen_active(istat_cen)%statistics%siga(islot)=sqrt(sum2a/REAL(ndata,JPRB))
         stat_cen_active(istat_cen)%statistics%sigo(islot)=sqrt(sum2o/REAL(ndata,JPRB))

       ELSE

         stat_cen_active(istat_cen)%statistics%mean_hxb_hxa (islot)=9999.
         stat_cen_active(istat_cen)%statistics%stdev_hxb_hxa(islot)=9999.
         stat_cen_active(istat_cen)%statistics%mean_y0_hxb  (islot)=9999.
         stat_cen_active(istat_cen)%statistics%stdev_y0_hxb (islot)=9999.
         stat_cen_active(istat_cen)%statistics%mean_y0_hxa  (islot)=9999.
         stat_cen_active(istat_cen)%statistics%stdev_y0_hxa (islot)=9999.
         stat_cen_active(istat_cen)%statistics%sigb(islot)=9999.
         stat_cen_active(istat_cen)%statistics%siga(islot)=9999.
         stat_cen_active(istat_cen)%statistics%sigo(islot)=9999.

       ENDIF
       
     ENDDO

  ENDDO

  DO islot=0,nslot_max
     IF (statistics_active%ndata(islot) .ne. jdata(islot)) THEN
        print *,'WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *,'statistics_active%ndata(islot) .ne. jdata(islot) for islot=',islot
        print *,'statistics_active%ndata(islot)',statistics_active%ndata(islot)
        print *,'jdata(islot)',jdata(islot)
     ENDIF

     ndata_all=statistics_active%ndata(islot)
     tab2_hxb_hxa_all=STATS(ndata_all,tab_hxb_hxa_all(islot,1:ndata_all))
     tab2_y0_hxb_all =STATS(ndata_all,tab_y0_hxb_all (islot,1:ndata_all))
     tab2_y0_hxa_all =STATS(ndata_all,tab_y0_hxa_all (islot,1:ndata_all))
     statistics_active%mean_hxb_hxa (islot)= tab2_hxb_hxa_all(1)
     statistics_active%stdev_hxb_hxa(islot)= tab2_hxb_hxa_all(2)
     statistics_active%mean_y0_hxb  (islot)= tab2_y0_hxb_all (1)
     statistics_active%stdev_y0_hxb (islot)= tab2_y0_hxb_all (2)
     statistics_active%mean_y0_hxa  (islot)= tab2_y0_hxa_all (1)
     statistics_active%stdev_y0_hxa (islot)= tab2_y0_hxa_all (2)
     sum2a=0.
     sum2b=0.
     sum2o=0.
     DO idata=1,ndata_all
        sum2a=sum2a+abs( (tab_hxb_hxa_all(islot,idata)-tab2_hxb_hxa_all(1))* &
                  & (tab_y0_hxa_all (islot,idata)-tab2_y0_hxa_all (1)) )
        sum2b=sum2b+abs( (tab_hxb_hxa_all(islot,idata)-tab2_hxb_hxa_all(1))* &
                  & (tab_y0_hxb_all (islot,idata)-tab2_y0_hxb_all (1)) )
        sum2o=sum2o+abs( (tab_y0_hxb_all (islot,idata)-tab2_y0_hxb_all (1))* &
                  & (tab_y0_hxa_all (islot,idata)-tab2_y0_hxa_all (1)) )
     ENDDO
     statistics_active%sigb(islot)=sqrt(sum2b/REAL(ndata_all,JPRB))
     statistics_active%siga(islot)=sqrt(sum2a/REAL(ndata_all,JPRB))
     statistics_active%sigo(islot)=sqrt(sum2o/REAL(ndata_all,JPRB))
  ENDDO

 END SUBROUTINE CALC_STATISTICS_ACTIVE_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION IS_DISTRIBUTION_GAUSSIAN(ndata_in,tab,mean,stdev) result(ll_gauss)
  USE parkind1, only: JPRB, JPIM
  USE gpssol_mod, only: ndata_min_in_class_gauss, nclass_gauss, nconf_gauss, &
                        density_method, kde_P
  implicit none

  ! I/O variables
  integer(KIND=JPIM), intent(in) :: ndata_in
  real(KIND=JPRB), dimension(1:ndata_in), intent(in) :: tab
  real(KIND=JPRB), intent(in) :: mean, stdev
  logical :: ll_gauss
  
  ! local variables
  integer(KIND=JPIM) :: nstep, i, j, threshold, ngauss, nd, nu
  real(KIND=JPRB)    :: xmin, xmax, dstep
  real(KIND=JPRB), dimension(:), allocatable :: xb, gauss, y2
  real(KIND=JPRB), dimension(:), allocatable :: hh_histo, hh_kde, kde_input,kde_output, hh_in
  real(KIND=JPRB), dimension(1:50,1:3)    :: CHI2_0
  real(KIND=JPRB), dimension(1:3) :: B=(/1.,1./6.,.1/3./)
  real(KIND=JPRB) :: P_loc, B_loc, H, ndata
  
  CHI2_0( 1,1:3)=(/3.84,	6.64,	10.83/)
  CHI2_0( 2,1:3)=(/5.99 ,9.21,	13.82/)
  CHI2_0( 3,1:3)=(/7.82 ,11.35,	16.27/)
  CHI2_0( 4,1:3)=(/9.49 ,13.28,	18.47/)
  CHI2_0( 5,1:3)=(/11.07 ,15.09,	20.52/)
  CHI2_0( 6,1:3)=(/12.59 ,16.81,	22.46/)
  CHI2_0( 7,1:3)=(/14.07 ,18.48,	24.32/)
  CHI2_0( 8,1:3)=(/15.51 ,20.09,	26.13/)
  CHI2_0( 9,1:3)=(/16.92 ,21.67,	27.88/)
  CHI2_0(10,1:3)=(/18.31 ,23.21,	29.59/)
  CHI2_0(11,1:3)=(/19.68 ,24.73,	31.26/)
  CHI2_0(12,1:3)=(/21.03 ,26.22,	32.91/)
  CHI2_0(13,1:3)=(/22.36 ,27.69,	34.53/)
  CHI2_0(14,1:3)=(/23.69 ,29.14,	36.12/)
  CHI2_0(15,1:3)=(/25.00 ,30.58,	37.70/)
  CHI2_0(16,1:3)=(/26.30 ,32.00,	39.25/)
  CHI2_0(17,1:3)=(/27.59 ,33.41,	40.79/)
  CHI2_0(18,1:3)=(/28.87 ,34.81,	42.31/)
  CHI2_0(19,1:3)=(/30.14 ,36.19,	43.82/)
  CHI2_0(20,1:3)=(/31.41 ,37.57,	45.32/)
  CHI2_0(21,1:3)=(/32.67 ,38.93,	46.80/)
  CHI2_0(22,1:3)=(/33.92 ,40.29,	48.27/)
  CHI2_0(23,1:3)=(/35.17 ,41.64,	49.73/)
  CHI2_0(24,1:3)=(/36.42 ,42.98,	51.18/)
  CHI2_0(25,1:3)=(/37.65 ,44.31,	52.62/)
  CHI2_0(26,1:3)=(/38.89 ,45.64,	54.05/)
  CHI2_0(27,1:3)=(/40.11 ,46.96,	55.48/)
  CHI2_0(28,1:3)=(/41.34 ,48.28,	56.89/)
  CHI2_0(29,1:3)=(/42.56 ,49.59,	58.30/)
  CHI2_0(30,1:3)=(/43.77 ,50.89,	59.70/)
  CHI2_0(31,1:3)=(/44.99,	52.19,	61.10/)
  CHI2_0(32,1:3)=(/46.19,	53.49,	62.49/)
  CHI2_0(33,1:3)=(/47.40,	54.78,	63.87/)
  CHI2_0(34,1:3)=(/48.60,	56.06,	65.25/)
  CHI2_0(35,1:3)=(/49.80,	57.34,	66.62/)
  CHI2_0(36,1:3)=(/51.00,	58.62,	67.99/)
  CHI2_0(37,1:3)=(/52.19,	59.89,	69.35/)
  CHI2_0(38,1:3)=(/53.38,	61.16,	70.71/)
  CHI2_0(39,1:3)=(/54.57,	62.43,	72.06/)
  CHI2_0(40,1:3)=(/55.76,	63.69,	73.41/)
  CHI2_0(41,1:3)=(/56.94,	64.95,	74.75/)
  CHI2_0(42,1:3)=(/58.12,	66.21,	76.09/)
  CHI2_0(43,1:3)=(/59.30,	67.46,	77.42/)
  CHI2_0(44,1:3)=(/60.48,	68.71,	78.75/)
  CHI2_0(45,1:3)=(/61.66,	69.96,	80.08/)
  CHI2_0(46,1:3)=(/62.83,	71.20,	81.40/)
  CHI2_0(47,1:3)=(/64.00,	72.44,	82.72/)
  CHI2_0(48,1:3)=(/65.17,	73.68,	84.03/)
  CHI2_0(49,1:3)=(/66.34,	74.92,	85.35/)
  CHI2_0(50,1:3)=(/67.51,	76.15,	86.66/)

  threshold=ndata_min_in_class_gauss
  nstep=nclass_gauss

  xmin=minval(tab)
  xmax=maxval(tab)
  dstep=(xmax-xmin)/REAL(nstep,JPRB)
  allocate(xb(nstep),gauss(nstep),hh_histo(nstep))
  DO i=1,nstep
     xb(i)=REAL(i-1,JPRB)*dstep+xmin+dstep/2.
  ENDDO
!  print *,'ndata_in',ndata_in
!  print *,'ndata ',ndata
!  ndata=ndata_in

  ! create histogram
  DO i=1,nstep
     xmin=xb(i)-dstep/2.
     xmax=xb(i)+dstep/2.
     hh_histo(i)=REAL(COUNT(tab .ge. xmin .and. tab .lt. xmax),JPRB)
  ENDDO
!  print *,'hh_histo'
!  print '(10I4)',INT(hh_histo)
!  print *,'estimated ndata from hh_histo',SUM(hh_histo)
  ! use Kernel Density Estimator
  P_loc=REAL(kde_P)
  H=0.003
  B_loc=B(INT(P_loc+1))
  allocate(kde_input(nstep),kde_output(nstep),hh_kde(nstep))
  hh_kde(:)=0.
  DO i=1,ndata_in
     kde_input=(xb(:)-tab(i))/H
     kde_output(:)=0.
     DO j=1,nstep
        IF (kde_input(j) .ge. -1. .and. kde_input(j) .lt. 1.) THEN
           kde_output(j)=((1.-kde_input(j)**2.)**P_loc)/(2.**(2.*P_loc+1.))/B_loc
        ENDIF
     ENDDO
!     print *,'kde_output(1)',kde_output(j)*dstep/H
     hh_kde(:)=hh_kde(:)+kde_output(:)*dstep/H
  ENDDO
  deallocate(kde_input,kde_output)
!  print *,'hh_kde'
!  print '(10F7.2)',hh_kde
!  print *,'estimated ndata from hh_kde',SUM(hh_kde)

  ! create expected gauss distribution
  allocate(hh_in(nstep))
  IF (trim(density_method) .eq. 'kde') THEN
     hh_in=hh_kde(:)
  ELSEIF (trim(density_method) .eq. 'histo') THEN
     hh_in=hh_histo(:)
  ELSE
     print *,'density_method unknown',trim(density_method)
     stop
  ENDIF
  ndata=SUM(hh_in(:))
!  print *,'ndata',ndata
  ngauss=nstep
  gauss=exp(-(xb-mean)**2./(2.*stdev**2.))/sqrt(2.*3.14159265)/stdev*dstep*ndata

!  print *,'xb',xb
!  print *,'hh_histo :'
!  print '(10I4)',INT(hh_histo)
!  print *,'mean',mean
!  print *,'stdev',stdev
!  print '(10I4)', INT(gauss)
!  print *,'merge'

! merges the cells of arrays NDATA and TAB so that all cells
! in NDATA contain at least <<THRESHOLD>> elements

  i=0
  DO WHILE (i .lt. ngauss-1)
    i=i+1
    DO WHILE (gauss(i) .lt. threshold .and. i .lt. ngauss)
      nd=ngauss
      gauss(i)=gauss(i)+gauss(i+1)
      IF (i .ne. nd-1) gauss(i+1:nd-1)=gauss(i+2:nd)
      ngauss=ngauss-1
      hh_in(i)=hh_in(i)+hh_in(i+1)
      IF (i .ne. nd) hh_in(i+1:nd-1)=hh_in(i+2:nd)
    END DO
  END DO

!  print *,'in-merge hh_in'
!  print '(10I4)',INT(hh_in(1:ngauss))

  i=ngauss+1
  DO WHILE (i .gt. 2)
    i=i-1
    DO WHILE (gauss(i) .lt. threshold .and. i .gt. 1)
      nd=ngauss
      gauss(i)=gauss(i)+gauss(i-1)
      IF (i .ne. 2) gauss(2:i-1)=gauss(1:i-2)
      gauss(1:nd-1)=gauss(2:nd)
      ngauss=ngauss-1
      hh_in(i)=hh_in(i)+hh_in(i-1)
      IF (i .ne. 2) hh_in(2:i-1)=hh_in(1:i-2)
      hh_in(1:nd-1)=hh_in(2:nd)
      i=i-1
    END DO
  END DO

!  print *,'post-merge hh_in'
!  print '(10I4)',INT(hh_in(1:ngauss))

!  print '(10I4)', INT(gauss(1:ngauss))

  allocate(y2(ngauss))
  y2=(hh_in(1:ngauss)-gauss(1:ngauss))**2./gauss(1:ngauss)
  nu=ngauss-1-2
  IF (nu .le. 1) THEN
!  IF (nu .le. 5) THEN
!    print *,'not enough data ',ngauss,nu
    ll_gauss=.false.
  ELSE
!    print *,'nu ',nu,'sum hh_in',SUM(hh_in(1:ngauss))
!    print *,'Chi 2 = ',SUM(y2)
!    print *,'Threshold = ',CHI2_0(nu,nconf_gauss)
!    print *,'H0 = Distribution drawn from a Gaussian distribution'
!    print *,'Reject H0 at 95%   ?',SUM(y2) .ge. CHI2_0(nu,nconf_gauss)
    ll_gauss=SUM(y2) .lt. CHI2_0(nu,nconf_gauss)
  ENDIF
!  print *,'over'
  deallocate(xb,gauss,y2,hh_histo,hh_kde,hh_in)

END FUNCTION IS_DISTRIBUTION_GAUSSIAN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SELECT_GPSSOL_STATIONS(stat_cen,error_flag)
  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod
  implicit none

  ! I/O variables
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max), intent(inout) :: stat_cen
  integer(KIND=JPIM), intent(out) :: error_flag

  ! local variables
  logical :: flagok
  character(len=4) :: stat4
  integer(KIND=JPIM), dimension(0:15) :: i_id
  integer(KIND=JPIM) :: istat_cen, icount, jstat_cen, ipbs, npbs, &
                     stats_pbs(1:nstat_cen_max), idt, dt1, &
                     counts(13), icen
  real(KIND=JPRB) :: stdevmin, diffmin, mean_stdevmin, &
                     dist1, pi=3.14159265358979323846264338328, &
                     deg2rad, rad2deg, &
                     lat_1, lon_1, lat_2, lon_2, cost, &
                     dists(1:nstat_cen_max,1:nstat_cen_max)
  character(len=2) :: fmt2
  character(len=150) :: chain_stat

  counts(:)=0
  DO istat_cen=1,nstat_cen
     idt=1
     DO WHILE (stat_cen(istat_cen)%ndts(idt) .ne. &
       & maxval(stat_cen(istat_cen)%ndts) .and. &
       & idt .lt. stat_cen(istat_cen)%ndt)
       idt=idt+1
     ENDDO
     dt1=stat_cen(istat_cen)%dts(idt)
     stat_cen(istat_cen)%ok = &
!      & (stat_cen(istat_cen)%statistics%time_coverage7_10days .gt. min_time_coverage) .and. &
!      & (stat_cen(istat_cen)%statistics%time_coverage_10days .gt. min_time_coverage) .and. &
      & (stat_cen(istat_cen)%statistics%time_coverage7 .gt. min_time_coverage) .and. &
!      & (stat_cen(istat_cen)%statistics%time_coverage .gt. min_time_coverage) .and. &
      & (stat_cen(istat_cen)%statistics%gauss) .and. &
      & (dt1 .lt. max_dt) .and. &
!      & (stat_cen(istat_cen)%namestr(5:7) .ne. 'SGN') .and. &
!      & (stat_cen(istat_cen)%namestr(5:7) .ne. 'GOP') .and. &
!      & (stat_cen(istat_cen)%namestr(5:7) .ne. 'ROB') .and. &
!      & (stat_cen(istat_cen)%namestr(5:7) .ne. 'BKG') .and. &
      & (abs(stat_cen(istat_cen)%alt_mod - stat_cen(istat_cen)%alt) .le. diff_alt_max .and. &
      &  stat_cen(istat_cen)%alt .le. alt_max) .and. &
      & (stat_cen(istat_cen)%lon .ge.  lonmin .and. &  ! conver from degree to arc
      &  stat_cen(istat_cen)%lon .le.  lonmax )          .and. &
      & (stat_cen(istat_cen)%lat .ge.  latmin .and. &
      &  stat_cen(istat_cen)%lat .le.  latmax )           .and. &
!      & (stat_cen(istat_cen)%ndt .eq. 1) .and. &
!      & (stat_cen(istat_cen)%dts(1) .ne. 60) .and. &
      & abs(stat_cen(istat_cen)%statistics%mean_10days) < biasmax .and. &
      & abs(stat_cen(istat_cen)%statistics%mean) < biasmax .and. &
      & stat_cen(istat_cen)%statistics%stdev_10days < stdevmax .and. &
      & stat_cen(istat_cen)%statistics%stdev < stdevmax .and. &
      & (stat_cen(istat_cen)%cst_alt) !.or. stat_cen(istat_cen)%namestr(5:7) .eq.  'ACR')
!      IF (.not.(stat_cen(istat_cen)%statistics%time_coverage7_10days .gt. min_time_coverage)) counts(1)=counts(1)+1
!      IF (.not.(stat_cen(istat_cen)%statistics%time_coverage_10days .gt. min_time_coverage)) counts(2)=counts(2)+1
      IF (n_cenexclude >= 1) THEN
         DO icen=1,n_cenexclude
            IF (cen_exclude(icen) == stat_cen(istat_cen)%namestr(5:8)) THEN
                    stat_cen(istat_cen)%ok = .FALSE.
                    counts(10)=counts(10)+1
                    print*,'GREP EXCLUDE',stat_cen(istat_cen)%namestr
            ENDIF
         ENDDO
      ENDIF
      IF (.not.(stat_cen(istat_cen)%statistics%time_coverage7 .gt. min_time_coverage)) then
       counts(3)=counts(3)+1
       PRINT*,'GREP RARE:',stat_cen(istat_cen)%namestr
      ENDIF
!      IF (.not.(stat_cen(istat_cen)%statistics%time_coverage .gt. min_time_coverage)) counts(4)=counts(4)+1
      IF (.not.(stat_cen(istat_cen)%statistics%gauss)) THEN
        counts(5)=counts(5)+1
        PRINT*,"GREP NOT GAUSS",stat_cen(istat_cen)%namestr
      ENDIF
      IF (.not.(dt1 .lt. max_dt)) THEN
        counts(6)=counts(6)+1
        PRINT*,'GREP timeshift',stat_cen(istat_cen)%namestr
      ENDIF
      IF (.not.(abs(stat_cen(istat_cen)%alt_mod - stat_cen(istat_cen)%alt) .le. diff_alt_max .and. &
      &         stat_cen(istat_cen)%alt .le. alt_max)) THEN
        counts(7)=counts(7)+1
        PRINT*,'GREP altitude',stat_cen(istat_cen)%namestr
      ENDIF
      ! lat, lon from odb results are in arc, not in degree! so convert with 2pi/360=0.0174
      IF (.not.((stat_cen(istat_cen)%lon .ge. lonmin .and. &
      &  stat_cen(istat_cen)%lon .le.  lonmax )          .and. &
      & (stat_cen(istat_cen)%lat .ge.  latmin .and. &
      &  stat_cen(istat_cen)%lat .le.  latmax))) THEN
        counts(8)=counts(8)+1
        PRINT*,'GREP out of domain:',stat_cen(istat_cen)%namestr
      ENDIF
      ! print *,stat_cen(istat_cen)%namestr,stat_cen(istat_cen)%lon,stat_cen(istat_cen)%lat, lonmin, lonmax,latmin,latmax
      IF (.not.(stat_cen(istat_cen)%cst_alt)) THEN
        counts(9)=counts(9)+1
         PRINT*,"GREP ALTMOVING",stat_cen(istat_cen)%namestr
      ENDIF
      IF (.not.(abs(stat_cen(istat_cen)%statistics%mean_10days) < biasmax .and. &
              & abs(stat_cen(istat_cen)%statistics%mean       ) < biasmax))THEN
        counts(11)=counts(11)+1
        PRINT*," GREP BIASTOBIG",stat_cen(istat_cen)%namestr,stat_cen(istat_cen)%statistics%mean_10days
        PRINT*,stat_cen(istat_cen)%statistics%mean,biasmax
      ENDIF
      IF (.not.(stat_cen(istat_cen)%statistics%stdev_10days < stdevmax .and. &
              & stat_cen(istat_cen)%statistics%stdev        < stdevmax)) then
          counts(12)=counts(12)+1
          PRINT*," GREP STDVTOBIG",stat_cen(istat_cen)%namestr
      ENDIF
      IF (stat_cen(istat_cen)%ok) THEN
        counts(13)=counts(13)+1
        !PRINT*,"GREP OK",stat_cen(istat_cen)%namestr
      ENDIF
  ENDDO

  print *,'Rejected counts: out of ',nstat_cen
  print *,'Time_coverage7_10days :',counts(1)
  print *,'Time_coverage_10days  :',counts(2)
  print *,'Time_coverage7        :',counts(3)
  print *,'Time_coverage         :',counts(4)
  print *,'Gauss                 :',counts(5)
  print *,'Time averaging        :',counts(6)
  print *,'Altitude              :',counts(7)
  print *,'Domain                :',counts(8)
  print *,'Altitude not constant :',counts(9)
  print *,'Analysis center excl  :',counts(10)
  print *,'Bias above max        :',counts(11)
  print *,'Stdev above max       :',counts(12)
  print *,'-------------------------------'
  print *,'Total                 :',counts(13)

  call PRINT_STATISTICS(stat_cen,'stats_1.txt'//trim(ficoutput_ext), &
                        okonly=.TRUE.)

  print *,'select_gpssol_stations: Select one processing center per station'
  DO istat_cen=1,nstat_cen
     ! look if there are several stations which are OK but with different
     ! processings
     stat4=stat_cen(istat_cen)%namestr(1:4)
     i_id=FIND_CHAIN(nstat_cen,4,nstat_cen_max,stat_cen%namestr(1:4),stat4, &
          MASK=stat_cen%ok)
     ! select a unique station from a pack
     IF (i_id(0) .gt. 1) THEN
        write(fmt2,'(I2)') i_id(0)
        chain_stat(:)=' '
        write(chain_stat,'('//fmt2//'A10)') stat_cen(i_id(1:i_id(0)))%namestr
        ! select using stdev
        IF ((maxval(stat_cen(i_id(1:i_id(0)))%statistics%stdev)- &
          & minval(stat_cen(i_id(1:i_id(0)))%statistics%stdev))*1.e3 .gt. 0.2) &
        THEN
           print *,'stdev method:'//trim(chain_stat)
           stdevmin=minval(stat_cen(i_id(1:i_id(0)))%statistics%stdev)
           DO icount=1,i_id(0)
              IF (stat_cen(i_id(icount))%statistics%stdev .eq. stdevmin) THEN
                 stat_cen(i_id(icount))%ok=.TRUE.
              ELSE
                 stat_cen(i_id(icount))%ok=.FALSE.
              ENDIF
           ENDDO
        ELSE
        ! select using difference mean_max-mean_min
           IF ((maxval( &
   & stat_cen(i_id(1:i_id(0)))%mean_max-stat_cen(i_id(1:i_id(0)))%mean_min) - &
               minval( &
   & stat_cen(i_id(1:i_id(0)))%mean_max-stat_cen(i_id(1:i_id(0)))%mean_min))* &
   & 1.e3 .gt. 5.) THEN
              print *,'diff method:'//trim(chain_stat)
              diffmin=minval( &
   & stat_cen(i_id(1:i_id(0)))%mean_max-stat_cen(i_id(1:i_id(0)))%mean_min)
 ! -------- Debut modif Paul poli Avr 2008
              !DO icount=1,i_id(0)
              !IF (stat_cen(i_id(icount))%mean_max- &
              !  & stat_cen(i_id(icount))%mean_min .eq. diffmin) THEN
              !   stat_cen(i_id(icount))%ok=.TRUE.
              !ELSE
              !   stat_cen(i_id(icount))%ok=.FALSE.
              !ENDIF
              !ENDDO
              print *,'series of diff:', &
               & stat_cen(i_id(1:i_id(0)))%mean_max - stat_cen(i_id(1:i_id(0)))%mean_min
             print *,'minimum diff:', diffmin
             DO icount=1,i_id(0)
             IF (stat_cen(i_id(icount))%mean_max- &
               & stat_cen(i_id(icount))%mean_min .le. diffmin+1.e-7) THEN
                stat_cen(i_id(icount))%ok=.TRUE.
             ELSE
                stat_cen(i_id(icount))%ok=.FALSE.
             ENDIF
             print *,'icount #',icount,'found diff=', stat_cen(i_id(icount))%mean_max- &
               & stat_cen(i_id(icount))%mean_min, 'compare with min=', diffmin, 'result:',stat_cen(i_id(icount))%ok
             ENDDO
 ! -------- Fin modif Paul poli Avr 2008
           ELSE
           ! select using stdev of mean
              print *,'stdev of mean method:'//trim(chain_stat)
              mean_stdevmin=minval(stat_cen(i_id(1:i_id(0)))%mean_stdev)
              DO icount=1,i_id(0)
              IF (stat_cen(i_id(icount))%mean_stdev .eq. mean_stdevmin) THEN
                 stat_cen(i_id(icount))%ok=.TRUE.
              ELSE
                 stat_cen(i_id(icount))%ok=.FALSE.
              ENDIF
              ENDDO
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  call PRINT_STATISTICS(stat_cen,'stats_2.txt'//trim(ficoutput_ext), &
                        okonly=.TRUE.)

  ! calculate stations inter-distances
  print *,'select_gpssol_stations: Horizontal thinning to ',INT(min_distance),' km'
  deg2rad=pi/180.
  rad2deg=180./pi
  DO istat_cen=1,nstat_cen-1
     DO jstat_cen=istat_cen+1,nstat_cen
        lat_1=stat_cen(istat_cen)%lat  !*deg2rad  in rad already, so remove conversion
        lon_1=stat_cen(istat_cen)%lon  !*deg2rad
        lat_2=stat_cen(jstat_cen)%lat  !*deg2rad
        lon_2=stat_cen(jstat_cen)%lon  !*deg2rad
        cost=sin(lat_2)*sin(lat_1)+cos(lat_2)*cos(lat_1)*cos(lon_2-lon_1)
        dist1=6378.*acos(cost)
        dists(istat_cen,jstat_cen)=dist1
     ENDDO
  ENDDO

  ! select problematic pairs (stations closer than xx km)
  npbs=0
  DO istat_cen=1,nstat_cen-1
     DO jstat_cen=istat_cen+1,nstat_cen
        IF (stat_cen(istat_cen)%ok .and. &
          & stat_cen(jstat_cen)%ok .and. &
          & dists(istat_cen,jstat_cen) .lt. min_distance) THEN
            npbs=npbs+1
            stats_pbs((npbs-1)*2+1)=istat_cen
            stats_pbs(npbs*2)=jstat_cen
            print *,'stations ', &
              stat_cen(istat_cen)%namestr,' ', &
              stat_cen(jstat_cen)%namestr, dists(istat_cen,jstat_cen)
        ENDIF
     ENDDO
  ENDDO

  ! remove stations which are present more than once from the list
  DO ipbs=1,npbs*2
     i_id(0)=COUNT(stats_pbs .eq. stats_pbs(ipbs))
     IF (i_id(0) .gt. 1) &
     stat_cen(stats_pbs(ipbs))%ok=.FALSE.
  ENDDO

  call PRINT_STATISTICS(stat_cen,'stats_3.txt'//trim(ficoutput_ext), &
                        okonly=.TRUE.)

  ! select again problematic pairs
  npbs=0
  DO istat_cen=1,nstat_cen-1
     DO jstat_cen=istat_cen+1,nstat_cen
        IF (stat_cen(istat_cen)%ok .and. &
          & stat_cen(jstat_cen)%ok .and. &
          & dists(istat_cen,jstat_cen) .lt. min_distance) THEN
            npbs=npbs+1
            stats_pbs((npbs-1)*2+1)=istat_cen
            stats_pbs(npbs*2)=jstat_cen
            print *,'stations ', &
              stat_cen(istat_cen)%namestr,' ', &
              stat_cen(jstat_cen)%namestr, dists(istat_cen,jstat_cen)
        ENDIF
     ENDDO
  ENDDO

  ! remove the problematic station in each pair. Use the stdev method.
  DO ipbs=1,npbs
     istat_cen=stats_pbs((ipbs-1)*2+1)
     jstat_cen=stats_pbs(ipbs*2)
     IF (stat_cen(istat_cen)%statistics%stdev .gt. &
       & stat_cen(jstat_cen)%statistics%stdev) THEN
        stat_cen(istat_cen)%ok=.FALSE.
     ELSE
        stat_cen(jstat_cen)%ok=.FALSE.
     ENDIF
  ENDDO

  call PRINT_STATISTICS(stat_cen,'stats_4.txt'//trim(ficoutput_ext), &
                        okonly=.TRUE.)

END SUBROUTINE SELECT_GPSSOL_STATIONS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_LIST_GPSSOL(list_gpssol,stat_cen,error_flag)
  USE parkind1, ONLY:JPRB,JPIM
  USE gpssol_mod
  implicit none

  ! I/O variables
  character(len=*), intent(in)   :: list_gpssol
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max), intent(inout) :: stat_cen
  integer(KIND=JPIM), intent(out) :: error_flag

  ! local variables
  integer(KIND=JPIM) :: iostatus, istat_cen, ns, idt
  real(KIND=JPRB)    :: dt1, sigma_o

  error_flag = no_error
  print '(A)','write_list_gpssol: Opening file '
  print '(A)',trim(list_gpssol)
  OPEN(unit=55,file=trim(list_gpssol),action='WRITE',iostat=iostatus)
  ns=0
  IF (iostatus == 0) THEN
    DO istat_cen=1,nstat_cen
       IF (stat_cen(istat_cen)%ok) THEN
          idt=1
          DO WHILE (stat_cen(istat_cen)%ndts(idt) .ne. &
             & maxval(stat_cen(istat_cen)%ndts) .and. &
             & idt .lt. stat_cen(istat_cen)%ndt)
             idt=idt+1
          ENDDO
          dt1=REAL(stat_cen(istat_cen)%dts(idt),JPRB)
          IF (dt1 .lt. 1.) THEN
             print *,'dt1 not good!!!!!!! ',dt1,stat_cen(istat_cen)%namestr
             stop
          ENDIF
          IF (obsstdeverror .lt. 1.e-7) THEN
             !sigma_o=stat_cen(istat_cen)%statistics%stdev_10days
             sigma_o=stat_cen(istat_cen)%statistics%stdev
          ELSE
             sigma_o=obsstdeverror
          ENDIF
          WRITE(unit=55,fmt='(A8,5X,F9.4,5X,F9.4,9X,F6.0,1X,F3.0,1X,F12.9,1X,F5.2,F9.4)', &
          iostat=iostatus) &
          stat_cen(istat_cen)%namestr, &
          stat_cen(istat_cen)%lat, &  !convert back to degree in output list
          stat_cen(istat_cen)%lon, &
          stat_cen(istat_cen)%alt, &
          dt1, &
          !stat_cen(istat_cen)%statistics%mean_10days, &
          !stat_cen(istat_cen)%statistics%mean, &
          0., &
          !stat_cen(istat_cen)%statistics%stdev_10days*1.e3, & ! junk
          stat_cen(istat_cen)%statistics%stdev*1.e3, & ! junk
          sigma_o                                             ! obs error stdev. in mm
          IF (iostatus /= 0) EXIT
          ns=ns+1
       ENDIF
    ENDDO

    readif: IF (iostatus > 0) THEN
              print *,'write_list_gpssol: An error has occurred while writing.'
              error_flag=error_writing_file
            ELSE
              print *,'write_list_gpssol: Finished writing. Found',ns,'records.'
            ENDIF readif
  ELSE
    print *,'read_write_gpssol: Could not open file ',trim(list_gpssol)
    error_flag=error_writing_file
  ENDIF
  CLOSE(unit=55)
  
END SUBROUTINE WRITE_LIST_GPSSOL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_LIST_GPSSOL(list_gpssol,nmaxs, &
    ns,statids,lats,lons,alts,dts,bias,obserrs)
  USE parkind1, ONLY:JPRB,JPIM
  implicit none

  ! I/O variables
  character(len=100), intent(in)   :: list_gpssol
  integer(KIND=JPIM), intent(in)   :: nmaxs
  integer(KIND=JPIM), intent(out)  :: ns
  character(len=8), dimension(nmaxs), intent(out):: statids
  real(KIND=JPRB),  dimension(nmaxs), intent(out):: lats,lons,alts,dts, &
                                                    bias,obserrs

  ! local variables
  integer(KIND=JPIM) :: iostatus
  character(len=100) :: string100
  character(len=91)  :: string90
  real(KIND=JPRB)    :: junk1

  print *,'read_list_gpssol: Opening file ',trim(list_gpssol)
  OPEN(unit=55,file=trim(list_gpssol),status='OLD',action='READ', &
       iostat=iostatus)
  ns=0
  IF (iostatus == 0) THEN
    readloop: DO WHILE (iostatus == 0)
      READ(unit=55,fmt='(A100)',iostat=iostatus) string100
      IF (iostatus /= 0) EXIT
      ns=ns+1
      statids(ns)=string100(1:8)
      string90=string100(10:len_trim(string100))
      READ(string90,*) lats(ns),lons(ns),alts(ns),dts(ns),bias(ns), &
                       junk1,obserrs(ns)
    ENDDO readloop

    readif: IF (iostatus > 0) THEN
              print *,'read_list_gpssol: An error has occurred while reading.'
            ELSE
              print *,'read_list_gpssol: Finished reading. Found',ns,'records.'
            ENDIF readif
  ELSE
    print *,'read_list_gpssol: Could not open file ',trim(list_gpssol)
  ENDIF
  CLOSE(unit=55)
  
END SUBROUTINE READ_LIST_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_OBSOUL_GPSSOL(fic_obsoul_gpssol,nmax, &
  n,statid,date,hhmmss,lat,lon,alt,dt,obsvalue,dat0,tim0)
  USE parkind1, ONLY:JPRB,JPIM
  implicit none

  ! I/O variables
  character(len=100), intent(in)   :: fic_obsoul_gpssol
  integer(KIND=JPIM), intent(in)   :: nmax
  integer(KIND=JPIM), intent(out)  :: n,dat0,tim0
  character(len=8), dimension(nmax), intent(out) :: statid
  integer(KIND=JPIM),dimension(nmax),intent(out) :: date,hhmmss
  real(KIND=JPRB),  dimension(nmax), intent(out) :: lat,lon,alt,dt, &
                                                    obsvalue

  ! local variables
  integer(KIND=JPIM) :: iostatus,ls,begstatid,endstatid,endhdr,i0,i1
  character(len=150) :: string100
  real(KIND=JPRB)    :: obserr0

  print *,'read_obsoul_gpssol: Opening file ',trim(fic_obsoul_gpssol)
  OPEN(unit=55,file=trim(fic_obsoul_gpssol),status='OLD',action='READ', &
       iostat=iostatus)
  n=0
  IF (iostatus == 0) THEN
    READ(unit=55,fmt=*,iostat=iostatus) dat0,tim0
    IF (iostatus == 0) THEN
      print *,'read_obsoul_gpssol: date',dat0,'time',tim0
    ELSE
      GOTO 100
    ENDIF
    readloop: DO WHILE (iostatus == 0)
      READ(unit=55,fmt='(A150)',iostat=iostatus) string100
      IF (iostatus /= 0) EXIT
      n=n+1
      ls=len_trim(string100)
      !print *,'string100 ','BEG'//string100(1:ls)//'END'
      IF (string100(1:11) .ne. ' 17 1 3110 ') stop
      string100=string100(12:ls)
      ls=len_trim(string100)
      begstatid=index(string100(1:ls),'''',BACK=.FALSE.)
      endstatid=index(string100(1:ls),'''',BACK=.TRUE.)
      IF (endstatid-begstatid+1 .ne. 10) stop
      endhdr=index(string100(1:ls),'1 11111 0 128')
      statid(n)=string100(begstatid+1:endstatid-1)
      READ(string100(1:begstatid-1),*) lat(n),lon(n)
      READ(string100(endstatid+1:endhdr-1),*) date(n),hhmmss(n),alt(n)
      string100=string100(endhdr+13:ls)
      READ(string100(1:len_trim(string100)),*) dt(n),obserr0,obsvalue(n)
      !print *,'statid/date/hhmmss/lat/lon/alt/dt/obsvalue'
      !print *,statid(n),date(n),hhmmss(n),lat(n),lon(n),alt(n),dt(n),obsvalue(n)
    ENDDO readloop

100 readif: IF (iostatus > 0) THEN
              print *,'read_obsoul_gpssol: An error has occurred while reading.'
            ELSE
              print *,'read_obsoul_gpssol: Finished reading. Found',n,'records.'
            ENDIF readif
  ELSE
    print *,'read_obsoul_gpssol: Could not open file ',trim(fic_obsoul_gpssol)
  ENDIF
  CLOSE(unit=55)

END SUBROUTINE READ_OBSOUL_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_TSLOT_GPSSOL(dat0,tim0,nmax,n,date,hhmmss,tslot,tslot_cent)
  USE parkind1, ONLY:JPRB,JPIM
  implicit none

  ! I/O variables
  integer(KIND=JPIM), intent(in)   :: dat0,tim0,nmax,n
  integer(KIND=JPIM), dimension(nmax), intent(in)  :: date
  integer(KIND=JPIM), dimension(nmax), intent(in)  :: hhmmss
  integer(KIND=JPIM), dimension(nmax), intent(out) :: tslot
  integer(KIND=JPIM), dimension(7),    intent(out) :: tslot_cent

  ! local variables
  integer(KIND=JPIM), dimension(8) :: tslots_00,tslots_06,tslots_12,tslots_18,&
                                      tslots
  integer(KIND=JPIM) :: i,j,obstime

  tslots_00=(/2100,2130,2230,2330,2430,2530,2630,2700/)
  tslots_06=(/300,330,430,530,630,730,830,900/)
  tslots_12=(/900,930,1030,1130,1230,1330,1430,1500/)
  tslots_18=(/1500,1530,1630,1730,1830,1930,2030,2100/)

  if (tim0 .eq. 0) tslots=100*tslots_00
  if (tim0 .eq. 6) tslots=100*tslots_06
  if (tim0 .eq.12) tslots=100*tslots_12
  if (tim0 .eq.18) tslots=100*tslots_18

  if (tim0 .eq. 0) tslot_cent=(/2115,2200,2300,0,100,200,245/)
  if (tim0 .eq. 6) tslot_cent=(/ 315, 400, 500,600,700,800,845/)
  if (tim0 .eq.12) tslot_cent=(/ 915,1000,1100,1200,1300,1400,1445/)
  if (tim0 .eq.18) tslot_cent=(/1515,1600,1700,1800,1900,2000,2045/)
  tslot_cent=tslot_cent*100

  !print *,'get_tslots_gpssol: tslots=',tslots

  tslot(1:n)=-1

  obsloop: DO i=1,n
    IF (date(i) > 20040101 .and. date(i) < 21000000 .and. &
        hhmmss(i) >= 0 .and. hhmmss(i) <= 235959) THEN
      IF (tim0 .eq. 0 .and. hhmmss(i) <= 30000) THEN
        obstime=hhmmss(i)+240000
      ELSE
        obstime=hhmmss(i)
      ENDIF
      IF (obstime >= tslots(1) .and. obstime <= tslots(2)) tslot(i)=1
      timloop: DO j=2,6
        IF (obstime > tslots(j) .and. obstime <= tslots(j+1)) tslot(i)=j
      ENDDO timloop
      IF (obstime > tslots(7) .and. obstime < tslots(8)) tslot(i)=7
    ELSE
      tslot(i)=-1
    ENDIF
  ENDDO obsloop
  !print *,'get_tslots_gpssol: tslot=',tslot(1:n)
  
END SUBROUTINE GET_TSLOT_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FILTER_GPSSOL(nmaxs,ns,statids,lats,lons,alts,dts, &
    nmax,n,statid,hhmmss,tslot,lat,lon,alt,dt,obsvalue, &
    nobs,mask)
  USE parkind1, ONLY:JPRB,JPIM
  implicit none

  ! I/O variables
  integer(KIND=JPIM), intent(in)   :: nmax,nmaxs,n,ns
  character(len=8),  dimension(nmaxs),intent(in) :: statids
  character(len=8),  dimension(nmax), intent(in) :: statid
  real(KIND=JPRB),   dimension(nmaxs),intent(in) :: lats,lons,alts,dts
  real(KIND=JPRB),   dimension(nmax), intent(in) :: lat,lon,alt,dt,obsvalue
  integer(KIND=JPIM),dimension(nmax), intent(in) :: hhmmss,tslot
  integer(KIND=JPIM), dimension(nmaxs,7), intent(out) :: nobs
  integer(KIND=JPIM), dimension(nmaxs,7,12), intent(out) :: mask

  ! local variables
  integer(KIND=JPIM) :: i,j,k,l,m,nobstot,hh0,mn0
  character(len=8)   :: statid0
  character(len=3)   :: cent0
  real(KIND=JPRB)    :: timecov

  nobs(:,:)=0
  k=0
  loopobs: DO i=1,n
    ! find whether the station is in the list
    l=-1
    DO j=1,ns
      IF (statid(i) .eq. statids(j)) l=j
    ENDDO
    IF (l == -1) GOTO 200
    statid0=statid(i)
    cent0=statid0(6:8)
    mn0=INT(hhmmss(i)/100)
    hh0=INT(hhmmss(i)/10000)
    mn0=mn0-hh0*100
    ! check that lat,lon,dt is the same as that specified in the list
    ! check that alt        is the same as that specified in the list (for 
    !         all centers EXCEPT ACR)
    ! check that minutes are either 00,15,30,or 45 for MET
    IF (abs(lat(i)-lats(l)) .gt. 0.1 .or. abs(lon(i)-lons(l)) .gt. 0.1 .or. &
      & (cent0 .ne. 'ACR' .and. abs(alt(i)-alts(l)) .gt. 0.5) .or. &
      & (cent0 .eq. 'MET' .and. (mn0.ne.0.and.mn0.ne.15.and.mn0.ne.30.and.mn0.ne.45)) .or. &
      & dt(i) .ne. dts(l)) GOTO 200
    ! check obsvalue physically possible and tslot defined
    IF (obsvalue(i) < 1.8 .or. obsvalue(i) > 3. .or. tslot(i) == -1) GOTO 200
    ! add observation into mask
    m=tslot(i)
    nobs(l,m)=nobs(l,m)+1
    mask(l,m,nobs(l,m))=i
200  k=k+1
  ENDDO loopobs

  print *,'filter_gpssol: Result of filtering operation:',ns,'stations'
  DO j=1,ns
    nobstot=SUM(nobs(j,:))
    timecov=FLOAT(nobstot)*dts(j)
    IF (SUM(nobs(j,:)) .ne. 0) THEN
      print '(A8,1X,7I3,1X,I3,A1)',statids(j),nobs(j,:),INT(timecov/360.*100.),'%'
    ELSE
      print '(A8,1X,A6)',statids(j),'no obs'
    ENDIF
  ENDDO

END SUBROUTINE FILTER_GPSSOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_OBSOUL_GPSSOL(fic_obsoul_gpssol_out, &
    nmax,nmaxs,ns,dat0,tim0,statids,date,hhmmss,tslot_cent, &
    lats,lons,alts,dts,obsvalue,bias,obserrs,nobs,mask)
  USE parkind1, ONLY:JPRB,JPIM
  implicit none
  
  ! I/O variables
  character(len=100), intent(in)   :: fic_obsoul_gpssol_out
  integer(KIND=JPIM), intent(in)   :: nmax,nmaxs,ns,dat0,tim0
  character(len=8),  dimension(nmaxs),intent(in) :: statids
  integer(KIND=JPIM),dimension(nmax), intent(in) :: date,hhmmss
  integer(KIND=JPIM),dimension(7),    intent(in) :: tslot_cent
  real(KIND=JPRB),   dimension(nmax), intent(in) :: obsvalue
  real(KIND=JPRB),   dimension(nmaxs),intent(in) :: lats,lons,alts,dts, &
                                                    bias,obserrs
  integer(KIND=JPIM), dimension(nmaxs,7), intent(in) :: nobs
  integer(KIND=JPIM), dimension(nmaxs,7,12), intent(in) :: mask

  ! local variables
  integer(KIND=JPIM) :: i,j,k,l,iostatus,nrecords,date_0
  real(KIND=JPRB)    :: obsv
  character(len=150) :: string100

  print *,'write_obsoul_gpssol: Opening file for writing ', &
    trim(fic_obsoul_gpssol_out)

  nrecords=0
  OPEN(unit=55,file=fic_obsoul_gpssol_out,action='WRITE',iostat=iostatus)
  IF (iostatus == 0) THEN
    WRITE(unit=55,fmt='(1X,I8,1X,I2)') dat0,tim0
    stationloop: DO i=1,ns
      timeloop: DO j=1,7
        IF ( (dts(i) >  5. .and. nobs(i,j) >  0) .or. &
           & (dts(i) == 5. .and. nobs(i,j) >  2) ) THEN
          obsv=0.
          obsloop: DO k=1,nobs(i,j)
            obsv=obsv+obsvalue(mask(i,j,k))
          ENDDO obsloop
          obsv=obsv/FLOAT(nobs(i,j))
          obsv=obsv-bias(i)
          l=mask(i,j,1)
          IF (tslot_cent(j) .eq. 0) THEN
            date_0=dat0
          ELSE
            date_0=date(l)
          ENDIF
          WRITE(string100,'(A10,1X,F10.5,2X,F10.5,2X,A1,A8,A1,2X,I8,1X,I6,1X,F6.1,1X,A13,1X,F5.1,1X,G11.5,1X,F7.5,1X,A10)') &
          '17 1 3110 ',lats(i),lons(i),'''',statids(i),'''',date_0,tslot_cent(j), &
          alts(i),'1 11111 0 128',dts(i)*float(nobs(i,j)),obserrs(i),obsv,'2147483647'
          WRITE(unit=55,fmt=*) trim(string100)
          nrecords=nrecords+1
        ENDIF
      ENDDO timeloop
    ENDDO stationloop
  ELSE
    print *,'write_obsoul_gpssol: Could not open file ', &
      trim(fic_obsoul_gpssol_out),' for writing.'
  ENDIF
  CLOSE(unit=55)
  print *,'write_obsoul_gpssol: Finished writing. Wrote',nrecords,'records.'

END SUBROUTINE WRITE_OBSOUL_GPSSOL

