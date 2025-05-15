! MODULE FOR GROUND-BASED SELECTION PACKAGE
!
! - defines structures
! - contains settings
!
! P. Poli original code 2006

MODULE GPSSOL_MOD

  USE parkind1, only: JPRB, JPIM

  implicit none

  integer(KIND=JPIM), parameter :: no_error = 0
  integer(KIND=JPIM), parameter :: error_reading_file = 1
  integer(KIND=JPIM), parameter :: error_allocating_array= 2
  integer(KIND=JPIM), parameter :: error_programming= 3
  integer(KIND=JPIM), parameter :: error_writing_file = 4

  integer(KIND=JPIM) :: nstat_cen, nstation, ncenter, &
                        nperiod, nperiod7, nstat_cen_ok, nstat_cen_active
  integer(KIND=JPIM), parameter :: nmaxs=5000,nmax=31000 !5000 31000
  integer(KIND=JPIM), parameter :: nstat_cen_max=2000 !2000 !nmaxs
  integer(KIND=JPIM), parameter :: nstation_max =900 !600 !nmaxs
  integer(KIND=JPIM), parameter :: ncenter_max  =15 ! old value 11, modified 2010/02/01
  integer(KIND=JPIM), parameter :: ndata_max    =50000
  integer(KIND=JPIM), parameter :: ndata_period_max=240 !120 ! 2 months 
  integer(KIND=JPIM), parameter :: nperiod_10days=80 !60 ! 1 months
  integer(KIND=JPIM), parameter :: nperiod_max  =240 !240  !120
  integer(KIND=JPIM), parameter :: n_cenexclude=1 ! set to 0 for investigation
  character(LEN=4),   parameter, dimension(n_cenexclude) :: cen_exclude=(/"TEST"/)
  !integer(KIND=JPIM), parameter :: n_cenexclude=4 ! set to 0 for investigation
  !character(LEN=4),   parameter, dimension(n_cenexclude) :: cen_exclude=(/"KNM1","BKG_","IEEC","SGN1"/)
  !real(KIND=JPRB),    parameter :: biasmax =200.E-3 ! set to 99 to retain all stations, otherwise set in meters
  !real(KIND=JPRB),    parameter :: biasmax =40.E-3 ! set to 99 to retain all stations, otherwise set in meters
  real(KIND=JPRB),    parameter :: biasmax=40.E-3 ! set to 99 to retain all stations, otherwise set in meters
  !real(KIND=JPRB),    parameter :: stdevmax=200.E-3 ! set to 99 to retain all stations, otherwise set in meters
  !real(KIND=JPRB),    parameter :: stdevmax=40.E-3 ! set to 99 to retain all stations, otherwise set in meters
  real(KIND=JPRB),    parameter :: stdevmax=30.E-3 ! set to 99 to retain all stations, otherwise set in meters
  real(KIND=JPRB),    parameter :: obsstdeverror=0.  ! 0. => use stat_cen(istat_cen)%statistics%stdev_10days,
                                                     ! otherwise set in meters
  logical,            parameter :: flag_separate_59mins = .true.
  logical,            parameter :: flag_exclude_59mins = .true.
  ! note: to activate flag_exclude_59mins you need to activate also flag_separate_59mins
  logical,            parameter :: flag_retain_only_select_mins = .false.
  integer(KIND=JPIM), parameter :: n_retains=1
  character(LEN=2),   parameter, dimension(n_retains) :: retain_mins=(/"00"/) ! force to retain obs at 00 for 'SGN_'
  character(LEN=4),   parameter, dimension(n_retains) :: cen_mins=(/"SGN_"/)
  logical,            parameter :: flag_retain_only_central_time_obs = .true.
  ! flag_retain_only_central_time_obs only works for passive obs
  ! for 3-hour windows:
  integer(KIND=JPIM), parameter, dimension(9) :: &
                   screening_window_bounds=(/223000,013000,043000,073000,103000,133000,163000,193000,223000/)
  integer(KIND=JPIM), parameter, dimension(8) :: &
              screening_window_centraltime=(/000000,030000,060000,090000,120000,150000,180000,210000/)
  integer(KIND=JPIM), parameter :: n_screening_windows=8
  integer(KIND=JPIM), parameter :: nslot_max=1 ! set to 1 for 3dvar or 7 for 4dvar
  real(KIND=JPRB),    parameter :: min_distance =15. ! km old value 10km for aladin/arome
  integer(KIND=JPIM), parameter :: max_dt=99 ! set to 99 to retain all stations, 16 to retain only data on 15 min
  !character(len=100) :: fic_list_files='/home/poli/run_odb/gpssol/WORK/list_files'
  character(len=100) :: fic_list_files='list_files'
  character(len=100) :: fic_list_gpssol='list_gpssol'
  character(len=5),   parameter :: density_method='histo' ! choose 'histo' [histogramme] or 'kde' [kernel density estimator]
  integer(KIND=JPIM), parameter :: ndata_min_in_class_gauss=10, &
                                   nclass_gauss=70, &
                                   nconf_gauss=3, &  ! 1,2, or 3 . Level of confidence
                                   kde_P=0           ! 0,1, or 2 . Only matters for 'kde'
  real(KIND=JPRB) :: lonmin=2., lonmax=39., latmin=33., latmax=55.
  real(KIND=JPRB),    parameter :: diff_alt_max=600., alt_max=2500.
  real(KIND=JPRB),    parameter :: min_time_coverage=40. ! %
  integer(KIND=JPIM), parameter :: ndt_max=10
  character(len=*),   parameter :: ficoutput_ext='_aust'  
  integer(KIND=JPIM), parameter :: nstat_cen_active_max=900 !250

  TYPE DATA_TYPE
    real(KIND=JPRB)    :: fg_depar !, ZWD_calc !obsvalue, ZHD_calc
  END TYPE

  TYPE DATA_ACTIVE_TYPE
    integer(KIND=JPIM) :: tslot
    real(KIND=JPRB)    :: fg_depar, ZWD_calc, an_depar !ZHD_calc, obsvalue
  END TYPE

  TYPE STATISTICS_PERIOD_TYPE
    integer(KIND=JPIM) :: ndata, abs_time_diff
    integer(KIND=JPIM), dimension(1:ndata_period_max) :: mask
    real(KIND=JPRB)    :: mean, stdev
    integer(KIND=JPIM), dimension(:), allocatable :: data_present7
    real(KIND=JPRB)    :: time_coverage7
  END TYPE
  
  TYPE STATISTICS_TYPE
    integer(KIND=JPIM) :: ndata
    integer(KIND=JPIM) :: ndata_10days
    real(KIND=JPRB)    :: mean, stdev
    real(KIND=JPRB)    :: mean_10days, stdev_10days
    logical            :: gauss
    integer(KIND=JPIM), dimension(:), allocatable :: data_present7
    integer(KIND=JPIM), dimension(:), allocatable :: data_present
    real(KIND=JPRB)    :: time_coverage7
    real(KIND=JPRB)    :: time_coverage7_10days
    real(KIND=JPRB)    :: time_coverage
    real(KIND=JPRB)    :: time_coverage_10days
  END TYPE

  TYPE STATISTICS_ACTIVE_TYPE
    integer(KIND=JPIM),dimension(0:nslot_max) :: ndata, ndata_rejected
    real(KIND=JPRB),dimension(0:nslot_max) :: mean_hxb_hxa,mean_y0_hxb,mean_y0_hxa, &
                                   stdev_hxb_hxa,stdev_y0_hxb,stdev_y0_hxa, &
                                      sigb, siga, sigo
  END TYPE

  TYPE STAT_CEN_TYPE
    character(len=8)    :: namestr
    real(KIND=JPRB)     :: lat, lon, alt, diff_alt, alt_mod
    integer(KIND=JPIM)  :: ndt, dts(1:ndt_max), ndts(1:ndt_max)
    logical             :: cst_lat, cst_lon, cst_alt, cst_diff_alt
    TYPE(DATA_TYPE), dimension(1:ndata_max) :: data
    TYPE(STATISTICS_TYPE) :: statistics
    TYPE(STATISTICS_PERIOD_TYPE), dimension(1:nperiod_max) :: tperiod_statistics
    logical             :: ok
    real(KIND=JPRB)     :: mean_min,mean_max,mean_mean,mean_stdev, &
                           stdev_min,stdev_max,stdev_mean,stdev_stdev
    integer(KIND=JPIM)  :: nperiod_ok
  END TYPE

  TYPE STAT_CEN_ACTIVE_TYPE
    character(len=8)    :: namestr
    TYPE(DATA_ACTIVE_TYPE), dimension(1:ndata_max) :: data
    TYPE(STATISTICS_ACTIVE_TYPE) :: statistics
  END TYPE

  TYPE CENTER_TYPE
    character(len=4)   :: namestr
    TYPE(STAT_CEN_TYPE), dimension(:), pointer :: stat_cen
    TYPE(STATISTICS_TYPE) :: statistics
  END TYPE

  TYPE STATION_TYPE
    character(len=4)   :: namestr
    TYPE(STAT_CEN_TYPE), dimension(:), pointer :: stat_cen
    TYPE(STATISTICS_TYPE) :: statistics
  END TYPE

  public :: find_statid, find_station, find_center

CONTAINS

FUNCTION FIND_STATID(nexisting,nids,ids,id1) result(i_id)
  ! I/O variables
  integer(KIND=JPIM),                     intent(in) :: nexisting, nids
  TYPE(STAT_CEN_TYPE), dimension(1:nids), intent(in) :: ids
  character(len=8),                       intent(in) :: id1
  integer(KIND=JPIM) :: i_id

  ! local variables
  logical            :: ll_found
  integer(KIND=JPIM) :: imax

  i_id=-1
  IF (nexisting .eq. 0) THEN
    return
  ENDIF
  ll_found=.FALSE.
  imax=minval((/nexisting,nids/))
  i_id=0
  DO WHILE (i_id .lt. imax .and. (.not.(ll_found)))
    i_id = i_id+1
    IF (ids(i_id)%namestr .eq. id1) ll_found=.TRUE.
  ENDDO
  IF (.not.(ll_found)) i_id=-1
END FUNCTION FIND_STATID

FUNCTION FIND_CHAIN(nexisting,lenstr,nids,ids,id1,mask) result(i_id)
  ! I/O variables
  integer(KIND=JPIM),                     intent(in) :: nexisting, lenstr, nids
  character(len=lenstr),dimension(1:nids),intent(in) :: ids
  character(len=lenstr),                  intent(in) :: id1
  logical, dimension(1:nids), optional,   intent(in) :: mask
  integer(KIND=JPIM), dimension(0:15) :: i_id

  ! local variables
  integer(KIND=JPIM) :: imax, icount, i0
  logical :: ll_found

  i_id(:)=-1
  IF (nexisting .eq. 0) THEN
    return
  ENDIF
  imax=minval((/nexisting,nids/))
  i0=0
  DO icount=1,imax
    ll_found = ids(icount) .eq. id1
    IF (present(mask)) ll_found=ll_found .and. mask(icount)
    IF (ll_found) THEN
       i0=i0+1
       i_id(i0)=icount
    ENDIF
  ENDDO
  i_id(0)=i0
END FUNCTION FIND_CHAIN

FUNCTION FIND_STATION(nexisting,nids,ids,id1) result(i_id)
  ! I/O variables
  integer(KIND=JPIM),                     intent(in) :: nexisting, nids
  TYPE(STATION_TYPE),  dimension(1:nids), intent(in) :: ids
  character(len=4),                       intent(in) :: id1
  integer(KIND=JPIM) :: i_id

  ! local variables
  logical            :: ll_found
  integer(KIND=JPIM) :: imax

  i_id=-1
  IF (nexisting .eq. 0) THEN
    return
  ENDIF
  ll_found=.FALSE.
  imax=minval((/nexisting,nids/))
  i_id=0
  DO WHILE (i_id .lt. imax .and. (.not.(ll_found)))
    i_id = i_id+1
    IF (ids(i_id)%namestr .eq. id1) ll_found=.TRUE.
  ENDDO
  IF (.not.(ll_found)) i_id=-1
END FUNCTION FIND_STATION

FUNCTION FIND_CENTER (nexisting,nids,ids,id1) result(i_id)
  ! I/O variables
  integer(KIND=JPIM),                     intent(in) :: nexisting, nids
  TYPE(CENTER_TYPE),   dimension(1:nids), intent(in) :: ids
  character(len=4),                       intent(in) :: id1
  integer(KIND=JPIM) :: i_id

  ! local variables
  logical            :: ll_found
  integer(KIND=JPIM) :: imax

  i_id=-1
  IF (nexisting .eq. 0) THEN
    return
  ENDIF
  ll_found=.FALSE.
  imax=minval((/nexisting,nids/))
  i_id=0
  DO WHILE (i_id .lt. imax .and. (.not.(ll_found)))
    i_id = i_id+1
    IF (ids(i_id)%namestr .eq. id1) ll_found=.TRUE.
  ENDDO
  IF (.not.(ll_found)) i_id=-1
END FUNCTION FIND_CENTER 

FUNCTION STATS(ndata,tab) result(tab2)

  ! I/O variables
  integer(KIND=JPIM), intent(in) :: ndata
  real(KIND=JPRB), dimension(1:ndata), intent(in) :: tab
  real(KIND=JPRB), dimension(2) :: tab2

  ! local variables
  integer(KIND=JPIM) :: i
  real(KIND=JPRB) :: sum1, sum2, mean, stdev

  sum1=0.
  DO i=1,ndata
     sum1=sum1+tab(i)
  ENDDO
  mean=sum1/REAL(ndata,JPRB)

  sum2=0.
  DO i=1,ndata
     sum2=sum2+(tab(i)-mean)**2.
!     sum2=sum2+tab(i)**2.
  ENDDO
  stdev=sqrt(sum2/REAL(ndata-1,JPRB))
!  stdev=sqrt(sum2/REAL(ndata,JPRB)-mean**2.)
  tab2=(/mean,stdev/)
  
END FUNCTION STATS

SUBROUTINE PRINT_STATISTICS(stat_cen,fic_stats,okonly)

  ! I/O variables
  TYPE(STAT_CEN_TYPE), dimension(1:nstat_cen_max), &
                     intent(in) :: stat_cen
  character(len=*),  intent(in) :: fic_stats
  logical, optional, intent(in) :: okonly

  ! local variables
  integer(KIND=JPIM) :: istat_cen, iperiod, idt
  character(len=3) :: strfmt, strdtadd, str3, txt
  character(len=150) :: strformat, strdt
  logical :: ll_okonly

  ll_okonly=.FALSE.
  IF (present(okonly)) ll_okonly=okonly
 
  print '(A)', 'print_statistics: Opening file ',trim(fic_stats)
  open (54,file=fic_stats)
  write (54,'(A)') &
  & "0000-000 Nb Tsltcv|bias: MIN AVG STD MAX |std: MIN TOT STD MAX |GaussCstAltOrog LaLo DT ALT"
  DO istat_cen=1,nstat_cen
     IF (stat_cen(istat_cen)%statistics%gauss) THEN
        str3(1:1)='T'
     ELSE
        str3(1:1)='F'
     ENDIF
     IF (stat_cen(istat_cen)%cst_alt) THEN
        str3(2:2)='T'
     ELSE
        str3(2:2)='F'
     ENDIF
     IF (ABS(stat_cen(istat_cen)%alt_mod - stat_cen(istat_cen)%alt ) .le. diff_alt_max .and. &
       & stat_cen(istat_cen)%alt .le. alt_max) THEN
        str3(3:3)='T'
     ELSE
        str3(3:3)='F'
     ENDIF
     IF (((ll_okonly .and. stat_cen(istat_cen)%ok) .or. &
        & .not.(ll_okonly)) .and. &
        (stat_cen(istat_cen)%statistics%ndata .gt. 2)) THEN
     strdt(:)=' '
     DO idt=1,stat_cen(istat_cen)%ndt
        strdtadd(:)=' '
        IF (stat_cen(istat_cen)%dts(idt) .lt. 10) THEN
           write (strdtadd,'(I2)') stat_cen(istat_cen)%dts(idt)
        ELSE
           write (strdtadd,'(I3)') stat_cen(istat_cen)%dts(idt)
        ENDIF
        strdt=trim(strdt)//trim(strdtadd)
     ENDDO
     IF (stat_cen(istat_cen)%ndt .gt. 1) THEN
        IF (stat_cen(istat_cen)%ndt .gt. 2) THEN
           txt="++"
           print *,'found more than 2 dts :', &
           stat_cen(istat_cen)%namestr, stat_cen(istat_cen)%ndt, &
           stat_cen(istat_cen)%dts(1:stat_cen(istat_cen)%ndt)
           print *,'ERROR!!!!!'
!          stop
        ELSE
           txt=""
        ENDIF
        write (strdtadd,'(I2,"%")') INT(REAL(stat_cen(istat_cen)%ndts(2))/ &
         & (REAL(stat_cen(istat_cen)%ndts(1))+ &
         &  REAL(stat_cen(istat_cen)%ndts(2)))*100.,JPIM)
        strdt=trim(strdt)//'->'//trim(strdtadd)//trim(txt)
     ENDIF
     write (strfmt,'("A",I2.2)') len_trim(strdt)
     strformat='(A8,I5,1X,I3,"%|", &
            & F5.1,F5.1,"+/-",F4.1,F5.1,"|", &
            & F4.1,1X,F4.1,"+/-",F4.1,1X,F4.1,"|", &
            & A3,F5.2,F6.2,'//trim(strfmt)//',1X,F6.1)'
     write (54,strformat) &
      stat_cen(istat_cen)%namestr, &
      stat_cen(istat_cen)%statistics%ndata, &
      INT(stat_cen(istat_cen)%statistics%time_coverage7,JPIM), &
!      INT(minval( &
!      (/stat_cen(istat_cen)%statistics%time_coverage7_10days, &
!        stat_cen(istat_cen)%statistics%time_coverage7       , &
!        stat_cen(istat_cen)%statistics%time_coverage_10days , &
!        stat_cen(istat_cen)%statistics%time_coverage        /)), &
!        JPIM), &
!      stat_cen(istat_cen)%nperiod_ok, &
      stat_cen(istat_cen)%mean_min  *1.e3, &
      stat_cen(istat_cen)%statistics%mean*1.e3, &
!      stat_cen(istat_cen)%mean_mean *1.e3, &
      stat_cen(istat_cen)%mean_stdev*1.e3, &
      stat_cen(istat_cen)%mean_max  *1.e3, &
      stat_cen(istat_cen)%stdev_min  *1.e3, &
      stat_cen(istat_cen)%statistics%stdev*1.e3, &
!      stat_cen(istat_cen)%stdev_mean*1.e3, &
      stat_cen(istat_cen)%stdev_stdev*1.e3, &
      stat_cen(istat_cen)%stdev_max  *1.e3, &
      str3, &
      stat_cen(istat_cen)%lat, &
      stat_cen(istat_cen)%lon, &
      strdt, &
      stat_cen(istat_cen)%alt
!      stat_cen(istat_cen)%statistics%gauss, &
!      stat_cen(istat_cen)%cst_lat.and.&
!       & stat_cen(istat_cen)%cst_lon,
!      stat_cen(istat_cen)%cst_alt
!      print *,stat_cen(istat_cen)%namestr, &
!              stat_cen(istat_cen)%dts(1:stat_cen(istat_cen)%ndt), &
!              stat_cen(istat_cen)%ndts(1:stat_cen(istat_cen)%ndt)
      ENDIF
  ENDDO

  IF (.not.(ll_okonly)) THEN
  DO istat_cen=1,nstat_cen
     DO iperiod=1,nperiod
        IF (stat_cen(istat_cen)%tperiod_statistics(iperiod)%ndata .gt. 2) &
        write (54,'(I3,1X,A8,1X,I6,1X,I3,"%",1X,F5.1,1X,F4.1)') &
         iperiod, &
         stat_cen(istat_cen)%namestr, &
         stat_cen(istat_cen)%tperiod_statistics(iperiod)%ndata, &
         INT(stat_cen(istat_cen)%tperiod_statistics(iperiod)%time_coverage7,JPIM), &
         stat_cen(istat_cen)%tperiod_statistics(iperiod)%mean*1.e3, &
         stat_cen(istat_cen)%tperiod_statistics(iperiod)%stdev*1.e3
     ENDDO
  ENDDO
  ENDIF

  close(54)
  
END SUBROUTINE PRINT_STATISTICS

SUBROUTINE PRINT_STATISTICS_ACTIVE(stat_cen_active,statistics_active,fic_stats)

  ! I/O variables
  TYPE(STAT_CEN_ACTIVE_TYPE), dimension(1:nstat_cen_active_max), &
                     intent(in) :: stat_cen_active
  TYPE(STATISTICS_ACTIVE_TYPE), intent(in) :: statistics_active
  character(len=*),  intent(in) :: fic_stats

  ! local variables
  integer(KIND=JPIM) :: istat_cen, islot, ndata, ndata_rejected
  character(len=1000) :: str_slot_add0,str_slot_add,str

  print '(A)', 'print_statistics_active: Opening file '
  print '(A)', trim(fic_stats)
  open (54,file=fic_stats)
  write (54,'(A,I3,A)') &
  & "0000-000 ",nstat_cen_active, &
  & "   Nb %rej   hxb_hxa, y0_hxb, y0_hxa :mean stdev|  sigo  sigb  siga"
  DO istat_cen=1,nstat_cen_active
     DO islot=0,0
        str=''
        ndata=stat_cen_active(istat_cen)%statistics%ndata(islot)
        ndata_rejected=stat_cen_active(istat_cen)%statistics%ndata_rejected(islot)
        write (str,'(A8,1X,I1,": ",I5,1X,I3,"%")') &
          stat_cen_active(istat_cen)%namestr, islot, ndata, &
          INT(REAL(ndata_rejected,JPRB)/(REAL(ndata,JPRB)+REAL(ndata_rejected,JPRB))*100.)
        str_slot_add(:)=' '
        IF (max(ndata,ndata_rejected) .gt. 0) THEN
           write (str_slot_add,'(3F6.1)') &
           stat_cen_active(istat_cen)%statistics%sigo(islot)*1.e3, &
           stat_cen_active(istat_cen)%statistics%sigb(islot)*1.e3, &
           stat_cen_active(istat_cen)%statistics%siga(islot)*1.e3
!           IF (islot.eq.0) THEN
              write (str_slot_add0,'(6F6.1)') &
              stat_cen_active(istat_cen)%statistics%mean_hxb_hxa(islot)*1.e3, &
              stat_cen_active(istat_cen)%statistics%mean_y0_hxb(islot)*1.e3, &
              stat_cen_active(istat_cen)%statistics%mean_y0_hxa(islot)*1.e3, &
              stat_cen_active(istat_cen)%statistics%stdev_hxb_hxa(islot)*1.e3, &
              stat_cen_active(istat_cen)%statistics%stdev_y0_hxb(islot)*1.e3, &
              stat_cen_active(istat_cen)%statistics%stdev_y0_hxa(islot)*1.e3
              str_slot_add=str_slot_add0(1:36)//' |'//str_slot_add(1:18)
!           ENDIF
           str=trim(str)//' '//trim(str_slot_add)
           write (54,'(A)') trim(str)
        ENDIF
     ENDDO
  ENDDO

  DO islot=0,nslot_max
     str=''
     ndata=statistics_active%ndata(islot)
     ndata_rejected=statistics_active%ndata_rejected(islot)
     write (str,'(A8,1X,I1,": ",I5,1X,I3,"%")') &
       'STAT-CEN', islot, ndata, &
       INT(REAL(ndata_rejected,JPRB)/(REAL(ndata,JPRB)+REAL(ndata_rejected,JPRB))*100.)
     str_slot_add(:)=' '
     IF (ndata .gt. 2) THEN
        write (str_slot_add,'(3F6.1)') &
        statistics_active%sigo(islot)*1.e3, &
        statistics_active%sigb(islot)*1.e3, &
        statistics_active%siga(islot)*1.e3
        write (str_slot_add0,'(6F6.1)') &
        statistics_active%mean_hxb_hxa(islot)*1.e3, &
        statistics_active%mean_y0_hxb(islot)*1.e3, &
        statistics_active%mean_y0_hxa(islot)*1.e3, &
        statistics_active%stdev_hxb_hxa(islot)*1.e3, &
        statistics_active%stdev_y0_hxb(islot)*1.e3, &
        statistics_active%stdev_y0_hxa(islot)*1.e3
        str_slot_add=str_slot_add0(1:36)//' |'//str_slot_add(1:18)
        str=trim(str)//' '//trim(str_slot_add)
        write (54,'(A)') trim(str)
     ENDIF
  ENDDO

  close(54)
  
END SUBROUTINE PRINT_STATISTICS_ACTIVE

SUBROUTINE GET_ABS_TIMEDIFF(hhmmss,abs_time_diff)
  ! I/O variables
  integer(KIND=JPIM), intent(in)  :: hhmmss
  integer(KIND=JPIM), intent(out) :: abs_time_diff

  ! local variables
  integer(KIND=JPIM) :: iwindow,jwindow,obstime,centraltime
  integer(KIND=JPIM) :: obstime_hh,obstime_mm !,obstime_ss
  integer(KIND=JPIM) :: centraltime_hh,centraltime_mm !,centraltime_ss

  jwindow=-1
  IF (hhmmss>=screening_window_bounds(1) .or. hhmmss < screening_window_bounds(2)) THEN
     jwindow=1
  ELSE
     DO iwindow=2,n_screening_windows
        IF (hhmmss >= screening_window_bounds(iwindow) .and. &
          & hhmmss <  screening_window_bounds(iwindow+1)) THEN
           jwindow=iwindow
        ENDIF
     ENDDO
  ENDIF
  IF (jwindow==-1) THEN
     print *,'ERROR!!! could not find time difference for obs ', hhmmss
     STOP
  ENDIF

  IF (jwindow == 1) THEN
     obstime=hhmmss+240000
     centraltime=screening_window_centraltime(jwindow)+240000
  ELSE
     obstime=hhmmss
     centraltime=screening_window_centraltime(jwindow)
  ENDIF

  obstime_hh=INT(obstime/10000)
  obstime_mm=INT(obstime-obstime_hh*10000)/100
!  obstime_ss=obstime-obstime_hh*10000-obstime_mm*100
  centraltime_hh=INT(centraltime/10000)
  centraltime_mm=INT(centraltime-centraltime_hh*10000)/100
!  centraltime_ss=centraltime-centraltime_hh*10000-centraltime_mm*100

  abs_time_diff=abs((obstime_hh*60+obstime_mm)-(centraltime_hh*60+centraltime_mm))
!  print *,'hhmmss,centraltime,jwindow,abs_time_diff',hhmmss,obstime_hh,obstime_mm,centraltime,centraltime_hh,centraltime_mm,jwindow,abs_time_diff

END SUBROUTINE GET_ABS_TIMEDIFF

END MODULE GPSSOL_MOD
     
