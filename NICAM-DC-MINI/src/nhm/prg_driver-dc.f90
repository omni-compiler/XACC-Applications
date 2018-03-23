!-------------------------------------------------------------------------------
!
!+  Program driver
!
!-------------------------------------------------------------------------------
program prg_driver
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This program is a driver of non-hydrostatic model based on an 
  !       icosahedral grid system.
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !      0.03      04-05-31   Change by addtion of mod[onestep].
  !                05-12-01   M.Satoh add history_setup
  !                05-12-19   S.Iga moved output_timeinfo after output_all
  !                06-04-18   T.Mitsui add sfc_restart
  !                06-04-21   H.Tomita:  remove output_timeinfo due to 
  !                                      computational efficeincy.
  !                                      Instead, this process is move to 
  !                                      mod[mod_output].
  !                06-08-07   W.Yanase add history_vars
  !                06-09-27   S.Iga add history_vars_cfmip
  !                07-03-23   Y.Niwa add ndg_setup, ndg_do, FLAG_NUDGING
  !                07-06-27   Y.Niwa add history_vars_setup
  !                07-07-24   K.Suzuki: implementing SPRINTARS aerosol model
  !                07-08-06   Y.Niwa: add history_obs, history_vars_obs
  !                07-11-07   T.Mitsui: add option to omit output_all
  !                08-03-10   T.Mitsui: add intermediate output of restart file
  !                08-05-24   T.Mitsui: trivial fix
  !                08-09-09   Y.Niwa : modfied for nudging
  !                09-01-23   H.Tomita: a) abolish mod_output, mod_extdata
  !                                     mod_history_vars_cfmip, mod_o3var.
  !                                     b) introduce mod_extdata2.
  !                09-04-14   T.Mitsui: arrange initialization of aerosols
  !                09-07-10   H.Tomita: Add the module [mod_embudget].
  !                09-08-05   S.Iga: remove latlon_setup (suggested by T.Mitsui)
  !                09-08-05   T.Mitsui: add conditioning by ADM_myprc_is_run
  !                                     to keep out extra-processes from main routines.
  !                09-08-18   T.Mitsui: change of 09-08-05 is not enough.
  !                10-03-08   C.Kodama: Modify for overwrite_restart option
  !                10-04-30   M.Satoh: move diagvar_setup
  !                11-09-03   H.Yashiro : New I/O
  !                11-11-28   Y.Yamada : merge Terai-san timer
  !                12-06-07   T.Seiki  : Application to Multi-job System
  !                12-10-12   R.Yoshida  : Modify for Dynamical Core test
  !                12-10-22   R.Yoshida  : add papi instructions
  !      ----------------------------------------------------------------------- 
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_MULTI_PRC,      &
     ADM_LOG_FID,        &
     ADM_prc_me,         &
     ADM_prc_run_master, &
     ADM_proc_init,      &
     ADM_proc_stop,      &
     ADM_setup, &
     ADM_prc_tab,ADM_rgn_vnum,ADM_IopJop  ! openacc
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup, &
     sendlist,sendlist_pl, &  ! openacc
     sendinfo,sendinfo_pl, &
     recvlist,recvlist_pl, &
     recvinfo,recvinfo_pl, &
     recvlist_r2r,sendlist_r2r, &
     recvlist_r2p,sendlist_r2p, &
     recvlist_p2r,sendlist_p2r, &
     recvlist_sgp,sendlist_sgp, &
     copyinfo_r2r,copyinfo_r2p,copyinfo_p2r,copyinfo_sgp, &
     nsmax,nsmax_pl,ncmax_r2r,ncmax_r2p,ncmax_p2r,nrmax,nrmax_pl,ncmax_sgp
  use mod_cnst, only: &
     CNST_setup
  use mod_calendar, only: &
     calendar_setup
  use mod_time, only: &
     TIME_setup,          &
     TIME_report,         &
     TIME_advance,        &
     TIME_LSTEP_MAX,      &
     cstep => TIME_CSTEP, &
     ctime => TIME_CTIME, &
     dtime => TIME_DTL
  use mod_grd, only: &
     GRD_setup, &
     GRD_afac,GRD_bfac,GRD_cfac,GRD_dfac,GRD_e,GRD_rdgz,GRD_rdgzh,GRD_x,GRD_xt,GRD_vz,GRD_zs  ! openacc
  use mod_gmtr, only: &
     GMTR_setup, &
     GMTR_T_var,GMTR_P_var,GMTR_A_var,GMTR_lat    ! openacc
  use mod_oprt, only: &
     OPRT_setup, &
     cdiv,cgrad,clap, cmdif_P,cmdif_T,cmdif_AH,cmdif_AT  ! openacc
  use mod_vmtr, only: &
     VMTR_setup, &
     VMTR_C2Wfact,VMTR_RGAMH,VMTR_RGAM,VMTR_RGSH, &  ! openaccc
     VMTR_RGSGAM2,VMTR_GSGAM2H,VMTR_GSGAM2, &
     VMTR_RGAM2H,VMTR_RGAM2,VMTR_RGSGAM2H,VMTR_GAM2H, &
     VMTR_GSGAMH,VMTR_GZXH,VMTR_GZYH,VMTR_GZZH,VMTR_PHI
  use mod_runconf, only: &
     runconf_setup, &
     FLAG_NUDGING, &
     CVW  ! openacc
  use mod_prgvar, only: &
     prgvar_setup,            &
     restart_input_basename,  &
     restart_output_basename, &
     restart_input,           &
     restart_output, &
     PRG_var, &    ! openacc
     PRG_var1, &   !
     DIAG_var      !
  use mod_diagvar, only: &
       diagvar_setup, &
       diagvar_restart_output
  use mod_sfcvar, only: &
       sfcvar_setup, &
       sfcvar, KSTR    ! openacc
  use mod_bsstate, only: &
       bsstate_setup, &
       pre_bs,tem_bs,rho_bs,phi  ! openacc
  use mod_bndcnd, only   :   &
       bndcnd_setup
  use mod_numfilter, only  : &
       numfilter_setup, &
       divdamp_coef,Kh_coef,Kh_coef_lap1  ! openacc
  use mod_forcing_driver, only : &
       forcing_init, &
       forcing
  use mod_ndg, only: &
       ndg_setup
  use mod_dynstep, only : &
       dynstep
  use mod_history, only: &
       history_setup, &
       history_out,   &
       HIST_output_step0, &
       v_save, &
       ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2    ! openacc
  use mod_history_vars, only: &
       history_vars_setup, &
       history_vars
  use mod_embudget, only: &
       embudget_setup, &
       embudget_monitor
  use mod_vi, only : &
       vi_alloc, &
       A2_o,CooCip,CooCim,D2,Mc,Mu,Ml  ! openacc
  use mod_trcadv_thuburn, only: &
       OPRT_divergence2_alloc, &  ! openacc
       local_t_var
  implicit none

  character(len=14) :: cdate

  integer :: n
  !-----------------------------------------------------------------------------

  call ADM_proc_init(ADM_MULTI_PRC)

  !---< admin module setup >---
  call ADM_setup('nhm_driver.cnf')

  !#############################################################################

  write(ADM_LOG_FID,*) '##### start  setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  setup     #####'
  endif

  call DEBUG_rapstart('Total')
  call DEBUG_rapstart('Setup ALL')

  !---< cnst module setup >---
  call CNST_setup

  !---< calendar module setup >---
  call calendar_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< grid module setup >---
  call GRD_setup

  !---< geometrics module setup >---
  call GMTR_setup

  !---< operator module setup >---
  call OPRT_setup

  !---< vertical metrics module setup >---
  call VMTR_setup

  !---< time module setup >---
  call TIME_setup


  !---< nhm_runconf module setup >---
  call runconf_setup

  !---< prognostic variable module setup >---
  call prgvar_setup
  call restart_input( restart_input_basename )

  !---< diagnostic variable module setup >---
  call diagvar_setup

  !---< surface variable module setup >---
  call sfcvar_setup


  !---< boundary condition module setup >---
  call bndcnd_setup

  !---< basic state module setup >---
  call bsstate_setup

  !---< numerical filter module setup >---
  call numfilter_setup

  !---< forcing module setup >---
  call forcing_init

  !---< energy&mass budget module setup >---
  call embudget_setup

  !---< history module setup >---
  call history_setup

  !---< history variable module setup >---
  call history_vars_setup

  !---< nudging module setup >---
  if( FLAG_NUDGING ) call ndg_setup( ctime, dtime )

  write(ADM_LOG_FID,*) '##### finish setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish setup     #####'
  endif

  call DEBUG_rapend('Setup ALL')

  !#############################################################################
#ifdef _FIPP_
  call fipp_start()
#endif
  call DEBUG_rapstart('Main ALL')

  write(ADM_LOG_FID,*) '##### start  main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  main loop #####'
  endif

  call TIME_report

  call vi_alloc()
  call OPRT_divergence2_alloc()

  !---- !$acc& pcopy(v_save) &

  !$acc data &
  !$acc& pcopyin(VMTR_C2Wfact,VMTR_RGAMH,VMTR_RGAM,VMTR_RGSH) &
  !$acc& pcopyin(VMTR_RGSGAM2,VMTR_GSGAM2H,VMTR_GSGAM2) &
  !$acc& pcopyin(VMTR_RGAM2H,VMTR_RGAM2,VMTR_RGSGAM2H,VMTR_GAM2H) &
  !$acc& pcopyin(VMTR_GSGAMH, VMTR_GZXH,VMTR_GZYH,VMTR_GZZH,VMTR_PHI) &
  !$acc& pcopyin(GRD_afac,GRD_bfac,GRD_cfac,GRD_dfac,GRD_e,GRD_rdgz,GRD_rdgzh,GRD_x,GRD_xt,GRD_vz,GRD_zs) &
  !$acc& pcopyin(GMTR_T_var,GMTR_P_var,GMTR_A_var,GMTR_lat) &
  !$acc& pcopyin(ADM_prc_tab,ADM_rgn_vnum,ADM_IopJop) &
  !$acc& pcopyin(cdiv,cgrad,clap, cmdif_P,cmdif_T,cmdif_AH,cmdif_AT) &
  !$acc& pcopyin(pre_bs,tem_bs,rho_bs,phi) &
  !$acc& pcopyin(divdamp_coef,Kh_coef,Kh_coef_lap1) &
  !$acc& pcopy(A2_o,CooCip,CooCim,D2,Mc,Mu,Ml) &
  !$acc& pcopyin(CVW) &
  !$acc& pcopyin(local_t_var) &
  !$acc& pcopyin(sendlist,sendlist_pl) &
  !$acc& pcopyin(sendinfo,sendinfo_pl) &
  !$acc& pcopyin(recvlist,recvlist_pl) &
  !$acc& pcopyin(recvinfo,recvinfo_pl) &
  !$acc& pcopyin(recvlist_r2r,sendlist_r2r) &
  !$acc& pcopyin(recvlist_r2p,sendlist_r2p) &
  !$acc& pcopyin(recvlist_p2r,sendlist_p2r) &
  !$acc& pcopyin(recvlist_sgp,sendlist_sgp) &
  !$acc& pcopyin(copyinfo_r2r,copyinfo_r2p,copyinfo_p2r,copyinfo_sgp) &
  !$acc& pcopyin(nsmax,nsmax_pl,ncmax_r2r,ncmax_r2p,ncmax_p2r,nrmax,nrmax_pl,ncmax_sgp) &
  !$acc& pcopyin(sfcvar,KSTR) &
  !$acc& pcopyin(ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2) &
  !$acc& pcopy(PRG_var,PRG_var1,DIAG_var)

  !--- history output at initial time
  if ( HIST_output_step0 ) then
     cstep = -1
     ctime = -dtime
     call history_vars
     call TIME_advance
     call history_out
  endif

  !
  do n = 1, TIME_LSTEP_MAX

     call DEBUG_rapstart('+Atmos')
     call dynstep
     call forcing
     call DEBUG_rapend  ('+Atmos')

     call DEBUG_rapstart('+History')
     call history_vars
     call TIME_advance

     !--- budget monitor
     call embudget_monitor
     call history_out

     ! if (n == TIME_LSTEP_MAX) then
     !    cdate = ""
     !    call restart_output( restart_output_basename )
     !    call diagvar_restart_output ( cdate )
     ! endif
     call DEBUG_rapend  ('+History')

  enddo

  !$acc end data

  cdate = ""
  call restart_output( restart_output_basename )
  call diagvar_restart_output ( cdate )

  !

  write(ADM_LOG_FID,*) '##### finish main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish main loop #####'
  endif

  call DEBUG_rapend('Main ALL')
#ifdef _FIPP_
  call fipp_stop()
#endif
  !#############################################################################

  call DEBUG_rapend('Total')
  call DEBUG_rapreport

  !--- all processes stop
  call ADM_proc_stop

  stop
end program prg_driver
!-------------------------------------------------------------------------------
