!-------------------------------------------------------------------------------
!>
!! History history variables
!!
!! @par Description
!!         This module prepares diagnostic/prognostic variavles for histroy output
!!
!! @author W.Yanase, S.Iga, H.Tomita, Y.Niwa
!!
!! @par History
!! @li      2007-06-27   Y.Niwa  : Imported from mod_history_vars_cfmip
!! @li      2007-07-31   K.Suzuki: add SPRINTARS variables
!! @li      2007-08-06   Y.Niwa  : bug fix  'th' => 'ml_th' 'rh' => 'ml_rh'
!! @li      2007-08-06   K.Suzuki: move SPRINTARS variables from mod_postphystep
!! @li      2007-08-20   W.Yanase: debug in "call getr1_in_mat"
!! @li      2007-11-22   T.Mitsui: Effective Radius is calculated in rd_driver => mp_driver
!!                                 and some variables are moved from aerovar 
!! @li      2007-11-30   T.Mitsui: del double use(trivial)
!! @li      2007-12-05   T.Mitsui: add radiative flux categolized as ISCCP-D series
!! @li      2008-02-19   T.Mitsui: Add output for slice-cloud normalized by cot 
!!                                 and trivial fix for output option
!! @li      2008-03-06   H.Tomita: add mp_wsm3_wcdiag.
!! @li      2008-05-24   T.Mitsui: add nc, nr, ni, ns, ng for 2moment option
!! @li      2008-06-13   T.Mitsui: change adiabat_diag => adiabat_diag2
!!                                 and add arguments (positive only cape, and all cin)
!! @li      2008-07-30   W.Yanase: add sl_u10m, sl_v10m, sl_tauu, sl_tauv
!! @li      2009-01-23   H.Tomita: add ml_omg, ml_hgt
!! @li      2009-04-14   H.Tomita: bugfix : zero clear of lwpqc,lwpqr.
!! @li      2009-04-14   T.Mitsui: sl_nc,r,i,s,g and sl_tauc,r,i,s,g
!!                                 sl_qi, sl_qg
!! @li      2009-07-10   H.Tomita: 1. The energy conservation terms are added.
!!                                 2. Some tag names were changed to avoid the confusion.
!!                                   'sl_swi' -> 'sl_swd_sfc'
!!                                   'sl_sswr' -> 'sl_swu_sfc'
!!                                   'sl_slwd' -> 'sl_lwd_sfc'
!!                                   'sl_slwu' -> 'sl_lwu_sfc'
!!                                   'sl_slh' ->  'sl_lh_sfc'
!!                                   'sl_ssh' ->  'sl_sh_sfc'
!!                                   'sl_sw_toai' ->  'sl_swd_toa'
!!                                   'sl_sw_toar' ->  'sl_swu_toa'
!!                                   'sl_lw_toa' ->  'sl_lwu_toa'
!! @li      2009-07-13   S.Iga   : Bug fix 'sl_lw_toa' ->  'sl_lwu_toa'
!! @li      2010-06-19   A.T.Noda: Allow to use a convection parameterization
!!                                 with an advanced microphysics schemes, such as G98, NSW?,
!! @li      2010-08-20   C.Kodama: land model output is filled with undef values over the ocean.
!! @li      2012-11-05 (H.Yashiro)  NICAM milestone project (Phase I:cleanup of shared module)
!!
!<
module mod_history_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_MAXFNAME, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: history_vars_setup
  public :: history_vars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save :: out_ucos_vcos = .false.
  logical, private, save :: out_omg       = .false.
  logical, private, save :: out_hgt       = .false.
  logical, private, save :: out_th        = .false.
  logical, private, save :: out_rh        = .false.

  logical, private, save :: out_cld_frac  = .false.
  logical, private, save :: out_cldw      = .false.
  logical, private, save :: out_qr        = .false.
  logical, private, save :: out_qs        = .false.
  logical, private, save :: out_qi        = .false.
  logical, private, save :: out_qg        = .false.
  logical, private, save :: out_cldi      = .false.
  logical, private, save :: out_slh       = .false.
  logical, private, save :: out_tau_uv    = .false.
  logical, private, save :: out_10m_uv    = .false.
  logical, private, save :: out_vap_atm   = .false.
  logical, private, save :: out_tem_atm   = .false.
  logical, private, save :: out_albedo    = .false.
  logical, private, save :: out_slp       = .false.
  logical, private, save :: out_850hPa    = .false.
  logical, private, save :: out_cape_cin  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine history_vars_setup
    use mod_history, only: &
       HIST_req_nmax, &
       item_save
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    do n = 1, HIST_req_nmax
       if(      item_save(n) == 'ml_ucos'            &
           .OR. item_save(n) == 'ml_vcos'            ) out_ucos_vcos = .true.
       if(      item_save(n) == 'ml_omg'             ) out_omg       = .true.
       if(      item_save(n) == 'ml_hgt'             ) out_hgt       = .true.
       if(      item_save(n) == 'ml_th'              ) out_th        = .true.
       if(      item_save(n) == 'ml_rh'              ) out_rh        = .true.

       if(      item_save(n) == 'sl_cld_frac'        ) out_cld_frac  = .true.
       if(      item_save(n) == 'sl_cldw'            ) out_cldw      = .true.
       if(      item_save(n) == 'sl_cldi'            ) out_cldi      = .true.
       if(      item_save(n) == 'sl_qr'              ) out_qr        = .true.
       if(      item_save(n) == 'sl_qs'              ) out_qs        = .true.
       if(      item_save(n) == 'sl_qi'              ) out_qi        = .true.
       if(      item_save(n) == 'sl_qg'              ) out_qg        = .true.

       if(      item_save(n) == 'sl_albedo'          ) out_albedo    = .true.

       if(      item_save(n) == 'sl_tauu'            &
           .OR. item_save(n) == 'sl_tauv'            &
           .OR. item_save(n) == 'sl_tauucos'         &
           .OR. item_save(n) == 'sl_tauvcos'         ) out_tau_uv    = .true.

       if(      item_save(n) == 'sl_u10m'            &
           .OR. item_save(n) == 'sl_v10m'            &
           .OR. item_save(n) == 'sl_ucos10m'         &
           .OR. item_save(n) == 'sl_vcos10m'         ) out_10m_uv    = .true.

       if(      item_save(n) == 'sl_u850'            &
           .OR. item_save(n) == 'sl_v850'            &
           .OR. item_save(n) == 'sl_w850'            &
           .OR. item_save(n) == 'sl_t850'            ) out_850hPa    = .true.

       if(      item_save(n) == 'sl_slh'             ) out_slh       = .true.
       if(      item_save(n) == 'sl_vap_atm'         ) out_vap_atm   = .true.
       if(      item_save(n) == 'sl_tem_atm'         ) out_tem_atm   = .true.
       if(      item_save(n) == 'sl_slp'             ) out_slp       = .true.
       if(      item_save(n) == 'sl_cape'            &
           .OR. item_save(n) == 'sl_cin'             ) out_cape_cin  = .true.
    enddo

    return
  end subroutine history_vars_setup

  !----------------------------------------------------------------------
  subroutine history_vars
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_KNONE,   &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EGRAV,     &
       CNST_RAIR,      &
       CNST_RVAP,      &
       CNST_TMELT,     &
       CNST_EPS_ZERO,  &
       CNST_UNDEF,     &
       CNST_DWATR
    use mod_gmtr, only :        &
         GMTR_area,             &
         GMTR_area_pl,          &
         GMTR_P_var,            &
         GMTR_P_var_pl,         &
         GMTR_P_IX,             &
         GMTR_P_IY,             &
         GMTR_P_IZ,             &
         GMTR_P_JX,             &
         GMTR_P_JY,             &
         GMTR_P_JZ,             &
         GMTR_P_LAT
    use mod_gtl, only :         &
         GTL_output_var2,       &
         GTL_output_var2_da,    &
         GTL_generate_uv,       &
         GTL_clip_region_1layer , &  ! [add] 2010.08.20 C.Kodama
       GTL_max, &
       GTL_min
    use mod_vmtr, only :        &
         VMTR_GSGAM2,           &
         VMTR_GSGAM2_pl,        &
         VMTR_GSGAM2H,          &
         VMTR_GSGAM2H_pl,       &
         VMTR_VOLUME,           &
         VMTR_VOLUME_pl
    use mod_prgvar, only :     &
         prgvar_get_withdiag
    use mod_sfcvar, only :    &
         sfcvar_get,          &
         sfcvar_get1,         &
         sfcvar_get2,         &
         I_ALBEDO_SFC,        &
         I_PRECIP_TOT,        &
         I_PRECIP_CP,         &
         I_EVAP_SFC,          &
         I_PRE_SFC,           &
         I_TEM_SFC,           &
         I_SH_FLUX_SFC,       &
         I_LH_FLUX_SFC,       &
         I_SFCRAD_ENERGY,     &
         I_TOARAD_ENERGY,     &
         I_EVAP_ENERGY,       &
         I_PRECIP_ENERGY,     &
         I_TAUX_SFC,          &
         I_TAUY_SFC,          &
         I_TAUZ_SFC,          &
         I_RFLUXS_SD,        &
         I_RFLUXS_SU,        &
         I_RFLUXS_LD,        &
         I_RFLUXS_LU,        &
         I_RFLUX_TOA_SD,     &
         I_RFLUX_TOA_SU,     &
         I_RFLUX_TOA_LD,     &
         I_RFLUX_TOA_LU,     &
         I_RFLUX_TOA_SD_C,   &
         I_RFLUX_TOA_SU_C,   &
         I_RFLUX_TOA_LD_C,   &
         I_RFLUX_TOA_LU_C,   &
         I_CUMFRC,           &
         I_VX10,             &
         I_VY10,             &
         I_VZ10,             &
         I_T2,               &
         I_Q2,               &
         INDEX_SEA,          &
         sfcvar, KSTR    ! openacc
    use mod_runconf, only : &
         LHV,LHS,           &
         RAIN_TYPE,         &
         CP_TYPE,           & ! 10/06/10 A.T.Noda
         NQW_MAX,           &
         TRC_VMAX,          &
         TRC_name,          &
         MP_TYPE,           & ! 07/05/08 H.tomita
         AE_TYPE,           & ! 07/07/30 K.Suzuki
         opt_2moment_water, & ! 08/05/24 T.Mitsui
         opt_incloud_aerosol,& ! 12/02/01 T.Seiki
         I_NC, I_NR, I_NI,  & ! 08/05/24 T.Mitsui
         I_NS, I_NG,        & ! 08/05/24 T.Mitsui
         NNW_STR, NNW_END,  & ! 08/05/24 T.Mitsui
         I_QV,              &
         I_QC,              &
         I_QR,              &
         I_QI,              &
         I_QS,              &
         I_QG,              &
         I_QH,              &
         NQW_STR,NQW_END,   &
         opt_carb_on, opt_dust_on, &! 12/02/01 T.Seiki
         opt_salt_on, opt_sulf_on, &! 12/02/01 T.Seiki
         NQDU_STR, NQDU_END,&       ! 12/02/01 T.Seiki
         NQDUIN_STR, NQDUIN_END, &  ! 12/02/01 T.Seiki
         NQCB_STR, NQCB_END, &      ! 12/02/01 T.Seiki
         NQCBIN_STR, NQCBIN_END, &  ! 12/02/01 T.Seiki
         NQSU_SO4,          &       ! 12/02/01 T.Seiki
         NQSU_SO2,          &       ! 12/02/01 T.Seiki
         NQSU_DMS,          &       ! 12/02/01 T.Seiki
         NQSO4IN,           &       ! 12/02/01 T.Seiki
         NQSA_STR, NQSA_END, &      ! 12/02/01 T.Seiki
         NQSAIN_STR, NQSAIN_END, &  ! 12/02/01 T.Seiki
         NRDIR,             &
         NRBND,             &
         NRBND_VIS,         &
         NRBND_NIR,         &
         NRBND_IR,          &
         NRDIR_DIRECT,      &
         NRDIR_DIFFUSE,     &
         NCRF,              &
         RAD_TYPE,          &  ! 07/12/05 T.Mitsui
         NTAU_ISCCP,        &  ! 07/12/05 T.Mitsui
         NPRES_ISCCP,       &  ! 07/12/05 T.Mitsui
         LAND_TYPE,         &  ! Y.Niwa add 070627
         OCEAN_TYPE            ! Y.Niwa add 070627
    use mod_bsstate, only : &
         phi, phi_pl
    use mod_cnvvar, only :  &
         cnvvar_p2d
    use mod_thrmdyn, only : &
         thrmdyn_th, &
         thrmdyn_th_ldim
    use mod_bndcnd, only : &
         bndcnd_thermo, &
         bndcnd_thermo_lall
    use mod_diagvar, only : &
         diagvar_get,       &
         I_GDCFRC,          &
         I_CUMCLW,          &
         I_UNCCN                 ! 08/02/19 T.Mitsui
    use mod_grd, only:          & 
         grd_zsfc,              &
         grd_zs,                &
         grd_zs_pl,             &
         GRD_vz,                &
         GRD_Z,                 &
         GRD_VEGINDX              ! [add] 2010.08.20 C.Kodama
    use mod_history, only : &
       history_in, &
       v_save,ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2
    implicit none

    real(8) :: rhog     (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhoge    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogq    (ADM_gall,   ADM_kall,ADM_lall   ,TRC_vmax)
    real(8) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(8) :: rho      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: pre      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: tem      (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vx       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vy       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vz       (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: w        (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: q        (ADM_gall,   ADM_kall,ADM_lall,   TRC_vmax)
    real(8) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(8) :: u   (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: v   (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: ucos(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: vcos(ADM_gall,ADM_kall,ADM_lall)
    real(8) :: wc  (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: omg (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: hgt (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: th  (ADM_gall,ADM_kall,ADM_lall)
!    real(8) :: rh  (ADM_gall,ADM_kall,ADM_lall)

    real(8) :: q_clw     (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: q_cli     (ADM_gall,ADM_kall,ADM_lall)
    real(8) :: qtot      (ADM_gall,ADM_kall,ADM_lall)
    ! real(8) :: cloud_frac(ADM_gall,ADM_kall,ADM_lall)

    real(8) :: rfluxs_sd(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_su(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_ld(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rfluxs_lu(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: rflux_toa_sd(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_su(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_ld(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_lu(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_su_c(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rflux_toa_lu_c(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: precip(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: sh_flux_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: lh_flux_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: sfcrad(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: toarad(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: precip_energy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: evap_energy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: taux(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tauy(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tauz(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: evap(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: rho_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: pre_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: tem_sfc(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vx10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vy10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: vz10m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: t2m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: q2m(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: u_850(ADM_gall,ADM_KNONE,ADM_lall)   ! [add] 20130705 R.Yoshida
    real(8) :: v_850(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: w_850(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: t_850(ADM_gall,ADM_KNONE,ADM_lall)

    real(8) :: albedo_sfc(ADM_gall,ADM_KNONE,ADM_lall,1:NRDIR,1:NRBND)
    real(8) :: slp(ADM_gall,ADM_KNONE,ADM_lall)
    real(8) :: albedo_sfc_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1:NRDIR,1:NRBND) ! 05/11/13 S.Iga

    integer, parameter :: cot_binmax=11

    real(8) :: dummy2D   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(8) :: dummy2D_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(8) :: dummy3D1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dummy3D2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dummy3D3_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dummy3D4_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dummy3D5_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rw(ADM_gall,ADM_lall)

    integer :: g, k, l, nq, K0
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE

    !$acc data &
    !$acc& create(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq) &
    !$acc& create(rho,pre,tem,vx,vy,vz,w,q) &
    !$acc& create(u,v,ucos,vcos,wc,omg,hgt,th) &
    !$acc& create(q_clw,q_cli,qtot) &
    !$acc& create(rfluxs_sd,rfluxs_su,rfluxs_ld,rfluxs_lu) &
    !$acc& create(rflux_toa_sd,rflux_toa_su,rflux_toa_ld,rflux_toa_lu,rflux_toa_su_c,rflux_toa_lu_c) &
    !$acc& create(precip) &
    !$acc& create(sh_flux_sfc,lh_flux_sfc,sfcrad,toarad,precip_energy,evap_energy) &
    !$acc& create(taux,tauy,tauz,evap,rho_sfc,pre_sfc,tem_sfc) &
    !$acc& create(vx10m,vy10m,vz10m,t2m,q2m,u_850,v_850,w_850,t_850) &
    !$acc& create(albedo_sfc,slp,dummy2D,rw) &
    !$acc& pcopyin(GRD_vz,GRD_zs) &
    !$acc& pcopyin(sfcvar,KSTR) &
    !$acc& pcopyin(ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2) &
    !$acc& pcopy(v_save)

    !--- get variables
    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    call bndcnd_thermo_lall( ADM_gall,   & ! [IN]
                             tem(:,:,:), & ! [INOUT]
                             rho(:,:,:), & ! [INOUT]
                             pre(:,:,:), & ! [INOUT]
                             phi(:,:,:)  ) ! [IN]
    !$acc wait

    !--- density, temperature & pressure
    do l = 1, ADM_lall
       call history_in( 'ml_rho',  rho(:,:,l) ) 
       call history_in( 'ml_tem',  tem(:,:,l) ) 
       call history_in( 'ml_pres', pre(:,:,l) ) 
    enddo

    !--- wind on cartesian
    do l = 1, ADM_lall
       call history_in( 'vx', vx(:,:,l) )
       call history_in( 'vy', vy(:,:,l) )
       call history_in( 'vz', vz(:,:,l) )
    enddo

    !--- zonal and meridonal wind
    call GTL_generate_uv( u,  dummy3D1_pl, & ! [OUT]
                          v,  dummy3D2_pl, & ! [OUT]
                          vx, dummy3D3_pl, & ! [IN]
                          vy, dummy3D4_pl, & ! [IN]
                          vz, dummy3D5_pl, & ! [IN]
                          icos = 0         ) ! [IN]

    do l = 1, ADM_lall
       call history_in( 'ml_u', u(:,:,l) )
       call history_in( 'ml_v', v(:,:,l) )
    enddo

    !--- vertical wind (at cell center)
    !$acc kernels pcopy(wc) pcopyin(w)
    do l = 1, ADM_lall
    do k = ADM_kmin-1, ADM_kmax+1
    do g = 1, ADM_gall
       if (( k >= ADM_kmin ) .and. ( k <= ADM_kmax )) then
          wc(g,k,l) = 0.5D0 * ( w(g,k,l) + w(g,k+1,l) )
       else
          wc(g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    do l = 1, ADM_lall
       call history_in( 'ml_w', wc(:,:,l) ) 
    enddo

    !--- zonal and meridonal wind with cos(phi)
    if (out_ucos_vcos) then
       call GTL_generate_uv( ucos, dummy3D1_pl, & ! [OUT]
                             vcos, dummy3D2_pl, & ! [OUT]
                             vx,   dummy3D3_pl, & ! [IN]
                             vy,   dummy3D4_pl, & ! [IN]
                             vz,   dummy3D5_pl, & ! [IN]
                             icos = 1           ) ! [IN]

       do l = 1, ADM_lall
          call history_in( 'ml_ucos', ucos(:,:,l) )  
          call history_in( 'ml_vcos', vcos(:,:,l) ) 
       enddo
    endif

    !--- omega
    if (out_omg) then
       !$acc kernels pcopy(omg) pcopyin(rho,wc)
       do l = 1, ADM_lall
          omg(:,:,l) = -CNST_EGRAV * rho(:,:,l) * wc(:,:,l)
       enddo
       !$acc end kernels
       do l = 1, ADM_lall
          call history_in( 'ml_omg', omg(:,:,l) ) 
       enddo
    endif

    !--- geopotential height
    ! [NOTE] H.Tomita Hydrostatic assumption
    if (out_hgt) then
       !$acc kernels pcopy(hgt) pcopyin(GRD_vz)
       do l = 1, ADM_lall
          hgt(:,:,l) = GRD_vz(:,:,l,GRD_Z)
       enddo
       !$acc end kernels
       do l = 1, ADM_lall
          call history_in( 'ml_hgt', hgt(:,:,l) ) 
       enddo
    endif

    !--- potential temperature
    if (out_th) then
       call thrmdyn_th_ldim( ADM_gall,   & ! [IN]
                             ADM_lall,   &
                             th (:,:,:), & ! [OUT]
                             tem(:,:,:), & ! [IN]
                             pre(:,:,:)  ) ! [IN]
       !$acc wait
       do l = 1, ADM_lall
          call history_in( 'ml_th', th(:,:,l) )
       enddo
    endif

    !--- relative humidity
!    if (out_rh) then
!       call moist_relative_humidity( rh(:,:,l),    & ! [OUT]
!                                     rho(:,:,l),   & ! [IN]
!                                     tem(:,:,l),   & ! [IN]
!                                     q(:,:,l,I_QV) ) ! [IN]
!
!       call history_in( 'ml_rh', rh(:,:,l) )
!    endif

    !--- tracers
    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
       call history_in( 'ml_'//TRC_name(nq), q(:,:,l,nq) )
    enddo
    enddo

    !--- hydrometeors
    !$acc kernels pcopy(q_clw,q_cli) pcopyin(q)
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       q_clw(g,k,l) = 0.D0
       q_cli(g,k,l) = 0.D0
       !$acc loop seq
       do nq = NQW_STR, NQW_END
          if ( nq == I_QC ) then
             q_clw(g,k,l) = q_clw(g,k,l) + q(g,k,l,nq)
          elseif( nq == I_QR ) then
             q_clw(g,k,l) = q_clw(g,k,l) + q(g,k,l,nq)
          elseif( nq == I_QI ) then
             q_cli(g,k,l) = q_cli(g,k,l) + q(g,k,l,nq)
          elseif( nq == I_QS ) then
             q_cli(g,k,l) = q_cli(g,k,l) + q(g,k,l,nq)
          elseif( nq == I_QG ) then
             q_cli(g,k,l) = q_cli(g,k,l) + q(g,k,l,nq)
          endif
       enddo
!       !$acc end loop
       qtot(g,k,l) = q_clw(g,k,l) + q_cli(g,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels
    do l = 1, ADM_lall
       call history_in( 'ml_qtot', qtot(:,:,l) )
    enddo

    ! do l = 1, ADM_lall
    !    cloud_frac(:,:,l) = 0.0D0
    !    where( qtot(:,:,l) >= 0.005D-3 ) !--- tompkins & craig
    !       cloud_frac(:,:,l) = 1.D0
    !    end where
    ! enddo

    if (out_cldw) then
       !$acc kernels pcopy(dummy2D) pcopyin(rho,q_clw, VMTR_VOLUME,GMTR_area)
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          dummy2D(g,K0,l) = 0.D0
          !$acc loop seq
          do k = ADM_kmin, ADM_kmax
             dummy2D(g,K0,l) = dummy2D(g,K0,l) &
                             + rho(g,k,l) * q_clw(g,k,l) * VMTR_VOLUME(g,k,l) / GMTR_area(g,l)
          enddo
!          !$acc end loop
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_cldw', dummy2D(:,:,l) )
       enddo
    endif

    if (out_cldi) then
       !$acc kernels pcopy(dummy2D) pcopyin(rho,q_cli, VMTR_VOLUME,GMTR_area)
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          dummy2D(g,K0,l) = 0.D0
          do k = ADM_kmin, ADM_kmax
             dummy2D(g,K0,l) = dummy2D(g,K0,l) &
                             + rho(g,k,l) * q_cli(g,k,l) * VMTR_VOLUME(g,k,l) / GMTR_area(g,l)
          enddo
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_cldi', dummy2D(:,:,l) )
       enddo
    endif

    !--- mass budget term
    call sfcvar_get( precip, dummy2D_pl, vid = I_PRECIP_TOT )
    call sfcvar_get( evap,   dummy2D_pl, vid = I_EVAP_SFC   )
    do l = 1, ADM_lall
       call history_in( 'sl_tppn', precip(:,:,l) )
       call history_in( 'sl_evap', evap  (:,:,l) )
    enddo

    !--- energy budget term
    call sfcvar_get( sfcrad,        dummy2D_pl, vid = I_SFCRAD_ENERGY )
    call sfcvar_get( toarad,        dummy2D_pl, vid = I_TOARAD_ENERGY )
    call sfcvar_get( evap_energy,   dummy2D_pl, vid = I_EVAP_ENERGY   )
    call sfcvar_get( precip_energy, dummy2D_pl, vid = I_PRECIP_ENERGY )
    call sfcvar_get( sh_flux_sfc,   dummy2D_pl, vid = I_SH_FLUX_SFC   )
    call sfcvar_get( lh_flux_sfc,   dummy2D_pl, vid = I_LH_FLUX_SFC   )
    do l = 1, ADM_lall
       call history_in( 'sl_rad_toa',     toarad       (:,:,l) )
       call history_in( 'sl_rad_sfc',     sfcrad       (:,:,l) )
       call history_in( 'sl_evap_energy', evap_energy  (:,:,l) )
       call history_in( 'sl_tppn_energy', precip_energy(:,:,l) )
       call history_in( 'sl_sh_sfc',      sh_flux_sfc  (:,:,l) )
       call history_in( 'sl_lh_sfc',      lh_flux_sfc  (:,:,l) )
    enddo

    !--- surface
    call sfcvar_get( tem_sfc, dummy2D_pl, vid = I_TEM_SFC  )
    call sfcvar_get( t2m,     dummy2D_pl, vid = I_T2       )
    call sfcvar_get( q2m,     dummy2D_pl, vid = I_Q2       )
    call sfcvar_get( taux,    dummy2D_pl, vid = I_TAUX_SFC )
    call sfcvar_get( tauy,    dummy2D_pl, vid = I_TAUY_SFC )
    call sfcvar_get( tauz,    dummy2D_pl, vid = I_TAUZ_SFC )
    call sfcvar_get( vx10m,   dummy2D_pl, vid = I_VX10     )
    call sfcvar_get( vy10m,   dummy2D_pl, vid = I_VY10     )
    call sfcvar_get( vz10m,   dummy2D_pl, vid = I_VZ10     )

    call sfcvar_get2( albedo_sfc, albedo_sfc_pl, I_ALBEDO_SFC, NRDIR, NRBND )

    call sv_pre_sfc_lall( ADM_gall,                 & ! [IN]
                          rho    (:,:,:),           & ! [IN]
                          pre    (:,:,:),           & ! [IN]
                          GRD_vz (:,:,:,GRD_Z),     & ! [IN]
                          GRD_zs (:,:,:,GRD_ZSFC), & ! [IN]
                          rho_sfc(:,:,:),          & ! [OUT]
                          pre_sfc(:,:,:)           ) ! [OUT]

    do l = 1, ADM_lall
       call history_in( 'sl_ps',      pre_sfc(:,:,l) )
       call history_in( 'sl_tem_sfc', tem_sfc(:,:,l) )
       call history_in( 'sl_t2m',     t2m    (:,:,l) )
       call history_in( 'sl_q2m',     q2m    (:,:,l) )

       call history_in( 'sl_albedo_sfc', albedo_sfc(:,:,l,NRDIR_DIRECT,NRBND_VIS) )
    enddo

    if (out_slp) then
       !$acc kernels pcopy(slp) pcopyin(pre_sfc,t2m, GRD_zs)
       do l = 1, ADM_lall
          slp(:,K0,l) = pre_sfc(:,K0,l) &
                      * exp( CNST_EGRAV * GRD_zs(:,K0,l,GRD_ZSFC) &
                      / ( CNST_RAIR * ( t2m(:,K0,l) + 0.5D0 * 0.005D0 * ( GRD_zs(:,K0,l,GRD_ZSFC) + 2.D0 ) ) ) )
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_slp', slp(:,:,l) )
       enddo
    endif

    if (out_850hPa) then   ! [add] 20130705 R.Yoshida
       call sv_plev_uvwt_lall( ADM_gall,        & ! [IN]
                               pre    (:,:,:),  & ! [IN]
                               u      (:,:,:),  & ! [IN]
                               v      (:,:,:),  & ! [IN]
                               w      (:,:,:),  & ! [IN]
                               tem    (:,:,:),  & ! [IN]
                               85000.d0,        & ! [IN]
                               u_850  (:,:,:), & ! [OUT]
                               v_850  (:,:,:), & ! [OUT]
                               w_850  (:,:,:), & ! [OUT]
                               t_850  (:,:,:)  ) ! [OUT]

       do l = 1, ADM_lall
          call history_in( 'sl_u850', u_850(:,:,l) )
          call history_in( 'sl_v850', v_850(:,:,l) )
          call history_in( 'sl_w850', w_850(:,:,l) )
          call history_in( 'sl_t850', t_850(:,:,l) )
       enddo
    endif

    if (out_tau_uv) then
       !$acc kernels pcopy(dummy2D) pcopyin(taux,tauy,tauz, GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = taux(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IX) &
                          + tauy(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IY) &
                          + tauz(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IZ)
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_tauu', dummy2D(:,:,l) )
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = dummy2D(:,K0,l) * cos( GMTR_P_var(:,K0,l,GMTR_P_LAT) )
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_tauucos', dummy2D(:,:,l) )
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(taux,tauy,tauz, GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = taux(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JX) &
                          + tauy(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JY) &
                          + tauz(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JZ)
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_tauv', dummy2D(:,:,l) )
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = dummy2D(:,K0,l) * cos( GMTR_P_var(:,K0,l,GMTR_P_LAT) )
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_tauvcos', dummy2D(:,:,l) )
       enddo
    endif

    if (out_10m_uv) then
       !$acc kernels pcopy(dummy2D) pcopyin(vx10m,vy10m,vz10m, GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = vx10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IX) &
                          + vy10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IY) &
                          + vz10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_IZ)
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_u10m', dummy2D(:,:,l) )  
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = dummy2D(:,K0,l) * cos( GMTR_P_var(:,K0,l,GMTR_P_LAT) )
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_ucos10m', dummy2D(:,:,l) )
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(vx10m,vy10m,vz10m, GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = vx10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JX) &
                          + vy10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JY) &
                          + vz10m(:,K0,l) * GMTR_P_var(:,K0,l,GMTR_P_JZ)
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_v10m', dummy2D(:,:,l) )
       enddo

       !$acc kernels pcopy(dummy2D) pcopyin(GMTR_P_var)
       do l = 1, ADM_lall
          dummy2D(:,K0,l) = dummy2D(:,K0,l) * cos( GMTR_P_var(:,K0,l,GMTR_P_LAT) )
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_vcos10m', dummy2D(:,:,l) )
       enddo
    endif

    if (out_vap_atm) then
       !$acc kernels pcopy(dummy2D) pcopyin(rho,q, VMTR_VOLUME,GMTR_area)
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          !$acc loop seq
          do k = ADM_kmin, ADM_kmax
             dummy2D(g,K0,l) = dummy2D(g,K0,l) &
                             + rho(g,k,l) * q(g,k,l,I_QV) * VMTR_VOLUME(g,k,l) / GMTR_area(g,l)
          enddo
!          !$acc end loop
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_vap_atm', dummy2D(:,:,l) )
       enddo
    endif

    if (out_tem_atm) then
       !$acc kernels pcopy(dummy2D) pcreate(rw) pcopyin(rho,tem, VMTR_VOLUME)
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          dummy2D(g,K0,l) = 0.D0
          rw(g,l)         = 0.D0 
          !$acc loop seq
          do k = ADM_kmin, ADM_kmax
             dummy2D(g,K0,l) = dummy2D(g,K0,l) + rho(g,k,l) * tem(g,k,l) * VMTR_VOLUME(g,k,l)
             rw     (g,l)    = rw     (g,l)    + rho(g,k,l)              * VMTR_VOLUME(g,k,l)
          enddo
!          !$acc end loop
          dummy2D(g,K0,l) = dummy2D(g,K0,l) / rw(g,l)
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_tem_atm', dummy2D(:,:,l) )
       enddo
    endif

    !--- radiation flux
    call sfcvar_get( rfluxs_sd,      dummy2D_pl, vid = I_RFLUXS_SD      )
    call sfcvar_get( rfluxs_su,      dummy2D_pl, vid = I_RFLUXS_SU      )
    call sfcvar_get( rfluxs_ld,      dummy2D_pl, vid = I_RFLUXS_LD      )
    call sfcvar_get( rfluxs_lu,      dummy2D_pl, vid = I_RFLUXS_LU      )
    call sfcvar_get( rflux_toa_sd,   dummy2D_pl, vid = I_RFLUX_TOA_SD   )
    call sfcvar_get( rflux_toa_su,   dummy2D_pl, vid = I_RFLUX_TOA_SU   )
    call sfcvar_get( rflux_toa_ld,   dummy2D_pl, vid = I_RFLUX_TOA_LD   )
    call sfcvar_get( rflux_toa_lu,   dummy2D_pl, vid = I_RFLUX_TOA_LU   )
    call sfcvar_get( rflux_toa_su_c, dummy2D_pl, vid = I_RFLUX_TOA_SU_C )
    call sfcvar_get( rflux_toa_lu_c, dummy2D_pl, vid = I_RFLUX_TOA_LU_C )
    do l = 1, ADM_lall
       call history_in( 'sl_swd_sfc',   rfluxs_sd     (:,:,l) )
       call history_in( 'sl_swu_sfc',   rfluxs_su     (:,:,l) )
       call history_in( 'sl_lwd_sfc',   rfluxs_ld     (:,:,l) )
       call history_in( 'sl_lwu_sfc',   rfluxs_lu     (:,:,l) )
       call history_in( 'sl_swd_toa',   rflux_toa_sd  (:,:,l) )
       call history_in( 'sl_swu_toa',   rflux_toa_su  (:,:,l) )
       call history_in( 'sl_lwu_toa',   rflux_toa_lu  (:,:,l) )
       call history_in( 'sl_lwd_toa',   rflux_toa_ld  (:,:,l) )
       call history_in( 'sl_lwu_toa_c', rflux_toa_lu_c(:,:,l) )
       call history_in( 'sl_swu_toa_c', rflux_toa_su_c(:,:,l) )
    enddo

    if (out_albedo) then
       !$acc kernels pcopy(dummy2D) pcopyin(rflux_toa_su,rflux_toa_sd)
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          dummy2D(g,K0,l) = min( rflux_toa_su(g,K0,l) / max( rflux_toa_sd(g,K0,l), CNST_EPS_ZERO ), 1.D0 )
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          call history_in( 'sl_albedo', dummy2D(:,:,l) )
       enddo
    endif

    !$acc end data

    return
  end subroutine history_vars

  subroutine sv_pre_sfc( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use mod_adm, only :  &
       kdim => ADM_kall,    &
       kmin => ADM_kmin
    use mod_cnst, only :  &
       CNST_EGRAV
    implicit none

    integer, intent(in)  :: ijdim
    real(8), intent(in)  :: rho    (ijdim,kdim)
    real(8), intent(in)  :: pre    (ijdim,kdim)
    real(8), intent(in)  :: z      (ijdim,kdim)
    real(8), intent(in)  :: z_srf  (ijdim)
    real(8), intent(out) :: rho_srf(ijdim)
    real(8), intent(out) :: pre_srf(ijdim)

    integer :: ij
    !---------------------------------------------------------------------------

    !--- surface density ( extrapolation )
    do ij = 1, ijdim
       rho_srf(ij) = rho(ij,kmin) &
                   - ( rho(ij,kmin+1)-rho(ij,kmin) ) * ( z(ij,kmin)-z_srf(ij) ) / ( z(ij,kmin+1)-z(ij,kmin) )
    enddo

    !--- surface pressure ( hydrostatic balance )
    do ij = 1, ijdim
       pre_srf(ij) = pre(ij,kmin) + 0.5D0 * ( rho_srf(ij)+rho(ij,kmin) ) * CNST_EGRAV * ( z(ij,kmin)-z_srf(ij) )
    enddo

    return
  end subroutine sv_pre_sfc

  !
  subroutine sv_pre_sfc_lall( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use mod_adm, only :  &
       ldim => ADM_lall,    &
       kdim => ADM_kall,    &
       kmin => ADM_kmin, &
       k0 => ADM_KNONE
    use mod_cnst, only :  &
       CNST_EGRAV
    implicit none

    integer, intent(in)  :: ijdim
    real(8), intent(in)  :: rho    (ijdim,kdim,ldim)
    real(8), intent(in)  :: pre    (ijdim,kdim,ldim)
    real(8), intent(in)  :: z      (ijdim,kdim,ldim)
    real(8), intent(in)  :: z_srf  (ijdim,k0,ldim)
    real(8), intent(out) :: rho_srf(ijdim,k0,ldim)
    real(8), intent(out) :: pre_srf(ijdim,k0,ldim)

    integer :: ij, l
    !---------------------------------------------------------------------------

    !$acc data pcopy(rho_srf,pre_srf) pcopyin(rho,pre,z,z_srf)

    !--- surface density ( extrapolation )
    !$acc kernels pcopy(rho_srf) pcopyin(rho,z,z_srf)
    do l = 1, ldim
    do ij = 1, ijdim
       rho_srf(ij,k0,l) = rho(ij,kmin,l) &
                   - ( rho(ij,kmin+1,l)-rho(ij,kmin,l) ) * ( z(ij,kmin,l)-z_srf(ij,k0,l) ) / ( z(ij,kmin+1,l)-z(ij,kmin,l) )
    enddo
    enddo
    !$acc end kernels

    !--- surface pressure ( hydrostatic balance )
    !$acc kernels pcopy(pre_srf) pcopyin(pre,rho_srf,rho,z,z_srf)
    do l = 1, ldim
    do ij = 1, ijdim
       pre_srf(ij,k0,l) = pre(ij,kmin,l) + 0.5D0 * ( rho_srf(ij,k0,l)+rho(ij,kmin,l) ) * CNST_EGRAV * ( z(ij,kmin,l)-z_srf(ij,k0,l) )
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine sv_pre_sfc_lall

  !-----------------------------------------------------------------------------
  ! [add] 20130705 R.Yoshida
  subroutine sv_plev_uvwt( &
       ijdim, &
       pre,   &
       u_z,   &
       v_z,   &
       w_z,   &
       t_z,   &
       plev,  &
       u_p,   &
       v_p,   &
       w_p,   &
       t_p    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    implicit none

    integer, intent(in)  :: ijdim
    real(8), intent(in)  :: pre(ijdim,kdim)
    real(8), intent(in)  :: u_z(ijdim,kdim)
    real(8), intent(in)  :: v_z(ijdim,kdim)
    real(8), intent(in)  :: w_z(ijdim,kdim)
    real(8), intent(in)  :: t_z(ijdim,kdim)
    real(8), intent(in)  :: plev
    real(8), intent(out) :: u_p(ijdim)
    real(8), intent(out) :: v_p(ijdim)
    real(8), intent(out) :: w_p(ijdim)
    real(8), intent(out) :: t_p(ijdim)

    integer :: kl(ijdim)
    integer :: ku(ijdim)

    real(8) :: wght_l, wght_u

    integer :: ij, k
    !---------------------------------------------------------------------------

    ! search z-level
    do ij = 1, ijdim
       do k = kmin, kdim
          if( pre(ij,k) < plev ) exit
       enddo
       if( k >= kdim ) stop "internal error! [sv_uvwp_850/mod_history_vars]"

       ku(ij) = k
       kl(ij) = k - 1
    enddo

    ! interpolate
    do ij = 1, ijdim
       wght_l = ( log(plev)           - log(pre(ij,ku(ij))) ) &
              / ( log(pre(ij,kl(ij))) - log(pre(ij,ku(ij))) )

       wght_u = ( log(pre(ij,kl(ij))) - log(plev)           ) &
              / ( log(pre(ij,kl(ij))) - log(pre(ij,ku(ij))) )

       u_p(ij) = wght_l * u_z(ij,kl(ij)) + wght_u * u_z(ij,ku(ij))
       v_p(ij) = wght_l * v_z(ij,kl(ij)) + wght_u * v_z(ij,ku(ij))
       w_p(ij) = wght_l * w_z(ij,kl(ij)) + wght_u * w_z(ij,ku(ij))
       t_p(ij) = wght_l * t_z(ij,kl(ij)) + wght_u * t_z(ij,ku(ij))
    enddo

    return
  end subroutine sv_plev_uvwt

  !-----------------------------------------------------------------------------
  subroutine sv_plev_uvwt_lall( &
       ijdim, &
       pre,   &
       u_z,   &
       v_z,   &
       w_z,   &
       t_z,   &
       plev,  &
       u_p,   &
       v_p,   &
       w_p,   &
       t_p    )
    use mod_adm, only: &
       ldim => ADM_lall, &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       k0 => ADM_KNONE
    implicit none

    integer, intent(in)  :: ijdim
    real(8), intent(in)  :: pre(ijdim,kdim,ldim)
    real(8), intent(in)  :: u_z(ijdim,kdim,ldim)
    real(8), intent(in)  :: v_z(ijdim,kdim,ldim)
    real(8), intent(in)  :: w_z(ijdim,kdim,ldim)
    real(8), intent(in)  :: t_z(ijdim,kdim,ldim)
    real(8), intent(in)  :: plev
    real(8), intent(out) :: u_p(ijdim,k0,ldim)
    real(8), intent(out) :: v_p(ijdim,k0,ldim)
    real(8), intent(out) :: w_p(ijdim,k0,ldim)
    real(8), intent(out) :: t_p(ijdim,k0,ldim)

    integer :: ku
    integer :: kl
    integer :: max_ku

    real(8) :: wght_l, wght_u

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    max_ku = kmin

    !$acc parallel pcopy(u_p,v_p,w_p,t_p) pcopyin(pre,u_z,v_z,w_z,t_z) reduction(max:max_ku)
    do l = 1, ldim
    do ij = 1, ijdim
       ! search z-level
       ku = kdim
       !$acc loop seq
       do k = kmin, kdim
          if( (k < ku) .and. (pre(ij,k,l) < plev) ) then
             ku = k
          endif
       enddo
!       !$acc end loop
       kl = ku - 1

       if ( max_ku < ku ) then
          max_ku = ku
       endif

       ! interpolate
       wght_l = ( log(plev)         - log(pre(ij,ku,l)) ) &
              / ( log(pre(ij,kl,l)) - log(pre(ij,ku,l)) )

       wght_u = ( log(pre(ij,kl,l)) - log(plev)         ) &
              / ( log(pre(ij,kl,l)) - log(pre(ij,ku,l)) )

       u_p(ij,k0,l) = wght_l * u_z(ij,kl,l) + wght_u * u_z(ij,ku,l)
       v_p(ij,k0,l) = wght_l * v_z(ij,kl,l) + wght_u * v_z(ij,ku,l)
       w_p(ij,k0,l) = wght_l * w_z(ij,kl,l) + wght_u * w_z(ij,ku,l)
       t_p(ij,k0,l) = wght_l * t_z(ij,kl,l) + wght_u * t_z(ij,ku,l)
    enddo
    enddo
    !$acc end parallel

    if( max_ku >= kdim ) stop "internal error! [sv_uvwp_850/mod_history_vars]"

    return
  end subroutine sv_plev_uvwt_lall

end module mod_history_vars
!-------------------------------------------------------------------------------
