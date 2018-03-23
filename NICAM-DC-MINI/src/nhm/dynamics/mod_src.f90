!-------------------------------------------------------------------------------
!
!+  Source calculation module
!
!-------------------------------------------------------------------------------
module mod_src
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module is for the caluculation of source terms
  !       in the nhm-model.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-08-11   Add sub[src_update_tracer] for tracer advection.
  !                07-04-25   K.Suzuki: branch for save_memory in src_update_tracer
  !                07-11-13   H.Tomita: Change the vertical advection limiter.
  !                                     apply the limiter from rho q to q.
  !                07-11-14   Y.Niwa: bug fix /rho_h => /rho_h_pl
  !                08-01-24   Y.Niwa: add src_update_tracer
  !                                   save_memory => TRC_SAVE_MEMORY
  !                08-04-28   Y.Niwa: add initialization and avoid zero-dividing. 
  !                09-05-26   Y.Yamada: add directive, derected by T.Asano
  !                11-09-27   T.Seiki: merge optimized routines for K by RIST and M.Terai
  !                12-03-29   T.Yamaura: optimized for K
  !                12-05-30   T.Yashiro: Change arguments from character to index/switch
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: src_advection_convergence_momentum
  public :: src_advection_convergence
  public :: src_flux_convergence

  public :: src_gradient
  public :: src_buoyancy

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_SRC_horizontal = 1 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_vertical   = 2 ! [add] H.Yashiro 20120530
  integer, public, parameter :: I_SRC_default    = 3 ! [add] H.Yashiro 20120530

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, parameter :: first_layer_remedy = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Advection convergence for momentum
  subroutine src_advection_convergence_momentum( &  
       vx,      vx_pl,      &
       vy,      vy_pl,      &
       vz,      vz_pl,      &
       w,       w_pl,       &
       rhog,    rhog_pl,    &
       rhogvx,  rhogvx_pl,  &
       rhogvy,  rhogvy_pl,  &
       rhogvz,  rhogvz_pl,  &
       rhogw,   rhogw_pl,   &
       grhogvx, grhogvx_pl, &
       grhogvy, grhogvy_pl, &
       grhogvz, grhogvz_pl, &
       grhogw,  grhogw_pl   )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_KNONE,   &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EOHM
    use mod_grd, only: &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_rscale, &
       GRD_e,      &
       GRD_e_pl,   &
       GRD_x,      &
       GRD_x_pl,   &
       GRD_afac,   &
       GRD_bfac,   &
       GRD_cfac,   &
       GRD_dfac
    use mod_runconf, only : &
       NON_HYDRO_ALPHA
    implicit none

    real(8), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: w    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: w_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(in)  :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8), intent(out) :: grhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(out) :: grhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: vvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: vvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: vvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
 
    real(8) :: dvvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dvvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: dvvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: dvvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: grhogwc   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: grhogwc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: prd, wc

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++src_advection_convergence_m')

    !$acc data &
    !$acc& pcopy(grhogvx,grhogvy,grhogvz,grhogw) &
    !$acc& pcopyin(vx,vy,vz,w,rhog,rhogvx,rhogvy,rhogvz,rhogw) &
    !$acc& create(vvx,vvy,vvz,dvvx,dvvy,dvvz,grhogwc) &
    !$acc& pcopyin(GRD_afac,GRD_bfac,GRD_cfac,GRD_dfac,GRD_e)

    !---< merge horizontal velocity & vertical velocity >

    !$acc kernels pcopy(vvx,vvy,vvz) pcopyin(w,vx,vy,vz,GRD_cfac,GRD_dfac,GRD_e) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin-1,ADM_kmax+1
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
       if ((k .ge. ADM_kmin).and.(k .le. ADM_kmax)) then
          wc = 0.5D0 * ( GRD_cfac(k) * w(g,k+1,l) &
                       + GRD_dfac(k) * w(g,k  ,l) )

          vvx(g,k,l) = vx(g,k,l) + wc * GRD_e(g,l,GRD_XDIR)
          vvy(g,k,l) = vy(g,k,l) + wc * GRD_e(g,l,GRD_YDIR)
          vvz(g,k,l) = vz(g,k,l) + wc * GRD_e(g,l,GRD_ZDIR)
       else
          vvx(g,k,l) = 0.D0
          vvy(g,k,l) = 0.D0
          vvz(g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             wc = 0.5D0 * ( GRD_cfac(k) * w_pl(g,k+1,l) &
                          + GRD_dfac(k) * w_pl(g,k  ,l) )

             vvx_pl(g,k,l) = vx_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_XDIR)
             vvy_pl(g,k,l) = vy_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_YDIR)
             vvz_pl(g,k,l) = vz_pl(g,k,l) + wc * GRD_e_pl(g,l,GRD_ZDIR)
          enddo
          enddo

          vvx_pl(:,ADM_kmin-1,l) = 0.D0
          vvx_pl(:,ADM_kmax+1,l) = 0.D0
          vvy_pl(:,ADM_kmin-1,l) = 0.D0
          vvy_pl(:,ADM_kmax+1,l) = 0.D0
          vvz_pl(:,ADM_kmin-1,l) = 0.D0
          vvz_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    !---< advection term for momentum >

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvx,    vvx_pl,    & ! [IN]
                                    dvvx,   dvvx_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvy,    vvy_pl,    & ! [IN]
                                    dvvy,   dvvy_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvz,    vvz_pl,    & ! [IN]
                                    dvvz,   dvvz_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]


    !---< coriolis force >
    !$acc kernels pcopy(dvvx,dvvy) pcopyin(vvx,vvy,rhog) 
    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax
    do g = 1, ADM_gall
       dvvx(g,k,l) = dvvx(g,k,l) - 2.D0 * rhog(g,k,l) * ( -CNST_EOHM * vvy(g,k,l) )
       dvvy(g,k,l) = dvvy(g,k,l) - 2.D0 * rhog(g,k,l) * (  CNST_EOHM * vvx(g,k,l) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !---< horizontalize & separate vertical velocity >
    !$acc kernels pcopy(grhogvx,grhogvy,grhogvz,grhogwc) pcopyin(dvvx,dvvy,dvvz,GRD_e) 
    do l = 1, ADM_lall
    do k = ADM_kmin-1, ADM_kmax+1
    do g = 1, ADM_gall
       if ( (k .ge. ADM_kmin) .and. (k .le. ADM_kmax) ) then
          prd = dvvx(g,k,l) * GRD_e(g,l,GRD_XDIR) &
               + dvvy(g,k,l) * GRD_e(g,l,GRD_YDIR) &
               + dvvz(g,k,l) * GRD_e(g,l,GRD_ZDIR)

          grhogvx(g,k,l) = dvvx(g,k,l) - prd * GRD_e(g,l,GRD_XDIR)
          grhogvy(g,k,l) = dvvy(g,k,l) - prd * GRD_e(g,l,GRD_YDIR)
          grhogvz(g,k,l) = dvvz(g,k,l) - prd * GRD_e(g,l,GRD_ZDIR)

          grhogwc(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=8)
       else
          grhogvx(g,k,l) = 0.D0
          grhogvy(g,k,l) = 0.D0
          grhogvz(g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(grhogw) pcopyin(grhogwc,GRD_afac,GRD_bfac) 
    !$acc loop
    do l = 1, ADM_lall
    !$acc loop vector(4)
    do k = ADM_kmin-1, ADM_kmax+1
    !$acc loop vector(64)
    do g = 1, ADM_gall
       if ( (k .ge. ADM_kmin+1) .and. (k .le. ADM_kmax) ) then
          grhogw(g,k,l) = 0.5D0 * ( GRD_afac(k) * grhogwc(g,k  ,l) &
                                  + GRD_bfac(k) * grhogwc(g,k-1,l) )
       else 
          grhogw (g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels
 
    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          !---< coriolis force >
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             dvvx_pl(g,k,l) = dvvx_pl(g,k,l) - 2.D0 * rhog_pl(g,k,l) * ( -CNST_EOHM * vvy_pl(g,k,l) )
             dvvy_pl(g,k,l) = dvvy_pl(g,k,l) - 2.D0 * rhog_pl(g,k,l) * (  CNST_EOHM * vvx_pl(g,k,l) )
          enddo
          enddo

          !---< horizontalize & separate vertical velocity >
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             prd = dvvx_pl(g,k,l) * GRD_e_pl(g,l,GRD_XDIR)  &
                 + dvvy_pl(g,k,l) * GRD_e_pl(g,l,GRD_YDIR)  &
                 + dvvz_pl(g,k,l) * GRD_e_pl(g,l,GRD_ZDIR)

             grhogvx_pl(g,k,l) = dvvx_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_XDIR)
             grhogvy_pl(g,k,l) = dvvy_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_YDIR)
             grhogvz_pl(g,k,l) = dvvz_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_ZDIR)

             grhogwc_pl(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=8)
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             grhogw_pl(g,k,l) = 0.5D0 * ( GRD_afac(k) * grhogwc_pl(g,k  ,l) &
                                        + GRD_bfac(k) * grhogwc_pl(g,k-1,l) )
          enddo
          enddo

          grhogvx_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvx_pl(:,ADM_kmax+1,l) = 0.D0
          grhogvy_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvy_pl(:,ADM_kmax+1,l) = 0.D0
          grhogvz_pl(:,ADM_kmin-1,l) = 0.D0
          grhogvz_pl(:,ADM_kmax+1,l) = 0.D0
          grhogw_pl (:,ADM_kmin-1,l) = 0.D0
          grhogw_pl (:,ADM_kmin  ,l) = 0.D0
          grhogw_pl (:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++src_advection_convergence_m')

    return
  end subroutine src_advection_convergence_momentum

  !-----------------------------------------------------------------------------
  !> Advection convergence
  subroutine src_advection_convergence( &
       rhogvx,   rhogvx_pl,   &
       rhogvy,   rhogvy_pl,   &
       rhogvz,   rhogvz_pl,   &
       rhogw,    rhogw_pl,    &
       scl,      scl_pl,      &
       grhogscl, grhogscl_pl, &
       fluxtype               )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    implicit none

    real(8), intent(in)  :: rhogvx     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz     (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw      (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: scl        (ADM_gall,   ADM_kall,ADM_lall   ) ! scalar
    real(8), intent(in)  :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhogscl   (ADM_gall,   ADM_kall,ADM_lall   ) ! scalar tendency
    real(8), intent(out) :: grhogscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in)  :: fluxtype                                    ! scheme type
                                                                        ! I_SRC_horizontal : horizontal convergence
                                                                        ! I_SRC_vertical   : vertical convergence
                                                                        ! I_SRC_default    : both of them

    real(8) :: rhogvxscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvxscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvyscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvyscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvzscl   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvzscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogwscl    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogwscl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++src_advection_convergence')

    !$acc data &
    !$acc& pcopy(grhogscl) &
    !$acc& pcopyin(rhogvx,rhogvy,rhogvz,rhogw,scl) &
    !$acc& create(rhogvxscl,rhogvyscl,rhogvzscl,rhogwscl) &
    !$acc& pcopyin(GRD_afac,GRD_bfac)

    !$acc kernels pcopy(rhogvxscl,rhogvyscl,rhogvzscl) pcopyin(rhogvx,rhogvy,rhogvz,scl) 
    rhogvxscl(:,:,:) = rhogvx(:,:,:) * scl(:,:,:)
    rhogvyscl(:,:,:) = rhogvy(:,:,:) * scl(:,:,:)
    rhogvzscl(:,:,:) = rhogvz(:,:,:) * scl(:,:,:)
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       rhogvxscl_pl(:,:,:) = rhogvx_pl(:,:,:) * scl_pl(:,:,:)
       rhogvyscl_pl(:,:,:) = rhogvy_pl(:,:,:) * scl_pl(:,:,:)
       rhogvzscl_pl(:,:,:) = rhogvz_pl(:,:,:) * scl_pl(:,:,:)
    endif

    ! rhogwscl = rhow * e_w ( at half level ).
    if ( fluxtype == I_SRC_horizontal ) then

       !$acc kernels pcopy(rhogwscl) 
       rhogwscl(:,:,:) = 0.D0
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          rhogwscl_pl(:,:,:) = 0.D0
       endif

    else

       !$acc kernels pcopy(rhogwscl) pcopyin(rhogw,scl,GRD_afac,GRD_bfac) 
       !$acc loop gang
       do l = 1, ADM_lall
       !$acc loop gang vector(4)
       do k = ADM_kmin, ADM_kmax+1
       !$acc loop gang vector(64)
       do g = 1, ADM_gall
          rhogwscl(g,k,l) = rhogw(g,k,l) * 0.5D0 * ( GRD_afac(k) * scl(g,k,  l) &
                                                   + GRD_bfac(k) * scl(g,k-1,l) )
       enddo
       enddo
       enddo
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             rhogwscl_pl(:,k,l) = rhogw_pl(:,k,l) * 0.5D0 * ( GRD_afac(k) * scl_pl(:,k  ,l) &
                                                            + GRD_bfac(k) * scl_pl(:,k-1,l) )
          enddo
          enddo
       endif

    endif

    !--- flux convergence step
    call src_flux_convergence( rhogvxscl, rhogvxscl_pl, & !--- [IN]
                               rhogvyscl, rhogvyscl_pl, & !--- [IN]
                               rhogvzscl, rhogvzscl_pl, & !--- [IN]
                               rhogwscl,  rhogwscl_pl,  & !--- [IN]
                               grhogscl,  grhogscl_pl,  & !--- [OUT]
                               fluxtype                 ) !--- [IN]

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++src_advection_convergence')

    return
  end subroutine src_advection_convergence

  !-----------------------------------------------------------------------------
  !> Flux convergence calculation
  !! 1. Horizontal flux convergence is calculated by using rhovx, rhovy, and 
  !!    rhovz which are defined at cell center (vertical) and A-grid (horizontal).
  !! 2. Vertical flux convergence is calculated by using rhovx, rhovy, rhovz, and rhow.
  !! 3. rhovx, rhovy, and rhovz can be replaced by rhovx*h, rhovy*h, and rhovz*h, respectively.
  !! 4. fluxtype can be set as below.
  subroutine src_flux_convergence( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       grhog,  grhog_pl,  &
       fluxtype           )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz
    use mod_vmtr, only: &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSH,       &
       VMTR_RGSH_pl,    &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_oprt, only: &
       OPRT_divergence
    implicit none

    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vz ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grhog    (ADM_gall,   ADM_kall,ADM_lall   ) ! source
    real(8), intent(out) :: grhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in)  :: fluxtype                                    ! scheme type
                                                                        ! I_SRC_horizontal : horizontal convergence
                                                                        ! I_SRC_vertical   : vertical convergence
                                                                        ! I_SRC_default    : both of them

    real(8) :: div_rhogvh   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: div_rhogvh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: rhogvx_vm   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvx_vm_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvy_vm   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvy_vm_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogvz_vm   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogvz_vm_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: rhogw_vm    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8) :: rhogw_vm_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++src_flux_convergence')

    !$acc data &
    !$acc& pcopy(grhog) &
    !$acc& pcopyin(rhogvx,rhogvy,rhogvz,rhogw) &
    !$acc& create(div_rhogvh,rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm) &
    !$acc& pcopyin(VMTR_C2Wfact,VMTR_RGAMH,VMTR_RGAM,VMTR_RGSH,GRD_rdgz)

    ! boundary condition
    !$acc kernels pcopy(rhogw_vm) 
    do l = 1, ADM_lall
    do g = 1, ADM_gall
       rhogw_vm(g,ADM_kmin  ,l) = 0.D0
       rhogw_vm(g,ADM_kmax+1,l) = 0.D0
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          rhogw_vm_pl(:,ADM_kmin  ,l) = 0.D0
          rhogw_vm_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    if ( fluxtype == I_SRC_horizontal ) then ! Horizontal

       !$acc kernels pcopy(rhogw_vm) pcopyin(rhogvx,rhogvy,rhogvz,VMTR_C2Wfact,VMTR_RGAMH) 
       !$acc loop gang vector(1)
       do l = 1, ADM_lall
       !$acc loop gang vector(4)
       do k = ADM_kmin+1, ADM_kmax
       !$acc loop gang vector(64)
       do g = 1, ADM_gall
          rhogw_vm(g,k,l) = ( VMTR_C2Wfact(g,1,k,l) * rhogvx(g,k  ,l) &
                            + VMTR_C2Wfact(g,2,k,l) * rhogvx(g,k-1,l) &
                            + VMTR_C2Wfact(g,3,k,l) * rhogvy(g,k  ,l) &
                            + VMTR_C2Wfact(g,4,k,l) * rhogvy(g,k-1,l) &
                            + VMTR_C2Wfact(g,5,k,l) * rhogvz(g,k  ,l) &
                            + VMTR_C2Wfact(g,6,k,l) * rhogvz(g,k-1,l) &
                            ) * VMTR_RGAMH(g,k,l) ! horizontal contribution
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(rhogvx_vm,rhogvy_vm,rhogvz_vm) pcopyin(rhogvx,rhogvy,rhogvz,VMTR_RGAM) 
       !$acc loop gang vector(1)
       do l = 1, ADM_lall
       !$acc loop gang vector(2)
       do k = 1, ADM_kall
       !$acc loop gang vector(128)
       do g = 1, ADM_gall
          rhogvx_vm(g,k,l) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvy_vm(g,k,l) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvz_vm(g,k,l) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
       enddo
       enddo
       enddo
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,1,k,l) * rhogvx_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,2,k,l) * rhogvx_pl(g,k-1,l) &
                                  + VMTR_C2Wfact_pl(g,3,k,l) * rhogvy_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,4,k,l) * rhogvy_pl(g,k-1,l) &
                                  + VMTR_C2Wfact_pl(g,5,k,l) * rhogvz_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,6,k,l) * rhogvz_pl(g,k-1,l) &
                                  ) * VMTR_RGAMH_pl(g,k,l) ! horizontal contribution
          enddo
          enddo
          enddo

          rhogvx_vm_pl(:,:,:) = rhogvx_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          rhogvy_vm_pl(:,:,:) = rhogvy_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          rhogvz_vm_pl(:,:,:) = rhogvz_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       endif

       !--- Horizontal flux convergence
       call OPRT_divergence( div_rhogvh, div_rhogvh_pl, & !--- [OUT]
                             rhogvx_vm,  rhogvx_vm_pl,  & !--- [IN]
                             rhogvy_vm,  rhogvy_vm_pl,  & !--- [IN]
                             rhogvz_vm,  rhogvz_vm_pl   ) !--- [IN]

    elseif( fluxtype == I_SRC_vertical ) then ! Vertical

       !$acc kernels pcopy(rhogw_vm) pcopyin(rhogw,VMTR_RGSH) 
       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          rhogw_vm(g,k,l) = rhogw(g,k,l) * VMTR_RGSH(g,k,l) ! vertical contribution
       enddo
       enddo
       enddo
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,k,l) = rhogw_pl(g,k,l) * VMTR_RGSH_pl(g,k,l) ! vertical contribution
          enddo
          enddo
          enddo
       endif

       !--- Horizontal flux convergence
       !$acc kernels pcopy(div_rhogvh) 
       div_rhogvh   (:,:,:) = 0.D0
       !$acc end kernels
       div_rhogvh_pl(:,:,:) = 0.D0

    else ! Default

       !--- horizontal + vertical contribution to rhogw_vm
       !$acc kernels pcopy(rhogw_vm) &
       !$acc& pcopyin(rhogvx,rhogvy,rhogvz,rhogw,vmtr_c2wfact,vmtr_rgamh,vmtr_rgsh) 
       !$acc loop gang vector(1)
       do l = 1, ADM_lall
       !$acc loop gang vector(4)
       do k = ADM_kmin+1, ADM_kmax
       !$acc loop gang vector(64)
       do g = 1, ADM_gall
          rhogw_vm(g,k,l) = ( VMTR_C2Wfact(g,1,k,l) * rhogvx(g,k  ,l) &
                            + VMTR_C2Wfact(g,2,k,l) * rhogvx(g,k-1,l) &
                            + VMTR_C2Wfact(g,3,k,l) * rhogvy(g,k  ,l) &
                            + VMTR_C2Wfact(g,4,k,l) * rhogvy(g,k-1,l) &
                            + VMTR_C2Wfact(g,5,k,l) * rhogvz(g,k  ,l) &
                            + VMTR_C2Wfact(g,6,k,l) * rhogvz(g,k-1,l) &
                            ) * VMTR_RGAMH(g,k,l)                     & ! horizontal contribution
                          + rhogw(g,k,l) * VMTR_RGSH(g,k,l)             ! vertical   contribution
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(rhogvx_vm,rhogvy_vm,rhogvz_vm) pcopyin(rhogvx,rhogvy,rhogvz,VMTR_RGAM) 
       !$acc loop gang vector(1)
       do l = 1, ADM_lall
       !$acc loop gang vector(2)
       do k = 1, ADM_kall
       !$acc loop gang vector(128)
       do g = 1, ADM_gall
          rhogvx_vm(g,k,l) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvy_vm(g,k,l) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvz_vm(g,k,l) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
       enddo
       enddo
       enddo
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,1,k,l) * rhogvx_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,2,k,l) * rhogvx_pl(g,k-1,l) &
                                  + VMTR_C2Wfact_pl(g,3,k,l) * rhogvy_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,4,k,l) * rhogvy_pl(g,k-1,l) &
                                  + VMTR_C2Wfact_pl(g,5,k,l) * rhogvz_pl(g,k  ,l) &
                                  + VMTR_C2Wfact_pl(g,6,k,l) * rhogvz_pl(g,k-1,l) &
                                  ) * VMTR_RGAMH_pl(g,k,l)                        & ! horizontal contribution
                                + rhogw_pl(g,k,l) * VMTR_RGSH_pl(g,k,l)             ! vertical   contribution
          enddo
          enddo
          enddo

          rhogvx_vm_pl(:,:,:) = rhogvx_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          rhogvy_vm_pl(:,:,:) = rhogvy_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
          rhogvz_vm_pl(:,:,:) = rhogvz_pl(:,:,:) * VMTR_RGAM_pl(:,:,:)
       endif

       !--- Horizontal flux convergence
       call OPRT_divergence( div_rhogvh, div_rhogvh_pl, & !--- [OUT]
                             rhogvx_vm,  rhogvx_vm_pl,  & !--- [IN]
                             rhogvy_vm,  rhogvy_vm_pl,  & !--- [IN]
                             rhogvz_vm,  rhogvz_vm_pl   ) !--- [IN]

    endif

    !--- Total flux convergence
    !$acc kernels pcopy(grhog) pcopyin(div_rhogvh,rhogw_vm,GRD_rdgz) 
    !$acc loop gang vector(1)
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin-1, ADM_kmax+1
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
       if (( k .ge. ADM_kmin ) .and. ( k .le. ADM_kmax )) then
          grhog(g,k,l) = - div_rhogvh(g,k,l) &
               - ( rhogw_vm(g,k+1,l)-rhogw_vm(g,k,l) ) * GRD_rdgz(k)
       else
          grhog(g,k,l) = 0.D0
       endif
    enddo
!     !$acc end loop
    enddo
!     !$acc end loop
    enddo
!     !$acc end loop
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             grhog_pl(:,k,l) = - div_rhogvh_pl(:,k,l) &
                               - ( rhogw_vm_pl(:,k+1,l)-rhogw_vm_pl(:,k,l) ) * GRD_rdgz(k)
          enddo

          grhog_pl(:,ADM_kmin-1,l) = 0.D0
          grhog_pl(:,ADM_kmax+1,l) = 0.D0
       enddo
    endif

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++src_flux_convergence')

    return
  end subroutine src_flux_convergence

  !-----------------------------------------------------------------------------
  !> Gradient operator
  subroutine src_gradient( &
       e,    e_pl,    &
       gex,  gex_pl,  &
       gey,  gey_pl,  &
       gez,  gez_pl,  &
       gevz, gevz_pl, &
       grad_type      )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz, &
       GRD_rdgzh
    use mod_vmtr, only: &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GAM2H,      &
       VMTR_GAM2H_pl,   &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_oprt, only: &
       OPRT_gradient,          &
       OPRT_horizontalize_vec
    implicit none

    real(8), intent(in)  :: e      (ADM_gall   ,ADM_kall,ADM_lall   ) ! phi * G^{1/2} * gamma^2
    real(8), intent(in)  :: e_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gex    (ADM_gall   ,ADM_kall,ADM_lall   ) ! horizontal grad. ( x-comp. )
    real(8), intent(out) :: gex_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gey    (ADM_gall   ,ADM_kall,ADM_lall   ) ! horizontal grad. ( y-comp. )
    real(8), intent(out) :: gey_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gez    (ADM_gall   ,ADM_kall,ADM_lall   ) ! horizontal grad. ( z-comp. )
    real(8), intent(out) :: gez_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gevz   (ADM_gall   ,ADM_kall,ADM_lall   ) ! vertical grad.   ( half point )
    real(8), intent(out) :: gevz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in)  :: grad_type                                 ! scheme type
                                                                      ! I_SRC_horizontal : horizontal gradient
                                                                      ! I_SRC_vertical   : vertical gradient ( no use )
                                                                      ! I_SRC_default    : both of them

    real(8) :: fex_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fex_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: fey_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fey_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8) :: fez_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: fez_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: gee   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8) :: gee_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++src_gradient')

    !$acc data &
    !$acc& pcopy(gex,gey,gez,gevz) &
    !$acc& pcopyin(e) &
    !$acc& create(fex_h,fey_h,fez_h,gee) &
    !$acc& pcopyin(VMTR_C2Wfact,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_GAM2H) &
    !$acc& pcopyin(GRD_rdgz,GRD_rdgzh)

    !------ horizontal gradient without mountain
    !$acc kernels pcopy(gee) pcopyin(e, VMTR_RGAM) 
    gee(:,:,:) = e(:,:,:) * VMTR_RGAM(:,:,:)
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl) then
       gee_pl(:,:,:) = e_pl(:,:,:) * VMTR_RGAM_pl(:,:,:) 
    endif

    call OPRT_gradient( gex, gex_pl, & ! [OUT]
                        gey, gey_pl, & ! [OUT]
                        gez, gez_pl, & ! [OUT]
                        gee, gee_pl  ) ! [IN]

    !--- horizontal gradient
    !$acc kernels pcopy(fex_h,fey_h,fez_h) pcopyin(e, VMTR_C2Wfact,VMTR_RGAMH) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin, ADM_kmax+1
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
       fex_h(g,k,l) = ( VMTR_C2Wfact(g,1,k,l) * e(g,k  ,l) &
                      + VMTR_C2Wfact(g,2,k,l) * e(g,k-1,l) ) * VMTR_RGAMH(g,k,l)
       fey_h(g,k,l) = ( VMTR_C2Wfact(g,3,k,l) * e(g,k  ,l) &
                      + VMTR_C2Wfact(g,4,k,l) * e(g,k-1,l) ) * VMTR_RGAMH(g,k,l)
       fez_h(g,k,l) = ( VMTR_C2Wfact(g,5,k,l) * e(g,k  ,l) &
                      + VMTR_C2Wfact(g,6,k,l) * e(g,k-1,l) ) * VMTR_RGAMH(g,k,l)  
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(gex,gey,gez) pcopyin(fex_h,fey_h,fez_h, GRD_rdgz) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
          gex(g,k,l) = gex(g,k,l) + ( fex_h(g,k+1,l) - fex_h(g,k,l) ) * GRD_rdgz(k)
          gey(g,k,l) = gey(g,k,l) + ( fey_h(g,k+1,l) - fey_h(g,k,l) ) * GRD_rdgz(k)
          gez(g,k,l) = gez(g,k,l) + ( fez_h(g,k+1,l) - fez_h(g,k,l) ) * GRD_rdgz(k)
    enddo
    enddo
    enddo
    !$acc end kernels

    !--- At the lowest layer, do not use the extrapolation value!
    if ( first_layer_remedy ) then
       !$acc kernels pcopy(gex,gey,gez) 
       do l = 1, ADM_lall
       do g = 1, ADM_gall
          gex(g,ADM_kmin,l) = gex(g,ADM_kmin+1,l)
          gey(g,ADM_kmin,l) = gey(g,ADM_kmin+1,l)
          gez(g,ADM_kmin,l) = gez(g,ADM_kmin+1,l)
       enddo
       enddo
       !$acc end kernels
    endif

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall_pl
             fex_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,1,k,l) * e_pl(g,k  ,l) &
                               + VMTR_C2Wfact_pl(g,2,k,l) * e_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)
             fey_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,3,k,l) * e_pl(g,k  ,l) &
                               + VMTR_C2Wfact_pl(g,4,k,l) * e_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)
             fez_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,5,k,l) * e_pl(g,k  ,l) &
                               + VMTR_C2Wfact_pl(g,6,k,l) * e_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)  
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax
             gex_pl(:,k,l) = gex_pl(:,k,l) + ( fex_h_pl(:,k+1,l) - fex_h_pl(:,k,l) ) * GRD_rdgz(k)
             gey_pl(:,k,l) = gey_pl(:,k,l) + ( fey_h_pl(:,k+1,l) - fey_h_pl(:,k,l) ) * GRD_rdgz(k)
             gez_pl(:,k,l) = gez_pl(:,k,l) + ( fez_h_pl(:,k+1,l) - fez_h_pl(:,k,l) ) * GRD_rdgz(k)
          enddo

          !--- At the lowest layer, do not use the extrapolation value!
          if ( first_layer_remedy ) then
              gex_pl(:,ADM_kmin,l) = gex_pl(:,ADM_kmin+1,l)
              gey_pl(:,ADM_kmin,l) = gey_pl(:,ADM_kmin+1,l)
              gez_pl(:,ADM_kmin,l) = gez_pl(:,ADM_kmin+1,l)
          endif
       enddo
    endif

    !--- horizontalize
    call OPRT_horizontalize_vec( gex, gex_pl, & ! [INOUT]
                                 gey, gey_pl, & ! [INOUT]
                                 gez, gez_pl  ) ! [INOUT]

    !--- vertical gradient ( note : half points )
    if ( grad_type == I_SRC_horizontal ) then

       !$acc kernels pcopy(gevz) 
       gevz(:,:,:) = 0.D0
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          gevz_pl(:,:,:) = 0.D0
       endif

    else

       !$acc kernels pcopy(gevz) pcopyin(e, VMTR_RGSGAM2,VMTR_GAM2H,GRD_rdgzh) 
       !$acc loop gang
       do l = 1, ADM_lall
       !$acc loop gang vector(4)
       do k = ADM_kmin-1, ADM_kmax+1
       !$acc loop gang vector(64)
       do g = 1, ADM_gall
          if (( k >= ADM_kmin ) .and. ( k <= ADM_kmax+1 )) then
             gevz(g,k,l) = ( e(g,k  ,l) * VMTR_RGSGAM2(g,k  ,l) &
                           - e(g,k-1,l) * VMTR_RGSGAM2(g,k-1,l) &
                           ) * GRD_rdgzh(k) * VMTR_GAM2H(g,k,l)
          else 
             gevz(g,k,l) = 0.D0
          endif
       enddo
       enddo
       enddo
       !$acc end kernels

       if ( ADM_prc_me == ADM_prc_pl ) then
          do l = 1, ADM_lall_pl
            do k = ADM_kmin, ADM_kmax+1
            do g = 1, ADM_gall_pl
               gevz_pl(g,k,l) = ( e_pl(g,k  ,l) * VMTR_RGSGAM2_pl(g,k  ,l) &
                                - e_pl(g,k-1,l) * VMTR_RGSGAM2_pl(g,k-1,l) &
                                ) * GRD_rdgzh(k) * VMTR_GAM2H_pl(g,k,l)
            enddo
            enddo
            gevz_pl(:,ADM_kmin-1,l) = 0.D0
          enddo
       endif

    endif

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++src_gradient')

    return
  end subroutine src_gradient

  !-----------------------------------------------------------------------------
  !> Calculation of buoyacy force
  !> NOTICE : Upward direction is positive for gbz.
  subroutine src_buoyancy( &
       rhog, rhog_pl, &
       gbz,  gbz_pl   )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_EGRAV
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl
    implicit none

    real(8), intent(in)  :: rhog   (ADM_gall   ,ADM_kall,ADM_lall   ) ! density perturbation ( gam2 X G^{1/2} )
    real(8), intent(in)  :: rhog_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: gbz    (ADM_gall   ,ADM_kall,ADM_lall   ) ! buoyancy force  at half level
    real(8), intent(out) :: gbz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++src_buoyancy')

    !$acc kernels pcopy(gbz) pcopyin(rhog, VMTR_RGSGAM2,VMTR_GSGAM2H, GRD_afac,GRD_bfac) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin-1, ADM_kmax+1
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
       if (( k >= ADM_kmin ) .and. ( k <= ADM_kmax+1 )) then
          gbz(g,k,l) = -0.5D0 * ( GRD_afac(k) * rhog(g,k  ,l) * VMTR_RGSGAM2(g,k  ,l) &
                                + GRD_bfac(k) * rhog(g,k-1,l) * VMTR_RGSGAM2(g,k-1,l) &
                                ) * VMTR_GSGAM2H(g,k,l) * CNST_EGRAV
       else
          gbz(g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             gbz_pl(g,k,l) = -0.5D0 * ( GRD_afac(k) * rhog_pl(g,k  ,l) * VMTR_RGSGAM2_pl(g,k  ,l) &
                                      + GRD_bfac(k) * rhog_pl(g,k-1,l) * VMTR_RGSGAM2_pl(g,k-1,l) &
                                      ) * VMTR_GSGAM2H_pl(g,k,l) * CNST_EGRAV
          enddo
          enddo
          gbz_pl(:,ADM_kmin-1,l) = 0.D0
       enddo
    endif

    !--- !$acc wait

    call DEBUG_rapend('++++src_buoyancy')

    return
  end subroutine src_buoyancy

end module mod_src
