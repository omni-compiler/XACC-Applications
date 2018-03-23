!-------------------------------------------------------------------------------
!>
!! 3D Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators using vertical metrics.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)    Imported from igdc-4.33
!! @li      2011-09-27 (T.Seiki)     merge optimization by RIST and M.Terai 
!!
!<
module mod_oprt3d
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
  public :: OPRT3D_divdamp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OPRT3D_divdamp( &
       grdx,   grdx_pl,   &
       grdy,   grdy_pl,   &
       grdz,   grdz_pl,   &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl   )
    use mod_adm, only: &
       ADM_W,        &
       ADM_TI,       &
       ADM_TJ,       &
       ADM_AI,       &
       ADM_AIJ,      &
       ADM_AJ,       &
       ADM_prc_tab,  &
       ADM_prc_me,   &
       ADM_prc_pl,   &
       ADM_rgn_vnum, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl,  &
       ADM_KNONE,    &
       ADM_kmin,     &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz
    use mod_gmtr, only: &
       GMTR_P_RAREA,  &
       GMTR_T_RAREA,  &
       GMTR_A_HNX,    &
       GMTR_A_HNY,    &
       GMTR_A_HNZ,    &
       GMTR_A_TNX,    &
       GMTR_A_TNY,    &
       GMTR_A_TNZ,    &
       GMTR_A_TN2X,   &
       GMTR_A_TN2Y,   &
       GMTR_A_TN2Z,   &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    use mod_vmtr, only: &
       VMTR_RGAM,       &
       VMTR_RGAM_pl,    &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSH,       &
       VMTR_RGSH_pl,    &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    implicit none

    real(8), intent(out) :: grdx    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdy    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(out) :: grdz    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(out) :: grdz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: sclt         (ADM_gall,   ADM_kall, ADM_lall, ADM_TI:ADM_TJ)  ! for openacc
    real(8) :: sclt_rhogw   (ADM_gall,   ADM_kall, ADM_lall, ADM_TI:ADM_TJ)  !
    real(8) :: sclt_pl      (ADM_gall_pl,ADM_kall)
    real(8) :: sclt_rhogw_pl(ADM_gall_pl,ADM_kall)

    real(8) :: rhogvx_vm   (ADM_gall, ADM_kall, ADM_lall)  ! for openacc
    real(8) :: rhogvy_vm   (ADM_gall, ADM_kall, ADM_lall)  !
    real(8) :: rhogvz_vm   (ADM_gall, ADM_kall, ADM_lall)  !
    real(8) :: rhogw_vm    (ADM_gall, ADM_kall, ADM_lall)  !
    real(8) :: rhogvx_vm_pl(ADM_gall_pl)
    real(8) :: rhogvy_vm_pl(ADM_gall_pl)
    real(8) :: rhogvz_vm_pl(ADM_gall_pl)
    real(8) :: rhogw_vm_pl (ADM_gall_pl,ADM_kall)

    real(8) :: sclt1, sclt2, sclt3, sclt4, sclt5, sclt6

    integer :: rgnid
    integer :: nstart, nend

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: g, k, l, n, v, k0

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,TNX,TNY,TNZ,TN2X,TN2Y,TN2Z,HNX,HNY,HNZ
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT3D_divdamp')

    !$acc data &
    !$acc& pcopy(grdx,grdy,grdz) &
    !$acc& pcopyin(rhogvx,rhogvy,rhogvz,rhogw) &
    !$acc& create(sclt,sclt_rhogw,rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm) &
    !$acc& pcopyin(VMTR_C2Wfact,VMTR_RGAMH,VMTR_RGSH,VMTR_RGAM) &
    !$acc& pcopyin(GRD_rdgz) &
    !$acc& pcopyin(GMTR_T_var,GMTR_P_var,GMTR_A_var) &
    !$acc& pcopyin(ADM_prc_tab,ADM_rgn_vnum)

    k0   = ADM_KNONE
    TI   = ADM_TI
    TJ   = ADM_TJ
    AI   = ADM_AI
    AIJ  = ADM_AIJ
    AJ   = ADM_AJ
    TNX  = GMTR_A_TNX
    TNY  = GMTR_A_TNY
    TNZ  = GMTR_A_TNZ
    HNX  = GMTR_A_HNX
    HNY  = GMTR_A_HNY
    HNZ  = GMTR_A_HNZ
    TN2X = GMTR_A_TN2X
    TN2Y = GMTR_A_TN2Y
    TN2Z = GMTR_A_TN2Z

    !$acc kernels pcopy(rhogw_vm) pcopyin(rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_RGAMH,VMTR_RGSH) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(2)
    do k = ADM_kmin, ADM_kmax+1
    !$acc loop gang vector(128)
    do g = 1, ADM_gall
       if (( k >= ADM_kmin+1 ) .and. ( k <= ADM_kmax )) then
          rhogw_vm(g,k,l) = ( VMTR_C2Wfact(g,1,k,l) * rhogvx(g,k  ,l) &
                            + VMTR_C2Wfact(g,2,k,l) * rhogvx(g,k-1,l) &
                            + VMTR_C2Wfact(g,3,k,l) * rhogvy(g,k  ,l) &
                            + VMTR_C2Wfact(g,4,k,l) * rhogvy(g,k-1,l) &
                            + VMTR_C2Wfact(g,5,k,l) * rhogvz(g,k  ,l) &
                            + VMTR_C2Wfact(g,6,k,l) * rhogvz(g,k-1,l) &
                            ) * VMTR_RGAMH(g,k,l)                     & ! horizontal contribution
                            + rhogw(g,k,l) * VMTR_RGSH(g,k,l)             ! vertical   contribution
       else
          rhogw_vm(g,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(rhogvx_vm,rhogvy_vm,rhogvz_vm) pcopyin(rhogvx,rhogvy,rhogvz,VMTR_RGAM) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(2)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(128)
    do g = 1, ADM_gall
       rhogvx_vm(g,k,l) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
       rhogvy_vm(g,k,l) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
       rhogvz_vm(g,k,l) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(sclt_rhogw) pcopyin(rhogw_vm,GRD_rdgz) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(2)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(128)
    do n = nstart, nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d

       sclt_rhogw(n,k,l,TI) = &
            ( ( rhogw_vm(ij,k+1,l) + rhogw_vm(ip1j,k+1,l) + rhogw_vm(ip1jp1,k+1,l) ) &
            - ( rhogw_vm(ij,k  ,l) + rhogw_vm(ip1j,k  ,l) + rhogw_vm(ip1jp1,k  ,l) ) &
            ) / 3.D0 * GRD_rdgz(k)

       sclt_rhogw(n,k,l,TJ) = &
            ( ( rhogw_vm(ij,k+1,l) + rhogw_vm(ijp1,k+1,l) + rhogw_vm(ip1jp1,k+1,l) ) &
            - ( rhogw_vm(ij,k  ,l) + rhogw_vm(ijp1,k  ,l) + rhogw_vm(ip1jp1,k  ,l) ) &
            ) / 3.D0 * GRD_rdgz(k)
    enddo
    enddo
    enddo
    !$acc end kernels

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(sclt) pcopyin(rhogvx_vm,rhogvy_vm,rhogvz_vm,sclt_rhogw,GMTR_A_var,GMTR_T_var) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(2)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(128)
    do n = nstart, nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d

       sclt(n,k,l,ADM_TI) = &
            ( - ( rhogvx_vm(ij    ,k,l) + rhogvx_vm(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNX) &
              - ( rhogvy_vm(ij    ,k,l) + rhogvy_vm(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNY) &
              - ( rhogvz_vm(ij    ,k,l) + rhogvz_vm(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNZ) &
              - ( rhogvx_vm(ip1j  ,k,l) + rhogvx_vm(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNX) &
              - ( rhogvy_vm(ip1j  ,k,l) + rhogvy_vm(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNY) &
              - ( rhogvz_vm(ip1j  ,k,l) + rhogvz_vm(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNZ) &
              + ( rhogvx_vm(ip1jp1,k,l) + rhogvx_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
              + ( rhogvy_vm(ip1jp1,k,l) + rhogvy_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
              + ( rhogvz_vm(ip1jp1,k,l) + rhogvz_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
              ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TI,GMTR_T_RAREA) &
              + sclt_rhogw(n,k,l,TI)

       sclt(n,k,l,ADM_TJ) = &
            ( - ( rhogvx_vm(ij    ,k,l) + rhogvx_vm(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
              - ( rhogvy_vm(ij    ,k,l) + rhogvy_vm(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
              - ( rhogvz_vm(ij    ,k,l) + rhogvz_vm(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
              + ( rhogvx_vm(ip1jp1,k,l) + rhogvx_vm(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNX) &
              + ( rhogvy_vm(ip1jp1,k,l) + rhogvy_vm(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNY) &
              + ( rhogvz_vm(ip1jp1,k,l) + rhogvz_vm(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNZ) &
              + ( rhogvx_vm(ijp1  ,k,l) + rhogvx_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNX) &
              + ( rhogvy_vm(ijp1  ,k,l) + rhogvy_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNY) &
              + ( rhogvz_vm(ijp1  ,k,l) + rhogvz_vm(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNZ) &
              ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TJ,GMTR_T_RAREA) &
              + sclt_rhogw(n,k,l,TJ)
    enddo
    enddo
    enddo
    !$acc end kernels

    nstart = suf(ADM_gmin,  ADM_gmin  )
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(grdx,grdy,grdz) pcopyin(sclt,GMTR_A_var,GMTR_P_var,ADM_prc_tab,ADM_rgn_vnum) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(2)
    do k = ADM_kmin-1, ADM_kmax+1
    !$acc loop gang vector(128)
    do n = nstart, nend
       if ( ( k >= ADM_kmin ) .and. ( k <= ADM_kmax ) ) then
          ij     = n
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          sclt1 = sclt(ijm1,  k,l,TJ) + sclt(ij,    k,l,TI)
          sclt2 = sclt(ij,    k,l,TI) + sclt(ij,    k,l,TJ)
          sclt3 = sclt(ij,    k,l,TJ) + sclt(im1j,  k,l,TI)
          sclt4 = sclt(im1jm1,k,l,TJ) + sclt(im1j,  k,l,TI)
          sclt5 = sclt(im1jm1,k,l,TI) + sclt(im1jm1,k,l,TJ)
          sclt6 = sclt(ijm1  ,k,l,TJ) + sclt(im1jm1,k,l,TI)
          if ( n == nstart ) then
             rgnid = ADM_prc_tab(l,ADM_prc_me)
             if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then
                sclt5 = sclt(ijm1  ,k,l,TJ) + sclt(im1jm1,k,l,TJ)
                sclt6 = 0.0
             endif
          endif

          grdx(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNX) &
                          + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNX) &
                          + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNX) &
                          - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNX) &
                          - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
                          - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNX) &
                          ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)

          grdy(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNY) &
                          + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNY) &
                          + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNY) &
                          - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNY) &
                          - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
                          - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNY) &
                          ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)

          grdz(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNZ) &
                          + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNZ) &
                          + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNZ) &
                          - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNZ) &
                          - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
                          - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNZ) &
                          ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA)
       else
          grdx(n,k,l) = 0.D0
          grdy(n,k,l) = 0.D0
          grdz(n,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    ! boundary condition
    if ( ADM_prc_me == ADM_prc_pl ) then
       rhogw_vm_pl(:,ADM_kmin  ) = 0.D0
       rhogw_vm_pl(:,ADM_kmax+1) = 0.D0
    endif

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
             do g = 1, ADM_gall_pl
                rhogw_vm_pl(g,k) = ( VMTR_C2Wfact_pl(g,1,k,l) * rhogvx_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(g,2,k,l) * rhogvx_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(g,3,k,l) * rhogvy_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(g,4,k,l) * rhogvy_pl(g,k-1,l) &
                                   + VMTR_C2Wfact_pl(g,5,k,l) * rhogvz_pl(g,k  ,l) &
                                   + VMTR_C2Wfact_pl(g,6,k,l) * rhogvz_pl(g,k-1,l) &
                                   ) * VMTR_RGAMH_pl(g,k,l)                        & ! horizontal contribution
                                 + rhogw_pl(g,k,l) * VMTR_RGSH_pl(g,k,l)             ! vertical   contribution
             enddo
          enddo

          n = ADM_GSLF_PL
          do k = ADM_kmin, ADM_kmax
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                sclt_rhogw_pl(v,k) = ( ( rhogw_vm_pl(n,k+1) + rhogw_vm_pl(ij,k+1) + rhogw_vm_pl(ijp1,k+1) ) &
                                     - ( rhogw_vm_pl(n,k  ) + rhogw_vm_pl(ij,k  ) + rhogw_vm_pl(ijp1,k  ) ) &
                                     ) / 3.D0 * GRD_rdgz(k)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax
             rhogvx_vm_pl(:) = rhogvx_pl(:,k,l) * VMTR_RGAM_pl(:,k,l)
             rhogvy_vm_pl(:) = rhogvy_pl(:,k,l) * VMTR_RGAM_pl(:,k,l)
             rhogvz_vm_pl(:) = rhogvz_pl(:,k,l) * VMTR_RGAM_pl(:,k,l)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                sclt_pl(v,k) = ( + ( rhogvx_vm_pl(n   ) + rhogvx_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNX ) &
                                 + ( rhogvy_vm_pl(n   ) + rhogvy_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNY ) &
                                 + ( rhogvz_vm_pl(n   ) + rhogvz_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNZ ) &
                                 + ( rhogvx_vm_pl(ij  ) + rhogvx_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2X) &
                                 + ( rhogvy_vm_pl(ij  ) + rhogvy_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2Y) &
                                 + ( rhogvz_vm_pl(ij  ) + rhogvz_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2Z) &
                                 - ( rhogvx_vm_pl(ijp1) + rhogvx_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNX ) &
                                 - ( rhogvy_vm_pl(ijp1) + rhogvy_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNY ) &
                                 - ( rhogvz_vm_pl(ijp1) + rhogvz_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNZ ) &
                               ) * 0.5D0 * GMTR_T_var_pl(ij,k0,l,GMTR_T_RAREA) &
                             + sclt_rhogw_pl(v,k)
             enddo
          enddo

          do k = ADM_kmin, ADM_kmax
             grdx_pl(n,k,l) = 0.D0
             grdy_pl(n,k,l) = 0.D0
             grdz_pl(n,k,l) = 0.D0

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 < ADM_gmin_pl ) ijm1 = ADM_gmax_pl ! cyclic condition

                grdx_pl(n,k,l) = grdx_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNX)
                grdy_pl(n,k,l) = grdy_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNY)
                grdz_pl(n,k,l) = grdz_pl(n,k,l) + ( sclt_pl(ijm1,k) + sclt_pl(ij,k) ) * GMTR_A_var_pl(ij,k0,l,HNZ)
             enddo

             grdx_pl(n,k,l) = grdx_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA)
             grdy_pl(n,k,l) = grdy_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA)
             grdz_pl(n,k,l) = grdz_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA)
          enddo

       enddo
    endif

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++OPRT3D_divdamp')

    return
  end subroutine OPRT3D_divdamp

end module mod_oprt3d
!-------------------------------------------------------------------------------
