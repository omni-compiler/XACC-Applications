!-------------------------------------------------------------------------------
!>
!! Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)    Imported from igdc-4.33
!! @li      2006-04-17 (H.Tomita)    Add sub[OPRT_divergence]
!! @li      2006-08-11 (H.Tomita)    Implementatio of miura scheme(2004) with Thurbuen(1996)'s limiter.
!!                                   Add sub[OPRT_divergence2]
!!                                       sub[OPRT_divergence2_prep] 
!!                                       sub[OPRT_divergence2_all].
!! @li      2006-08-22 (Y.Niwa)      divide the rows for the calc. of clap and clap_pl due to the rule of XLF.
!! @li      2007-01-26 (H.Tomita)    Optimization of sub[oprt_diffusion].
!! @li      2007-11-28 (T.Mitsui)    bugfix in oprt_divergence2, _all
!! @li      2008-01-24 (Y.Niwa)      add OPRT_divergence2{_prep,,_all}_rev	
!! @li      2008-04-28 (T.Mitsui)    bug fix in OPRT_divergence2{_all}_rev
!! @li      2009-09-04 (H.Taniguchi) bug fix in OPRT_divergence2{_all}_rev Zero clear of wrk[_pl] is needed.
!! @li      2010-06-08 (S.Iga)       new grid is implemented (see, string XTMS)
!! @li      2011-09-27 (T.Seiki)     merge optimization by RIST and M.Terai 
!! @li      2012-01-17 (M.Terai)     update optimization(case6) in div2rev
!! @li      2012-05-01 (T.Yamaura)   bug fix in div2rev
!! @li      2012-06-28 (M.Terai)     Removed wrapper subroutine to invoked the operators directly macro var.
!!
!<
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT_setup
  public :: OPRT_divergence
  public :: OPRT_gradient
  public :: OPRT_laplacian
  public :: OPRT_diffusion
  public :: OPRT_horizontalize_vec
  public :: OPRT_vorticity
  public :: OPRT_divdamp

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
  integer, private, save :: OPRT_nstart
  integer, private, save :: OPRT_nend

  ! < for divergence operator >
  real(8), public,  allocatable, save :: cdiv   (:,:,:,:)   ! for openacc
  real(8), private, allocatable, save :: cdiv_pl(:,:,:,:)

  ! < for gradient operator >
  real(8), public,  allocatable, save :: cgrad   (:,:,:,:)   ! for openacc
  real(8), private, allocatable, save :: cgrad_pl(:,:,:,:)

  ! < for laplacian operator >
  real(8), public,  allocatable, save :: clap   (:,:,:)    ! for openacc
  real(8), private, allocatable, save :: clap_pl(:,:,:)

  ! < for diffusion operator >   ! for openacc
  real(8), public, allocatable, save :: cmdif_P (:,:)     !(n,l) <- GMTR_P_VAR(n,1,l,P_RAREA)
  real(8), public, allocatable, save :: cmdif_T (:,:,:)   !(TI:TJ,n,l) <- GMTR_T_VAR(n,1,l,TI:TJ,T_RAREA)
  real(8), public, allocatable, save :: cmdif_AH(:,:,:,:) !(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,HN[XYZ])
  real(8), public, allocatable, save :: cmdif_AT(:,:,:,:) !(AI:AJ,1:3,n,l) <-GMTR_A_VAR(n,1,l,TI:TJ,TN[XYZ])

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OPRT_setup
    use mod_adm, only: &
       ADM_W,          &
       ADM_TI,         &
       ADM_TJ,         &
       ADM_AI,         &
       ADM_AIJ,        &
       ADM_AJ,         &
       ADM_prc_tab,    &
       ADM_prc_me,     &
       ADM_prc_pl,     &
       ADM_rgn_vnum,   &
       ADM_vlink_nmax, &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_KNONE,      &
       ADM_kmin,       &
       ADM_kmax
    use mod_gmtr, only: &
       GMTR_P_rarea,  &
       GMTR_T_W1,     &
       GMTR_T_W2,     &
       GMTR_T_W3,     &
       GMTR_T_rarea,  &
       GMTR_A_hnx,    &
       GMTR_A_hny,    &
       GMTR_A_hnz,    &
       GMTR_A_tnx,    &
       GMTR_A_tny,    &
       GMTR_A_tnz,    &
       GMTR_A_tn2x,   &
       GMTR_A_tn2y,   &
       GMTR_A_tn2z,   &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    implicit none

    integer :: n0,n1,n2,n3,n4

    integer :: ij
    integer :: im1j,ijm1,im1jm1
    integer :: ip1j,ijp1,ip1jp1

    integer :: rgnid
    integer :: n, l, m, md, v

    integer :: suf, i, j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,W1,W2,W3
    integer :: k0,a0
    integer :: tx1,ty1,tz1
    integer :: tx2,ty2,tz2
    integer :: hx1,hy1,hz1
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[oprt]/Category[common share]'

    k0  = ADM_KNONE

    TI  = ADM_TI
    TJ  = ADM_TJ
    AI  = ADM_AI
    AIJ = ADM_AIJ
    AJ  = ADM_AJ
    W1  = GMTR_T_W1
    W2  = GMTR_T_W2
    W3  = GMTR_T_W3

    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)

    !---< setup coefficient of divergence operator >
    write(ADM_LOG_FID,*) '*** setup coefficient of divergence operator'

    allocate( cdiv   (ADM_gall   ,ADM_lall   ,0:6,             3) )
    allocate( cdiv_pl(ADM_gall_pl,ADM_lall_pl,0:ADM_vlink_nmax,3) )

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do m = 1, 3
          md = m + GMTR_A_HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             ! ij
             cdiv(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0*GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(n,l,6,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             ! ij
             cdiv(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1j
             cdiv(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cdiv(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijp1
             cdiv(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1j
             cdiv(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! im1jm1
             cdiv(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
             ! ijm1
             cdiv(n,l,6,m) = ( + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5D0 * GMTR_P_var(n,k0,l,GMTR_P_rarea)
          enddo
       endif

    enddo ! loop l

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             cdiv_pl(n,l,0,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cdiv_pl(n,l,0,m) = cdiv_pl(n,l,0,m) + ( GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                                      + GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ijp1,k0,l,md) )
             enddo
             cdiv_pl(n,l,0,m) = cdiv_pl(n,l,0,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cdiv_pl(n,l,v-1,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                       + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                     ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif



    !---< setup coefficient of gradient operator >

    write(ADM_LOG_FID,*) '*** setup coefficient of gradient operator'

    allocate( cgrad   (ADM_gall   ,ADM_lall   ,0:6,             3) )
    allocate( cgrad_pl(ADM_gall_pl,ADM_lall_pl,0:ADM_vlink_nmax,3) )

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do m = 1, 3
          md = m + GMTR_A_HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             ! ij
             cgrad(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.D0 * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.D0 * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                                + 2.D0 * GMTR_A_var(ijm1  ,k0,l,AJ ,md)                          &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cgrad(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cgrad(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cgrad(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cgrad(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cgrad(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cgrad(n,l,6,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + ADM_gall_1d
             ijp1   = n     + ADM_gall_1d
             im1j   = n - 1
             im1jm1 = n - 1 - ADM_gall_1d
             ijm1   = n     - ADM_gall_1d

             ! ij
             cgrad(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.D0 * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.D0 * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.D0 * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1j
             cgrad(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ip1jp1
             cgrad(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijp1
             cgrad(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1j
             cgrad(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! im1jm1
             cgrad(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
             ! ijm1
             cgrad(n,l,6,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_rarea)
          enddo
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + GMTR_A_hnx - 1

             cgrad_pl(n,l,0,m) = 0.D0
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cgrad_pl(n,l,0,m) = cgrad_pl(n,l,0,m) &
                                  + 2.D0 * ( GMTR_T_var_pl(ij,k0,l,W1) - 1.D0 ) * GMTR_A_var_pl(ijp1,k0,l,md)
             enddo
             cgrad_pl(n,l,0,m) = cdiv_pl(n,l,0,m) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cgrad_pl(n,l,v-1,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                        + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                      ) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif

    ! ---- setup coefficient of laplacian operator

    write(ADM_LOG_FID,*) '*** setup coefficient of laplacian operator'

    allocate( clap   (ADM_gall   ,ADM_lall   ,0:6             ) )
    allocate( clap_pl(ADM_gall_pl,ADM_lall_pl,0:ADM_vlink_nmax) )

    a0  = GMTR_T_rarea
    tx1 = GMTR_A_tnx
    ty1 = GMTR_A_tny
    tz1 = GMTR_A_tnz
    hx1 = GMTR_A_hnx
    hy1 = GMTR_A_hny
    hz1 = GMTR_A_hnz

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          ! ij
          clap(ij,l,0) = ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          clap(ij,l,0) = clap(ij,l,0) + ( &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1j
          clap(ij,l,1) = ( &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1jp1
          clap(ij,l,2) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijp1
          clap(ij,l,3) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &     
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1j
          clap(ij,l,4) = ( &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1jm1
          clap(ij,l,5) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &  
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijm1
          clap(ij,l,6) = ( &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TI,a0)*GMTR_A_var(ijm1,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

       enddo

       if ( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          ! 0: ij
          clap(ij,l,0) = ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0  ! Y.Niwa add 06/08/22

          clap(ij,l,0) = clap(ij,l,0) + ( &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1j
          clap(ij,l,1) = ( &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -2*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ij    ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(ijm1  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1  ,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ip1jp1
          clap(ij,l,2) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          -1*GMTR_A_var(ip1j,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijp1
          clap(ij,l,3) = ( &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AIJ,hz1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(ijp1,k0,l,ADM_AI ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ij  ,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          -1*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &     
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1j
          clap(ij,l,4) = ( &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hx1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          +1*GMTR_A_var(im1j,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hy1) &
          -1*GMTR_A_var(im1j,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
          +2*GMTR_A_var(ij  ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j,k0,l,ADM_TI,a0)*GMTR_A_var(ij,k0,l,ADM_AJ,hz1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -2*GMTR_A_var(ij    ,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1j  ,k0,l,ADM_TI,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -2*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! im1jm1
          clap(ij,l,5) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hx1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hy1) &
       +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1j,k0,l,ADM_AI,hz1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(im1jm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(im1j  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(im1jm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

          ! ijm1
          clap(ij,l,6) = ( &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &  
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hx1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hy1) &
      +2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(im1jm1,k0,l,ADM_AIJ,hz1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tx1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hx1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
      -1*GMTR_A_var(ijm1,k0,l,ADM_AIJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      +1*GMTR_A_var(ijm1,k0,l,ADM_AJ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,tz1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hz1) &
      -2*GMTR_A_var(ij  ,k0,l,ADM_AI ,ty1)*GMTR_T_var(ijm1,k0,l,ADM_TJ,a0)*GMTR_A_var(ij,k0,l,ADM_AI,hy1) &
         )*GMTR_P_var(ij,k0,l,GMTR_P_rarea)/12.0d0

      endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then

      n =ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4
      k0=ADM_KNONE
      a0=GMTR_T_rarea
      tx1=gMtr_a_tnx
      tx2=GMTR_A_tn2x
      ty1=GMTR_A_tny
      ty2=GMTR_A_tn2y
      tz1=GMTR_A_tnz
      tz2=GMTR_A_tn2z
      hx1=GMTR_A_hnx
      hy1=GMTR_A_hny
      hz1=GMTR_A_hnz

      do l = 1,ADM_lall_pl
          clap_pl(n,l,0)=( &  
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0   ! Y.Niwa add 060822

          clap_pl(n,l,0)= clap_pl(n,l,0) + ( &                        ! Y.Niwa add 060822 
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n0
          clap_pl(n,l,1)=( &  
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n1
          clap_pl(n,l,2)=( &  
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tx2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n0,k0,l,ty2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n0,k0,l,tz2)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n0,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n2
          clap_pl(n,l,3)=( &  
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n1,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tx2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n1,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n1,k0,l,ty2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n1,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n1,k0,l,tz2)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n1,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n1,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n3
          clap_pl(n,l,4)=( &  
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n2,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tx2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n2,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n2,k0,l,ty2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n2,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n2,k0,l,tz2)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n2,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n2,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) &
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
          !
          ! n4
          clap_pl(n,l,5)=( &  
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n0,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n3,k0,l,hz1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tx2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -2*GMTR_A_var_pl(n3,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       -1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n4,k0,l,tx2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +2*GMTR_A_var_pl(n0,k0,l,tx1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hx1) &
                       +1*GMTR_A_var_pl(n3,k0,l,ty2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -2*GMTR_A_var_pl(n3,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       -1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n4,k0,l,ty2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +2*GMTR_A_var_pl(n0,k0,l,ty1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hy1) &
                       +1*GMTR_A_var_pl(n3,k0,l,tz2)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -2*GMTR_A_var_pl(n3,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       -1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n3,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +1*GMTR_A_var_pl(n4,k0,l,tz2)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                       +2*GMTR_A_var_pl(n0,k0,l,tz1)*GMTR_T_var_pl(n4,k0,l,a0)*GMTR_A_var_pl(n4,k0,l,hz1) & 
                      )*GMTR_P_var_pl(n,k0,l,GMTR_P_rarea)/12.0d0
      enddo
    endif



    allocate( cmdif_P (                  ADM_gall,ADM_lall) )
    allocate( cmdif_T (    ADM_gall,ADM_lall,ADM_TI:ADM_TJ) )    ! cnvtd cmdif_T
    allocate( cmdif_AT(ADM_gall,ADM_lall,ADM_AI:ADM_AJ,1:3) )    ! cnvtd cmdif_AT
    allocate( cmdif_AH(ADM_gall,ADM_lall,ADM_AI:ADM_AJ,1:3) )    ! cnvtd cmdif_AH

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cmdif_P(n,l) = GMTR_P_var(n,k0,l,GMTR_P_RAREA)
    enddo
    enddo

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cmdif_T(n,l,TI) = GMTR_T_var(n,k0,l,TI,GMTR_T_RAREA)    ! cnvtd cmdif_T
       cmdif_T(n,l,TJ) = GMTR_T_var(n,k0,l,TJ,GMTR_T_RAREA)    ! cnvtd cmdif_T
    enddo
    enddo

    do m=ADM_AI,ADM_AJ
       do l=1,ADM_lall
          do n=1,ADM_gall
             cmdif_AT( n, l,m, 1 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNX)    ! cnvtd cmdif_AT
             cmdif_AT( n, l,m, 2 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNY)    ! cnvtd cmdif_AT
             cmdif_AT( n, l,m, 3 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_TNZ)    ! cnvtd cmdif_AT

             cmdif_AH( n, l,m, 1 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNX)    ! cnvtd cmdif_AH
             cmdif_AH( n, l,m, 2 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNY)    ! cnvtd cmdif_AH
             cmdif_AH( n, l,m, 3 ) =  GMTR_A_var(n, ADM_KNONE, l, m, GMTR_A_HNZ)    ! cnvtd cmdif_AH
          enddo
       enddo
    enddo

    return
  end subroutine OPRT_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence( &
       scl, scl_pl, &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       mfact        )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(inout) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divergence')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    !$acc kernels pcopy(scl) &
    !$acc& pcopyin(vx,vy,vz,cdiv) 
    !$acc loop gang vector(1)
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(64)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ip1jp1 = n + 1 + ADM_gall_1d
       ijp1   = n     + ADM_gall_1d
       im1j   = n - 1
       im1jm1 = n - 1 - ADM_gall_1d
       ijm1   = n     - ADM_gall_1d
       scl(n,k,l) = ( cdiv(n,l,0,1) * vx(ij    ,k,l) &
                    + cdiv(n,l,1,1) * vx(ip1j  ,k,l) &
                    + cdiv(n,l,2,1) * vx(ip1jp1,k,l) &
                    + cdiv(n,l,3,1) * vx(ijp1  ,k,l) &
                    + cdiv(n,l,4,1) * vx(im1j  ,k,l) &
                    + cdiv(n,l,5,1) * vx(im1jm1,k,l) &
                    + cdiv(n,l,6,1) * vx(ijm1  ,k,l) &
                    + cdiv(n,l,0,2) * vy(ij    ,k,l) &
                    + cdiv(n,l,1,2) * vy(ip1j  ,k,l) &
                    + cdiv(n,l,2,2) * vy(ip1jp1,k,l) &
                    + cdiv(n,l,3,2) * vy(ijp1  ,k,l) &
                    + cdiv(n,l,4,2) * vy(im1j  ,k,l) &
                    + cdiv(n,l,5,2) * vy(im1jm1,k,l) &
                    + cdiv(n,l,6,2) * vy(ijm1  ,k,l) &
                    + cdiv(n,l,0,3) * vz(ij    ,k,l) &
                    + cdiv(n,l,1,3) * vz(ip1j  ,k,l) &
                    + cdiv(n,l,2,3) * vz(ip1jp1,k,l) &
                    + cdiv(n,l,3,3) * vz(ijp1  ,k,l) &
                    + cdiv(n,l,4,3) * vz(im1j  ,k,l) &
                    + cdiv(n,l,5,3) * vz(im1jm1,k,l) &
                    + cdiv(n,l,6,3) * vz(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          scl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( cdiv_pl(n,l,v-1,1) * vx_pl(v,k,l) &
                                             + cdiv_pl(n,l,v-1,2) * vy_pl(v,k,l) &
                                             + cdiv_pl(n,l,v-1,3) * vz_pl(v,k,l) )
          enddo

          scl_pl(n,k,l) = scl_pl(n,k,l) * fact
       enddo
       enddo
    endif

    !--- !$acc wait

    call DEBUG_rapend('++++OPRT_divergence')

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient( &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl,  &
       scl, scl_pl, &
       mfact        )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(in)    :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_gradient')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    !$acc kernels pcopy(vx,vy,vz) pcopyin(cgrad,scl) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = 1, ADM_kall
    !$acc loop gang vector(64)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ip1jp1 = n + 1 + ADM_gall_1d
       ijp1   = n     + ADM_gall_1d
       im1j   = n - 1
       im1jm1 = n - 1 - ADM_gall_1d
       ijm1   = n     - ADM_gall_1d

       vx(n,k,l) = ( cgrad(ij,l,0,1) * scl(ij    ,k,l) &
                   + cgrad(ij,l,1,1) * scl(ip1j  ,k,l) &
                   + cgrad(ij,l,2,1) * scl(ip1jp1,k,l) &
                   + cgrad(ij,l,3,1) * scl(ijp1  ,k,l) &
                   + cgrad(ij,l,4,1) * scl(im1j  ,k,l) &
                   + cgrad(ij,l,5,1) * scl(im1jm1,k,l) &
                   + cgrad(ij,l,6,1) * scl(ijm1  ,k,l) ) * fact

       vy(n,k,l) = ( cgrad(ij,l,0,2) * scl(ij    ,k,l) &
                   + cgrad(ij,l,1,2) * scl(ip1j  ,k,l) &
                   + cgrad(ij,l,2,2) * scl(ip1jp1,k,l) &
                   + cgrad(ij,l,3,2) * scl(ijp1  ,k,l) &
                   + cgrad(ij,l,4,2) * scl(im1j  ,k,l) &
                   + cgrad(ij,l,5,2) * scl(im1jm1,k,l) &
                   + cgrad(ij,l,6,2) * scl(ijm1  ,k,l) ) * fact

       vz(n,k,l) = ( cgrad(ij,l,0,3) * scl(ij    ,k,l) &
                   + cgrad(ij,l,1,3) * scl(ip1j  ,k,l) &
                   + cgrad(ij,l,2,3) * scl(ip1jp1,k,l) &
                   + cgrad(ij,l,3,3) * scl(ijp1  ,k,l) &
                   + cgrad(ij,l,4,3) * scl(im1j  ,k,l) &
                   + cgrad(ij,l,5,3) * scl(im1jm1,k,l) &
                   + cgrad(ij,l,6,3) * scl(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          vx_pl(n,k,l) = 0.D0
          vy_pl(n,k,l) = 0.D0
          vz_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             vx_pl(n,k,l) = vx_pl(n,k,l) + cgrad_pl(n,l,v-1,1) * scl_pl(v,k,l)
             vy_pl(n,k,l) = vy_pl(n,k,l) + cgrad_pl(n,l,v-1,2) * scl_pl(v,k,l)
             vz_pl(n,k,l) = vz_pl(n,k,l) + cgrad_pl(n,l,v-1,3) * scl_pl(v,k,l)
          enddo

          vx_pl(n,k,l) = vx_pl(n,k,l) * fact
          vy_pl(n,k,l) = vy_pl(n,k,l) * fact
          vz_pl(n,k,l) = vz_pl(n,k,l) * fact
       enddo
       enddo
    endif

    !--- !$acc wait

    call DEBUG_rapend('++++OPRT_gradient')

    return
  end subroutine OPRT_gradient

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       mfact          )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    real(8), intent(inout) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: fact

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_laplacian')

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    !$acc kernels pcopy(dscl) pcopyin(clap,scl) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = 1, ADM_kall
    !$acc loop gang vector(64)
    do n = OPRT_nstart, OPRT_nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d
       im1j   = n - 1
       ijm1   = n     - ADM_gall_1d
       im1jm1 = n - 1 - ADM_gall_1d

       dscl(n,k,l) = ( clap(ij,l,0) * scl(ij    ,k,l) &
                     + clap(ij,l,1) * scl(ip1j  ,k,l) &
                     + clap(ij,l,2) * scl(ip1jp1,k,l) &
                     + clap(ij,l,3) * scl(ijp1  ,k,l) &
                     + clap(ij,l,4) * scl(im1j  ,k,l) &
                     + clap(ij,l,5) * scl(im1jm1,k,l) &
                     + clap(ij,l,6) * scl(ijm1  ,k,l) ) * fact
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          dscl_pl(n,k,l) = 0.D0

          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + clap_pl(n,l,v-1) * scl_pl(v,k,l)
          enddo

          dscl_pl(n,k,l) = dscl_pl(n,k,l) * fact
       enddo
       enddo
    endif

    !--- !$acc wait

    call DEBUG_rapend('++++OPRT_laplacian')

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------  
  subroutine OPRT_diffusion( &
       dscl, dscl_pl,        &
       scl,scl_pl,           &
       kh, kh_pl,            &
       mfact )
    use mod_adm, only :   &
         !--- public parameters
         ADM_prc_pl,      &
         ADM_KNONE,       &
         ADM_TI,          &
         ADM_TJ,          &
         adm_ai,          &
         adm_aij,         &
         adm_aj,          &
         adm_w,           &
         ADM_lall_pl,     &
         ADM_GMAX_PL,     &
         ADM_gmin_pl,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         ADM_rgn_vnum,    &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
    use mod_gmtr, only : &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_T_rarea,    &
         GMTR_A_hnx,      &
         GMTR_A_hny,      &
         GMTR_A_hnz,      &
         GMTR_A_htx,      &
         GMTR_A_hty,      &
         GMTR_A_htz,      &
         GMTR_A_tnx,      &
         GMTR_A_tny,      &
         GMTR_A_tnz,      &
         GMTR_A_tn2x,     &
         GMTR_A_tn2y,     &
         GMTR_A_tn2z,     &
         GMTR_P_rarea,    &
         !--- public variables
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    implicit none

    real(8), intent(inout) :: dscl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in) :: kh(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in) :: kh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in), optional :: mfact
    !
    real(8)  :: vxt(ADM_gall, ADM_kall, ADM_lall, ADM_TI:ADM_TJ)
    real(8)  :: vyt(ADM_gall, ADM_kall, ADM_lall, ADM_TI:ADM_TJ)
    real(8)  :: vzt(ADM_gall, ADM_kall, ADM_lall, ADM_TI:ADM_TJ)
    real(8)  :: flux(ADM_gall,ADM_kall, ADM_lall, ADM_AI:ADM_AJ)    ! cnvtd flux

    real(8)  :: vxt_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl)
    real(8)  :: vyt_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl)
    real(8)  :: vzt_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl)
    real(8)  :: flux_pl(ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: rgnid
    real(8) :: fact
    real(8) :: u1,u2,u3
    real(8) :: smean
    !
    integer :: nstart, nend
    integer :: nstart1, nstart2, nstart3 !2011/02/25 for K by RIST 
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_diffusion')

    !$acc data &
    !$acc& pcopy(dscl) &
    !$acc& pcopyin(scl,kh) &
    !$acc& create(vxt,vyt,vzt,flux) &
    !$acc& pcopyin(cmdif_T,cmdif_AT,cmdif_AH,cmdif_P) &
    !$acc& pcopyin(ADM_prc_tab,ADM_rgn_vnum) &
    !$acc& pcopyin(GMTR_P_var)

    if ( present(mfact) ) then
       fact = mfact
    else
       fact=1.0D0
    end if

    !
    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax  ,ADM_gmax  )
    !$acc kernels pcopy(vxt,vyt,vzt) pcopyin(scl, cmdif_T,cmdif_AT) 
    do l=1,ADM_lall
       do k= 1, ADM_kall
          do n = nstart, nend

             smean = (scl(n,k,l)+scl(n+1,k,l)+scl(n+1+ADM_gall_1d,k,l))/3.0D0
             u1=-cmdif_T( n, l,ADM_TI) &
                    *(0.5D0*(scl(n,k,l)+scl(n+1,k,l))-smean)
             u2=-cmdif_T( n, l,ADM_TI) &
                    *(0.5D0*(scl(n+1,k,l)+scl(n+1+ADM_gall_1d,k,l))-smean)
             u3=+cmdif_T( n, l,ADM_TI) &
                    *(0.5D0*(scl(n+1+ADM_gall_1d,k,l)+scl(n,k,l))-smean)
             vxt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT( n  , l,ADM_AI , 1 )&
                  +u2*cmdif_AT( n+1, l,ADM_AJ , 1 )&
                  +u3*cmdif_AT( n  , l,ADM_AIJ, 1 ) 
             vyt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT( n  , l,ADM_AI , 2 )&
                  +u2*cmdif_AT( n+1, l,ADM_AJ , 2 )&
                  +u3*cmdif_AT( n  , l,ADM_AIJ, 2 ) 
             vzt( n, k, l, ADM_TI ) &
                  =u1*cmdif_AT( n  , l,ADM_AI , 3 )&
                  +u2*cmdif_AT( n+1, l,ADM_AJ , 3 )&
                  +u3*cmdif_AT( n  , l,ADM_AIJ, 3 ) 

             smean = (scl(n,k,l)+scl(n+1+ADM_gall_1d,k,l)+scl(n+ADM_gall_1d,k,l))/3.0D0
             u1=-cmdif_T( n, l,ADM_TJ) &
                    *(0.5D0*(scl(n              ,k,l)+scl(n+1+ADM_gall_1d,k,l))-smean)
             u2=+cmdif_T( n, l,ADM_TJ) &
                    *(0.5D0*(scl(n+1+ADM_gall_1d,k,l)+scl(n  +ADM_gall_1d,k,l))-smean)
             u3=+cmdif_T( n, l,ADM_TJ) &
                    *(0.5D0*(scl(n  +ADM_gall_1d,k,l)+scl(n              ,k,l))-smean)
             vxt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT( n            , l,ADM_AIJ, 1 )&
                  +u2*cmdif_AT( n+ADM_gall_1d, l,ADM_AI , 1 )&
                  +u3*cmdif_AT( n            , l,ADM_AJ , 1 ) 
             vyt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT( n            , l,ADM_AIJ, 2 )&
                  +u2*cmdif_AT( n+ADM_gall_1d, l,ADM_AI , 2 )&
                  +u3*cmdif_AT( n            , l,ADM_AJ , 2 ) 
             vzt( n, k, l, ADM_TJ ) &
                  =u1*cmdif_AT( n            , l,ADM_AIJ, 3 )&
                  +u2*cmdif_AT( n+ADM_gall_1d, l,ADM_AI , 3 )&
                  +u3*cmdif_AT( n            , l,ADM_AJ , 3 ) 
             !
          end do
       enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(vxt,vyt,vzt) pcopyin(ADM_prc_tab,ADM_rgn_vnum) 
    do l=1,ADM_lall
       do k= 1, ADM_kall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI) = vxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI) = vyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI) = vzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
          !
       end do
    end do
    !$acc end kernels

    if(ADM_prc_me==ADM_prc_pl) then
       do k= 1, ADM_kall
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                smean = (scl_pl(ADM_gslf_pl,k,l)+scl_pl(n,k,l)+scl_pl(n+1,k,l))/3.0D0
                u1=+GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(ADM_gslf_pl,k,l)+scl_pl(n,k,l))-smean)
                u2=+GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(n,k,l)+scl_pl(n+1,k,l))-smean)
                u3=-GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_RAREA)&
                     *(0.5D0*(scl_pl(n+1,k,l)+scl_pl(ADM_gslf_pl,k,l))-smean)
                vxt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNX)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2X)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNX)
                vyt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNY)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Y)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNY)
                vzt_pl(n,k,l)&
                     =u1*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TNZ)&
                     +u2*GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_TN2Z)&
                     +u3*GMTR_A_var_pl(n+1,ADM_KNONE,l,GMTR_A_TNZ)
             enddo
             smean = (scl_pl(ADM_gslf_pl,k,l)+scl_pl(ADM_GMAX_PL,k,l)+scl_pl(ADM_gmin_pl,k,l))/3.0D0
             u1=+GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_gslf_pl,k,l)+scl_pl(ADM_GMAX_PL,k,l))-smean)
             u2=+GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_GMAX_PL,k,l)+scl_pl(ADM_gmin_pl,k,l))-smean)
             u3=-GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_RAREA)&
                  *(0.5D0*(scl_pl(ADM_gmin_pl,k,l)+scl_pl(ADM_gslf_pl,k,l))-smean)
             vxt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNX)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2X)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNX)
             vyt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNY)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Y)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNY)
             vzt_pl(ADM_GMAX_PL,k,l)&
                  =u1*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TNZ)&
                  +u2*GMTR_A_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_A_TN2Z)&
                  +u3*GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_TNZ)
          enddo
       enddo
    end if

    !
    nstart1 = suf(ADM_gmin-1,ADM_gmin  )
    nstart2 = suf(ADM_gmin-1,ADM_gmin-1)
    nstart3 = suf(ADM_gmin  ,ADM_gmin-1)
    nstart = min( nstart1, nstart2, nstart3 )
    nend   = suf(ADM_gmax  ,ADM_gmax  )
    !$acc kernels pcopy(flux) pcopyin(vxt,vyt,vzt,kh, cmdif_AH) 
    do l=1,ADM_lall
       do k= 1, ADM_kall
          do n = nstart,nend
             if ( n >= nstart1 ) then
                flux( n,k,l, ADM_AI )  =(&
                     &(vxt( n-ADM_gall_1d, k, l, ADM_TJ ) +vxt( n, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AI, 1 )&
                     +(vyt( n-ADM_gall_1d, k, l, ADM_TJ ) +vyt( n, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AI, 2 )&
                     +(vzt( n-ADM_gall_1d, k, l, ADM_TJ ) +vzt( n, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AI, 3 )&
                     &)*0.5D0
                flux( n,k,l, ADM_AI )  = flux( n,k,l, ADM_AI )  * (kh(n,k,l)+kh(n+1,k,l))*0.5D0
             endif

             if ( n >= nstart2 ) then
                flux( n,k,l, ADM_AIJ )  =(&
                     &(vxt( n, k, l, ADM_TI ) +vxt( n, k, l, ADM_TJ ) ) &
                     *cmdif_AH( n, l,ADM_AIJ, 1 )&
                     +(vyt( n, k, l, ADM_TI ) +vyt( n, k, l, ADM_TJ ) ) &
                     *cmdif_AH( n, l,ADM_AIJ, 2 )&
                     +(vzt( n, k, l, ADM_TI ) +vzt( n, k, l, ADM_TJ ) ) &
                     *cmdif_AH( n, l,ADM_AIJ, 3 )&
                     &)*0.5D0
                flux( n,k,l, ADM_AIJ )  = flux( n,k,l, ADM_AIJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d+1,k,l))*0.5D0
             endif

             if ( n >= nstart3 ) then
                flux( n,k,l, ADM_AJ )  =(&
                     &(vxt( n, k, l, ADM_TJ ) +vxt( n-1, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AJ, 1 )&
                     +(vyt( n, k, l, ADM_TJ ) +vyt( n-1, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AJ, 2 )&
                     +(vzt( n, k, l, ADM_TJ ) +vzt( n-1, k, l, ADM_TI ) ) &
                     *cmdif_AH( n, l,ADM_AJ, 3 )&
                     &)*0.5D0
                flux( n,k,l, ADM_AJ )  = flux( n,k,l, ADM_AJ )  * (kh(n,k,l)+kh(n+ADM_gall_1d,k,l))*0.5D0
             endif
          end do
       enddo
    enddo
    !$acc end kernels

    !
    nstart = suf(ADM_gmin  ,ADM_gmin  )
    nend   = suf(ADM_gmax  ,ADM_gmax  )
    !$acc kernels pcopy(dscl) pcopyin(flux, cmdif_P) 
    do l=1,ADM_lall
       do k= 1, ADM_kall
          do n = nstart,nend
             dscl(n,k,l)=(&
                  +flux( n              ,k,l, ADM_AI  ) &
                  +flux( n              ,k,l, ADM_AIJ ) &
                  +flux( n              ,k,l, ADM_AJ  ) &
                  -flux( n-1            ,k,l, ADM_AI  ) &
                  -flux( n-1-ADM_gall_1d,k,l, ADM_AIJ ) &
                  -flux( n  -ADM_gall_1d,k,l, ADM_AJ ) & 
                  ) * cmdif_P( n, l)&
                  * fact
          end do
       enddo
    enddo
    !$acc end kernels

    !
    !$acc kernels pcopy(dscl) pcopyin(flux, GMTR_P_var, ADM_prc_tab,ADM_rgn_vnum) 
    do l=1,ADM_lall
       do k= 1, ADM_kall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             dscl(suf(ADM_gmin,ADM_gmin),k,l)=(&
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AI)& 
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AIJ)&
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AJ)& 
                  -flux(suf(ADM_gmin-1,ADM_gmin),k,l,ADM_AI)&
                  -flux(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_AIJ)&
                  ) * GMTR_P_var(suf(ADM_gmin,ADM_gmin),ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end if
       end do
    end do
    !$acc end kernels

    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k= 1, ADM_kall
          do l=1,ADM_lall_pl
             flux_pl(ADM_gmin_pl)&
                  =((vxt_pl(ADM_GMAX_PL,k,l)+vxt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNX)&
                  +(vyt_pl(ADM_GMAX_PL,k,l)+vyt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNY)&
                  +(vzt_pl(ADM_GMAX_PL,k,l)+vzt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HNZ) )*0.5D0
             flux_pl(ADM_gmin_pl) = flux_pl(ADM_gmin_pl)&
                  * ( kh_pl(ADM_gslf_pl,k,l)+kh_pl(ADM_gmin_pl,k,l) )*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =((vxt_pl(n-1,k,l)+vxt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNX)&
                     +(vyt_pl(n-1,k,l)+vyt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNY)&
                     +(vzt_pl(n-1,k,l)+vzt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HNZ) )*0.5D0
                flux_pl(n) = flux_pl(n)&
                     * ( kh_pl(ADM_gslf_pl,k,l)+kh_pl(n,k,l) )*0.5D0
             enddo
             !
             dscl_pl(ADM_gslf_pl,k,l)=(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
       enddo
    end if

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++OPRT_diffusion')

  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------  
  subroutine OPRT_horizontalize_vec( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl,  &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_e,    &
       GRD_e_pl
    implicit none

    real(8), intent(inout) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(8), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(8) :: prd
    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_horizontalize_vec')

    !$acc kernels pcopy(vx,vy,vz) pcopyin(GRD_e) 
    !$acc loop gang
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = 1, ADM_kall
    !$acc loop gang vector(64)
    do g = 1, ADM_gall
       prd = ( vx(g,k,l) * GRD_e(g,l,GRD_XDIR) &
             + vy(g,k,l) * GRD_e(g,l,GRD_YDIR) &
             + vz(g,k,l) * GRD_e(g,l,GRD_ZDIR) )

       vx(g,k,l) = vx(g,k,l) - prd * GRD_e(g,l,GRD_XDIR)
       vy(g,k,l) = vy(g,k,l) - prd * GRD_e(g,l,GRD_YDIR)
       vz(g,k,l) = vz(g,k,l) - prd * GRD_e(g,l,GRD_ZDIR)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prd = ( vx_pl(g,k,l) * GRD_e_pl(g,l,GRD_XDIR) &
                + vy_pl(g,k,l) * GRD_e_pl(g,l,GRD_YDIR) &
                + vz_pl(g,k,l) * GRD_e_pl(g,l,GRD_ZDIR) )

          vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_XDIR)
          vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_YDIR)
          vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_e_pl(g,l,GRD_ZDIR)
       enddo
       enddo
       enddo
    endif

    !--- !$acc wait

    call DEBUG_rapend('++++OPRT_horizontalize_vec')

    return
  end subroutine OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  subroutine OPRT_vorticity(       &
       scl, scl_pl,                &
       vx, vx_pl,                  &
       vy, vy_pl,                  &
       vz, vz_pl,                  &
       mfact )
    !
    use mod_adm, only :   &
         !--- public parameters
         ADM_W,           &
         ADM_prc_pl,      &
         ADM_TI,          &
         ADM_TJ,          &
         ADM_AI,          &
         ADM_AIJ,         &
         ADM_AJ,          &
         ADM_KNONE,       &
         ADM_lall_pl,     &
         ADM_gmin_pl,     &
         ADM_GMAX_PL,     &
         ADM_gslf_pl,     &
         ADM_gall_pl,     &
         !--- public variables
         ADM_prc_me,      &
         ADM_prc_tab,     &
         ADM_rgn_vnum,    &
         ADM_gall_1d,     &
         ADM_kall,        &
         ADM_lall,        &
         ADM_kmin,        &
         ADM_kmax,        &
         ADM_gmin,        &
         ADM_gmax,        &
         ADM_gall
    use mod_gmtr, only :  &
         !--- public parameters
         GMTR_T_W1,       &
         GMTR_T_W2,       &
         GMTR_T_W3,       &
         GMTR_A_htx,      &
         GMTR_A_hty,      &
         GMTR_A_htz,      &
         GMTR_P_rarea,    &
         !--- public variables
         GMTR_T_var,      &
         GMTR_T_var_pl,   &
         GMTR_P_var,      &
         GMTR_P_var_pl,   &
         GMTR_A_var,      &
         GMTR_A_var_pl
    !
    implicit none
    !
    real(8), intent(inout) :: scl(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in)  :: vx(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vy(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)  :: vz(ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(8), intent(in), optional :: mfact
    !
    real(8)  :: vxt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vxt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vyt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vyt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: vzt(         &
         ADM_gall,           &
         ADM_kall,           &
         ADM_lall,           &
         ADM_TI:ADM_TJ)
    real(8)  :: vzt_pl(      &
         ADM_gall_pl,        &
         ADM_kall,           &
         ADM_lall_pl)
    !
    real(8)  :: flux(ADM_gall, ADM_kall, ADM_lall, ADM_AI:ADM_AJ)    ! for openacc
    real(8)  :: flux_pl(     &
         ADM_gall_pl)
    !
    integer :: l,n,k
    integer :: rgnid
    real(8) :: fact
    !
    integer :: nstart,nend
    !
    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    call DEBUG_rapstart('++++OPRT_vorticity')

    if(present(mfact)) then
       fact=mfact
    else
       fact=1.0D0
    end if

    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       do k=ADM_kmin,ADM_kmax
          !
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend = suf(ADM_gmax,ADM_gmax)
          do n = nstart,nend
             vxt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vx(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vx(n+1+ADM_gall_1d,k,l)
             vxt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vx(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vx(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vx(n+ADM_gall_1d,k,l)
             vyt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vy(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vy(n+1+ADM_gall_1d,k,l)
             vyt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vy(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vy(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vy(n+ADM_gall_1d,k,l)
             vzt(n,k,l,ADM_TI)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W1)&
                  *vz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W2)&
                  *vz(n+1,k,l)                       &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TI,GMTR_T_W3)&
                  *vz(n+1+ADM_gall_1d,k,l)
             vzt(n,k,l,ADM_TJ)                            &
                  =GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W1)&
                  *vz(n,k,l)                         &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W2)&
                  *vz(n+1+ADM_gall_1d,k,l)                     &
                  +GMTR_T_var(n,ADM_KNONE,l,ADM_TJ,GMTR_T_W3)&
                  *vz(n+ADM_gall_1d,k,l)
          enddo
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vxt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vyt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_TI)     &
                  =vzt(suf(ADM_gmin,ADM_gmin-1),k,l,ADM_TJ)
          end if
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       !
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             do n=ADM_gmin_pl,ADM_GMAX_PL-1
                vxt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vx_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vx_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vx_pl(n+1,k,l)
                vyt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vy_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vy_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vy_pl(n+1,k,l)
                vzt_pl(n,k,l)                            &
                     =GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W1)&
                     *vz_pl(ADM_gslf_pl,k,l)              &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W2)&
                     *vz_pl(n,k,l)                        &
                     +GMTR_T_var_pl(n,ADM_KNONE,l,GMTR_T_W3)&
                     *vz_pl(n+1,k,l)
             enddo
             vxt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vx_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vx_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vx_pl(ADM_gmin_pl,k,l)
             vyt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vy_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vy_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vy_pl(ADM_gmin_pl,k,l)
             vzt_pl(ADM_GMAX_PL,k,l)&
                  =GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W1)&
                  *vz_pl(ADM_gslf_pl,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W2)&
                  *vz_pl(ADM_GMAX_PL,k,l)                        &
                  +GMTR_T_var_pl(ADM_GMAX_PL,ADM_KNONE,l,GMTR_T_W3)&
                  *vz_pl(ADM_gmin_pl,k,l)
          enddo
       enddo
    end if
    !
    do l=1,ADM_lall
       rgnid=ADM_prc_tab(l,ADM_prc_me)
       !
       do k=ADM_kmin,ADM_kmax
          !
          nstart = suf(ADM_gmin-1,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,k,l,ADM_AI)&    ! cnvtd flux
                  =((vxt(n-ADM_gall_1d,k,l,ADM_TJ)+vxt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTX)&
                  +(vyt(n-ADM_gall_1d,k,l,ADM_TJ)+vyt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTY)&
                  +(vzt(n-ADM_gall_1d,k,l,ADM_TJ)+vzt(n,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AI,GMTR_A_HTZ))*0.5D0
          enddo
          nstart = suf(ADM_gmin-1,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,k,l,ADM_AIJ)&    ! cnvtd flux
                  =((vxt(n,k,l,ADM_TI)+vxt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTX)&
                  +(vyt(n,k,l,ADM_TI)+vyt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTY)&
                  +(vzt(n,k,l,ADM_TI)+vzt(n,k,l,ADM_TJ))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AIJ,GMTR_A_HTZ))*0.5D0
          enddo
          nstart = suf(ADM_gmin  ,ADM_gmin-1)
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             flux(n,k,l,ADM_AJ)&    ! cnvtd flux
                  =((vxt(n,k,l,ADM_TJ)+vxt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTX)&
                  +(vyt(n,k,l,ADM_TJ)+vyt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTY)&
                  +(vzt(n,k,l,ADM_TJ)+vzt(n-1,k,l,ADM_TI))&
                  *GMTR_A_var(n,ADM_KNONE,l,ADM_AJ,GMTR_A_HTZ))*0.5D0
          enddo
          !
          nstart = suf(ADM_gmin  ,ADM_gmin  )
          nend   = suf(ADM_gmax  ,ADM_gmax  )
          do n = nstart,nend
             scl(n,k,l)=-(&
                  +flux(n,k,l,ADM_AI)&    ! cnvtd flux
                  +flux(n,k,l,ADM_AIJ)&    ! cnvtd flux
                  +flux(n,k,l,ADM_AJ)&    ! cnvtd flux
                  -flux(n-1,k,l,ADM_AI)&    ! cnvtd flux
                  -flux(n-1-ADM_gall_1d,k,l,ADM_AIJ)&    ! cnvtd flux
                  -flux(n-ADM_gall_1d,k,l,ADM_AJ)&    ! cnvtd flux
                  ) * GMTR_P_var(n,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
          !
          if(ADM_rgn_vnum(ADM_W,rgnid)==3) then
             scl(suf(ADM_gmin,ADM_gmin),k,l)=-(&
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AI)&    ! cnvtd flux
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AIJ)&    ! cnvtd flux
                  +flux(suf(ADM_gmin,ADM_gmin),k,l,ADM_AJ)&    ! cnvtd flux
                  -flux(suf(ADM_gmin-1,ADM_gmin),k,l,ADM_AI)&    ! cnvtd flux
                  -flux(suf(ADM_gmin-1,ADM_gmin-1),k,l,ADM_AIJ)&    ! cnvtd flux
                  ) * GMTR_P_var(suf(ADM_gmin,ADM_gmin),ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          end if
       enddo
    enddo
    !
    if(ADM_prc_me==ADM_prc_pl) then
       do k=ADM_kmin,ADM_kmax
          do l=1,ADM_lall_pl
             !
             flux_pl(ADM_gmin_pl)&
                  =((vxt_pl(ADM_GMAX_PL,k,l)+vxt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTX)&
                  +(vyt_pl(ADM_GMAX_PL,k,l)+vyt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTY)&
                  +(vzt_pl(ADM_GMAX_PL,k,l)+vzt_pl(ADM_gmin_pl,k,l))&
                  *GMTR_A_var_pl(ADM_gmin_pl,ADM_KNONE,l,GMTR_A_HTZ) )*0.5D0
             do n=ADM_gmin_pl+1,ADM_GMAX_PL
                flux_pl(n)&
                     =((vxt_pl(n-1,k,l)+vxt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTX)&
                     +(vyt_pl(n-1,k,l)+vyt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTY)&
                     +(vzt_pl(n-1,k,l)+vzt_pl(n,k,l))&
                     *GMTR_A_var_pl(n,ADM_KNONE,l,GMTR_A_HTZ) )*0.5D0
             enddo
             !
             scl_pl(ADM_gslf_pl,k,l)=-(&
                  +flux_pl(ADM_gmin_pl  )&
                  +flux_pl(ADM_gmin_pl+1)&
                  +flux_pl(ADM_gmin_pl+2)&
                  +flux_pl(ADM_gmin_pl+3)&
                  +flux_pl(ADM_gmin_pl+4)&
                  ) * GMTR_P_var_pl(ADM_gslf_pl,ADM_KNONE,l,GMTR_P_RAREA)&
                  * fact
          enddo
       enddo
    end if

    call DEBUG_rapend('++++OPRT_vorticity')

  end subroutine OPRT_vorticity

  !-----------------------------------------------------------------------------
  subroutine OPRT_divdamp( &
       grdx, grdx_pl, &
       grdy, grdy_pl, &
       grdz, grdz_pl, &
       vx,   vx_pl,   &
       vy,   vy_pl,   &
       vz,   vz_pl,   &
       mfact          )
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
    implicit none

    real(8), intent(inout) :: grdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(inout) :: grdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(inout) :: grdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in)    :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(8), intent(in)    :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(8), intent(in), optional :: mfact

    real(8) :: sclt   (ADM_gall,   ADM_kall, ADM_lall, ADM_TI:ADM_TJ)  ! for openacc
    real(8) :: sclt_pl(ADM_gall_pl,ADM_kall)
    real(8) :: fact

    real(8) :: sclt1, sclt2, sclt3, sclt4, sclt5, sclt6

    integer :: rgnid
    integer :: nstart, nend

    integer :: ij
    integer :: im1j, ijm1, im1jm1
    integer :: ip1j, ijp1, ip1jp1

    integer :: k, l, n, v, k0

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)

    integer :: TI,TJ,AI,AIJ,AJ,TNX,TNY,TNZ,TN2X,TN2Y,TN2Z,HNX,HNY,HNZ
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('++++OPRT_divdamp')

    !$acc data &
    !$acc& pcopy(grdx,grdy,grdz) &
    !$acc& pcopyin(vx,vy,vz) &
    !$acc& create(sclt) &
    !$acc& pcopyin(GMTR_A_var,GMTR_T_var,GMTR_P_var) &
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

    if ( present(mfact) ) then
       fact = mfact
    else
       fact = 1.D0
    endif

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(sclt) pcopyin(vx,vy,vz, GMTR_A_var,GMTR_T_var) 
    !$acc loop gang vector(1)
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin, ADM_kmax
    !$acc loop gang vector(64)
    do n = nstart, nend
       ij     = n
       ip1j   = n + 1
       ijp1   = n     + ADM_gall_1d
       ip1jp1 = n + 1 + ADM_gall_1d

       sclt(n,k,l,ADM_TI) = ( - ( vx(ij    ,k,l) + vx(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNX) &
            - ( vy(ij    ,k,l) + vy(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNY) &
            - ( vz(ij    ,k,l) + vz(ip1j  ,k,l) ) * GMTR_A_var(ij,  k0,l,AI, TNZ) &
            - ( vx(ip1j  ,k,l) + vx(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNX) &
            - ( vy(ip1j  ,k,l) + vy(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNY) &
            - ( vz(ip1j  ,k,l) + vz(ip1jp1,k,l) ) * GMTR_A_var(ip1j,k0,l,AJ, TNZ) &
            + ( vx(ip1jp1,k,l) + vx(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
            + ( vy(ip1jp1,k,l) + vy(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
            + ( vz(ip1jp1,k,l) + vz(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
            ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TI,GMTR_T_RAREA)

       sclt(n,k,l,ADM_TJ) = ( - ( vx(ij    ,k,l) + vx(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNX) &
            - ( vy(ij    ,k,l) + vy(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNY) &
            - ( vz(ij    ,k,l) + vz(ip1jp1,k,l) ) * GMTR_A_var(ij,  k0,l,AIJ,TNZ) &
            + ( vx(ip1jp1,k,l) + vx(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNX) &
            + ( vy(ip1jp1,k,l) + vy(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNY) &
            + ( vz(ip1jp1,k,l) + vz(ijp1  ,k,l) ) * GMTR_A_var(ijp1,k0,l,AI, TNZ) &
            + ( vx(ijp1  ,k,l) + vx(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNX) &
            + ( vy(ijp1  ,k,l) + vy(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNY) &
            + ( vz(ijp1  ,k,l) + vz(ij    ,k,l) ) * GMTR_A_var(ij,  k0,l,AJ, TNZ) &
            ) * 0.5D0 * GMTR_T_var(ij,k0,l,ADM_TJ,GMTR_T_RAREA)
    enddo
    enddo
    enddo
    !$acc end kernels

    nstart = suf(ADM_gmin,  ADM_gmin  )
    nend   = suf(ADM_gmax,  ADM_gmax  )

    !$acc kernels pcopy(grdx,grdy,grdz) pcopyin(sclt, GMTR_A_var,GMTR_P_var, ADM_prc_tab,ADM_rgn_vnum) 
    !$acc loop gang vector(1)
    do l = 1, ADM_lall
    !$acc loop gang vector(4)
    do k = ADM_kmin-1, ADM_kmax+1
    !$acc loop gang vector(64)
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
          rgnid = ADM_prc_tab(l,ADM_prc_me)
          if (( ADM_rgn_vnum(ADM_W,rgnid) == 3 ) .and. ( n == nstart )) then
             sclt5 = sclt(ijm1  ,k,l,TJ) + sclt(im1jm1,k,l,TJ)
             sclt6 = 0.0
          endif

          grdx(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNX) &
               + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNX) &
               + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNX) &
               - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNX) &
               - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
               - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNX) &
               ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

          grdy(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNY) &
               + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNY) &
               + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNY) &
               - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNY) &
               - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
               - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNY) &
               ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact

          grdz(n,k,l) = ( + ( sclt1 ) * GMTR_A_var(ij,    k0,l,AI, HNZ) &
               + ( sclt2 ) * GMTR_A_var(ij,    k0,l,AIJ,HNZ) &
               + ( sclt3 ) * GMTR_A_var(ij,    k0,l,AJ, HNZ) &
               - ( sclt4 ) * GMTR_A_var(im1j,  k0,l,AI, HNZ) &
               - ( sclt5 ) * GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
               - ( sclt6 ) * GMTR_A_var(ijm1,  k0,l,AJ, HNZ) &
               ) * 0.5D0 * GMTR_P_var(ij,k0,l,GMTR_P_RAREA) * fact
       else 
          grdx(n,k,l) = 0.D0
          grdy(n,k,l) = 0.D0
          grdz(n,k,l) = 0.D0
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ADM_prc_me == ADM_prc_pl ) then
       n = ADM_GSLF_PL
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                sclt_pl(v,k) = ( + ( vx_pl(n   ,k,l) + vx_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNX ) &
                                 + ( vy_pl(n   ,k,l) + vy_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNY ) &
                                 + ( vz_pl(n   ,k,l) + vz_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNZ ) &
                                 + ( vx_pl(ij  ,k,l) + vx_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2X) &
                                 + ( vy_pl(ij  ,k,l) + vy_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Y) &
                                 + ( vz_pl(ij  ,k,l) + vz_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Z) &
                                 - ( vx_pl(ijp1,k,l) + vx_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNX ) &
                                 - ( vy_pl(ijp1,k,l) + vy_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNY ) &
                                 - ( vz_pl(ijp1,k,l) + vz_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNZ ) &
                               ) * 0.5D0 * GMTR_T_var_pl(ij,k0,l,GMTR_T_RAREA)
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

             grdx_pl(n,k,l) = grdx_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
             grdy_pl(n,k,l) = grdy_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
             grdz_pl(n,k,l) = grdz_pl(n,k,l) * 0.5D0 * GMTR_P_var_pl(n,k0,l,GMTR_P_RAREA) * fact
          enddo

       enddo
    endif

    !--- !$acc wait

    !$acc end data

    call DEBUG_rapend('++++OPRT_divdamp')

    return
  end subroutine OPRT_divdamp

end module mod_oprt
!-------------------------------------------------------------------------------
