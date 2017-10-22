!==========================================================================
! This is the FORTRAN implementation of the 2-FAST algorithm.
! ref: Gebhardt and Jeong 
!	"2-FAST: Fast and Accurate Computation of Projected Two-point Functions"
! 
! 3 Autust 2017
! Donghui Jeong
!==========================================================================
PROGRAM twoFAST_driver
!--------------------------------------------------------------------------
 USE header, only:double,single
 USE interpolation, only:open_file,close_file,func,xmin,xmax

 Implicit none
 include "fftw3.f"

 real(kind=double), parameter :: pi = 3.1415926535897932385D0
 real(kind=double), parameter :: twopi = 6.2831853071795864769D0

 integer, parameter :: ell = 42
 real(kind=double), parameter :: nu = 0

! power spectrum file
 character(len=128) :: pkfname
 character(len=128) :: phiMfname
 character(len=128) :: outfname
! size of the Fourier grid
 integer, parameter :: nmesh = 32!1024
 integer, parameter :: nnmesh = nmesh-1
 integer, parameter :: cnmesh = nmesh/2
! for the FFTlog array
 real(kind=double), parameter :: kmin = 1.d-5
 real(kind=double), parameter :: kmax = 1.d3
 real(kind=double) :: Deltalnk, dlnk, kval, pkval
 integer :: kindx
 real(kind=double) :: qbias !, parameter :: qbias = 0.9d0
 character(len=4) :: c_qbias
 integer(kind=double) :: plan_FFT, plan_iFFT
 real(kind=double),allocatable,dimension(:) :: biasedPk
 real(kind=double),allocatable,dimension(:) :: biasedPk_check
 complex(kind=double),allocatable,dimension(:) :: phiq
 integer :: status
! For the inverseFFTlog
 real(kind=double) :: tFund, tval, tTotal
 integer :: tindx
 real(kind=double),parameter :: k0 = kmin
 real(kind=double),parameter :: r0 = 1.d-3
 real(kind=double),parameter :: alpha = k0*r0

 integer :: rindx
 real(kind=double) :: rho,wellell,dlgamma
 complex(kind=double) :: cdlgamma, elltest
 integer :: qindx 

 pkfname = '../planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat'
 ! open_file(filename, skip, number of lines, number of columns, x column, ycolumn, log_interpolation?, log_interpolation?)
 CALL open_file(pkfname,1,988,2,1,2,.true.,.true.)

 print*, xmin,xmax

! allocate the array
 allocate(biasedPk(0:nnmesh),STAT=status)
 if(status /= 0) then
   print *,'biasedPk: failed to allocate'
 endif
 allocate(phiq(0:cnmesh),STAT=status)
 if(status /= 0) then
   print *,'biasedPk: failed to allocate'
 endif

do qindx=9,19
! initialize FFTW
 print *,'initialize the FFTW'
 CALL dfftw_plan_dft_r2c_1d(plan_FFT,nmesh,biasedPk,phiq,FFTW_ESTIMATE)
 CALL dfftw_plan_dft_c2r_1d(plan_iFFT,nmesh,phiq,biasedPk,FFTW_ESTIMATE)
 print *,'Done!'


 qbias = dble(qindx)/10.d0

! Fill the array
 Deltalnk = dlog(kmax) - dlog(kmin)
 dlnk     = Deltalnk/dble(nmesh)
 do kindx=0,nnmesh
   kval = k0 * dexp(dlnk*dble(kindx))
   pkval = func(kval)
   biasedPk(kindx) = pkval*(kval/k0)**(3.d0-qbias)
 enddo

! FFTlog
 print *,'FFT!'
 CALL dfftw_execute(plan_FFT)
 print *,'Done!'

! fundamental frequency
 tFund = twopi/Deltalnk
 tTotal = tFund * dble(nmesh)
! modifying the phiq array 
 write(c_qbias,'(f4.1)')qbias
 phiMfname = 'phiqMq_q='//trim(adjustl(c_qbias))//'.dat'
 open(13,file=phiMfname,form='formatted')

 do tindx=0,cnmesh
    tval = tFund * dble(tindx)

    write(13,*) tval,cdabs(phiq(tindx))/Deltalnk,cdabs(Mellq(ell,qbias,tval))
!    write(13,*) tval,phiq(tindx)

!    print *, tval, cdabs(phiq(tindx))
    phiq(tindx) = conjg(phiq(tindx))/Deltalnk * Mellq(ell,qbias,tval)
!    print *, Mellq(ell,qbias,tval)
 enddo
 close(13)

 CALL dfftw_execute(plan_iFFT)
 
 outfname = 'wllp_q='//trim(adjustl(c_qbias))//'.dat'
 open(13,file=outfname,form='formatted')
 do rindx=0,nnmesh
    rho = dble(rindx)*twopi/tTotal
    wellell = 4.d0*k0**3./tTotal*dexp(-qbias*rho)*biasedPk(rindx)
    write(13,*)r0*dexp(rho),wellell
 enddo
 close(13)


! Wrapping up: destroy plan, deallocate arrays
 CALL dfftw_destroy_plan(plan_FFT)
 CALL dfftw_destroy_plan(plan_iFFT)

enddo
 deallocate(biasedPk,phiq)
 CALL close_file

CONTAINS
 complex(kind=double) Function Mellq(ell,q,t)
  implicit none
  integer, intent(in) :: ell
  real(kind=double), intent(in) :: q,t
  real(kind=double) :: dell
  complex(kind=double) :: n,cdlgamma,Ullp,lnG2F1

  dell = dble(ell)
  n = dcmplx(q-1.d0, -t)

  lnG2F1  = cdlgamma(dcmplx(1.5d0+dell,0.d0)) + cdlgamma(1.d0-n)&
          - cdlgamma((3.d0+2.d0*dell-n)/2.d0)                  &
          - cdlgamma((2.d0-n)/2.d0)

  Ullp  = 2.d0**(n-2.d0)*pi                             &
        * cdexp( cdlgamma((1.d0+2.d0*dell+n)/2.d0)      &
               - cdlgamma((2.d0-n)/2.d0)                &
               - cdlgamma(dcmplx(1.5d0+dell,0.d0)) + lnG2F1)        

  Mellq = alpha**dcmplx(-q,t)*Ullp

  return
 End Function Mellq
END PROGRAM twoFAST_driver
!==========================================================================
