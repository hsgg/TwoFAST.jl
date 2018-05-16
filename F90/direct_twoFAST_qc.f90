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
! include "fftw3.f"

 integer, parameter :: quad = kind(0.q0)

 real(kind=quad), parameter :: pi = 3.1415926535897932384626433832795028842q0
 real(kind=quad), parameter :: twopi = 6.2831853071795864769252867665590057684q0
 integer, parameter :: ell = 42
 real(kind=quad), parameter :: nu = 0.q0

! power spectrum file
 character(len=128) :: pkfname
 character(len=128) :: phiMfname
 character(len=128) :: outfname
! size of the Fourier grid
 integer, parameter :: nmesh = 129
 real(kind=quad) :: qnmesh
 integer, parameter :: nnmesh = nmesh-1
 integer, parameter :: cnmesh = nmesh/2
! for the FFTlog array
 real(kind=quad), parameter :: kmin = 1.q-5
 real(kind=quad), parameter :: kmax = 1.q3
 real(kind=quad) :: Deltalnk, dlnk, kval, pkval
 integer :: kindx
 real(kind=quad) :: qbias !, parameter :: qbias = 0.9d0
 character(len=4) :: c_qbias
 integer(kind=double) :: plan_phiq, plan_iresult

 real(kind=quad),allocatable,dimension(:) :: biasedPk
! complex(kind=quad),allocatable,dimension(:) :: biasedPk
 complex(kind=quad),allocatable,dimension(:) :: phiq

 complex(kind=quad),allocatable,dimension(:) :: phiq_all
 complex(kind=quad),allocatable,dimension(:) :: iresult

 integer :: status

! For the inverseFFTlog
 real(kind=quad) :: tFund, tval, tTotal
 integer :: tindx
 real(kind=quad),parameter :: k0 = kmin
 real(kind=quad),parameter :: r0 = 1.d-3
 real(kind=quad),parameter :: alpha = k0*r0

 integer :: rindx
 real(kind=quad) :: rho,wellell
 complex(kind=quad) :: cqlgamma
 integer :: qindx 

 complex(kind=quad) :: phiq_this, Mellq_this

 qnmesh = real(nmesh,quad)

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
   print *,'phiq: failed to allocate'
 endif
 allocate(phiq_all(0:nnmesh),STAT=status)
 if(status /= 0) then
   print *,'phiq_all: failed to allocate'
 endif
 allocate(iresult(0:nnmesh),STAT=status)
 if(status /= 0) then
   print *,'iresult: failed to allocate'
 endif

! initialize FFTW
 print *,'initialize the FFTW'
! CALL dfftw_plan_dft_r2c_1d(plan_phiq,nmesh,biasedPk,phiq,FFTW_ESTIMATE)
! CALL dfftw_plan_dft_1d(plan_iresult,nmesh,phiq_all,iresult,FFTW_BACKWARD,FFTW_ESTIMATE)
 do kindx = 0, nnmesh
    biasedPk(kindx) = 0.q0
    if(kindx .le. cnmesh) phiq(kindx) = qcmplx(0.q0,0.q0)
    phiq_all(kindx) = qcmplx(0.q0,0.q0)
    iresult(kindx) = qcmplx(0.q0,0.q0)
 enddo
 print *,'Done!'

do qindx=9,19
 qbias = real(qindx,quad)/10.q0

! Fill the array
 Deltalnk = qlog(kmax) - qlog(kmin)
 dlnk     = Deltalnk/real(nmesh,quad)
 do kindx=0,nnmesh
   kval = k0 * qexp(dlnk*real(kindx,quad))
   pkval = real(func(dble(kval)),quad)
   biasedPk(kindx) = qcmplx(pkval*(kval/k0)**(3.q0-qbias),0.q0)
 enddo

! FFTlog
 print *,'FFT!'
! CALL dfftw_execute(plan_phiq)
 do tindx = 0, cnmesh
    phiq(tindx) = sum((/(biasedPk(kindx)*cqexp(qcmplx(0.d0,twopi*real(tindx*kindx,quad)/qnmesh)),kindx=0,nnmesh)/))
 enddo
 print *,'Done!'

! fundamental frequency
 tFund = twopi/Deltalnk
 tTotal = tFund * qnmesh
! modifying the phiq array 
 write(c_qbias,'(f4.1)')qbias
 phiMfname = 'direct_qc_phiqMq_q='//trim(adjustl(c_qbias))//'.dat'
 open(13,file=phiMfname,form='formatted')

 do tindx=0,nnmesh
    if(tindx .le. cnmesh) then
      tval = tFund * real(tindx,quad)
      phiq_all(tindx) = phiq(tindx)/Deltalnk
    else
      tval = tFund * real(tindx-nmesh,quad)
      phiq_all(tindx) = qconjg(phiq(nmesh-tindx))/Deltalnk
    endif

    phiq_this = phiq_all(tindx)
    Mellq_this = Mellq(ell,qbias,tval)

!    write(13,'(3e40.30)') tval,cqabs(phiq_all(tindx)),cqabs(Mellq(ell,qbias,tval))
    phiq_all(tindx) = phiq_this * Mellq_this
    if((tindx .eq. cnmesh) .or. (tindx .eq. 0)) phiq_all(tindx) = qcmplx(qreal(phiq_this*Mellq_this),0.q0)
    write(13,'(7e30.20)') tval,qreal(phiq_this),qimag(phiq_this),qreal(Mellq_this),qimag(Mellq_this),qreal(phiq_all(tindx)),qimag(phiq_all(tindx))
 enddo
 
 close(13)

! CALL dfftw_execute_dft(plan_iresult,phiq_all,iresult)
 do rindx = 0, nnmesh
    iresult(rindx) = sum((/(phiq_all(tindx)*cqexp(qcmplx(0.d0,twopi*real(rindx*tindx,quad)/qnmesh)),tindx=0,nnmesh)/))
 enddo

 outfname = 'direct_qc_wllp_q='//trim(adjustl(c_qbias))//'.dat'
 open(13,file=outfname,form='formatted')
 do rindx=0,nnmesh
    rho = real(rindx,quad)*twopi/tTotal
    wellell = 4.d0*k0**3./tTotal*qexp(-qbias*rho)*qreal(iresult(rindx))
    write(13,'(3e40.30)')r0*qexp(rho),wellell,qimag(iresult(rindx))
 enddo
 close(13)

enddo

! Wrapping up: destroy plan, deallocate arrays
! CALL dfftw_destroy_plan(plan_phiq)
! CALL dfftw_destroy_plan(plan_iresult)
 deallocate(biasedPk,phiq,phiq_all,iresult)
 CALL close_file

CONTAINS
 complex(kind=quad) Function Mellq(ell,q,t)
  implicit none
  integer, intent(in) :: ell
  real(kind=quad), intent(in) :: q,t
  real(kind=quad) :: dell
  complex(kind=quad) :: n,cqlgamma,Ullp,lnG2F1

  dell = dble(ell)
  n = qcmplx(q-1.q0, -t)

  lnG2F1  = cqlgamma(qcmplx(1.5q0+dell,0.q0)) + cqlgamma(qcmplx(1.q0-n))&
          - cqlgamma(qcmplx((3.q0+2.q0*dell-n)/2.q0))                 &
          - cqlgamma(qcmplx((2.q0-n)/2.q0))

  Ullp  = 2.q0**(n-2.q0)*pi                             &
        * cqexp( cqlgamma(qcmplx((1.q0+2.q0*dell+n)/2.q0))      &
               - cqlgamma(qcmplx((2.q0-n)/2.q0))                &
               - cqlgamma(qcmplx(1.5q0+dell,0.q0)) + lnG2F1)        

  Mellq = alpha**qcmplx(-q,t)*qcmplx(Ullp)

  return
 End Function Mellq
END PROGRAM twoFAST_driver
!==========================================================================
