! ##################################################
! ### sample initial conditions with a Gaussian ###
! ### random number generator.                  ###
! ##################################################
subroutine LSCsample(R,P,MWig,SucC,MoveC,PScoord,OccN)
use parameters, only : Ndof, Nnuc, Nel, Zi, pZi, GammaI, NDCorr, NelNO, Ek, Nbath, pi, GMMST
implicit none

integer,intent(inout)   ::      SucC, MoveC
real*8, intent(inout)   ::      R, P
real*8, intent(inout)   ::      MWig
real*8, intent(out)     ::      PScoord(2,Ndof)
integer, intent(out)    ::      OccN(Nbath)
real*8                  ::      gauss, rnd
integer                 ::      i

!Perform MC move to sample R,P
do i = 1, NDCorr
   call MCmove(R,P,MWig,SucC,MoveC)
end do 
PScoord(1,1) = R
PScoord(2,1) = P

!Z is sampled from W.T. of coherent state aroung Zi, Pzi
PScoord(1,2) = gauss(dsqrt(0.5d0/GammaI)) + Zi
PScoord(2,2) = gauss(dsqrt(0.5d0*GammaI)) + pZi

!For all elec dofs, we use focused sampling
!Occupied states have 0.5*(x^2+p^2-gamma) = 1
!Unoccupied states have 0.5*(x^2+p^2-gamma) = 0.

! NO starts in uncharged NO state. NO- is unpopulated.
do i = 1, NelNO
   call random_number(rnd)
   PScoord(1,Nnuc+i) = dsqrt(GMMST)*dcos(2.d0*pi*rnd)
   PScoord(2,Nnuc+i) = dsqrt(GMMST)*dsin(2.d0*pi*rnd)
end do

!Get occupation numbers for metal dofs and then sample
! Total number of occ states is variable.
! It averages to M/2
do i = 1, Nbath
  !Determine whether ith mode is occupied or not
  call random_number(rnd)
  if (rnd.gt.fermi(Ek(i))) then
     OccN(i) = 0
     call random_number(rnd)
     PScoord(1,Nnuc+NelNO+i) = dsqrt(GMMST)*dcos(2.d0*pi*rnd)
     PScoord(2,Nnuc+NelNO+i) = dsqrt(GMMST)*dsin(2.d0*pi*rnd)
  else 
     OccN(i) = 1
     call random_number(rnd)
     PScoord(1,Nnuc+NelNO+i) = dsqrt(GMMST+2.d0)*dcos(2.d0*pi*rnd)
     PScoord(2,Nnuc+NelNO+i) = dsqrt(GMMST+2.d0)*dsin(2.d0*pi*rnd)
  end if 
end do

contains 

!define fermi distribution
function fermi(ee)
use parameters, only : BetaT, Mu
implicit none
real*8 :: fermi
real*8 :: ee
fermi = 1.d0/(dexp(BetaT*(ee-Mu))+1.d0)
end function fermi

end subroutine
!##############################################
!     Gaussian random number generator
!##############################################
double precision function gauss(sigma)

implicit none
real*8, intent(in)      :: sigma
real*8                  :: w, v1, v2, l
real*8                  :: s1, s2

w = 2.d0

do
        call random_number(s1)
        call random_number(s2)

        v1 = 2.d0*s1 - 1.d0
        v2 = 2.d0*s2 - 1.d0
        w = v1*v1 + v2*v2

        if (w.lt.1.d0) exit
end do

l     = v1*sqrt(-2.d0*log(w)/(w))
l     = sigma*l
gauss = l

end function gauss
