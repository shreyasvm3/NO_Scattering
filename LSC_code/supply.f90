! MPI data division
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in) :: n1,n2,nprocs,irank
integer, intent(out) :: ista,iend
integer :: iwork1,iwork2
iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank*iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 -1
if (iwork2 .gt. irank) then
iend = iend + 1
endif
return
end subroutine para_range
!#############################################################
! initiate seed
subroutine init_random_seed(myrank)
integer :: myrank
integer :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) + myrank*clock

call random_seed(put = seed)

deallocate(seed)
end subroutine init_random_seed
!#############################################################
!       Calculate the Wigner transform of projection   
!       onto |NuI><NuI| using a DVR grid 
!       for all states 0 to NuI+2. (not normalized)
!#############################################################
subroutine ProjWigner(R,P,PWT)
use parameters, only : NR, NuI, VibProjMat 
implicit none

real*8, intent(in)      :: R, P
real*8, intent(out)     :: PWT(0:NuI+2)
complex*16              :: DMat(NR,NR), cf
integer                 :: i

call DeltaMatrix(R,P,DMat)
do i = 0, NuI+2
  call corrf(VibProjMat(i,:,:),DMat,cf)
  PWT(i) = real(cf)
end do

end subroutine ProjWigner
!#############################################################
!       Calculate the Wigner transform of projection   
!       onto |NuI><NuI| using a DVR grid
!#############################################################
subroutine Morse_Wigner(R,P,MWT)
use parameters, only : NR, VibRho 
implicit none

real*8, intent(in)      :: R, P
real*8, intent(out)     :: MWT
complex*16              :: DMat(NR,NR), cf

call DeltaMatrix(R,P,DMat)
call corrf(VibRho,DMat,cf)

MWT = real(cf)

end subroutine Morse_Wigner
!###############################################################
!       Delta matrix used for wigner transform 
!###############################################################
subroutine DeltaMatrix(q,p,Delta)
use parameters, only : NR, gridR, dR, pi, Iu, Rmin, Rmax
implicit none

real*8,intent(in)       :: q, p
complex*16,intent(out)  :: Delta(NR,NR)
real*8                  :: xbar, pB, tr
integer                 :: i, j

Delta = 0.d0
pB    = pi/dR
tr    = 0.d0

do i = 1, NR
   do j = 1, NR
      xbar = 0.5d0*(GridR(i)+GridR(j))
      if (xbar.eq.q) then
         Delta(i,j) = dR/pi*2.d0*(pB-dabs(p))*cdexp(Iu*p*(GridR(i)-GridR(j)))
      else 
         Delta(i,j) = dR/pi*dsin(2.d0*(pB-dabs(p))*(xbar-q))/(xbar-q)&
                      *cdexp(Iu*p*(GridR(i)-GridR(j)))
      end if
   end do
end do

if (q.lt.Rmin.or.q.gt.Rmax) then
   print*, "R outside [Rmin,Rmax]. R =", q
elseif (dabs(p).gt.pB) then
   print*, "P outside [Pmin,Pmax]. P =", p
end if

end subroutine DeltaMatrix
!##############################################################
! Calculate the correlation function from all the matrices
!###############################################################
subroutine Corrf(Amat,Delta,cf)
use parameters, only : NR, pi
implicit none

real*8, intent(in)      :: Amat(NR,NR)
complex*16,intent(in)   :: Delta(NR,NR)
complex*16              :: temp(NR,NR)
integer                 :: i
complex*16,intent(out)  :: cf

cf      = 0.d0
temp    = matmul(Delta,Amat)

do i = 1, NR
   cf  = cf + temp(i,i)
end do

end subroutine Corrf 
!###############################################################################
!       Perform 1 MC move
!###############################################################################
subroutine MCmove(R,P,MWig,SucCount,MoveCount)
use parameters, only : MCstep, A0, Rmin, Rmax, dR, pi
implicit none

real*8, intent(inout)   :: R, P, Mwig
real*8                  :: Rp, Pp, MwigP,rnd, pacc, pB
integer                 :: flag
integer, intent(inout)  :: SucCount, MoveCount

! Perform steps in scaled dimensionless coordinates:
! Rp = A0*R and Pp = P/A0/hbar
! Hopefully, step size can be the same for these scaled coordinates

Rp      = R*A0
Pp      = P/A0
MwigP   = 0.d0
pB      = pi/dR

call random_number(rnd)
rnd = 2.d0*rnd - 1.d0
Rp = Rp + rnd*MCstep

call random_number(rnd)
rnd = 2.d0*rnd - 1.d0
Pp = Pp + rnd*MCstep
!Step size is chosen to be the same since we are using scaled coordinates

!Transform back to unscaled coordinates
Rp = Rp/A0
Pp = Pp*A0
 
111 continue
! R needs to b.w Rmin and Rmax 
if (Rp.lt.Rmin) then
   Rp = Rp + Rmax - Rmin !Move to increase R to be in range again
   goto 111              ! repeat if still outside box
elseif (Rp.gt.Rmax) then
   Rp = Rp - (Rmax -Rmin)!Move to decrease R to be in range again
   goto 111              ! repeat if still outside box
end if 

! P needs to be b/w -pi/dR to + pi/dR
112 continue
if (Pp.lt.-pB) then
   Pp = Pp + 2.d0*pB !Move to increase P to be in range again
   goto 112              ! repeat if still outside box
elseif (Pp.gt.pB) then
   Pp = Pp - 2.d0*pB !Move to decrease R to be in range again
   goto 112              ! repeat if still outside box
end if 

!Calculate WT at new point
call Morse_Wigner(Rp,Pp,MWigp)

pacc = dmin1(1.d0,dabs(MWigp)/dabs(Mwig))
call random_number(rnd)
if (rnd.le.pacc) then
  R     = Rp
  P     = Pp
  Mwig  = Mwigp
  flag  = 1
else 
  flag = 0
end if
SucCount        = SucCount  + flag
MoveCount       = MoveCount + 1

end subroutine MCmove
!#############################################################
!###############################################################################
!       Subroutine that calculates the optimum MC step size, given a guess.
!###############################################################################
subroutine mcstepsize(R,P,MWT)
use parameters, only : MCstepGuess, MCstep
implicit none

integer                 :: UBFix, LBFix, i, j, SucCount, MoveCount
real*8,intent(inout)    :: R, P, MWT
real*8                  :: UpB, LowB, AccRate

MCstep  = MCstepGuess
UBFix   = 0
LBFix   = 0
UpB     = 2.d0*MCstep
LowB    = 0.5d0*MCstep

SucCount  = 0
MoveCount = 0
do j = 1, 1000
 call MCmove(R,P,MWT,SucCount,MoveCount)
end do
Accrate = dble(SucCount)/dble(MoveCount)
 
i = 1
do while ( (i.le.100) .and. ((AccRate.lt.4.d-1).or.(AccRate.gt.6.d-1)) )
print*, i, "Step size", MCstep, "Acc. Rate", AccRate
   if (AccRate.lt.5.d-1) then
      if (i.eq.1) then
         UBFix = 1       ! Ub can only be dec when this is true.
      end if
      UpB     = MCstep
      MCstep  = 0.5d0*(MCstep + LowB)
      if (LBFix.ne.1) then
         LowB = 0.5d0*LowB
      end if 
   else if (AccRate.gt.6.d-1) then
      if (i.eq.1) then
         LBFix = 1       ! Lb can only be inc when this true
      end if
      LowB    = MCstep
      MCstep  = 0.5d0*(MCstep + UpB)
      if (UBFix.ne.1) then
         UpB = 2.d0*UpB
      end if 
   end if
   SucCount  = 0
   MoveCount = 0
   do j = 1, 1000
      call MCmove(R,P,MWT,SucCount,MoveCount)
   end do
   Accrate = dble(SucCount)/dble(MoveCount)

   i = i + 1
end do
print*, i, "Step size", MCstep, "Acc. Rate", AccRate
if (i.ge.101) then
   print*, "Not enough loops to get optimum step size"
   call exit(0)
end if
end subroutine mcstepsize

