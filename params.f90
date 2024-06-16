! ######################### !
! ### Global parameters ### !
! ######################### !
! last edited 04/27/2023 S. Malpathak
module parameters
implicit none

! pi and the imaginary unit
real*8, parameter     :: pi = dacos(-1.d0)
complex*16, parameter :: Iu = (0.d0,1.d0)
real*8,parameter      :: EVtoau=3.67493d-2, Angtoau=1.8897d0, fstoau=4.134d1
real*8,parameter      :: kB=1.380649d-23, JtoEV=6.242d18
real*8                :: autoEV, autoAng, autofs

! System variables
real*8             :: MassN, MassO
real*8,allocatable :: Mass(:,:), MassInv(:,:)
!Parameters for potentials
real*8  :: D0, A0, R0, b0, Z0, C0       ! For NO
real*8  :: D1, A1, R1, D2, A2, Z1, C1   ! For NO-
real*8  :: Gammab, Vbar, Abar, kC       ! Coupling

! Coherent state parameters for Z
real*8  :: Zi, Ei, GammaI, pZi
! Initial vibrational quantum #
integer :: NuI

! Matrix for initial vib. projection matrix in DVR basis
real*8,allocatable  :: VibProjMat(:,:,:), VibRho(:,:)

integer :: NumberMCsteps ! number of trajectory pairs
integer :: Ncheck        ! check convergence after using 'NumberMCsteps'
                         ! trajectories, repeat 'Ncheck' times
! propagator parameters
integer :: Ntime                     ! number of timesteps
real*8  :: TimeStep, EnergyTolerance ! timestep, energy tolerance

integer :: Ndof ! full system dimensionality
integer :: Nnuc, NelNO, Nbath, Nel

!Parameters for the metal surface
real*8  :: BetaT, BW, Mu        ! Temp. (K), width of condunction band and fermi energy in EV
real*8  :: Emax, Emin,dE        ! Range of Es: Emin to +Emax, dE is energy spacing 
real*8,allocatable  :: Ek(:)    ! Energies of bath modes

! DVR grid for Wigner Transform
integer                 :: NR
real*8                  :: Rmin, Rmax
real*8                  :: dR
real*8, allocatable     :: gridR(:)

!Step size and its guess for MC for R dof in a.u.
real*8                  :: MCStepGuess, MCstep
integer                 :: NDCorr !Decorrelation length for MC sampling
!Gamma for MMST sampling
real*8                  :: GMMST
!Occupation numbers for metal states
integer,allocatable :: OccN(:)

end module

! ########################### !
! ### Read the input file ### !
! ########################### !
subroutine input
use parameters
implicit none

integer         :: i, j, k
character*75    :: infostr

!unit conversion
autoEV  = 1.d0/EVtoau
autoAng = 1.d0/Angtoau
autofs  = 1.d0/fstoau

open(555,file='input',status='old')
read(555,'(a75)') infostr
read(555,*) Nnuc, NelNO, Nbath

Nel  = NelNO + Nbath
Ndof = Nnuc  + Nel

allocate(Mass(Nnuc,Nnuc),MassInv(Nnuc,Nnuc))
allocate(Ek(Nbath),OccN(Nbath))

Mass    = 0.d0
MassInv = 0.d0
Ek      = 0.d0
OccN    = 0

read(555,'(a75)') infostr
read(555,*)  MassN, MassO

Mass(2,2) = MassN + MassO !Mass of translational dof
Mass(1,1) = 1.d0/(1.d0/MassN + 1.d0/MassO) ! Reduced mass for vib. dof

!Parameters for R morse oscillators for NO and NO-
read(555,'(a75)') infostr
read(555,*)  D0, A0, R0
read(555,*)  D1, A1, R1 
!Convert to a.u.
D0      = D0*EVtoau
D1      = D1*EVtoau
A0      = A0/Angtoau
A1      = A1/Angtoau
R0      = R0*Angtoau
R1      = R1*Angtoau

!Parameters for Z potentials
read(555,'(a75)') infostr
read(555,*) B0, Z0
read(555,*) D2, A2, Z1
!Convert to a.u.
B0      = B0/Angtoau
Z0      = Z0*Angtoau
D2      = D2*EVtoau
A2      = A2/Angtoau
Z1      = Z1*Angtoau

!Constants for both surfaces
read(555,'(a75)') infostr
read(555,*) C0, C1
C0      = C0*EVtoau
C1      = C1*EVtoau

!Parameters for metal-bath coupling 
read(555,'(a75)') infostr
read(555,*)       GammaB, Abar
GammaB  = GammaB*EVtoau
ABar    = Abar*Angtoau
Vbar    = dsqrt(GammaB/2.d0/pi)

!Scaling factor for metal-molecule coupling
read(555,'(a75)') infostr
read(555,*)       kC
!kc is 1.d0
!Parameters for metal-bath coupling and temp. of surface
read(555,'(a75)') infostr
read(555,*)     BetaT, BW, Mu
!BetaT is read in as T in K
BetaT = 1.d0/kB/BetaT           !Convert T to Beta in joules**-1
BetaT = BetaT/JtoEV             !Convert BetaT to ev**-1
BetaT = BetaT/EVtoau            !Convert BetaT to a.u.^-1
BW    = BW*EVtoau
Mu    = Mu*EVtoau

read(555,'(a75)') infostr
read(555,*)  Zi, Ei, GammaI
!Convert to a.u. and calculate PzI
Zi      = Zi*Angtoau
GammaI  = GammaI/Angtoau**2
pZi     = -dsqrt(2.d0*Mass(2,2)*Ei*EVtoau) ! Energy converted to a.u. first

read(555,'(a75)') infostr
read(555,*) NuI

read(555,'(a75)') infostr
read(555,*)       Rmin, Rmax, NR
Rmin = Rmin*Angtoau
Rmax = Rmax*Angtoau

read(555,'(a75)') infostr
read(555,*)       MCstepGuess, NDCorr

read(555,'(a75)') infostr
read(555,*)       GMMST

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance
!Convert to a.u.
Timestep= Timestep*fstoau

read(555,'(a75)') infostr
read(555,*) NumberMCsteps

read(555,'(a75)') infostr
read(555,*) Ncheck

close(555)

! Spacing in R grid
dR = (Rmax-Rmin)/dble(NR-1)
allocate(gridR(NR),VibRho(NR,NR),VibProjMat(0:NuI+2,NR,NR))
gridR           = 0.d0
VibRho          = 0.d0
VibProjMat      = 0.d0
!Build R grid
do i = 0, NR-1
  gridR(i+1) = Rmin + dble(i)*dR
enddo

!Define Inverse Mass matrix
do i = 1, Nnuc
   MassInv(i,i) = 1.d0/Mass(i,i)
end do
! Calculate the Proj. operator onto initial vib.state in DVR basis
call Amatrix(VibRho)
call ProjMatrix(VibProjMat)


!Get energy levels and their couplings
!Energy levels are b.w Emin and +Emax centered around E=mu
Emin = Mu - BW/2.d0
Emax = Mu + BW/2.d0
dE   = (Emax-Emin)/dble(Nbath) 
do i = 0, Nbath-1
   Ek(i+1)        = Emin + dble(i)*dE 
   !write(*,*), i+1, Ek(i+1)*autoEV
end do

end subroutine input
!########################################
!       A matrix in DVR basis
!########################################
subroutine Amatrix(Amat)
use parameters, only : NR, pi, NuI, Mass, D0, A0, R0, gridR
implicit none

real*8,intent(out)  :: Amat(NR,NR)
real*8              :: tr, lambda, w, wp, lg, lgp
integer             :: i, ip

Amat    = 0.d0
tr      = 0.d0
lambda  = dsqrt(2.d0*Mass(1,1)*D0)/A0
do i = 1, NR
   w = 2.d0*lambda*dexp(-A0*(gridR(i)-R0))
   call LaguerreP(NuI,2.d0*lambda-2.d0*dble(NuI)-1.d0,w,lg)
   do ip = 1, NR
       wp = 2.d0*lambda*dexp(-A0*(gridR(ip)-R0))
       call LaguerreP(NuI,2.d0*lambda-2.d0*dble(NuI)-1.d0,wp,lgp)
       Amat(i,ip)= (w**(lambda-dble(NuI)-0.5d0))*dexp(-0.5d0*w)*lg & 
                 * (wp**(lambda-dble(NuI)-0.5d0))*dexp(-0.5d0*wp)*lgp &
                 * 2.d0*gamma(dble(NuI+1))*A0*(lambda-0.5d0-dble(NuI)) &
                 / gamma(2.d0*lambda-dble(NuI))
   end do
end do

do i = 1, NR
   tr = tr + Amat(i,i)
end do

Amat = Amat/tr   ! Normalize initial density.
           
end subroutine Amatrix
!########################################
!       Proj onto vib. state Nu matrix in DVR basis
!########################################
subroutine ProjMatrix(Pmat)
use parameters, only : NR, pi, Nui, Mass, D0, A0, R0, gridR
implicit none

real*8,intent(out)  :: Pmat(0:NuI+2,NR,NR)
real*8              :: tr, lambda, w, wp, lg, lgp
integer             :: i, ip, Nu

Pmat    = 0.d0
tr      = 0.d0
lambda  = dsqrt(2.d0*Mass(1,1)*D0)/A0
do Nu = 0, NuI+2
 do i = 1, NR
   w = 2.d0*lambda*dexp(-A0*(gridR(i)-R0))
   call LaguerreP(Nu,2.d0*lambda-2.d0*dble(Nu)-1.d0,w,lg)
   do ip = 1, NR
       wp = 2.d0*lambda*dexp(-A0*(gridR(ip)-R0))
       call LaguerreP(Nu,2.d0*lambda-2.d0*dble(Nu)-1.d0,wp,lgp)
       Pmat(Nu,i,ip) = (w**(lambda-dble(Nu)-0.5d0))*dexp(-0.5d0*w)*lg & 
                     * (wp**(lambda-dble(Nu)-0.5d0))*dexp(-0.5d0*wp)*lgp &
                     * 2.d0*gamma(dble(Nu+1))*A0*(lambda-0.5d0-dble(Nu)) &
                     / gamma(2.d0*lambda-dble(Nu))
   end do
 end do
end do
        
end subroutine Projmatrix
!##############################################################
!       Calculates the generalized Laguerre polynomial L_n^(a)(x)
!###############################################################
subroutine LaguerreP(n,a,x,LP)
implicit none

integer, intent(in) :: n
real*8, intent(in)  :: a, x
integer             :: i
real*8, intent(out) :: LP

LP = 0.d0

do i = 0, n
  LP = LP + ((-x)**i)*gamma(dble(n+a+1))/gamma(dble(n-i+1))/gamma(dble(a+i+1))/gamma(dble(i+1))
end do

end subroutine LaguerreP

