! ################################################### !
! ### this program computes observables for       ### !
! ### Eric's model of NO scattering using LSC-IVr ### !
! ################################################### !
! last edited 06/26/2023 S.V.Malpathak
program CF
use parameters
implicit none

include 'mpif.h'

integer             :: i, j, k, l, m
integer             :: flagE, flagR
integer             :: countE, brokenE, countR, brokenR
integer*8           :: brokenEF, trajused, brokenRF
integer             :: ierr, myrank, nprocs, istart, iend, length, imc
integer             :: MCMoveC, MCSucC, MCMoveT, MCSucT
real*8, allocatable :: PScoord(:,:), qpfwd(:,:,:), pop11(:,:) 
real*8, allocatable :: Initialq(:), Initialp(:)
real*8, allocatable :: p2(:,:), pop2(:,:,:), Rt(:), Zt(:), WigNormT(:)
real*8, allocatable :: rTCF(:), iTCF(:), rTCF2(:,:), iTCF2(:,:), elem11(:)
real*8, allocatable :: AverR(:), AverZ(:), AverRt(:,:), AverZt(:,:)
real*8, allocatable :: TProj(:,:), TProjT(:,:,:), TP(:,:), TProjW(:)
real*8              :: TCF, dist0

real*8              :: normC, summ, MWT
real*8              :: WigNorm, WigNormC

call input

allocate(PScoord(2,Ndof),qpfwd(2,Ndof,0:Ntime),Zt(0:Ntime))
allocate(pop11(Nel,0:Ntime),Rt(0:Ntime),elem11(Nel))
allocate(AverR(0:Ntime),AverZ(0:Ntime),WigNormT(Ncheck))
allocate(AverRt(Ncheck,0:Ntime),AverZt(Ncheck,0:Ntime))
allocate(Initialq(Ndof),Initialp(Ndof))
allocate(p2(Nel,0:Ntime),pop2(Ncheck,Nel,0:Ntime))
allocate(TProjW(0:Nui+2),TP(0:NuI+2,0:100))
allocate(TProj(0:Nui+2,0:100),TProjT(Ncheck,0:NuI+2,0:100))


rTCF2   = 0.d0
iTCF2   = 0.d0
p2      = 0.d0
pop11   = 0.d0
pop2    = 0.d0
elem11  = 0.d0
AverR   = 0.d0
AverZ   = 0.d0
AverRt  = 0.d0
AverZt  = 0.d0
WigNormC= 0.d0
WigNormT= 0.d0
MCSucT  = 0
MCMoveT = 0
Tproj   = 0.d0
TprojT  = 0.d0
TprojW  = 0.d0

! accumulate broken trajectory info
brokenEF = 0    ! broken energy
brokenRF = 0    ! broken b/c inv. pot.
trajused = 0    ! total trajectories used

! set up parallelization
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,myrank,ierr)
call mpi_barrier(mpi_comm_world,ierr)

call init_random_seed(myrank)
!Get MC step size to sample Morse dof
if (myrank.eq.0) then
   !Set initial (R,P) to (R0,0) and get initial WT
   Initialq(1) = R0
   Initialp(1) = 0.d0
   call Morse_Wigner(Initialq(1),Initialp(1),MWT)
   !Get MC step size
   call MCstepsize(Initialq(1),Initialp(1),MWT)   

   !write percentage of broken trajectories
   open(111,file='Traj_Info.out',status='unknown')
end if

call mpi_barrier(mpi_comm_world,ierr)
! Broadcast rest of the variables.
call mpi_bcast(MCStep,1,mpi_double_precision,0,mpi_comm_world,ierr)

do m = 1, Ncheck ! write correlation function at each m

 countE  = 0     ! count broken energy
 countR  = 0     ! count inv. pot trajs
 p2      = 0.d0
 Rt      = 0.d0
 Zt      = 0.d0
 WigNorm = 0.d0
 TP      = 0.d0

 call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)

 ! Initialize MC with a little bit of sampling
 Initialq(1) = R0
 Initialp(1) = 0.d0
 call Morse_Wigner(Initialq(1),Initialp(1),MWT)
 do j = 1, 100
   call MCmove(Initialq(1),Initialp(1),MWT,MCSucC,MCMoveC)     
 end do    
 MCMoveC = 0
 MCSucC  = 0
 do imc = istart, iend  ! loop over trajectories
   ! sample scheme for LSC
   call LSCsample(Initialq(1),Initialp(1),MWT,MCSucC,MCMoveC,PScoord,OccN)
   flagE  = 0
   flagR  = 0
   ! Initial conditions with LSCsample sampling
   Initialq = PScoord(1,:)
   Initialp = PScoord(2,:)
   ! propagate trajectory
   call PropagateFwd(Initialq,Initialp,qpfwd,flagE,flagR,myrank)
   ! check for conservation of ...
   if (flagE.eq.1) then                 ! energy
      countE = countE + 1
      goto 222
   endif
   if (flagR.eq.1) then                 ! Inv. pot
      countR = countR + 1
      goto 222
   endif

   flagE = 0 
   flagR = 0
   ! The initial rho_wig/LSCsample function
   dist0 = dsign(1.d0,MWT) ! Estimator from WT of Morse dof
   ! Keep track of normalization for Morse wigner MC sampling
   WigNorm = WigNorm + dsign(1.d0,MWT)

   do i = 0, Ntime
      TCF = dist0
      !Calculate average R and Z and vib. transition probs at time t
      Rt(i) = Rt(i) + TCF*qpfwd(1,1,i)
      Zt(i) = Zt(i) + TCF*qpfwd(1,2,i)
      if (modulo(i,Ntime/100).eq.0) then
        ! Calculate the wig.trans of proj onto vib. states 0 to NuI+2
        call ProjWigner(qpfwd(1,1,i),qpfwd(2,1,i),TProjW)
        ! Calculate the trans. prob to go all states 0 to Nui+2
        TP(:,i/(Ntime/100)) = TP(:,i/(Ntime/100)) + TCF*TProjW(:)
      end if

      !This part calculates the time-dependednt population of NO- state
      !and all the metal states using B = 1/2*(x^2 + p^2 - GMMST)
      do j = 1, Nel
         elem11(j) = 0.5d0*(qpfwd(1,Nnuc+j,i)**2+qpfwd(2,Nnuc+j,i)**2-GMMST)
      end do
      ! Populations of all elec states at time t stored.
      p2(:,i) = p2(:,i) + TCF*elem11(:)
      !}
   enddo
  222 continue

 enddo
 call mpi_barrier(mpi_comm_world,ierr)
! Gather electronic population
 call mpi_reduce(p2,pop11,Nel*(Ntime+1),mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
! Gather nuclear averages
 call mpi_reduce(Rt,AverR,Ntime+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(Zt,AverZ,Ntime+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(TP,Tproj,(NuI+3)*101,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
! gather conservation counters and wigner normalization
 call mpi_reduce(countE,brokenE,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(countR,brokenR,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(MCSucC,MCSucT,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(MCMoveC,MCMoveT,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(WigNorm,WigNormC,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

 if (myrank.eq.0) then
   pop2(m,:,:)= pop11(:,:)
   AverRt(m,:)= AverR(:)
   AverZt(m,:)= AverZ(:)
   TprojT(m,:,:)= Tproj(:,:)
   WigNormT(m)= WigNormC

   brokenEF   = brokenEF + brokenE
   brokenRF   = brokenRF + brokenR
   trajused   = trajused + (NumberMCsteps - brokenE - brokenR)
   WigNormC    = sum(WigNormT(:))/dble(trajused)
   do i = 0, Ntime
     do j = 1, Nel
        pop11(j,i) = sum(pop2(:,j,i))/dble(trajused)/WigNormC
     end do
     AverR(i)    = sum(AverRt(:,i))/dble(trajused)/WigNormC
     AverZ(i)    = sum(AverZt(:,i))/dble(trajused)/WigNormC

     write(1000+m,'(E15.6)',advance='no'), i*TimeStep*autofs
     write(1000+m,'(E15.6)',advance='no'), AverR(i)*autoAng
     write(1000+m,'(E15.6)'), AverZ(i)*autoAng

     write(3000+m,'(E15.6)',advance='no'), i*TimeStep*autofs
     do j = 1, Nel-1
        write(3000+m,'(E15.6)',advance='no'), pop11(j,i)
     end do 
     write(3000+m,'(E15.6)'), pop11(Nel,i)

     if (modulo(i,Ntime/100).eq.0) then
        do j = 0, NuI+2
          TProj(j,i/(Ntime/100)) = sum(TprojT(:,j,i/(Ntime/100)))/dble(trajused)/WigNormC
        end do

        summ = sum(TProj(:,i/(Ntime/100)))
        write(2000+m,'(E15.6)',advance='no'), dble(i)*TimeStep*autofs
        do j = 0, NuI+1
           write(2000+m,'(E15.6)',advance='no'), TProj(j,i/(Ntime/100))/summ
        end do 
        write(2000+m,'(E15.6)'), TProj(Nui+2,i/(Ntime/100))/summ
     end if
   enddo
   
   ! write broken trajectory information
   write(111,*) m
   write(111,*) 'Broken energy      :', brokenEF
   write(111,*) 'Broken b/c inv.pot :', brokenRF
   write(111,*) 'Total traj used    :', trajused
   write(111,*) 'MC Success rate    :', dble(MCSucT)/dble(MCMoveT)*1.d2
   write(111,*) ''
 endif

enddo

close(111) 

call mpi_finalize(ierr)

stop

end program CF
