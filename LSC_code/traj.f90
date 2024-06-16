! ############################################################
! ############################################################
!       this subroutine propagates a classical traj using the
!       ABM predictor-corrector algorithm.
! #############################################################
! #############################################################

subroutine PropagateFwd(Initialq,Initialp,coord,flagE,flagR,myrank)
use parameters, only: TimeStep, Ntime, EnergyTolerance, Ndof, Nnuc, Nel, MassInv, autofs, autoAng
implicit none

real*8, intent(in)     :: Initialq(Ndof), Initialp(Ndof)
real*8, intent(out)    :: coord(2,Ndof,0:Ntime)
integer, intent(inout) :: flagE, flagR
integer, intent(in)    :: myrank
integer                :: i, j, k
real*8                 :: InitialEnergy, Energy, H2
real*8                 :: U, dU(Nnuc), fq(4,Ndof), fp(4,Ndof)
real*8                 :: V(Nel,Nel), dV(Nel,Nel,Nnuc)
real*8, dimension(Ndof):: q, p, qp, pp, fqp, fpp
real*8                 :: qc(4,Ndof), pc(4,Ndof)

  flagE  = 0
  flagR  = 0
  ! initial conditions
  q = Initialq
  p = Initialp
  fp = 0.d0 !force for p
  fq = 0.d0 !force for q
  fqp= 0.d0
  fpp= 0.d0

  ! total initial energy
  call Hamiltonian2(q(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),H2)
  InitialEnergy = 0.5d0*dot_product(p(1:Nnuc),matmul(MassInv,p(1:Nnuc))) + H2

  ! initial conditions
  coord(1,:,0) = Initialq
  coord(2,:,0) = Initialp

  ! get V(R) matrix elements with new nuclear position
  call potbits(q(1:Nnuc),U,V,dU,dV)  
  !Store forces at t=0  in f1
  call Forces(q,p,dU,V,dV,fq(1,:),fp(1,:))
  qc(1,:) = q
  pc(1,:) = p

  !Calculate forces at t = -dt, -2dt and -3dt to initialize ABM integrator
  do j = 1, 3
     call ME_predictor(qc(j,:),pc(j,:),fq(j,:),fp(j,:),qp,pp) !predict 
     call potbits(qp(1:Nnuc),U,V,dU,dV)          
     call Forces(qp,pp,dU,V,dV,fqp,fpp)           ! get forces at predicted coords
     call ME_corrector(qc(j,:),pc(j,:),fq(j,:),fp(j,:),fqp,fpp,qc(j+1,:),pc(j+1,:)) !correct prediction
     call potbits(qc(j+1,1:Nnuc),U,V,dU,dV)           
     call Forces(qc(j+1,:),pc(j+1,:),dU,V,dV,fq(j+1,:),fp(j+1,:))   ! get forces at corrected coords
     ! Inverted potential can happen here too
     if (qc(j+1,1).ge.1.d1.or.qc(j+1,1).le.0.5d0) then   ! Throw inverted pot trajs
        flagR = 1
        goto 112
     elseif (any(isnan(qc)).or.any(isnan(pc)).or.any(isnan(fq)).or.any(isnan(fp))) then
        flagR = 1
        goto 112
     end if
  end do

  do i = 1, Ntime

     !Predict next time step using ABM_predictor
     call ABM_predictor(q,p,fq,fp,qp,pp)
     !Get forces at predicted coords
     call potbits(qp(1:Nnuc),U,V,dU,dV)           
     call Forces(qp,pp,dU,V,dV,fqp,fpp)          
     ! Rearrange fs. throw f4, f3 -->f4, f2--> f3, f1 -->f2.
     ! f1 are forces at predicted coords
     fq(2:4,:) = fq(1:3,:)
     fp(2:4,:) = fp(1:3,:)
     fq(1,:)   = fqp
     fp(1,:)   = fpp

     !Correct corrdinates
     call ABM_corrector(q,p,fq,fp)
     !Get forces at corrected coords and store in f1
     call potbits(q(1:Nnuc),U,V,dU,dV)           
     call Forces(q,p,dU,V,dV,fq(1,:),fp(1,:))          
  
     ! energy conservation check
     call Hamiltonian2(q(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),H2)
     Energy = 0.5d0*dot_product(p(1:Nnuc),matmul(MassInv,p(1:Nnuc))) + H2

     ! store trajectory 
      coord(1,:,i)       = q
      coord(2,:,i)       = p

     ! check energy conservation
     if (q(1).ge.1.d1.or.q(1).le.0.5d0) then   ! Throw inverted pot trajs
        flagR = 1
        goto 112
     elseif (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance) then
        flagE = 1
        goto 112
     endif
     ! proceed to drop trajectory of not conserved

  end do
  
  
112 continue

end subroutine PropagateFwd
! ############################################################
! Calculate forces for all q and p from the eoms
! ############################################################
subroutine Forces(q,p,dU,V,dV,fq,fp)
use parameters, only    : Ndof, Nnuc, Mass, Nel
implicit none

real*8, intent(in)      :: q(Ndof), p(Ndof), dU(Nnuc)
real*8, intent(in)      :: V(Nel,Nel), dV(Nel,Nel,NNuc)
real*8, intent(out)     :: fq(Ndof), fp(Ndof)
integer                 :: j

do j = 1, Nnuc
   fq(j) = p(j)/Mass(j,j)
   fp(j) = -dU(j) - 0.5d0*(dot_product(p(Nnuc+1:Ndof),matmul(dV(:,:,j),p(Nnuc+1:Ndof)))& 
           + dot_product(q(Nnuc+1:Ndof),matmul(dV(:,:,j),q(Nnuc+1:Ndof))))
end do

fq(Nnuc+1:Ndof) = matmul(V,p(Nnuc+1:Ndof))
fp(Nnuc+1:Ndof) =-matmul(V,q(Nnuc+1:Ndof))

end subroutine Forces
! ############################################################
! Modified Euler predcitor
! ############################################################
subroutine ME_predictor(q,p,fq,fp,qp,pp)
use parameters, only    : Ndof, TimeStep, Nel, Mass
implicit none

integer                 :: i 
real*8, intent(in)      :: q(Ndof), p(Ndof)
real*8, intent(in)      :: fq(Ndof), fp(Ndof)
real*8, intent(out)     :: qp(Ndof), pp(Ndof)
real*8                  :: dt

dt = -Timestep ! propagate backwards

qp = q + fq*dt
pp = p + fp*dt
 
end subroutine ME_predictor
! ############################################################
! Modified Euler corrector
! ############################################################
subroutine ME_corrector(q,p,fq1,fp1,fq2,fp2,qp,pp)
use parameters, only    : Ndof, TimeStep, Nel, Mass
implicit none

integer                 :: i 
real*8, intent(in)      :: q(Ndof), p(Ndof)
real*8, intent(in)      :: fq1(Ndof), fp1(Ndof), fq2(Ndof), fp2(Ndof)
real*8, intent(out)     :: qp(Ndof), pp(Ndof)
real*8                  :: dt

dt = -Timestep ! propagate backwards

qp = q + 0.5d0*(fq1+fq2)*dt
pp = p + 0.5d0*(fp1+fp2)*dt
 
end subroutine ME_corrector
! ############################################################
! ABM  predcitor
! ############################################################
subroutine ABM_predictor(q,p,fq,fp,qp,pp)
use parameters, only    : Ndof, TimeStep
implicit none

real*8, intent(in)      :: q(Ndof), p(Ndof)
real*8, intent(in)      :: fq(4,Ndof), fp(4,Ndof)
real*8, intent(out)     :: qp(Ndof), pp(Ndof)

qp = q + (5.5d1*fq(1,:) - 5.9d1*fq(2,:) + 3.7d1*fq(3,:) - 9.d0*fq(4,:))/2.4d1*Timestep
pp = p + (5.5d1*fp(1,:) - 5.9d1*fp(2,:) + 3.7d1*fp(3,:) - 9.d0*fp(4,:))/2.4d1*Timestep
 
end subroutine ABM_predictor
! ############################################################
! ABM  corrector
! ############################################################
subroutine ABM_corrector(q,p,fq,fp)
use parameters, only    : Ndof, TimeStep
implicit none

real*8, intent(inout)   :: q(Ndof), p(Ndof)
real*8, intent(in)      :: fq(4,Ndof), fp(4,Ndof)

q = q + (9.d0*fq(1,:) + 1.9d1*fq(2,:) - 5.d0*fq(3,:) + 1.d0*fq(4,:))/2.4d1*Timestep
p = p + (9.d0*fp(1,:) + 1.9d1*fp(2,:) - 5.d0*fp(3,:) + 1.d0*fp(4,:))/2.4d1*Timestep
 
end subroutine ABM_corrector
! ############################################################
! ######################################## !
! ### compute elements of V as well    ### !
! ### the derivatives w.r.t. nuclear R ### ! 
! ######################################## !
subroutine potbits(R,U,V,dU,dV)
use parameters, only   : Nnuc, Nel, D0, D1, A0, A1, R0, R1, B0, Z0, D2, A2, Z1, &
                         C0, C1, Ek, NelNO, Mass, Nbath, OccN, BW, kC, Vbar, Abar, &
                         EVtoau
implicit none

real*8, intent(in)  :: R(Nnuc)
real*8, intent(out) :: U, dU(Nnuc)
real*8, intent(out) :: V(Nel,Nel), dV(Nel,Nel,Nnuc) 
real*8              :: f1, f2
real*8              :: V1, V2, dV1(Nnuc), dV2(Nnuc)
real*8              :: V12, dV12(Nnuc)
integer             :: i, No

U       = 0.d0
V       = 0.d0
dU      = 0.d0
dV      = 0.d0
No      = sum(OccN) !Total number of occupied states

V1      = 0.d0
V2      = 0.d0
V12     = 0.d0
dV1     = 0.d0
dV2     = 0.d0
dV12    = 0.d0

!Morse dof for R: V1 and V2
f1      = dexp(-A0*(R(1)-R0))
f2      = dexp(-A1*(R(1)-R1))
!Potential
V1      = D0*(f1**2-2.d0*f1)
V2      = D1*(f2**2-2.d0*f2)
!Derivatives w.r.t R
dV1(1)  = 2.d0*A0*D0*f1*(1.d0-f1)
dV2(1)  = 2.d0*A1*D1*f2*(1.d0-f2)
!Derivatives w.r.t Z
dV1(2)  = 0.d0
dV2(2)  = 0.d0

!For Z dof
f1      = dexp(-B0*(R(2)-Z0))
f2      = dexp(-A2*(R(2)-Z1)) 
! Potential
V1      = V1 + f1
V2      = V2 + D2*(f2**2-2.d0*f2)
!Derivatives w.r.t R
dV1(1)  = dV1(1) + 0.d0
dV2(1)  = dV2(1) + 0.d0
!Derivatives w.r.t Z
dV1(2)  = dV1(2) - B0*f1
dV2(2)  = dV2(2) + 2.d0*A2*D2*f2*(1.d0-f2)

!Constants for both potentials
V1      = V1 + C0
V2      = V2 + C1

! V12 coupling elements
! kC = 1.d0
V12     = Vbar*(1.d0-dtanh(R(2)/Abar))*dsqrt(BW/dble(Nbath))*kC
dV12(1) = 0.d0
dV12(2) =-Vbar/Abar/(dcosh(R(2)/Abar)**2)*dsqrt(BW/dble(Nbath))*kC

!Potential matrix elements of usual H_MMST
!1st dof is the NO- elec dof, next Nbath dofs are metal ones.
V(1,1)          = V2 - V1
V(1,NelNO+1:Nel)= V12
V(NelNO+1:Nel,1)= V(1,NelNO+1:Nel)
do i = 1, Nbath  ! sum over bath dofs
   V(NelNO+i,NelNO+i)  = Ek(i) 
end do

!Derivatives w.r.t R 
dV(1,1,1)                       = dV2(1) - dV1(1)
dV(1,NelNO+1:Nel,1)             = dV12(1)
dV(NelNO+1:Nel,1,1)             = dV(1,NelNO+1:Nel,1)
dV(NelNO+1:Nel,NelNO+1:Nel,1)   = 0.d0

!Derivative w.r.t Z
dV(1,1,2)                       = dV2(2) - dV1(2)
dV(1,NelNO+1:Nel,2)             = dV12(2)
dV(NelNO+1:Nel,1,2)             = dV(1,NelNO+1:Nel,2)
dV(NelNO+1:Nel,NelNO+1:Nel,2)   = 0.d0

! Construct symmetrized H_MMST
!Construct U and dU
do i = 1, Nel
   U     = U + V(i,i)
   dU(:) = dU(:) + dV(i,i,:) 
end do
!At this point U = Tr[V] and dU = Tr[dV]

!Construct Vbar and return in variable V
do i = 1, Nel
   V(i,i)    = V(i,i) - U/dble(Nel)
   dV(i,i,:) = dV(i,i,:) - dU(:)/dble(Nel)
end do

!Construct complete U and dU
U   = U*dble(No)/dble(Nel)
dU  = dU*dble(No)/dble(Nel)

!Ground state NO potential goes in U now
U       = U  + V1
dU(:)   = dU(:) + dV1(:)
return

end subroutine potbits
! ######################################## !
! ### compute Hamiltonian              ### !
! ###   H2                             ### !
! ######################################## !
subroutine Hamiltonian2(Rnuc,x,p,H2)
use parameters, only : Nnuc, Nel
implicit none

real*8, intent(in)  :: Rnuc(Nnuc), x(Nel), p(Nel)
real*8              :: V(Nel,Nel), dV(Nel,Nel,Nnuc)
real*8              :: U, dU(Nnuc)
real*8, intent(out) :: H2
integer             :: i

call potbits(Rnuc,U,V,dU,dV)
! H2 for symmetrized H_MMST = U(R) + 0.5*(p.Vbar.p + x.Vbar.x)

H2 = U + 0.5d0*(dot_product(p,matmul(V,p)) + dot_product(x,matmul(V,x))) 

!This term is not present in H_sym
!do i = 1, Nel
!   H2 = H2 - 0.5d0*V(i,i)
!end do

end subroutine Hamiltonian2

