!--------------------------------------------------------------------- 
!---------------------------------------------------------------------
!    Spatiotemporal stability of periodic travelling waves
!         in a heteroclinic-cycle model 
! --------------------------------------------------------------------
! TWstability: Essential spectrum and Eckhaus stability boundary of
! large-wavelength travelling waves
!
! --------------------------------------------------------------------
! Ref.: Hasan, Osinga, Postlethwaite and Rucklidge,
!       arXiv 1911.10447 (revised version from November 2020)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

  SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION  Duf(3,3),Rlam,Ilam,phi,mu,Dlam,lamp,lampp
      DOUBLE PRECISION  Az, Bz, Cz, rho, gam, sig, zeta, L
      DOUBLE PRECISION  Uvar(3), Uz(3), ReV(3), ReVz(3), ImV(3), ImVz(3)
      DOUBLE PRECISION  Vp(3), Vpz(3), Vpp(3), Vppz(3)
      DOUBLE PRECISION  a, b, c
      INTEGER k

      !Variables in original coordinates
      a=EXP(U(1))
      b=EXP(U(2))
      c=EXP(U(3))
      rho=a+b+c
      !Derivatives of variables in logarithmic coordinates     
      Az=U(4)
      Bz=U(5)
      Cz=U(6)
      !wavespeed
      gam=PAR(1)
      !system parameters
      sig=PAR(2)
      zeta=PAR(3)

      !The imaginary part of the spatial Floquet exponent
      phi=PAR(4)
      !The real part of the spatial Floquet exponent (always kept at mu=0)
      mu=PAR(5)
      !Real and imaginary parts of lambda
      Rlam=PAR(6) 
      Ilam=PAR(7)
      !wavelength
      L=PAR(13)
      !first and second derivatives of lambda with respect to nu=i*phi
      lamp=PAR(16)
      lampp=PAR(17)

      !Variable vector Uvar (can't use U) and its derivative Uz
      Uvar=U(1:3)
      Uz=U(4:6)
      !Eigenfunction V and its derivative Vz
      ReV=U(7:9)
      ReVz=U(10:12)
      ImV=U(13:15)
      ImVz=U(16:18)
      !First and second derivatives of the eigenfunction V with respect to nu=i*phi
      Vp=U(19:21)
      Vpz=U(22:24)
      Vpp=U(25:27)
      Vppz=U(28:30)

      !Travelling-frame equations
      F(1:3)=L*Uz
      F(4)=  L*(gam*Az- (1-rho-(sig+zeta)*b+ zeta*c) - Az**2)
      F(5)=  L*(gam*Bz- (1-rho-(sig+zeta)*c+ zeta*a) - Bz**2 )
      F(6)=  L*(gam*Cz- (1-rho-(sig+zeta)*a+ zeta*b) - Cz**2 )
      IF (NDIM==6) RETURN
      
      !Jacobian matrix
      Duf(1,1) = -a
      Duf(1,2) = -(1+sig+zeta)*b
      Duf(1,3) = -(1-zeta)*c

      Duf(2,1) = -(1-zeta)*a
      Duf(2,2) = -b
      Duf(2,3) = -(1+sig+zeta)*c
      
      Duf(3,1) = -(1+sig+zeta)*a
      Duf(3,2) = -(1-zeta)*b
      Duf(3,3) = -c
 
      !Eigenvalue problem
      !Equations for the real part of the eigenfunction    
      F(7:9)   = L*(ReVz-(mu*ReV-phi*ImV))
      F(10:12) = L*(Rlam*ReV -Ilam*ImV +gam*ReVz -MATMUL(Duf,ReV) -2*Uz*ReVz-(mu*ReVz-phi*ImVz))
      !Equations for the imaginary part of the eigenfunction    
      F(13:15) = L*(ImVz-(mu*ImV+phi*ReV))
      F(16:18) = L*(Rlam*ImV +Ilam*ReV +gam*ImVz -MATMUL(Duf,ImV) -2*Uz*ImVz-(mu*ImVz+phi*ReVz))
      IF (NDIM==18) RETURN

      !Equations for the first and second derivtives wrt to phi
      F(19:21)= L*(Vpz-ReV) 
      F(22:24)= L*((gam-2*Uz)*Vpz - MATMUL(Duf,Vp)- ReVz + lamp*ReV)
      F(25:27)= L*(Vppz-2*Vp)  
      F(28:30)= L*((gam-2*Uz)*Vppz- MATMUL(Duf,Vpp)-2*Vpz+ lampp*ReV+ 2*lamp*Vp )    
      IF (NDIM==30) RETURN
      
      RETURN
      END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      
      SUBROUTINE STPNT(NDIM,U,PAR)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION  sig

      !Initialising parameters
      PAR(1)= 0.1d00 !gamma
      PAR(2)= 3.2d00 !sigma
      PAR(3)= 1.0d00 !zeta
      PAR(4)= 0.0d00 !phi
      PAR(5)= 0.0d00 !mu
      PAR(6)= 0.0d00 !Rlam
      PAR(7)= 0.0d00 !Ilam
      PAR(12)= 0.0d00 !h
      PAR(13)= 1.0d00 !L
      PAR(16)= 0.0d00!lamp
      PAR(17)= 0.0d00!lampp

      sig=PAR(2)
      !Starting from the coexistence steady-state solution
      U(1:3)=LOG(1/(3+sig))
      U(4:6)=0.0d0
      IF (NDIM==6) RETURN
      !Starting form a nonzero constant eigenfunction
      U(7)= 1.0d0
      U(8:18)=0.0d0
      IF (NDIM==18) RETURN
      !Starting from a trivial solution for the derivatives of the eigenfunction
      U(19:30)=0.0d0
      IF (NDIM==30) RETURN
      
      RETURN
      END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

      !Imposing periodic boundary conditions  
      FB(1:6) = U1(1:6) - U0(1:6)
      IF (NBC==6) RETURN
      FB(7:18)= U1(7:18)- U0(7:18) 
      IF (NBC==18) RETURN
      FB(19:30) = U1(19:30) - U0(19:30)

      RETURN 
      END 

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
     
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
      DOUBLE PRECISION h

      !Norm of the eigenfunctions
      h=PAR(12)

      !Phase condition for the periodic travelling wave
      FI(1)=DOT_PRODUCT(UPOLD(1:6),(UOLD(1:6)-U(1:6)))
      !Normalizing condition for the eigenfunction
      FI(2)=DOT_PRODUCT(U(7:18),U(7:18)) - h
      IF (NINT==2) RETURN
      !Phase condition for the eigenfunction
      FI(3) =  DOT_PRODUCT(UOLD(7:12),U(13:18)) - DOT_PRODUCT(UOLD(13:18),U(7:12))
      IF (NINT==3) RETURN
      !Integral conditions necassry for computing the instability curve
      FI(4)=DOT_PRODUCT(U(7:9),U(19:21))
      FI(5)=DOT_PRODUCT(U(7:9),U(25:27))
     
      RETURN 
      END 

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FOPT 
      RETURN 
      END 

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PVLS
      RETURN 
      END 

! ----------------------------------------------------------------------
      
