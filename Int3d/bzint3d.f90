module BZInt3D
  use OMP_LIB
  use SplitMesh
  use LinTetra
  use TetraSupply
  implicit none
contains
! Linear tetrahedron method for weight function W(k) = Θ(eF-E(k))
! with Blöchl correction
  subroutine Lin3DRuleTheta(Emesh,eF,Wmesh)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, allocatable, dimension(:,:,:) :: Tetras
    real(wp), allocatable, dimension(:,:) :: ETetras, WTetras
    integer :: dim, i, j

    call Mesh2Tetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),Tetras)

    dim=size(Tetras,dim=3)
    allocate(ETetras(4,dim))
    do i=1,dim
      do j=1,4
        ETetras(j,i)=Emesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))
      end do
    end do

    allocate(WTetras(4,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WTetras(:,i)=LinTetraTheta_Blochl(ETetras(:,i),eF)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,4
        Wmesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))&
          =Wmesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))+WTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6
  end subroutine Lin3DRuleTheta

! Linear tetrahedron method for weight function W(k) = δ(eF-E(k))
  subroutine Lin3DRuleDelta(Emesh,eF,Wmesh)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, allocatable, dimension(:,:,:) :: Tetras
    real(wp), allocatable, dimension(:,:) :: ETetras, WTetras
    integer :: dim, i, j

    call Mesh2Tetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),Tetras)

    dim=size(Tetras,dim=3)
    allocate(ETetras(4,dim))
    do i=1,dim
      do j=1,4
        ETetras(j,i)=Emesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))
      end do
    end do

    allocate(WTetras(4,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WTetras(:,i)=LinTetraDelta(ETetras(:,i),eF)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,4
        Wmesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))&
          =Wmesh(Tetras(1,j,i),Tetras(2,j,i),Tetras(3,j,i))+WTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6
  end subroutine Lin3DRuleDelta

! Recursive tetrahedron method for weight function W(k) = Θ(eF-E(k))
  subroutine Quad3DRuleTheta(Emesh,eF,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: EQTetras, WQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(EQTetras(10,dim))
    do i=1,dim
      do j=1,10
        EQTetras(j,i)=Emesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraTheta(EQTetras(:,i),eF,iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,EQTetras,WQTetras)
  end subroutine Quad3DRuleTheta

! Recursive tetrahedron method for weight function W(k) = W(k) = 1/D(k) Θ(eF-E(k))
  subroutine Quad3DRuleThetaFrakD(Emesh,eF,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh,Dmesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: EQTetras, WQTetras, DQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(EQTetras(10,dim))
    do i=1,dim
      do j=1,10
        EQTetras(j,i)=Emesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(DQTetras(10,dim))
    do i=1,dim
      do j=1,10
        DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraThetaFrakD(EQTetras(:,i),eF,DQTetras(:,i),iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,EQTetras,DQTetras,WQTetras)
  end subroutine Quad3DRuleThetaFrakD

! Recursive tetrahedron method for weight function W(k) = δ(eF-E(k))
  subroutine Quad3DRuleDelta(Emesh,eF,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: EQTetras, WQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(EQTetras(10,dim))
    do i=1,dim
      do j=1,10
        EQTetras(j,i)=Emesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraDelta(EQTetras(:,i),eF,iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,EQTetras,WQTetras)
  end subroutine Quad3DRuleDelta

! Recursive tetrahedron method for weight function W(k) = δ(D(k)) Θ(eF-E(k))
  subroutine Quad3DRuleThetaDelta(Emesh,eF,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: Emesh,Dmesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: EQTetras, WQTetras, DQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(EQTetras(10,dim))
    do i=1,dim
      do j=1,10
        EQTetras(j,i)=Emesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(DQTetras(10,dim))
    do i=1,dim
      do j=1,10
        DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraThetaDelta(EQTetras(:,i),eF,DQTetras(:,i),iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2),size(Emesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,EQTetras,DQTetras,WQTetras)
  end subroutine Quad3DRuleThetaDelta

! Recursive tetrahedron method for weight function W(k) = Θ(x1(k))*Θ(x2(k))
  subroutine Quad3DRuleThetaTheta(X1mesh,X2mesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, WQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(X1QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(X2QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraThetaTheta(X1QTetras(:,i),X2QTetras(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,X1QTetras,X2QTetras,WQTetras)
  end subroutine Quad3DRuleThetaTheta

! Recursive tetrahedron method for weight function W(k) = 1/D(k) Θ(x1(k))*Θ(x2(k))
  subroutine Quad3DRuleThetaThetaFrakD(X1mesh,X2mesh,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh, Dmesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, DQtetras, WQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(X1QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(X2QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(DQTetras(10,dim))
    do i=1,dim
      do j=1,10
        DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraThetaThetaFrakD(X1QTetras(:,i),X2QTetras(:,i),DQTetras(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,X1QTetras,X2QTetras,DQTetras,WQTetras)
  end subroutine Quad3DRuleThetaThetaFrakD

! Recursive tetrahedron method for weight function W(k) = δ(D(k)) Θ(x1(k))*Θ(x2(k))
  subroutine Quad3DRuleThetaThetaDelta(X1mesh,X2mesh,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh, Dmesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTetras
    real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, DQtetras, WQTetras
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)

    dim=size(QTetras,dim=3)
    allocate(X1QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(X2QTetras(10,dim))
    do i=1,dim
      do j=1,10
        X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(DQTetras(10,dim))
    do i=1,dim
      do j=1,10
        DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
      end do
    end do

    allocate(WQTetras(10,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTetras(:,i)=QuadTetraThetaThetaDelta(X1QTetras(:,i),X2QTetras(:,i),DQTetras(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
    Wmesh=0
    do i=1,dim
      do j=1,10
        Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
          =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
      end do
    end do
    Wmesh=Wmesh/dim*6

    deallocate(QTetras,X1QTetras,X2QTetras,DQTetras,WQTetras)
  end subroutine Quad3DRuleThetaThetaDelta

! Recursive tetrahedron method for weight function W(k) = δ(x1(k))*δ(x2(k))
  subroutine Quad3DRuleDeltaDelta(X1mesh,X2mesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh
    real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    real(wp), allocatable, dimension(:,:,:) :: Wmesh0, Wmeshd
    real(wp) :: dx1

    if (.not. present(iter)) then
      iter=2
    end if

    dx1=maxval(abs(X1mesh))*0.0005_wp
    call Quad3DRuleThetaDelta(X1mesh,0.0_wp,X2mesh,Wmesh0,iter)
    call Quad3DRuleThetaDelta(X1mesh,dx1,X2mesh,Wmeshd,iter)

    allocate(Wmesh(size(Wmesh0,dim=1),size(Wmesh0,dim=2),size(Wmesh0,dim=3)))
    Wmesh=(Wmeshd-Wmesh0)/dx1

    deallocate(Wmesh0,Wmeshd)
  end subroutine Quad3DRuleDeltaDelta

! Recursive tetrahedron method for weight function W(k) = 𝒲(x1(k))
! In order to use the functions below, PLEASE DEFINE THE FUNCTIONS ScrW, etc. FIRST
!   subroutine Quad3DRuleSrcW(X1mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh
!     real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTetras
!     real(wp), allocatable, dimension(:,:) :: X1QTetras, WQTetras
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)
!
!     dim=size(QTetras,dim=3)
!     allocate(X1QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(WQTetras(10,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTetras(:,i)=QuadTetraSrcW(X1QTetras(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,10
!         Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
!           =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
!       end do
!     end do
!     Wmesh=Wmesh/dim*6
!     deallocate(QTetras,X1QTetras,WQTetras)
!   end subroutine Quad3DRuleSrcW
!
! Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))
!   subroutine Quad3DRuleSrcWSrcW(X1mesh,X2mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh
!     real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTetras
!     real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, WQTetras
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)
!
!     dim=size(QTetras,dim=3)
!     allocate(X1QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(X2QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(WQTetras(10,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTetras(:,i)=QuadTetraSrcWSrcW(X1QTetras(:,i),X2QTetras(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,10
!         Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
!           =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
!       end do
!     end do
!     Wmesh=Wmesh/dim*6
!     deallocate(QTetras,X1QTetras,X2QTetras,WQTetras)
!   end subroutine Quad3DRuleSrcWSrcW

! Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))*𝒲3(x3(k))
!   subroutine Quad3DRuleSrcWSrcWSrcW(X1mesh,X2mesh,X3mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh, X3mesh
!     real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTetras
!     real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, X3QTetras, WQTetras
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)
!
!     dim=size(QTetras,dim=3)
!     allocate(X1QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(X2QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(X3QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X3QTetras(j,i)=X3mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(WQTetras(10,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTetras(:,i)=QuadTetraSrcWSrcWSrcW(X1QTetras(:,i),X2QTetras(:,i),X3QTetras(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,10
!         Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
!           =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
!       end do
!     end do
!     Wmesh=Wmesh/dim*6
!     deallocate(QTetras,X1QTetras,X2QTetras,X3QTetras,WQTetras)
!   end subroutine Quad3DRuleSrcWSrcWSrcW

! Recursive tetrahedron method for weight function W(k) = 𝒲(x1(k))/D(k)
!   subroutine Quad3DRuleSrcWFrakD(X1mesh,Dmesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, Dmesh
!     real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTetras
!     real(wp), allocatable, dimension(:,:) :: X1QTetras, DQTetras, WQTetras
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)
!
!     dim=size(QTetras,dim=3)
!     allocate(X1QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(DQTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(WQTetras(10,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTetras(:,i)=QuadTetraSrcWFrakD(X1QTetras(:,i),DQTetras(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,10
!         Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
!           =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
!       end do
!     end do
!     Wmesh=Wmesh/dim*6
!     deallocate(QTetras,X1QTetras,DQTetras,WQTetras)
!   end subroutine Quad3DRuleSrcWFrakD
!
! Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))/D(k)
!   subroutine Quad3DRuleSrcWSrcWFrakD(X1mesh,X2mesh,Dmesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:,:), intent(in) :: X1mesh, X2mesh, Dmesh
!     real(wp), allocatable, dimension(:,:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTetras
!     real(wp), allocatable, dimension(:,:) :: X1QTetras, X2QTetras, DQTetras, WQTetras
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTetra(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3),QTetras)
!
!     dim=size(QTetras,dim=3)
!     allocate(X1QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X1QTetras(j,i)=X1mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(X2QTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         X2QTetras(j,i)=X2mesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(DQTetras(10,dim))
!     do i=1,dim
!       do j=1,10
!         DQTetras(j,i)=Dmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))
!       end do
!     end do
!
!     allocate(WQTetras(10,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTetras(:,i)=QuadTetraSrcWSrcWFrakD(X1QTetras(:,i),X2QTetras(:,i),DQTetras(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2),size(X1mesh,dim=3)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,10
!         Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))&
!           =Wmesh(QTetras(1,j,i),QTetras(2,j,i),QTetras(3,j,i))+WQTetras(j,i)
!       end do
!     end do
!     Wmesh=Wmesh/dim*6
!     deallocate(QTetras,X1QTetras,X2QTetras,DQTetras,WQTetras)
!   end subroutine Quad3DRuleSrcWSrcWFrakD
end module BZInt3D
