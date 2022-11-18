module BZInt2D
  use OMP_LIB
  use SplitMesh
  use LinTrig
  use TrigSupply
  implicit none
contains
! Linear triangle method for weight function W(k) = Θ(eF-E(k))
! with Blöchl correction
  subroutine Lin2DRuleTheta(Emesh,eF,Wmesh)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, allocatable, dimension(:,:,:) :: Trigs
    real(wp), allocatable, dimension(:,:) :: ETrigs, WTrigs
    integer :: dim, i, j

    call Mesh2Trig(size(Emesh,dim=1),size(Emesh,dim=2),Trigs)

    dim=size(Trigs,dim=3)
    allocate(ETrigs(3,dim))
    do i=1,dim
      do j=1,3
        ETrigs(j,i)=Emesh(Trigs(1,j,i),Trigs(2,j,i))
      end do
    end do

    allocate(WTrigs(3,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WTrigs(:,i)=LinTrigTheta_Blochl(ETrigs(:,i),eF)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,3
        Wmesh(Trigs(1,j,i),Trigs(2,j,i))=WTrigs(j,i)+Wmesh(Trigs(1,j,i),Trigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(Trigs,ETrigs,WTrigs)
  end subroutine Lin2DRuleTheta

! Linear triangle method for weight function W(k) = δ(eF-E(k))
  subroutine Lin2DRuleDelta(Emesh,eF,Wmesh)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, allocatable, dimension(:,:,:) :: Trigs
    real(wp), allocatable, dimension(:,:) :: ETrigs, WTrigs
    integer :: dim, i, j

    call Mesh2Trig(size(Emesh,dim=1),size(Emesh,dim=2),Trigs)

    dim=size(Trigs,dim=3)
    allocate(ETrigs(3,dim))
    do i=1,dim
      do j=1,3
        ETrigs(j,i)=Emesh(Trigs(1,j,i),Trigs(2,j,i))
      end do
    end do

    allocate(WTrigs(3,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WTrigs(:,i)=LinTrigDelta(ETrigs(:,i),eF)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,3
        Wmesh(Trigs(1,j,i),Trigs(2,j,i))=WTrigs(j,i)+Wmesh(Trigs(1,j,i),Trigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(Trigs,ETrigs,WTrigs)
  end subroutine Lin2DRuleDelta

! Recursive triangle method for weight function W(k) = Θ(eF-E(k))
  subroutine Quad2DRuleTheta(Emesh,eF,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: EQTrigs, WQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(Emesh,dim=1),size(Emesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(EQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        EQTrigs(j,i)=Emesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigTheta(EQTrigs(:,i),eF,iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,EQTrigs,WQTrigs)
  end subroutine Quad2DRuleTheta

! Recursive triangle method for weight function W(k) = W(k) = 1/D(k) Θ(eF-E(k))
  subroutine Quad2DRuleThetaFrakD(Emesh,eF,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh,Dmesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: EQTrigs, WQTrigs, DQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(Emesh,dim=1),size(Emesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(EQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        EQTrigs(j,i)=Emesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(DQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigThetaFrakD(EQTrigs(:,i),eF,DQTrigs(:,i),iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,EQTrigs,DQTrigs,WQTrigs)
  end subroutine Quad2DRuleThetaFrakD

! Recursive triangle method for weight function W(k) = δ(eF-E(k))
  subroutine Quad2DRuleDelta(Emesh,eF,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: EQTrigs, WQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(Emesh,dim=1),size(Emesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(EQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        EQTrigs(j,i)=Emesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigDelta(EQTrigs(:,i),eF,iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,EQTrigs,WQTrigs)
  end subroutine Quad2DRuleDelta

! Recursive triangle method for weight function W(k) = δ(D(k)) Θ(eF-E(k))
  subroutine Quad2DRuleThetaDelta(Emesh,eF,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: Emesh,Dmesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    real(wp), intent(in) :: eF
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: EQTrigs, WQTrigs, DQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(Emesh,dim=1),size(Emesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(EQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        EQTrigs(j,i)=Emesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(DQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigThetaDelta(EQTrigs(:,i),eF,DQTrigs(:,i),iter)
    end do

    allocate(Wmesh(size(Emesh,dim=1),size(Emesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,EQTrigs,DQTrigs,WQTrigs)
  end subroutine Quad2DRuleThetaDelta

! Recursive triangle method for weight function W(k) = Θ(x1(k))*Θ(x2(k))
  subroutine Quad2DRuleThetaTheta(X1mesh,X2mesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh,X2mesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, WQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(X1QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(X2QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigThetaTheta(X1QTrigs(:,i),X2QTrigs(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,X1QTrigs,X2QTrigs,WQTrigs)
  end subroutine Quad2DRuleThetaTheta

! Recursive triangle method for weight function W(k) = 1/D(k) Θ(x1(k))*Θ(x2(k))
  subroutine Quad2DRuleThetaThetaFrakD(X1mesh,X2mesh,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh,X2mesh,Dmesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, Dqtrigs, WQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(X1QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(X2QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(DQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigThetaThetaFrakD(X1QTrigs(:,i),X2QTrigs(:,i),Dqtrigs(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,X1QTrigs,X2QTrigs,WQTrigs,Dqtrigs)
  end subroutine Quad2DRuleThetaThetaFrakD

! Recursive triangle method for weight function W(k) = δ(D(k)) Θ(x1(k))*Θ(x2(k))
  subroutine Quad2DRuleThetaThetaDelta(X1mesh,X2mesh,Dmesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh,X2mesh,Dmesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    integer, allocatable, dimension(:,:,:) :: QTrigs
    real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, Dqtrigs, WQTrigs
    integer :: dim, i, j

    if (.not. present(iter)) then
      iter=2
    end if

    call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)

    dim=size(QTrigs,dim=3)
    allocate(X1QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(X2QTrigs(6,dim))
    do i=1,dim
      do j=1,6
        X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(DQTrigs(6,dim))
    do i=1,dim
      do j=1,6
        DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do

    allocate(WQTrigs(6,dim))
!$OMP PARALLEL DO
    do i=1,dim
      WQTrigs(:,i)=QuadTrigThetaThetaDelta(X1QTrigs(:,i),X2QTrigs(:,i),Dqtrigs(:,i),iter)
    end do

    allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
    Wmesh=0
    do i=1,dim
      do j=1,6
        Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
      end do
    end do
    Wmesh=Wmesh/dim*2

    deallocate(QTrigs,X1QTrigs,X2QTrigs,WQTrigs,Dqtrigs)
  end subroutine Quad2DRuleThetaThetaDelta

! Recursive triangle method for weight function W(k) = δ(x1(k))*δ(x2(k))
  subroutine Quad2DRuleDeltaDelta(X1mesh,X2mesh,Wmesh,iter)
    implicit none
    real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh, X2mesh
    real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
    integer, intent(inout), optional :: iter
    real(wp), allocatable, dimension(:,:) :: Wmesh0, Wmeshd
    real(wp) :: dx1

    if (.not. present(iter)) then
      iter=2
    end if

    dx1=maxval(abs(X1mesh))*0.0005_wp
    call Quad2DRuleThetaDelta(X1mesh,0.0_wp,X2mesh,Wmesh0,iter)
    call Quad2DRuleThetaDelta(X1mesh,dx1,X2mesh,Wmeshd,iter)

    allocate(Wmesh(size(Wmesh0,dim=1),size(Wmesh0,dim=2)))
    Wmesh=(Wmeshd-Wmesh0)/dx1

    deallocate(Wmesh0,Wmeshd)
  end subroutine Quad2DRuleDeltaDelta

! Recursive triangle method for weight function W(k) = 𝒲(x1(k))
! In order to use the following functions, PLEASE DEFINE THE FUNCTIONS SrcW, etc. FIRST
!   subroutine Quad2DRuleSrcW(X1mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh
!     real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTrigs
!     real(wp), allocatable, dimension(:,:) :: X1QTrigs, WQTrigs
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)
!
!     dim=size(QTrigs,dim=3)
!     allocate(X1QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(WQTrigs(6,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTrigs(:,i)=QuadTrigSrcW(X1QTrigs(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,6
!         Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!     Wmesh=Wmesh/dim*2
!
!     deallocate(QTrigs,X1QTrigs,WQTrigs)
!   end subroutine Quad2DRuleSrcW
!
! Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))
!   subroutine Quad2DRuleSrcWSrcW(X1mesh,X2mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh, X2mesh
!     real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTrigs
!     real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, WQTrigs
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)
!
!     dim=size(QTrigs,dim=3)
!     allocate(X1QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(X2QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(WQTrigs(6,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTrigs(:,i)=QuadTrigSrcWSrcW(X1QTrigs(:,i),X2QTrigs(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,6
!         Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!     Wmesh=Wmesh/dim*2
!
!     deallocate(QTrigs,X1QTrigs,X2QTrigs,WQTrigs)
!   end subroutine Quad2DRuleSrcWSrcW
!
! Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))*𝒲3(x3(k))
!   subroutine Quad2DRuleSrcWSrcWSrcW(X1mesh,X2mesh,X3mesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh, X2mesh, X3mesh
!     real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTrigs
!     real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, X3QTrigs, WQTrigs
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)
!
!     dim=size(QTrigs,dim=3)
!     allocate(X1QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(X2QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(X3QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X3QTrigs(j,i)=X3mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(WQTrigs(6,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTrigs(:,i)=QuadTrigSrcWSrcWSrcW(X1QTrigs(:,i),X2QTrigs(:,i),X3QTrigs(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,6
!         Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!     Wmesh=Wmesh/dim*2
!
!     deallocate(QTrigs,X1QTrigs,X2QTrigs,X3QTrigs,WQTrigs)
!   end subroutine Quad2DRuleSrcWSrcWSrcW
!
! Recursive triangle method for weight function W(k) = 𝒲(x1(k))/ D(k)
!   subroutine Quad2DRuleSrcWFrakD(X1mesh,Dmesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh, Dmesh
!     real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTrigs
!     real(wp), allocatable, dimension(:,:) :: X1QTrigs, DQTrigs, WQTrigs
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)
!
!     dim=size(QTrigs,dim=3)
!     allocate(X1QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(DQTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(WQTrigs(6,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTrigs(:,i)=QuadTrigSrcWFrakD(X1QTrigs(:,i),DQTrigs(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,6
!         Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!     Wmesh=Wmesh/dim*2
!
!     deallocate(QTrigs,X1QTrigs,DQTrigs,WQTrigs)
!   end subroutine Quad2DRuleSrcWFrakD
!
! Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))/ D(k)
!   subroutine Quad2DRuleSrcWSrcWFrakD(X1mesh,X2mesh,Dmesh,Wmesh,iter)
!     implicit none
!     real(wp), allocatable, dimension(:,:), intent(in) :: X1mesh, X2mesh, Dmesh
!     real(wp), allocatable, dimension(:,:), intent(out) :: Wmesh
!     integer, intent(inout), optional :: iter
!     integer, allocatable, dimension(:,:,:) :: QTrigs
!     real(wp), allocatable, dimension(:,:) :: X1QTrigs, X2QTrigs, DQTrigs, WQTrigs
!     integer :: dim, i, j
!
!     if (.not. present(iter)) then
!       iter=2
!     end if
!
!     call Mesh2QTrig(size(X1mesh,dim=1),size(X1mesh,dim=2),QTrigs)
!
!     dim=size(QTrigs,dim=3)
!     allocate(X1QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X1QTrigs(j,i)=X1mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(X2QTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         X2QTrigs(j,i)=X2mesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(DQTrigs(6,dim))
!     do i=1,dim
!       do j=1,6
!         DQTrigs(j,i)=Dmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!
!     allocate(WQTrigs(6,dim))
! !$OMP PARALLEL DO
!     do i=1,dim
!       WQTrigs(:,i)=QuadTrigSrcWSrcWFrakD(X1QTrigs(:,i),X2QTrigs(:,i),DQTrigs(:,i),iter)
!     end do
!
!     allocate(Wmesh(size(X1mesh,dim=1),size(X1mesh,dim=2)))
!     Wmesh=0
!     do i=1,dim
!       do j=1,6
!         Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))=WQTrigs(j,i)+Wmesh(QTrigs(1,j,i),QTrigs(2,j,i))
!       end do
!     end do
!     Wmesh=Wmesh/dim*2
!
!     deallocate(QTrigs,X1QTrigs,X2QTrigs,DQTrigs,WQTrigs)
!   end subroutine Quad2DRuleSrcWSrcWFrakD
end module BZInt2D
