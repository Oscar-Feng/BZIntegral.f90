module TetraSupply
  use QuadTetra
  implicit none
contains
! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))*Θ(X2(k))
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  recursive function QuadTetraThetaTheta(X1qtetra,X2qtetra,iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: qweights, x1qtetras, x2qtetras
    real(wp), dimension(10) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtetra)>0 .and. minval(X2qtetra)>0) then
        qw=QuadTetraTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
           0.0_wp,0.0_wp,0.0_wp/),1.0_wp,0)
      else if (maxval(X1qtetra)<0 .and. maxval(X2qtetra)<0) then
        qw=0
      else
        qw=(-QuadTetraTheta(X1qtetra*X2qtetra,0.0_wp,0)+QuadTetraTheta(&
            -X1qtetra,0.0_wp,0)+QuadTetraTheta(-X2qtetra,0.0_wp,0))/2
      end if
    else
      x1qtetras=QTetraInterpolation(X1qtetra)
      x2qtetras=QTetraInterpolation(X2qtetra)
      do i=1,8
        qweights(:,i)=QuadTetraThetaTheta(x1qtetras(:,i),x2qtetras(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraThetaTheta

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅1/D(k)
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  recursive function QuadTetraThetaThetaFrakD(X1qtetra,X2qtetra,Dqtetra,iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra, Dqtetra
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: x1qtetras, x2qtetras, dqtetras, qweights
    real(wp), dimension(10) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtetra)>0 .and. minval(X2qtetra)>0) then
        qw=QuadTetraThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
           0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtetra,0)
      else if (maxval(X1qtetra)<0 .and. maxval(X2qtetra)<0) then
        qw=0
      else
        qw=(-QuadTetraThetaFrakD(X1qtetra*X2qtetra,0.0_wp,Dqtetra,0)&
           +QuadTetraThetaFrakD(-X1qtetra,0.0_wp,Dqtetra,0)+QuadTetraThetaFrakD&
           (-X2qtetra,0.0_wp,Dqtetra,0))/2
      end if
    else
      x1qtetras=QTetraInterpolation(X1qtetra)
      x2qtetras=QTetraInterpolation(X2qtetra)
      dqtetras=QTetraInterpolation(Dqtetra)
      do i=1,8
        qweights(:,i)=QuadTetraThetaThetaFrakD(x1qtetras(:,i),x2qtetras(:,i),&
                      dqtetras(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraThetaThetaFrakD

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅δ(D(k))
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  recursive function QuadTetraThetaThetaDelta(X1qtetra,X2qtetra,Dqtetra,iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra, Dqtetra
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: x1qtetras, x2qtetras, dqtetras, qweights
    real(wp), dimension(10) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtetra)>0 .and. minval(X2qtetra)>0) then
        qw=QuadTetraThetaDelta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
           0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtetra,0)
      else if (maxval(X1qtetra)<0 .and. maxval(X2qtetra)<0) then
        qw=0
      else
        qw=(-QuadTetraThetaDelta(X1qtetra*X2qtetra,0.0_wp,Dqtetra,0)&
           +QuadTetraThetaDelta(-X1qtetra,0.0_wp,Dqtetra,0)+QuadTetraThetaDelta&
           (-X2qtetra,0.0_wp,Dqtetra,0))/2
      end if
    else
      x1qtetras=QTetraInterpolation(X1qtetra)
      x2qtetras=QTetraInterpolation(X2qtetra)
      dqtetras=QTetraInterpolation(Dqtetra)
      do i=1,8
        qweights(:,i)=QuadTetraThetaThetaDelta(x1qtetras(:,i),x2qtetras(:,i),&
                      dqtetras(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraThetaThetaDelta

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = 𝒲(X1(k))
! I think ScrW, etc. could be complex functions that take in complex input
! In order to use the functions below, PLEASE DEFINE THE FUNCTIONS ScrW, etc. FIRST
  ! recursive function QuadTetraScrW(X1qtetra,iter) result(qw)
  !   implicit none
  !   real(wp), dimension(10), intent(in) :: X1qtetra
  !   integer, intent(in) :: iter
  !   real(wp), dimension(10) :: Wqtetra
  !   real(wp), dimension(10,8) :: x1qtetras, qweights
  !   real(wp), dimension(10) :: qw
  !   integer :: i
  !
  !   if (iter==0) then
  !     do i=1,10
  !       Wqtetra(i)=ScrW(X1qtetra(i))
  !     end do
  !     if (maxval(abs(Wqtetra))==0) then
  !       qw=0
  !     else
  !       qw=QuadTetraTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
  !          0.0_wp,0.0_wp,0.0_wp/),1.0_wp,0)*Wqtetra
  !     end if
  !   else
  !     x1qtetras=QTetraInterpolation(X1qtetra)
  !     do i=1,8
  !       qweights(:,i)=QuadTetraScrW(x1qtetras(:,i),iter-1)
  !     end do
  !     qw=CollectQWeights(qweights)
  !   end if
  ! end function QuadTetraScrW

! Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))
  ! recursive function QuadTetraScrWScrW(X1qtetra,X2qtetra,iter) result(qw)
  !   implicit none
  !   real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra
  !   integer, intent(in) :: iter
  !   real(wp), dimension(10) :: W1qtetra, W2qtetra
  !   real(wp), dimension(10,8) :: x1qtetras, x2qtetras, qweights
  !   real(wp), dimension(10) :: qw
  !   integer :: i
  !
  !   if (iter==0) then
  !     do i=1,10
  !       W1qtetra(i)=ScrW1(X1qtetra(i))
  !       W2qtetra(i)=ScrW2(X2qtetra(i))
  !     end do
  !     if (maxval(abs(W1qtetra))==0 .or. maxval(abs(W2qtetra))==0) then
  !       qw=0
  !     else
  !       qw=QuadTetraTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
  !          0.0_wp,0.0_wp,0.0_wp/),1.0_wp,0)*W1qtetra*W2qtetra
  !     end if
  !   else
  !     x1qtetras=QTetraInterpolation(X1qtetra)
  !     x2qtetras=QTetraInterpolation(X2qtetra)
  !     do i=1,8
  !       qweights(:,i)=QuadTetraScrWScrW(x1qtetras(:,i),x2qtetras(:,i),iter-1)
  !     end do
  !     qw=CollectQWeights(qweights)
  !   end if
  ! end function QuadTetraScrWScrW

! Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))*𝒲3(X3(k))
  ! recursive function QuadTetraScrWScrWScrW(X1qtetra,X2qtetra,X3qtetra,iter) result(qw)
  !   implicit none
  !   real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra, X3qtetra
  !   integer, intent(in) :: iter
  !   real(wp), dimension(10) :: W1qtetra, W2qtetra, W3qtetra
  !   real(wp), dimension(10,8) :: x1qtetras, x2qtetras, x3qtetras, qweights
  !   real(wp), dimension(10) :: qw
  !   integer :: i
  !
  !   if (iter==0) then
  !     do i=1,10
  !       W1qtetra(i)=ScrW1(X1qtetra(i))
  !       W2qtetra(i)=ScrW2(X2qtetra(i))
  !       W3qtetra(i)=ScrW3(X3qtetra(i))
  !     end do
  !     if (maxval(abs(W1qtetra))==0 .or. maxval(abs(W2qtetra))==0 .or. maxval(abs(W3qtetra))==0) then
  !       qw=0
  !     else
  !       qw=QuadTetraTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
  !          0.0_wp,0.0_wp,0.0_wp/),1.0_wp,0)*W1qtetra*W2qtetra*W3qtetra
  !     end if
  !   else
  !     x1qtetras=QTetraInterpolation(X1qtetra)
  !     x2qtetras=QTetraInterpolation(X2qtetra)
  !     x3qtetras=QTetraInterpolation(X3qtetra)
  !     do i=1,8
  !       qweights(:,i)=QuadTetraScrWScrWScrW(x1qtetras(:,i),x2qtetras(:,i),x3qtetras(:,i),iter-1)
  !     end do
  !     qw=CollectQWeights(qweights)
  !   end if
  ! end function QuadTetraScrWScrWScrW

! Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is  W(k) = 𝒲(X1(k))/D(k)
  ! recursive function QuadTetraScrWFrakD(X1qtetra,Dqtetra,iter) result(qw)
  !   implicit none
  !   real(wp), dimension(10), intent(in) :: X1qtetra, Dqtetra
  !   integer, intent(in) :: iter
  !   real(wp), dimension(10) :: Wqtetra
  !   real(wp), dimension(10,8) :: x1qtetras, dqtetras, qweights
  !   real(wp), dimension(10) :: qw
  !   integer :: i
  !
  !   if (iter==0) then
  !     do i=1,10
  !       Wqtetra(i)=ScrW(X1qtetra(i))
  !     end do
  !     if (maxval(abs(Wqtetra))==0) then
  !       qw=0
  !     else
  !       qw=QuadTetraThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
  !          0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtetra,0)*Wqtetra
  !     end if
  !   else
  !     x1qtetras=QTetraInterpolation(X1qtetra)
  !     dqtetras=QTetraInterpolation(Dqtetra)
  !     do i=1,8
  !       qweights(:,i)=QuadTetraScrWFrakD(x1qtetras(:,i),dqtetras(:,i),iter-1)
  !     end do
  !     qw=CollectQWeights(qweights)
  !   end if
  ! end function QuadTetraScrWFrakD

! Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))𝒲2(X2(k))/D(k)
  ! recursive function QuadTetraScrWScrWFrakD(X1qtetra,X2qtetra,Dqtetra,iter) result(qw)
  !   implicit none
  !   real(wp), dimension(10), intent(in) :: X1qtetra, X2qtetra, Dqtetra
  !   integer, intent(in) :: iter
  !   real(wp), dimension(10) :: W1qtetra, W2qtetra
  !   real(wp), dimension(10,8) :: x1qtetras, x2qtetras, dqtetras, qweights
  !   real(wp), dimension(10) :: qw
  !   integer :: i
  !
  !   if (iter==0) then
  !     do i=1,10
  !       W1qtetra(i)=ScrW1(X1qtetra(i))
  !       W2qtetra(i)=ScrW2(X2qtetra(i))
  !     end do
  !     if (maxval(abs(W1qtetra))==0 .or. maxval(abs(W2qtetra))==0) then
  !       qw=0
  !     else
  !       qw=QuadTetraThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
  !          0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtetra,0)*W1qtetra*W2qtetra
  !     end if
  !   else
  !     x1qtetras=QTetraInterpolation(X1qtetra)
  !     x2qtetras=QTetraInterpolation(X2qtetra)
  !     dqtetras=QTetraInterpolation(Dqtetra)
  !     do i=1,8
  !       qweights(:,i)=QuadTetraScrWScrWFrakD(x1qtetras(:,i),x2qtetras(:,i),dqtetras(:,i),iter-1)
  !     end do
  !     qw=CollectQWeights(qweights)
  !   end if
  ! end function QuadTetraScrWScrWFrakD
end module TetraSupply
