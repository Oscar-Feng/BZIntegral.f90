module TrigSupply
  use QuadTrig
  implicit none
contains
! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))*Θ(X2(k))
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  recursive function QuadTrigThetaTheta(X1qtrig,X2qtrig,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig
    integer, intent(in) :: iter
    real(wp), dimension(6,4) :: qweights, x1qtrigs, x2qtrigs
    real(wp), dimension(6) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtrig)>0 .and. minval(X2qtrig)>0) then
        qw=QuadTrigTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,0)
      else if (maxval(X1qtrig)<0 .and. maxval(X2qtrig)<0) then
        qw=0
      else
        qw=(-QuadTrigTheta(X1qtrig*X2qtrig,0.0_wp,0)+QuadTrigTheta(&
            -X1qtrig,0.0_wp,0)+QuadTrigTheta(-X2qtrig,0.0_wp,0))/2
      end if
    else
      x1qtrigs=QTrigInterpolation(X1qtrig)
      x2qtrigs=QTrigInterpolation(X2qtrig)
      do i=1,4
        qweights(:,i)=QuadTrigThetaTheta(x1qtrigs(:,i),x2qtrigs(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTrigThetaTheta

! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅1/D(k)
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  recursive function QuadTrigThetaThetaFrakD(X1qtrig,X2qtrig,Dqtrig,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig, Dqtrig
    integer, intent(in) :: iter
    real(wp), dimension(6,4) :: qweights, x1qtrigs, x2qtrigs, dqtrigs
    real(wp), dimension(6) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtrig)>0 .and. minval(X2qtrig)>0) then
        qw=QuadTrigThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtrig,0)
      else if (maxval(X1qtrig)<0 .and. maxval(X2qtrig)<0) then
        qw=0
      else
        qw=(-QuadTrigThetaFrakD(X1qtrig*X2qtrig,0.0_wp,Dqtrig,0)+&
           QuadTrigThetaFrakD(-X1qtrig,0.0_wp,Dqtrig,0)+QuadTrigThetaFrakD&
           (-X2qtrig,0.0_wp,Dqtrig,0))/2
      end if
    else
      x1qtrigs=QTrigInterpolation(X1qtrig)
      x2qtrigs=QTrigInterpolation(X2qtrig)
      dqtrigs=QTrigInterpolation(Dqtrig)
      do i=1,4
        qweights(:,i)=QuadTrigThetaThetaFrakD(x1qtrigs(:,i),x2qtrigs(:,i),dqtrigs(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTrigThetaThetaFrakD

! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅δ(D(k))
! Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
  function QuadTrigThetaThetaDelta(X1qtrig,X2qtrig,Dqtrig,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig, Dqtrig
    integer, intent(in) :: iter
    real(wp), dimension(6,4) :: qweights, x1qtrigs, x2qtrigs, dqtrigs
    real(wp), dimension(6) :: qw
    integer :: i

    if (iter==0) then
      if (minval(X1qtrig)>0 .and. minval(X2qtrig)>0) then
        qw=QuadTrigThetaDelta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,Dqtrig,0)
      else if (maxval(X1qtrig)<0 .and. maxval(X2qtrig)<0) then
        qw=0
      else
        qw=(-QuadTrigThetaDelta(X1qtrig*X2qtrig,0.0_wp,Dqtrig,0)+&
           QuadTrigThetaDelta(-X1qtrig,0.0_wp,Dqtrig,0)+QuadTrigThetaDelta&
           (-X2qtrig,0.0_wp,Dqtrig,0))/2
      end if
    else
      x1qtrigs=QTrigInterpolation(X1qtrig)
      x2qtrigs=QTrigInterpolation(X2qtrig)
      dqtrigs=QTrigInterpolation(Dqtrig)
      do i=1,4
        qweights(:,i)=QuadTrigThetaThetaDelta(x1qtrigs(:,i),x2qtrigs(:,i),dqtrigs(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTrigThetaThetaDelta

! Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))
! I think SrcW, etc. could be complex functions that take in complex input
! In order to use the following functions, PLEASE DEFINE THE FUNCTIONS SrcW, etc. FIRST
!   recursive function QuadTrigSrcW(X1qtrig,iter) result(qw)
!     implicit none
!     real(wp), dimension(6), intent(in) :: X1qtrig
!     real(wp), dimension(6) :: Wqtrig
!     integer, intent(in) :: iter
!     real(wp), dimension(6) :: qw
!     real(wp), dimension(6,4) :: x1qtrigs, qweights
!     integer :: i
!
!     if (iter==0) then
!       do i=1,6
!         Wqtrig(i)=SrcW(X1qtrig(i))
!       end do
!       if (maxval(abs(Wqtrig))==0) then
!         qw=0
!       else
!         qw=QuadTrigTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,&
!            0)*Wqtrig
!       end if
!     else
!       x1qtrigs=QTrigInterpolation(X1qtrig)
!       do i=1,4
!         qweights(:,i)=QuadTrigSrcW(x1qtrigs(:,i),iter-1)
!       end do
!       qw=CollectQWeights(qweights)
!     end if
!   end function QuadtrigSrcW
!
! Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))
!   recursive function QuadTrigSrcWSrcW(X1qtrig,X2qtrig,iter) result(qw)
!     implicit none
!     real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig
!     real(wp), dimension(6) :: W1qtrig, W2qtrig
!     integer, intent(in) :: iter
!     real(wp), dimension(6) :: qw
!     real(wp), dimension(6,4) :: x1qtrigs, x2qtrigs, qweights
!     integer :: i
!
!     if (iter==0) then
!       do i=1,6
!         W1qtrig(i)=SrcW1(X1qtrig(i))
!         W2qtrig(i)=SrcW2(X2qtrig(i))
!       end do
!       if (maxval(abs(W1qtrig))==0 .or. maxval(abs(W2qtrig))==0) then
!         qw=0
!       else
!         qw=QuadTrigTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,&
!            0)*W1qtrig*W2qtrig
!       end if
!     else
!       x1qtrigs=QTrigInterpolation(X1qtrig)
!       x2qtrigs=QTrigInterpolation(X2qtrig)
!       do i=1,4
!         qweights(:,i)=QuadTrigSrcWSrcW(x1qtrigs(:,i),x2qtrigs(:,i),iter-1)
!       end do
!       qw=CollectQWeights(qweights)
!     end if
!   end function QuadtrigSrcWSrcW
!
! Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))*𝒲3(X3(k))
!   recursive function QuadTrigSrcWSrcWSrcW(X1qtrig,X2qtrig,X3qtrig,iter) result(qw)
!     implicit none
!     real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig, X3qtrig
!     real(wp), dimension(6) :: W1qtrig, W2qtrig, W3qtrig
!     integer, intent(in) :: iter
!     real(wp), dimension(6) :: qw
!     real(wp), dimension(6,4) :: x1qtrigs, x2qtrigs, x3qtrigs, qweights
!     integer :: i
!
!     if (iter==0) then
!       do i=1,6
!         W1qtrig(i)=SrcW1(X1qtrig(i))
!         W2qtrig(i)=SrcW2(X2qtrig(i))
!         W3qtrig(i)=SrcW3(X3qtrig(i))
!       end do
!       if (maxval(abs(W1qtrig))==0 .or. maxval(abs(W2qtrig))==0 .or. maxval(abs(W3qtrig))==0) then
!         qw=0
!       else
!         qw=QuadTrigTheta((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),1.0_wp,&
!            0)*W1qtrig*W2qtrig*W3qtrig
!       end if
!     else
!       x1qtrigs=QTrigInterpolation(X1qtrig)
!       x2qtrigs=QTrigInterpolation(X2qtrig)
!       x3qtrigs=QTrigInterpolation(X3qtrig)
!       do i=1,4
!         qweights(:,i)=QuadTrigSrcWSrcWSrcW(x1qtrigs(:,i),x2qtrigs(:,i),x3qtrigs(:,i),iter-1)
!       end do
!       qw=CollectQWeights(qweights)
!     end if
!   end function QuadtrigSrcWSrcWSrcW
!
! Recursive triangle rule in a single quadratic triangle,weight function is  W(k) = 𝒲(X1(k))/D(k)
!   recursive function QuadTrigSrcWFrakD(X1qtrig,Dqtrig,iter) result(qw)
!     implicit none
!     real(wp), dimension(6), intent(in) :: X1qtrig, Dqtrig
!     real(wp), dimension(6) :: Wqtrig
!     integer, intent(in) :: iter
!     real(wp), dimension(6) :: qw
!     real(wp), dimension(6,4) :: x1qtrigs, dqtrigs, qweights
!     integer :: i
!
!     if (iter==0) then
!       do i=1,6
!         Wqtrig(i)=SrcW(X1qtrig(i))
!       end do
!       if (maxval(abs(Wqtrig))==0) then
!         qw=0
!       else
!         qw=QuadTrigThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),&
!            1.0_wp,Dqtrig,0)*Wqtrig
!       end if
!     else
!       x1qtrigs=QTrigInterpolation(X1qtrig)
!       dqtrigs=QTrigInterpolation(Dqtrig)
!       do i=1,4
!         qweights(:,i)=QuadTrigSrcWFrakD(x1qtrigs(:,i),dqtrigs(:,i),iter-1)
!       end do
!       qw=CollectQWeights(qweights)
!     end if
!   end function QuadtrigSrcWFrakD
!
! Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))𝒲2(X2(k))/D(k)
!   recursive function QuadTrigSrcWSrcWFrakD(X1qtrig,X2qtrig,Dqtrig,iter) result(qw)
!     implicit none
!     real(wp), dimension(6), intent(in) :: X1qtrig, X2qtrig, Dqtrig
!     real(wp), dimension(6) :: W1qtrig, W2qtrig
!     integer, intent(in) :: iter
!     real(wp), dimension(6) :: qw
!     real(wp), dimension(6,4) :: x1qtrigs, x2qtrigs, dqtrigs, qweights
!     integer :: i
!
!     if (iter==0) then
!       do i=1,6
!         W1qtrig(i)=SrcW1(X1qtrig(i))
!         W2qtrig(i)=SrcW2(X2qtrig(i))
!       end do
!       if (maxval(abs(W1qtrig))==0 .or. maxval(abs(W2qtrig))==0) then
!         qw=0
!       else
!         qw=QuadTrigThetaFrakD((/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/),&
!            1.0_wp,Dqtrig,0)*W1qtrig*W2qtrig
!       end if
!     else
!       x1qtrigs=QTrigInterpolation(X1qtrig)
!       x2qtrigs=QTrigInterpolation(X2qtrig)
!       dqtrigs=QTrigInterpolation(Dqtrig)
!       do i=1,4
!         qweights(:,i)=QuadTrigSrcWSrcWFrakD(x1qtrigs(:,i),x2qtrigs(:,i),dqtrigs(:,i),iter-1)
!       end do
!       qw=CollectQWeights(qweights)
!     end if
!   end function QuadtrigSrcWSrcWFrakD
end module TrigSupply
