module QuadTrig
  use LinTrig
  implicit none

  integer, dimension(3,6), parameter :: InterpEq1=reshape((/4, 1, 2,  4, 2, 1,&
                    6, 2, 3,  6, 3, 2,  5, 3, 1,  5, 1, 3/),(/3,6/))
  integer, dimension(5,3), parameter :: InterpEq2=reshape((/4, 5, 6, 2, 3,&
                    4, 6, 5, 1, 3,  5, 6, 4, 1, 2/),(/5,3/))
  integer, dimension(6,4), parameter :: subqtrig2qtrig=reshape((/1,4,5, 7,12,13,&
                    5,6,3, 15,11,10,  4,2,6, 8,14,9,  6,5,4, 15,14,13/),(/6,4/))
  integer, dimension(3,4), parameter :: qtrig2trig=reshape((/1, 4, 5,  4, 2, 6,&
                    5, 6, 3,  6, 5, 4/),(/3,4/))
  real(wp), dimension(15,6), parameter :: subqW2qW=reshape((/0.25_wp, 0.0_wp,&
   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.09375_wp, -0.03125_wp, 0.0_wp, 0.0_wp, &
   -0.03125_wp, 0.09375_wp, 0.0_wp, -0.03125_wp, -0.03125_wp, 0.0_wp, 0.25_wp,&
   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -0.03125_wp, 0.09375_wp, 0.09375_wp,&
   -0.03125_wp, 0.0_wp, 0.0_wp, -0.03125_wp, 0.0_wp, -0.03125_wp, 0.0_wp,&
   0.0_wp, 0.25_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -0.03125_wp, &
   0.09375_wp, 0.09375_wp, -0.03125_wp, -0.03125_wp, -0.03125_wp, 0.0_wp, &
   0.0_wp, 0.0_wp, 0.0_wp, 0.25_wp, 0.0_wp, 0.0_wp, 0.1875_wp, 0.1875_wp, &
   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.125_wp, 0.125_wp, 0.0625_wp, 0.0_wp, &
   0.0_wp, 0.0_wp, 0.0_wp, 0.25_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,&
   0.1875_wp, 0.1875_wp, 0.125_wp, 0.0625_wp, 0.125_wp, 0.0_wp, 0.0_wp, 0.0_wp,&
   0.0_wp, 0.0_wp, 0.25_wp, 0.0_wp, 0.0_wp, 0.1875_wp, 0.1875_wp, 0.0_wp, &
   0.0_wp, 0.0625_wp, 0.125_wp, 0.125_wp/),(/15,6/))
contains
  function QTrigInterpolation(qdata) result(subqdata_qtrig)
    implicit none
    real(wp), dimension(6), intent(in) :: qdata
    real(wp), dimension(6,4) :: subqdata_qtrig
    real(wp), dimension(15) :: rawdata
    integer :: i

    rawdata=0
    do i=1,6
      rawdata(i)=qdata(i)
    end do

    do i=7,12
      rawdata(i)=(6*qdata(InterpEq1(1,i-6))+3*qdata(InterpEq1(2,i-6))&
        -qdata(InterpEq1(3,i-6)))/8
    end do

    do i=13,15
      rawdata(i)=(4*qdata(InterpEq2(1,i-12))+4*qdata(InterpEq2(2,i-12))&
                 +2*qdata(InterpEq2(3,i-12))-qdata(InterpEq2(4,i-12))&
                 -qdata(InterpEq2(5,i-12)))/8
    end do

    do i=1,4
      subqdata_qtrig(:,i)=rawdata(subqtrig2qtrig(:,i))
    end do
  end function QTrigInterpolation

  function CollectQWeights(qweights)
    implicit none
    real(wp), dimension(15) :: subqweight
    real(wp), dimension(6,4), intent(in) :: qweights
    real(wp), dimension(6) :: CollectQWeights
    integer :: i

    subqweight=0
    do i=1,4
      subqweight(subqtrig2qtrig(:,i))=qweights(:,i)+subqweight(subqtrig2qtrig(:,i))
    end do

    CollectQWeights=matmul(subqweight,subqW2qW)
  end function CollectQWeights
    
! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k)).
  recursive function QuadTrigTheta(Eqtrig,eF,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: Eqtrig
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(6) :: qw
    real(wp), dimension(3) :: res
    real(wp), dimension(6,4) :: eqtrigs, qweights
    integer :: i,j

    if (iter==0) then
      if (maxval(Eqtrig)<=eF) then
        qw=(/0.041666666666666664_wp, 0.041666666666666664_wp,&
             0.041666666666666664_wp, 0.125_wp, 0.125_wp, 0.125_wp/)
      else
        qw=0
        if (minval(Eqtrig)<eF) then
          do i=1,4
            res = LinTrigTheta(Eqtrig(qtrig2trig(:,i)),eF)/4
            do j=1,3
              qw(qtrig2trig(j,i))=res(j)+qw(qtrig2trig(j,i))
            end do
          end do
        end if
      end if
    else
      eqtrigs=QTrigInterpolation(Eqtrig)
      do i=1,4
        qweights(:,i)=QuadTrigTheta(eqtrigs(:,i),eF,iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadtrigTheta

! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = δ(eF-E(k)).
  recursive function QuadTrigDelta(Eqtrig,eF,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: Eqtrig
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(3) :: res
    real(wp), dimension(6) :: qw
    real(wp), dimension(6,4) :: eqtrigs, qweights
    integer :: i,j

    if (iter==0) then
      if (minval(Eqtrig)<eF .and. eF<maxval(Eqtrig)) then
        qw=0
        do i=1,4
          res=LinTrigDelta(Eqtrig(qtrig2trig(:,i)),eF)/4
          do j=1,3
            qw(qtrig2trig(j,i))=res(j)+qw(qtrig2trig(j,i))
          end do
        end do
      else
        qw=0
      end if
    else
      eqtrigs=QTrigInterpolation(Eqtrig)
      qweights=0
      do i=1,4
        qweights(:,i)=QuadTrigDelta(eqtrigs(:,i),eF,iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTrigDelta

! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k))δ(D(k)).
  recursive function QuadTrigThetaDelta(Eqtrig,eF,Dqtrig,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: Eqtrig, Dqtrig
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(3) :: res
    real(wp), dimension(6) :: qw
    real(wp), dimension(6,4) :: eqtrigs,dqtrigs,qweights
    integer :: i,j

    if (iter==0) then
      if (minval(Eqtrig)<eF) then
        qw=0
        do i=1,4
          res=LinTrigThetaDelta(Eqtrig(qtrig2trig(:,i)),eF,Dqtrig(qtrig2trig(:,i)))/4
          do j=1,3
            qw(qtrig2trig(j,i))=res(j)+qw(qtrig2trig(j,i))
          end do
        end do
      else
        qw=0
      end if
    else
      eqtrigs=QTrigInterpolation(Eqtrig)
      dqtrigs=QTrigInterpolation(Dqtrig)
      qweights=0
      do i=1,4
        qweights(:,i)=QuadTrigThetaDelta(eqtrigs(:,i),eF,dqtrigs(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTrigThetaDelta

! Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k))⋅1/D(k)
  recursive function QuadTrigThetaFrakD(Eqtrig,eF,Dqtrig,iter) result(qw)
    implicit none
    real(wp), dimension(6), intent(in) :: Eqtrig, Dqtrig
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp) :: dmax, dmin, rtol
    real(wp), dimension(6) :: qw
    real(wp), dimension(3) :: res
    real(wp), dimension(6,4) :: eqtrigs,dqtrigs,qweights
    integer :: i,j

    dmax=maxval(Dqtrig)
    dmin=minval(Dqtrig)
    rtol=0.25_wp
    if (dmax*dmin>0 .and. isapprox(dmax,dmin,rtol)) then
      qw=QuadTrigTheta(Eqtrig,eF,iter)/Dqtrig
    else
      if (iter==0) then
        qw=0
        if (minval(Eqtrig)<eF) then
          do i=1,4
            res=LinTrigThetaFrakD(Eqtrig(qtrig2trig(:,i)),eF,Dqtrig(qtrig2trig(:,i)))/4
            do j=1,3
              qw(qtrig2trig(j,i))=res(j)+qw(qtrig2trig(j,i))
            end do
          end do
        end if
      else
        eqtrigs=QTrigInterpolation(Eqtrig)
        dqtrigs=QTrigInterpolation(Dqtrig)
        do i=1,4
          qweights(:,i)=QuadTrigThetaFrakD(eqtrigs(:,i),eF,dqtrigs(:,i),iter-1)
        end do
        qw=CollectQWeights(qweights)
      end if
    end if
  end function QuadTrigThetaFrakD
end module QuadTrig
