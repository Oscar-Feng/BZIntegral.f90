module QuadTetra
  use LinTetra
  implicit none

  integer, dimension(3,12), parameter :: InterpEq1=reshape((/7,1,4, 7,4,1,&
                            9,4,3, 9,3,4, 8,3,2, 8,2,3, 5,2,1,&
                            5,1,2, 6,1,3, 6,3,1, 10,2,4, 10,4,2/),(/3,12/))
  integer, dimension(5,12), parameter :: InterpEq2=reshape((/5,7,10,2,4,&
                      5,10,7,1,4, 7,10,5,1,2, 10,9,8,3,2, 10,8,9,4,3,&
                      9,8,10,2,4, 5,8,6,1,3, 5,6,8,2,3, 6,8,5,1,2,&
                      6,7,9,3,4, 6,9,7,1,4, 7,9,6,1,3/),(/5,12/))
  integer, dimension(10), parameter :: InterpEq3=(/5,6,7,8,9,10,1,2,3,4/)
  integer, dimension(10,8), parameter :: subqtetra2qtetra=reshape(&
                (/1,5,6,7, 18,19,11, 30,32,23, 5,2,8,10, 17,29,24, 16,27,21,&
                  6,8,3,9, 31,20,33, 15,14,28, 7,10,9,4, 25,34,12, 26,13,22,&
                  8,5,6,7, 29,31,35, 30,32,23, 5,7,8,10, 23,29,24, 35,27,25,&
                  9,8,7,6, 28,34,33, 35,32,31, 7,10,9,8, 25,34,35, 26,28,27/),&
                  (/10,8/))
  integer, dimension(4,8), parameter :: qtetra2tetra=reshape((/1,5,6,7,&
                                      5,2,8,10, 6,8,3,9, 7,10,9,4, 8,5,6,7,&
                                      5,7,8,10, 9,8,7,6, 7,10,9,8/),(/4,8/))
  real(wp), dimension(35,10), parameter :: subqW2qW=reshape((/0.125_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.046875_wp,&
    -0.015625_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,0.046875_wp,&
    0.046875_wp,-0.015625_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,-0.015625_wp,&
    0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,0.0_wp,-0.015625_wp,0.0_wp,-0.015625_wp,&
    -0.015625_wp,-0.015625_wp, 0.0_wp,0.125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,&
    0.046875_wp,0.046875_wp,-0.015625_wp,0.0_wp,0.0_wp,0.046875_wp,&
    -0.015625_wp,-0.015625_wp,0.0_wp,-0.015625_wp,-0.015625_wp,0.0_wp,&
    -0.015625_wp,0.0_wp,-0.015625_wp,-0.015625_wp,0.0_wp,0.0_wp,0.0_wp,&
    -0.015625_wp,0.0_wp,0.0_wp,0.125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,0.046875_wp,0.046875_wp,&
    -0.015625_wp,0.0_wp,0.0_wp,-0.015625_wp,0.046875_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,-0.015625_wp,-0.015625_wp,0.0_wp,-0.015625_wp,-0.015625_wp,&
    0.0_wp,-0.015625_wp,0.0_wp,-0.015625_wp,-0.015625_wp, 0.0_wp,0.0_wp,0.0_wp,&
    0.125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,&
    0.046875_wp,0.046875_wp,-0.015625_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,-0.015625_wp,0.046875_wp,-0.015625_wp,-0.015625_wp,0.0_wp,0.0_wp,&
    -0.015625_wp,-0.015625_wp,0.0_wp,0.0_wp,0.0_wp,-0.015625_wp,-0.015625_wp,&
    0.0_wp,-0.015625_wp, 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.125_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.09375_wp,&
    0.09375_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0625_wp,0.0625_wp,0.03125_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0625_wp,0.0625_wp,0.03125_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.03125_wp, 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.125_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.09375_wp,0.09375_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.03125_wp,0.0625_wp,0.0625_wp,0.0625_wp,0.0625_wp,0.03125_wp,&
    0.03125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.125_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.09375_wp,0.09375_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0625_wp,0.03125_wp,0.0625_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0625_wp,0.03125_wp,0.0625_wp,&
    0.03125_wp, 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.125_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.09375_wp,0.09375_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.03125_wp,&
    0.0625_wp,0.0625_wp,0.0625_wp,0.03125_wp,0.0625_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.03125_wp, 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.125_wp,0.0_wp,0.0_wp,0.0_wp,0.09375_wp,0.09375_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0625_wp,&
    0.03125_wp,0.0625_wp,0.0_wp,0.0_wp,0.0_wp,0.03125_wp,0.0625_wp,0.0625_wp,&
    0.03125_wp, 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.0_wp,0.09375_wp,0.09375_wp,0.03125_wp,0.0625_wp,0.0625_wp,0.0625_wp,&
    0.0625_wp,0.03125_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,&
    0.03125_wp/),(/35,10/))
contains
! Given data on vertex of a quadratic tetrahedron, compute quadratic interpolated value at vertex  !   with it recognized as a subqtetrahedron.
! results are given in terms of list of data on each of its 8 subordinate quadratic tetrahedrons
  function QTetraInterpolation(qdata) result(subqdata_qtetra)
    implicit none
    real(wp), dimension(10), intent(in) :: qdata
    real(wp), dimension(10,8) :: subqdata_qtetra
    real(wp),dimension(35):: rawdata
    integer::i

    rawdata=0
    do i=1,10
      rawdata(i)=qdata(i)
    end do

    do i=11,22
      rawdata(i)=(6*qdata(InterpEq1(1,i-10))+3*qdata(InterpEq1(2,i-10))&
      -qdata(InterpEq1(3,i-10)))/8
    end do

    do i=23,34
      rawdata(i)=(4*qdata(InterpEq2(1,i-22))+4*qdata(InterpEq2(2,i-22))&
                   +2*qdata(InterpEq2(3,i-22))-qdata(InterpEq2(4,i-22))&
                   -qdata(InterpEq2(5,i-22)))/8
    end do

    rawdata(35)=(2*sum(qdata(InterpEq3(1:6)))-sum(qdata(InterpEq3(7:10))))/8

    do i=1,8
        subqdata_qtetra(:,i)=rawdata(subqtetra2qtetra(:,i))
    end do
  end function QTetraInterpolation

! Given weights of 8 subordinate quadratic tetrahedrons (data interpolated from master quadratic tetrahedrons),
!   compute weights for master quadratic tetrahedrons
! sum(qweights[i,:]) should have same normalization convention as sum(result)
  function CollectQWeights(qweights)
    implicit none
    real(wp), dimension(10) :: CollectQWeights
    real(wp), dimension(35) :: subqweight
    real(wp), dimension(10,8), intent(in):: qweights
    integer :: i

    subqweight=0
    do i=1,8
      subqweight(subqtetra2qtetra(:,i))=qweights(:,i)+subqweight(subqtetra2qtetra(:,i))
    end do

    CollectQWeights=matmul(subqweight,subqW2qW)
  end function CollectQWeights

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k)).
  recursive function QuadTetraTheta(Eqtetra, eF, iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: Eqtetra
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: qweights, eqtetras
    real(wp), dimension(10) :: qw
    real(wp), dimension(4) :: res
    integer :: i,j

    if (iter==0) then
      if (maxval(Eqtetra)<=eF) then
        qw=(/0.005208333333333333_wp,0.005208333333333333_wp,&
             0.005208333333333333_wp,0.005208333333333333_wp,&
             0.020833333333333332_wp,0.020833333333333332_wp,&
             0.031249999999999997_wp,0.031249999999999997_wp,&
             0.020833333333333332_wp,0.020833333333333332_wp/)
      else
        qw=0
        if (minval(Eqtetra)<eF) then
          do i=1,8
            res=LinTetraTheta(Eqtetra(qtetra2tetra(:,i)), eF)/8
            do j=1,4
              qw(qtetra2tetra(j,i))=res(j)+qw(qtetra2tetra(j,i))
            end do
          end do
        end if
      end if
    else
      eqtetras=QTetraInterpolation(Eqtetra)
      do i=1,8
        qweights(:,i)=QuadTetraTheta(eqtetras(:,i),eF,iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraTheta

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = δ(eF-E(k)).
  recursive function QuadTetraDelta(Eqtetra, eF, iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: Eqtetra
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: qweights, eqtetras
    real(wp), dimension(10) :: qw
    real(wp), dimension(4) :: res
    integer :: i,j

    if (iter==0) then
      qw=0
      if (minval(Eqtetra)<eF .and. maxval(Eqtetra)>eF) then
        do i=1,8
          res=LinTetraDelta(Eqtetra(qtetra2tetra(:,i)), eF)/8
          do j=1,4
            qw(qtetra2tetra(j,i))=res(j)+qw(qtetra2tetra(j,i))
          end do
        end do
      end if
    else
      eqtetras=QTetraInterpolation(Eqtetra)
      do i=1,8
        qweights(:,i)=QuadTetraDelta(eqtetras(:,i),eF,iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraDelta

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k))δ(D(k)).
  recursive function QuadTetraThetaDelta(Eqtetra,eF,Dqtetra,iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: Eqtetra, Dqtetra
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: qweights, eqtetras, dqtetras
    real(wp), dimension(10) :: qw
    real(wp), dimension(4) :: res
    integer :: i,j

    if (iter==0) then
      qw=0
      if (minval(Eqtetra)<eF) then
        do i=1,8
          res=LinTetraThetaDelta(Eqtetra(qtetra2tetra(:,i)),eF,Dqtetra(&
          qtetra2tetra(:,i)))/8
          do j=1,4
            qw(qtetra2tetra(j,i))=res(j)+qw(qtetra2tetra(j,i))
          end do
        end do
      end if
    else
      eqtetras=QTetraInterpolation(Eqtetra)
      dqtetras=QTetraInterpolation(Dqtetra)
      do i=1,8
        qweights(:,i)=QuadTetraThetaDelta(eqtetras(:,i),eF,dqtetras(:,i),iter-1)
      end do
      qw=CollectQWeights(qweights)
    end if
  end function QuadTetraThetaDelta

! Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅1/D(k)
  recursive function QuadTetraThetaFrakD(Eqtetra,eF,Dqtetra,iter) result(qw)
    implicit none
    real(wp), dimension(10), intent(in) :: Eqtetra, Dqtetra
    real(wp), intent(in) :: eF
    integer, intent(in) :: iter
    real(wp), dimension(10,8) :: qweights, eqtetras, dqtetras
    real(wp), dimension(10) :: qw
    real(wp), dimension(4) :: res
    real(wp) ::  dmax, dmin, rtol
    integer :: i,j

    dmax=maxval(Dqtetra)
    dmin=minval(Dqtetra)
    rtol=0.25_wp
    if (dmax*dmin>0 .and. isapprox(dmax, dmin, rtol)) then
      qw=QuadTetraTheta(Eqtetra,eF,iter)/Dqtetra
    else
      if (iter==0) then
        qw=0
        if (minval(Eqtetra)<eF) then
          do i=1,8
            res=LinTetraThetaFrakD(Eqtetra(qtetra2tetra(:,i)),eF,Dqtetra(qtetra2tetra(:,i)))/8
            do j=1,4
              qw(qtetra2tetra(j,i))=res(j)+qw(qtetra2tetra(j,i))
            end do
          end do
        end if
      else
        eqtetras=QTetraInterpolation(Eqtetra)
        dqtetras=QTetraInterpolation(Dqtetra)
        do i=1,8
          qweights(:,i)=QuadTetraThetaFrakD(eqtetras(:,i),eF,dqtetras(:,i),iter-1)
        end do
        qw=CollectQWeights(qweights)
      end if
    end if
  end function QuadTetraThetaFrakD
end module QuadTetra
