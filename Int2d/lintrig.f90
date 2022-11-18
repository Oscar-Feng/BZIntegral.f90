module LinTrig
  implicit none

  INTRINSIC KIND, SELECTED_REAL_KIND

!     INTEGER, PARAMETER :: skind = SELECTED_REAL_KIND(p=6, r=37)
!     INTEGER, PARAMETER :: dkind = SELECTED_REAL_KIND(p=15, r=307)
  INTEGER, PARAMETER :: skind = KIND(0.0E0)
  INTEGER, PARAMETER :: dkind = KIND(0.0D0)

  INTEGER, PARAMETER :: wp = dkind
contains
  subroutine subtrig1(E,eF,lininterp1,vol1)
    implicit none
    real(wp), intent(in), dimension(3) :: E
    real(wp), intent(in) :: eF
    real(wp), intent(out), dimension(3,3) :: lininterp1
    real(wp), intent(out) :: vol1

    lininterp1=0

    lininterp1(1,1) = 1
    lininterp1(2,1) = (E(2)-eF)/(E(2)-E(1))
    lininterp1(2,2) = (eF-E(1))/(E(2)-E(1))
    lininterp1(3,1) = (E(3)-eF)/(E(3)-E(1))
    lininterp1(3,3) = (eF-E(1))/(E(3)-E(1))

    vol1 = (eF-E(1))**2/(E(2)-E(1))/(E(3)-E(1))
  end subroutine subtrig1

  subroutine subtrig2(E,eF,lininterp2,vol2)
    implicit none
    real(wp), intent(in),dimension(3)::E
    real(wp), intent(in)::eF
    real(wp), intent(out), dimension(3,3) :: lininterp2
    real(wp), intent(out)::vol2

    lininterp2=0

    lininterp2(1,1) = (E(3)-eF)/(E(3)-E(1))
    lininterp2(1,3) = (eF-E(1))/(E(3)-E(1))
    lininterp2(2,2) = (E(3)-eF)/(E(3)-E(2))
    lininterp2(2,3) = (eF-E(2))/(E(3)-E(2))
    lininterp2(3,3) = 1

    vol2 = (E(3)-eF)**2/(E(3)-E(1))/(E(3)-E(2))
  end subroutine subtrig2
  
  function sortperm(array)
    implicit none
    real(wp), dimension(:), intent(in) :: array
    integer :: leng, count, index, upper_bound, temp
    integer, dimension(lbound(array, dim=1):ubound(array, dim=1)) :: sortperm

    leng=size(array)
    upper_bound=leng-1

    do count=1,leng
      sortperm(count)=count
    end do

    do count=1, leng
      do index=1, upper_bound
        if (array(sortperm(index))>array(sortperm(index+1))) then
          temp=sortperm(index)
          sortperm(index)=sortperm(index+1)
          sortperm(index+1)=temp
        end if
      end do
      upper_bound=upper_bound-1
    end do
  end function sortperm

! Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k)).
! With Blöchl correction in Blöchl et al. PhysRevB.49.16223 (1994)
  function LinTrigTheta_Blochl(Eraw,eF) result(res)
    implicit none
    real(wp), dimension(3), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(3) :: ind
    real(wp), dimension(3) :: res, w, E
    real(wp) :: dos, C

    ind=sortperm(Eraw)
    E=Eraw(ind)

    if (E(1)>=eF) then
      w=0
      dos=0
    else if (E(3)<eF) then
      w=1.0_wp/6.0_wp
      dos=0
    else if (E(1)<eF .and. eF<=E(2)) then
      C=(eF-E(1))**2/((E(2)-E(1))*(E(3)-E(1)))/6
      w=(/3-(eF-E(1))*(1/(E(2)-E(1))+1/(E(3)-E(1))), (eF-E(1))/(E(2)-E(1)), &
         (eF-E(1))/(E(3)-E(1))/)*C
      dos=2*(eF-E(1))/((E(2)-E(1))*(E(3)-E(1)))/2
    else if (E(2)<eF .and. eF<=E(3)) then
      C=(E(3)-eF)**2/((E(3)-E(1))*(E(3)-E(2)))/6
      w=1.0_wp/6.0_wp-(/(E(3)-eF)/(E(3)-E(1)),(E(3)-eF)/(E(3)-E(2)),&
        (3-(E(3)-eF)*(1/(E(3)-E(1))+1/(E(3)-E(2))))/)*C
      dos=2*(E(3)-eF)/((E(3)-E(1))*(E(3)-E(2)))/2
    end if

    res(ind)=w+dos*(sum(E)-3*E)/24
  end function LinTrigTheta_Blochl

! Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k)).
! By splitting triangle by the Fermi surface
  function LinTrigTheta(Eraw,eF) result(res)
    implicit none
    real(wp), dimension(3), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(3) :: ind
    real(wp) :: vol1, vol2
    real(wp), dimension(3) :: res, w, E
    real(wp), dimension(3,3) :: lininterp1,lininterp2

    ind=sortperm(Eraw)
    E=Eraw(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(3)<eF) then
      w=1.0_wp/6.0_wp
    else if (E(1)<eF .and. E(2)>=eF) then
      call subtrig1(E,eF,lininterp1,vol1)
      w=matmul(transpose(lininterp1),(/1.0_wp/6.0_wp,1.0_wp/6.0_wp,1.0_wp/6.0_wp/))*vol1
    else if (E(2)<eF .and. E(3)>=eF) then
      call subtrig2(E,eF,lininterp2,vol2)
      w=(/1.0_wp/6.0_wp,1.0_wp/6.0_wp,1.0_wp/6.0_wp/)-matmul(transpose(&
         lininterp2),(/1.0_wp/6.0_wp,1.0_wp/6.0_wp,1.0_wp/6.0_wp/))*vol2
    end if

    res(ind)=w
  end function LinTrigTheta

! Linear triangle rule in a single triangle, weight function is W(k) = δ(eF-E(k)).
! By differentiating The closed form rule for Θ(eF-E(k))
  function LinTrigDelta(Eraw,eF) result(res)
    implicit none
    real(wp), dimension(3), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(3) :: ind
    real(wp), dimension(3) :: E,res,w
    real(wp) :: C,dC

    ind=sortperm(Eraw)
    E=Eraw(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(3)<eF) then
      w=0
    else if (E(1)<eF .and. E(2)>=eF) then
      C=(eF-E(1))**2/((E(2)-E(1))*(E(3)-E(1)))/6
      dC=2*(eF-E(1))/((E(2)-E(1))*(E(3)-E(1)))/6
      w=(/(3-(eF-E(1))*(1/(E(2)-E(1))+1/(E(3)-E(1)))),(eF-E(1))/(E(2)-E(1)),&
         (eF-E(1))/(E(3)-E(1))/)*dC+(/-(1/(E(2)-E(1))+1/(E(3)-E(1))),&
         1/(E(2)-E(1)),1/(E(3)-E(1))/)*C
    else if (E(2)<eF .and. E(3)>=eF) then
      C=(E(3)-eF)**2/((E(3)-E(1))*(E(3)-E(2)))/6
      dC=-2*(E(3)-eF)/((E(3)-E(1))*(E(3)-E(2)))/6
      w=-dC*(/(E(3)-eF)/(E(3)-E(1)),(E(3)-eF)/(E(3)-E(2)),&
        (3-(E(3)-eF)*(1/(E(3)-E(1))+1/(E(3)-E(2))))/)+(/1/(E(3)-E(1)),&
        1/(E(3)-E(2)),-(1/(E(3)-E(1))+1/(E(3)-E(2)))/)*C
    end if

    res(ind)=w
  end function LinTrigDelta

! Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k))⋅δ(D(k)).
! By splitting triangle by the Fermi surface
  function LinTrigThetaDelta(Eraw,eF,Draw) result(res)
    implicit none
    real(wp), dimension(3), intent(in) :: Eraw,Draw
    real(wp), intent(in) :: eF
    integer, dimension(3) :: ind
    real(wp) :: dF,vol1, vol2
    real(wp), dimension(3) :: E,D,w,res
    real(wp), dimension(3,3) :: lininterp1,lininterp2

    ind=sortperm(Eraw)
    E=Eraw(ind)
    D=Draw(ind)
    dF=0.0_wp

    if (E(1)>=eF) then
      w=0
    else if (E(3)<eF) then
      w=LinTrigDelta(D,dF)
    else if (E(1)<eF .and. E(2)>=eF) then
      call subtrig1(E,eF,lininterp1,vol1)
      w=matmul(transpose(lininterp1),LinTrigDelta(matmul(lininterp1,D),dF))*vol1
    else if (E(2)<eF .and. E(3)>=eF) then
      call subtrig2(E,eF,lininterp2,vol2)
      w=LinTrigDelta(D,dF)-matmul(transpose(lininterp2),LinTrigDelta(matmul(&
        lininterp2,D),dF))*vol2
    end if

    res(ind)=w
  end function LinTrigThetaDelta

  function isapprox(x,y,rtol)
    implicit none
    real(wp), intent(in) :: x,y,rtol
    logical :: isapprox

    if (abs(x-y)<=rtol*max(abs(x),abs(y))) then
      isapprox=.true.
    else
      isapprox=.false.
    end if
  end function isapprox

! weight for x: RawExp(x,y,z)+RawExp(x,z,y)
! singular when x->y , y->z
  function RawExp(x,y,z)
    implicit none
    real(wp), intent(in) :: x,y,z
    real(kind=16) :: a,b,c,RawExp !It is said that such kind declaration is not good...
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(x,y,rtol) .or. isapprox(y,z,rtol)) then
      a=x
      b=y
      c=z
      RawExp = -(b/(2*(a - b)*(-b + c))) - (b**2*(-log(abs(a)) + log(abs(b))))/&
               (2*(a - b)**2*(-b + c))
    else
      RawExp=real(-(y/(2*(x - y)*(-y + z))) - (y**2*(-log(abs(x)) + log(abs(y))))/&
             (2*(x - y)**2*(-y + z)),16)
    end if
  end function RawExp

! limit case: RawExp(x,y,z), x->y
  function ExpXeqY(x,y,z)
    implicit none
    real(wp), intent(in) :: x,y,z
    real(kind=16) :: b,c,ExpXeqY
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(y,z,rtol)) then
      b=y
      c=z
      ExpXeqY = 1/(4*b - 4*c)
    else
      ExpXeqY =real(1/(4*y - 4*z),16)
    end if
  end function ExpXeqY

! limit case: RawExp(x,y,z)+RawExp(x,z,y), y->z
  function ExpYeqZ(x,y,z)
    implicit none
    real(wp), intent(in) :: x,y,z
    real(kind=16) :: a,c,ExpYeqZ
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(x,y,rtol)) then
      a=x
      c=z
      ExpYeqZ=(a**2 - c**2 - 2*a*c*(log(abs(a))-log(abs(c))))/(2*(a - c)**3)
    else
      ExpYeqZ=real((x**2-z**2-2*x*z*(log(abs(x))-log(abs(z))))/(2*(x-z)**3),16)
    end if
  end function ExpYeqZ

! weight for x: RawExp(x,y,z)+RawExp(x,z,w)
  function FracTrigWeight(x,y,z) result(res)
    real(wp), intent(in) :: x,y,z
    real(wp), dimension(2) :: yz
    real(wp), dimension(3) :: v
    integer, dimension(2) :: ind
    real(wp) :: rtol,res

    rtol=0.00000001_wp
    yz=(/y,z/)
    ind=sortperm(yz)
    v=(/x,yz(ind(1)),yz(ind(2))/)

    if (isapprox(v(1),v(2),rtol)) then
      if (isapprox(v(2),v(3),rtol)) then
        res=1/v(1)/6
      else
        res = ExpXeqY(v(1),v(2),v(3))+RawExp(v(1),v(3),v(2))
      end if
    else if (isapprox(v(1),v(3),rtol)) then
        res = ExpXeqY(v(1),v(3),v(2))+RawExp(v(1),v(2),v(3))
    else if (isapprox(v(2),v(3),rtol)) then
        res = ExpYeqZ(v(1),v(2),v(3))
    else
        res = RawExp(v(1),v(2),v(3))+RawExp(v(1),v(3),v(2))
    end if
  end function FracTrigWeight

! Linear triangle rule in a single triangle, weight function is W(k) = 1/𝔇(k)
  function LinTrigFrakD(frakD)
    implicit none
    real(wp), dimension(3), intent(in) :: frakD
    real(wp), dimension(3) :: LinTrigFrakD

    LinTrigFrakD=(/FracTrigWeight(frakD(1),frakD(2),frakD(3)),&
                 FracTrigWeight(frakD(2),frakD(3),frakD(1)),&
                 FracTrigWeight(frakD(3),frakD(1),frakD(2))/)
  end function LinTrigFrakD

! Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k))⋅ 1/D(k)
  function LinTrigThetaFrakD(Eraw,eF,Dnom) result(res)
    implicit none
    real(wp), intent(in), dimension(3) :: Eraw, Dnom
    real(wp), intent(in) :: eF
    integer, dimension(3) :: ind
    real(wp), dimension(3) :: E,D,w,res
    real(wp) :: vol1,vol2
    real(wp), dimension(3,3) :: lininterp1,lininterp2

    ind=sortperm(Eraw)
    E=Eraw(ind)
    D=Dnom(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(3)<eF) then
      w=LinTrigFrakD(D)
    else if (E(1)<eF .and. E(2)>=eF) then
      call subtrig1(E,eF,lininterp1,vol1)
      w=matmul(transpose(lininterp1),LinTrigFrakD(matmul(lininterp1,D)))*vol1
    else if (E(2)<eF .and. E(3)>=eF) then
      call subtrig2(E,eF,lininterp2,vol2)
      w=LinTrigFrakD(D)-matmul(transpose(lininterp2),LinTrigFrakD(matmul(lininterp2,D)))*vol2
    end if
    res(ind)=w
  end function LinTrigThetaFrakD
end module LinTrig
