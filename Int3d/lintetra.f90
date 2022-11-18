module LinTetra
  implicit none

  INTRINSIC KIND

!     INTEGER, PARAMETER :: skind = SELECTED_REAL_KIND(p=6, r=37)
!     INTEGER, PARAMETER :: dkind = SELECTED_REAL_KIND(p=15, r=307)
  INTEGER, PARAMETER :: skind = KIND(0.0E0)
  INTEGER, PARAMETER :: dkind = KIND(0.0D0)

  INTEGER, PARAMETER :: wp = dkind
contains
  subroutine subtetra1(E,eF,lininterp1,vol1)
    implicit none
    real(wp), intent(in), dimension(4)::E
    real(wp), intent(in) ::eF
    real(wp), intent(out), dimension(4,4) :: lininterp1
    real(wp), intent(out)::vol1

    lininterp1=0

    lininterp1(1,1)=1
    lininterp1(2,1)=(E(2)-eF)/(E(2)-E(1))
    lininterp1(2,2)=(eF-E(1))/(E(2)-E(1))
    lininterp1(3,1)=(E(3)-eF)/(E(3)-E(1))
    lininterp1(3,3)=(eF-E(1))/(E(3)-E(1))
    lininterp1(4,1)=(E(4)-eF)/(E(4)-E(1))
    lininterp1(4,4)=(eF-E(1))/(E(4)-E(1))

    vol1=(eF-E(1))**3/(E(2)-E(1))/(E(3)-E(1))/(E(4)-E(1))
  end subroutine subtetra1

  subroutine subtetra2(E, eF, lininterp21, lininterp22, lininterp23, vol2)
    implicit none
    real(wp), intent(in),dimension(4)::E
    real(wp), intent(in)::eF
    real(wp), intent(out), dimension(4,4) :: lininterp21, lininterp22,&
                                               lininterp23
    real(wp), intent(out), dimension(3)::vol2

    lininterp21=0

    lininterp21(1,1)=1
    lininterp21(2,2)=(E(3)-eF)/(E(3)-E(2))
    lininterp21(2,3)=(eF-E(2))/(E(3)-E(2))
    lininterp21(3,1)=(E(3)-eF)/(E(3)-E(1))
    lininterp21(3,3)=(eF-E(1))/(E(3)-E(1))
    lininterp21(4,1)=(E(4)-eF)/(E(4)-E(1))
    lininterp21(4,4)=(eF-E(1))/(E(4)-E(1))

    vol2(1)=(eF-E(1))**2/(E(3)-E(1))/(E(4)-E(1))*(eF-E(3))/(E(2)-E(3))

    lininterp22=0

    lininterp22(1,1)=1
    lininterp22(2,2)=(E(4)-eF)/(E(4)-E(2))
    lininterp22(2,4)=(eF-E(2))/(E(4)-E(2))
    lininterp22(3,2)=(E(3)-eF)/(E(3)-E(2))
    lininterp22(3,3)=(eF-E(2))/(E(3)-E(2))
    lininterp22(4,1)=(E(4)-eF)/(E(4)-E(1))
    lininterp22(4,4)=(eF-E(1))/(E(4)-E(1))

    vol2(2)=(eF-E(2))/(E(3)-E(2))*(eF-E(4))/(E(2)-E(4))*(eF-E(1))/(E(4)-E(1))

    lininterp23=0

    lininterp23(1,1)=1
    lininterp23(2,2)=1
    lininterp23(3,2)=(E(3)-eF)/(E(3)-E(2))
    lininterp23(3,3)=(eF-E(2))/(E(3)-E(2))
    lininterp23(4,2)=(E(4)-eF)/(E(4)-E(2))
    lininterp23(4,4)=(eF-E(2))/(E(4)-E(2))

    vol2(3)=(eF-E(2))**2/(E(4)-E(2))/(E(3)-E(2))
  end subroutine subtetra2

  subroutine subtetra3(E, eF, lininterp3, vol3)
    implicit none
    real(wp), intent(in),dimension(4)::E
    real(wp), intent(in)::eF
    real(wp), intent(out), dimension(4,4) :: lininterp3
    real(wp), intent(out)::vol3

    lininterp3=0

    lininterp3(1,1)=(E(4)-eF)/(E(4)-E(1))
    lininterp3(1,4)=(eF-E(1))/(E(4)-E(1))
    lininterp3(2,2)=(E(4)-eF)/(E(4)-E(2))
    lininterp3(2,4)=(eF-E(2))/(E(4)-E(2))
    lininterp3(3,3)=(E(4)-eF)/(E(4)-E(3))
    lininterp3(3,4)=(eF-E(3))/(E(4)-E(3))
    lininterp3(4,4)=1

    vol3=(E(4)-eF)**3/(E(4)-E(1))/(E(4)-E(3))/(E(4)-E(2))
  end subroutine subtetra3
    
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

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k)).
! With Blöchl correction in Blöchl et al. PhysRevB.49.16223 (1994)
  function LinTetraTheta_Blochl(Eraw, eF) result(res)
    implicit none
    real(wp), dimension(4), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(4) :: ind
    real(wp), dimension(4) :: res, w, E
    real(wp) :: dos, C, C1, C2, C3

    ind=sortperm(Eraw)
    E=Eraw(ind)
    if (E(1)>=eF) then
      w=0
      dos=0
    else if (E(4)<eF) then
      w=1.0_wp/24.0_wp
      dos=0
    else if (E(1)<eF .and. eF<=E(2)) then
      C=(eF-E(1))**3/(E(2)-E(1))/(E(3)-E(1))/(E(4)-E(1))/24
      w=(/(4-(eF-E(1))*(1/(E(2)-E(1))+1/(E(3)-E(1))+1/(E(4)-E(1)))),&
        (eF-E(1))/(E(2)-E(1)), (eF-E(1))/(E(3)-E(1)), (eF-E(1))/(E(4)-E(1))/)*C
      dos=3*(eF-E(1))**2/((E(2)-E(1))*(E(3)-E(1))*(E(4)-E(1)))/6
    else if (E(2)<eF .and. eF<=E(3)) then
      C1=(eF-E(1))**2/((E(4)-E(1))*(E(3)-E(1)))/24
      C2=(eF-E(1))*(eF-E(2))*(E(3)-eF)/((E(4)-E(1))*(E(3)-E(2))*(E(3)-E(1)))/24
      C3=(eF-E(2))**2*(E(4)-eF)/((E(4)-E(2))*(E(3)-E(2))*(E(4)-E(1)))/24
      w=(/C1+(C1+C2)*(E(3)-eF)/(E(3)-E(1))+(C1+C2+C3)*(E(4)-eF)/(E(4)-E(1)),&
        C1+C2+C3+(C2+C3)*(E(3)-eF)/(E(3)-E(2))+C3*(E(4)-eF)/(E(4)-E(2)),&
        (C1+C2)*(eF-E(1))/(E(3)-E(1))+(C2+C3)*(eF-E(2))/(E(3)-E(2)),&
        (C1+C2+C3)*(eF-E(1))/(E(4)-E(1))+C3*(eF-E(2))/(E(4)-E(2))/)
      dos=1/((E(3)-E(1))*(E(4)-E(1)))/6*(3*(E(2)-E(1))+6*(eF-E(2))-3*(E(3)-E(1)&
          +E(4)-E(2))*(eF-E(2))**2/((E(3)-E(2))*(E(4)-E(2))))
    else if (E(3)<eF .and. eF<=E(4)) then
      C=(E(4)-eF)**3/((E(4)-E(1))*(E(4)-E(2))*(E(4)-E(3)))/24
      w=1.0_wp/24.0_wp-(/(E(4)-eF)/(E(4)-E(1)), (E(4)-eF)/(E(4)-E(2)),&
        (E(4)-eF)/(E(4)-E(3)), (4-(E(4)-eF)*(1/(E(4)-E(1))+1/(E(4)-E(2))+1/(E(4)-E(3))))/)*C
      dos=3*(E(4)-eF)**2/((E(4)-E(1))*(E(4)-E(2))*(E(4)-E(3)))/6
    end if

    res(ind)=w+dos*(sum(E)-4*E)/40
  end function LinTetraTheta_Blochl

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k)).
! By splitting tetrahedron by the Fermi surface
  function LinTetraTheta(Eraw, eF) result(res)
    implicit none
    real(wp), dimension(4), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(4) :: ind
    real(wp) :: vol1, vol3
    real(wp), dimension(3)::vol2
    real(wp), dimension(4):: res, w, E
    real(wp), dimension(4,4) :: lininterp1, lininterp21, lininterp22,&
                                lininterp23, lininterp3

    ind=sortperm(Eraw)
    E=Eraw(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(4)<eF) then
      w=1.0_wp/24.0_wp
    else if (E(1)<eF .and. eF<=E(2)) then
      call subtetra1(E, eF, lininterp1, vol1)
      w=matmul(transpose(lininterp1),(/1.0_wp/24.0_wp,1.0_wp/24.0_wp,&
      1.0_wp/24.0_wp,1.0_wp/24.0_wp/))*vol1
    else if (E(2)<eF .and. eF<=E(3)) then
      call subtetra2(E, eF, lininterp21, lininterp22, lininterp23, vol2)
      w=matmul(transpose(lininterp21)*vol2(1)+transpose(lininterp22)*vol2(2)+&
        transpose(lininterp23)*vol2(3),(/1.0_wp/24.0_wp,1.0_wp/24.0_wp,&
        1.0_wp/24.0_wp,1.0_wp/24.0_wp/))
    else if (E(3)<eF .and. eF<=E(4)) then
      call subtetra3(E, eF, lininterp3, vol3)
      w=(/1.0_wp/24.0_wp,1.0_wp/24.0_wp,1.0_wp/24.0_wp,1.0_wp/24.0_wp/)-matmul(&
      transpose(lininterp3),(/1.0_wp/24.0_wp,1.0_wp/24.0_wp,1.0_wp/24.0_wp,&
      1.0_wp/24.0_wp/))*vol3
    end if

    res(ind)=w
  end function LinTetraTheta

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = δ(eF-E(k)).
! By differentiating The closed form rule for Θ(eF-E(k))
  function LinTetraDelta(Eraw,eF) result(res)
    implicit none
    real(wp), dimension(4), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    integer, dimension(4) :: ind
    real(wp), dimension(4) :: res, w, E
    real(wp) :: C, C1, C2, C3, dC, dC1, dC2, dC3

    ind=sortperm(Eraw)
    E=Eraw(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(4)<eF) then
      w=0
    else if (E(1)<eF .and. E(2)>=eF) then
      C=(eF-E(1))**3/(E(2)-E(1))/(E(3)-E(1))/(E(4)-E(1))/24
      dC=3*(eF-E(1))**2/(E(2)-E(1))/(E(3)-E(1))/(E(4)-E(1))/24
      w=(/4-(eF-E(1))*(1/(E(2)-E(1))+1/(E(3)-E(1))+1/(E(4)-E(1))),&
          (eF-E(1))/(E(2)-E(1)),(eF-E(1))/(E(3)-E(1)),(eF-E(1))/(E(4)-E(1))/)*dC+&
          (/-(1/(E(2)-E(1))+1/(E(3)-E(1))+1/(E(4)-E(1))),1/(E(2)-E(1)),&
          1/(E(3)-E(1)),1/(E(4)-E(1))/)*C
    else if (E(2)<eF .and. E(3)>=eF) then
      C1=(eF-E(1))**2/((E(4)-E(1))*(E(3)-E(1)))/24
      dC1=2*(eF-E(1))/((E(4)-E(1))*(E(3)-E(1)))/24
      C2=(eF-E(1))*(eF-E(2))*(E(3)-eF)/((E(4)-E(1))*(E(3)-E(2))*(E(3)-E(1)))/24
      dC2=((eF-E(2))*(E(3)-eF)+(eF-E(1))*(E(3)-eF)-(eF-E(1))*(eF-E(2)))/&
            ((E(4)-E(1))*(E(3)-E(2))*(E(3)-E(1)))/24
      C3=(eF-E(2))**2*(E(4)-eF)/((E(4)-E(2))*(E(3)-E(2))*(E(4)-E(1)))/24
      dC3=(2*(eF-E(2))*(E(4)-eF)-(eF-E(2))**2)/((E(4)-E(2))*(E(3)-E(2))*&
            (E(4)-E(1)))/24
      w=(/dC1+(dC1+dC2)*(E(3)-eF)/(E(3)-E(1))+(dC1+dC2+dC3)*(E(4)-eF)/(E(4)-E(1)),&
           dC1+dC2+dC3+(dC2+dC3)*(E(3)-eF)/(E(3)-E(2))+dC3*(E(4)-eF)/(E(4)-E(2)),&
           (dC1+dC2)*(eF-E(1))/(E(3)-E(1))+(dC2+dC3)*(eF-E(2))/(E(3)-E(2)),&
           (dC1+dC2+dC3)*(eF-E(1))/(E(4)-E(1))+dC3*(eF-E(2))/(E(4)-E(2))/)+&
           (/-(C1+C2)/(E(3)-E(1))-(C1+C2+C3)/(E(4)-E(1)),&
           -(C2+C3)/(E(3)-E(2))-C3/(E(4)-E(2)),(C1+C2)/(E(3)-E(1))+(C2+C3)/(E(3)-E(2)),&
           (C1+C2+C3)/(E(4)-E(1))+C3/(E(4)-E(2))/)
    else if (E(3)<eF .and. E(4)>=eF) then
      C=(E(4)-eF)**3/((E(4)-E(1))*(E(4)-E(2))*(E(4)-E(3)))/24
      dC=-3*(E(4)-eF)**2/((E(4)-E(1))*(E(4)-E(2))*(E(4)-E(3)))/24
      w=-dC*(/(E(4)-eF)/(E(4)-E(1)),(E(4)-eF)/(E(4)-E(2)),(E(4)-eF)/(E(4)-E(3)),&
          (4-(E(4)-eF)*(1/(E(4)-E(1))+1/(E(4)-E(2))+1/(E(4)-E(3))))/)+&
          C*(/1/(E(4)-E(1)),1/(E(4)-E(2)),1/(E(4)-E(3)),&
          -(1/(E(4)-E(1))+1/(E(4)-E(2))+1/(E(4)-E(3)))/)
    end if

    res(ind)=w
  end function LinTetraDelta

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅δ(D(k)).
! By splitting tetrahedron by the Fermi surface
  function LinTetraThetaDelta(Eraw,eF,Draw) result(res)
    implicit none
    real(wp), dimension(4), intent(in) :: Eraw
    real(wp), intent(in) :: eF
    real(wp), dimension(4), intent(in) :: Draw
    integer, dimension(4) :: ind
    real(wp), dimension(4):: res, w, E, D
    real(wp) :: vol1, vol3, dF
    real(wp), dimension(3)::vol2
    real(wp), dimension(4,4) :: lininterp1, lininterp21, lininterp22,&
                                lininterp23, lininterp3

    ind=sortperm(Eraw)
    E=Eraw(ind)
    D=Draw(ind)
    dF=0.0_wp

    if (E(1)>=eF) then
      w=0
    else if (E(4)<eF) then
      w=LinTetraDelta(D,dF)
    else if (E(1)<eF .and. eF<=E(2)) then
      call subtetra1(E,eF,lininterp1,vol1)
      w=matmul(transpose(lininterp1),LinTetraDelta(matmul(lininterp1,D),dF))*vol1
    else if (E(2)<eF .and. eF<=E(3)) then
      call subtetra2(E,eF,lininterp21,lininterp22,lininterp23,vol2)
      w=matmul(transpose(lininterp21),LinTetraDelta(matmul(lininterp21,D),dF))*vol2(1)+&
        matmul(transpose(lininterp22),LinTetraDelta(matmul(lininterp22,D),dF))*vol2(2)+&
        matmul(transpose(lininterp23),LinTetraDelta(matmul(lininterp23,D),dF))*vol2(3)
    else if (E(3)<eF .and. eF<=E(4)) then
      call subtetra3(E,eF,lininterp3,vol3)
      w=LinTetraDelta(D,dF)-matmul(transpose(lininterp3),LinTetraDelta(&
        matmul(lininterp3,D),dF))*vol3
    end if

    res(ind)=w
  end function LinTetraThetaDelta

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

! weight for x: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z)
! satisfied property: RawExp(x,y,z,w) == RawExp(x,y,w,z)
! singular when x->y , y->z , y->w
  function RawExp(x,y,z,w)
    implicit none
    real(wp), intent(in) :: x,y,z,w
    real(kind=16) :: a,b,c,d,RawExp
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(x,y,rtol) .or. isapprox(y,z,rtol) .or. isapprox(y,w,rtol)) then
        a=x
        b=y
        c=z
        d=w
        RawExp = ((-a/9 + b/4)*a**2 - (5*b**3)/36 + (a/3- b/2)*a**2*log(abs(a))&
             + b**3/6*log(abs(b)))/((a - b)**2*(b - c)*(b - d))
    else
        Rawexp=real(((-x/9 + y/4)*x**2 - (5*y**3)/36 + (x/3- y/2)*x**2*log(abs(x))+&
             y**3/6*log(abs(y)))/((x - y)**2*(y - z)*(y - w)),16)
    end if
  end function Rawexp

! limit case: RawExp(x,y,z,w), x->y
  function ExpXeqY(x,y,z,w)
    implicit none
    real(wp), intent(in) :: x,y,z,w
    real(kind=16) :: b,c,d,ExpXeqY
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(y,z,rtol) .or. isapprox(y,w,rtol)) then
        b=y
        c=z
        d=w
        ExpXeqY=(b*log(abs(b)))/(2*(b - c)*(b - d))
    else
        ExpXeqY=real((y*log(abs(y)))/(2*(y - z)*(y - w)),16)
    end if
  end function ExpXeqY

! limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y), y->z
  function ExpYeqZ(x,y,z,w)
    implicit none
    real(wp), intent(in) :: x,y,z,w
    real(kind=16) :: a,c,d,ExpYeqZ
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(x,z,rtol) .or. isapprox(z,w,rtol)) then
        a=x
        c=z
        d=w
        ExpYeqZ=((a-c)*(4*a**3+10*a*c*(c-d)+c**2*(6*c-d)-a**2*(8*c+d))-&
                6*a**2*(2*a**2-6*a*c+6*c**2+a*d-3*c*d)*log(abs(a))+&
                6*c**2*(2*a*c-3*a*d+c*d)*log(abs(c)))/(36*(a-c)**3*(c-d)**2)
    else
        ExpYeqZ=real(((x-z)*(4*x**3+10*x*z*(z-w)+z**2*(6*z-w)-x**2*(8*z+w))-6*x**2*&
              (2*x**2-6*x*z+6*z**2+x*w-3*z*w)*log(abs(x))+6*z**2*&
              (2*x*z-3*x*w+z*w)*log(abs(z)))/(36*(x-z)**3*(z-w)**2),16)
    end if
  end function ExpYeqZ

! limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z), y,z->w
  function ExpYZeqW(x,y,z,w)
    implicit none
    real(wp), intent(in) :: x,y,z,w
    real(kind=16) :: a,d,ExpYZeqW
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(x,w,rtol)) then
        a=x
        d=w
        ExpYZeqW=(2*a**3 + 3*a**2*d - 6*a*d**2 + d**3 - 6*a**2*d*log(abs(a)) +&
                  6*a**2*d*log(abs(d)))/(12*(a - d)**4)
    else
        ExpYZeqW=real((2*x**3+3*x**2*w-6*x*w**2+w**3-6*x**2*w*log(abs(x))+6*x**2*w*&
               log(abs(w)))/(12*(x-w)**4),16)
    end if
  end function ExpYZeqW

! limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y), x,y->z
  function ExpXYeqZ(x,y,z,w)
    implicit none
    real(wp), intent(in) :: x,y,z,w
    real(kind=16) :: c,d,ExpXYeqZ
    real(wp) :: rtol

    rtol=0.00001_wp
    if (isapprox(z,w,rtol)) then
        c=z
        d=w
        ExpXYeqZ=-((-c + d + (2*c + d)*log(abs(c)))/(6*(c - d)**2))
    else
        ExpXYeqZ=real(-((-z+w+(2*z+w)*log(abs(z)))/(6*(z-w)**2)),16)
    end if
  end function ExpXYeqZ

! weight for x: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z)
  function FracTetraWeight(x,y,z,w) result(res)
    real(wp), intent(in) :: x,y,z,w
    real(wp), dimension(3) :: yzw
    real(wp), dimension(4) :: v
    integer, dimension(3) :: ind
    real(wp) :: rtol,res

    rtol=0.00000001_wp
    yzw=(/y,z,w/)
    ind=sortperm(yzw)
    v=(/x,yzw(ind(1)),yzw(ind(2)),yzw(ind(3))/)

    if (isapprox(v(1),v(2),rtol)) then
      if (isapprox(v(2),v(3),rtol)) then
        if (isapprox(v(3),v(4),rtol)) then
            res=1/v(1)/24
        else
            res = ExpXYeqZ(v(1),v(2),v(3),v(4))+RawExp(v(1),v(4),v(2),v(3))
        end if
      else
        if (isapprox(v(3),v(4),rtol)) then
            res = ExpXeqY(v(1),v(2),v(3),v(4))+ExpYeqZ(v(1),v(3),v(4),v(2))
        else
            res = ExpXeqY(v(1),v(2),v(3),v(4))+RawExp(v(1),v(3),v(4),v(2))+RawExp(v(1),v(4),v(2),v(3))
        end if
      end if
    else if (isapprox(v(1),v(3),rtol)) then
      if (isapprox(v(3),v(4),rtol)) then
          res = ExpXYeqZ(v(1),v(3),v(4),v(2))+RawExp(v(1),v(2),v(3),v(4))
      else
          res = ExpXeqY(v(1),v(3),v(4),v(2))+RawExp(v(1),v(2),v(3),v(4))+RawExp(v(1),v(4),v(2),v(3))
      end if
    else if (isapprox(v(1),v(4),rtol)) then
      if (isapprox(v(2),v(3),rtol)) then
          res = ExpXeqY(v(1),v(4),v(2),v(3))+ExpYeqZ(v(1),v(2),v(3),v(4))
      else
          res = ExpXeqY(v(1),v(4),v(2),v(3))+RawExp(v(1),v(2),v(3),v(4))+RawExp(v(1),v(3),v(4),v(2))
      end if
    else if (isapprox(v(2),v(3),rtol)) then
      if (isapprox(v(3),v(4),rtol)) then
          res = ExpYZeqW(v(1),v(2),v(3),v(4))
      else
          res = ExpYeqZ(v(1),v(2),v(3),v(4)) + RawExp(v(1),v(4),v(2),v(3))
      end if
    else if (isapprox(v(3),v(4),rtol)) then
        res = ExpYeqZ(v(1),v(3),v(4),v(2)) + RawExp(v(1),v(2),v(3),v(4))
    else
        res = RawExp(v(1),v(2),v(3),v(4))+RawExp(v(1),v(3),v(4),v(2))+RawExp(v(1),v(4),v(2),v(3))
    end if
  end function FracTetraWeight

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = 1/𝔇(k).
  function LinTetraFrakD(FrakD)
    implicit none
    real(wp), dimension(4), intent(in) :: FrakD
    real(wp), dimension(4) :: LinTetraFrakD

    LinTetraFrakD=(/FracTetraWeight(FrakD(1),FrakD(2),FrakD(3),FrakD(4)),&
                    FracTetraWeight(FrakD(2),FrakD(3),FrakD(4),FrakD(1)),&
                    FracTetraWeight(FrakD(3),FrakD(4),FrakD(1),FrakD(2)),&
                    FracTetraWeight(FrakD(4),FrakD(1),FrakD(2),FrakD(3))/)
  end function LinTetraFrakD

! Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅ 1/D(k).
  function LinTetraThetaFrakD(Eraw,eF,Dnom) result(res)
    implicit none
    real(wp), dimension(4), intent(in) :: Eraw, Dnom
    real(wp), intent(in) :: eF
    integer, dimension(4) :: ind
    real(wp), dimension(4) :: E,D,w,res
    real(wp) :: vol1, vol3
    real(wp), dimension(3)::vol2
    real(wp), dimension(4,4) :: lininterp1, lininterp21, lininterp22,&
                                lininterp23, lininterp3

    ind=sortperm(Eraw)
    E=Eraw(ind)
    D=Dnom(ind)

    if (E(1)>=eF) then
      w=0
    else if (E(4)<eF) then
      w=LinTetraFrakD(D)
    else if (E(1)<eF .and. E(2)>=eF) then
      call subtetra1(E,eF,lininterp1,vol1)
      w=matmul(transpose(lininterp1),LinTetraFrakD(matmul(lininterp1,D)))*vol1
    else if (E(2)<eF .and. E(3)>=eF) then
      call subtetra2(E, eF, lininterp21, lininterp22, lininterp23, vol2)
      w=matmul(transpose(lininterp21),LinTetraFrakD(matmul(lininterp21,D)))*vol2(1)+&
          matmul(transpose(lininterp22),LinTetraFrakD(matmul(lininterp22,D)))*vol2(2)+&
          matmul(transpose(lininterp23),LinTetraFrakD(matmul(lininterp23,D)))*vol2(3)
    else if (E(3)<eF .and. E(4)>=eF) then
      call subtetra3(E,eF,lininterp3,vol3)
      w=LinTetraFrakD(D)-matmul(transpose(lininterp3),LinTetraFrakD(matmul(lininterp3,D)))*vol3
    end if

    res(ind)=w
  end function LinTetraThetaFrakD
end module LinTetra
