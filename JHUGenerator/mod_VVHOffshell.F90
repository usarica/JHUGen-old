module ModVVHOffshell
use ModParameters
implicit none

private

   real(dp), private, parameter :: tol = 0.00000010_dp
   real(8), parameter :: T3_nulep(1:2)= (/ 0.5d0, -0.5d0 /)
   real(8), parameter :: T3_ud(1:2)= (/ 0.5d0, -0.5d0 /)
   real(8), parameter :: Q_nulep(1:2)= (/ 0d0, -1d0 /)
   real(8), parameter :: Q_ud(1:2)= (/ 0.66666666666666666666666666667d0, -0.33333333333333333333333333333d0 /)
   real(8), parameter :: xw = sitW**2
   real(8), parameter :: cw = dsqrt(1d0-xw)
   real(8), parameter :: eA = dsqrt(4d0*pi*alpha_QED)
   real(8), parameter :: eW = eA/sitW ! g=e/sintW
   real(8), parameter :: eZ = eW*cw ! g*costw
   real(8), parameter :: eW_cW = eW/cw ! g/costw
   integer, parameter :: isW=1,isZ=2,isA=0 ! W+:1,W-:-1,Z:2,A:0


!----- notation for subroutines
!  public :: EvalAmp_VVHOffshell

contains



!  subroutine EvalAmp_VHiggs(id,helicity,MomExt,inv_mass,mass,me2)

! Use Weyl basis to make Psi=(Psi_R,Psi_L)
function gammaMatrix()
   implicit none
   complex(dp) :: gammaMatrix(1:4,1:4,1:4)

   gammaMatrix(:,:,:)=0_dp
   gammaMatrix(1,1,3)=1_dp ! gamma^0
   gammaMatrix(1,2,4)=1_dp
   gammaMatrix(1,3,1)=1_dp
   gammaMatrix(1,4,2)=1_dp
   gammaMatrix(2,1,4)=-1_dp ! gamma^x
   gammaMatrix(2,2,3)=-1_dp
   gammaMatrix(2,3,2)=1_dp
   gammaMatrix(2,4,1)=1_dp
   gammaMatrix(3,1,4)=ci ! gamma^y
   gammaMatrix(3,2,3)=-ci
   gammaMatrix(3,3,2)=-ci
   gammaMatrix(3,4,1)=ci
   gammaMatrix(3,1,3)=-1_dp ! gamma^z
   gammaMatrix(3,2,4)=1_dp
   gammaMatrix(3,3,1)=1_dp
   gammaMatrix(3,4,2)=-1_dp
   return
end function gammaMatrix

function MultiplySquareMatrices(ND,m1,m2)
   implicit none
   integer, intent(in) :: ND
   complex(dp), intent(in) :: m1(1:ND,1:ND),m2(1:ND,1:ND)
   complex(dp) :: MultiplySquareMatrices(1:ND,1:ND)
   integer :: i,mu,nu

   MultiplySquareMatrices=czero
   if(ND.lt.1) return
   do mu=1,ND
   do nu=1,ND
   do i=1,ND
      MultiplySquareMatrices(mu,nu) = MultiplySquareMatrices(mu,nu) + m1(mu,i)*m2(i,nu)
   enddo
   enddo
   enddo

end function MultiplySquareMatrices

function gamma_gamma_mu(N,p,ubar,v)
   implicit none
   integer, intent(in) :: N
   complex(dp), intent(in) :: p(1:4,2:N)
   complex(dp), intent(in) :: ubar(1:4),v(1:4)
   complex(dp) :: gamma_gamma_mu(1:4)
   integer :: a,b,i,mu
   complex(dp) :: gamma(1:4,1:4,1:4,1:N),gamma_tmp(1:4,1:4,2:N),gtmp(1:4,1:4,1:4),gacc(1:4,1:4,1:4)

   gamma_gamma_mu=czero
   if(N.lt.2) return

   gamma(:,:,:,1)=gammaMatrix()
   do i=2,N
      gamma(:,:,:,i)=gamma(:,:,:,1)
   enddo
   do i=2,N
      gamma_tmp(:,:,i) = gamma(1,:,:,i)*p(1,i) - gamma(2,:,:,i)*p(2,i) - gamma(3,:,:,i)*p(3,i) - gamma(4,:,:,i)*p(4,i) ! gama^mu p_mu
   enddo

   gacc=gamma(:,:,:,1)
   do i=2,N
   do mu=1,4
      gtmp(mu,:,:) = MultiplySquareMatrices(4,gacc(mu,:,:),gamma_tmp(:,:,i))
      gacc=gtmp
   enddo
   enddo

   do a=1,4
   do b=1,4
      gamma_gamma_mu(:) = gamma_gamma_mu(:) + ubar(a)*gacc(:,a,b)*v(b)
   enddo
   enddo
   return
end function gamma_gamma_mu

!     ubar spinor, massless
function ubar0(p,i)
   implicit none
   complex(dp), intent(in) :: p(4)
   integer, intent(in) :: i
   complex(dp) :: ubar0(4)
   complex(dp) :: fc, fc2
   real(dp)   :: p0,px,py,pz

   p0=real(p(1),dp)
   px=real(p(2),dp)
   py=real(p(3),dp)
   pz=real(p(4),dp)

   fc2 = p0 + pz
   fc=sqrt(fc2)

   if (abs(fc2).gt. tol) then
      if (i.eq.1) then
        ubar0(1)=czero
        ubar0(2)=czero
        ubar0(3)=fc
        ubar0(4)=(px-ci*py)/fc
      elseif (i.eq.-1) then
        ubar0(1)=(px+ci*py)/fc
        ubar0(2)=-fc
        ubar0(3)=czero
        ubar0(4)=czero
      else
        stop 'ubar0: i out of range'
      endif
   else
      if (i.eq.1) then
        ubar0(1) = czero
        ubar0(2) = czero
        ubar0(3) = czero
        ubar0(4) = sqrt(cone*two*p0)
      elseif (i.eq.-1) then
        ubar0(1) = sqrt(cone*two*p0)
        ubar0(2) = czero
        ubar0(3) = czero
        ubar0(4) = czero
      else
        stop 'ubar0: i out of range'
      endif
   endif
   return
end function ubar0


!     ubar spinor, massive
function ubarm(cp,i)
   implicit none
   complex(dp), intent(in) :: cp(4)
   integer, intent(in) :: i
   complex(dp) :: ubarm(4)
   complex(dp) :: fc2
   real(dp)   :: E,x,y,z,m,p,fac

   E=real(cp(1),dp)
   x=real(cp(2),dp)
   y=real(cp(3),dp)
   z=real(cp(4),dp)
   m=dsqrt(dabs(E**2-x**2-y**2-z**2))
   p=dsqrt(dabs(x**2+y**2+z**2))

   fc2 = p + z

   if (abs(fc2).gt. tol) then
      fac = cone/sqrt(2_dp*p*fc2)
      if (i.eq.1) then
         ubarm(1)=sqrt(E-p)*fac*fc2
         ubarm(2)=sqrt(E-p)*fac*(x-ci*y)
         ubarm(3)=sqrt(E+p)*fac*fc2
         ubarm(4)=sqrt(E+p)*fac*(x-ci*y)
      elseif (i.eq.-1) then
         ubarm(1)=sqrt(E+p)*fac*(x+ci*y)
         ubarm(2)=-sqrt(E+p)*fac*fc2
         ubarm(3)=sqrt(E-p)*fac*(x+ci*y)
         ubarm(4)=-sqrt(E-p)*fac*fc2
      else
        stop 'ubarm: i out of range'
      endif
   else
      if (i.eq.1) then
        ubarm(1) = czero
        ubarm(2) = sqrt(cone*(E-p))
        ubarm(3) = czero
        ubarm(4) = sqrt(cone*(E+p))
      elseif (i.eq.-1) then
        ubarm(1) = sqrt(cone*(E+p))
        ubarm(2) = czero
        ubarm(3) = sqrt(cone*(E-p))
        ubarm(4) = czero
      else
        stop 'ubarm: i out of range'
      endif
   endif
   return
end function ubarm


  ! -- v0  spinor, massless
function v0(p,i)
   implicit none
   complex(dp), intent(in) :: p(4)
   integer, intent(in)      :: i
   complex(dp) :: v0(4)
   complex(dp) :: fc2, fc
   real(dp)   :: p0,px,py,pz

   p0=real(p(1),dp)
   px=real(p(2),dp)
   py=real(p(3),dp)
   pz=real(p(4),dp)

   fc2 = p0 + pz
   fc=sqrt(fc2)

   if (abs(fc2).gt. tol) then
      if (i.eq.1) then
        v0(1)=czero
        v0(2)=czero
        v0(3)=(px-ci*py)/fc
        v0(4)=-fc
      elseif (i.eq.-1) then
        v0(1)=fc
        v0(2)=(px+ci*py)/fc
        v0(3)=czero
        v0(4)=czero
      else
        stop 'v0: i out of range'
      endif
   else
      if (i.eq.1) then
        v0(1)=czero
        v0(2)=czero
        v0(3)=sqrt(cone*two*p0)
        v0(4)=czero
      elseif (i.eq.-1) then
        v0(1)=czero
        v0(2)=sqrt(cone*two*p0)
        v0(3)=czero
        v0(4)=czero
      else
        stop 'v0: i out of range'
      endif
   endif
   return
end function v0


!     v spinor, massive
function vm(cp,i)
   implicit none
   complex(dp), intent(in) :: cp(4)
   integer, intent(in) :: i
   complex(dp) :: vm(4)
   complex(dp) :: fc2
   real(dp)   :: E,x,y,z,m,p,fac

   E=real(cp(1),dp)
   x=real(cp(2),dp)
   y=real(cp(3),dp)
   z=real(cp(4),dp)
   m=dsqrt(dabs(E**2-x**2-y**2-z**2))
   p=dsqrt(dabs(x**2+y**2+z**2))

   fc2 = p + z

   if (abs(fc2).gt. tol) then
      fac = cone/sqrt(2_dp*p*fc2)
      if (i.eq.1) then
         vm(1)=-sqrt(E-p)*fac*(x-ci*y)
         vm(2)=sqrt(E-p)*fac*fc2
         vm(3)=sqrt(E+p)*fac*(x-ci*y)
         vm(4)=-sqrt(E+p)*fac*fc2
      elseif (i.eq.-1) then
         vm(1)=sqrt(E+p)*fac*fc2
         vm(2)=sqrt(E+p)*fac*(x+ci*y)
         vm(3)=-sqrt(E-p)*fac*fc2
         vm(4)=-sqrt(E-p)*fac*(x+ci*y)
      else
        stop 'vm: i out of range'
      endif
   else
      if (i.eq.1) then
        vm(1) = -sqrt(cone*(E-p))
        vm(2) = czero
        vm(3) = sqrt(cone*(E+p))
        vm(4) = czero
      elseif (i.eq.-1) then
        vm(1) = czero
        vm(2) = sqrt(cone*(E+p))
        vm(3) = czero
        vm(4) = -sqrt(cone*(E-p))
      else
        stop 'vm: i out of range'
      endif
   endif
   return
end function vm


function CurrentPrefactor(id,hel,useA)
   implicit none
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   integer :: idV
   logical, optional :: useA
   complex(dp) :: CurrentPrefactor
   integer :: testAcoupl
   integer :: idc(1:2)

   CurrentPrefactor = czero
   if(hel.eq.0) return
   testAcoupl=0
   if(present(useA)) then
      if(useA) testAcoupl=1
   endif
   idV = CoupledVertex(id,hel,useAHcoupl=testAcoupl)
   if(idV.eq.Not_a_particle_) return
   idc(1)=convertLHE(id(1))
   idc(2)=convertLHE(id(2))

   if( idV.eq.Wp_ .or. idV.eq.Wm_ ) then ! W->ffb'
      CurrentPrefactor = ci*eW/sqrt2*CKM(idc(1),idc(2))/ScaleFactor(idc(1),idc(2))
   elseif( ((abs(idc(1)).eq.12) .or. (abs(idc(1)).eq.14) .or. (abs(idc(1)).eq.16)) .and. testAcoupl.eq.0 ) then ! Z->nunu
      CurrentPrefactor = ci*eW_cW*(T3_nulep(1)-Q_nulep(1)*xw)
   elseif( (abs(idc(1)).eq.11) .or. (abs(idc(1)).eq.13) .or. (abs(idc(1)).eq.15) ) then ! A/Z->ll
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_nulep(2)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_nulep(2)-Q_nulep(2)*xw)
      else
         CurrentPrefactor = ci*eW_cW*(-Q_nulep(2)*xw)
      endif
   elseif( (abs(idc(1)).eq.pdfUp_) .or. (abs(idc(1)).eq.pdfChm_) .or. (abs(idc(1)).eq.pdfTop_) ) then ! Z->uu
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_ud(1)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_ud(1)-Q_ud(1)*xw)
      else
         CurrentPrefactor = ci*eW_cW*(-Q_ud(1)*xw)
      endif
   elseif( (abs(idc(1)).eq.pdfDn_) .or. (abs(idc(1)).eq.pdfStr_) .or. (abs(idc(1)).eq.pdfBot_) ) then ! Z->uu
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_ud(2)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_ud(2)-Q_ud(2)*xw)
      else
         CurrentPrefactor = ci*eW_cW*(-Q_ud(2)*xw)
      endif
   endif
   return
end function CurrentPrefactor

! Returns i*charge*J^mu_L/R
function Vcurrent(p,id,hel,idV,useAcoupl)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2)
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   integer, intent(out) :: idV
   logical, optional :: useAcoupl
   integer :: idc(1:2)
   complex(dp) :: Vcurrent(4),Ub(4),V(4),prefactor
   integer :: testAcoupl

   Vcurrent(:)=czero

   if(hel.eq.0) return
   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1
   endif
   idV = CoupledVertex(id,hel,useAHcoupl=testAcoupl)
   if(idV.eq.Not_a_particle_) return

   idc(1)=convertLHE(id(1))
   idc(2)=convertLHE(id(2))
   if(idc(1).gt.0 .and. idc(2).lt.0) then
      Ub(:)=ubar0(cmplx(p(1:4,1),kind=dp),hel)
      V(:)=v0(cmplx(p(1:4,2),kind=dp),-hel)
   else
      Ub(:)=ubar0(cmplx(p(1:4,2),kind=dp),hel)
      V(:)=v0(cmplx(p(1:4,1),kind=dp),-hel)
   endif
   ! 1=E,2=px,3=py,4=pz
   ! This is an expression for Ub(+/-)) Gamma^\mu V(-/+)
   Vcurrent(1)=(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
   Vcurrent(2)=(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
   Vcurrent(3)=ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
   Vcurrent(4)=(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))

   if(present(useAcoupl)) then
      prefactor=CurrentPrefactor(id,hel,useA=useAcoupl)
   else
      prefactor=CurrentPrefactor(id,hel)
   endif
   Vcurrent(:)=Vcurrent(:)*prefactor
   return
end function Vcurrent

function pol_massive(p,hel)
   implicit none
   real(dp), intent(in) :: p(1:4)
   integer, intent(in) :: hel
   real(dp) :: sincos(1:4), abs3p
   complex(dp) :: inv_mass
   complex(dp) :: pol_massive(1:4)
   real(dp) :: epsilon = 1d-14*GeV ! A small quantity ~ 1e-5 eV
      
   pol_massive(:)=czero
   sincos = angles(p)

   ! Transverse:
   ! pol_massive(1) ! 0
   ! pol_massive(2) ! 1/sqrt(2) * { cos(theta)*cos(phi) - lambda*i*sin(phi) }
   ! pol_massive(3) ! 1/sqrt(2) * { cos(theta)*sin(phi) + lambda*i*cos(phi) }
   ! pol_massive(4) ! -1/sqrt(2) * sin(theta)
   ! Longitudinal:
   ! pol_massive(1) ! |p|/m
   ! pol_massive(2) ! E/m*px/|p|
   ! pol_massive(3) ! E/m*py/|p|
   ! pol_massive(4) ! E/m*pz/|p|

   if(hel.ne.0) then ! lambda = +-1
      pol_massive(1)= 0d0
      pol_massive(2)= (sincos(3)*sincos(1)-ci*hel*sincos(4))/sqrt2
      pol_massive(3)= (sincos(4)*sincos(1)+ci*hel*sincos(3))/sqrt2
      pol_massive(4)= -sincos(2)/sqrt2
   else ! lambda = L
      abs3p = sqrt(abs(p(2)**2+p(3)**2+p(4)**2))
      inv_mass = (p(1)**2-abs3p**2)
      inv_mass = sqrt(inv_mass) ! make sure E**2-p**2=m**2 is satisfied by extending to the complex notation since sum_h{eps^mu_h eps*^nu_h}=-(g^munu-p^mup^nu/p**2) has to be true.
      if(abs(inv_mass).gt.epsilon) then
         pol_massive(1)= abs3p/inv_mass
         pol_massive(2)= sincos(3)*sincos(2)*p(1)/inv_mass
         pol_massive(3)= sincos(4)*sincos(2)*p(1)/inv_mass
         pol_massive(4)= sincos(1)*p(1)/inv_mass
      endif
   endif
   return
end function pol_massive

function pol_massless(p,hel)
   implicit none
   real(dp), intent(in) :: p(1:4)
   integer, intent(in) :: hel
   real(dp) :: sincos(1:4)
   complex(dp) :: pol_massless(1:4)
      
   pol_massless(:)=czero
   sincos = angles(p)

   ! pol_massless(1) ! 0
   ! pol_massless(2) ! 1/sqrt(2) * { cos(theta)*cos(phi) - lambda*i*sin(phi) }
   ! pol_massless(3) ! 1/sqrt(2) * { cos(theta)*sin(phi) + lambda*i*cos(phi) }
   ! pol_massless(4) ! -1/sqrt(2) * sin(theta)

   if(hel.ne.0) then ! lambda = +-1
      pol_massless(1)= 0d0
      pol_massless(2)= (sincos(3)*sincos(1)-ci*hel*sincos(4))/sqrt2
      pol_massless(3)= (sincos(4)*sincos(1)+ci*hel*sincos(3))/sqrt2
      pol_massless(4)= -sincos(2)/sqrt(2d0)
   endif
   return
end function pol_massless

function angles(vector, halfangles)
   implicit none
   real(dp), intent(in) :: vector(4)
   logical, optional :: halfangles
   logical :: usehalftheta
   real(dp) :: angles(4)
   real(dp) :: epsilon = 1d-14*GeV ! A small quantity ~ 1e-5 eV
   real(dp) :: abs3p, abspt

   ! angles(1)=cos(theta)
   ! angles(2)=sin(theta)
   ! angles(3)=cos(phi)
   ! angles(4)=sin(phi)

   usehalftheta=.false.
   if(present(halfangles)) then
      usehalftheta=halfangles
   endif

   angles(:)=0_dp

   abs3p = sqrt(abs(vector(2)**2+vector(3)**2+vector(4)**2))
   if(abs3p.lt.epsilon) then
     angles(1)=1_dp ! theta=0
   else
     angles(1)=vector(4)/abs3p
     angles(2)=sqrt(abs((1d0+angles(1))*(1d0-angles(1))))
     if(usehalftheta) then
        angles(2)=sqrt(abs((1_dp-angles(1))/2_dp))
        angles(1)=sqrt(abs((1_dp+angles(1))/2_dp))
     endif
   endif

   abspt = sqrt(abs(vector(2)**2+vector(3)**2))
   if(abspt.lt.epsilon) then
      angles(3)=1_dp
   else
      angles(3)=vector(2)/abspt
      angles(4)=vector(3)/abspt
   endif
   return
end function angles


function dgmunuab_ewk(p) ! 2 g_munu g_ab - g_mua g_nub -g_mub g_nua
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4) ! Lorentz, current index
   complex(dp) :: dgmunuab_ewk

   dgmunuab_ewk = 2_dp*(p(:,1).dot.p(:,2))*(p(:,3).dot.p(:,4))-(p(:,1).dot.p(:,3))*(p(:,2).dot.p(:,4))-(p(:,1).dot.p(:,4))*(p(:,2).dot.p(:,3))
   return
end function dgmunuab_ewk

function dp_gmunu_ewk(p,current) ! (p+-p-)_mu g_ab + (p- -q)_a g_mub + (q-p+)_b g_mua
   use ModMisc
   implicit none
   complex(dp), intent(in) :: current(1:4,1:3) ! Lorentz, current Z/A W+-
   real(dp), intent(in) :: p(1:4,1:2) ! Lorentz, mom. W+-
   real(dp) :: q(1:4)
   complex(dp) :: dpm(1:4),dmq(1:4),dqp(1:4)
   complex(dp) :: dp_gmunu_ewk

   q(:) = -p(:,1)-p(:,2)
   dpm(:) = dcmplx(p(:,2)-p(:,3))
   dmq(:) = dcmplx(-q(:)+p(:,3))
   dqp(:) = dcmplx(q(:)-p(:,2))

   dp_gmunu_ewk = (dpm(:).dot.current(:,1))*(current(:,2).dot.current(:,3)) &
                + (dmq(:).dot.current(:,2))*(current(:,3).dot.current(:,1)) &
                + (dqp(:).dot.current(:,3))*(current(:,1).dot.current(:,2))
   return
end function dp_gmunu_ewk

function epsilon_munuab(e1,e2,e3,e4)
   implicit none
   complex(dp), intent(in) :: e1(4), e2(4), e3(4), e4(4)
   complex(dp) :: epsilon_munuab

   epsilon_munuab =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3)  &
                     -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
                     +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
                     -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
                     +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
                     -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
                     +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
                     -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
                     +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
                     -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
                     +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
                     -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)

   return
end function epsilon_munuab

function epsilon_munuab_single(mu,nu,a,b)
   implicit none
   integer mu,nu,a,b
   real(dp) epsilon_munuab_single
   epsilon_munuab_single=real((mu-nu)*(mu-a)*(mu-b)*(nu-a)*(nu-b)*(a-b),kind=dp)/12_dp
   return
end function epsilon_munuab_single

function Id_Order(Npart,idV,idT)
   implicit none
   integer, intent(in) :: Npart
   integer, intent(in) :: idV(1:NPart),idT(1:NPart)
   integer :: Id_Order(1:Npart),target,part

   Id_Order(:)=0
   
   do part=1,NPart
      do target=1,NPart
         if(idV(part).eq.idT(target) .and. Id_Order(target).ne.0) then
            Id_Order(target)=part
            exit
         endif
      enddo
   enddo

   return
end function Id_Order


function CoupledVertex(id,hel,useAHcoupl)
   implicit none
   integer, optional :: useAHcoupl
   integer, intent(in) :: id(1:2),hel
   integer :: testAHcoupl
   integer :: CoupledVertex

   testAHcoupl = 0
   if(present(useAHcoupl)) then
      testAHcoupl = useAHcoupl
   endif
   if( (&
   (id(1).eq.ElP_ .and. id(2).eq.NuE_) .or. (id(2).eq.ElP_ .and. id(1).eq.NuE_) .or. &
   (id(1).eq.MuP_ .and. id(2).eq.NuM_) .or. (id(2).eq.MuP_ .and. id(1).eq.NuM_) .or. &
   (id(1).eq.TaP_ .and. id(2).eq.NuT_) .or. (id(2).eq.TaP_ .and. id(1).eq.NuT_) .or. &
   (id(1).eq.Up_  .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Up_  .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Chm_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Chm_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Top_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Top_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wp_
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ANuE_) .or. (id(2).eq.ElM_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.MuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.ANuT_) .or. (id(2).eq.TaM_ .and. id(1).eq.ANuT_) .or. &
   (id(1).eq.AUp_  .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AUp_  .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.AChm_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AChm_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.ATop_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.ATop_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wm_
   elseif( &
   (id(1).eq.ElM_ .and. id(2).eq.ElP_) .or. (id(2).eq.ElM_ .and. id(1).eq.ElP_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.MuP_) .or. (id(2).eq.MuM_ .and. id(1).eq.MuP_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.TaP_) .or. (id(2).eq.TaM_ .and. id(1).eq.TaP_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.AUp_) .or. (id(2).eq.Up_  .and. id(1).eq.AUp_) .or. &
   (id(1).eq.Dn_  .and. id(2).eq.ADn_) .or. (id(2).eq.Dn_  .and. id(1).eq.ADn_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AChm_) .or. (id(2).eq.Chm_ .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.AStr_) .or. (id(2).eq.Str_ .and. id(1).eq.Astr_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.ATop_) .or. (id(2).eq.Top_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ABot_) .or. (id(2).eq.Bot_ .and. id(1).eq.ABot_)      &
   ) then
      if(testAHcoupl.eq.1) then
         CoupledVertex=Pho_
      elseif(testAHcoupl.eq.2) then
         CoupledVertex=Hig_
      else
         CoupledVertex=Z0_
      endif
   elseif( (&
   (id(1).eq.NuE_ .and. id(2).eq.ANuE_) .or. (id(2).eq.NuE_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.NuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.NuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.NuT_ .and. id(2).eq.ANuT_) .or. (id(2).eq.NuT_ .and. id(1).eq.ANuT_)      &
   ) .and. hel.lt.0) then
      CoupledVertex=Z0_ ! Only Z coupling to nuL-nubR
   else
      CoupledVertex=Not_a_particle_
   endif

   return
end function CoupledVertex

function CoupledVertex_FlavorViolating(id,hel,useAHcoupl)
   implicit none
   integer, optional :: useAHcoupl
   integer, intent(in) :: id(1:2),hel
   integer :: testAHcoupl
   integer :: CoupledVertex_FlavorViolating

   testAHcoupl = 0
   if(present(useAHcoupl)) then
      testAHcoupl = useAHcoupl
   endif
   if( (&
   (id(1).eq.Up_  .and. id(2).eq.AChm_) .or. (id(2).eq.Up_  .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.ATop_) .or. (id(2).eq.Up_  .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AUp_ ) .or. (id(2).eq.Chm_ .and. id(1).eq.AUp_ ) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.ATop_) .or. (id(2).eq.Chm_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.AUp_ ) .or. (id(2).eq.Top_ .and. id(1).eq.AUp_ ) .or. &
   (id(1).eq.Top_ .and. id(2).eq.AChm_) .or. (id(2).eq.Top_ .and. id(1).eq.AChm_) .or. &

   (id(1).eq.Dn_  .and. id(2).eq.AStr_) .or. (id(2).eq.Dn_  .and. id(1).eq.AStr_) .or. &
   (id(1).eq.Dn_  .and. id(2).eq.ABot_) .or. (id(2).eq.Dn_  .and. id(1).eq.ABot_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.ADn_ ) .or. (id(2).eq.Str_ .and. id(1).eq.ADn_ ) .or. &
   (id(1).eq.Str_ .and. id(2).eq.ABot_) .or. (id(2).eq.Str_ .and. id(1).eq.ABot_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ADn_ ) .or. (id(2).eq.Bot_ .and. id(1).eq.ADn_ ) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.AStr_) .or. (id(2).eq.Bot_ .and. id(1).eq.AStr_)      &
   ) .and. hel.lt.0) then
      if(testAHcoupl.eq.1) then
         CoupledVertex_FlavorViolating=Pho_
      elseif(testAHcoupl.eq.2) then
         CoupledVertex_FlavorViolating=Hig_
      else
         CoupledVertex_FlavorViolating=Z0_
      endif
   else
      CoupledVertex_FlavorViolating=Not_a_particle_
   endif

   return
end function CoupledVertex_FlavorViolating

function WDaughterPair(id)
   implicit none
   integer, intent(in) :: id
   integer :: WDaughterPair(1:3)

   WDaughterPair(:)=Not_a_particle_
   if(id.eq.ElP_) then
      WDaughterPair(1)=NuE_
   elseif(id.eq.MuP_) then
      WDaughterPair(1)=NuM_
   elseif(id.eq.TaP_) then
      WDaughterPair(1)=NuT_
   elseif(id.eq.NuE_) then
      WDaughterPair(1)=ElP_
   elseif(id.eq.NuM_) then
      WDaughterPair(1)=MuP_
   elseif(id.eq.NuT_) then
      WDaughterPair(1)=TaP_

   elseif(id.eq.ElM_) then
      WDaughterPair(1)=ANuE_
   elseif(id.eq.MuM_) then
      WDaughterPair(1)=ANuM_
   elseif(id.eq.TaM_) then
      WDaughterPair(1)=ANuT_
   elseif(id.eq.ANuE_) then
      WDaughterPair(1)=ElM_
   elseif(id.eq.ANuM_) then
      WDaughterPair(1)=MuM_
   elseif(id.eq.ANuT_) then
      WDaughterPair(1)=TaM_

   elseif(id.eq.Up_ .or. id.eq.Chm_ .or. id.eq.Top_) then
      WDaughterPair(1)=ADn_
      WDaughterPair(2)=AStr_
      WDaughterPair(3)=ABot_
   elseif(id.eq.AUp_ .or. id.eq.AChm_ .or. id.eq.ATop_) then
      WDaughterPair(1)=Dn_
      WDaughterPair(2)=Str_
      WDaughterPair(3)=Bot_

   elseif(id.eq.Dn_ .or. id.eq.Str_ .or. id.eq.Bot_) then
      WDaughterPair(1)=AUp_
      WDaughterPair(2)=AChm_
!      WDaughterPair(3)=ATop_
   elseif(id.eq.ADn_ .or. id.eq.AStr_ .or. id.eq.ABot_) then
      WDaughterPair(1)=Up_
      WDaughterPair(2)=Chm_
!      WDaughterPair(3)=Top_
   endif

   return
end function WDaughterPair



! WWZZ,WWZA,WWAA,WWWW
function QuarticEWKVertex(p,idV,outType) ! p(:,1:4): Currents; order: (incoming W+,incoming W-, outgoing Z or A or W, outgoing Z or A or W); outType(1:2) for Z, A, W+-
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4) ! Lorentz, current index
   integer, intent(in) :: idV(1:4) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:4)
   complex(dp) :: QuarticEWKVertex

   QuarticEWKVertex = czero
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))
   elseif( (outType(1).eq.Pho_ .and. outType(2).eq.Z0_) .or. (outType(2).eq.Pho_ .and. outType(1).eq.Z0_) ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Pho_ /))
   elseif( outType(1).eq.Pho_ .and. outType(2).eq.Pho_ ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Pho_,Pho_ /))
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Z0_ /))
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
      print *,"QuarticEWKVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   if( outType(1).eq.Wp_ .or. outType(1).eq.Wm_ ) then ! W+W-W+W-
      QuarticEWKVertex = ci*dgmunuab_ewk( (/ p(1:4,order(1)), p(1:4,order(3)), p(1:4,order(2)), p(1:4,order(4)) /) )
      QuarticEWKVertex = QuarticEWKVertex*eW**2
   else ! W+W-(ZZ/ZA/AA); symmetric in munu and ab indices
      QuarticEWKVertex = -ci*dgmunuab_ewk( (/ p(1:4,order(1)), p(1:4,order(2)), p(1:4,order(3)), p(1:4,order(4)) /) )
      if(outType(1).eq.Z0_) then
         QuarticEWKVertex = QuarticEWKVertex*eZ
      else
         QuarticEWKVertex = QuarticEWKVertex*eA
      endif
      if(outType(2).eq.Z0_) then
         QuarticEWKVertex = QuarticEWKVertex*eZ
      else
         QuarticEWKVertex = QuarticEWKVertex*eA
      endif
   endif
   return
end function QuarticEWKVertex

function TripleEWKVertex(p,current,idV,useAcoupl)
   implicit none
   real(dp), intent(in) :: p(1:4,1:3)
   complex(dp), intent(in) :: current(1:4,1:3) ! Lorentz, particle index; p and current are for Z/A and W+W-
   integer, intent(in) :: idV(1:3)
   logical, optional :: useAcoupl
   logical :: forceAcoupl
   integer :: order(1:3)
   complex(dp) :: TripleEWKVertex

   TripleEWKVertex = czero
   forceAcoupl=.false.
   if(present(useAcoupl)) forceAcoupl = useAcoupl
   if(forceAcoupl) then
      order(:)=Id_Order(3,idV,(/ Pho_,Wp_,Wm_ /))
   else
      order(:)=Id_Order(3,idV,(/ Z0_ ,Wp_,Wm_ /))
   endif

   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
      print *,"TripleEWKVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   TripleEWKVertex = ci*dp_gmunu_ewk( (/ p(1:4,order(2)), p(1:4,order(3)) /), (/ current(1:4,order(1)), current(1:4,order(2)), current(1:4,order(3)) /) )

   if(forceAcoupl) then
      TripleEWKVertex = TripleEWKVertex*eA
   else
      TripleEWKVertex = TripleEWKVertex*eZ
   endif

   return
end function TripleEWKVertex

! A/Zff,f->fA/Zf'f'
function ZAf_fZAfpfp(p,id,hel,Ub,V,useAcoupl)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:2) ! Lorentz, 1,i=ubar(i); 2,i=v(i)
   integer, intent(in) :: id(1:2,1:2)
   integer, intent(in) :: hel(1:2)
   complex(dp), intent(in) :: Ub(1:4,1:2),V(1:4,1:2)
   logical, optional :: useAcoupl
   integer :: idV(1:2) ! V->f1 fb1 or V->f2 fb2
   integer :: idc(1:2,1:2)
   complex(dp) :: ZAf_fZAfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: testAcoupl,i,j
   real(dp) :: q(1:4),q1(1:4),q2(1:4)
   complex(dp) :: innercurrent_tmp(1:4),q1scprop,prefactor(1:2),innerprefactor
   integer :: idV_tmp

   if(hel(1).eq.0 .or. hel(2).eq.0) return
   do i=1,2;do j=1,2
   idc(i,j)=convertLHE(id(i,j))
   enddo;enddo
   if( idc(1,1).ne.-idc(2,1) .or. idc(1,2).ne.-idc(2,2) ) return

   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1 ! For the root of the diagram
      prefactor(1) = CurrentPrefactor(id(:,1),hel(1),useA=useAcoupl) ! A->f1f1 V->f2f2
      prefactor(2) = CurrentPrefactor(id(:,2),hel(2),useA=useAcoupl) ! A->f2f2 V->f1f1
   else
      prefactor(1) = CurrentPrefactor(id(:,1),hel(1)) ! Z->f1f1 V->f2f2
      prefactor(2) = CurrentPrefactor(id(:,2),hel(2)) ! Z->f2f2 V->f1f1
   endif

   ZAf_fZAfpfp(:) = czero
   q(:)=0_dp
   do i=1,2
      q(:)=q(:)+p(:,1,i)+p(:,2,i)
   enddo

   idV(1) = CoupledVertex(id(:,1),hel(1),useAHcoupl=testAcoupl)
   idV(2) = CoupledVertex(id(:,2),hel(2),useAHcoupl=testAcoupl)
   if(idV(1).ne.Not_a_particle_) then
      ! Case 1: V->f1 fb1, f1->f1 Z -> f2 fb2
      q2(:)=p(:,1,2)+p(:,2,2)
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1))
      else
         q1(:)=q(:)-p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1))
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_tmp,useAcoupl=.false.) ! idV_tmp=Z0_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      q1scprop = ScalarPropagator(id(1,1),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(1)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)

      ! Case 2: V->f1 fb1, fb1->fb1 Z -> f2 fb2
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1))
      else
         q1(:)=-q(:)+p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1))
      endif
      q1scprop = ScalarPropagator(id(1,1),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(1)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)
   endif

   if(idV(2).ne.Not_a_particle_) then
      ! Case 3: V->f2 fb2, f2->f2 Z -> f1 f1b
      q2(:)=p(:,1,1)+p(:,2,1)
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=q(:)-p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2))
      else
         q1(:)=q(:)-p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2))
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,1),id(1:2,1),hel(1),idV_tmp,useAcoupl=.false.) ! idV_tmp=Z0_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      q1scprop = ScalarPropagator(id(1,2),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(2)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)

      ! Case 4: V->f2 fb2, fb2->fb2 Z -> f1 fb1
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=-q(:)+p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2))
      else
         q1(:)=-q(:)+p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2))
      endif
      q1scprop = ScalarPropagator(id(1,2),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(2)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)
   endif

   ! Pho_
   ! Case 1: V->f1 fb1, f1->f1 A -> f2 f2b
   if( .not.((abs(idc(1,2)).eq.12) .or. (abs(idc(1,2)).eq.14) .or. (abs(idc(1,2)).eq.16)) .and. idV(1).ne.Not_a_particle_ ) then
      q2(:)=p(:,1,2)+p(:,2,2)
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1),useA=useAcoupl)
      else
         q1(:)=q(:)-p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1),useA=useAcoupl)
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_tmp,useAcoupl=.true.) ! idV_tmp=Pho_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      q1scprop = ScalarPropagator(id(1,1),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(1)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)

      ! Case 2: V->f1 fb1, fb1->fb1 A -> f2 f2b
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1),useA=useAcoupl)
      else
         q1(:)=-q(:)+p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1),useA=useAcoupl)
      endif
      q1scprop = ScalarPropagator(id(1,1),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(1)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)
   endif

   if( .not.((abs(idc(1,1)).eq.12) .or. (abs(idc(1,1)).eq.14) .or. (abs(idc(1,1)).eq.16)) .and. idV(2).ne.Not_a_particle_ ) then
      ! Case 3: V->f2 fb2, f2->f2 A -> f1 fb1
      q2(:)=p(:,1,1)+p(:,2,1)
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=q(:)-p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2),useA=useAcoupl)
      else
         q1(:)=q(:)-p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2),useA=useAcoupl)
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,1),id(1:2,1),hel(1),idV_tmp,useAcoupl=.true.) ! idV_tmp=Pho_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      q1scprop = ScalarPropagator(id(1,2),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(2)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)

      ! Case 4: V->f2 fb2, fb2->fb2 A -> f1 fb1
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=-q(:)+p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2),useA=useAcoupl)
      else
         q1(:)=-q(:)+p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2),useA=useAcoupl)
      endif
      q1scprop = ScalarPropagator(id(1,2),q1)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor(2)*innerprefactor
      ZAf_fZAfpfp(:) = ZAf_fZAfpfp(:) + restmp(:)
   endif

   return
end function ZAf_fZAfpfp


! A/Zff,f->fWf'f'
function ZAf_fWfpfp(p,id,hel,Ub,V,useAcoupl)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:2) ! Lorentz, 1,i=ubar(i); 2,i=v(i)
   integer, intent(in) :: id(1:2,1:2)
   integer, intent(in) :: hel(1:2)
   complex(dp), intent(in) :: Ub(1:4,1:2),V(1:4,1:2)
   logical, optional :: useAcoupl
   integer :: idV(1:2) ! V->f1 fb1 only
   integer :: idc(1:2,1:2)
   complex(dp) :: ZAf_fWfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: testAcoupl,i,j
   real(dp) :: q(1:4),q1(1:4),q2(1:4)
   complex(dp) :: innercurrent_tmp(1:4),q1scprop,prefactor,innerprefactor
   integer :: idV_tmp

   if(hel(1).eq.0 .or. hel(2).ne.-1) return
   do i=1,2;do j=1,2
   idc(i,j)=convertLHE(id(i,j))
   enddo;enddo
   if( idc(1,1).eq.-idc(2,1) .or. idc(1,1).ne.-sign(idc(1,1),idc(2,1)) .or. idc(1,2).eq.-idc(2,2) .or. idc(1,2).ne.-sign(idc(1,2),idc(2,2)) ) return ! The presence of W+- inside changes the ids

   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1 ! For the root of the diagram
   endif

   ZAf_fWfpfp(:) = czero
   q(:)=0_dp
   do i=1,2
      q(:)=q(:)+p(:,1,i)+p(:,2,i)
   enddo

   if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
      idV(1) = CoupledVertex((/-id(2,1),id(2,1)/),hel(1),useAHcoupl=testAcoupl)
      idV(2) = CoupledVertex((/id(1,1),-id(1,1)/),hel(1),useAHcoupl=testAcoupl)
   else
      idV(1) = CoupledVertex((/-id(1,1),id(1,1)/),hel(1),useAHcoupl=testAcoupl)
      idV(2) = CoupledVertex((/id(2,1),-id(2,1)/),hel(1),useAHcoupl=testAcoupl)
   endif

   if(idV(1).ne.Not_a_particle_) then
      ! Case 1: V->f1 fb1, f1->f1p W -> f2 fb2
      q2(:)=p(:,1,2)+p(:,2,2)
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),id(2,1)/),hel(1))
         q1scprop = ScalarPropagator(id(2,1),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/-id(2,1),id(2,1)/),hel(1),useA=useAcoupl) ! A->f1f1 V->f2f2
         else
            prefactor = CurrentPrefactor((/-id(2,1),id(2,1)/),hel(1)) ! Z->f1f1 V->f2f2
         endif
      else
         q1(:)=q(:)-p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),id(1,1)/),hel(1))
         q1scprop = ScalarPropagator(id(1,1),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1),useA=useAcoupl) ! A->f1f1 V->f2f2
         else
            prefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1)) ! Z->f1f1 V->f2f2
         endif
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_tmp,useAcoupl=.false.) ! idV_tmp=Wpm_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      ZAf_fWfpfp(:) = ZAf_fWfpfp(:) + restmp(:)
   endif
   if(idV(2).ne.Not_a_particle_) then
      ! Case 2: V->f1 fb1, fb1->fb1p W -> f2 fb2
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),id(1,1)/),hel(1))
         q1scprop = ScalarPropagator(id(1,1),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/-id(1,1),id(1,1)/),hel(1),useA=useAcoupl) ! A->f1f1 V->f2f2
         else
            prefactor = CurrentPrefactor((/-id(1,1),id(1,1)/),hel(1)) ! Z->f1f1 V->f2f2
         endif
      else
         q1(:)=-q(:)+p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),id(2,1)/),hel(1))
         q1scprop = ScalarPropagator(id(2,1),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1),useA=useAcoupl) ! A->f1f1 V->f2f2
         else
            prefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1)) ! Z->f1f1 V->f2f2
         endif
      endif
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      ZAf_fWfpfp(:) = ZAf_fWfpfp(:) + restmp(:)
   endif

   
   if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
      idV(1) = CoupledVertex((/-id(2,2),id(2,2)/),hel(2),useAHcoupl=testAcoupl)
      idV(2) = CoupledVertex((/id(1,2),-id(1,2)/),hel(2),useAHcoupl=testAcoupl)
   else
      idV(1) = CoupledVertex((/-id(1,2),id(1,2)/),hel(2),useAHcoupl=testAcoupl)
      idV(2) = CoupledVertex((/id(2,2),-id(2,2)/),hel(2),useAHcoupl=testAcoupl)
   endif
   if(idV(1).ne.Not_a_particle_) then
      ! Case 1: V->f2 fb2, f2->f2p W -> f1 fb1
      q2(:)=p(:,1,1)+p(:,2,1)
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=q(:)-p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),id(2,2)/),hel(2))
         q1scprop = ScalarPropagator(id(2,2),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/-id(2,2),id(2,2)/),hel(2),useA=useAcoupl) ! A->f2f2 W->f1f1
         else
            prefactor = CurrentPrefactor((/-id(2,2),id(2,2)/),hel(2)) ! Z->f2f2 W->f1f1
         endif
      else
         q1(:)=q(:)-p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),id(1,2)/),hel(2))
         q1scprop = ScalarPropagator(id(1,2),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2),useA=useAcoupl) ! A->f2f2 V->f1f1
         else
            prefactor = CurrentPrefactor((/id(1,2),-id(1,2)/),hel(2)) ! Z->f2f2 V->f1f1
         endif
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,1),id(1:2,1),hel(1),idV_tmp,useAcoupl=.false.) ! idV_tmp=Wpm_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(1),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      ZAf_fWfpfp(:) = ZAf_fWfpfp(:) + restmp(:)
   endif
   if(idV(2).ne.Not_a_particle_) then
      ! Case 2: V->f2 fb2, fb2->fbp2 W -> f1 fb1
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=-q(:)+p(:,1,2)
         innerprefactor = CurrentPrefactor((/id(2,2),id(1,2)/),hel(2))
         q1scprop = ScalarPropagator(id(1,2),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/-id(1,2),id(1,2)/),hel(2),useA=useAcoupl) ! A->f2f2 V->f1f1
         else
            prefactor = CurrentPrefactor((/-id(1,2),id(1,2)/),hel(2)) ! Z->f2f2 V->f1f1
         endif
      else
         q1(:)=-q(:)+p(:,2,2)
         innerprefactor = CurrentPrefactor((/id(1,2),id(2,2)/),hel(2))
         q1scprop = ScalarPropagator(id(2,2),q1)
         if(present(useAcoupl)) then
            prefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2),useA=useAcoupl) ! A->f2f2 V->f1f1
         else
            prefactor = CurrentPrefactor((/id(2,2),-id(2,2)/),hel(2)) ! Z->f2f2 V->f1f1
         endif
      endif
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))*q1scprop
      restmp = VectorPropagator(idV(2),q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      ZAf_fWfpfp(:) = ZAf_fWfpfp(:) + restmp(:)
   endif

   return
end function ZAf_fWfpfp


! Wff,f->fA/Zf'f'
function Wf_fZAfpfp(p,id,hel,Ub,V)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:2) ! Lorentz, 1,i=ubar(i); 2,i=v(i)
   integer, intent(in) :: id(1:2,1:2)
   integer, intent(in) :: hel(1:2)
   complex(dp), intent(in) :: Ub(1:4,1:2),V(1:4,1:2)
   integer :: idV ! V->f1 fb1 only
   integer :: idc(1:2,1:2)
   complex(dp) :: Wf_fZAfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: i,j
   real(dp) :: q(1:4),q1(1:4),q2(1:4)
   complex(dp) :: innercurrent_tmp(1:4),q1scprop,prefactor,innerprefactor
   integer :: idV_tmp

   if(hel(1).ne.-1 .or. hel(2).eq.0) return
   do i=1,2;do j=1,2
   idc(i,j)=convertLHE(id(i,j))
   enddo;enddo
   if( idc(1,1).eq.-idc(2,1) .or. idc(1,1).ne.-sign(idc(1,1),idc(2,1)) .or. idc(1,2).ne.-idc(2,2) ) return ! The presence of Z or A inside does not change the ids

   prefactor = CurrentPrefactor(id(:,1),hel(1)) ! W->f1f1 V->f2f2
   idV = CoupledVertex(id(:,1),hel(1))

   Wf_fZAfpfp(:) = czero
   q(:)=0_dp
   do i=1,2
      q(:)=q(:)+p(:,1,i)+p(:,2,i)
   enddo

   if(idV.ne.Not_a_particle_) then
      ! Case 1: W->f1 fb1, f1->f1 Z -> f2 fb2
      q2(:)=p(:,1,2)+p(:,2,2)
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1))
         q1scprop = ScalarPropagator(-id(1,1),q1)
      else
         q1(:)=q(:)-p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1))
         q1scprop = ScalarPropagator(-id(2,1),q1)
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_tmp,useAcoupl=.false.) ! idV_tmp=Z0_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV,q,restmp) ! Add -igmunu/[W] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      Wf_fZAfpfp(:) = Wf_fZAfpfp(:) + restmp(:)

      ! Case 2: W->f1 fb1, fb1->fb1 Z -> f2 fb2
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         innerprefactor = CurrentPrefactor((/-id(2,1),id(2,1)/),hel(1))
         q1scprop = ScalarPropagator(id(2,1),q1)
      else
         q1(:)=-q(:)+p(:,2,1)
         innerprefactor = CurrentPrefactor((/-id(1,1),id(1,1)/),hel(1))
         q1scprop = ScalarPropagator(id(1,1),q1)
      endif
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV,q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      Wf_fZAfpfp(:) = Wf_fZAfpfp(:) + restmp(:)

      ! Case 3: W->f1 fb1, f1->f1 A -> f2 fb2
      q2(:)=p(:,1,2)+p(:,2,2)
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         innerprefactor = CurrentPrefactor((/id(1,1),-id(1,1)/),hel(1),useA=.true.)
         q1scprop = ScalarPropagator(-id(1,1),q1)
      else
         q1(:)=q(:)-p(:,1,1)
         innerprefactor = CurrentPrefactor((/id(2,1),-id(2,1)/),hel(1),useA=.true.)
         q1scprop = ScalarPropagator(-id(2,1),q1)
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_tmp,useAcoupl=.true.) ! idV_tmp=Pho_
      innercurrent_tmp = VectorPropagator(idV_tmp,q2,innercurrent_tmp)
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV,q,restmp) ! Add -igmunu/[W] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      Wf_fZAfpfp(:) = Wf_fZAfpfp(:) + restmp(:)

      ! Case 4: W->f1 fb1, fb1->fb1 A -> f2 fb2
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         innerprefactor = CurrentPrefactor((/-id(2,1),id(2,1)/),hel(1),useA=.true.)
         q1scprop = ScalarPropagator(id(2,1),q1)
      else
         q1(:)=-q(:)+p(:,2,1)
         innerprefactor = CurrentPrefactor((/-id(1,1),id(1,1)/),hel(1),useA=.true.)
         q1scprop = ScalarPropagator(id(1,1),q1)
      endif
      restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))*q1scprop
      restmp = VectorPropagator(idV,q,restmp) ! Add -igmunu/[V] piece
      restmp(:) = restmp(:)*prefactor*innerprefactor
      Wf_fZAfpfp(:) = Wf_fZAfpfp(:) + restmp(:)
   endif

   return
end function Wf_fZAfpfp


! Wff,f->fpWf'f'
function Wf_fWfpfp(p,id,hel,Ub,V)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:2) ! Lorentz, 1,i=ubar(i); 2,i=v(i)
   integer, intent(in) :: id(1:2,1:2)
   integer, intent(in) :: hel(1:2)
   complex(dp), intent(in) :: Ub(1:4,1:2),V(1:4,1:2)
   integer :: idV(2) ! V->f1 fb1 only
   integer :: idc(1:2,1:2)
   complex(dp) :: Wf_fWfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: i,j
   real(dp) :: q(1:4),q1(1:4),q2(1:4)
   complex(dp) :: innercurrent_tmp(1:4),q1scprop,prefactor,innerprefactor,allprefactor
   integer :: idV_tmp, idV_dum
   logical :: isHadronicNonW(1:2)
   integer :: rootPair(1:3),rp

   isHadronicNonW(1:2)=.false.
   if(hel(1).ne.-1 .or. hel(2).ne.-1) return
   do i=1,2;do j=1,2
   idc(i,j)=convertLHE(id(i,j))
   enddo;enddo

   idV(1) = CoupledVertex(id(:,1),hel(1))
   idV(2) = CoupledVertex(id(:,2),hel(2))
   if( .not.( &
       ( ( (idV(1).eq.Wm_ .or. idV(2).eq.Z0_) .or. (idV(1).eq.Wp_ .or. idV(2).eq.Z0_) ) .or. ( (idV(2).eq.Wm_ .or. idV(1).eq.Z0_) .or. (idV(2).eq.Wp_ .or. idV(1).eq.Z0_) ) ) .or. &
       ( ( (idV(1).eq.Wm_ .or. idV(2).eq.Not_a_particle_) .or. (idV(1).eq.Wp_ .or. idV(2).eq.Not_a_particle_) ) .or. ( (idV(2).eq.Wm_ .or. idV(1).eq.Not_a_particle_) .or. (idV(2).eq.Wp_ .or. idV(1).eq.Not_a_particle_) ) ) &
            ) &
   ) return
   if(idV(1).eq.Not_a_particle_) then
      idV(1)=CoupledVertex_FlavorViolating(id(:,1),hel(1))
      if(idV(1).eq.Not_a_particle_) return
      isHadronicNonW(1)=.true.
   endif
   if(idV(2).eq.Not_a_particle_) then
      idV(2)=CoupledVertex_FlavorViolating(id(:,2),hel(2))
      if(idV(2).eq.Not_a_particle_) return
      isHadronicNonW(2)=.true.
   endif
   if(isHadronicNonW(1).and.isHadronicNonW(2)) return ! This means the flavor is somehow unidentified

   ! For this combination, there is only one possible diagram. The other one with the currents swapped is either a W--Z (leptonic or identical q1qb1) or not a valid f1-fb1 combination.
   Wf_fWfpfp(:) = czero
   q(:)=0_dp
   do i=1,2
      q(:)=q(:)+p(:,1,i)+p(:,2,i)
   enddo

   if(idV(1).eq.Z0_) then ! Remember, the flavor-violating combination is also a Z0_ now.
      ! Case 1: W->u1 db1, u1->d1 W -> f2 fb2
      q2(:)=p(:,1,2)+p(:,2,2)

      allprefactor = czero
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=q(:)-p(:,2,1)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(2,1))
         do rp=1,3
            idV_tmp = CoupledVertex((/ rootPair(rp), id(2,1) /),hel(1))
            if(idV_tmp .eq. idV(2)) then
               prefactor = CurrentPrefactor((/ rootPair(rp), id(2,1) /),hel(1)) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/id(1,1),-rootPair(rp)/),hel(1))
               q1scprop = ScalarPropagator(-rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      else
         q1(:)=q(:)-p(:,1,1)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,1))
         do rp=1,3
            idV_tmp = CoupledVertex((/ rootPair(rp), id(1,1) /),hel(1))
            if(idV_tmp .eq. idV(2)) then
               prefactor = CurrentPrefactor((/ rootPair(rp), id(1,1) /),hel(1)) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/id(2,1),-rootPair(rp)/),hel(1))
               q1scprop = ScalarPropagator(-rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,2),id(1:2,2),hel(2),idV_dum) ! idV_dum=Wpm_, includes W->f2 fb2
      if(idV_dum .eq. idV(2)) then
         innercurrent_tmp = VectorPropagator(idV_dum,q2,innercurrent_tmp)
         restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))
         restmp = VectorPropagator(idV_tmp,q,restmp) ! Add -igmunu/[W] piece
         restmp(:) = restmp(:)*allprefactor
         Wf_fWfpfp(:) = Wf_fWfpfp(:) + restmp(:)
      endif

      ! Case 2: W->f1 fb1, fb1->fb1 W -> f2 fb2
      allprefactor = czero
      if(idc(1,1).gt.0 .and. idc(2,1).lt.0) then
         q1(:)=-q(:)+p(:,1,1)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,1))
         do rp=1,3
            idV_tmp = CoupledVertex((/ id(1,1), rootPair(rp) /),hel(1))
            if(idV_tmp .eq. idV(2)) then
               prefactor = CurrentPrefactor((/ id(1,1), rootPair(rp) /),hel(1)) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/ -rootPair(rp), id(2,1) /),hel(1))
               q1scprop = ScalarPropagator(rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      else
         q1(:)=-q(:)+p(:,2,1)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,1))
         do rp=1,3
            idV_tmp = CoupledVertex((/ id(2,1), rootPair(rp) /),hel(1))
            if(idV_tmp .eq. idV(2)) then
               prefactor = CurrentPrefactor((/ id(2,1), rootPair(rp) /),hel(1)) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/ -rootPair(rp), id(1,1) /),hel(1))
               q1scprop = ScalarPropagator(rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      if(idV_dum .eq. idV(2)) then
         restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,1),V(:,1))
         restmp = VectorPropagator(idV_tmp,q,restmp) ! Add -igmunu/[W] piece
         restmp(:) = restmp(:)*allprefactor
         Wf_fWfpfp(:) = Wf_fWfpfp(:) + restmp(:)
      endif
   endif

   if(idV(2).eq.Z0_) then ! Remember, the flavor-violating combination is also a Z0_ now.
      ! Case 1: W->f2 fb2, f2->fp2 W -> f1 fb1
      q2(:)=p(:,1,1)+p(:,2,1)

      allprefactor = czero
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=q(:)-p(:,2,2)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(2,2))
         do rp=1,3
            idV_tmp = CoupledVertex((/ rootPair(rp), id(2,2) /),hel(2))
            if(idV_tmp .eq. idV(1)) then
               prefactor = CurrentPrefactor((/ rootPair(rp), id(2,2) /),hel(2)) ! W->f2 fb2
               innerprefactor = CurrentPrefactor((/id(1,2),-rootPair(rp)/),hel(2)) ! W->f2 fp2 
               q1scprop = ScalarPropagator(-rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      else
         q1(:)=q(:)-p(:,1,2)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,2))
         do rp=1,3
            idV_tmp = CoupledVertex((/ rootPair(rp), id(1,2) /),hel(2))
            if(idV_tmp .eq. idV(1)) then
               prefactor = CurrentPrefactor((/ rootPair(rp), id(1,2) /),hel(2)) ! W->f2 fb2
               innerprefactor = CurrentPrefactor((/id(2,2),-rootPair(rp)/),hel(2)) ! W->f2 fp2 
               q1scprop = ScalarPropagator(-rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      innercurrent_tmp = Vcurrent(p(1:4,1:2,1),id(1:2,1),hel(1),idV_dum) ! idV_dum=Wpm_, includes W->f1 fb1
      if(idV_dum .eq. idV(1)) then
         innercurrent_tmp = VectorPropagator(idV_dum,q2,innercurrent_tmp)
         restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))
         restmp = VectorPropagator(idV_tmp,q,restmp) ! Add -igmunu/[W] piece
         restmp(:) = restmp(:)*allprefactor
         Wf_fWfpfp(:) = Wf_fWfpfp(:) + restmp(:)
      endif

      ! Case 2: W->f2 fb2, fb2->fbp2 W -> f1 fb1
      allprefactor = czero
      if(idc(1,2).gt.0 .and. idc(2,2).lt.0) then
         q1(:)=-q(:)+p(:,1,2)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,2))
         do rp=1,3
            idV_tmp = CoupledVertex((/ id(1,2), rootPair(rp) /),hel(2))
            if(idV_tmp .eq. idV(1)) then
               prefactor = CurrentPrefactor((/ id(1,2), rootPair(rp) /),hel(2)) ! W->f2 fb2
               innerprefactor = CurrentPrefactor((/ -rootPair(rp), id(2,2) /),hel(2)) ! W->fb2 fbp2 
               q1scprop = ScalarPropagator(rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      else
         q1(:)=-q(:)+p(:,2,2)
         rootPair(:)=Not_a_particle_
         rootPair = WDaughterPair(id(1,2))
         do rp=1,3
            idV_tmp = CoupledVertex((/ id(2,2), rootPair(rp) /),hel(2))
            if(idV_tmp .eq. idV(1)) then
               prefactor = CurrentPrefactor((/ id(2,2), rootPair(rp) /),hel(2)) ! W->f2 fb2
               innerprefactor = CurrentPrefactor((/ -rootPair(rp), id(1,2) /),hel(2)) ! W->fb2 fbp2 
               q1scprop = ScalarPropagator(rootPair(rp),q1)
               allprefactor = allprefactor + prefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      if(idV_dum .eq. idV(1)) then
         restmp = gamma_gamma_mu(3,(/ cmplx(q1(:),kind=dp), innercurrent_tmp(1:4) /),Ub(:,2),V(:,2))
         restmp = VectorPropagator(idV_tmp,q,restmp) ! Add -igmunu/[W] piece
         restmp(:) = restmp(:)*allprefactor
         Wf_fWfpfp(:) = Wf_fWfpfp(:) + restmp(:)
      endif
   endif

   return
end function Wf_fWfpfp



! 4H
function QuarticHVertex(idV)
   implicit none
   integer, intent(in) :: idV(1:4) ! current ids
   integer :: order(1:4)
   complex(dp) :: QuarticHVertex

   QuarticHVertex = czero
   order(:)=0
   order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Hig_,Hig_ /))
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
      print *,"QuarticHVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   QuarticHVertex = -ci*3_dp*M_Reso**2/vev
   return
end function QuarticHVertex

! 3H
function TripleHVertex(idV)
   implicit none
   integer, intent(in) :: idV(1:3) ! current ids
   integer :: order(1:3)
   complex(dp) :: TripleHVertex

   TripleHVertex = czero
   order(:)=0
   order(:)=Id_Order(3,idV,(/ Hig_,Hig_,Hig_ /))
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
      print *,"QuarticHVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   TripleHVertex = -ci*3_dp*(M_Reso/vev)**2
   return
end function TripleHVertex

! VVHH
function DoubleVHVertex(p,current,idV,outType)
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4),current(1:4,1:4) ! Lorentz, current index
   integer, intent(in) :: idV(1:4) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:4)
   complex(dp) :: DoubleVHVertex

   DoubleVHVertex = czero
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Wp_,Wm_ /))
      M_V=M_W
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Z0_,Z0_ /))
      M_V=M_Z
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
      print *,"DoubleVHVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   DoubleVHVertex = ci*2_dp*(M_V/vev)**2*( current(:,order(3)).dot.current(:,order(4)) )
   return
end function DoubleVHVertex

! VVH
function VVHVertex(p,current,idV,outType)
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:3),current(1:4,1:3) ! Lorentz, current index
   integer, intent(in) :: idV(1:3) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:3)
   complex(dp) :: VVHVertex
   complex(dp) :: aa(1:3)
   complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
   complex(dp) :: e3_e4,q3_q4,e3_q4,e4_q3
   real(dp) :: q_q,q3_q3,q4_q4

   VVHVertex = czero
   aa(:)=czero
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Wp_,Wm_ /))
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Z0_,Z0_ /))
   elseif( (outType(1).eq.Z0_ .and. outType(2).eq.Pho_) .or. (outType(1).eq.Pho_ .and. outType(2).eq.Z0_) ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Z0_,Pho_ /))
   elseif( outType(1).eq.Pho_ .and. outType(2).eq.Pho_ ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Pho_,Pho_ /))
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
      print *,"VVHVertex::Target mismatch (ids, order)",idV,order
      return
   endif

   q_q   = p(:,order(1)).dot.p(:,order(1))
   q3_q3 = p(:,order(2)).dot.p(:,order(2))
   q4_q4 = p(:,order(3)).dot.p(:,order(3))
   q3_q4 = p(:,order(2)).dot.p(:,order(3))
   e3_e4 = current(:,order(2)).dot.current(:,order(3))
   e3_q4 = current(:,order(2)).dot.p(:,order(3))
   e4_q3 = p(:,order(2)).dot.current(:,order(3))

   if( outType(1).ne.Pho_ .and. outType(2).ne.Pho_ ) then
      ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q_q)
      ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q_q)
      if(outType(1).eq.Z0_) then
         M_V=M_Z
      else
         M_V=M_W
      endif
      aa(1) = ghz1_dyn*M_V**2                             &
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   elseif( outType(1).ne.Pho_ .or. outType(2).ne.Pho_ ) then ! Parameterization goes as Zgs in ModHiggs
      ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q_q)
      ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q_q)

      aa(1) = ghz1_dyn*q_q                                & ! Notice mG**2 instead
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   else
      ghz1_dyn = czero
      ghz2_dyn = HVVSpinZeroDynamicCoupling(9, 0d0,0d0,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(10,0d0,0d0,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(11,0d0,0d0,q_q)

      aa(1) =                                             &
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   endif

   aa(1) = aa(1)*e3_e4
   aa(2) = aa(2)*e3_q4*e4_q3
   aa(3) = aa(3)*epsilon_munuab(current(:,order(2)),current(:,order(3)),p(:,order(2)),p(:,order(3)))

   VVHVertex = ci*( aa(1) + aa(2) + aa(3) )/vev
   return
end function VVHVertex


function ScalarPropagator(idV,p,cmasssq)
   use ModMisc
   implicit none
   integer, intent(in) :: idV
   real(dp), intent(in) :: p(1:4)
   complex(dp), optional :: cmasssq
   complex(dp) :: ScalarPropagator
   real(dp) :: s

   ScalarPropagator = czero
   s = p(:).dot.p(:)
   if(idV .eq. Wm_ .or. idV.eq.Wp_) then
      ScalarPropagator = -ci/( s - M_W**2 + ci*M_W*Ga_W)
   elseif(idV .eq. Z0_) then
      ScalarPropagator = -ci/( s - M_Z**2 + ci*M_Z*Ga_Z)
   elseif(idV .eq. Pho_) then
      ScalarPropagator = -ci/s
   elseif(idV .eq. Hig_) then
      ScalarPropagator = ci/( s - M_Reso**2 + ci*M_Reso*Ga_Reso)
   elseif(present(cmasssq)) then
      ScalarPropagator = ci/( s - cmasssq)
   else
      ScalarPropagator = ci/s
   endif

   return
end function ScalarPropagator


function VectorPropagator(idV,p,current)
   use ModMisc
   implicit none
   integer, intent(in) :: idV
   real(dp), intent(in) :: p(1:4)
   complex(dp), intent(in) :: current(1:4)
   complex(dp) :: VectorPropagator(1:4)
   complex(dp) :: prefactor

   VectorPropagator(:) = czero
   prefactor = ScalarPropagator(idV,p)

   if(idV .eq. Wm_ .or. idV.eq.Wp_ .or. idV .eq. Z0_) then
      VectorPropagator(:) = prefactor*(current(:)-p(:)*((cmplx(p(:),kind=dp)).dot.current(:))/(p(:).dot.p(:)))
   elseif(idV .eq. Pho_) then
      VectorPropagator(:) = prefactor*current(:)
   endif

   return
end function VectorPropagator

function EuclideanMagnitude_Complex(ND,p)
   implicit none
   integer, intent(in) :: ND
   complex(dp), intent(in) :: p(1:ND)
   real(dp) :: EuclideanMagnitude_Complex
   integer :: i

   EuclideanMagnitude_Complex=0_dp
   if(ND.gt.0) then
      do i=1,ND
         EuclideanMagnitude_Complex = EuclideanMagnitude_Complex + p(i)*dconjg(p(i))
      enddo
      EuclideanMagnitude_Complex = dsqrt(EuclideanMagnitude_Complex)
   endif
   return
end function EuclideanMagnitude_Complex


function getCurrents_VA(p,id,hel,tryA,isWBF,current,qV,idV)
   implicit none
   real(dp), intent(in) :: p(1:4,1:8) ! 17,28 - 34,56
   integer, intent(in) :: id(1:8)
   integer, intent(in) :: hel(1:4)
   integer, intent(in) :: tryA(1:2) ! First and second A current index
   logical, intent(in) :: isWBF
   complex(dp), intent(out) :: current(1:4,1:4) ! Lorentz, current
   real(dp), intent(out) :: qV(1:4,1:4)
   integer, intent(out) :: idV(1:4)
   real(dp) :: pout(1:4,1:8),eucmagtmp
   integer :: idout(1:8)
   logical :: getCurrents_VA
   integer :: i
   integer :: order(1:2,1:4)

   getCurrents_VA=.true.
   eucmagtmp=0_dp

   pout(:,3:8)=p(:,3:8)
   pout(:,1)=-p(:,1)
   pout(:,2)=-p(:,2)

   idout(3:8)=id(3:8)
   idout(1)=-id(1)
   idout(2)=-id(2)

   current(:,:)=czero
   if(isWBF) then
      order(:,1) = (/ 1, 7 /)
      order(:,2) = (/ 2, 8 /)
   else
      order(:,1) = (/ 1, 2 /)
      order(:,2) = (/ 7, 8 /)
   endif
   order(:,3) = (/ 3, 4 /)
   order(:,4) = (/ 5, 6 /)

   do i=1,4
      if(i.eq.tryA(1) .or. i.eq.tryA(2)) then
         current(:,i)=Vcurrent((/pout(:,order(1,i)),pout(:,order(2,i))/),(/idout(order(1,i)),idout(order(2,i))/),hel(i),idV(i),useAcoupl=.true.)
      else
         current(:,i)=Vcurrent((/pout(:,order(1,i)),pout(:,order(2,i))/),(/idout(order(1,i)),idout(order(2,i))/),hel(i),idV(i))
      endif
      qV(:,i) = pout(:,order(1,i))+pout(:,order(2,i))
      if( (idV(i).eq.Not_a_particle_).or.(EuclideanMagnitude_Complex(4,current(:,i)).eq.0_dp) ) then
         getCurrents_VA=.false.
         exit
      else
         current(:,i)=VectorPropagator(idV(i),qV(:,i),current(:,i))
      endif
   enddo
   return
end function getCurrents_VA


subroutine amp_WWZZ(current,idV,ime)
   implicit none
   complex(dp), intent(in) :: current(1:4,1:4) ! Lorentz, current
   integer, intent(in) :: idV(1:4)
   complex(dp), intent(out) :: ime
   integer :: outType(1:2)

   ime = czero

   outType(:) = (/ Z0_, Z0_ /)
   ime = QuarticEWKVertex(current,idV,outType)
   return
end subroutine amp_WWZZ

subroutine amp_WWAA_4to4(current,idV,ime)
   implicit none
   complex(dp), intent(in) :: current(1:4,1:4) ! Lorentz, current
   integer, intent(in) :: idV(1:4)
   complex(dp), intent(out) :: ime
   integer :: outType(1:2)

   ime = czero

   outType(:) = (/ Pho_, Pho_ /)
   ime = QuarticEWKVertex(current,idV,outType)
   return
end subroutine amp_WWAA_4to4

subroutine amp_WWZA_4to4(current,idV,ime)
   implicit none
   complex(dp), intent(in) :: current(1:4,1:4) ! Lorentz, current
   integer, intent(in) :: idV(1:4)
   complex(dp), intent(out) :: ime
   integer :: outType(1:2)

   ime = czero

   outType(:) = (/ Z0_, Pho_ /)
   ime = QuarticEWKVertex(current,idV,outType)
   return
end subroutine amp_WWZA_4to4

subroutine amp_WWWW(current,qV,idV,ime)
   implicit none
   complex(dp), intent(in) :: current(1:4,1:4) ! Lorentz, current
   real(dp), intent(in) :: qV(1:4,1:4)
   integer, intent(in) :: idV(1:4)
   complex(dp), intent(out) :: ime
   complex(dp) :: ime_WW_Z_WW,ime_WW_A_WW
   integer :: outType(1:2)
   real(dp) :: qZA(1:4,1:2)
   complex(dp) :: cZA(1:4,1:2)
   integer :: order(1:4),h,mu

   qZA(:,:) = 0_dp
   cZA(:,:) = czero
   ime = czero
   ime_WW_Z_WW = czero
   ime_WW_A_WW = czero

   outType(:) = (/ Wp_, Wm_ /)
   order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))

   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      ime = QuarticEWKVertex(current,idV,outType)

      qZA(:,1) = -qV(:,order(1)) - qV(:,order(2))
      qZA(:,2) = -qV(:,order(3)) - qV(:,order(4))
      do h=-1,1
         cZA(:,1) = pol_massive(qZA(:,1),h)
         do mu=1,4
            cZA(mu,2) = -conjg(cZA(mu,1)) ! cZA(:,2) = -pol_massive(qZA(:,2),h)
         enddo

         ime_WW_Z_WW = ime_WW_Z_WW + &
                       TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(2)) /), (/ cZA(:,1), current(:,order(1)), current(:,order(2)) /), (/ Z0_, idV(order(1)), idV(order(2)) /) ) * &
                       TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(4)) /), (/ cZA(:,2), current(:,order(3)), current(:,order(4)) /), (/ Z0_, idV(order(3)), idV(order(4)) /) )

         if(h.ne.0) then
            cZA(:,1) = pol_massless(qZA(:,1),h)
            do mu=1,4
               cZA(mu,2) = -conjg(cZA(mu,1))
            enddo
            ime_WW_A_WW = ime_WW_A_WW + &
                          TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(2)) /), (/ cZA(:,1), current(:,order(1)), current(:,order(2)) /), (/ Pho_, idV(order(1)), idV(order(2)) /), useAcoupl=.true. ) * &
                          TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(4)) /), (/ cZA(:,2), current(:,order(3)), current(:,order(4)) /), (/ Pho_, idV(order(3)), idV(order(4)) /), useAcoupl=.true. )
         endif
      enddo
      ime_WW_Z_WW = ime_WW_Z_WW * ScalarPropagator(Z0_,qZA(:,1))
      ime_WW_A_WW = ime_WW_A_WW * ScalarPropagator(Pho_,qZA(:,1))

      qZA(:,1) = -qV(:,order(1)) - qV(:,order(4))
      qZA(:,2) = -qV(:,order(3)) - qV(:,order(2))
      do h=-1,1
         cZA(:,1) = pol_massive(qZA(:,1),h)
         do mu=1,4
            cZA(mu,2) = -conjg(cZA(mu,1))
         enddo

         ime_WW_Z_WW = ime_WW_Z_WW + &
                       TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(4)) /), (/ cZA(:,1), current(:,order(1)), current(:,order(4)) /), (/ Z0_, idV(order(1)), idV(order(4)) /) ) * &
                       TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(2)) /), (/ cZA(:,2), current(:,order(3)), current(:,order(2)) /), (/ Z0_, idV(order(3)), idV(order(2)) /) )

         if(h.ne.0) then
            cZA(:,1) = pol_massless(qZA(:,1),h)
            do mu=1,4
               cZA(mu,2) = -conjg(cZA(mu,1))
            enddo
            ime_WW_A_WW = ime_WW_A_WW + &
                          TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(4)) /), (/ cZA(:,1), current(:,order(1)), current(:,order(4)) /), (/ Pho_, idV(order(1)), idV(order(4)) /), useAcoupl=.true. ) * &
                          TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(2)) /), (/ cZA(:,2), current(:,order(3)), current(:,order(2)) /), (/ Pho_, idV(order(3)), idV(order(2)) /), useAcoupl=.true. )
         endif
      enddo
      ime_WW_Z_WW = ime_WW_Z_WW * ScalarPropagator(Z0_,qZA(:,1))
      ime_WW_A_WW = ime_WW_A_WW * ScalarPropagator(Pho_,qZA(:,1))

      ime = ime + ime_WW_Z_WW + ime_WW_A_WW
   endif
   return
end subroutine amp_WWWW

subroutine amp_WW_V_4f(p,id,hel,Ub,V,current,ime)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:4) ! Lorentz, particle 1/2, root vector boson
   integer, intent(in) :: id(1:2,1:4)
   integer, intent(in) :: hel(1:4)
   complex(dp), intent(in) :: Ub(1:4,1:4),V(1:4,1:4)
   complex(dp), intent(in) :: current(1:4,1:4) ! Pass the currents as well since they are already computed
   complex(dp), intent(out) :: ime
   integer :: idc(1:2,1:4),order(1:4)
   complex(dp) :: composite(1:4,1:2)
   integer :: idV(1:4),i
   real(dp) :: pordered(1:4,1:2,1:2), qV(1:4,1:4),qZA(1:4)
   integer :: idordered(1:2,1:2),mu

   ime = czero
   composite(:,:)=czero

   do i=1,4
      idV(i) = CoupledVertex(id(1:2,i),hel(i))
      if(idV(i).eq.Not_a_particle_) return
      do mu=1,4
         qV(mu,i) = p(mu,1,i)+p(mu,2,i)
      enddo
      idc(1,i)=convertLHE(id(1,i))
      idc(2,i)=convertLHE(id(2,i))
   enddo

   order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Z0_ /)) ! Z0_,Z0_ is dummy here, just to ensure identical particles are present. There is no need to call it a second time for Pho_,Pho_; it is checked inside ZAf_fZAfpfp in a more correct way.
   if(order(1).eq.0 .or.order(2).eq.0 .or.order(3).eq.0 .or.order(4).eq.0) return

   do i=3,4
      if(idc(1,order(i)).gt.idc(2,order(i))) then
         idordered(1,i-2)=id(1,order(i))
         idordered(2,i-2)=id(2,order(i))
         do mu=1,4
            pordered(mu,1,i-2)=p(mu,1,order(i))
            pordered(mu,2,i-2)=p(mu,2,order(i))
         enddo
      else
         idordered(1,i-2)=id(2,order(i))
         idordered(2,i-2)=id(1,order(i))
         do mu=1,4
            pordered(mu,1,i-2)=p(mu,2,order(i))
            pordered(mu,2,i-2)=p(mu,1,order(i))
         enddo
      endif
   enddo
   qZA(:) = -qV(:,order(1))-qV(:,order(2))

   composite(:,1) = ZAf_fZAfpfp( pordered, idordered, (/hel(order(3)),hel(order(4))/), (/Ub(1:4,order(3)),Ub(1:4,order(4))/), (/V(1:4,order(3)),V(1:4,order(4))/) ) ! Z
   composite(:,2) = ZAf_fZAfpfp( pordered, idordered, (/hel(order(3)),hel(order(4))/), (/Ub(1:4,order(3)),Ub(1:4,order(4))/), (/V(1:4,order(3)),V(1:4,order(4))/), useAcoupl=.true. ) ! A
   ime = &
         TripleEWKVertex( (/ qZA(:), qV(:,order(1)), qV(:,order(2)) /), (/ composite(:,1), current(:,order(1)), current(:,order(2)) /), (/ Z0_ , idV(order(1)), idV(order(2)) /) ) &
       + TripleEWKVertex( (/ qZA(:), qV(:,order(1)), qV(:,order(2)) /), (/ composite(:,2), current(:,order(1)), current(:,order(2)) /), (/ Pho_, idV(order(1)), idV(order(2)) /) )

   return
end subroutine amp_WW_V_4f

subroutine amp_tchannelV_VffVfpfp(p,id,hel,Ub,V,current,ime)
   use ModMisc
   implicit none
   real(dp), intent(in) :: p(1:4,1:2,1:4) ! Lorentz, particle 1/2, root vector boson
   integer, intent(in) :: id(1:2,1:4)
   integer, intent(in) :: hel(1:4)
   complex(dp), intent(in) :: Ub(1:4,1:4),V(1:4,1:4)
   complex(dp), intent(in) :: current(1:4,1:4) ! Pass the currents as well since they are already computed
   complex(dp), intent(out) :: ime
   complex(dp) :: composite(1:4,1:2)
   logical :: isFlavorChangingW(1:4)
   integer :: idV(1:4),i,order(1:4)
   real(dp) :: pordered(1:4,1:2,1:4), qV(1:4,1:4), qtchannel(1:4)
   integer :: idordered(1:2,1:4),mu,lineorder(1:4)

   ime = czero
   composite(:,:)=czero
   isFlavorChangingW(:) = .false.

   do i=1,4
      idV(i) = CoupledVertex(id(1:2,i),hel(i))
      if(idV(i).eq.Not_a_particle_) idV(i) = CoupledVertex_FlavorViolating(id(1:2,i),hel(i)) ! Will return Z0_ if has valid pairing
      if(idV(i).eq.Not_a_particle_) then
         return
      else
         isFlavorChangingW(i)=.true.
      endif
      do mu=1,4
         qV(mu,i) = p(mu,1,i)+p(mu,2,i)
      enddo
   enddo

   ! Case 1: 12.Z.34.Z.56.Z.78
   order(:)=Id_Order(4,idV,(/ Z0_,Z0_,Z0_,Z0_ /)) ! Z0_,Z0_ is dummy here, just to ensure identical particles are present. There is no need to call it a second time for Pho_,Pho_; it is checked inside ZAf_fZAfpfp in a more correct way.
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then

      lineorder(:)=(/1,2,3,4/) ! 1+2 -- 3+4
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the Z propagator
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator

      lineorder(:)=(/1,3,2,4/) ! 1+3 -- 2+4
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the Z propagator
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator

      lineorder(:)=(/1,4,2,3/) ! 1+4 -- 2+3
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the Z propagator
      composite(:,1) = ZAf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator
   endif

   ! Case 2: 12.Z.34.W.56.Z.78 / 12.W.34.W.56.W.78 / 12.W.34.Z.56.Z.78 : Second one looks like a Z in the root!
   order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Z0_ /))
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      lineorder(:)=(/1,3,2,4/) ! W1Z3 -- W2Z4
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = Wf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      if(isFlavorChangingW(order(lineorder(2)))) then
         composite(:,1) = composite(:,1) &
                        + Wf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 2 and 1
      endif
      composite(:,2) = Wf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      if(isFlavorChangingW(order(lineorder(4)))) then
         composite(:,2) = composite(:,2) &
                        + Wf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      endif
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Wm_,qtchannel) ! Avoid double-counting the W propagator

      lineorder(:)=(/1,4,2,3/) ! W1Z4 -- W2Z3
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = Wf_fZAfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      if(isFlavorChangingW(order(lineorder(2)))) then
         composite(:,1) = composite(:,1) &
                        + Wf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 2 and 1
      endif
      composite(:,2) = Wf_fZAfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      if(isFlavorChangingW(order(lineorder(4)))) then
         composite(:,2) = composite(:,2) &
                        + Wf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      endif
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Wm_,qtchannel) ! Avoid double-counting the W propagator

      lineorder(:)=(/1,2,3,4/) ! W1W2 -- Z3Z4
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the W propagator
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator
   endif


   ! Case 3: 12.W+.34.Z.56.W-.78
   order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      lineorder(:)=(/1,2,3,4/) ! Wp1Wm1 -- Wp2Wm2
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the Z propagator
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator

      lineorder(:)=(/1,4,3,2/) ! Wp1Wm2 -- Wp2Wm1
      do i=1,4
         if(id(1,order(lineorder(i))).gt.id(2,order(lineorder(i)))) then
            idordered(1,i)=id(1,order(lineorder(i)))
            idordered(2,i)=id(2,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,1,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,2,order(lineorder(i)))
            enddo
         else
            idordered(1,i)=id(2,order(lineorder(i)))
            idordered(2,i)=id(1,order(lineorder(i)))
            do mu=1,4
               pordered(mu,1,i)=p(mu,2,order(lineorder(i)))
               pordered(mu,2,i)=p(mu,1,order(lineorder(i)))
            enddo
         endif
      enddo
      qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/) ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/) ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Z0_,qtchannel) ! Avoid double-counting the Z propagator
      composite(:,1) = ZAf_fWfpfp( pordered(:,:,1:2), idordered(:,1:2), (/hel(order(lineorder(1))),hel(order(lineorder(2)))/), (/Ub(1:4,order(lineorder(1))),Ub(1:4,order(lineorder(2)))/), (/V(1:4,order(lineorder(1))),V(1:4,order(lineorder(2)))/), useAcoupl=.true. ) ! Lines 1 and 2
      composite(:,2) = ZAf_fWfpfp( pordered(:,:,3:4), idordered(:,3:4), (/hel(order(lineorder(3))),hel(order(lineorder(4)))/), (/Ub(1:4,order(lineorder(3))),Ub(1:4,order(lineorder(4)))/), (/V(1:4,order(lineorder(3))),V(1:4,order(lineorder(4)))/), useAcoupl=.true. ) ! Lines 3 and 4
      ime = ime + (composite(:,1).dot.composite(:,2))/ScalarPropagator(Pho_,qtchannel) ! Avoid double-counting the A propagator
   endif

   return
end subroutine amp_tchannelV_VffVfpfp




end module ModVVHOffshell
