! Candidate continuous spectrum initializer
! No hypergeometric functions, unscaled Legendre values, or complex sine.
module sphere_riccati_init_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter :: dp=real64
  real(dp), parameter :: pi=acos(-1.0_dp)
  public :: riccati_mode_covariances, riccati_mode_filter, calibrate_spherical_c2, spherical_mode_cutoff
contains
  ! State is [g_m, (1-z*z)*dg_m/dz/kappa], kappa=sqrt(m*m+sqrt(c0/c2)).
  ! Return C(:,:,i)=Cov(state_i,state_i), X(:,:,i)=Cov(state_(i+1),state_i).
  ! rc: 1 invalid inputs; 2 nonfinite integration/covariance; 3 non-SPD
  ! covariance/innovation; 4 calibration did not converge within limits.
  subroutine riccati_mode_covariances(z,m,c0,c2,cov,cross,rc,refinement)
    real(dp), intent(in) :: z(:),c0,c2
    integer, intent(in) :: m
    real(dp), intent(out) :: cov(2,2,size(z)),cross(2,2,size(z)-1)
    integer, intent(out) :: rc
    integer, optional, intent(in) :: refinement
    real(dp) :: tt(size(z)), ratio, kap, fac, cur_t, target, hlim, delta
    complex(dp) :: bp(size(z)),bm(size(z)),jp(size(z)),d,jint,aa,bb,s,ss
    complex(dp) :: pole, ee
    integer :: i,j,k,n,ref
    cov=0; cross=0; rc=1; n=size(z); ref=1
    if(present(refinement)) ref=refinement
    if(n<2.or.m<0.or.c0<=0.or.c2<=0.or.ref<1) return
    if(.not.all(ieee_is_finite(z)).or..not.ieee_is_finite(c0).or..not.ieee_is_finite(c2)) return
    if(any(abs(z)>=1).or.any(z(2:)<=z(:n-1))) return
    ratio=sqrt(c0/c2); pole=cmplx(0.0_dp,ratio,dp)
    kap=sqrt(real(m,dp)**2+ratio); fac=1/(2*pi*c2*ratio)
    tt=atanh(z)
    cur_t=max(20.0_dp,tt(n)+2,-tt(1)+2)
    ! D=B+m. Unlike m*m-B*B, D*(2*m-D) retains small corrections at large m.
    d=pole/(2*(real(m,dp)+1)*cosh(cur_t)**2)
    jint=cmplx(0.0_dp,0.0_dp,dp)
    ! RK4 step cap: accuracy for sech(t)^2 and stability of the backward ODE.
    hlim=min(0.00125_dp,0.125_dp/max(1.0_dp,kap))/real(ref,dp)
    ! Merge the two descending lists t(n:1:-1) and -t(1:n).
    i=n; j=1
    do while(i>=1.or.j<=n)
      if(i<1) then
        target=-tt(j)
      else if(j>n) then
        target=tt(i)
      else
        target=max(tt(i),-tt(j))
      end if
      call integrate_to(target,cur_t,d,jint,m,pole,hlim)
      if(i>=1) then
        if(abs(tt(i)-target)<64*epsilon(1.0_dp)*max(1.0_dp,abs(target))) then
          bp(i)=-real(m,dp)+d; jp(i)=jint; i=i-1
        end if
      end if
      if(j<=n) then
        if(abs(-tt(j)-target)<64*epsilon(1.0_dp)*max(1.0_dp,abs(target))) then
          bm(j)=-real(m,dp)+d; j=j+1
        end if
      end if
    end do
    rc=2
    if(.not.all(ieee_is_finite(real(bp))).or..not.all(ieee_is_finite(aimag(bp)))) return
    if(.not.all(ieee_is_finite(real(bm))).or..not.all(ieee_is_finite(aimag(bm)))) return
    do k=1,n
      aa=-bm(k); bb=bp(k); s=1/(aa-bb)
      cov(1,1,k)=fac*aimag(s)
      cov(1,2,k)=fac*aimag(bb*s)/kap
      cov(2,1,k)=cov(1,2,k)
      cov(2,2,k)=fac*aimag(aa*bb*s)/kap**2
      if(k<n) then
        ! Logarithmic ratio: v(t_next)/v(t), never an enormous v(t).
        delta=tt(k+1)-tt(k)
        ee=cmplx(-real(m,dp)*delta,0.0_dp,dp)+jp(k+1)-jp(k)
        if(real(ee)<log(tiny(1.0_dp))+2.0_dp) then
          ss=cmplx(0.0_dp,0.0_dp,dp) ! a negligible cross-ring covariance
        else
          ss=exp(ee)*s
        end if
        cross(1,1,k)=fac*aimag(ss)
        cross(1,2,k)=fac*aimag(aa*ss)/kap
        cross(2,1,k)=fac*aimag(bp(k+1)*ss)/kap
        cross(2,2,k)=fac*aimag(aa*bp(k+1)*ss)/kap**2
      end if
    end do
    if(.not.all(ieee_is_finite(cov)).or..not.all(ieee_is_finite(cross))) return
    rc=0
  end subroutine

  subroutine integrate_to(target,t,d,jint,m,pole,hlim)
    real(dp), intent(in) :: target,hlim
    real(dp), intent(inout) :: t
    complex(dp), intent(inout) :: d,jint
    complex(dp), intent(in) :: pole
    integer, intent(in) :: m
    complex(dp) :: k1,k2,k3,k4,d2,d3,d4
    real(dp) :: h
    do while(t-target>8*epsilon(t)*max(1.0_dp,abs(t)))
      h=-min(hlim,t-target)
      k1=rhs(t,d,m,pole)
      d2=d+0.5_dp*h*k1; k2=rhs(t+0.5_dp*h,d2,m,pole)
      d3=d+0.5_dp*h*k2; k3=rhs(t+0.5_dp*h,d3,m,pole)
      d4=d+h*k3; k4=rhs(t+h,d4,m,pole)
      jint=jint+h*(d+2*d2+2*d3+d4)/6
      d=d+h*(k1+2*k2+2*k3+k4)/6
      t=t+h
    end do
    t=target
  end subroutine

  pure complex(dp) function rhs(t,d,m,pole) result(f)
    real(dp), intent(in) :: t
    complex(dp), intent(in) :: d,pole
    integer, intent(in) :: m
    f=d*(2*real(m,dp)-d)-pole/cosh(t)**2
  end function

  ! Output uses transposed storage convention: the walk accesses
  ! trans(input_component,output_component,ring), likewise innov.
  subroutine riccati_mode_filter(z,m,c0,c2,innov,trans,rc,refinement)
    real(dp), intent(in) :: z(:),c0,c2
    integer, intent(in) :: m
    real(dp), intent(out) :: innov(2,2,size(z)),trans(2,2,size(z)-1)
    integer, intent(out) :: rc
    integer, optional, intent(in) :: refinement
    real(dp) :: cov(2,2,size(z)),cross(2,2,size(z)-1),chol(2,2),at(2,2),q(2,2),v(2)
    integer :: k,i
    innov=0; trans=0
    call riccati_mode_covariances(z,m,c0,c2,cov,cross,rc,refinement)
    if(rc/=0) return
    call chol2(cov(:,:,1),chol,rc)
    if(rc/=0) return
    innov(:,:,1)=transpose(chol)
    do k=1,size(z)-1
      call chol2(cov(:,:,k),chol,rc)
      if(rc/=0) return
      do i=1,2
        v(1)=cross(i,1,k)/chol(1,1)
        v(2)=(cross(i,2,k)-chol(2,1)*v(1))/chol(2,2)
        at(i,2)=v(2)/chol(2,2)
        at(i,1)=(v(1)-chol(2,1)*at(i,2))/chol(1,1)
      end do
      q=cov(:,:,k+1)-matmul(at,transpose(cross(:,:,k)))
      q=0.5_dp*(q+transpose(q))
      call chol2(q,chol,rc)
      if(rc/=0) return
      innov(:,:,k+1)=transpose(chol); trans(:,:,k)=transpose(at)
    end do
  end subroutine

  subroutine chol2(a,l,rc)
    real(dp), intent(in) :: a(2,2)
    real(dp), intent(out) :: l(2,2)
    integer, intent(out) :: rc
    real(dp) :: rem
    l=0; rc=3
    if(.not.all(ieee_is_finite(a))) return
    if(a(1,1)<=0.or.a(2,2)<=0) return
    l(1,1)=sqrt(a(1,1)); l(2,1)=a(2,1)/l(1,1)
    rem=a(2,2)-l(2,1)**2
    ! Never hide a failed covariance by adding jitter or clipping eigenvalues.
    if(rem<=0) return
    l(2,2)=sqrt(rem); rc=0
  end subroutine

  ! A uniform bound on omitted point variance. Modes m>M are a subset of
  ! spherical harmonics l>M. Their variance is bounded by
  ! 1/(4*pi*c2*(M+1)^2). A finite positive sum gives a lower bound on
  ! the full point variance, hence a conservative relative-error criterion.
  subroutine spherical_mode_cutoff(c0,c2,relative_tol,mmax,rc)
    real(dp),intent(in)::c0,c2,relative_tol
    integer,intent(out)::mmax,rc
    real(dp)::total,ell,needed
    integer::l,last
    mmax=0;rc=1
    if(c0<=0.or.c2<=0.or.relative_tol<=0.or.relative_tol>=1) return
    if(.not.ieee_is_finite(c0).or..not.ieee_is_finite(c2)) return
    last=max(256,ceiling(10*(c0/c2)**0.25_dp)); total=0
    do l=0,last
      ell=real(l,dp)
      total=total+(2*ell+1)/(c0+c2*(ell*(ell+1))**2)
    end do
    needed=sqrt(1/(c2*total*relative_tol))
    if(.not.ieee_is_finite(needed).or.needed>1000000) return
    mmax=max(1,ceiling(needed));rc=0
  end subroutine

  ! Grid-independent calibration. Bound the unsummed positive spectrum using
  ! sum_{l>L}(2l+1)/(c0+c2*l^2*(l+1)^2) <= 1/(c2*(L+1)^2).
  ! Hence the infinite-spectrum correlation differs by <=2*tail/sum0.
  subroutine calibrate_spherical_c2(xcorr_deg,c0,c2,rc,tail_tol)
    real(dp), intent(in) :: xcorr_deg,c0
    real(dp), intent(out) :: c2
    integer, intent(out) :: rc
    real(dp), optional, intent(in) :: tail_tol
    real(dp), allocatable :: w(:),wtheta(:),u2(:)
    real(dp) :: angle,ct,p0,p1,p2,lo,hi,mid,cor,total,bound,tol,wgt,goal
    integer :: l,lmax,newmax,iter,outer
    c2=0; rc=1; tol=1.e-6_dp
    if(present(tail_tol)) tol=tail_tol
    if(.not.ieee_is_finite(xcorr_deg).or..not.ieee_is_finite(c0)) return
    if(xcorr_deg<=0.or.xcorr_deg>=180.or.c0<=0.or.tol<=0) return
    angle=xcorr_deg*pi/180; ct=cos(angle); goal=exp(-0.5_dp)
    lmax=max(512,ceiling(10.0_dp/angle)); rc=4
    do outer=1,8
      if(lmax>2000000) return
      allocate(w(0:lmax),wtheta(0:lmax),u2(0:lmax))
      p0=1; p1=ct
      do l=0,lmax
        w(l)=real(2*l+1,dp); u2(l)=(real(l,dp)*real(l+1,dp))**2
        if(l==0) then
          wtheta(l)=w(l)
        else if(l==1) then
          wtheta(l)=w(l)*ct
        else
          p2=(real(2*l-1,dp)*ct*p1-real(l-1,dp)*p0)/real(l,dp)
          wtheta(l)=w(l)*p2; p0=p1; p1=p2
        end if
      end do
      lo=log(c0)+4*log(angle)-log(1.e4_dp)
      hi=log(c0)+4*log(angle)+log(1.e4_dp)
      call sums(exp(lo),cor,total)
      if(cor>goal) then
        deallocate(w,wtheta,u2); lmax=2*lmax; cycle
      end if
      call sums(exp(hi),cor,total)
      if(cor<goal) then
        deallocate(w,wtheta,u2); return
      end if
      do iter=1,80
        mid=0.5_dp*(lo+hi); c2=exp(mid)
        call sums(c2,cor,total)
        if(abs(cor-goal)<1.e-11_dp) exit
        if(cor>goal) then
          hi=mid
        else
          lo=mid
        end if
      end do
      bound=2/(c2*real(lmax+1,dp)**2*total)
      if(bound<=tol) then
        deallocate(w,wtheta,u2); rc=0; return
      end if
      newmax=max(2*lmax,ceiling(sqrt(2/(c2*total*tol))*1.02_dp))
      deallocate(w,wtheta,u2); lmax=newmax
    end do
  contains
    subroutine sums(c,rat,total)
      real(dp),intent(in)::c
      real(dp),intent(out)::rat,total
      real(dp)::num,den
      integer::j
      total=0; num=0
      do j=0,lmax
        den=c0+c*u2(j)
        total=total+w(j)/den
        num=num+wtheta(j)/den
      end do
      rat=num/total
    end subroutine
  end subroutine
end module
