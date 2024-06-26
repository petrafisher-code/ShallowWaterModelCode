      !>@author
      !>Netlib 2019
      !>@copyright Public Domain
      !>@brief
      !> Find a minimum between bounds (single precision version)
      !>Downloaded from Netlib, 2019, 
      !>@param[in] ax: lower bound
      !>@param[in] bx: upper bound
      !>@param[in] f: function to minimize
      !>@param[in] tol: tolerance
      !>@return sfmin: location of minimum
c  To get r1mach, mail netlib
c       send r1mach from core
      real function sfmin(ax,bx,f,tol)
      real ax,bx,f,tol
      external f
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result (.ge.0.)
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a
c        minimum
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(sfmin)+(tol/3), where eps is  approximately  the  square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated  by  at least  eps*abs(x)+(tol/3), then  sfmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(sfmin)+tol.  if   f   is  not  unimodal,
c  then sfmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
c
      real  a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w,fu,fv,fw,
     2    fx,x,tol3
      real   abs, sqrt,r1mach
c
c  c is the squared inverse of the golden ratio
      c=0.5e0*(3.0e0- sqrt(5.0e0))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
   10 eps=r1mach(4)
      tol1=eps+1.0e0
      eps= sqrt(eps)
c
      a=ax
      b=bx
      v=a+c*(b-a)
      w=v
      x=v
      e=0.0e0
      fx=f(x)
      fv=fx
      fw=fx
      tol3=tol/3.0e0
c
c  main loop starts here
c
   20 xm=0.5e0*(a+b)
      tol1=eps* abs(x)+tol3
      t2=2.0e0*tol1
c
c  check stopping criterion
c
      if ( abs(x-xm).le.(t2-0.5e0*(b-a))) go to 190
      p=0.0e0
      q=0.0e0
      r=0.0e0
      if ( abs(e).le.tol1) go to 50
c
c  fit parabola
c
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0e0*(q-r)
      if (q.le.0.0e0) go to 30
      p=-p
      go to 40
   30 q=-q
   40 r=e
      e=d
   50 if (( abs(p).ge. abs(0.5e0*q*r)).or.(p.le.q*(a-x))
     2          .or.(p.ge.q*(b-x))) go to 60
c
c  a parabolic-interpolation step
c
      d=p/q
      u=x+d
c
c  f must not be evaluated too close to ax or bx
c
      if (((u-a).ge.t2).and.((b-u).ge.t2)) go to 90
      d=tol1
      if (x.ge.xm) d=-d
      go to 90
c
c  a golden-section step
c
   60 if (x.ge.xm) go to 70
      e=b-x
      go to 80
   70 e=a-x
   80 d=c*e
c
c  f must not be evaluated too close to x
c
   90 if ( abs(d).lt.tol1) go to 100
      u=x+d
      go to 120
  100 if (d.le.0.0e0) go to 110
      u=x+tol1
      go to 120
  110 u=x-tol1
  120 fu=f(u)
c
c  update  a, b, v, w, and x
c
      if (fx.gt.fu) go to 140
      if (u.ge.x) go to 130
      a=u
      go to 140
  130 b=u
  140 if (fu.gt.fx) go to 170
      if (u.ge.x) go to 150
      b=x
      go to 160
  150 a=x
  160 v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
      go to 20
  170 if ((fu.gt.fw).and.(w.ne.x)) go to 180
      v=w
      fv=fw
      w=u
      fw=fu
      go to 20
  180 if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) go to 20
      v=u
      fv=fu
      go to 20
c
c  end of main loop
c
  190 sfmin=x
      return
      end