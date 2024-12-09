c------------------------------------------------------------------------------
c      program readspline.f
c
c---------------------------------------------------------------------
c
c     Outputs CHAOS4 field values from spherical harmonic coefficients
c     (intermediate output) interpolated from an order 6 spline
c
c     adapted from f77 code for CHAOS (order 4 spline) from DTU Space
c     by SMAC Sep 2010
c     Extended to calculate derivatives by RTH 11/5/2011, and simplified
c     slightly to enable easy conversion to other spline degrees.
c
c---------------------------------------------------------------------
c
c     uses code by: Jeremy Bloxham, Andrew Jackson,
c     Rick O'Connell, and Carl de Boor
c
c---------------------------------------------------------------------
c     for details of B-splines see: 
c            Carl de Boor "A Practical Guide to Splines"
c            Springer-Verlag, 1978., 2nd edition 2000
c---------------------------------------------------------------------
c
c     lmax   is maximum degree of spherical harmonics
c     nspl   is number of B-splines
c
c     gt     is the full array of coefficients in the expansion of
c               the geomagnetic potential in spherical harmonics
c               and cubic B-splines
c
c     g      is an array of geomagnetic main field coefficients at a
c               particular point in dt
c
c     gd     is an array of geomagnetic secular variation coefficients 
c               at a particular point in dt. Higher derivatives follow
c               the same structure
c
c     tknts  is an array of knot points
c
c---------------------------------------------------------------------
c
c     jord  = 6 (order of B-splines)
c
c---------------------------------------------------------------------
c
c     CALLS:    interv   - calculates which knot lies immediately left
c                          of the current dt point
c
c               bspline  - calculates B-splines at current dt point
c
c               bspline1 - calculates first dt derivative of 
c                          B-splines at current dt point
c
c               b0, b1   - functions required by bspline1
c
c---------------------------------------------------------------------- 
      implicit none
      
      integer lmax,nspl,n,np,nl,jord
      integer it,ns1,dermax

      parameter (lmax=14)
      parameter (nspl=95)
      parameter (n=lmax*(lmax+2))
      parameter (np=n*nspl)
      parameter (nl=(lmax+1)*(lmax+2)/2)
      parameter (jord=4)
      parameter (dermax = 2)

      dimension gt(n,nspl),gtd(n,nspl)
      dimension spl(nspl),tknts(nspl+6)      
      real*8 g(n,0:dermax),gta(n,nspl,0:dermax)

      integer jordin,lm,ns,k,nleft,j,der
      
      real*8 gt,gtd,spl,tknts,dt
      
      CHARACTER(LEN=8) :: date
      integer lw,mw,sw
c*********************************************************************
c----------------------
c    input spline file
c
c change the path to the correct model and the name of the output file at the end of the program
c----------------------
      open(1,file=
     &'/nfs/see-fs-01_users/earsmaf/geomag/models
     &/COVOBS.x1/COV-OBS.x1-int')

     
c-----
c     input current model
      
      read(1,*) 
      read(1,*) lm,ns,jordin
      read(1,*) (tknts(k),k=1,ns+jordin)
      read(1,*) ((gt(k,ns1), k = 1,lm*(lm+2)), ns1 = 1,ns)
            
      if (ns .ne. nspl .or. lm .ne. lmax .or. jord .ne. jordin)
     &    stop 'Model and code do not match'

c     copy these coefficients to the more general array
c     for coefficients of field and derivatives

      do j=1,nspl
        do k=1,n
          gta(k,j,0) = gt(k,j)
        end do
      end do

c-----
c     calculate secular variation coefficients and higher
c     derivatives (to degree dermax) if required at dt
c     first calculate the derivative of the spline coefficients
c     by differencing following de Boor, p116
c     obtain higher degrees by iteration
c     der gives the order of differentiation
      do der=1,dermax
        call derivcoeff(n,nspl,jord,der,tknts,gt,gtd)
        do j=1,nspl
          do k=1,n
            gta(k,j,der) = gtd(k,j)
c           reuse gt for ease of iteration
            gt(k,j) = gtd(k,j)
          end do
        end do
      end do


c-----     open(2,file=)

c     calculate main field and time deriv coefficients at dt
c   
c      do it = 120,120
      do it = 0,0 !96,96
        PRINT*,'Enter date'
        read*, dt
    !    dt = 2009 ! +it/12.
        call interv(tknts,dt,ns,nleft,jord)

        do der = 0,dermax
c         calculate the B-splines, reduced for the derivatives
          call bspline(tknts,dt,ns,jord-der,nleft,
     &                         spl(nleft-jord+1+der))

          do  k=1,n
            g(k,der)=0.0
            do j=1,jord-der
              g(k,der) = g(k,der)
     &           + spl(j+nleft-jord+der)*gta(k,j+nleft-jord+der,der)
            enddo 
          enddo 
        end do
         PRINT*, dt
         
c----------------
c open output file
c----------------
	write(date , '(f8.3)') dt
	open(2,file='COVOBS.x1_'//date//'.dat')
	write(2,*) 'COV-OBS.x1-int Gauss coefficients for '//date
	lw=1
	mw=0
	sw=0

	do  k=1,n
          write (6,'(8f18.6)') (g(k,j),j=0,dermax)
          
	  if (sw==0) then
	    if (mw==0) then
	      write (2,'(2i5,2f18.6)') lw, mw, g(k,0), 0
c	don't skip the next writing cycle	      
	      sw = 0
	    else
	      write (2,'(2i5,2f18.6)') lw, mw, g(k,0), g(k+1,0)
c	I wrote both g and h: skip the next writing cycle	      
	      sw = 1
	    endif
            mw=mw+1;
            if (mw.gt.lw) then
              lw = lw +1;
              mw=0;
            endif
          elseif (sw==1) then
c	I skipped a writing cycle. Don't skip the next
            sw = 0
          endif
          
          
        enddo
       

        close(2)
        
      end do

      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine interv(tknts,dt,nspl,nleft,jord)

      implicit real*8 (a-h,o-z)
 
      dimension tknts(nspl+jord)
       
c  knots equally spaced (no stripping of repeated knots needed)
c  interval starts at tknts(jord) and ends at tknts(nspl+1)

c-----
c    calculate nleft: 
c                  tknts(nleft) < tknts(nleft+1)
c                  tknts(nleft) <= dt <= tknts(nleft+1)
c 
  
c  check we are in-range
      if(dt.lt.tknts(jord).or.dt.gt.tknts(nspl+1)) return
      
      do 200 n=jord+1,nspl+1
       if(dt.le.tknts(n)) then
        nleft=n-1
        goto 210
       endif
200   continue
210   continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bspline(tknts,t,nspl,jorder,nleft,spl)

c calculate splines of order jorder
       
      implicit real*8 (a-h,o-z)
      dimension tknts(nspl+jorder)
      dimension spl(jorder)
       
      dimension deltal(jorder),deltar(jorder)
       
      spl(1)=1.0
      
      do 200 j=1,jorder-1
      
      deltar(j) = tknts(nleft+j) - t
      deltal(j) = t - tknts(nleft+1-j)
      saved=0.0
       
      do 100 i=1,j
        term = spl(i)/(deltar(i)+deltal(j+1-i))
        spl(i) = saved + deltar(i)*term
        saved = deltal(j+1-i)*term
100   continue

      spl(j+1) = saved
       
200   continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivcoeff(n,nspl,jordin,der,tknts,gt,gtd)

c     calculate first derivative by deriving a new set of
c     coefficients

      implicit real*8 (a-h,o-z)
      dimension tknts(nspl+jordin)
      dimension gt(n,nspl),gtd(n,nspl)
      dimension fact(nspl)

      integer der

c     der defines the level of derivative sought, jordin
c     the order of the original spline

      jord = jordin - der + 1

c     first calculate a vector for the time differentials
      delt0 = tknts(nspl)-tknts(1)
      do j=1,nspl
        delt = ( tknts(j+jord-1) - tknts(j) )
        if (abs(delt/delt0) .lt. 1.0e-6) then
c         Following deBoor p117 that anything multiplied by
c         nothing is nothing, getting rid of the singularity
          fact(j) = 0.0
        else
          fact(j) = 1.0/delt
        endif
      end do
      do k=1,n
c       Need a special case for j=1 (assume 0 bspline is 0)
        gtd(k,1) = (jord-1)*gt(k,1)*fact(1)
        do j=2,nspl
          gtd(k,j) = (jord-1)*(gt(k,j)-gt(k,j-1))*fact(j)
        end do
      end do
      
      return

      end
