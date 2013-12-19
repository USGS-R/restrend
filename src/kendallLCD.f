      subroutine kendalllcd(xx, cx, yy, cy, n, results, slopel, 
     $     slopeh, nsl)
c     Kendall test for trend of left-censored data.
c
c     Arguments:
c        xx is the adjusted x-values, with NAs removed
c        cx is the censor codes for x
c        yy is the adjusted y-values, with NAs removed
c        cy is the censor codes for y
c        n is the length xx, cx, yy, and cy
c        nsl indicates whether slopes are desired 0 -> yes, <0 -> no
c
c     Output:
c        results is the output vector containing
c           results(1) = s
c           results(2) = n*(n-1) = 2 * J
c           results(3) = tt
c           results(4) = uu
c           results(5) = vars
c           results(6) = adjustment to vars
c        slopel is a vector of the lowest possible slopes
c        slopeh is a vector of the upper possible slopes
c               both must be dimension(max(n,n*(n/nseas-1)/2))
c        nsl is the number of slopes
c
c     Limitations:
c        There must be at least two years of data.
c
c     References:
c     1) A Nonparametric Trend Test for Seasonal Data with Serial
c          Dependence, R.M.Hirsch & J.R.Slack, Water Resources Research,
c          Vol.20, No.6, Pages 727-732, June 1984.
c     2) Techniques of Trend Analysis for Monthly Water Quality Data,
c          R.M.Hirsch, J.R.Slack, & R.A.Smith, Water Resources Research,
c          Vol.18, No.1, Pages 107-121, February 1982.
c
c     Variables names correspond to the original code in S, which is based
c        on the MiniTab code by Dennis Helsel and Ed Gilroy.
c
c     Code history:
c	 22 Nov 05  DLLorenz Initial Conversion to left-censored data
c
      double precision xx(*), cx(*), yy(*), cy(*), results(6), slopel(*)
      double precision slopeh(*)
      double precision, dimension (n,n) :: diffx, diffcx, xplus, diffy, 
     $     diffcy, yplus, signyx
      double precision, dimension (n) :: dxx, dcx, dyy, dcy
      double precision, dimension (4) :: xties, yties
      double precision tt, uu, cix, ciy, tplus, uplus, tot, xpl, ypl
      double precision xdif, y1l, y2l, y1h, y2h
      integer, dimension (n) :: dorder
      double precision signum

      call outersub(xx, n, xx, n, diffx)
      call outersub(cx, n, cx, n, diffcx)
      call outeradd(cx, n, cx, n, xplus)
      call outersub(yy, n, yy, n, diffy)
      call outersub(cy, n, cy, n, diffcy)
      call outeradd(cy, n, cy, n, yplus)
      call vectmult(diffy, diffx, n, n, signyx)
      call vectsign(signyx, n, n, signyx)
c
c     compute # ties in x (tt) and # ties in y (uu)
c
      tt = 0.d0
      uu = 0.d0
      do i=1,n
         do j=1,n
            tt = tt + 1.d0 - abs(signum(diffx(i,j)))
            uu = uu + 1.d0 - abs(signum(diffy(i,j)))
         enddo
      enddo
      tt = (tt - n)/2.d0
      uu = (uu - n)/2.d0
c
c     adjust tt and uu
c
      cix = 0.d0
      ciy = 0.d0
      do i=1,n
         do j=1,n
            if(diffcx(i,j) * diffx(i,j) .gt. 0.d0) then
               cix = cix + 1.d0
               signyx(i,j) = 0.d0
            endif
            if(diffcy(i,j) * diffy(i,j) .gt. 0.d0) then
               ciy = ciy + 1.d0
               signyx(i,j) = 0.d0
            endif
         enddo
      enddo
      tt = tt + cix/2.d0
      uu = uu + ciy/2.d0
c
c     one final adjustment to tt and uu; compute concordances
c
      tplus = 0.d0
      uplus = 0.d0
      tot = 0.d0
      do i=1,n
         do j=1,n
            xpl = 0.d0
            if(xplus(i,j) .gt. 1.d0) xpl = 1.d0
            ypl = 0.d0
            if(yplus(i,j) .gt. 1.d0) ypl = 1.d0
            tot = tot + signyx(i,j) * (1.d0 - xpl) * (1.d0 - ypl)
            if(diffx(i,j) .eq. 0.d0) xpl = 0.d0
            if(diffy(i,j) .eq. 0.d0) ypl = 0.d0
            tplus = tplus + xpl
            uplus = uplus + ypl
         enddo
      enddo
      tt = tt + tplus/2.d0
      uu = uu + uplus/2.d0
c
c     first part of results
c
      results(1) = tot/2.d0
      results(2) = n * (n - 1)
      results(3) = tt
      results(4) = uu
      results(5) = results(2) * (2.d0 * n + 5.d0) / 18.d0
c
c     adjust vars for ties, from tiesboth.mac by Ed Gilroy
c
c     These code easily if the data are sorted
c
      call porder(xx, n, dorder)
      do i=1,n
         dxx(i) = xx(dorder(i))
         dcx(i) = cx(dorder(i))
      enddo
      call porder(yy, n, dorder)
      do i=1,n
         dyy(i) = yy(dorder(i))
         dcy(i) = cy(dorder(i))
      enddo
c
c     adjust for censored ties
c
      call ccties(dxx, dcx, n, xties, 1.d0)
      call ccties(dyy, dcy, n, yties, 1.d0)
      results(6) = (xties(1) + yties(1)) / 18.d0
      results(6) = results(6) - xties(2) * yties(2) / (9 * n * (n - 1)
     $     * (n - 2)) - xties(3) * yties(3) / (2 * n * (n - 1))
CCC DLL removed 05/04/06     $     - xties(4) - yties(4)
c
c     adjust for censored data > uncensored data
c
      call cucties(dxx, dcx, n, xties)
      call cucties(dyy, dcy, n, yties)
      results(6) = results(6) + (xties(1) + yties(1)) / 18.d0
      results(6) = results(6) - xties(2) * yties(2) / (9 * n * (n - 1)
     $     * (n - 2)) - xties(3) * yties(3) / (2 * n * (n - 1))
c
c     adjust for uncensored ties
c
      call ccties(dxx, dcx, n, xties, 0.d0)
      call ccties(dyy, dcy, n, yties, 0.d0)
      results(6) = results(6) + (xties(1) + yties(1)) / 18.d0
      results(6) = results(6) - xties(2) * yties(2) / (9 * n * (n - 1)
     $     * (n - 2)) - xties(3) * yties(3) / (2 * n * (n - 1))
      if(nsl .lt. 0) return
c
c     compute slopes, requires that x is uncensored and untied
c
      nsl = 0
      if(cx(1) .gt. 0.9d0) return
      do i=2,n
         if(cx(i) .gt. 0.9d0) return
         if(dxx(i) .eq. dxx(i-1)) return
      enddo
      do i=1,n
         do 10 j=(i+1),n
C*            if(cy(i) .gt. 0.9d0 .and. cy(j) .gt. 0.9d0) goto 10
C*            if(cy(i) .gt. 0.9d0 .and. yy(i) .gt. yy(j)) goto 10
C*            if(cy(j) .gt. 0.9d0 .and. yy(j) .gt. yy(i)) goto 10
            xdif = xx(j) - xx(i)
            y1h = yy(j)
            y2h = yy(i)
            if(cy(i) .gt. 0.9d0) then
               y2l = 0.d0
            else
               y2l = y2h
            endif
            if(cy(j) .gt. 0.9d0) then
               y1l = 0.d0
            else
               y1l = y1h
            endif
            nsl = nsl + 1
            slopel(nsl) = (y1l - y2h) / xdif
            slopeh(nsl) = (y1h - y2l) / xdif
C$            if(cy(i) .gt. 0.9d0 .and. cy(j) .gt. 0.9d0) then
C$               slopel(nsl) = 0.d0
C$               slopeh(nsl) = 0.d0
C$            endif
C$            if(cy(i) .gt. 0.9d0 .and. yy(i) .gt. yy(j)) then
C$               slopel(nsl) = 0.d0
C$               slopeh(nsl) = 0.d0
C$            endif
C$            if(cy(j) .gt. 0.9d0 .and. yy(j) .gt. yy(i)) then
C$               slopel(nsl) = 0.d0
C$               slopeh(nsl) = 0.d0
C$            endif
 10      continue
      enddo
      return
      end
c
c     end of main code
c
      subroutine cucties(xx, cx, n, ties)
      double precision xx(*), cx(*), ties(*)
      double precision xsum, base
      do i=1,4
         ties(i) = 0.d0
      enddo
c
c     step through the data counting index of censored data
c
      xsum = 0.d0
      base = 0.d0
      do i=1,n
         if(base .eq. 0.d0) then
            if(cx(i) .eq. 0.d0) base = i
         else
            if(cx(i) .eq. 1.d0) xsum=xsum + i - base
         endif
      enddo
      ties(1) = xsum * 18.d0
      ties(3) = xsum * 2.d0
      return
      end
c
      subroutine ccties(xx, cx, n, ties, test)
      double precision xx(*), cx(*), ties(*), test
      double precision xlast

      do i=1,4
         ties(i) = 0.d0
      enddo
c
c     step through the data counting ties of censored data
c
      nties = 1
      xlast = xx(1)
      do i=2,n
         if(cx(i) .ne. test) then ! force accumulation
            ties(1) = ties(1) + nties * (nties - 1) * (nties * 2 + 5)
            ties(2) = ties(2) + nties * (nties - 1) * (nties - 2)
            ties(3) = ties(3) + nties * (nties - 1)
            ties(4) = ties(4) + nties - 1
            nties = 1
            xlast = xx(i)
         else if(xx(i) .ne. xlast) then ! force accum
            ties(1) = ties(1) + nties * (nties - 1) * (nties * 2 + 5)
            ties(2) = ties(2) + nties * (nties - 1) * (nties - 2)
            ties(3) = ties(3) + nties * (nties - 1)
            ties(4) = ties(4) + nties - 1
            nties = 1
            xlast = xx(i)
         else ! must be tie
            nties = nties + 1
         endif
      enddo
      ties(1) = ties(1) + nties * (nties - 1) * (nties * 2 + 5)
      ties(2) = ties(2) + nties * (nties - 1) * (nties - 2)
      ties(3) = ties(3) + nties * (nties - 1)
      ties(4) = ties(4) + nties - 1
      return
      end
c
      subroutine outersub(x, n, y, m, diff)
      double precision x(n), y(m), diff(n,m)
      do i = 1,n
         do j = 1,m
            diff(i,j) = x(i)-y(j)
         enddo
      enddo
      return
      end
c
      subroutine outeradd(x, n, y, m, sum)
      double precision x(n), y(m), sum(n,m)
      do i = 1,n
         do j = 1,m
            sum(i,j) = x(i)+y(j)
         enddo
      enddo
      return
      end
c
      subroutine vectmult(x, y, n, m, out)
      integer n, m
      double precision x(n,m), y(n,m), out(n,m)
      do i = 1,n
         do j = 1,m
            out(i,j) = x(i,j) * y(i,j)
         enddo
      enddo
      return
      end
c
      subroutine vectsign(x, n, m, out)
      double precision x(n,m), out(n,m)
      double precision signum
      do i = 1,n
         do j = 1,m
            out(i,j) = signum(x(i,j))
         enddo
      enddo
      return
      end
c
c     sorted x is
c     xs(I) = x(ix(i))
c
c
      double precision function signum(x)
      double precision x
      if(x .gt. 0.d0) then 
         signum = 1.d0
      else if(x .lt. 0.d0) then
         signum = -1.d0
      else
         signum = 0.d0
      endif
      return
      end
