      subroutine seakenlcd(xx, cx, yy, cy, n, results)  
c     seasonal adjustment to Kendall test for trend of left-censored data.
c
c     Arguments:
c        xx is the adjusted x-values, with NAs removed
c        cx is the censor codes for x
c        yy is the adjusted y-values, with NAs removed
c        cy is the censor codes for y
c        n is the length xx, cx, yy, and cy
c        the xs are from one season, the ys from the other
c
c     Output:
c        results is the output 
c           results = sigmagh
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
      double precision xx(*), yy(*), results(2)
      integer cx(*),  cy(*)
      double precision kgh, rgh, sg, sh, signum
c
c     accumulate kgh and rgh
c
      kgh = 0.d0
      rgh = 0.d0
      do j=2,n
c
c     kgh part
c
         if(cx(j) .eq. -1) goto 20
         if(cy(j) .eq. -1) goto 20
         do i=1,j
            if(cx(i) .eq. -1) goto 10
            if(cy(i) .eq. -1) goto 10
            if(cx(j) .eq. 0) then
               if(cx(i) .eq. 0) then
                  sg = signum(xx(j) - xx(i)) ! eq 2 (both uncensored)
               else
                  if(xx(j) .ge. xx(i)) then  ! eq 5 (j un, i cen)
                     sg = 1.d0
                  else
                     sg = 0.d0
                  endif
               endif
            else
               if(cx(i) .eq. 1) then 
                  sg = 0.d0                  ! eq 3 (both censored)
               else
                  if(xx(j) .le. xx(i)) then  ! eq 4 (j cen, i un)
                     sg = -1.d0
                  else
                     sg = 0.d0
                  endif
               endif
            endif
            if(cy(j) .eq. 0) then
               if(cy(i) .eq. 0) then
                  sh = signum(yy(j) - yy(i)) ! eq 2
               else
                  if(yy(j) .ge. yy(i)) then  ! eq 5
                     sh = 1.d0
                  else
                     sh = 0.d0
                  endif
               endif
            else
               if(cy(i) .eq. 1) then        ! eq 3
                  sh = 0.d0
               else
                  if(yy(j) .le. yy(i)) then ! eq 4
                     sh = -1.d0
                  else
                     sh = 0.d0
                  endif
               endif
            endif
            kgh = kgh + sg * sh
 10         continue
         enddo
 20      continue
      enddo
c
c     rgh part
c
      do j=1,n
         if(cx(j) .eq. -1) goto 40 ! skip missing values
         if(cy(j) .eq. -1) goto 40
         do i=1,n
            if(i .eq. j) goto 30
            if(cx(i) .eq. -1) goto 30 ! skip missing values
            if(cx(j) .eq. 0) then ! always sg = 0
               if(cx(i) .eq. 0) then
                  sg = signum(xx(j) - xx(i)) ! eq 2
               else
                  if(xx(j) .ge. xx(i)) then  ! eq 5
                     sg = 1.d0
                  else
                     sg = 0.d0
                  endif
               endif
            else
               if(cx(i) .eq. 1) then
                  sg = 0.d0                  ! eq 3
               else
                  if(xx(j) .le. xx(i)) then  ! eq 4
                     sg = -1.d0
                  else
                     sg = 0.d0
                  endif
               endif
            endif
            if(sg .ne. 0.d0) then ! no reason to check 
               do k=1,n
                  if(j .eq. k) goto 25 ! always sh = 0
                  if(cy(k) .eq. -1) goto 25 ! skip missing values
                  if(cy(j) .eq. 0) then
                     if(cy(k) .eq. 0) then
                        sh = signum(yy(j) - yy(k)) ! eq 2
                     else
                        if(yy(j) .ge. yy(k)) then ! eq 5
                           sh = 1.d0
                        else
                           sh = 0.d0
                        endif
                     endif
                  else
                     if(cy(k) .eq. 1) then       ! eq 3
                        sh = 0.d0
                     else
                        if(yy(j) .le. yy(k)) then ! eq 4
                           sh = -1.d0
                        else
                           sh = 0.d0
                        endif
                     endif
                  endif
                  rgh = rgh + sg * sh
 25               continue
               enddo
            endif
 30         continue
         enddo
 40      continue
      enddo
      results(1) = kgh 
      results(2) = rgh
      return
      end
