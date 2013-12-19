      subroutine regionsk(x, Nsite, Nseas, Nyear, results, scomp, slope)
c This subroutine computes the covariance adjustment to S for regional
c seasonal Kendall test. 
c The argument x is an matrix composed of rows of the individual 
c observations at each site and columns of sites.
c It calls seakenI so that routine computes the regional covariance for 
c each season. The results are returned in results:
c      1 S, 2 var(S), 3 serialcov(S), 4 spatialcov(S)
      dimension x(Nyear*Nseas, Nsite), results(4)
      dimension partres(5)
      real, dimension (Nyear*Nsite) :: xtosk
      real, dimension (Nsite) :: slopes
      do i = 1, 4
         results(i) = 0.0
      enddo
      scomp = 0.0
c
c Compute the seasonal Kendall test on each site, accumulating the
c S statistics
c
      do k=1,Nsite
         call seakeni(x(1, k), Nyear*Nseas, Nseas, partres)
         results(1) = results(1) + partres(3)
         results(2) = results(2) + partres(4)
         results(3) = results(3) + partres(5)
         slopes(k) = partres(2)
         scomp = scomp + partres(1)
      enddo
c
c Compute the seasonal Kendall test on each season, giving the
c spatial covariance.
c
      do j=1,Nseas
         iout = 0
         do i=1,Nyear
            do k=1,Nsite
               iout = iout + 1
               xtosk(iout) = x(j + (i-1)*Nseas, k)
            enddo
         enddo
         call seakeni(xtosk, Nyear*Nsite, Nsite, partres)
         results(4) = results(4) + partres(5)
      enddo
c
c compute the median of the median slopes
c
      call vssort(slopes, Nsite)
      slope = (slopes((Nsite + 1)/2) + slopes((Nsite + 2)/2)) / 2.0
      return
      end
