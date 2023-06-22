cc -*- mode: fortran; kept-new-versions: 25; kept-old-versions: 20 -*-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  rrcov : Scalable Robust Estimators with High Breakdown Point
cc
cc  This program is free software; you can redistribute it and/or modify
cc  it under the terms of the GNU General Public License as published by
cc  the Free Software Foundation; either version 2 of the License, or
cc  (at your option) any later version.
cc
cc  This program is distributed in the hope that it will be useful,
cc  but WITHOUT ANY WARRANTY; without even the implied warranty of
cc  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
cc  GNU General Public License for more details.
cc
cc  You should have received a copy of the GNU General Public License
cc  along with this program; if not, write to the Free Software
cc  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
cc
cc  I would like to thank Peter Rousseeuw and Katrien van Driessen for
cc  providing the initial code of this function.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   Computes the MCD estimator of multivariate location and scatter.
cc   This estimator is given by the subset of h observations for which
cc   the determinant of their covariance matrix is minimal. The MCD
cc   location estimate is then the mean of those h points, and the MCD
cc   scatter estimate is their covariance matrix. This value of h may be
cc   chosen by the user; its default value is roughly n/2.
cc
cc   The MCD estimator was first introduced in:
cc
cc	Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
cc	Journal of the American Statistical Association, Vol. 79,
cc	pp. 871-881. [See page 877.]
cc
cc   The MCD is a robust estimator in the sense that the estimates are
cc   not unduly influenced by outliers in the data, even if there
cc   are many outliers. Its robustness was proved in:
cc
cc	Rousseeuw, P.J. (1985), "Multivariate Estimation with High
cc	Breakdown Point," in Mathematical Statistics and Applications,
cc	edited by  W. Grossmann, G. Pflug, I. Vincze, and W. Wertz.
cc	Dordrecht: Reidel Publishing Company, pp. 283-297.
cc
cc	Rousseeuw, P.J. and Leroy, A.M. (1987), Robust Regression and
cc	Outlier Detection, Wiley-Interscience, New York. [Chapter 7]
cc
cc   The program also computes the distance of each observation
cc   from the center (location) of the data, relative to the shape
cc   (scatter) of the data:
cc
cc   * Using the classical estimates yields the Mahalanobis distance
cc     MD(i). Often, outlying points fail to have a large Mahalanobis
cc     distance because of the masking effect.
cc
cc   * Using the MCD estimates yields a robust distance RD(i).
cc     These distances allow us to easily identify the outliers.
cc
cc   For applications of robust distances in a regression context see:
cc
cc	Rousseeuw, P.J. and van Zomeren, B.C. (1990), "Unmasking
cc	Multivariate Outliers and Leverage Points," Journal of the
cc	American Statistical Association, Vol. 85, 633-639.
cc
cc   There also a diagnostic plot is given to distinguish between
cc   regular observations, vertical outliers, good leverage points,
cc   and bad leverage points.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   The new FAST_MCD algorithm introduced here is due to
cc
cc	Rousseeuw, P.J. and Van Driessen, K. (1997), "A Fast
cc	Algorithm for the Minimum Covariance Determinant
cc	Estimator," in preparation.
cc
cc   The algorithm works as follows:
cc
cc	 The dataset contains n cases, and nvar variables are used.
cc  Let n_0 := 2 * nmini (== 600).
cc	 When n <  n_0, the algorithm will analyze the dataset as a whole.
cc	 When n >= n_0, the algorithm will use several subdatasets.
cc
cc	 When the dataset is analyzed as a whole, a trial
cc	 subsample of nvar+1 cases is taken, of which the mean and
cc	 covariance matrix is calculated. The h cases with smallest
cc	 relative distances are used to calculate the next mean and
cc	 covariance matrix, and this cycle is repeated k1 times.
cc	 [For small n we can consider all subsets of nvar+1 out of n, else
cc	 the algorithm draws 500 random subsets.]
cc	 Afterwards, the best 10 solutions (covariance matrices and
cc	 corresponding means) are used as starting values for the final
cc	 iterations. These iterations stop when two subsequent determinants
cc	 become equal. (At most k3 iteration steps are taken.)
cc	 The solution with smallest determinant is retained.
cc
cc	 When the dataset contains more than 2*nmini cases, the algorithm
cc	 does part of the calculations on (at most) kmini nonoverlapping
cc	 subdatasets, of (roughly) nmini cases.
cc
cc	 Stage 1: For each trial subsample in each subdataset,
cc	 k1 iterations are carried out in that subdataset.
cc	 For each subdataset, the 10 best solutions are stored.
cc
cc	 Stage 2 considers the union of the subdatasets, called the
cc	 merged set. (If n is large, the merged set is a proper subset of
cc	 the entire dataset.) In this merged set, each of the 'best
cc	 solutions' of stage 1 are used as starting values for k2
cc	 iterations. Also here, the 10 best solutions are stored.
cc
cc	 Stage 3 depends on n, the total number of cases in the
cc	 dataset. If n <= 5000, all 10 preliminary solutions are iterated
cc	 k3 times. If n > 5000, only the best preliminary
cc	 solution is iterated, and the number of iterations decreases to 1
cc	 according to n*nvar. (If n*nvar <= 100,000 we iterate k3 times,
cc	 whereas for n*nvar > 1,000,000 we take only one iteration step.)
cc
cc   An important advantage of the algorithm FAST_MCD is that it allows
cc   for exact fit situations, where more than h observations lie on
cc   a hyperplane. Then the program still yields the MCD location and
cc   scatter matrix, the latter being singular (as it should be), as
cc   well as the equation of the hyperplane.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE FNANMCD(dat,n,nvar,nhalff,krep,initcov,initmean,
     *    inbest,det,fit,
     *    kount,temp,index1,index2,nmahad,ndist,am,am2,sd,
     *    means,bmeans,rec,sscp1,cova1,corr1,cinv1,cova2,cinv2,cstock,
     *    mstock,c1stock,m1stock,dath,xdat,ddd,xint,nlen,xdble,ndble,
     *    mvcode,xindex1,xindex2)

        implicit integer(i-n), double precision(a-h,o-z)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  ALGORITHM PARAMETERS:
cc
cc	To change the number of subdatasets and their size, the values of
cc	kmini and nmini can be changed.
cc
        parameter (kmini=5)
        parameter (nmini=300)
cc
cc	The number of iteration steps in stages 1,2 and 3 can be changed
cc	by adapting the parameters k1, k2, and k3.
cc
        parameter (k1=2)
        parameter (k2=2)
        parameter (k3=100)

c MM: below, '10' ("the ten best solutions") is also hardcoded in many places,
c	 related to "stock" !

cc
cc  krep := the total number of trial subsamples
cc          to be drawn when n exceeds 2*nmini;
cc  	was hardcoded krep := 500; now an *argument*
cc

cc  The following lines need not be modified.
cc
        parameter (km10=10*kmini)
        parameter (nmaxi=nmini*kmini)
cc
        integer rfncomb
        integer ierr,matz,tottimes,step
        integer pnsel
        integer flag(km10)
        integer mini(kmini)
        integer subdat(2,nmaxi)
        double precision mcdndex(10,2,kmini)
        integer subndex(450)
c FIXME                 ^^^ how does this depend on (kmini,nmini,...) ???
        integer replow
        integer fit
        double precision medi2

        integer inbest(nhalff)
        double precision dat(n,nvar)
        double precision initcov(nvar*nvar)
        double precision initmean(nvar)

        integer temp(n)
        integer index1(n)
        integer index2(n)
        double precision nmahad(n)
        double precision ndist(n)
        double precision am(n),am2(n)

        double precision sd(nvar)
        double precision means(nvar)
        double precision bmeans(nvar)

        double precision rec(nvar+1)
        double precision sscp1((nvar+1)*(nvar+1))
        double precision cova1(nvar*nvar)
        double precision corr1(nvar*nvar)
        double precision cinv1(nvar*nvar)
        double precision cova2(nvar*nvar)
        double precision cinv2(nvar*nvar)

        double precision cstock(10,nvar*nvar)
        double precision mstock(10,nvar)
        double precision c1stock(km10,nvar*nvar)
        double precision m1stock(km10,nvar*nvar)
        double precision dath(nmaxi,nvar)

C       Added for handling missing values - parameters and work
C       for calling the EM algo - EMNINT().
C       ddd = (2 + 3 * p + p^2)/2 - dimension of the packed storage
C       nint  = (p+1)*(p+1) + n*p + 3*n + 3*p - dim of the integer
C               srorage
C       ndble = 4*d + 2*p + n*p - dimension of the double storage
C       mvcode = missing value code (defaults to -99999)
C       xindex1, xilen1 = indexes and number of complete observations in
C                       the data matrix dat()
C       xindex2, xilen2 = same as above in the group data matrix dath()

        integer ddd, nlen, ndble, xint(nlen)
        integer xindex1(n), xindex2(n), xilen1, xilen2
        double precision xdble(ndble), mvcode
        double precision xdat(n*nvar)

C       Local variables
        double precision percen
        logical all,part,fine,final,rfodd,class
        logical complete

C       CALL INTPR('Entering RFFASTMCD - KREP: ',-1,KREP,1)

        call rndstart
C            -------- == GetRNGstate() in C
        nmax = n
        nvmax = nvar
        nvmax1=nvmax+1
        nvmax2=nvmax*nvmax
        nvm12=nvmax1*nvmax1

        nrep = krep
        part=.false.
        fine=.false.
        final=.false.
        all=.true.
        kstep=k1
        medi2=0
cc
cc  From here on, the sample size n is known.
cc  Some initializations can be made. First of all, h (= the number of
cc  observations on which the MCD is based) is given by the integer variable
cc  nhalff.
cc  If nhalff equals n, the MCD is the classical covariance matrix.
cc  The logical value class indicates this situation.
cc  The variable jbreak is the breakdown point of the MCD estimator
cc  based on nhalff observations, whereas jdefaul = (n+nvar+1)/2
cc  would be the optimal value of nhalff, with maximal breakdown point.
cc  The variable percen is the corresponding percentage (MM: rather "fraction").
cc
        percen = (1.D0*nhalff)/(1.D0*n)

        if(nvar.lt.5) then
          eps=1.0D-12
        else
          if(nvar.ge.5.and.nvar.le.8) then
            eps=1.0D-14
          else
            eps=1.0D-16
          endif
        endif

        class= .false.
        if(nhalff.ge.n) then
c	  compute *only* the classical estimate
          class= .true.
          goto 9500
        endif

cc  p >= 2   in the following
cc  ------
cc  Some initializations:
cc    seed = starting value for random generator
cc    matz = auxiliary variable for the subroutine rs, indicating whether
cc	     or not eigenvectors are calculated
cc    nsel = number of variables + 1
cc    ngroup = number of subdatasets
cc    part = logical value, true if the dataset is split up
cc    fine = logical value, becomes true when the subsets are merged
cc    final = logical value, to indicate the final stage of the algorithm
cc    all = logical value, true for small n, if all (p+1)-subsets out of
cc	    n can be drawn
cc    subdat = matrix with a first row containing indices of observations
cc	       and a second row indicating the corresponding subdataset
cc
        matz=1
        nsel=nvar+1
        ngroup=1
        part=.false.
        fine=.false.
        final=.false.
        all=.true.
        do 21,i=1,nmaxi
          subdat(1,i)=1000000
          subdat(2,i)=1000000
21      continue

C  VT::: Determine the number of complete observations in the data set -
C  this is necessary because the (p+1)-subsamples for initiating the
C  iterations are drawn from complete observations only.
C  This means also that for determining the number of trials nrep will
C  be based on the number of the complete observations, not on the
C  number of all observations n.
C  If there are not enough complete observations (<= 2*nvar) then allow to
C  use also incomplete observations for computing the initial estimates,
C  but with no more than one missing item. if still not sufficient set
C  fit=3 and exit.
        xilen1 = 0
        xilen2 = 0
        do 24 i=1,n
          complete = .TRUE.
          nx = 0
          do 25 j=1,nvar
            IF(dat(i,j) .EQ. mvcode) THEN
              complete = .FALSE.
              nx = nx + 1
            ENDIF
25        continue
          IF(complete) THEN
            xilen1 = xilen1 + 1
            xindex1(xilen1) = i
          ENDIF
C       Try observations with one missing item
C          IF(nx <= 1) THEN
C            xilen2 = xilen2 + 1
C            xindex2(xilen2) = i
C          ENDIF
24      continue

        IF(xilen1 .LE. 2*nsel) THEN
            CALL INTPR('NOT ENOUGH COMPLETE OBS.',-1,0,0)
            fit = 3
            goto 9999

C       not enough complete observations for drawing the
C       (p+1)-subsamples - take also these with only one missing item.
C          IF(xilen2 .LE. 2*nvar) THEN
C            CALL INTPR('NOT ENOUGH COMPLETE OBS.',-1, xilen2,1)
C            fit = 3
C            goto 9999
C          ENDIF
C          xilen1 = xilen2
C          DO 26 i=1,xilen2
C26          xindex1(i) = xindex2(i)

        ENDIF
cc
cc  Determine whether the dataset needs to be divided into subdatasets
cc  or can be treated as a whole. The subroutine rfrdraw constructs
cc  nonoverlapping subdatasets, with uniform distribution of the case numbers.
cc  For small n, the number of trial subsamples is determined.
cc
c MM(FIXME):  The following code depends crucially on  kmini == 5

        do 22 i=1,kmini
          mini(i)=0
22      continue
        if(n.gt.(2*nmini-1)) then
          kstep=k1
          part=.true.
          ngroup=int(n/(nmini*1.D0))
          if(n.ge.(2*nmini) .and. n.le.(3*nmini-1)) then
            if(rfodd(n)) then
              mini(1)=int(n/2)
              mini(2)=int(n/2)+1
            else
              mini(1)=n/2
              mini(2)=n/2
            endif
          else
            if(n.ge.(3*nmini) .and. n.le.(4*nmini-1)) then
              if(3*(n/3) .eq. n) then
                mini(1)=n/3
                mini(2)=n/3
                mini(3)=n/3
              else
                mini(1)=int(n/3)
                mini(2)=int(n/3)+1
                if(3*(n/3) .eq. n-1) then
                  mini(3)=int(n/3)
                else
                  mini(3)=int(n/3)+1
                endif
              endif
            else
              if(n.ge.(4*nmini) .and. n.le.(5*nmini-1)) then
                if(4*(n/4) .eq. n) then
                  mini(1)=n/4
                  mini(2)=n/4
                  mini(3)=n/4
                  mini(4)=n/4
                else
                  mini(1)=int(n/4)
                  mini(2)=int(n/4)+1
                  if(4*(n/4) .eq. n-1) then
                    mini(3)=int(n/4)
                    mini(4)=int(n/4)
                  else
                    if(4*(n/4) .eq. n-2) then
                      mini(3)=int(n/4)+1
                      mini(4)=int(n/4)
                    else
                      mini(3)=int(n/4)+1
                      mini(4)=int(n/4)+1
                    endif
                  endif
                endif
              else
                mini(1)=nmini
                mini(2)=nmini
                mini(3)=nmini
                mini(4)=nmini
                mini(5)=nmini
              endif
            endif
          endif
          nhalf=int(mini(1)*percen)
          if(ngroup.gt.kmini) ngroup=kmini
          nrep=int((krep*1.D0)/ngroup)
          minigr=mini(1)+mini(2)+mini(3)+mini(4)+mini(5)
          call rfrdraw(subdat,n,minigr,mini,ngroup,kmini)
        else
          minigr=n
          nhalf=nhalff
          kstep=k1
          if(xilen1.le.replow(nsel)) then
c     		use all combinations; happens iff  nsel = nvar+1 = p+1 <= 6
             nrep=rfncomb(nsel,xilen1)
          else
            nrep=krep
            all=.false.
          endif
        endif

C        CALL INTPR('xilen1',-1,xilen1,1)
C        CALL INTPR('nrep',-1,nrep,1)

cc
cc  Some more initializations:
cc    m1stock = matrix containing the means of the ngroup*10 best estimates
cc		obtained in the subdatasets.
cc    c1stock = matrix containing the covariance matrices of the ngroup*10
cc		best estimates obtained in the subdatasets.
cc    mstock = matrix containing the means of the ten best estimates
cc	       obtained after merging the subdatasets and iterating from
cc	       their best estimates.
cc    cstock = matrix containing the covariance matrices of the ten best
cc	       estimates obtained after merging the subdatasets
cc	       and iterating from their best estimates.
cc    means = mean vector
cc    bmeans = initial MCD location estimate
cc    sd = standard deviation vector
cc    nmahad = vector of mahalanobis distances
cc    ndist = vector of general (possibly robust) distances
cc    inbest = best solution vector
cc    index1 = index vector of subsample observations
cc    index2 = index vector of ordered mahalanobis distances
cc    temp  = auxiliary vector
cc    flag = vector with components indicating the occurrence of a
cc	     singular intermediate MCD estimate.
cc
        do 31 j=1,nvmax
           do 33 k=1,10
              mstock(k,j)=1000000.D0
              do 35 kk=1,kmini
                m1stock((kk-1)*10+k,j)=1000000.D0
35            continue
              do 37 i=1,nvmax
                 do 39 kk=1,kmini
                   c1stock((kk-1)*10+k,(j-1)*nvmax+i)=1000000.D0
39               continue
                 cstock(k,(j-1)*nvmax+i)=1000000.D0
37            continue
33         continue
           means(j)=0.D0
           bmeans(j)=0.D0
           sd(j)=0.D0
31      continue

        do 41 j=1,nmax
           nmahad(j)=0.D0
           ndist(j)=0.D0
           index1(j)=1000000
           index2(j)=1000000
           temp(j)=1000000
41      continue
        do 43 j=1,km10
          flag(j)=1
43      continue


9500    continue
cc
cc ********* Compute the classical estimates **************
cc
c        CALL INTPR('Compute the classical estimates ... ',-1,det,0)
C-EM-MCD-START
c        call rfcovinit(sscp1,nvar+1,nvar+1)
c        do 51 i=1,n
c          do 53 j=1,nvar
c53          rec(j)=dat(i,j)
c          call rfadmit(rec,nvar,nvar+1,sscp1)
c51      continue
c        call rfcovar(n,nvar,nvar+1,sscp1,cova1,means,sd)
C-EM-MCD-END

        do 54 i=1,n
          do 55 j=1,nvar
            xdat(n*(j-1)+i) = dat(i,j)
55        continue
54      continue

        CALL EMNINT(xdat,n,nvar,
     *              ddd,xint,nlen,xdble,ndble,
     *              means,cova1,mvcode)

c        do 57 j=1,nvar
c          if(sd(j).eq.0.D0) goto 5001
c57      continue

        CALL MTXCPY(cova1,cinv1,nvar,nvar)
        CALL MTXINV(cinv1, nvar, det, eps, IERR)
        IF(IERR.NE.0) THEN
          CALL INTPR("SINGULAR CLASSICAL COVARIANCE!",-1,0,0)
          GOTO 5001
        ENDIF

c       if just classical estimate, we are done
        if(class) then
          CALL MTXCPY(means,initmean,nvar,1)
          CALL MTXCPY(cova1,initcov,nvar,nvar)
          goto 9999
        endif

        goto 5002

c  singularity '1' (exact fit := 1) :
5001    continue
C-EM-MCD    EXACT FIT
        fit=1
        goto 9999

5002    continue
cc
cc  Compute and store classical Mahalanobis distances.
cc      --- not used ---
c        CALL NAMAH(dat,n,nvar,means,cova1,cinv1,rec,nmahad,temp,am2,
c     *          temp,am,corr1,mvcode,EPS)
c

cc******** Compute the MCD estimates ************** ------------------------------

cc     Main loop: inspects the subsamples.
cc     Every time the sscp of the subsample is placed in sscp1,
cc     its covariance matrix in cova1, and its inverse in cinv1 .
cc     The minimum covariance determinant matrix is placed in cova2,
cc     and its inverse in cinv2.
cc     The robust distances are placed in ndist.
cc     In the main loop we count the total number of iteration steps
cc     with the variable tottimes.
cc
cc     The algorithm returns here twice when the dataset is divided
cc     at the beginning of the program. According to the situation,
cc     new initializations are made. The second stage, where the subdatasets
cc     are merged, is indicated by the logical value fine and
cc     the last stage, when the whole dataset is considered, by the logical
cc     variable final. In the last stage, the number of iterations nrep
cc     is determined according to the total number of observations
cc     and the dimension.
cc
        tottimes=0
5555    object=10.D25

C        CALL INTPR('MAIN LOOP - NUMBER of TRIALS NREP: ',-1,NREP,1)

        if(.not. part .or. final) then
          nn=n
        else
c	  (part .and. .not. final)
          if(fine) then
            nn=minigr
          endif
        endif

C
C       Determine the number of iterations depending on the size of the
C       data set:
C               nrep = 10 if n <= 5000, otherwise nrep=1
C               kstep = k3 or 10, 9, ..., 1
C
        if(fine .or.(.not.part.and.final)) then
          nrep=10
          nsel=nhalf
          kstep=k2
          if(final) then
            nhalf=nhalff
            ngroup=1
            if(n*nvar.le.100000) then
              kstep=k3
            else
              if(n*nvar .gt. 100000 .and. n*nvar .le. 200000) then
                kstep=10
              else
                if(n*nvar .gt.200000 .and. n*nvar
     *            .le.300000) then
                  kstep=9
                else
                  if(n*nvar .gt.300000 .and. n*nvar
     *               .le.400000) then
                    kstep=8
                  else
                     if (n*nvar .gt.400000 .and. n*nvar
     *                  .le.500000) then
                        kstep=7
            else
              if (n*nvar .gt.500000 .and. n*nvar
     *            .le.600000) then
              kstep=6
            else
            if (n*nvar .gt.600000 .and. n*nvar
     *            .le.700000) then
              kstep=5
            else
              if (n*nvar .gt.700000 .and. n*nvar
     *            .le.800000) then
                kstep=4
              else
                if (n*nvar .gt.800000 .and. n*nvar
     *               .le.900000) then
                  kstep=3
                else
                  if (n*nvar .gt.900000 .and. n*nvar
     *               .le.1000000) then
                kstep=2
                  else
                kstep =1
                  endif
                endif
              endif
            endif
              endif
            endif
                  endif
                endif
              endif
            endif
            if (n.gt.5000) then
              nrep=1
            endif
          else
            nhalf=int(minigr*percen)
          endif
        endif

        do 81 i=1,nsel-1
          index1(i)=i
81      continue
        index1(nsel)=nsel-1

cc
cc  Initialization of the matrices to store partial results. For the
cc  first stage of the algorithm, the currently best covariance matrices and
cc  means are stored in the matrices c1stock and m1stock initialized earlier.
cc  The corresponding objective values and the number of the trial subset
cc  are stored in the matrix mcdndex.
cc  For the second stage of the algorithm or for small datasets, only the
cc  currently best objective values are stored in the same matrix mcdndex
cc  and the corresponding covariance matrices and mean vectors are stored in
cc  the matrices cstock and mstock initialized earlier.
cc
        if(.not. final) then
          do 83 i=1,10
            do 85 j=1,ngroup
              mcdndex(i,1,j)=10.d25
              mcdndex(i,2,j)=10.d25
85          continue
83         continue
        endif

cc
cc  The matrix dath contains the observations to be used in the
cc  algorithm. In the first stage of the split-up procedure dath contains
cc  nmini objects, corresponding to the original observations, with the index
cc  of the processed group in the array subdat. For the second stage, the
cc  data points of all the subdatasets are merged in dath.
cc  The variable kount indicates the occurrence of a singular subsample leading
cc  to the corresponding plane. In some situations the variable kount counts
cc  the number of observations on that plane.
cc
        if(fine .and. .not. final) then
          do 91, j=1,minigr
            do 93, k=1,nvar
               dath(j,k)=dat(subdat(1,j),k)
93          continue
91        continue
        endif
        kount=0

c---- For-Loop over groups  - - - - - - - - - - - - - - - - - - - - -
C        CALL INTPR('LOOP OVER GROUPS - NGROUP ',-1,ngroup,1)
        do 1111 ii= 1,ngroup
          if(.not.fine) kount=0
          if(part .and. .not. fine) nn=mini(ii)
          do 101 i=1,nn
            index2(i)=i
101       continue
          if(part .and. .not. fine) then
            jndex=0
            do 103 j=1,minigr
              if(subdat(2,j).eq.ii) then
                jndex=jndex+1
                subndex(jndex)=subdat(1,j)
              endif
103         continue

            xilen2 = 0
            do 105 j=1,mini(ii)
              complete = .TRUE.
              do 107 k=1,nvar
                IF(dat(subndex(j),k) .EQ. mvcode) complete = .FALSE.
                dath(j,k)=dat(subndex(j),k)
107           continue
                IF(complete) THEN
                  xilen2 = xilen2 + 1
                  xindex2(xilen2) = j
                ENDIF
105         continue
C                CALL INTPR('xilen2',-1,xilen2,1)
C                CALL INTPR('xindex2',-1,xindex2,xilen2)
          endif

cc  The number of trial subsamples is represented by nrep, which depends
cc  on the data situation.
cc  When all (p+1)-subsets out of n can be drawn, the subroutine rfgenpn
cc  is used. Otherwise, random subsamples are drawn by the routine
cc  rfrangen. The trial subsamples are put in the array index1. The
cc  same thing happens for large datasets, except that the number of
cc  observations is nmini instead of n.
cc
cc  When a trial subsample is singular, the algorithm counts the number of
cc  observations that lie on the hyperplane corresponding to this sample.
cc  If, for small datasets, this number is larger than nhalff, the program
cc  stops (exact fit) and gives the mean and the covariance matrix
cc  of the observations on the hyperplane, together with the equation
cc  of the hyperplane.
cc  For large datasets, the algorithm first checks whether there are more
cc  than nhalff observations on the hyperplane. If this is the case, the
cc  program stops for the same reason of exact fit and gives the covariance
cc  matrix and mean of the observations on the hyperplane. If not, the
cc  algorithm counts the number of observations that lie on the hyperplane.
cc  When this number is smaller than the current nhalf in the subdataset, these
cc  observations are extended to nhalf observations by adding those
cc  observations that have smallest orthogonal distances to the hyperplane
cc  and the algorithm continues.
cc  When larger, the coefficients of the hyperplane are stored in the matrix
cc  m1stock for use as starting value in the next stage, and the flag of this
cc  estimate gets the value zero.
cc
cc  In the second stage of the algorithm, when the subdatasets are merged,
cc  the array index2 contains the indices of the observations
cc  corresponding to the nhalf observations with minimal relative distances
cc  with respect to the best estimates of the first stage.
cc  When the estimate of the first stage is a hyperplane, the algorithm
cc  investigates whether there are more than the current nhalf observations of
cc  the merged subdataset on that hyperplane. If so, the coefficients of the
cc  hyperplane are again stored, now in the matrix mstock, for the final
cc  stage of the algorithm.
cc  If not, the observations on the hyperplane are extended to nhalf
cc  observations by adding the observations in the merged dataset with
cc  smallest orthogonal distances to that hyperplane.
cc  For small datasets or for larger datasets with n <= nmini*kmini,
cc  the algorithm already stops when one solution becomes singular,
cc  since we then have an exact fit.
cc
cc  In the third stage, the covariance matrices and means of the best
cc  solutions of the second stage are used as starting values.
cc  Again, when a solution becomes singular, the subroutine 'exact'
cc  determines the hyperplane through at least nhalff observations and stops
cc  because of the exact fit.
cc
cc  When the program stops because of an exact fit, the covariance matrix and
cc  mean of the observations on the hyperplane will always be given.
cc

c        CALL INTPR('NREP: ',-1,nrep,1)
c        CALL INTPR('NSEL: ',-1,nsel,1)
c        CALL INTPR('XILEN1: ',-1,xilen1,1)
c        CALL INTPR('XILEN2: ',-1,xilen2,1)
        do 1000 i=1,nrep
          pnsel=nsel
          tottimes=tottimes+1
          deti=0.D0
          detimin1=0.D0
          step=0

          if((part.and..not.fine).or.(.not.part.and..not.final)) then
            if(part) then
              call rfrangen(xilen2,nsel,index1)
            else
              if(all) then
                call rfgenpn(xilen1,nsel,index1)
              else
                call rfrangen(xilen1,nsel,index1)
              endif
            endif
          endif
cc
cc  The covariance matrix and mean of the initial subsamples are
cc  calculated with the subroutine covar and represented by
cc  the variables cova1 and means.
cc
cc  In the following stages of the algorithm, the covariance matrices and means
cc  used as starting values are already stored in the matrices c1stock
cc  and m1stock (for the second stage), and in the matrices cstock and mstock
cc  (for the third stage).
cc
cc  The inverse cinv1 of the covariance matrix is calculated by the
cc  subroutine rfcovsweep, together with its determinant det.
cc
9550      continue
c        CALL INTPR('I=1,NREP ',-1,i,1)

C VT::: Copy the selected subsample from dat or dath respectively into
C       xdat and than calculate the mean and covariance matrix of xdat
C       Select only complete observations

          if(.not. final .and. .not. fine) then

            do 270 j=1,pnsel
              if(.not.fine.and.part) then
                do 271 m=1,nvar
                 xdat(pnsel*(m-1)+j) = dath(xindex2(index1(j)),m)
271             continue
              elseif(.not.part.and..not.final) THEN
                do 272 m=1,nvar
                 xdat(pnsel*(m-1)+j) = dat(xindex1(index1(j)),m)
272             continue
              else
                CALL INTPR("ERROR - should never come here!",-1,0,0)
              endif
270         continue

c            CALL INTPR('index1',-1,index1,pnsel)
            CALL EMNINT(xdat,pnsel,nvar,
     *                   ddd,xint,nlen,xdble,ndble,
     *                   means,cova1,mvcode)

c            CALL DBLEPR('MEANS',-1,means,nvar)
          elseif(final) then
            if(mstock(i,1).ne.1000000.D0) then
              do 125 jj=1,nvar
                means(jj)=mstock(i,jj)
                do 127 kk=1,nvar
                  cova1((jj-1)*nvar+kk)=cstock(i,(jj-1)*nvar+kk)
127             continue
125           continue
            else
              goto 1111
            endif
          elseif(fine) then
            if(m1stock((ii-1)*10+i,1).ne.1000000.D0) then
              do 131 jj=1,nvar
                means(jj)=m1stock((ii-1)*10+i,jj)
                do 133 kk=1,nvar
                  cova1((jj-1)*nvar+kk)=c1stock((ii-1)*10+i,
     *                  (jj-1)*nvar+kk)
133             continue
131           continue
            else
              goto 1111
            endif
          else
             CALL INTPR("ERROR - should never come here!",-1,0,0)
          endif

C VT::: Invert the cov matrix. If singular add one more observation and
C       retry
C
c          CALL INTPR('I=1,NREP ....... invert now',-1,i,1)
          CALL MTXCPY(cova1,cinv1,nvar,nvar)
          CALL MTXINV(cinv1, nvar, det, eps, IERR)
          IF(IERR.NE.0) THEN
            call rfishsort(index1,pnsel)
            call prdraw(index1,pnsel, nn)
            pnsel=pnsel+1
C            CALL INTPR('REDRAW',-1,index1,pnsel)
            IF(pnsel .GT. nhalf) THEN
               fit = 1
               goto 9999
            ENDIF
            goto 9550
          endif
cc
cc  Mahalanobis distances are computed with the subroutine rfmahad
cc  and stored in the array ndist.
cc  The k-th order statistic of the mahalanobis distances is stored
cc  in dist2. The array index2 containes the indices of the
cc  corresponding observations.
cc
c          CALL INTPR('I=1,NREP ....... compute distances',-1,i,1)
          do 151 j=1,nn
c            CALL INTPR('Distances, j=1,nn',-1,j,1)
            if(.not.part.or.final) then
            do 152 mm=1,nvar
              rec(mm)=dat(j,mm)
152         continue
            else
            do 153 mm=1,nvar
              rec(mm)=dath(j,mm)
153         continue
            endif
c            t=rfmahad(rec,nvar,means,cinv1)
c            ndist(j)=t

             CALL NAMDIST(rec,nvar,means,cova1,cinv1,
     *              ndist(j),temp(j),am2(j),
     *              xint,am,corr1,mvcode,EPS)

             am2(j) = ABS(am2(j))
151       continue
c        CALL DBLEPR('I=1,NREP ....... OK',-1,am2,nn)

          dist2=rffindq(am2,nn,nhalf,index2)
c        CALL INTPR('I=1,NREP ....... OK',-1,i,1)

cc
cc  The variable kstep represents the number of iterations. They depend on
cc  the situation of the program (k1, k2, or k3). Within each
cc  iteration the mean and covariance matrix of nhalf observations are
cc  calculated. The nhalf smallest corresponding mahalanobis distances
cc  determine the subset for the next iteration.
cc  The best subset for the whole data is stored in the array inbest.
cc  The iteration stops when two subsequent determinants become equal.
cc

         do 400 step=1,kstep
            tottimes=tottimes+1
c        CALL INTPR('step=1,kstep',-1,step,1)

C VT::: Compute the covariance matrix of the nhalf observations
C       with smallest D2
            do 155 j=1,nhalf
              temp(j)=index2(j)
155         continue
            call rfishsort(temp,nhalf)

C-EM-MCD-START
c            call rfcovinit(sscp1,nvar+1,nvar+1)
c            do 157 j=1,nhalf
c              if(.not.part.or.final) then
c                do 158 mm=1,nvar
c158               rec(mm)=dat(temp(j),mm)
c              else
c                do 159 mm=1,nvar
c159               rec(mm)=dath(temp(j),mm)
c              endif
c              call rfadmit(rec,nvar,nvar+1,sscp1)
c157         continue
c            call rfcovar(nhalf,nvar,nvar+1,sscp1,cova1,means,sd)
C-EM-MCD-END

            do 157 j=1,nhalf
              if(.not.part.or.final) then
                do 158 mm=1,nvar
                  xdat(nhalf*(mm-1)+j) = dat(temp(j),mm)
158             continue
              else
                do 159 mm=1,nvar
                  xdat(nhalf*(mm-1)+j) = dath(temp(j),mm)
159             continue
              endif
157         continue
            CALL EMNINT(xdat,nhalf,nvar,
     *                   ddd,xint,nlen,xdble,ndble,
     *                   means,cova1,mvcode)

C VT::: Invert the matrix and check for singularity
            CALL MTXCPY(cova1,cinv1,nvar,nvar)
            CALL MTXINV(cinv1, nvar, det, eps, IERR)

c            if(final) then
c            CALL INTPR('TMP BEST',-1,temp,nhalf)
c            CALL DBLEPR('TMP MEANS',-1,means,nvar)
c            CALL DBLEPR('TMP DET',-1,det,1)
c            endif

            IF(IERR.NE.0) THEN
              if(final .or. .not. part .or.
     *          (fine.and..not. final .and. n .le. nmini*kmini)) then
                  fit=1
                  goto 9999
              endif
              if(part.and..not.fine.and.kount.eq.0) then
                fit = 2
                goto 9999
              endif
            ENDIF

C VT::: Break the iterations if the determinant does not decrease
            if(step.ge.2 .and. det.eq.detimin1) then
              goto 5000
            endif

            detimin1=deti
            deti=det

C VT::: Compute the distances, standardize them using wilson-hilferty
C       transformation, take the ABS of them and sort - the result is
C       the first nhalf of index2
            do 171 j=1,nn
              if(.not.part.or.final) then
                do 172 mm=1,nvar
                  rec(mm)=dat(j,mm)
172           continue
              else
                do 173 mm=1,nvar
                  rec(mm)=dath(j,mm)
173             continue
              endif
              CALL NAMDIST(rec,nvar,means,cova1,cinv1,
     *              ndist(j),temp(j),am2(j),
     *              xint,am,corr1,mvcode,EPS)

              am2(j) = ABS(am2(j))
171         continue
            dist2=rffindq(am2,nn,nhalf,index2)

C VT::: Check if the new solution is better and if so, store it
            if(((i.eq.1.and.step.eq.1.and..not.fine)
     *         .or.det.lt.object).and.(final)) then
              medi2=rffindq(ndist,nn,int(n/2),index1)
              object=det
              do 175 jjj=1,nhalf
                inbest(jjj)=index2(jjj)
175           continue
              CALL MTXCPY(cova1,cova2,nvar,nvar)
              CALL MTXCPY(cinv1,cinv2,nvar,nvar)
              CALL MTXCPY(means,bmeans,nvar,1)

C       Best so far solution found - print it
c        CALL INTPR('BEST: ',-1,inbest,nhalf)
c        CALL DBLEPR('BEST MEANS: ',-1,means,nvar)
c        CALL DBLEPR('BEST DET: ',-1,det,1)

            endif
400       continue

cc  After each iteration, it has to be checked whether the new solution
cc  is better than some previous one and therefore needs to be stored. This
cc  isn't necessary in the third stage of the algorithm, where only the best
cc  solution is kept.

5000      if(.not. final) then
            if(part .and. .not. fine) then
              iii=ii
            else
              iii=1
cc	      At the end of the algorithm, only the ten
cc	      best solutions need to be stored.
            endif

cc  For each data group :
cc    If the objective function is lower than the largest value in the
cc    matrix mcdndex :
cc    A distinction is made between different stages of the algorithm:
cc	* At the first stage of the split-up situation:
cc	  -If the new objective value did not yet occur in mcdndex
cc	   its value and corresponding covariance matrix and mean are
cc	   stored at the right place in the matrices mcdndex, c1stock and
cc	   m1stock, and other values are shifted to their new position
cc	   in these arrays.
cc	  -If the new objective value already occurs in mcdndex, a
cc	   comparison is made between the new mean vector and covariance matrix
cc	   and those estimates with the same determinant.
cc	   When for an equal determinant, the mean vector or covariance matrix
cc	   do not correspond, both of them are kept in the matrices mcdndex
cc	   and nbest.
cc	* In the second stage of the algorithm, the covariances and means
cc	  are stored :
cc	  - If the new objective value did not yet occur
cc	    in the matrix mcdndex, it is inserted by shifting the greater
cc	    determinants upwards and doing the same in the arrays mstock
cc	    and cstock.
cc	  - If the new objective value already occurs in the array mcdndex,
cc	    it is compared with all solutions with the same determinant.
cc	    In the case of an equality, the means and covariances
cc	    are compared to determine whether or not to insert the
cc	    new solution.
cc    Otherwise nothing happens. When a singularity occurs,
cc    the determinant in the matrix mcdndex is zero and the
cc    corresponding flag is zero too, so the search in the arrays mcdndex,
cc    m1stock, c1stock, mstock and cstock is done on the rows with flag one.
cc

            if(flag((iii-1)*10+1).eq.1) then
              lll=1
            else
              lll=2
            endif
            do 201 j=lll,10
              if (det .le. mcdndex(j,2,iii)) then
                if(det.ne.mcdndex(j,2,iii)) then
                  if(.not.fine.and.part) goto 203
                  goto 205
                else
                  do 207 kkk=j,10
                    if(det.eq.mcdndex(kkk,2,iii)) then
                      do 209, jjj=1,nvar
                        if(part.and..not.fine) then
                          if(means(jjj).ne.m1stock((iii-1)*10+ kkk,jjj))
     *                                                          goto 203
                        else
                          if(means(jjj).ne.mstock(kkk,jjj))     goto 205
                        endif
209                   continue
                      do 211, jjj=1,nvar*nvar
                        if(part.and..not.fine) then
                          if(cova1(jjj).ne.c1stock((iii-1)*10+ kkk,jjj))
     *                                                          goto 203
                        else
                          if(cova1(jjj).ne.cstock(kkk,jjj))     goto 205
                        endif
211                   continue
                    endif
207               continue
                endif
                goto 1000

203             do 221 k=10,j+1,-1
                  do 223 kk=1,nvar*nvar
                    c1stock((iii-1)*10+k,kk)=
     *              c1stock((iii-1)*10+k-1,kk)
223               continue
                  do 225 kk=1,nvar
                    m1stock((iii-1)*10+k,kk)=
     *              m1stock((iii-1)*10+k-1,kk)
225               continue
                  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
221             continue
                do 227 kk=1,nvar
                  do 229 kkk=1,nvar
                    c1stock((iii-1)*10+j,(kk-1)*nvar+kkk)=
     *              cova1((kk-1)*nvar+kkk)
229               continue
                  m1stock((iii-1)*10+j,kk)=means(kk)
227             continue
                mcdndex(j,1,iii)=i
                mcdndex(j,2,iii)=det
                goto 1000

205             do 231 k=10,j+1,-1
                  do 233 kk=1,nvar*nvar
                    cstock(k,kk) = cstock(k-1,kk)
233               continue
                  do 235 kk=1,nvar
                    mstock(k,kk) = mstock(k-1,kk)
235               continue
                  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
231             continue
                do 237 kk=1,nvar
                  do 239 kkk=1,nvar
                    cstock(j,(kk-1)*nvar+kkk)=
     *                  cova1((kk-1)*nvar+kkk)
239               continue
                  mstock(j,kk)=means(kk)
237             continue
                mcdndex(j,1,iii)=i
                mcdndex(j,2,iii)=det
                goto 1000

              endif
201         continue

          endif
c		(not final)

1000    continue
1111    continue
c---- - - - - - end [ For-loop -- ii= 1, ngroup	 ]  - - - - - - - - -

cc  Determine whether the algorithm needs to be run again or not.
cc
        if(part .and. .not. fine) then
          fine= .true.
          goto 5555
        else if(.not. final .and. ((part.and.fine).or. .not.part)) then
          final= .true.
          goto 5555
        endif

cc******** end { Main Loop } ************** --------------------------------

C VT::: Prepare the output parameters. No need of distances, weights,
C       and adjusted cov matrix
C
C       initmean, initcov, det
C
        CALL MTXCPY(bmeans,initmean,nvar,1)
        CALL MTXCPY(cova2,cova1,nvar,nvar)
        CALL MTXCPY(cova1,initcov,nvar,nvar)
        det=object

c        CALL DBLEPR('FINAL MEANS: ',-1,means,nvar)
c        CALL DBLEPR('FINAL DET: ',-1,det,1)

        goto 9999
cc ******************************************************************

 9999   continue
        call rndend
C            ------ == PutRNGstate() in C
        return
        end
ccccc   end {rffastmcd}


        subroutine prdraw(a, pnsel, nn)
cc
        implicit none
        integer nn, a(nn), pnsel
c
        double precision unifrnd
        integer jndex, nrand, i,j
cc
        jndex=pnsel
        nrand=int(unifrnd() * (nn-jndex))+1
C         if(nrand .gt. nn-jndex) then
C            call intpr(
C      1          '** prdraw(): correcting nrand > nn-jndex; nrand=',
C      2          -1, nrand, 1)
C            nrand=nn-jndex
C         endif

        jndex=jndex+1
        a(jndex)=nrand+jndex-1
        do 5, i=1,jndex-1
           if(a(i).gt.nrand+i-1) then
              do 6,j=jndex,i+1,-1
                 a(j)=a(j-1)
 6            continue
              a(i)=nrand+i-1
              goto 10
c             ------- break
           endif
5       continue
10      continue
        return
        end
ccccc
ccccc
        function rfmahad(rec,nvar,means,sigma)
cc
cc  Computes a Mahalanobis-type distance.
cc
        double precision rec(nvar), means(nvar), sigma(nvar,nvar)
        double precision rfmahad, t
cc
        t=0
        do 100 j=1,nvar
          do 90 k=1,nvar
            t=t+(rec(j)-means(j))*(rec(k)-means(k))*sigma(j,k)
90        continue
100     continue
        rfmahad=t
        return
        end
ccccc
ccccc

      subroutine rfstore1(nvar,c1stock,m1stock,nvmax2,nvmax,
     * kmini,cova1,means,i,km10,ii,mcdndex,kount)
cc
      double precision c1stock(km10,nvmax2)
      double precision m1stock(km10,nvmax)
      double precision mcdndex(10,2,kmini)
      double precision cova1(nvar,nvar)
      double precision means(nvar)
cc
      do 10,k=10,2,-1
      do 20 kk=1,nvar*nvar
       c1stock((ii-1)*10+k,kk)=
     * c1stock((ii-1)*10+k-1,kk)
20    continue
      do 30 kk=1,nvar
       m1stock((ii-1)*10+k,kk)=
     *   m1stock((ii-1)*10+k-1,kk)
30    continue
      mcdndex(k,1,ii)=mcdndex(k-1,1,ii)
      mcdndex(k,2,ii)=mcdndex(k-1,2,ii)
10    continue
      do 40 kk=1,nvar
      m1stock((ii-1)*10+1,kk)=means(kk)
      do 50 jj=1,nvar
      c1stock((ii-1)*10+1,(kk-1)*nvar+jj)=
     *    cova1(kk,jj)
50    continue
40    continue
      mcdndex(1,1,ii)=i
      mcdndex(1,2,ii)=kount
      return
      end

ccccc
ccccc
c
c--   Routines common to
c--   fastLTS ( ./rfltsreg.f )  and
c--   fastMCD ( ./rffastmcd.f )
c
c
      subroutine rfrangen(n, nsel, index)
c
c     Randomly draws nsel cases out of n cases.
c     Here, index is the index set.
c
      implicit none
      integer n, nsel, index(nsel)
c     real uniran
      double precision unifrnd
      integer i,j, num
c
      do 100 i=1,nsel
10      num=int(unifrnd()*n)+1
C     if(num .gt. n) then
C       call intpr('** rfrangen(): num > n; num=', -1, num, 1)
C       num=n
C     endif
         if(i.gt.1) then
            do 50 j=1,i-1
               if(index(j).eq.num) goto 10
50          continue
         endif
         index(i)=num
100   continue
      return
      end

      subroutine rfgenpn(n,nsel,index)
cc
cc    Constructs all subsets of nsel cases out of n cases.
cc
        implicit none
        integer n,nsel,index(nsel)
cc
        integer k,i

      k=nsel
      index(k)=index(k)+1
c while
10    if(k.eq.1 .or. index(k).le.(n-(nsel-k))) goto 100
      k=k-1
      index(k)=index(k)+1
      do 50 i=k+1,nsel
       index(i)=index(i-1)+1
50    continue
      goto 10
c end{while}
100   return
      end
ccccc
ccccc
      subroutine rfshsort(a,n)
cc
cc    Sorts the array a of length n.
cc
        implicit none
        integer n
      double precision a(n)
c
      double precision t
      integer gap, i,j, nextj

      gap=n
 100    gap=gap/2
      if(gap.eq.0) goto 200
      do 180 i=1,n-gap
         j=i
 120     if(j.lt.1) goto 180
      nextj=j+gap
      if(a(j).gt.a(nextj)) then
        t=a(j)
        a(j)=a(nextj)
        a(nextj)=t
      else
        j=0
      endif
      j=j-gap
      goto 120
 180    continue
      goto 100
 200    return
      end
ccccc
ccccc
        subroutine rfishsort(a,kk)
cc
cc  Sorts the integer array a of length kk.
cc
        implicit none
        integer kk, a(kk)
c
        integer t, gap, i,j, nextj

        gap=kk
100     gap=gap/2
        if(gap.eq.0) goto 200
        do 180 i=1,kk-gap
          j=i
120       if(j.lt.1) goto 180
      nextj=j+gap
      if(a(j).gt.a(nextj)) then
        t=a(j)
        a(j)=a(nextj)
        a(nextj)=t
      else
        j=0
      endif
      j=j-gap
      goto 120
 180    continue
      goto 100
 200   return
      end
ccccc
ccccc
      function replow(k)
cc
cc    Find out which combinations of n and p are
cc    small enough in order to perform exaustive search
cc    Returns the maximal n for a given p, for which
cc    exhaustive search is to be done
cc
cc    k is the number of variables (p)
cc
      implicit none
      integer replow, k
c
      integer irep(6)
      data irep/500,50,22,17,15,14/
c
      if(k .le. 6) then
        replow = irep(k)
      else
        replow = 0
      endif
      return
      end
ccccc
ccccc
      function rfncomb(k,n)
cc
cc  Computes the number of combinations of k out of n.
cc  (To avoid integer overflow during the computation,
cc  ratios of reals are multiplied sequentially.)
cc  For comb > 1E+009 the resulting 'comb' may be too large
cc  to be put in the integer 'rfncomb', but the main program
cc  only calls this function for small enough n and k.
cc
      implicit none
      integer rfncomb,k,n
c
      double precision comb,fact
      integer j
c
      comb=dble(1.0)
      do 10 j=1,k
         fact=(dble(n-j+1.0))/(dble(k-j+1.0))
         comb=comb*fact
 10   continue
      rfncomb=int(comb+0.5D0)
      return
      end
ccccc
ccccc
        function rffindq(aw,ncas,k,index)
cc
cc  Finds the k-th order statistic of the array aw of length ncas.
cc
cc MM{FIXME}: "rather" use R's C API   rPsort (double* X, int N, int K)

        implicit none
        integer ncas,k,index(ncas)
        double precision rffindq, aw(ncas)
c
        double precision ax,wa
        integer i,j,l,lr,jnc
cc
        do 10 j=1,ncas
          index(j)=j
10      continue
        l=1
        lr=ncas
c    while(l < lr)
20      if(l.ge.lr) goto 90
        ax=aw(k)
        jnc=l
        j=lr

30      if(jnc.gt.j) goto 80
40      if(aw(jnc).ge.ax) goto 50
        jnc=jnc+1
        goto 40
50      if(aw(j).le.ax) goto 60
        j=j-1
        goto 50
60      if(jnc .le. j) then
           i=index(jnc)
           index(jnc)=index(j)
           index(j)=i
           wa=aw(jnc)
           aw(jnc)=aw(j)
           aw(j)=wa
           jnc=jnc+1
           j=j-1
        endif
        goto 30
80      if(j.lt.k) l=jnc
        if(k.lt.jnc) lr=j
        goto 20
90      rffindq=aw(k)
        return
        end
ccccc
ccccc
        subroutine rfrdraw(a,n,ntot,mini,ngroup,kmini)
cc
cc  Draws ngroup nonoverlapping subdatasets out of a dataset of size n,
cc  such that the selected case numbers are uniformly distributed from 1 to n.
cc
        implicit none
        integer ntot, kmini, a(2,ntot), n, mini(kmini), ngroup
c
        double precision unifrnd
c
        integer jndex, nrand, k,m,i,j
cc
        jndex=0
        do 10 k=1,ngroup
           do 20 m=1,mini(k)
            nrand=int(unifrnd()*(n-jndex))+1
C       if(nrand .gt. n-jndex) then
C          call intpr(
C      1         '** rfrdraw(): need to correct nrand > n-jndex; nrand=',
C      2	                  -1, nrand, 1)
C          nrand=n-jndex
C        endif

            jndex=jndex+1
            if(jndex.eq.1) then
              a(1,jndex)=nrand
              a(2,jndex)=k
            else
               a(1,jndex)=nrand+jndex-1
               a(2,jndex)=k
               do 5,i=1,jndex-1
                  if(a(1,i).gt.nrand+i-1) then
                     do 6, j=jndex,i+1,-1
                        a(1,j)=a(1,j-1)
                        a(2,j)=a(2,j-1)
6                    continue
                     a(1,i)=nrand+i-1
                     a(2,i)=k
                     goto 20
c                    ------- break
                 endif
5              continue
            endif
20         continue
10      continue
        return
        end
ccccc
ccccc
        function rfodd(n)
        logical rfodd

        rfodd=.true.
        if(2*(n/2).eq.n) rfodd=.false.
        return
        end
