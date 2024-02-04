C       EMNINT  - an interface to EMNCOVC providing linear
C               storage for easier call
C
C       The essenetial input parameter is x(n,p) and mvcode and
C       the output ones are mu and sigma.
C       The work arrays are xint(nint) and xdble(ndble)
C       Additinally 'd' is passed which is the dimension of the packed
C       storage for EMN()
C
C       d = d=(2+3*p+p^2)/2
C       nint = (p+1)*(p+1) + n*p + 3*n + 3*p
C       ndble = 4*d + 3*p + n*p
C
C       Below is given the mapping of the vectors.
C
C nint: (p+1)*(p+1) + n*p + 3*n + 3*p
C xint: psi           r     mdpst nmdp ro  nmis oc  mc
C       (p+1)*(p+1) + n*p + n +   n +  n + p +  p + p
C
C ndble: 4*d + 3*p + n*p
C xdble: old  theta  t   tobs  c    sdv  xbar  xtmp
C        d  + d    + d + d +   p +  p +  p  +  n*p
C
        SUBROUTINE EMNINT(x, n, p,
     1        d, xint,nlen,xdble, ndble,
     2        mu, sigma,
     3        mvcode)

        integer n,p,d,nlen,ndble
        integer xint(nlen)
        double precision x(n,p), mvcode
        double precision xdble(ndble)
        double precision mu(p), sigma(p,p)

        integer ipsi, ir, imdpst, inmdp, iro, inmis, ioc, imc
        integer iold, itheta, it, itobs, ic, isdv, ixbar,ixtmp
        integer inext, idnext

        ipsi   = 1
        ir     = ipsi + (p+1)*(p+1)
        imdpst = ir + n*p
        inmdp  = imdpst + n
        iro    = inmdp + n
        inmis  = iro + n
        ioc    = inmis + p
        imc    = ioc + p
        inext  = imc + p

        iold   = 1
        itheta = iold + d
        it     = itheta + d
        itobs  = it + d
        ic     = itobs + d
        isdv   = ic + p
        ixbar  = isdv + p
        ixtmp  = ixbar + p
        idnext = ixtmp + n*p

c        CALL INTPR("N:",-1,n,1)
c        CALL INTPR("P:",-1,p,1)
c        CALL INTPR("D:",-1,d,1)

c        CALL INTPR("ndble:",-1,ndble,1)
c        CALL INTPR("NLEN:",-1,nlen,1)

c        CALL INTPR(">>>>ENTER EMNCOV",-1,0,0)
        CALL EMNCOV(x, n, p,
     *        d, xint(ipsi),
     *        xdble(iold), xdble(itheta), xdble(it), xdble(itobs),
     1        xint(ir),xint(imdpst),xint(inmdp),xint(inmis),xint(iro),
     2        xint(ioc),xint(imc),
     3        xdble(ic), xdble(isdv), xdble(ixbar),
     4        mu, sigma,
     5        mvcode,
     6        xdble(ixtmp))

c       CALL INTPR("<<<<EXIT EMNCOV",-1,0,0)

        RETURN
        END

C       EMNCOV - calculates the EM mean and covariance of x(n,p).
C       At each EM step emn() from package norm is called. For each EM
C       step theta must be in sweep(0) condition. After execution,
C       the new parameter value is contained in t, and theta is left
C       swept on columns corresponding to observed variables in the
C       first missingness pattern.

C       x(n,p)  - input matrix
C       d       - integer, dimension of the packed storage,
C                 d=(2+3*p+p^2)/2
C       psi     integer (p+1,p+1), packed storage matrix
C
C       old, theta, t, tobs - double(d) - work arrays. The final result
C               as packed storage is in theta and can be retrieved by
C               calling GETPAR() - this is done authomatically at the
C               end and the par[mu, sigma] is returned too
C
C       r       - integer(n,p)
C       ndpst, nmdp, ro - integer(n)
C       nmis    - integer(p) - unmber of missing values in each variable
C       oc, mc  - integer(p), work
C       c       - double(p), work
C       sdv,xbar- double(p), work
C       mu      - double(p) OUTPUT
C       sigma   - double(p,p) OUTPUT
C       mvcode  - double, missingness code
C       xtmp    - double(n,p), work
C
        SUBROUTINE EMNCOV(x, n, p,
     *        d, psi,
     *        old, theta, t, tobs,
     1        r, mdpst, nmdp, nmis, ro,
     2        oc,mc,
     3        c, sdv, xbar,
     4        mu, sigma,
     5        mvcode,
     6        xtmp)

        integer n,p,d,psi(0:p,0:p),r(n,p),nmdp(n)
        integer mdpst(n), ro(n), nmis(p),oc(p), mc(p)
        double precision x(n,p), mvcode, xtmp(n,p)
        double precision old(d), theta(d), t(d), tobs(d)
        double precision mu(p), sigma(p,p)
        double precision c(p)
        double precision sdv(p), xbar(p)

        integer mle, npatt
        double precision tau,m
        double precision tmax, delta, criterion

C       actually these should be parameters
        maxits = 1000
        criterion = 1e-04

C       not used - no prior is specified, finds the mle
        mle = 1
        tau = 1
        m   = 1

C Form matrix of packed storage indices
C        d <- as.integer((2 + 3 * p + p^2)/2)
C        psi <- .Fortran("mkpsi", p, matrix(as.integer(0), p + 1, p + 1))[[2]]

        d = int(2 + 3 * p + p*p)/2
        CALL MKPSI(p, psi)

C Find missingness patterns
C    r <- 1 * is.na(x)
C    nmis <- as.integer(apply(r, 2, sum))

        DO 10 j=1,p
           nmis(j) = 0
10      CONTINUE

        DO 20 i=1,n
           DO 30 j=1,p
              r(i,j) = 0
              IF(x(i,j) .EQ. mvcode) THEN
                 r(i,j) = 1
                 nmis(j) = nmis(j) + 1
              ENDIF
30         CONTINUE
20      CONTINUE

C Index the missing data patterns
C   mdp <- as.integer((r %*% (2^((1:ncol(x)) - 1))) + 1)

        DO 40 j=1,p
           oc(j) = 2 ** (j - 1 )
40      CONTINUE
C        CALL INTPR('OC: ',-1,oc,p)

        DO 50 i=1,n
           nmdp(i) = 0
           DO 60 j=1,p
              nmdp(i) = nmdp(i) + r(i,j)*oc(j)
60         CONTINUE
           nmdp(i) = nmdp(i) + 1
50      CONTINUE

C Do row sort and rearrange the rows of x, r, nmdp
        CALL MYORD(nmdp, n, ro)
        DO 301 i=1,n
           do 300 j=1,p
              xtmp(i,j) = x(ro(i),j)
300        CONTINUE
301    CONTINUE
        DO 311 i=1,n
           do 310 j=1,p
              x(i,j) = xtmp(i,j)
310        CONTINUE
311    CONTINUE

        DO 341 i=1,n
           do 340 j=1,p
              xtmp(i,j) = r(ro(i),j)
340        CONTINUE
341    CONTINUE

        DO 351 i=1,n
           Do 350 j=1,p
              r(i,j) = xtmp(i,j)
350        CONTINUE
351    CONTINUE

        DO 360 i=1,n
           xtmp(i,1) = nmdp(ro(i))
360     CONTINUE
        DO 370 i=1,n
           nmdp(i) = xtmp(i,1)
370     CONTINUE


C Compress missing data patterns - only the non-duplicates
C       npatt is the number of non-dupl. patterns which
C       are stored in mdpst(npatt)
        CALL MYNDUPL(nmdp, n, mdpst, npatt)

C Compress the matrix r(n,p) -> r(npatt,p)
        DO 381 i=1,npatt
           do 380 j=1,p
              xtmp(i,j) = 1 - r(mdpst(i),j)
380        CONTINUE
381     CONTINUE
        CALL SETMAT(r, n, p, npatt, xtmp)

C Other bookkeeping quantities - nmdp will contain
C       the number of observations for each unique pattern
        IF(npatt .EQ. 1) THEN
           nmdp(1) = n
        ELSE
           DO 390 i=2,npatt
              nmdp(i-1) = mdpst(i) - mdpst(i-1)
390        CONTINUE
           nmdp(npatt) = n + 1 - mdpst(npatt)
        ENDIF
C        CALL INTPR('NMDP: ',-1,nmdp,npatt)

C Center and scale the columns of x
        CALL CTRSC(x, n, p, xbar, sdv, mvcode)
c        CALL DBLEPR("XBAR:",-1,xbar,p)
C Calculate initial values
        CALL STVALN(d, theta, p, psi)
        CALL TOBSN(d, tobs, p, psi, n, x, npatt, r, mdpst, nmdp, oc)

C Main loop call EMN from package norm
        DO 100 it=1,maxits
           DO 80 j=1,d
              t(j) = theta(j)
              old(j) = theta(j)
80         CONTINUE
           CALL EMN(d, t, theta, tobs, p, psi, n, x,
     1           npatt, r, mdpst, nmdp, oc, mc, c, mle, tau, m,
     2           mu, sigma)

c          CALL DBLEPR('old: ',-1,old,d)
c          CALL DBLEPR('theta: ',-1,theta,d)
          tmax = -1
          DO 90 J=1,d
                delta = abs(old(j)-theta(j))
                IF(delta .GT. tmax) THEN
                   tmax = delta
                ENDIF
90         continue
c           CALL DBLEPR('tmax: ',-1,tmax,1)
           if(tmax .LE. criterion) goto 9999
100     continue
9999    continue

        CALL GETPAR(p, d, theta, psi, sdv, xbar, mu, sigma)
        RETURN
        END

        SUBROUTINE SETMAT(r, n, p, npatt, xtmp)
        integer npatt, p, r(npatt,p)
        double precision xtmp(n,p)

        DO 391 i=1,npatt
           Do 390 j=1,p
              r(i,j) = xtmp(i,j)
390        CONTINUE
391     CONTINUE

C       CALL INTPR('R: ',-1,r,npatt*p)
        RETURN
        END
C
C Returns mu=means and sigma=covariance from theta (packed storage,
C   centered and standardised with xbar, sdev)
C
        SUBROUTINE GETPAR(p, d, theta, psi, sdv, xbar, mu, sigma)
        integer p, d, psi(0:p,0:p)
        double precision theta(d), sdv(p), xbar(p)
        double precision mu(p), sigma(p,p)

C       mu <- theta[s$psi[1, 2:(s$p + 1)]] * s$sdv + s$xbar
        do 10 i=1,p
           mu(i) = theta(psi(0, i)) * sdv(i) + xbar(i)
10      CONTINUE

C       sigma <- theta[s$psi[2:(s$p + 1), 2:(s$p + 1)]]
C       sigma <- matrix(sigma, s$p, s$p)
C       tmp <- matrix(s$sdv, s$p, s$p)
C       sigma <- sigma * tmp * t(tmp)

        do 30 i=1,p
           do 20 j=1,p
              sigma(i,j) = theta(psi(i, j)) * sdv(i) * sdv(j)
20         continue
30      continue
        return
        end

C***********************************************************************
        subroutine stvaln(d,theta,p,psi)
C Gets starting value of theta: mu=0 and sigma=I
        integer d,p,psi(0:p,0:p)
        double precision theta(d)
        call initn(d,theta)
        theta(1)=-1.
        do 5 j=1,p
           theta(psi(j,j))=1.
5       continue
        return
        end
C***********************************************************************
        subroutine tobsn(d,tobs,p,psi,n,x,npatt,r,mdpst,nmdp,oc)
C Tabulates the known part of the sscp matrix for all missingness
C patterns.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p)
        double precision tobs(d)

        call initn(d,tobs)
        do 40 patt=1,npatt
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 30 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 20 j=1,noc
                 tobs(psi(0,oc(j)))=tobs(psi(0,oc(j)))+x(i,oc(j))
                 do 10 k=j,noc
                    tobs(psi(oc(j),oc(k)))=tobs(psi(oc(j),oc(k)))+
     /                   x(i,oc(j))*x(i,oc(k))
10               continue
20            continue
30         continue
40      continue
        return
        end
C************************************************************************
        subroutine emn(d,theta,t,tobs,p,psi,n,x,npatt,r,mdpst,nmdp,
     /      oc,mc,c,mle,tau,m,mu,lmbinv)
C Performs one step of em. Theta must be in sweep(0) condition.
C After execution, the new parameter value is contained in t, and
C theta is left swept on columns corresponding to observed variables
C in the first missingness pattern.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,mc(p),nmc,patt,mle
        double precision x(n,p)
        double precision theta(d),t(d),tobs(d),c(p),tau,m,mu(p)
        double precision lmbinv(p,p)
        do 1 i=1,d
           t(i)=tobs(i)
1       continue
        do 200 patt=1,npatt
           call swpobs(d,theta,p,psi,npatt,r,patt)
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 150 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 50 j=1,nmc
                 c(mc(j))=theta(psi(0,mc(j)))
                 do 40 k=1,noc
                    c(mc(j))=c(mc(j))+theta(psi(oc(k),mc(j)))*x(i,oc(k))
40               continue
50            continue
              do 100 j=1,nmc
                 t(psi(0,mc(j)))=t(psi(0,mc(j)))+c(mc(j))
                 do 70 k=1,noc
                    t(psi(oc(k),mc(j)))=t(psi(oc(k),mc(j)))+
     /                   x(i,oc(k))*c(mc(j))
70               continue
                 do 80 k=j,nmc
                    t(psi(mc(k),mc(j)))=t(psi(mc(k),mc(j)))+
     /                c(mc(k))*c(mc(j))+theta(psi(mc(k),mc(j)))
80               continue
100           continue
150        continue
200     continue
        if(mle.eq.0) then
           call moden(d,t,p,psi,n,tau,m,mu,lmbinv)
        endif
        do 210 i=2,d
           t(i)=t(i)/dble(n)
210     continue
        call swp(d,t,0,p,psi,p,1)
        return
        end
C************************************************************************
        subroutine swpobs(d,theta,p,psi,npatt,r,patt)
C Sweeps theta to condition on the observed variables
        integer d,p,psi(0:p,0:p),npatt,r(npatt,p),patt
        double precision theta(d)
        do 10 j=1,p
           if((r(patt,j).eq.1).and.(theta(psi(j,j)).gt.0.))then
              call swp(d,theta,j,p,psi,p,1)
           elseif((r(patt,j).eq.0).and.(theta(psi(j,j)).lt.0.))then
              call swp(d,theta,j,p,psi,p,-1)
           endif
10      continue
        return
        end

C***********************************************************************
        subroutine swp(d,theta,pivot,p,psi,submat,dir)
C Performs sweep on a symmetric matrix in packed storage.
C Sweeps on pivot position. Sweeps only the (0:submat,0:submat)
C submatrix.
C If dir=1, performs ordinary sweep. If dir=-1, performs reverse sweep.
        integer d,p,pivot,psi(0:p,0:p),submat,dir
        double precision theta(d),a,b,c
        a=theta(psi(pivot,pivot))
        theta(psi(pivot,pivot))=-1./a
        do 10 j=0,submat
           if(j.ne.pivot) theta(psi(j,pivot))=theta(psi(j,pivot))/a*dir
10      continue
        do 30 i=0,submat
           do 20 j=i,submat
              if((i.ne.pivot).and.(j.ne.pivot))then
                 b=theta(psi(i,pivot))
                 c=theta(psi(j,pivot))
                 theta(psi(i,j))=theta(psi(i,j))-a*b*c
              endif
20         continue
30      continue
        return
        end
C***********************************************************************
        subroutine initn(d,theta)
C Initializes theta
        integer d
        double precision theta(d)
        theta(1)=1.
        do 1 i=2,d
           theta(i)=0.
1       continue
        return
        end
C************************************************************************
        subroutine gtmc(p,npatt,r,patt,mc,nmc,last)
C Finds the column numbers of the missing variables, and stores them
C in the first nmc elements of mc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,mc(p),nmc,last
        nmc=0
        do 10 j=1,last
           if(r(patt,j).eq.0)then
              nmc=nmc+1
              mc(nmc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtoc(p,npatt,r,patt,oc,noc,last)
C Finds the column numbers of the observed variables, and stores them
C in the first noc elements of oc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,oc(p),noc,last
        noc=0
        do 10 j=1,last
           if(r(patt,j).eq.1)then
              noc=noc+1
              oc(noc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine moden(d,t,p,psi,n,tau,m,mu,lmbinv)
C Alters the sufficient statistics to yield a posterior mode.
        integer d,p,psi(0:p,0:p),n
        double precision t(d),tau,m,mu(p),lmbinv(p,p),c,e
        do 5 j=1,p
           mu(j)=mu(j)*dble(n)
5       continue
        c=tau/(dble(n)*(tau+dble(n)))
        e=dble(n)/(dble(n)+m+dble(p)+2.)
        do 20 j=1,p
           do 10 k=j,p
              t(psi(j,k))=t(psi(j,k))+lmbinv(j,k)-(t(psi(0,j))*
     /           t(psi(0,k)))/dble(n)
              t(psi(j,k))=t(psi(j,k))+c*(t(psi(0,j))-mu(j))*
     /           (t(psi(0,k))-mu(k))
              t(psi(j,k))=t(psi(j,k))*e
10         continue
20      continue
        c=dble(n)/(tau+dble(n))
        e=1.-c
        do 30 j=1,p
           t(psi(0,j))=c*t(psi(0,j))+e*mu(j)
30      continue
        do 50 j=1,p
           do 40 k=j,p
              t(psi(j,k))=t(psi(j,k))+(t(psi(0,j))*
     /           t(psi(0,k)))/dble(n)
40         continue
50      continue
        return
        end

C***********************************************************************
        subroutine mkpsi(p, psi)
C Generates a symmetric matrix of integers indicating the linear
C position in packed storage of the matrix elements
        integer p, psi(0:p,0:p), posn
        posn=0
        do 10 j=0,p
           posn=posn+1
           psi(j,j)=posn
           do 5 k=j+1,p
              posn=posn+1
              psi(j,k)=posn
              psi(k,j)=posn
5          continue
10      continue
        return
        end

C***********************************************************************
        subroutine ctrsc(x,n,p,xbar,sdv,mvcode)
C Centers and scales the data matrix so that the observed data in every
C column have mean zero and variance one. If a column has zero variance
C or less than 2 observations then the data are centered (set equal
C to zero)
        integer n,p,count
C real x(n,p),mvcode was replaced by line below
        double precision x(n,p),mvcode
        double precision xbar(p),sdv(p),sum1,sum2
        do 10 j=1,p
           sum1=0
           sum2=0
           count=0
           do 5 i=1,n
              if(x(i,j).ne.mvcode) then
                 count=count+1
                 sum1=sum1+x(i,j)
                 sum2=sum2+x(i,j)**2.
              endif
5          continue
           if(count.gt.0) then
              xbar(j)=sum1/count
              sdv(j)=sqrt((sum2-(sum1**2.)/count)/count)
              do 7 i=1,n
                 if(x(i,j).ne.mvcode) x(i,j)=x(i,j)-xbar(j)
7             continue
              if(sdv(j).gt.0.) then
                 do 9 i=1,n
                    if(x(i,j).ne.mvcode) x(i,j)=x(i,j)/sdv(j)
9                continue
               else
                 sdv(j)=1.
              endif
           else
              sdv(j)=1.
           endif
10      continue
        return
        end

        subroutine mysort(a, kk)
        implicit none
        integer kk, a(kk)

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
180     continue
        goto 100
200     return
        end

        SUBROUTINE MYORD(a, kk, ind)
        implicit none
        integer kk, a(kk), ind(kk)

        integer t, gap, i,j, nextj

        do 10 i=1,kk
           ind(i) = i
10      CONTINUE
        gap=kk
100     gap=gap/2
        if(gap.eq.0) goto 200
        do 180 i=1,kk-gap
          j=i
120       if(j.lt.1) goto 180
          nextj=j+gap
          if(a(ind(j)) .GT. a(ind(nextj))) then
            t=ind(j)
            ind(j)=ind(nextj)
            ind(nextj)=t
          else
            j=0
          endif
          j=j-gap
          goto 120
180     continue
        goto 100
200     return
        end
C
C Return in ret the nret non-duplicated entries of the
C       integer array a(kk)
C
        SUBROUTINE MYNDUPL(a, kk, ret, nret)
        implicit none
        integer kk, a(kk), ret(kk), nret
        integer ndupl, i, j

        nret = 0
        do 10 i=1,kk
           ndupl = 1
           if(i .NE. 1) THEN
              do 20 j=1,i-1
                 if(a(j) .EQ. a(i)) THEN
                    ndupl = 0
                    goto 30
                 ENDIF
20            continue
           ENDIF
30         continue
           if(ndupl.eq.1) THEN
              nret = nret + 1
              ret(nret) = i
           ENDIF
10      continue
        RETURN
        END

C
C Computes a Mahalanobis-type distance of the observation in rec
C       relative to the mean 'mu' and the cov matrix (inverted)
C       'cinv'
C
        FUNCTION MDIST(rec, p, mu, cinv)
        INTEGER p
        DOUBLE PRECISION rec(p), mu(p), cinv(p,p)
        DOUBLE PRECISION MDIST
        DOUBLE PRECISION t

c        call dblepr("rec",-1,rec,p)
c        call dblepr("mu",-1,mu,p)
c        call dblepr("ctmp",-1,cinv,p*p)

        t=0
        do 100 j=1,p
          do 90 k=1,p
            t=t+(rec(j)-mu(j))*(rec(k)-mu(k))*cinv(j,k)
90        continue
100     continue
        MDIST=t
        RETURN
        END

C Returns an array of Mahalanobis-type distance of the  possibly
C       incomplete data observations in x(n,p) relative to the
C       mean 'mu' and the cov matrix 'sigma'
C
C       INPUT:
C       x(n,p) - data matrix
C       mu
C       sigma
C       cinv - the inverse of sigma, to be used in case of complete
C               observations
C       OUTPUT:
C       amah - array of mahalanobis distances
C       anov - array of degrees of freedom for chisq (dim of
C              observed values in each observation)
C       az   - array of standardized distances using Wilson-Hilferty
C               transformation
C
        SUBROUTINE NAMAH(x,n,p,mu,sigma,cinv,rec,amah,anov,az,
     *          ov,mutmp,ctmp,mvcode,EPS)
        INTEGER n, p, ov(p), anov(n)
        INTEGER nov
        DOUBLE PRECISION x(n,p), amah(n), az(n)
        DOUBLE PRECISION rec(p), mu(p), sigma(p,p), cinv(p,p)
        DOUBLE PRECISION mutmp(p), ctmp(p*p)
        DOUBLE PRECISION mvcode, EPS, mah, z

        do 10 i=1,n
           do 20 j=1,p
              rec(j) = x(i,j)
20         CONTINUE
           CALL NAMDIST(rec,p,mu,sigma,cinv,mah,nov,z,
     *                     ov,mutmp,ctmp,mvcode,EPS)

              amah(i) = mah
              az(i) = z
              anov(i) = nov
10        continue

        RETURN
        END
C
C Computes a Mahalanobis-type distance of the incomplete
C       observation in rec relative to the mean 'mu'
C       and the cov matrix 'sigma'
C
C       INPUT:
C       rec - data observation (possibly incomplete)
C       mu
C       sigma
C       cinv - the inverse of sigma, to be used in case of complete
C               observation
C       OUTPUT:
C       mah - mahalanobis distance
C       nov - degrees of freedom for chisq (dim of observed values in
C               rec)
C       z  - standardized distances using Wilson-Hilferty
C               transformation
C
        SUBROUTINE NAMDIST(rec,p,mu,sigma,cinv,mah,nov,z,
     *                     ov,mutmp,ctmp,mvcode,EPS)
        INTEGER p, ov(p), nov
        DOUBLE PRECISION rec(p), mu(p), sigma(p,p), cinv(p,p)
        DOUBLE PRECISION mutmp(p), ctmp(p*p)
        DOUBLE PRECISION mah, z, mvcode
        DOUBLE PRECISION MDIST, DET, EPS

C Collect the indexes of the observed values in 'ov'
        mah = 0D0
        z = 0D0
        nov=0
        DO 10 i=1,p
           if(rec(i) .NE. mvcode) THEN
              nov = nov + 1
              ov(nov) = i
           ENDIF
10      CONTINUE

        IF(nov .EQ. 0) THEN
C All items are missing - return 0.
           mah = 0.
           GOTO 9000
        ELSEIF(nov .EQ. p) THEN
C Complete observation - copy mu and cinv in mutmp and ctmp
          CALL MTXCPY(mu, mutmp, p, 1)
          CALL MTXCPY(cinv, ctmp, p, p)
        ELSE
C Compress the data, mu and sigma - keep only observded values
C       and invert sigma
          DO 20 i=1,nov
             rec(i) = rec(ov(i))
             mutmp(i) = mu(ov(i))
             DO 30 j=1,nov
               ctmp((i-1)*nov+j) = sigma(ov(i), ov(j))
               ctmp((j-1)*nov+i) = sigma(ov(j), ov(i))
30           CONTINUE
20        CONTINUE
          CALL MTXINV(CTMP, nov, DET, EPS, IER)
          IF(IER .NE. 0) THEN
            CALL DBLEPR("ERROR INVERTING COV",-1,rec,nov)
            GOTO 9000
          ENDIF
        ENDIF

        mah=MDIST(rec, nov, mutmp, ctmp)
        z = (mah/nov) ** (1D0/3D0)
c        z = (mah/nov) ** (1D0/3D0) - 1D0 + 2D0/(9D0*nov)
c        z = z / SQRT(2D0/(9D0*nov))

9000    RETURN
        END

        SUBROUTINE MTXCPY(A,B,N,M)
        DOUBLE PRECISION A(N,M), B(N,M)
        DO 100 I=1,N
           DO 90 J=1,M
              B(I,J) = A(I,J)
90         CONTINUE
100     CONTINUE
        RETURN
        END

        SUBROUTINE MTXINV(A,P,DET,EPS,IER)
        INTEGER P,IER
        DOUBLE PRECISION pivot, A(P*P), DET, EPS

        IER = 0
        DET = 1.D0

        DO 10 J=1,P
          pivot = A((J-1)*P + J)
          det = det * pivot
          IF(pivot.lt.eps) THEN
             IER = J
             GOTO 999
          ENDIF
          CALL MTXSWP(A,P,J)
10      CONTINUE
c        CALL DBLEPR('DET: ',-1,det,1)
c        CALL INTPR('IER: ',-1,IER,1)
        RETURN
999     CONTINUE
c        CALL DBLEPR('DET: ',-1,det,1)
c        CALL INTPR('IER: ',-1,IER,1)
        RETURN
        END

        SUBROUTINE MTXSWP(A, P, K)
        INTEGER P,K
        DOUBLE PRECISION A(P,P), B, D

        D = A(K,K)

        DO 100 j=1,P
          A(k,j)=A(k,j)/D
100     continue
        DO 1000 i=1,P
          IF(i.NE.k) THEN
            b=a(i,k)
            do 200 j=1,p
              a(i,j)=a(i,j)-b*a(k,j)
200         continue
            a(i,k)=-b/d
          endif
1000    continue
        a(k,k)=1/d
        RETURN
        END
