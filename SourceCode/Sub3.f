c
c-----------------------------------------------------------------------
c     V3  2026 07 : Vortex cylinder model
c
c     Section 3 of Li, Gaunaa, Pirrung, Horcas:
c     "A computationally efficient engineering aerodynamic model
c      for non-planar wind turbine rotors",
c     Wind Energ. Sci. 7, 75-104 (2022)
c
c     The converged axial induction factors abem(i) of the BEM
c     sweep define nsec+1 semi-infinite trailing vortex cylinders,
c     Eqs. (20)/(38):
c
c        gam_t(k) = 2 U0 ( a(k) - a(k-1) )
c
c     with ghost sections a(0) = a(nsec+1) = 0.
c     Cylinder k sits at the annulus boundary radius
c
c        Rk = rroot + (k-1) dr
c
c     and, for the non-planar (prebent) blade, starts at the
c     axial (downwind) position of the deflected blade at that
c     radius (Sect. 3.3: "the starting position of the cylindrical
c     vortex sheets follows the curved bound vortex line").
c
c     Radial induced velocity of one semi-infinite right cylinder,
c     Eqs. (10)+(11) - valid for evaluation points inside and
c     outside the cylinder, x is the axial distance between
c     evaluation point and cylinder start plane (u_r is symmetric
c     in x, only x**2 enters):
c
c        u_r(r,x) = -gam_t/(2 pi) sqrt(Rk/r)
c                   [ (2-k2)/sk K(k2) - 2/sk E(k2) ]
c
c        k2 = 4 r Rk / ( (Rk+r)**2 + x**2 ) ,  sk = sqrt(k2)
c
c     K, E: complete elliptic integrals of 1st/2nd kind.
c     Sign: the tip cylinder (gam_t < 0 for a wind turbine)
c     produces a positive = outboard radial velocity (wake
c     expansion), consistent with the Madsen correction of V2.
c     The total u_r at a blade section is the superposition
c     (sum) over all nsec+1 cylinders.
c-----------------------------------------------------------------------
c
      subroutine VCradial(nsc,dr,ursec)
c
      use mem
c
      integer nsc, j, k
      real dr
      real ursec(nsc)
      real rk, zk, rj, zj, gamt, a1, a2, x, ak2, sk, brack
      real ellK, ellE, dum
c
      do j = 1,nsc
         ursec(j) = 0.
      end do
c
      do j = 1,nsc
c
c        evaluation point: section midpoint, on the deflected blade
c
         rj = rroot + 0.5*dr + float(j-1)*dr
         call SPLINT(rsecsp,oopdefsp,oopdefS,ndes,rj,zj,dum)
c
c        superposition of all nsec+1 cylinders
c
         do k = 1,nsc+1
c
            if (k.eq.1) then
               a1 = 0.
            else
               a1 = abem(k-1)
            endif
            if (k.eq.nsc+1) then
               a2 = 0.
            else
               a2 = abem(k)
            endif
c
c           Eq. (20)/(38)
c
            gamt = 2.*vwind*(a2 - a1)
            if (abs(gamt).lt.1.e-10) goto 10
c
            rk = rroot + float(k-1)*dr
            if (rk.lt.1.e-3) goto 10
c
c           axial start position of cylinder k (non-planar blade)
c
            call SPLINT(rsecsp,oopdefsp,oopdefS,ndes,rk,zk,dum)
c
            x   = zj - zk
            ak2 = 4.*rj*rk/((rk+rj)**2 + x*x)
            if (ak2.gt.0.999999) ak2 = 0.999999
            if (ak2.lt.1.e-12  ) goto 10
            sk  = sqrt(ak2)
c
            brack = (2.-ak2)/sk * ellK(ak2) - 2./sk * ellE(ak2)
c
c           Eq. (10)
c
            ursec(j) = ursec(j)
     +               - gamt/(2.*pi) * sqrt(rk/rj) * brack
c
 10         continue
         end do
      end do
c
      return
      end
c
c-----------------------------------------------------------------------
c     complete elliptic integral of the FIRST kind K(m), 0 <= m < 1
c     polynomial approximation, Abramowitz & Stegun 17.3.34/35,
c     |error| < 2e-8
c-----------------------------------------------------------------------
c
      real function ellK(am)
c
      real am, am1, aa, bb
c
      am1 = 1. - am
      if (am1.lt.1.e-12) am1 = 1.e-12
c
      aa = 1.38629436112 + am1*(0.09666344259 + am1*(0.03590092383
     +     + am1*(0.03742563713 + am1*0.01451196212)))
      bb = 0.5 + am1*(0.12498593597 + am1*(0.06880248576
     +     + am1*(0.03328355346 + am1*0.00441787012)))
c
      ellK = aa - bb*log(am1)
c
      return
      end
c
c-----------------------------------------------------------------------
c     complete elliptic integral of the SECOND kind E(m), 0 <= m < 1
c     polynomial approximation, Abramowitz & Stegun 17.3.36/37,
c     |error| < 2e-8
c-----------------------------------------------------------------------
c
      real function ellE(am)
c
      real am, am1, aa, bb
c
      am1 = 1. - am
      if (am1.lt.1.e-12) am1 = 1.e-12
c
      aa = 1. + am1*(0.44325141463 + am1*(0.06260601220
     +     + am1*(0.04757383546 + am1*0.01736506451)))
      bb = am1*(0.24998368310 + am1*(0.09200180037
     +     + am1*(0.04069697526 + am1*0.00526449639)))
c
      ellE = aa - bb*log(am1)
c
      return
      end
c
