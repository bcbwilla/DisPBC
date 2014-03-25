*Deck DispF
      Subroutine DispF(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,FXYZ)
C
C     Compute the forces analytically due to Grimme DFT-D2 dispersion.
C
C     Arguments
C     In:
C     IOut   - Log file to print to
C     IPrint - Print flag.  Ranges 0-4. Higher prints more
C     IWhat  - 0 -> Energy, 1 -> Forces, 2 -> Trans. vector forces 
C              numerically (for consistency check) 
C     NAtoms - Number of atoms in unit cell
C     NPDir  - Number of periodic directions (translation vectors)
C     C      - Array of atom positions, including translation vectors
C              Atoms are stored first, followed by trans. vectors.
C     R0     - Array of atomic van der Waal radii.  Atom order is same
C              as in C array.
C     C6     - Arrary of atomic C6 coefficients.  Atom order is same as
C              in C array.
C     D6     - D6 parameter
C     S6     - Functional-dependent parameter
C     SR6    - Some other parameter related to van Der Waals radii
C     NUCExt - The number of unit cells to extend from the reference
C              unit cell in 1 direction.
C
C     Out:
C     FXYZ   - Array of forces on atoms and translation vectors.
C
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*), FXYZ(3,*)
      Dimension R0(*), C6(*)

 111  Format(' Atom',I5,I5)
 112  Format(' C6I=',F10.5,' C6J=', F10.5,' C6IJ=',F10.5,/,
     $       ' R0I=',F10.5,' R0J=',F10.5,' R0IJ=',F10.5,/,
     $       ' F=',E12.6,' DF=',E12.6,/,' RIJ=',F10.5,
     $       ' S6=',F10.5,' SR6=',F10.5,' D6=',F10.5)
 113  Format(' l=',I5,' m=',I5,' n=',I5)


      If(IPrint.gt.1) then
        Write(IOut,*)'Entering DispF.'
        endif


C     Compute forces on atoms first
C     Forces for interactions within the primary unit cell
      Epuc = 0
      Do 30 I=1,NAtoms
        Do 31 J=I+1,NAtoms

          C6I = C6(I)
          C6J = C6(J)
          C6IJ = Sqrt(C6I*C6J)
              
          R0I = R0(I)
          R0J = R0(J)
          R0IJ = R0I + R0J

C     For each coordinate
          Do 32 K=1,3 
C     kth component of atoms I and J
            CKI = C(K,I)
            CKJ = C(K,J) 
              
            RIJ = RSep(I,J,C,0,0,0,NAtoms,NPDir)
            F = FDmp(RIJ,D6,SR6,R0IJ)
            FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,K,
     $              C,0,0,0,0,NAtoms,NPDir,0)                

            Fk = -S6*C6IJ*(-6*F*(CKI-CKJ)/RIJ**8 + FDr/RIJ**6)

            FXYZ(K,I) = FXYZ(K,I) + Fk
            FXYZ(K,J) = FXYZ(K,J) - Fk

   32     Continue

          If(IPrint.gt.2) then
            Write(IOut,111)I,J
            Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,FDr,RIJ,
     $        S6,SR6,D6
            endif

   31   Continue
   30 Continue

      If(IPrint.gt.1) then
        Call PrintM(IOut,'Local forces on atoms:',FXYZ,NAtoms)
        endif

C     Forces for interactions between puc and neighbors
      If(NPDir.gt.0) then

C       Set loop limits based on dimension of periodicity
        NUC1 = NUCExt
        If(NPDir.ge.2) then
          NUC2 = NUCExt
        else
          NUC2 = 0
          endif
        
        If(NPDir.eq.3) then
          NUC3 = NUCExt
        else
          NUC3 = 0
          endif

        Do 40 I=1,NAtoms
          Do 41 J=1,NAtoms

            C6I = C6(I)
            C6J = C6(J)
            C6IJ = Sqrt(C6I*C6J)
          
            R0I = R0(I)
            R0J = R0(J)
            R0IJ = R0I + R0J

            Do 42 l=-NUC1,NUC1
              Do 43 m=-NUC2,NUC2
                Do 44 n=-NUC3,NUC3

C     Skip central unit cell
                  If(.not.((l.eq.0).and.(m.eq.0).and.(n.eq.0))) then

                    Do 45 K=1,3
C     kth component of atoms I and J
                      CKI = C(K,I)
                      CKJ = C(K,J) 
C     kth component of the TVs 
                      If(NPDir.ge.1) then
                        T1K = C(k,NAtoms+1)
                      else
                        T1K = 0
                        endif
                 
                      If(NPDir.ge.2) then
                        T2K = C(k,NAtoms+2)
                      else
                        T2K = 0
                        endif 

                      If(NPDir.ge.3) then
                        T3K = C(k,NAtoms+3)
                      else
                        T3K = 0
                        endif 
                      
                      RIJ = RSep(I,J,C,l,m,n,NAtoms,NPDir)
                      F = FDmp(RIJ,D6,SR6,R0IJ)
                      FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,K,
     $                            C,l,m,n,0,NAtoms,NPDir,0)                

                      Fk = 6*F*(CKI-CKJ-l*T1K-m*T2K-n*T3K)/RIJ**8
                      Fk = -Fk + FDr/RIJ**6
                      Fk = S6*C6IJ*Fk

                      FXYZ(K,I) = FXYZ(K,I) + Fk

   45                 Continue
                    
                    endif

                  If(IPrint.gt.2) then
                    Write(IOut,111)I,J
                    Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,RIJ,S6,SR6
                    endif

   44             Continue
   43           Continue
   42         Continue

   41       Continue
   40     Continue

        endif

      If(IPrint.gt.1) then
        Call PrintM(IOut,'All forces on atoms:',FXYZ,NAtoms)
        endif

C     Compute forces on translational vectors
      If(NPDir.gt.0) then

        If(IPrint.gt.1) then
          Write(IOut,*)'Forces on translational vectors.'
          endif

C     For each TV
        Do 50 IP=1,NPDir

          Do 51 I=1,NAtoms
            Do 52 J=1,NAtoms

              C6I = C6(I)
              C6J = C6(J)
              C6IJ = Sqrt(C6I*C6J)
          
              R0I = R0(I)
              R0J = R0(J)
              R0IJ = R0I + R0J

              Do 53 l=-NUC1,NUC1
                Do 54 m=-NUC2,NUC2
                  Do 55 n=-NUC3,NUC3

C     Skip central unit cell
                    If(.not.((l.eq.0).and.(m.eq.0).and.(n.eq.0))) then

                      Do 56 K=1,3

C     Kth component of atoms I and J
                        CKI = C(K,I)
                        CKJ = C(K,J) 

C     Get the k component of the TVs 
C     Translational vector 1
                        If(NPDir.ge.1) then
                          T1K = C(K,NAtoms+1)
                        else
                          T1K = 0
                          endif
                 
                        If(NPDir.ge.2) then
                          T2K = C(K,NAtoms+2)
                        else
                          T2K = 0
                          endif 

                        If(NPDir.ge.3) then
                          T3K = C(K,NAtoms+3)
                        else
                          T3K = 0
                          endif 

C     Get coefficient based on TV
                        If(IP.eq.1) then
                          Ia = l
                        else if(IP.eq.2) then
                          Ia = m
                        else if(IP.eq.3) then
                          Ia = n
                          endif
                      
                        RIJ = RSep(I,J,C,l,m,n,NAtoms,NPDir)
                        F = FDmp(RIJ,D6,SR6,R0IJ)
                        FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,k,C,
     $                               l,m,n,IP,NAtoms,NPDir,1)      

                        Fk = Ia*6*F*(CKI-CKJ-l*T1K-m*T2K-n*T3K)/RIJ**8
                        Fk = Fk + FDr/RIJ**6
                        Fk = S6*C6IJ*Fk
                        FXYZ(K,NAtoms+IP) = FXYZ(K,NAtoms+IP) + Fk

   56                   Continue

                      endif

   55               Continue
   54             Continue
   53           Continue

              If(IPrint.gt.2) then
                Write(IOut,111)I,J
                Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,RIJ,S6,SR6
                endif

   52         Continue
   51       Continue

   50     Continue
        endif
      
      If(IPrint.gt.1) then
        Call PrintM(IOut,'All forces:',FXYZ,NAtoms+NPDir)
        endif

      If(IPrint.gt.1) then
        Write(IOut,*)'Exiting DispF.'
        endif

      Return
      End
