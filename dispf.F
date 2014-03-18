*Deck DispF
      Subroutine DispF(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,FXYZ)
C
C     Compute the forces on atoms and translational vectors due to
C     Grimme DFT-D2 dispersion.
C
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*), FXYZ(3,*)
      Dimension R0(*), C6(*)

 111  Format(' Atom',I5,I5)
 112  Format(' C6I=',F6.4,' C6J=', F6.4,' C6IJ=',F6.4,/,
     $       ' R0I=',F6.4,' R0J=',F6.4,' R0IJ=',F6.4,/,
     $       ' F=',E12.6,/,' RIJ=',F6.4,' S6=',F6.4,' SR6=',F6.4)
 113  Format(' l=',I5,' m=',I5,' n=',I5)


      If(IPrint.gt.1) then
        Write(IOut,*)'Entering DispF.'
        endif


C     Compute forces on atoms first
C     Forces for interactions within the primary unit cell
      If(IPrint.gt.1) then
        Write(IOut,*)'Forces for interactions within unit cell.'
        endif

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
            KI = C(K,I)
            KJ = C(K,J) 
              
            RIJ = RSep(I,J,C,0,0,0,NAtoms,NPDir)
            F = FDmp(RIJ,D6,SR6,R0IJ)
            FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,K,
     $              C,0,0,0,0,NAtoms,NPDir,0)                

            Fk = 6*F*(KI - KJ)/RIJ**8 - FDr/RIJ**6

            FXYZ(K,I) = FXYZ(K,I) + Fk
            FXYZ(K,J) = FXYZ(K,J) - Fk

   32     Continue

          If(IPrint.gt.2) then
            Write(IOut,111)I,J
            Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,RIJ,S6,SR6
            endif

   31   Continue
   30 Continue

      If(IPrint.gt.1) then
        Call PrintM(IOut,'Local forces on atoms:',FXYZ,NAtoms+NPDir)
        endif

C     Forces for interactions between primary unit cell and neighbor unit cells
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

        If(IPrint.gt.1) then
          Write(IOut,*)'Forces for interactions out of unit cell.'
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
                      KI = C(K,I)
                      KJ = C(K,J) 
C     kth component of the TVs 
                      If(NPDir.ge.1) then
                        T1k = C(k,NAtoms+1)
                      else
                        T1k = 0
                        endif
                 
                      If(NPDir.ge.2) then
                        T2k = C(k,NAtoms+2)
                      else
                        T2k = 0
                        endif 

                      If(NPDir.ge.3) then
                        T3k = C(k,NAtoms+3)
                      else
                        T3k = 0
                        endif 
                      
                      RIJ = RSep(I,J,C,l,m,n,NAtoms,NPDir)
                      F = FDmp(RIJ,D6,SR6,R0IJ)
                      FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,K,
     $                            C,l,m,n,0,NAtoms,NPDir,0)                

                      Fk = 6*F*(KI - KJ - l*Tk1 - m*Tk2 - n*Tk3)/RIJ**8
                      Fk = Fk - FDr/RIJ**6

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
        Call PrintM(IOut,'All forces on atoms:',FXYZ,NAtoms+NPDir)
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
                        KI = C(K,I)
                        KJ = C(K,J) 

C     Get the k component of the TVs 
C     Translational vector 1
                        If(NPDir.ge.1) then
                          T1k = C(k,NAtoms+1)
                        else
                          T1k = 0
                          endif
                 
                        If(NPDir.ge.2) then
                          T2k = C(k,NAtoms+2)
                        else
                          T2k = 0
                          endif 

                        If(NPDir.ge.3) then
                          T3k = C(k,NAtoms+3)
                        else
                          T3k = 0
                          endif 

                      
                        RIJ = RSep(I,J,C,l,m,n,NAtoms,NPDir)
                        F = FDmp(RIJ,D6,SR6,R0IJ)
                        FDr = FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,k,C,
     $                               l,m,n,IP,NAtoms,NPDir,1)              

                        Fk = -6*F*(KI-KJ-l*T1K-m*T2K-n*T3K)/RIJ**8
                        Fk = Fk - FDr/RIJ**6

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
