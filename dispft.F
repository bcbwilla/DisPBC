*Deck DispFT
      Subroutine DispFT(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,FXYZ)
C
C     Compute the forces on translational vectors due to
C     Grimme DFT-D2 dispersion numerically
C
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*),C2(3,NAtoms+NPDir),FXYZ(3,*)
      Dimension FXYZ2(3,NAtoms+NPDir)
      Dimension R0(*), C6(*), FP(4), E(4)

 111  Format(' Atom',I5,I5)
 112  Format(' C6I=',F10.5,' C6J=', F10.5,' C6IJ=',F10.5,/,
     $       ' R0I=',F10.5,' R0J=',F10.5,' R0IJ=',F10.5,/,
     $       ' F=',E12.6,' DF=',E12.6,/,' RIJ=',F10.5,
     $       ' S6=',F10.5,' SR6=',F10.5,' D6=',F10.5)
 113  Format(' l=',I5,' m=',I5,' n=',I5)

      dT = 0.000001

      If(IPrint.gt.1) then
        Write(IOut,*)'Entering DispFT.'
        endif

C     Copy coordinates into new array

      Do JA=1,NAtoms+NPDir
        Do IA=1,3
          C2(IA,JA) = C(IA,JA)
          enddo
        enddo

C     Five point stencil function shift coefficients
      FP(1) = 2D0*dT
      FP(2) = dT
      FP(3) = -dT
      FP(4) = -2D0*dT

C     Each translation vector.    
      Do 10 N=1,NPDir
C     Each coordinate
        Do 20 I=1,3
C     Do each Five Point Stencil term.
          Do 30 J=1,4
C     Shift translation vectors according to five pt. term.
            C2(I,NAtoms+N) = C(I,NAtoms+N) + FP(J)

            Write(*,*)'TV: ',N,' Coord: ',I,' Term: ',J
            Write(*,*)'FPJ: ',FP(J)
        
C    Compute the energy with new TV
            Call DispE(IOut,0,NAtoms,NPDir,C2,R0,C6,D6,S6,SR6,
     $        NUCExt,En)

            Write(*,*)'Energy: ',En
            E(J) = En

C    Restore translation vector for next calc
            C2(I,NAtoms+N) = C(I,NAtoms+N)


   30       Continue

          FT = -FPtSt(E,dT)
          FXYZ2(I,NAtoms+N) = FT


   20     Continue

   10   Continue

      
      If(IPrint.gt.1) then
        Call PrintM(IOut,'Numerical Tv Forces:',FXYZ2,NAtoms+NPDir)
        endif

      If(IPrint.gt.1) then
        Write(IOut,*)'Exiting DispFT.'
        endif

            Write(*,*)'Numerical TV Forces:'
            Do JA=1,NAtoms+NPDir
              Write(*,*)(FXYZ2(IA,JA),IA=1,3)
              enddo

      Return
      End
