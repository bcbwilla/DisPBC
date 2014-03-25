*Deck DispFT
      Subroutine DispFT(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,FXYZ)
C
C     Compute the forces on translational vectors due to
C     Grimme DFT-D2 dispersion numerically.  This is mostly
C     for consistency checks.
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
      Parameter (dT=0.0000001)
      Real*8 C(3,*),C2(3,NAtoms+NPDir),FXYZ(3,*)
      Real*8 R0(*), C6(*), FP(4), E(4)


 111  Format(' Atom',I5,I5)
 112  Format(' C6I=',F10.5,' C6J=', F10.5,' C6IJ=',F10.5,/,
     $       ' R0I=',F10.5,' R0J=',F10.5,' R0IJ=',F10.5,/,
     $       ' F=',E12.6,' DF=',E12.6,/,' RIJ=',F10.5,
     $       ' S6=',F10.5,' SR6=',F10.5,' D6=',F10.5)
 113  Format(' l=',I5,' m=',I5,' n=',I5)


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
        
C    Compute the energy with new TV
            Call DispE(IOut,0,NAtoms,NPDir,C2,R0,C6,D6,S6,SR6,
     $        NUCExt,En)

            E(J) = En

C    Restore translation vector for next calc
            C2(I,NAtoms+N) = C(I,NAtoms+N)

   30       Continue

          FT = -FPtSt(E,dT)
          FXYZ(I,NAtoms+N) = FT

   20     Continue

   10   Continue

      If(IPrint.gt.1) then
        Write(IOut,*)'Exiting DispFT.'
        endif

      Return
      End
