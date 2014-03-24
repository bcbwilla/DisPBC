*Deck DisPBC
      Subroutine DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $  C,D6,S6,SR6,C6,R0,ED,FXYZ)
C
C     Main subroutine to compute the energy and forces due to Grimme
C     DFT-D2 dispersion correction.
C
C
C     Arguments
C     In:
C     IOut   - Log file to print to
C     IPrint - Print flag.  Ranges 0-4. Higher prints more
C     IWhat  - 0 -> Compute energy,  1 -> Compute forces 
C     NAtoms - Number of atoms in unit cell
C     NPDir  - Number of periodic directions (translation vectors)
C     IAn    - Array of atomic numbers of atoms
C     C      - Array of atom positions, including translation vectors
C              Atoms are stored first, followed by trans. vectors.
C
C     Out:
C     ED     - Dispersion Energycl
C     FXYZ   - Array of forces on atoms and translation vectors
C
      Implicit Real*8(A-H,O-Z)

      Integer IAn(NAtoms+NPDir)
      Dimension C(3,NAtoms+NPDir), FXYZ(3,NAtoms+NPDir)
      Dimension C6(NAtoms), R0(NAtoms), FXYZ2(3,NAtoms+NPDir)
      Logical First
      Save NUCExt
      Data NUCExt/0/

      Write(IOut,*)'Entering DisPBC.'

      If(NUCExt.eq.0) then
        Write(IOut,*)'Computing dispersion energy.'
C     Determine how many unit cells out to include disp. interactions
        Call DisCon(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,ED)
        endif

      If(IWhat.eq.1) then
        Write(IOut,*)'Computing dispersion correction to forces.'
C     Compute the forces
        Call DispF(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,FXYZ)

C     Compute the forces on TVs numerically
        Call DispFT(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,FXYZ2)
        endif

      Write(IOut,*)'Exiting DisPBC.'

      Return
      End      
