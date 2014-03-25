C
C     Compute the energy and forces due to van der Waals dispersion
C     energy in periodic systems using the DFT-D2 correction proposed
C     by Grimme.
C
C     Grimme, S. (2006), Semiempirical GGA-type density functional 
C     constructed with a long-range dispersion correction. 
C     J. Comput. Chem., 27: 1787–1799. doi: 10.1002/jcc.20495
C
C
*Deck DisPBC
      Subroutine DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $  C,R0,C6,D6,S6,SR6,ED,FXYZ)
C
C     Main subroutine to compute the energy and forces
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
C
C     Out:
C     ED     - Dispersion Energy
C     FXYZ   - Array of forces on atoms and translation vectors
C
      Implicit Real*8(A-H,O-Z)

      Integer IAn(NAtoms+NPDir)
      Real*8 C(3,*), FXYZ(3,*), C6(*), R0(*)
      Save NUCExt
      Data NUCExt/0/

      Write(IOut,*)'Entering DisPBC.'

C     NUCExt is the number of unit cells to add in one direction
      If(NUCExt.eq.0) then
        Write(IOut,*)'Computing dispersion energy.'
C     Determine how many unit cells out to include disp. interactions
C     and compute the energy
        Call DisCon(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,ED)
        endif

      If(IWhat.eq.1) then
        Write(IOut,*)'Computing dispersion correction to forces.'
C     Compute the forces
        Call DispF(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,FXYZ)

      else if(IWhat.eq.2) then
        Write(IOut,*)'Computing dispersion forces on translation',
     $   ' vectors numerically.'
        Call DispFT(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $    S6,SR6,NUCExt,FXYZ)
        endif

      Write(IOut,*)'Exiting DisPBC.'

      Return
      End      
