*Deck DispConv
      Subroutine DisCon(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,ECon,NUCExt,ED)
C
C     Determine how many neighboring cells to include in calculations.
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
C     ECon   - Change in dispersion energy that signifies convergence
C     NUCExt - The number of unit cells to extend from the reference
C              unit cell in 1 direction.
C
C     Out:
C     ED     - Dispersion Energy
C
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,NAtoms+NPDir), R0(NAtoms), C6(NAtoms)

111   Format(' Cycle',I3,/,'  NUCExt:',I3,' En:',F15.11,' Ediff:',
     $  F13.9,' NAt2:',I8)
222   Format('  Convergence failed after',I3,' cycles.',/,'   NUCExt:',
     $  I3,' En:',F15.11,'   Ediff:',F17.12,'   NAt2:',I8)
333   Format('  Convergence met after',I3,' cycles.',/,'   NUCExt:',
     $  I3,' En:',F15.11,'   Ediff:',F17.12,'   NAt2:',I8)

      If(IPrint.gt.0) then
        Write(IOut,*)'Entering DisCon.'
        endif

C     Determine energy convergence criteria
      If(ECon.eq.0) then
        ECon = 0.00001
        endif

C     Calculate the first dispersion energy without augmenting the UC

      Call DispE(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,S6,SR6,NUCExt,E1)

      If(IPrint.gt.1) then
        Write(IOut,*)' Initial dispersion energy: ',E1
        endif

C     Difference in energy for systems with different extension lengths
      EDiff = Abs(E1) + 4*ECon

C     Number of unit cells to add to the original cell in 1 dimension
      NUCExt = 0
C     Loop counter
      ik = 1
C     Limit on number of unit cell additions
      UCL = 50
      Do While ((EDiff.gt.ECon).and.(ik.le.UCL))
        NUCExt = NUCExt + 1

C     Number of atoms in extended system
        NAt2 = NAtoms*(2*NUCExt+1)**NPDir

C     Convergence failure
        If(ik.ge.UCL) then
          Write(IOut,222)ik,NUCExt,E1,EDiff,NAt2
          Write(IOut,*)'Out of cycles in DispConv.'
  999     Stop
          endif

        Call DispE(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,S6,SR6,
     $    NUCExt,E2)

        EDiff = ABS(E1 - E2)

        If(IPrint.gt.1) then
          Write(IOut,111)ik,NUCExt,E2,Ediff,NAt2
          endif

C     Successful convergence
        If((EDiff.lt.ECon).and.(IPrint.gt.0)) then
          Write(IOut,333)ik,NUCExt,E2,Ediff,NAt2
          endif

C     Update values for next iteration
        ik = ik + 1
        E1 = E2
        E2 = 0D0

        enddo

      ED = E1

      If(IPrint.gt.0) then
        Write(IOut,*)'Cell extended ',NUCExt,' translations.'
        Write(IOut,*)'Exiting DisCon.'
        endif

      End


