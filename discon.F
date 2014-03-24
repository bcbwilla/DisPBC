*Deck DispConv
C Description
C
C
      Subroutine DisCon(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,ED)

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

C     The energy cutoff
      ECut = 0.0001

C     Calculate the first dispersion energy without augmenting the UC

      Call DispE(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,S6,SR6,NUCExt,E1)

      If(IPrint.gt.1) then
        Write(IOut,*)' Initial dispersion energy: ',E1
        endif

C     Difference in energy for systems with different extension lengths
      EDiff = Abs(E1) + 4*ECut

C     Number of unit cells to add to the original cell in 1 dimension
      NUCExt = 0
C     Loop counter
      ik = 1
C     Limit on number of unit cell additions
      UCL = 50
      Do While ((EDiff.gt.ECut).and.(ik.le.UCL))
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
        If((EDiff.lt.ECut).and.(IPrint.gt.0)) then
          Write(IOut,333)ik,NUCExt,E2,Ediff,NAt2
          endif

C     Update values for next iteration
        ik = ik + 1
        E1 = E2
        E2 = 0D0

        enddo

      ED = E1

      If(IPrint.gt.0) then
        Write(IOut,*)'Exiting DisCon.'
        endif

      End


