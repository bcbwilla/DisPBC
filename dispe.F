*Deck DispE
      Subroutine DispE(IOut,IPrint,NAtoms,NPDir,C,R0,C6,D6,
     $  S6,SR6,NUCExt,ED)
C
C     Compute the energy due to Grimme DFT-D2 dispersion.
C
C
C     Arguments
C     In:
C     IOut   - Log file to print to
C     IPrint - Print flag.  Ranges 0-3
C              IPrint = 0 -> Nothing
C              IPrint = 1 -> Basic Information
C              IPrint = 2 -> Basic Debugging
C              IPrint = 3 -> Extreme Debugging
C     NAtoms - Number of atoms in unit cell
C     NPDir  - Number of periodic directions (translation vectors)
C     IAn    - Array of atomic numbers of atoms
C     C      - Array of atom positions, including translation vectors
C              Atoms are stored first, followed by trans. vectors.
C     XXJXDHUH
C     Out:
C     ED     - Dispersion Energy
C     FXYZ   - Array of forces on atoms and translation vectors
C
      Implicit Real*8(A-H,O-Z)

      Real*8 C(3,*), R0(*), C6(*)

 111  Format(' Atom',I5,I5)
 112  Format(' C6I=',F10.5,' C6J=', F10.5,' C6IJ=',F10.5,/,
     $       ' R0I=',F10.5,' R0J=',F10.5,' R0IJ=',F10.5,/,
     $       ' F=',F10.5,/,' RIJ=',F10.5,' S6=',F10.5,' SR6=',F10.5)
 113  Format(' l=',I5,' m=',I5,' n=',I5)


      If(IPrint.gt.1) then
        Write(IOut,*)'Entering DispE.'
        endif


C     Compute the dispersion energy of atoms inside unit cell only
      If(IPrint.gt.1) then
        Write(IOut,*)'Disp. energy for interactions in unit cell'
        endif

      Epuc = 0
      Do 10 I=1,NAtoms
        Do 11 J=I+1,NAtoms
          
          C6I = C6(I)
          C6J = C6(J)
          C6IJ = Sqrt(C6I*C6J)
          
          R0I = R0(I)
          R0J = R0(J)
          R0IJ = R0I + R0J
         
          RIJ = RSep(I,J,C,0,0,0,NAtoms,NPDir)
          F = FDmp(RIJ,D6,SR6,R0IJ)
 
          Epuc = Epuc - S6*(C6IJ/RIJ**6)*F 

          If(IPrint.gt.2) then
            Write(IOut,111)I,J
            Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,RIJ,S6,SR6
            endif

   11     Continue
   10   Continue

C     Compute the dispersion energy of atoms inside the unit cell
C     interacting with atoms in neighbor unit cells

      Enuc = 0
      If(NPDir.gt.0) then

        If(IPrint.gt.1) then
          Write(IOut,*)'Disp. energy for interactions out of unit cell'
          endif

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

        Do 20 I=1,NAtoms
          Do 21 J=1,NAtoms
            C6I = C6(I)
            C6J = C6(J)
            C6IJ = Sqrt(C6I*C6J)
          
            R0I = R0(I)
            R0J = R0(J)
            R0IJ = R0I + R0J

            Do 22 l=-NUC1,NUC1
              Do 23 m=-NUC2,NUC2
                Do 24 n=-NUC3,NUC3

C     Skip central unit cell
                  If(.not.((l.eq.0).and.(m.eq.0).and.(n.eq.0))) then

                    RIJ = RSep(I,J,C,l,m,n,NAtoms,NPDir)
                    F = FDmp(RIJ,D6,SR6,R0IJ)

                    Enuc = Enuc - S6*(C6IJ/RIJ**6)*F 

                    If(IPrint.gt.3) then
                      Write(IOut,113)l,m,n
                      Write(IOut,*)'RIJ: ',RIJ
                      endif

                    endif                  

   24             Continue
   23           Continue
   22         Continue

            If(IPrint.gt.2) then
              Write(IOut,111)I,J
              Write(IOut,112)C6I,C6J,C6IJ,R0I,R0J,R0IJ,F,RIJ,S6,SR6
              endif

   21       Continue
   20     Continue

        endif

      ED = Epuc + Enuc

      If(IPrint.gt.1) then
        Write(IOut,*)'EPuc: ',EPuc
        Write(IOut,*)'Enuc: ',ENuc
        endif

      If(IPrint.gt.1) then
        Write(IOut,*)'Exiting DispE.'
        endif

      Return
      End
