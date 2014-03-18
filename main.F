      Program Main
C
C     This program runs tests of DisPBC.
C
      Implicit Real*8(A-H,O-Z)
      Parameter(N=50)
      Dimension C(3,N),IAn(N),C6(N),R0(N),FXYZ(3,N)

C     Log file
      IOut = 42
      Open(IOut,file='iout.log')

C     Test cases
C -------------------------------------------------------------------
C     Test 1
      NAtoms = 2
      NPDir = 3

C     Atom List
      IAn(1) = 2
      IAn(2) = -2
      IAn(3) = -2
      IAn(4) = -2

C     Atoms
      C(1,1) = 0.0
      C(2,1) = 0.0
      C(3,1) = 0.0

      C(1,2) = 1.0
      C(2,2) = 0.0
      C(3,2) = 0.0

C     Translation Vectors
      C(1,3) = 3.0
      C(2,3) = 0.0
      C(3,3) = 0.0

      C(1,4) = 0.0
      C(2,4) = 3.0
      C(3,4) = 0.0

      C(1,5) = 0.0
      C(2,5) = 0.0
      C(3,5) = 3.0

C     Helium atom values, PBEPBE functional.
      D6 = 20.0
      S6 = 0.75
      SR6 = 1.1

C     Values for He atom
      Do I=1,NAtoms
        C6(I) = 1.387622157456241
        R0(I) = 3.477096069199714
        enddo

      Do J=1,NAtoms
        Write(*,*)'R1J ',R0(J)
        enddo


C -------------------------------------------------------------------

C     Build empty FXYZ
    
      Do j=1,NAtoms+NPDir
        Do i=1,3
          FXYZ(i,j) = 0D0
          enddo
        enddo

      Write(IOut,*)'Atoms initialized.'

      IWhat = 0
      IPrint = 2

      Call DisPBC(IOut,IPrint,IWhat,IEmpDi,NAtoms,NPDir,
     $  IAn,C,D6,S6,SR6,C6,R0,ED,FXYZ)

      Write(IOut,*)'Energy: ',ED

      IWhat = 1
      Call DisPBC(IOut,IPrint,IWhat,IEmpDi,NAtoms,NPDir,
     $  IAn,C,D6,S6,SR6,C6,R0,ED,FXYZ)




      Close(IOut)

      Write(*,*)'Test complete.'
      
      Return
      End

