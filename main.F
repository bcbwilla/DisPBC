      Program Main
C
C     This program runs tests of DisPBC.
C
      Implicit Real*8(A-H,O-Z)
      Parameter(N=50)
      Dimension C(3,N),C6(N),R0(N),FXYZ(3,N)

      Write(*,*)'Running tests.'

C     Log file
      IOut = 42
      Open(IOut,file='iout.log')

C     Print level
      IPrint = 2

      Do I=1,1
C     Initialize test values
        Call Test(I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)

        Write(IOut,'(" *  Test",I3,"  * ")')I
        Write(IOut,'(" NAtoms =",I2," NPDir =",I2)')NAtoms,NPDir
        Call PrintM(IOut,"C:",C,NAtoms+NPDir)
        Write(IOut,*)

C     Energy
        IWhat = 0
        Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $    C,D6,S6,SR6,C6,R0,ED,FXYZ)

        Write(IOut,*)'Energy: ',ED

C     Forces
        IWhat = 1
        Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $    C,D6,S6,SR6,C6,R0,ED,FXYZ)


        Call PrintM(IOut,"Forces:",FXYZ,NAtoms+NPDir)       

        Write(IOut,'(/,/)')
        enddo

      Close(IOut)

      Write(*,*)'Test complete.'
      
      End

*Deck Test
      Subroutine Test(I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)

      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*),C6(*),R0(*),FXYZ(3,*)

C     Test cases

      If(I.eq.1) then
C     Test 1
        NAtoms = 2
        NPDir = 1
C     Atoms
        C(1,1) =-3.779452249130
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = 3.779452249130
        C(2,2) = 0.0
        C(3,2) = 0.0
C     Translation Vectors
        C(1,3) = 94.48630622825
        C(2,3) = 0.0
        C(3,3) = 0.0

C     Helium atom values, PBEPBE functional.
        D6 = 20.000000
        S6 = 0.7500000
        SR6 = 1.1000000

C     Values for He atom
        Do J=1,NAtoms
          C6(J) = 1.387622157456241
          R0(J) = 1.0
          enddo
      
      else if(I.eq.1) then
C     Test 2
        NAtoms = 2
        NPDir = 3
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
        Do J=1,NAtoms
          C6(J) = 1.387622157456241
          R0(J) = 3.477096069199714
          enddo

        endif


C     Build empty FXYZ
      Do J=1,NAtoms+NPDir
        Do K=1,3
          FXYZ(K,J) = 0D0
          enddo
        enddo

      Return
      End
