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

      Do I=1,3

        Write(IOut,'(" *  Test",I3,"  * ")')I
C     Initialize test values
        Call Tests(IOut,I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)


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
