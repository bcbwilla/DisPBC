      Program Main
C
C     This program runs tests of DisPBC.
C
      Implicit Real*8(A-H,O-Z)
      Parameter(N=50)
      Dimension C(3,N),C6(N),R0(N),FXYZ(3,N),FXYZ2(3,N)

      Write(*,*)'Running tests.'

C     Log file
      IOut = 42
      Open(IOut,file='iout.log')

C     Print level
      IPrint = 1

C     Test
      I = 1

      Write(IOut,'(" *  Test",I3,"  * ")')I
C     Initialize test values
      Call Tests(IOut,I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)

      Write(IOut,'(" NAtoms =",I2," NPDir =",I2)')NAtoms,NPDir
      Call PrintM(IOut,"C:",C,NAtoms+NPDir)
      Write(IOut,*)

C     Energy
      IWhat = 0
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
 $      C,R0,C6,D6,S6,SR6,ED,FXYZ)

      Write(IOut,*)'Energy: ',ED

C     Forces
      IWhat = 1
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
 $      C,R0,C6,D6,S6,SR6,ED,FXYZ)

      Call PrintM(IOut,"Analytical Forces:",FXYZ,NAtoms+NPDir)

C      Numerical translational vector forces
      IWhat = 2
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
 $      C,R0,C6,D6,S6,SR6,ED,FXYZ2)

      Call PrintM(IOut,"Numerical Forces on TVs:",FXYZ2,NAtoms+NPDir)

      Close(IOut)

      Write(*,*)'Test complete.'
      
      End
