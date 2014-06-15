      Program Main
C
C     This program runs tests of DisPBC.
C
      Implicit Real*8(A-H,O-Z)
      Parameter(N=50)
      Logical DoCon,DoNum
      Dimension C(3,N),C6(N),R0(N),FXYZ(3,N),FXYZ2(3,N)

      Write(*,*)'Running tests.'

C     Log file
      IOut = 42
      Open(IOut,file='iout.log')

C     Print level
      IPrint = 2

C     Test
      I = 7

      ECon = 0
      DoCon = .true.
      NUCExt = 0
      DoNum = .false.

      Write(IOut,'(" *  Test",I3,"  * ")')I
C     Initialize test values
      Call Tests(IOut,I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)

      Write(IOut,'(" NAtoms =",I2," NPDir =",I2)')NAtoms,NPDir
      Call PrintM(IOut,"C:",C,NAtoms+NPDir)
      Write(IOut,*)


C     Energy
      IWhat = 0
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $  C,R0,C6,D6,S6,SR6,ECon,DoCon,DoNum,NUCExt,ED,FXYZ)

      Write(IOut,*)'Energy: ',ED

C     Forces
      IWhat = 1
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir,
     $  C,R0,C6,D6,S6,SR6,ECon,DoCon,DoNum,NUCExt,ED,FXYZ)

      Call PrintM(IOut,"Analytical Forces:",FXYZ,NAtoms+NPDir+1)

C      Numerical translational vector forces
      IWhat = 2
      DoNum = .true.
      Call DisPBC(IOut,IPrint,IWhat,NAtoms,NPDir, 
     $ C,R0,C6,D6,S6,SR6,ECon,DoCon,DoNum,NUCExt,ED,FXYZ2)

      Call PrintM(IOut,"Analytical TV Forces:",FXYZ,NAtoms+NPDir+1)
      Call PrintM(IOut,"Numerical TV Forces:",FXYZ2,NAtoms+NPDir+1)

      Close(IOut)

      Write(*,*)'Test complete.'
      
      End
