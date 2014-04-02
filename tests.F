*Deck Tests
      Subroutine Tests(IOut,I,C,C6,R0,FXYZ,D6,S6,SR6,NAtoms,NPDir)
C
C     Various test cases for DispPBC.
C

      Implicit Real*8(A-H,O-Z)
      Dimension C(3,*),C6(*),R0(*),FXYZ(3,*)
C     Standard values for tests.  In atomic units.
C     Parameters
C     Hyrdogen
      R0H = 1.719650773354206
      C6H = 2.428338775548422

C     Helium
      R0He = 1.738548034599857
      C6He = 1.387622157456241

C     Lithium
      R0Li = 1.417294593423796
      C6Li = 27.92589591880686

C     Carbon
      R0C = 2.494438484425882
      C6C = 30.35423469435528

C     Neon
      R0Ne = 2.135390520758520
      C6Ne = 10.92752448996790

C     Argon
      R0Ar = 2.740102880619339
      C6Ar = 79.96172682341592

C     Calcium
      R0Ca = 2.532233006917183
      C6Ca = 187.3289912565926

C     PBEPBE Functional
      D6 = 20.000000
      S6 = 0.7500000
      SR6 = 1.1000000

C     Distances in Bohr
C     1 Angstrom
      D1A = 1.889726124565 
C     2 Angstroms
      D2A = 3.779452249130
C     3 Angstroms
      D3A = 5.669178373695
C     10 Angstroms
      D10A = 18.89726124565
C     50 Angstroms
      D50A = 94.48630622825

C     Test cases
      If(I.eq.1) then
C     Test 1
        Write(IOut,*)'Two He atoms with large 1D unit cell.'
        Write(IOut,*)'Atoms are 4 ang. apart. 50 ang. TV.'
        NAtoms = 2
        NPDir = 1
C     Atoms
        C(1,1) =-D2A
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = D2A
        C(2,2) = 0.0
        C(3,2) = 0.0
C     Translation Vectors
        C(1,3) = D50A
        C(2,3) = 0.0
        C(3,3) = 0.0

C     Values for He atom
        Do J=1,NAtoms
          C6(J) = C6He
          R0(J) = R0He
          enddo
      
      else if(I.eq.2) then
C     Test 2
        Write(IOut,*)'Two He atoms with large 3D unit cell.'
        Write(IOut,*)'Atoms are 4 ang. apart. 50 ang. TVs.'
        NAtoms = 2
        NPDir = 3
C     Atoms
        C(1,1) =-D2A
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = D2A
        C(2,2) = 0.0
        C(3,2) = 0.0
C     Translation Vectors
        C(1,3) = D50A
        C(2,3) = 0.0
        C(3,3) = 0.0

        C(1,4) = 0.0
        C(2,4) = D50A
        C(3,4) = 0.0

        C(1,5) = 0.0
        C(2,5) = 0.0
        C(3,5) = D50A

C     Values for He atom
        Do J=1,NAtoms
          C6(J) = C6He
          R0(J) = R0He
          enddo

      else if(I.eq.3) then
C     Test 3
        Write(IOut,*)'Two He atoms with medium 3D unit cell.'
        Write(IOut,*)'Atoms are 4 ang. apart. 10 ang. TVs.'
        NAtoms = 2
        NPDir = 3
C     Atoms
        C(1,1) =-D2A
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = D2A
        C(2,2) = 0.0
        C(3,2) = 0.0
C     Translation Vectors
        C(1,3) = D10A
        C(2,3) = 0.0
        C(3,3) = 0.0

        C(1,4) = 0.0
        C(2,4) = D10A
        C(3,4) = 0.0

        C(1,5) = 0.0
        C(2,5) = 0.0
        C(3,5) = D10A

        Do J=1,NAtoms
          C6(J) = C6He
          R0(J) = R0He
          enddo

      else if(I.eq.4) then
C     Test 4
        Write(IOut,*)'Two He atoms with medium 3D unit cell.'
        Write(IOut,*)'Atoms are 4 ang. apart. 8 ang. TVs.'
        NAtoms = 2
        NPDir = 3
C     Atoms
        C(1,1) =-D2A
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = D2A
        C(2,2) = 0.0
        C(3,2) = 0.0
C     Translation Vectors
        C(1,3) = 4*D2A
        C(2,3) = 0.0
        C(3,3) = 0.0

        C(1,4) = 0.0
        C(2,4) = 4*D2A
        C(3,4) = 0.0

        C(1,5) = 0.0
        C(2,5) = 0.0
        C(3,5) = 4*D2A

        Do J=1,NAtoms
          C6(J) = C6He
          R0(J) = R0He
          enddo

      else if(I.eq.5) then
C     Test 5
        Write(IOut,*)'NeArHe triangle with large unit cell.'
        Write(IOut,*)'50 Ang TVs.'
        NAtoms = 3
        NPDir = 3
C     Atoms
        C(1,1) = 0.0
        C(2,1) = 0.0
        C(3,1) = 0.0

        C(1,2) = D2A
        C(2,2) = 0.0
        C(3,2) = D1A

        C(1,3) = D3A
        C(2,3) = D3A
        C(3,3) = 0.0
C     Translation Vectors
        C(1,4) = D50A
        C(2,4) = 0.0
        C(3,4) = 0.0

        C(1,5) = 0.0
        C(2,5) = D50A
        C(3,5) = 0.0

        C(1,6) = 0.0
        C(2,6) = 0.0
        C(3,6) = D50A

        
        C6(1) = C6Ne
        R0(1) = R0Ne
        C6(2) = C6Ar
        R0(2) = R0Ar
        C6(3) = C6He
        R0(3) = R0He

        endif


C     Build empty FXYZ
      Do J=1,NAtoms+NPDir
        Do K=1,3
          FXYZ(K,J) = 0D0
          enddo
        enddo

      Return
      End
