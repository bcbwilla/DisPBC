C
C     A collection of commonly used utility functions and subroutines
C

*Deck Rsep
      Real*8 Function RSep(I,J,C,l,m,n,NAtoms,NPDir)
C
C     Computes the separation between atoms I and J, with
C     J possibly being translated by translation vectors.
C
C     Arguments
C     I/J    - Position of Ith/Jth atom in C
C     C      - Array of atom and translation vector (TV) coordinates
C     l/m/n  - Number of translations of TV1/TV2/TV3
C     NAtoms - Number of atoms in primary unit cell  
C
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,NAtoms+NPDir)

C     Atoms I and J
      XI = C(1,I)
      YI = C(2,I)
      ZI = C(3,I)
      XJ = C(1,J)
      YJ = C(2,J)
      ZJ = C(3,J)

C     Translational vectors
      If(NPDir.ge.1) then
        T1X = C(1,NAtoms+1)
        T1Y = C(2,NAtoms+1)
        T1Z = C(3,NAtoms+1)
      else
        T1X = 0D0
        T1Y = 0D0
        T1Z = 0D0
        endif
     
      If(NPDir.ge.2) then
        T2X = C(1,NAtoms+2)
        T2Y = C(2,NAtoms+2)
        T2Z = C(3,NAtoms+2)
      else
        T2X = 0D0
        T2Y = 0D0
        T2Z = 0D0
        endif 

      If(NPDir.ge.3) then
        T3X = C(1,NAtoms+3)
        T3Y = C(2,NAtoms+3)
        T3Z = C(3,NAtoms+3)
      else
        T3X = 0D0
        T3Y = 0D0
        T3Z = 0D0
        endif 

C     Compute the separation distance
      XTrm = XI - (XJ + l*T1X + m*T2X + n*T3X)
      YTrm = YI - (YJ + l*T1Y + m*T2Y + n*T3Y)
      ZTrm = ZI - (ZJ + l*T1Z + m*T2Z + n*T3Z)

      Rsep = Sqrt(XTrm**2 + YTrm**2 + ZTrm**2)

      Return
      End


*Deck FDmp
      Real*8 Function FDmp(RIJ,D6,SR6,R0IJ)
C
C     Compute the damping function for a given atom separation RIJ
C
C     Arguments
C     RIJ    - Separation of atoms RI and RJ
C     D6     - D6 parameter
C     SR6    - SR6 parameter
C     R0IJ   - Sum of atomic van der Waal radii for atoms I,J 
C
      Implicit Real*8 (A-H,O-Z)

      FDmp = 1D0 / (1D0 + Exp(-D6*(RIJ/(SR6*R0IJ) - 1D0)))
  
      Return
      End


*Deck FDmpDr
      Real*8 Function FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,K,C,l,m,n,IP,
     $  NAtoms,NPDir,Iwrt)
C
C     Derivative of FDmp.  If WRT = "TV", it is with respect to TVk.
C     Otherwise, it is with respect to RIk.
C
C     Arguments
C     RIJ    - Separation of atoms RI and RJ
C     D6     - D6 parameter
C     SR6    - SR6 parameter
C     R0IJ   - Sum of atomic van der Waal radii for atoms I,J 
C     F      - Value of damping function
C     I/J    - Position of Ith/Jth atom in C
C     K      - Coordinate to differentiate w.r.t. K=1,2,3 -> x,y,z
C     C      - Array of atom and translation vector (TV) coordinates
C     l/m/n  - Number of translations of TV1/TV2/TV3
C     IP     - The translation vector to differentiate w.r.t. (1,2,3)
C     NAtoms - Number of atoms in primary unit cell  
C     NPDir  - Number of periodic directions (TVs)
C     Iwrt   - Do derivative w.r.t. coordinate or TV.
C              Iwrt = 0 -> coordinate, Iwrt = 1 -> TV
C
      Implicit Real*8 (A-H,O-Z)
      Dimension C(3,NAtoms+NPDir)

C     Kth component of Atoms I, J
      CKI = C(K,I)
      CKJ = C(K,J)
     
C     Get the Kth component of the TVs 
      If(NPDir.ge.1) then
        T1K = C(K,NAtoms+1)
      else
        T1K = 0D0
        endif
     
      If(NPDir.ge.2) then
        T2K = C(K,NAtoms+2)
      else
        T2K = 0D0
        endif 

      If(NPDir.ge.3) then
        T3K = C(K,NAtoms+3)
      else
        T3K = 0D0
        endif 

      FTmp = D6*(CKI - CKJ - l*T1K - m*T2K - n*T3K)
      FTmp = FTmp * Exp(-D6*(RIJ/(SR6*R0IJ)-1D0)) * (F**2)
      FTmp = FTmp / (RIJ*(SR6*R0IJ))

      If(Iwrt.eq.1) then
C     Differential w.r.t. T

C     See if you multiply the front by l, m, or n.  This depends on
C     what TV you're differentiating with respect to.
        If(IP.eq.1) then
          Ia = l
        else if(IP.eq.2) then
          Ia = m
        else if(IP.eq.3) then
          Ia = n
          endif
        
        FTmp = FTmp*(-Ia)
        endif

      FDmpDr = FTmp
    
      Return
      End

*Deck FPtSt
      Real*8 Function FPtSt(F,dT)
C
C     Compute numerical derivative using Five point stencil method
C
C     Arguments
C     F      - Array of function evaluations corresponding to terms
C              of five point method  
C     dT     - Size of differential
C
      Implicit Real*8(A-H,O-Z)
      Real*8 F(4)

      FPtSt = (-F(1) + 8*F(2) - 8*F(3) + F(4))/(12*dT)

      Return
      End


*Deck PrintM
      Subroutine PrintM(IOut,Msg,A,N)
C
C     Prints an (3,N) matrix in tidy form.
C
C     For debugging purposes only.
C
C     Arguments
C     IOut   - Output file number
C     Msg    - Message to print above matrix
C     A      - (3,N) matrix to print
C     N      - Number of rows in matrix
C
      Integer N
      Character*(*) Msg
      Real*8 A(3,N)

 555  Format("  ",F21.14,F21.14,F21.14)

      Write(IOut,*)Msg
      Do j=1,N
        Write(IOut,555)(A(i,j),i=1,3)
        enddo 

      Return
      End
