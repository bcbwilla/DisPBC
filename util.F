*Deck Rsep
      Real*8 Function RSep(I,J,C,l,m,n,NAtoms,NPDir)
C
C     Compute the separation between atoms I and J, with
C     J possibly being translated by translation vectors.
Cmak
      Implicit Real*8(A-H,O-Z)
      Dimension C(3,NAtoms+NPDir)

C     Atom I
      XI = C(1,I)
      YI = C(2,I)
      ZI = C(3,I)
C     Atom J
      XJ = C(1,J)
      YJ = C(2,J)
      ZJ = C(3,J)

C     Translational vector 1
      If(NPDir.ge.1) then
        T1X = C(1,NAtoms+1)
        T1Y = C(2,NAtoms+1)
        T1Z = C(3,NAtoms+1)
      else
        T1X = 0
        T1Y = 0
        T1Z = 0
        endif
     
      If(NPDir.ge.2) then
        T2X = C(1,NAtoms+2)
        T2Y = C(2,NAtoms+2)
        T2Z = C(3,NAtoms+2)
      else
        T2X = 0
        T2Y = 0
        T2Z = 0
        endif 

      If(NPDir.ge.3) then
        T3X = C(1,NAtoms+3)
        T3Y = C(2,NAtoms+3)
        T3Z = C(3,NAtoms+3)
      else
        T3X = 0
        T3Y = 0
        T3Z = 0
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
      Implicit Real*8 (A-H,O-Z)

      FDmp = 1.0 / (1 + Exp(-D6*(RIJ/(SR6*R0IJ) - 1)))
  
      Return
      End


*Deck FDmpDr
      Real*8 Function FDmpDr(RIJ,D6,SR6,R0IJ,F,I,J,k,C,l,m,n,IP,
     $  NAtoms,NPDir,Iwrt)
C
C     Derivative of FDmp.  If WRT = "TV", it is with respect to TVk.
C     Otherwise, it is with respect to RIk.
C
      Implicit Real*8 (A-H,O-Z)
      Dimension C(3,NAtoms+NPDir)

C      Write(*,*)'RIJ: ',RIJ,' D6: ',D6,' SR6: ',SR6
C      Write(*,*)'R0IJ: ',R0IJ,' F: ',F,' I: ',I
C      Write(*,*)'J: ',J,'k: ',k,' l: ',l
C      Write(*,*)'m: ',m,' n: ',n,'IP: ',IP

C     Atom I
      CKI = C(K,I)
C     Atom J
      CKJ = C(K,J)
     
C     Get the k component of the TVs 
C     Translational vector 1
      If(NPDir.ge.1) then
        T1K = C(K,NAtoms+1)
      else
        T1K = 0
        endif
     
      If(NPDir.ge.2) then
        T2K = C(K,NAtoms+2)
      else
        T2K = 0
        endif 

      If(NPDir.ge.3) then
        T3K = C(K,NAtoms+3)
      else
        T3K = 0
        endif 

      FTmp = D6*(CKI - CKJ - l*T1K - m*T2K - n*T3K)
      FTmp = FTmp * Exp(-D6*(RIJ/(SR6*R0IJ)-1)) * (F**2)
      FTmp = FTmp / (RIJ*(SR6*R0IJ))
C     Doing a derivative with respect to T or R?
      If(Iwrt.eq.1) then
C     Differential w.r.t. T

C     See if you multiply the front by l, m, or n.  This depends on
C     what TV you're differentiating with respect to.
        If(IP.eq.1) then
          Ia = m
        else if(IP.eq.2) then
          Ia = l
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
