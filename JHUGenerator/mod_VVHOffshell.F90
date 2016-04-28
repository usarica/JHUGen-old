module ModVVHOffshell
use ModParameters
use ModMisc
implicit none

#define fullH diag%HCurrent
#define fullV diag%VCurrent
#define fullA diag%ACurrent
#define VVHOffshell_debug 1


private

   integer, parameter :: kbsi=2,ksig=1,kbkg=0 ! These are used in V->4f functions
   real(dp), private, parameter :: tol = 0.00000010_dp
   real(dp), parameter :: T3_nulep(1:2) = (/ 0.5d0, -0.5d0 /)
   real(dp), parameter :: T3_ud(1:2) = (/ 0.5d0, -0.5d0 /)
   real(dp), parameter :: Q_nulep(1:2) = (/ 0d0, -1d0 /)
   real(dp), parameter :: Q_ud(1:2) = (/ 0.66666666666666666666666666667d0, -0.33333333333333333333333333333d0 /)
   real(dp), parameter :: xw = sitW**2
   real(dp), parameter :: cw = dsqrt(1d0-xw)
   real(dp), parameter :: eA = -dsqrt(4d0*pi*alpha_QED) ! eA=-e
   real(dp), parameter :: eW = -eA/sitW ! g=e/sintW
   real(dp), parameter :: eZ = eW*cw ! g*costw
   real(dp), parameter :: eW_cW = eW/cw ! g/costw

   ! Temporary parameters
   logical, parameter :: do1HVertex = .true. ! HVV or Hff
   logical, parameter :: do2HVertex = .true. ! HHVV
   logical, parameter :: do3HVertex = .true. ! HHH
   logical, parameter :: do4HVertex = .true. ! HHHH
   logical, parameter :: doSignal = (do1HVertex .or. do2HVertex .or. do3HVertex .or. do4HVertex)
   logical, parameter :: doBkg = .true.


   type :: ScalarCurrent
      ! Inputs to store
      real(dp) :: p(1:4,1:2)
      integer :: id(1:2)
      integer :: hel
      ! Outputs to use further
      integer :: idVertex ! =Not_a_particle_ if the vertex is strictly not valid
      real(dp) :: pVertex(1:4)
      complex(dp) :: Ub(1:4)
      complex(dp) :: V(1:4)

      ! Vector current
      complex(dp) :: scprefactor
      complex(dp) :: pscprefactor
      complex(dp) :: current
      complex(dp) :: prop ! ScalarPropagator
      complex(dp) :: propcurrent ! prop*current

      contains
         procedure :: init_ScalarCurrent
         procedure, nopass :: reset_ScalarCurrent
   end type

   type :: VACurrent
      ! Inputs to store
      real(dp) :: p(1:4,1:2)
      integer :: id(1:2)
      integer :: hel
      ! Outputs to use further
      integer :: idVertex ! =Not_a_particle_ if the vertex is strictly not valid, =Pho_ if onshell photon
      integer :: idVertex_flavorviolating ! =Not_a_particle_ unless a W(~Z/~A)
      logical :: isOnshellBoson
      real(dp) :: pVertex(1:4)
      complex(dp) :: Ub(1:4)
      complex(dp) :: V(1:4)

      ! Vector current
      complex(dp) :: prefactor
      complex(dp) :: current(1:4)
      complex(dp) :: prop ! ScalarPropagator
      complex(dp) :: propcurrent(1:4) ! VertexPropagator

      contains
         procedure :: init_VACurrent
         procedure, nopass :: reset_VACurrent
   end type

   type :: Single4VDiagram
      type(ScalarCurrent) :: HCurrent(1:4)
      type(VACurrent) :: VCurrent(1:4)
      type(VACurrent) :: ACurrent(1:4)
      real(dp) :: p(1:4,1:2,1:4)
      integer :: id(1:2,1:4)
      integer :: hel(1:4)
      real(dp) :: permutation_factor
      complex(dp) :: Amplitude
      logical :: isValid

      contains
         procedure :: init_Single4VDiagram
         procedure, nopass :: reset_Single4VDiagram
         procedure, nopass :: computeContributions

         procedure, nopass :: amp_VVHVV ! Signal VV.H.VV

         procedure, nopass :: amp_WWZA ! Background WW.ZA+WW.AZ
         procedure, nopass :: amp_WWWW ! Background WW.WW+W1W2.Z.W3W4+W1W2.A.W3W4+W1W4.Z.W3W2+W1W4.A.W3W2
         procedure, nopass :: amp_VV_VtVp
         procedure, nopass :: amp_tVVV
   end type

   type :: ProcessTree
      real(dp) :: mesq
      type(Single4VDiagram), dimension(:), allocatable :: Diagrams
      integer, dimension(:,:), allocatable :: fbperm
      integer :: nCombinations
      contains
         procedure :: init_ProcessTree
         procedure, nopass :: reset_ProcessTree
         procedure, nopass :: compute_ProcessTree
   end type

!----- List of  subroutines
   public :: amp_WWZA, amp_WWWW, amp_VV_VtVp, amp_tVVV
   public :: EvalAmp_VVHOffshell

contains

SUBROUTINE NEXPER(N, A, MTC, EVEN)
   INTEGER I,J,N,M
   INTEGER, DIMENSION(N) :: A
   INTEGER S, D, NM3, IA, I1, L
   LOGICAL MTC, EVEN

   IF (MTC) GOTO 10
   NM3 = N-3
   DO 1 I=1,N
1     A(I)=I
   MTC=.TRUE.
5  EVEN=.TRUE.
   IF(N .EQ. 1) GOTO 8
6  IF(A(N) .NE. 1 .OR. A(1) .NE. 2+MOD(N,2)) RETURN
   IF(N .LE. 3) GOTO 8
   DO 7 I=1,NM3
      IF(A(I+1) .NE. A(I)+1) RETURN
7     CONTINUE
8  MTC=.FALSE.
   RETURN
10 IF(N .EQ. 1) GOTO 27
   IF(.NOT. EVEN) GOTO 20
   IA=A(1)
   A(1)=A(2)
   A(2)=IA
   EVEN=.FALSE.
   GOTO 6
20 S=0
   DO 26 I1=2,N
25    IA=A(I1)
      I=I1-1
      D=0
      DO 30 J=1,I
30       IF(A(J) .GT. IA) D=D+1
      S=D+S
      IF(D .NE. I*MOD(S,2)) GOTO 35
26    CONTINUE
27 A(1)=0
   GOTO 8
35 M=MOD(S+1,2)*(N+1)
   DO 40 J=1,I
      IF(ISIGN(1,A(J)-IA) .EQ. ISIGN(1,A(J)-M)) GOTO 40
      M=A(J)
      L=J
40    CONTINUE
   A(L)=IA
   A(I1)=M
   EVEN=.TRUE.
   RETURN
END SUBROUTINE NEXPER


subroutine init_ScalarCurrent(cur, p, id, hel, Ub, V)
   implicit none
   class(ScalarCurrent) :: cur
   real(dp), intent(in) :: p(1:4,1:2)
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   complex(dp), intent(in) :: Ub(1:4)
   complex(dp), intent(in) :: V(1:4)
   integer :: idc(1:2)
   logical :: isValid

   isValid=.false.
   idc(1)=convertLHE(id(1))
   idc(2)=convertLHE(id(2))

   if(idc(1).gt.0 .and. idc(1).eq.-idc(2)) then ! f-fb
      cur%p=p
      cur%id=id
      isValid=.true.
   else if(idc(1).lt.0 .and. idc(1).eq.-idc(2)) then
      cur%p(:,1)=p(:,2)
      cur%id(1)=id(2)
      cur%p(:,2)=p(:,1)
      cur%id(2)=id(1)
      isValid=.true.
   endif
   if(isValid) then
      cur%hel=hel
      cur%Ub(:)=Ub
      cur%V(:)=V
      cur%pVertex(1:4)=p(1:4,1)+p(1:4,2)
      cur%idVertex=Hig_
      cur%prop = ScalarPropagator(cur%idVertex,cur%pVertex)
      cur%current = HFFCurrent(cur%id,cur%Ub,cur%V,prefactor_s_out=cur%scprefactor,prefactor_ps_out=cur%pscprefactor)
      cur%propcurrent = cur%current*cur%prop
      if(cur%scprefactor.eq.czero .and. cur%pscprefactor.eq.czero) then
         call cur%reset_ScalarCurrent(cur)
      endif
   else
      call cur%reset_ScalarCurrent(cur)
   endif
   return
end subroutine init_ScalarCurrent

subroutine reset_ScalarCurrent(cur)
   implicit none
   class(ScalarCurrent) :: cur

   cur%p(:,:)=0d0
   cur%id(:)=Not_a_particle_
   cur%hel=0
   cur%idVertex=Not_a_particle_
   cur%pVertex(:)=0d0
   cur%Ub(:)=czero
   cur%V(:)=czero
   cur%scprefactor=czero
   cur%pscprefactor=czero
   cur%current=czero
   cur%prop=czero
   cur%propcurrent=czero

   return
end subroutine reset_ScalarCurrent


subroutine init_VACurrent(cur, p, id, hel, useA)
   implicit none
   class(VACurrent) :: cur
   real(dp), intent(in) :: p(1:4,1:2)
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   logical, intent(in) :: useA
   integer :: idc(1:2)
   logical :: isValid

   isValid=.false.
   idc(1)=convertLHE(id(1))
   idc(2)=convertLHE(id(2))

   if(idc(1).gt.0 .and. idc(2).lt.0) then ! f-fb
      cur%p=p
      cur%id=id
      isValid=.true.
   else if(idc(1).lt.0 .and. idc(2).gt.0) then
      cur%p(:,1)=p(:,2)
      cur%id(1)=id(2)
      cur%p(:,2)=p(:,1)
      cur%id(2)=id(1)
      isValid=.true.
   endif
   if(isValid) then
      cur%hel=hel
      cur%pVertex(1:4)=p(1:4,1)+p(1:4,2)
      cur%current = Vcurrent(cur%p,cur%id,cur%hel,cur%idVertex,idV_flavorviolating=cur%idVertex_flavorviolating,useAcoupl=useA,Ub_out=cur%Ub,V_out=cur%V,prefactor_out=cur%prefactor, isOnshellV=cur%isOnshellBoson)
      if(cur%isOnshellBoson) then
         cur%current = pol_massless(cur%pVertex,cur%hel) ! current needs to be re-calculated
         cur%current(:) = conjg(cur%current(:)) ! Make polarization outgoing
         cur%prop = cone
         cur%propcurrent = cur%current
      else
         cur%propcurrent = VectorPropagator(cur%idVertex,cur%pVertex,cur%current,scprop=cur%prop)

         if(cur%idVertex.eq.Pho_ .and. .not.includeGammaStar) then ! Already passes off-shellness; make sure off-shell photons are not present if includeGammaStar is false.
            call reset_VACurrent(cur)
         endif
      endif

      if(cur%idVertex_flavorviolating.eq.Not_a_particle_ .and. cur%prefactor.eq.czero) then ! If it is flavor-violating, prefactor would be 0. Otherwise, it would be 1 for on-shell and something else for off-shell.
         call cur%reset_VACurrent(cur)
      endif
   else
      call cur%reset_VACurrent(cur)
   endif
   return
end subroutine init_VACurrent

subroutine reset_VACurrent(cur)
   implicit none
   class(VACurrent) :: cur

   cur%p(:,:)=0d0
   cur%id(:)=Not_a_particle_
   cur%hel=0
   cur%idVertex=Not_a_particle_
   cur%idVertex_flavorviolating=Not_a_particle_
   cur%pVertex(:)=0d0
   cur%Ub(:)=czero
   cur%V(:)=czero
   cur%prefactor=czero
   cur%current(:)=czero
   cur%prop=czero
   cur%propcurrent(:)=czero

   return
end subroutine reset_VACurrent


subroutine init_Single4VDiagram(diag, ptmp, idtmp, hel, even)
   implicit none
   class(Single4VDiagram) :: diag
   real(dp), intent(in) :: ptmp(1:4,1:8)
   integer, intent(in) :: idtmp(1:8)
   integer, intent(in) :: hel(1:4)
   logical, intent(in) :: even
   real(dp) :: pin(1:4,1:2)
   integer :: idin(1:2)
   integer :: iV, mu, ipair(1:2),j
   logical :: isValid

   diag%Amplitude=czero
   diag%isValid=.true.

   do iV=1,4
      diag%hel(iV) = hel(iV)

      ipair(1) = iV*2-1
      ipair(2) = iV*2
      do j=1,2
         pin(:,j)=ptmp(:,ipair(j))
         idin(j)=idtmp(ipair(j))
      enddo
      call diag%VCurrent(iV)%init_VACurrent(pin,idin,hel(iV),.false.)
      call diag%ACurrent(iV)%init_VACurrent(pin,idin,hel(iV),.true.)
      if( & ! Try gamma* first
         .not.diag%ACurrent(iV)%isOnshellBoson .and.                                                                        &
         (diag%ACurrent(iV)%idVertex.ne.Not_a_particle_ .or. diag%ACurrent(iV)%idVertex_flavorviolating.ne.Not_a_particle_) &
         ) then
         call diag%HCurrent(iV)%init_ScalarCurrent(pin,idin,hel(iV),diag%ACurrent(iV)%Ub, diag%ACurrent(iV)%V)
      else ! It is possible that diag%VCurrent(iV)%idVertex.ne.Z0_, but diag%HCurrent(iV) would then be initiated to a 0 current. Better to not leave it uninitialized.
         call diag%HCurrent(iV)%init_ScalarCurrent(pin,idin,hel(iV),diag%VCurrent(iV)%Ub, diag%VCurrent(iV)%V)
      endif
      if( &
         fullH(iV)%idVertex.eq.Not_a_particle_ .and.                 &
         fullV(iV)%idVertex.eq.Not_a_particle_ .and.                 &
         fullV(iV)%idVertex_flavorviolating.eq.Not_a_particle_ .and. &
         fullA(iV)%idVertex.eq.Not_a_particle_ .and.                 &
         fullA(iV)%idVertex_flavorviolating.eq.Not_a_particle_       &
         ) then
         diag%isValid=.false.
         exit
      endif
   enddo

   if(diag%isValid) then
      if (even) then
         diag%permutation_factor=1d0
      else
         diag%permutation_factor=-1d0
      endif
      call diag%computeContributions(diag)
   else
#if VVHOffshell_debug==1
      print *,"init_Single4VDiagram :: Diagram is not valid; resetting it."
#endif
      call diag%reset_Single4VDiagram(diag)
   endif
   return
end subroutine init_Single4VDiagram

subroutine reset_Single4VDiagram(diag)
   implicit none
   class(Single4VDiagram) :: diag
   integer :: iV

   do iV=1,4
      call fullH(iV)%reset_ScalarCurrent(fullH(iV))
      call fullA(iV)%reset_VACurrent(fullA(iV))
      call fullV(iV)%reset_VACurrent(fullV(iV))
      diag%p(:,:,iV)=0_dp
      diag%id(:,iV)=Not_a_particle_
      diag%hel(iV)=0
   enddo
   diag%permutation_factor=0_dp
   diag%Amplitude=czero
   diag%isValid=.false.

   return
end subroutine reset_Single4VDiagram

subroutine computeContributions(diag)
   implicit none
   class(Single4VDiagram) :: diag

   call diag%amp_VVHVV(diag)
   call diag%amp_WWZA(diag)
   call diag%amp_WWWW(diag)
   call diag%amp_VV_VtVp(diag)
   call diag%amp_tVVV(diag)

   return
end subroutine


! MomExt, MY_IDUP should be in the outgoing-momenta conventions
subroutine init_ProcessTree(tree,MomExt,MY_IDUP,hel,npart)
   implicit none
   class(ProcessTree) :: tree
   integer, parameter :: nhel = 4
   integer, intent(in) :: npart
   real(dp), intent(in) :: MomExt(1:4,1:npart) ! Mandatory input is a fermion-antifermion pair! Does not have to be ordered
   integer, intent(in) :: MY_IDUP(1:npart) ! Mandatory input is a fermion-antifermion pair! Does not have to be ordered
   integer, intent(in) :: hel(1:nhel)
   integer :: idf(1:4),idfb(1:4),idb(1:3)
   real(dp) :: pf(1:4,1:4),pfb(1:4,1:4),pb(1:4,1:3)
   integer :: ipart, ix, iy, ib, ic, nferm, nfermb, nbos, idV_tmp, ncomb
   real(dp) :: ptmp(1:4,1:8)
   integer :: idtmp(1:8)
   logical :: combPass
   integer, dimension(:), allocatable :: fermb_pairing ! These are input to get the group of permutations
   logical :: mtc ! These are input to get the group of permutations
   logical :: even ! These are input to get the group of permutations

   ! Convention for helicity assignment: Assign to f and fb in the order passed, then to the bosons

   tree%mesq = 0_dp

   idf(:)=Not_a_particle_
   idfb(:)=Not_a_particle_
   idb(:)=Not_a_particle_
   pf(:,:)=0_dp
   pfb(:,:)=0_dp
   pb(:,:)=0_dp
   nferm=1
   nfermb=1
   nbos=1

   if(npart.lt.5 .or. npart.gt.8) return ! At least a qqb->GGG and at most a 2q->6f state

   ! Group the ids
   do ipart=1,npart
      if(MY_IDUP(ipart).eq.Not_a_particle_) then
         continue
      else if(MY_IDUP(ipart).eq.Pho_) then ! Massless bosons
         if(nbos.eq.4) then ! This should never happen!
            print *,"init_ProcessTree :: Error: Number of bosons hit 4. Resetting the ProcessTree and aborting init_ProcessTree"
            call tree%reset_ProcessTree(tree)
            return
         endif
         idb(nbos)=MY_IDUP(ipart)
         pb(:,nbos)=MomExt(:,ipart)
         nbos=nbos+1
      else if( convertLHE(MY_IDUP(ipart)).gt.0 ) then
         idf(nferm)=MY_IDUP(ipart)
         pf(:,nferm)=MomExt(:,ipart)
         nferm=nferm+1
      else if( convertLHE(MY_IDUP(ipart)).lt.0 ) then
         idfb(nfermb)=MY_IDUP(ipart)
         pf(:,nfermb)=MomExt(:,ipart)
         nfermb=nfermb+1
      endif
   enddo
   if(nferm.ne.nfermb) return ! f-fb=0 has to be satisfied in the case where everything is outgoing
   if((nferm+nbos).ne.nhel) return ! We need the helicities to be defined
   allocate(fermb_pairing(nfermb))

   ! Count the combinations
   ncomb=0
   do ix=1,nfermb ! Identity
      fermb_pairing(ix)=ix
   enddo
   mtc = .false. ! Start with this statement
50 continue
   call NEXPER(nfermb, fermb_pairing(:), mtc, even)
   combPass=.true.
   do ix=1,nferm
      iy = fermb_pairing(ix)
      if(hel(ix).eq.hel(iy)) then
         idV_tmp=CoupledVertex((/idf(ix),idfb(iy)/),hel(ix))
         if(idV_tmp.eq.Not_a_particle_) idV_tmp = CoupledVertex_FlavorViolating((/idf(ix),idfb(iy)/),hel(ix))
         if(idV_tmp.eq.Not_a_particle_) combPass=.false.
      endif
      if(.not.combPass) exit ! No need to evaluate the others. The particular permutation is invalid.
#if VVHOffshell_debug==1
      print *,"init_ProcessTree :: Doing permutation (",ncomb,", ",iy,")"
#endif
   enddo
   if(combPass) then
      ncomb = ncomb+1
   endif
   if(mtc) GOTO 50

   if(ncomb.eq.0) then
      deallocate(fermb_pairing)
      return
   endif

   tree%nCombinations=ncomb
   allocate(tree%Diagrams(ncomb))
   allocate(tree%fbperm(4,ncomb))

   ! Start over, this time to initialize the combinations
   ic=0 ! Track the combination index
   do ix=1,nfermb ! First time fill, identity
      fermb_pairing(ix)=ix
   enddo
   mtc = .false. ! Start with this statement
   even = .true.
51 continue
   call NEXPER(nfermb, fermb_pairing(:), mtc, even)
   combPass=.true.
   do ix=1,nferm
      iy = fermb_pairing(ix)
      if(hel(ix).eq.hel(iy)) then
         idV_tmp=CoupledVertex((/idf(ix),idfb(iy)/),hel(ix))
         if(idV_tmp.eq.Not_a_particle_) idV_tmp = CoupledVertex_FlavorViolating((/idf(ix),idfb(iy)/),hel(ix))
         if(idV_tmp.eq.Not_a_particle_) combPass=.false.
      endif
   enddo
   if(combPass) then
      ic = ic + 1
      tree%fbperm(:,ic)=0 ! Positive only for fermion pairs
      do ix=1,nferm
         iy = fermb_pairing(ix)
         tree%fbperm(ix,ic)=iy
         ptmp(:,2*ix-1)=pf(:,ix)
         ptmp(:,2*ix)=pfb(:,iy)
         idtmp(2*ix-1)=idf(ix)
         idtmp(2*ix)=idfb(iy)
      enddo
      do ib=1,nbos
         ptmp(:,2*nferm+2*ib-1)=pb(:,ib)
         idtmp(2*nferm+2*ib-1)=idb(ib)
         ptmp(:,2*nferm+2*ib)=0d0
         idtmp(2*nferm+2*ib)=Not_a_particle_
      enddo
      call tree%Diagrams(ic)%init_Single4VDiagram(ptmp, idtmp, hel, even)
   endif
   if(mtc) GOTO 51

   deallocate(fermb_pairing) ! No longer needed

   if(ncomb .ne. ic) then
      print *,"init_ProcessTree :: ncomb .ne. ic!"
      call tree%reset_ProcessTree(tree)
      print *,"init_ProcessTree :: Diagrams are reset."
      stop
   endif
   return
end subroutine init_ProcessTree

subroutine reset_ProcessTree(tree)
   implicit none
   class(ProcessTree) :: tree

   tree%mesq=0
   if(allocated(tree%Diagrams)) then
      deallocate(tree%Diagrams)
   endif
   if(allocated(tree%fbperm)) then
      deallocate(tree%fbperm)
   endif
   tree%nCombinations=0

   return
end subroutine reset_ProcessTree

subroutine compute_ProcessTree(tree)
   implicit none
   class(ProcessTree) :: tree
   integer :: idiag,jdiag

   do idiag=1,tree%nCombinations
      do jdiag=idiag,tree%nCombinations
         call interfereDiagrams(tree,idiag,jdiag)
      enddo
   enddo

   return
end subroutine compute_ProcessTree

subroutine interfereDiagrams(tree,idiag,jdiag)
   implicit none
   class(ProcessTree) :: tree
   integer, intent(in) :: idiag,jdiag
   ! The following trackUnique_* variables track the unique
   logical :: isChecked(1:4)
   integer :: trackUnique_idf(1:4) ! Fermion id
   integer :: trackUnique_noccurance(1:4) ! Number of occurance of the particular current type
   integer :: tracker ! Iterator/number of unique currents
   ! End tracking variables
   integer :: iV, it1, it2
   real(dp) :: nquarksfactor
   real(dp) :: symmetryFactor
   real(dp) :: metmp

   if(idiag.gt.jdiag) then
      return
   endif

   nquarksfactor=1_dp
   symmetryFactor=1_dp

   isChecked(:)=.false.
   trackUnique_idf(:)=Not_a_particle_
   trackUnique_noccurance(:)=0
   tracker=0

   ! Count classes of swaps and how many leaps to make in each
   do iV=1,4
      if(isChecked(iV)) then
         continue
      else if(tree%fbperm(iV,idiag).eq.0) then
         isChecked(iV) = .true.
         continue
      else if(tree%fbperm(iV,idiag).gt.0) then
         isChecked(iV) = .true.

         tracker = tracker+1
         trackUnique_idf(tracker) = tree%Diagrams(idiag)%id(1,iV)
         trackUnique_noccurance(tracker) = trackUnique_noccurance(tracker)+1

         it1 = -1
         it2 = iV
107      continue
         call getIndexLocation(4, tree%fbperm(:,idiag), tree%fbperm(it2,jdiag), it1)
         if(it1.lt.1) then
            print *,"it1.lt.1. Aborting!"
            call abort
         endif
         if(it1.ne.iV) then
            it2=it1
            isChecked(it1) = .true.
            trackUnique_noccurance(tracker) = trackUnique_noccurance(tracker)+1
            goto 107
         endif
      endif
   enddo

   do iV=1,4
      if(trackUnique_idf(iV).eq.Not_a_particle_ .or. trackUnique_noccurance(iV).eq.0) then
         continue
      else if(IsAQuark(trackUnique_idf(iV)) .and. trackUnique_noccurance(iV).eq.1) then
         nquarksfactor = nquarksfactor*9_dp
      else if(IsAQuark(trackUnique_idf(iV)) .and. trackUnique_noccurance(iV).gt.1) then
         nquarksfactor = nquarksfactor*3_dp
         symmetryFactor = real(factorial(trackUnique_noccurance(iV)), kind=dp)
      else if(.not.IsAQuark(trackUnique_idf(iV)) .and. trackUnique_noccurance(iV).gt.1) then
         symmetryFactor = real(factorial(trackUnique_noccurance(iV)), kind=dp)
      endif
   enddo

   if(idiag.eq.jdiag) then
      metmp = abs((tree%Diagrams(idiag)%Amplitude)**2)
   else ! Notice the factor of 2!
      metmp = two * real(tree%Diagrams(idiag)%Amplitude*conjg(tree%Diagrams(jdiag)%Amplitude),kind=dp)
   endif
   metmp = metmp*nquarksfactor/symmetryFactor

   tree%mesq = tree%mesq + metmp
   return
end subroutine


subroutine getIndexLocation(nArray, checkArray, checkValue, iterator)
   integer, intent(in) :: nArray
   integer, intent(in) :: checkArray(1:nArray)
   integer, intent(in) :: checkValue
   integer :: iterator
   integer :: iV

   do iV=1,nArray
      if(checkArray(iV).eq.checkValue) then
         iterator = iV
         return
      endif
   enddo

   iterator = -1
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main subroutine: EvalAmp_VVHOffshell
! MomExt, MY_IDUP should be in the lab conventions. The first two always refer to the incoming particles, so their ids and momenta acquire a (-) sign.
subroutine EvalAmp_VVHOffshell(MomExt,MY_IDUP,MY_HEL,npart,mesq)
   integer, parameter :: nhel = 4
   type(ProcessTree) :: tree
   integer, intent(in) :: npart
   real(dp), intent(in) :: MomExt(1:4,1:npart)
   integer, intent(in) :: MY_IDUP(1:npart)
   integer, intent(in) :: MY_HEL(1:nhel)
   real(dp), intent(out) :: mesq
   real(dp) :: pin(1:4,1:npart)
   integer :: idin(1:npart)
   integer :: idiag, jdiag

   if(npart.lt.3) then
      return
   endif

   pin(:,1:2) = -MomExt(:,1:2); idin(1:2) = -MY_IDUP(1:2)
   pin(:,3:npart) = MomExt(:,3:npart); idin(3:npart) = -MY_IDUP(3:npart)

   ! The following directive initiates the process tree and its diagrams.
   ! init_Single4VDiagrams also calculates the contributions if the 4V diagram is valid.
   call tree%init_ProcessTree(pin,idin,MY_HEL,npart)
   ! Back to this subroutine to interfere the diagrams
   call tree%compute_ProcessTree(tree)
   ! Assign mesq
   mesq = tree%mesq
   ! Clean the process tree allocated variables
   call tree%reset_ProcessTree(tree)

   ! Done
   return
end subroutine EvalAmp_VVHOffshell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function factorial(factor)
   implicit none
   integer, intent(in) :: factor
   integer :: factorial
   integer :: j

   if(factor.lt.0) then
      factorial=0
   else
      factorial=1
      if(factor.gt.1) then
         do j=2,factor
            factorial = factorial*j
         enddo
      endif
   endif
   return
end function factorial

! Use Weyl basis to make Psi=(Psi_R,Psi_L)
! satisfying P_R,L=(1+-gamma5)/2
function gammaMatrix()
   implicit none
   complex(dp) :: gammaMatrix(1:4,1:4,1:4)

   gammaMatrix(:,:,:)=czero

   gammaMatrix(1,1,3)=1_dp ! gamma^0
   gammaMatrix(1,2,4)=1_dp
   gammaMatrix(1,3,1)=1_dp
   gammaMatrix(1,4,2)=1_dp
   gammaMatrix(2,1,4)=-1_dp ! gamma^x
   gammaMatrix(2,2,3)=-1_dp
   gammaMatrix(2,3,2)=1_dp
   gammaMatrix(2,4,1)=1_dp
   gammaMatrix(3,1,4)=ci ! gamma^y
   gammaMatrix(3,2,3)=-ci
   gammaMatrix(3,3,2)=-ci
   gammaMatrix(3,4,1)=ci
   gammaMatrix(4,1,3)=-1_dp ! gamma^z
   gammaMatrix(4,2,4)=1_dp
   gammaMatrix(4,3,1)=1_dp
   gammaMatrix(4,4,2)=-1_dp
   return
end function gammaMatrix
function gammaFive(scl)
   implicit none
   complex(dp) :: gammaFive(1:4,1:4)
   complex(dp), optional :: scl

   gammaFive(:,:)=czero

   if( present(scl) ) then
      gammaFive(1,1)=scl
      gammaFive(2,2)=scl
      gammaFive(3,3)=-scl
      gammaFive(4,4)=-scl
   else
      gammaFive(1,1)=1_dp
      gammaFive(2,2)=1_dp
      gammaFive(3,3)=-1_dp
      gammaFive(4,4)=-1_dp
   endif
   return
end function gammaFive
function gammaMatrixDotted(p)
   implicit none
   complex(dp) :: gammaMatrixDotted(1:4,1:4)
   complex(dp), intent(in) :: p(1:4)

   gammaMatrixDotted(:,:)=czero

   ! The coefficients of p(2-4) receive a (-) sign wrt. gammaMatrix() bc. of g_munu.
   gammaMatrixDotted(1,3)=p(1)+p(4)
   gammaMatrixDotted(2,4)=p(1)-p(4)
   gammaMatrixDotted(3,1)=p(1)-p(4)
   gammaMatrixDotted(4,2)=p(1)+p(4)

   gammaMatrixDotted(1,4)=p(2)-ci*p(3)
   gammaMatrixDotted(2,3)=p(2)+ci*p(3)
   gammaMatrixDotted(3,2)=-p(2)+ci*p(3)
   gammaMatrixDotted(4,1)=-p(2)-ci*p(3)

   return
end function gammaMatrixDotted
function identityMatrix(scl)
   implicit none
   complex(dp) :: identityMatrix(1:4,1:4)
   complex(dp), optional :: scl

   identityMatrix(:,:)=0_dp

   if( present(scl) ) then
      identityMatrix(1,1)=scl
      identityMatrix(2,2)=scl
      identityMatrix(3,3)=scl
      identityMatrix(4,4)=scl
   else
      identityMatrix(1,1)=1_dp
      identityMatrix(2,2)=1_dp
      identityMatrix(3,3)=1_dp
      identityMatrix(4,4)=1_dp
   endif
   return
end function identityMatrix

function MultiplySquareMatrices(ND,m1,m2)
   implicit none
   integer, intent(in) :: ND
   complex(dp), intent(in) :: m1(1:ND,1:ND),m2(1:ND,1:ND)
   complex(dp) :: MultiplySquareMatrices(1:ND,1:ND)
   integer :: i,a,b

   MultiplySquareMatrices=czero
   if(ND.lt.1) return
   do a=1,ND
   do b=1,ND
   do i=1,ND
      MultiplySquareMatrices(a,b) = MultiplySquareMatrices(a,b) + m1(a,i)*m2(i,b)
   enddo
   enddo
   enddo

end function MultiplySquareMatrices

function gamma_gamma_mu(N,p,ubar,v,whichUndotted,massMatrix,gammaFiveMatrix)
   implicit none
   integer, intent(in) :: N
   complex(dp), intent(in) :: p(1:4,2:N)
   complex(dp), intent(in) :: ubar(1:4),v(1:4)
   integer, optional :: whichUndotted
   complex(dp), optional :: massMatrix(2:N)
   complex(dp), optional :: gammaFiveMatrix(2:N)
   complex(dp) :: gamma_gamma_mu(1:4)
   integer :: iUndotted,a,b,i,j,mu
   complex(dp) :: gamma(1:4,1:4,1:4),gamma_tmp(1:4,1:4,2:N),gtmp(1:4,1:4,1:4),gacc(1:4,1:4,1:4)

   gamma_gamma_mu(:)=czero
   if(N.lt.2) return
   iUndotted=1
   if( present(whichUndotted) ) then
      if(whichUndotted.le.N .and. whichUndotted.gt.0) then
         iUndotted=whichUndotted
      endif
   endif

   gamma(:,:,:)=gammaMatrix()

   ! Now construct gamma_slash
   j=2
   do i=1,N
      if(i.ne.iUndotted) then
         gamma_tmp(:,:,j) = gammaMatrixDotted(p(:,j)) ! gamma(1,:,:)*p(1,j) - gamma(2,:,:)*p(2,j) - gamma(3,:,:)*p(3,j) - gamma(4,:,:)*p(4,j) ! (gamma^mu p_i_mu)_ab
         if( present(massMatrix) ) then
            if(massMatrix(j).ne.czero) then
               gamma_tmp(:,:,j) = gamma_tmp(:,:,i) + identityMatrix(massMatrix(j))
            endif
         endif
         if( present(gammaFiveMatrix) ) then
            if(gammaFiveMatrix(j).ne.czero) then
               gamma_tmp(:,:,j) = gamma_tmp(:,:,i) + gammaFive(gammaFiveMatrix(j))
            endif
         endif
         j=j+1
      else
         continue
      endif
   enddo

   gacc=gamma(:,:,:)
   if(iUndotted.lt.N) then
      do i=iUndotted+1,N
         do mu=1,4
            gtmp(mu,:,:) = MultiplySquareMatrices(4,gacc(mu,:,:),gamma_tmp(:,:,i)) ! (gamma_j^mu)_ac*( Prod_{i=j+1..N} (gamma_i^nu p_i_nu) )_cb
            gacc=gtmp
         enddo
      enddo
   endif
   if(iUndotted.gt.1) then
      do i=2,iUndotted
         do mu=1,4
            gtmp(mu,:,:) = MultiplySquareMatrices(4,gamma_tmp(:,:,i),gacc(mu,:,:)) ! ( Prod_{i=2..j} (gamma_i^nu p_i_nu) )_ac * (gamma_j^mu)_cb
            gacc=gtmp
         enddo
      enddo
   endif

   do a=1,4
      do b=1,4
         do mu=1,4
            gamma_gamma_mu(mu) = gamma_gamma_mu(mu) + ubar(a)*gacc(mu,a,b)*v(b) ! Returns a 4-vector in the end
         enddo
      enddo
   enddo
   return
end function gamma_gamma_mu

function gamma_dotted(N,p,massMatrix,gammaFiveMatrix)
   implicit none
   integer, intent(in) :: N
   complex(dp), intent(in) :: p(1:4,1:N)
   complex(dp), optional :: massMatrix(1:N)
   complex(dp), optional :: gammaFiveMatrix(1:N)
   complex(dp) :: gamma_dotted(1:4,1:4)
   integer :: a,b,i
   complex(dp) :: gamma_tmp(1:4,1:4,1:N)

   gamma_dotted=czero
   if(N.lt.1) then
      gamma_dotted=identityMatrix()
      return
   endif

   ! Now construct gamma_slash
   do i=1,N
      gamma_tmp(:,:,i) = gammaMatrixDotted(p(:,i)) ! gamma(1,:,:)*p(1,j) - gamma(2,:,:)*p(2,j) - gamma(3,:,:)*p(3,j) - gamma(4,:,:)*p(4,j) ! (gamma^mu p_i_mu)_ab
      if( present(massMatrix) ) then
         if(massMatrix(i).ne.czero) then
            gamma_tmp(:,:,i) = gamma_tmp(:,:,i) + identityMatrix(massMatrix(i))
         endif
      endif
      if( present(gammaFiveMatrix) ) then
         if(gammaFiveMatrix(i).ne.czero) then
            gamma_tmp(:,:,i) = gamma_tmp(:,:,i) + gammaFive(gammaFiveMatrix(i))
         endif
      endif
   enddo

   gamma_dotted=gamma_tmp(:,:,1)
   if(N.gt.1) then
      do i=2,N
         gamma_dotted(:,:) = MultiplySquareMatrices(4,gamma_dotted(:,:),gamma_tmp(:,:,i))
      enddo
   endif
   return
end function gamma_dotted

function gamma_dotted_contracted(N,p,ubar,v,massM,gFiveM)
   implicit none
   integer, intent(in) :: N
   complex(dp), intent(in) :: p(1:4,1:N)
   complex(dp), intent(in) :: ubar(1:4),v(1:4)
   complex(dp), optional :: massM(1:N)
   complex(dp), optional :: gFiveM(1:N)
   complex(dp) :: gamma_dotted_contracted
   integer :: a,b,i
   complex(dp) :: gamma_tmp(1:4,1:4,1:N),gacc(1:4,1:4)

   gamma_dotted_contracted=czero
   if( present(massM) .and. present(gFiveM) ) then
      gacc=gamma_dotted(N,p,massMatrix=massM,gammaFiveMatrix=gFiveM)
   else if( present(gFiveM) ) then
      gacc=gamma_dotted(N,p,gammaFiveMatrix=gFiveM)
   else if( present(massM) ) then
      gacc=gamma_dotted(N,p,massMatrix=massM)
   else
      gacc=gamma_dotted(N,p)
   endif

   do a=1,4
      do b=1,4
         gamma_dotted_contracted = gamma_dotted_contracted + ubar(a)*gacc(a,b)*v(b)
      enddo
   enddo
   return
end function gamma_dotted_contracted


function CoupledVertex(id,hel,useAHcoupl)
   implicit none
   integer, optional :: useAHcoupl
   integer, intent(in) :: id(1:2),hel
   integer :: testAHcoupl
   integer :: CoupledVertex

   testAHcoupl = 0
   if(present(useAHcoupl)) then
      testAHcoupl = useAHcoupl
   endif
   if( (&
   (id(1).eq.ElP_ .and. id(2).eq.NuE_) .or. (id(2).eq.ElP_ .and. id(1).eq.NuE_) .or. &
   (id(1).eq.MuP_ .and. id(2).eq.NuM_) .or. (id(2).eq.MuP_ .and. id(1).eq.NuM_) .or. &
   (id(1).eq.TaP_ .and. id(2).eq.NuT_) .or. (id(2).eq.TaP_ .and. id(1).eq.NuT_) .or. &
   (id(1).eq.Up_  .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Up_  .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Chm_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Chm_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_)) .or. &
   (id(1).eq.Top_ .and. (id(2).eq.ADn_ .or. id(2).eq.AStr_ .or. id(2).eq.ABot_)) .or. (id(2).eq.Top_ .and. (id(1).eq.ADn_ .or. id(1).eq.AStr_ .or. id(1).eq.ABot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wp_
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ANuE_) .or. (id(2).eq.ElM_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.MuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.ANuT_) .or. (id(2).eq.TaM_ .and. id(1).eq.ANuT_) .or. &
   (id(1).eq.AUp_  .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AUp_  .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.AChm_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.AChm_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_)) .or. &
   (id(1).eq.ATop_ .and. (id(2).eq.Dn_ .or. id(2).eq.Str_ .or. id(2).eq.Bot_)) .or. (id(2).eq.ATop_ .and. (id(1).eq.Dn_ .or. id(1).eq.Str_ .or. id(1).eq.Bot_))      &
   ) .and. hel.lt.0) then
      CoupledVertex=Wm_
   elseif( (&
   (id(1).eq.ElM_ .and. id(2).eq.ElP_) .or. (id(2).eq.ElM_ .and. id(1).eq.ElP_) .or. &
   (id(1).eq.MuM_ .and. id(2).eq.MuP_) .or. (id(2).eq.MuM_ .and. id(1).eq.MuP_) .or. &
   (id(1).eq.TaM_ .and. id(2).eq.TaP_) .or. (id(2).eq.TaM_ .and. id(1).eq.TaP_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.AUp_) .or. (id(2).eq.Up_  .and. id(1).eq.AUp_) .or. &
   (id(1).eq.Dn_  .and. id(2).eq.ADn_) .or. (id(2).eq.Dn_  .and. id(1).eq.ADn_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AChm_) .or. (id(2).eq.Chm_ .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.AStr_) .or. (id(2).eq.Str_ .and. id(1).eq.Astr_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.ATop_) .or. (id(2).eq.Top_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ABot_) .or. (id(2).eq.Bot_ .and. id(1).eq.ABot_)      &
   ) .and. hel.ne.0) then
      if(testAHcoupl.eq.1) then
         CoupledVertex=Pho_
      elseif(testAHcoupl.eq.2) then
         CoupledVertex=Hig_
      else
         CoupledVertex=Z0_
      endif
   elseif( (&
   (id(1).eq.NuE_ .and. id(2).eq.ANuE_) .or. (id(2).eq.NuE_ .and. id(1).eq.ANuE_) .or. &
   (id(1).eq.NuM_ .and. id(2).eq.ANuM_) .or. (id(2).eq.NuM_ .and. id(1).eq.ANuM_) .or. &
   (id(1).eq.NuT_ .and. id(2).eq.ANuT_) .or. (id(2).eq.NuT_ .and. id(1).eq.ANuT_)      &
   ) .and. hel.lt.0) then
      CoupledVertex=Z0_ ! Only Z coupling to nuL-nubR
   else
      CoupledVertex=Not_a_particle_
   endif

   return
end function CoupledVertex

function CoupledVertex_FlavorViolating(id,hel,useAHcoupl)
   implicit none
   integer, optional :: useAHcoupl
   integer, intent(in) :: id(1:2),hel
   integer :: testAHcoupl
   integer :: CoupledVertex_FlavorViolating

   testAHcoupl = 0
   if(present(useAHcoupl)) then
      testAHcoupl = useAHcoupl
   endif
   if( (&
   (id(1).eq.Up_  .and. id(2).eq.AChm_) .or. (id(2).eq.Up_  .and. id(1).eq.AChm_) .or. &
   (id(1).eq.Up_  .and. id(2).eq.ATop_) .or. (id(2).eq.Up_  .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.AUp_ ) .or. (id(2).eq.Chm_ .and. id(1).eq.AUp_ ) .or. &
   (id(1).eq.Chm_ .and. id(2).eq.ATop_) .or. (id(2).eq.Chm_ .and. id(1).eq.ATop_) .or. &
   (id(1).eq.Top_ .and. id(2).eq.AUp_ ) .or. (id(2).eq.Top_ .and. id(1).eq.AUp_ ) .or. &
   (id(1).eq.Top_ .and. id(2).eq.AChm_) .or. (id(2).eq.Top_ .and. id(1).eq.AChm_) .or. &

   (id(1).eq.Dn_  .and. id(2).eq.AStr_) .or. (id(2).eq.Dn_  .and. id(1).eq.AStr_) .or. &
   (id(1).eq.Dn_  .and. id(2).eq.ABot_) .or. (id(2).eq.Dn_  .and. id(1).eq.ABot_) .or. &
   (id(1).eq.Str_ .and. id(2).eq.ADn_ ) .or. (id(2).eq.Str_ .and. id(1).eq.ADn_ ) .or. &
   (id(1).eq.Str_ .and. id(2).eq.ABot_) .or. (id(2).eq.Str_ .and. id(1).eq.ABot_) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.ADn_ ) .or. (id(2).eq.Bot_ .and. id(1).eq.ADn_ ) .or. &
   (id(1).eq.Bot_ .and. id(2).eq.AStr_) .or. (id(2).eq.Bot_ .and. id(1).eq.AStr_)      &
   ) .and. hel.lt.0) then ! Tests W coupling that looks like a Z decay, so only allow left-handed decays
      if(testAHcoupl.eq.1) then
         CoupledVertex_FlavorViolating=Pho_
      elseif(testAHcoupl.eq.2) then
         CoupledVertex_FlavorViolating=Hig_
      else
         CoupledVertex_FlavorViolating=Z0_
      endif
   else
      CoupledVertex_FlavorViolating=Not_a_particle_
   endif

   return
end function CoupledVertex_FlavorViolating

function WDaughterPair(id, considerTops)
   implicit none
   integer, intent(in) :: id
   integer :: WDaughterPair(1:3)
   logical, optional :: considerTops

   WDaughterPair(:)=Not_a_particle_
   if(id.eq.ElP_) then
      WDaughterPair(1)=NuE_
   elseif(id.eq.MuP_) then
      WDaughterPair(1)=NuM_
   elseif(id.eq.TaP_) then
      WDaughterPair(1)=NuT_
   elseif(id.eq.NuE_) then
      WDaughterPair(1)=ElP_
   elseif(id.eq.NuM_) then
      WDaughterPair(1)=MuP_
   elseif(id.eq.NuT_) then
      WDaughterPair(1)=TaP_

   elseif(id.eq.ElM_) then
      WDaughterPair(1)=ANuE_
   elseif(id.eq.MuM_) then
      WDaughterPair(1)=ANuM_
   elseif(id.eq.TaM_) then
      WDaughterPair(1)=ANuT_
   elseif(id.eq.ANuE_) then
      WDaughterPair(1)=ElM_
   elseif(id.eq.ANuM_) then
      WDaughterPair(1)=MuM_
   elseif(id.eq.ANuT_) then
      WDaughterPair(1)=TaM_

   elseif(id.eq.Up_ .or. id.eq.Chm_ .or. id.eq.Top_) then
      WDaughterPair(1)=ADn_
      WDaughterPair(2)=AStr_
      WDaughterPair(3)=ABot_
   elseif(id.eq.AUp_ .or. id.eq.AChm_ .or. id.eq.ATop_) then
      WDaughterPair(1)=Dn_
      WDaughterPair(2)=Str_
      WDaughterPair(3)=Bot_

   elseif(id.eq.Dn_ .or. id.eq.Str_ .or. id.eq.Bot_) then
      WDaughterPair(1)=AUp_
      WDaughterPair(2)=AChm_
      if(present(considerTops)) then
         if(considerTops) then
            WDaughterPair(3)=ATop_
         endif
      endif
   elseif(id.eq.ADn_ .or. id.eq.AStr_ .or. id.eq.ABot_) then
      WDaughterPair(1)=Up_
      WDaughterPair(2)=Chm_
      if(present(considerTops)) then
         if(considerTops) then
            WDaughterPair(3)=Top_
         endif
      endif
   endif

   return
end function WDaughterPair

function ScalarPropagator(idV,p,compmass)
   use ModMisc
   implicit none
   integer, intent(in) :: idV
   real(dp), intent(in) :: p(1:4)
   complex(dp), optional :: compmass
   complex(dp) :: ScalarPropagator
   real(dp) :: s

   ScalarPropagator = czero
   s = p(:).dot.p(:)
   if( present(compmass) ) then
      compmass=czero
   endif
   if(idV .eq. Wm_ .or. idV.eq.Wp_) then
      ScalarPropagator = -ci/( s - M_W**2 + ci*M_W*Ga_W)
      if( present(compmass) ) then
         compmass=M_W
      endif
   elseif(idV .eq. Z0_) then
      ScalarPropagator = -ci/( s - M_Z**2 + ci*M_Z*Ga_Z)
      if( present(compmass) ) then
         compmass=M_Z
      endif
   elseif(idV .eq. Pho_) then
      ScalarPropagator = -ci/s
   elseif(idV .eq. Hig_) then
      ScalarPropagator = ci/( s - M_Reso**2 + ci*M_Reso*Ga_Reso)
      if( present(compmass) ) then
         compmass=M_Reso
      endif
   elseif(idV .eq. Top_) then
      ScalarPropagator = ci/( s - M_Top**2 + ci*M_Top*Ga_Top)
      if( present(compmass) ) then
         compmass=M_Top
      endif
   elseif(idV .eq. TaM_) then
      ScalarPropagator = ci/( s - m_tau**2 + ci*m_tau*Ga_tau)
      if( present(compmass) ) then
         compmass=m_tau
      endif
   else
      ScalarPropagator = ci/s
   endif

   return
end function ScalarPropagator

function VectorPropagator(idV,p,current,scprop)
   use ModMisc
   implicit none
   integer, intent(in) :: idV
   real(dp), intent(in) :: p(1:4)
   complex(dp), intent(in) :: current(1:4)
   complex(dp), optional :: scprop
   complex(dp) :: VectorPropagator(1:4)
   complex(dp) :: prefactor

   VectorPropagator(:) = czero
   prefactor = ScalarPropagator(idV,p)
   if(present(scprop)) then
      scprop = prefactor
   endif

   if(idV .eq. Wm_ .or. idV.eq.Wp_ .or. idV .eq. Z0_) then
      VectorPropagator(:) = prefactor*(current(:)-p(:)*((cmplx(p(:),kind=dp)).dot.current(:))/(p(:).dot.p(:)))
   elseif(idV .eq. Pho_) then
      VectorPropagator(:) = prefactor*current(:)
   endif

   return
end function VectorPropagator

function EuclideanMagnitude_Complex(ND,p)
   implicit none
   integer, intent(in) :: ND
   complex(dp), intent(in) :: p(1:ND)
   real(dp) :: EuclideanMagnitude_Complex
   integer :: i

   EuclideanMagnitude_Complex=0_dp
   if(ND.gt.0) then
      do i=1,ND
         EuclideanMagnitude_Complex = EuclideanMagnitude_Complex + p(i)*dconjg(p(i))
      enddo
      EuclideanMagnitude_Complex = dsqrt(EuclideanMagnitude_Complex)
   endif
   return
end function EuclideanMagnitude_Complex

!     ubar spinor, massless
function ubar0(p,i)
   implicit none
   complex(dp), intent(in) :: p(4)
   integer, intent(in) :: i
   complex(dp) :: ubar0(4)
   complex(dp) :: fc, fc2
   real(dp) :: p0,px,py,pz

   p0=real(p(1),dp)
   px=real(p(2),dp)
   py=real(p(3),dp)
   pz=real(p(4),dp)

   fc2 = p0 + pz
   fc=sqrt(fc2)

   if (abs(fc2).gt. tol) then
      if (i.eq.1) then
        ubar0(1)=czero
        ubar0(2)=czero
        ubar0(3)=fc
        ubar0(4)=(px-ci*py)/fc
      elseif (i.eq.-1) then
        ubar0(1)=(px+ci*py)/fc
        ubar0(2)=-fc
        ubar0(3)=czero
        ubar0(4)=czero
      else
        stop 'ubar0: i out of range'
      endif
   else
      if (i.eq.1) then
        ubar0(1) = czero
        ubar0(2) = czero
        ubar0(3) = czero
        ubar0(4) = sqrt(cone*two*p0)
      elseif (i.eq.-1) then
        ubar0(1) = sqrt(cone*two*p0)
        ubar0(2) = czero
        ubar0(3) = czero
        ubar0(4) = czero
      else
        stop 'ubar0: i out of range'
      endif
   endif
   return
end function ubar0

!     ubar spinor, massive (not used)
function ubarm(cp,i)
   implicit none
   complex(dp), intent(in) :: cp(4)
   integer, intent(in) :: i
   complex(dp) :: ubarm(4)
   complex(dp) :: fc2
   real(dp) :: E,x,y,z,m,p,fac

   E=real(cp(1),dp)
   x=real(cp(2),dp)
   y=real(cp(3),dp)
   z=real(cp(4),dp)
   m=dsqrt(dabs(E**2-x**2-y**2-z**2))
   p=dsqrt(dabs(x**2+y**2+z**2))

   fc2 = p + z

   if (abs(fc2).gt. tol) then
      fac = cone/sqrt(2_dp*p*fc2)
      if (i.eq.1) then
         ubarm(1)=sqrt(E-p)*fac*fc2
         ubarm(2)=sqrt(E-p)*fac*(x-ci*y)
         ubarm(3)=sqrt(E+p)*fac*fc2
         ubarm(4)=sqrt(E+p)*fac*(x-ci*y)
      elseif (i.eq.-1) then
         ubarm(1)=sqrt(E+p)*fac*(x+ci*y)
         ubarm(2)=-sqrt(E+p)*fac*fc2
         ubarm(3)=sqrt(E-p)*fac*(x+ci*y)
         ubarm(4)=-sqrt(E-p)*fac*fc2
      else
        stop 'ubarm: i out of range'
      endif
   else
      if (i.eq.1) then
        ubarm(1) = czero
        ubarm(2) = sqrt(cone*(E-p))
        ubarm(3) = czero
        ubarm(4) = sqrt(cone*(E+p))
      elseif (i.eq.-1) then
        ubarm(1) = sqrt(cone*(E+p))
        ubarm(2) = czero
        ubarm(3) = sqrt(cone*(E-p))
        ubarm(4) = czero
      else
        stop 'ubarm: i out of range'
      endif
   endif
   return
end function ubarm

!     v spinor, massless
function v0(p,i)
   implicit none
   complex(dp), intent(in) :: p(4)
   integer, intent(in) :: i
   complex(dp) :: v0(4)
   complex(dp) :: fc2, fc
   real(dp) :: p0,px,py,pz

   p0=real(p(1),dp)
   px=real(p(2),dp)
   py=real(p(3),dp)
   pz=real(p(4),dp)

   fc2 = p0 + pz
   fc=sqrt(fc2)

   if (abs(fc2).gt. tol) then
      if (i.eq.1) then
        v0(1)=czero
        v0(2)=czero
        v0(3)=(px-ci*py)/fc
        v0(4)=-fc
      elseif (i.eq.-1) then
        v0(1)=fc
        v0(2)=(px+ci*py)/fc
        v0(3)=czero
        v0(4)=czero
      else
        stop 'v0: i out of range'
      endif
   else
      if (i.eq.1) then
        v0(1)=czero
        v0(2)=czero
        v0(3)=sqrt(cone*two*p0)
        v0(4)=czero
      elseif (i.eq.-1) then
        v0(1)=czero
        v0(2)=sqrt(cone*two*p0)
        v0(3)=czero
        v0(4)=czero
      else
        stop 'v0: i out of range'
      endif
   endif
   return
end function v0

!     v spinor, massive (not used)
function vm(cp,i)
   implicit none
   complex(dp), intent(in) :: cp(4)
   integer, intent(in) :: i
   complex(dp) :: vm(4)
   complex(dp) :: fc2
   real(dp) :: E,x,y,z,m,p,fac

   E=real(cp(1),dp)
   x=real(cp(2),dp)
   y=real(cp(3),dp)
   z=real(cp(4),dp)
   m=dsqrt(dabs(E**2-x**2-y**2-z**2))
   p=dsqrt(dabs(x**2+y**2+z**2))

   fc2 = p + z

   if (abs(fc2).gt. tol) then
      fac = cone/sqrt(2_dp*p*fc2)
      if (i.eq.1) then
         vm(1)=-sqrt(E-p)*fac*(x-ci*y)
         vm(2)=sqrt(E-p)*fac*fc2
         vm(3)=sqrt(E+p)*fac*(x-ci*y)
         vm(4)=-sqrt(E+p)*fac*fc2
      elseif (i.eq.-1) then
         vm(1)=sqrt(E+p)*fac*fc2
         vm(2)=sqrt(E+p)*fac*(x+ci*y)
         vm(3)=-sqrt(E-p)*fac*fc2
         vm(4)=-sqrt(E-p)*fac*(x+ci*y)
      else
        stop 'vm: i out of range'
      endif
   else
      if (i.eq.1) then
        vm(1) = -sqrt(cone*(E-p))
        vm(2) = czero
        vm(3) = sqrt(cone*(E+p))
        vm(4) = czero
      elseif (i.eq.-1) then
        vm(1) = czero
        vm(2) = sqrt(cone*(E+p))
        vm(3) = czero
        vm(4) = -sqrt(cone*(E-p))
      else
        stop 'vm: i out of range'
      endif
   endif
   return
end function vm


function CurrentPrefactor(id,hel,useA)
   implicit none
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   integer :: idV
   logical, optional :: useA
   complex(dp) :: CurrentPrefactor
   integer :: testAcoupl
   integer :: idc(1:2)

   CurrentPrefactor = czero
   if(hel.eq.0) return
   testAcoupl=0
   if(present(useA)) then
      if(useA) testAcoupl=1
   endif
   idV = CoupledVertex(id,hel,useAHcoupl=testAcoupl)
   if(idV.eq.Not_a_particle_) return
   idc(1)=convertLHE(id(1))
   idc(2)=convertLHE(id(2))

   if( idV.eq.Wp_ .or. idV.eq.Wm_ ) then ! W->ffb'
      CurrentPrefactor = ci*eW/sqrt2*CKM(idc(1),idc(2))/ScaleFactor(idc(1),idc(2))
   elseif( ((abs(idc(1)).eq.12) .or. (abs(idc(1)).eq.14) .or. (abs(idc(1)).eq.16)) .and. testAcoupl.eq.0 ) then ! Z->nunu
      CurrentPrefactor = ci*eW_cW*(T3_nulep(1)-Q_nulep(1)*xw)
   elseif( (abs(idc(1)).eq.11) .or. (abs(idc(1)).eq.13) .or. (abs(idc(1)).eq.15) ) then ! A/Z->ll
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_nulep(2)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_nulep(2)-Q_nulep(2)*xw)
      elseif(hel.gt.0) then
         CurrentPrefactor = ci*eW_cW*(-Q_nulep(2)*xw)
      endif
   elseif( (abs(idc(1)).eq.pdfUp_) .or. (abs(idc(1)).eq.pdfChm_) .or. (abs(idc(1)).eq.pdfTop_) ) then ! Z->uu
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_ud(1)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_ud(1)-Q_ud(1)*xw)
      elseif(hel.gt.0) then
         CurrentPrefactor = ci*eW_cW*(-Q_ud(1)*xw)
      endif
   elseif( (abs(idc(1)).eq.pdfDn_) .or. (abs(idc(1)).eq.pdfStr_) .or. (abs(idc(1)).eq.pdfBot_) ) then ! Z->uu
      if(testAcoupl.eq.1) then
         CurrentPrefactor = ci*eA*Q_ud(2)
      elseif(hel.lt.0) then
         CurrentPrefactor = ci*eW_cW*(T3_ud(2)-Q_ud(2)*xw)
      elseif(hel.gt.0) then
         CurrentPrefactor = ci*eW_cW*(-Q_ud(2)*xw)
      endif
   endif
   return
end function CurrentPrefactor

function ScalarFFPrefactor(id)
   implicit none
   integer, intent(in) :: id(1:2)
   complex(dp) :: ScalarFFPrefactor

   ScalarFFPrefactor = czero
   if(id(1).eq.-id(2)) then
      ScalarFFPrefactor = -ci*getMass(id(1))/vev*kappa
   endif
   return
end function ScalarFFPrefactor

function PseudoScalarFFPrefactor(id)
   implicit none
   integer, intent(in) :: id(1:2)
   complex(dp) :: PseudoScalarFFPrefactor

   PseudoScalarFFPrefactor = czero
   if(id(1).eq.-id(2)) then
      PseudoScalarFFPrefactor = getMass(id(1))/vev*kappa_tilde ! = -ci*getMass(id(1))/vev*ci*kappa_tilde
   endif
   return
end function PseudoScalarFFPrefactor

! Returns i*charge*J^mu_L/R
function Vcurrent(p,id,hel, idV,idV_flavorviolating,useAcoupl,Ub_out,V_out,prefactor_out, isOnshellV)
   implicit none
   real(dp), intent(in) :: p(1:4,1:2)
   integer, intent(in) :: id(1:2)
   integer, intent(in) :: hel
   integer, intent(out) :: idV
   integer, optional :: idV_flavorviolating
   logical, optional :: useAcoupl
   complex(dp), optional :: Ub_out(4)
   complex(dp), optional :: V_out(4)
   complex(dp), optional :: prefactor_out
   logical, optional :: isOnshellV
   integer :: idc(1:2), idV_flvio
   complex(dp) :: Vcurrent(4),Ub(4),V(4),prefactor
   integer :: testAcoupl

   ! Initialize return values
   Ub(:)=czero
   V(:)=czero
   Vcurrent(:)=czero
   idV=Not_a_particle_
   idV_flvio=Not_a_particle_
   if( present(isOnshellV) ) then
      isOnshellV=.false.
   endif
   if( present(prefactor_out) ) then
      prefactor_out=czero
   endif

   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1
   endif

   ! Give idV regardless
   idV = CoupledVertex(id,hel,useAHcoupl=testAcoupl)
   idV_flvio = CoupledVertex_FlavorViolating(id,hel,useAHcoupl=testAcoupl)
   if(present(idV_flavorviolating)) then
      idV_flavorviolating = idV_flvio
   endif

   ! Compute ubar and v if there is either an actual association or a hidden W
   ! Hidden W: W->ud, u->W'd' or d->u'W', giving W->uu'W' or dd'W'. The fermions may violate flavor.
   ! Hidden Z, Z->uu, u->Wd', giving Z->ud'W, is also possible. In this case, the fermions would only be mis-identified, but idV would be valid albeit not te corect one.
   ! Hidden Z cases are tested separately in the corresponding amplitudes.
   ! The two hidden vertex cases are the reasons for returning Ub ans V as well.
   if(idV.ne.Not_a_particle_ .or. idV_flvio.ne.Not_a_particle_) then
      idc(1)=convertLHE(id(1))
      idc(2)=convertLHE(id(2))
      if(idc(1).gt.0 .and. idc(2).lt.0) then
         Ub(:)=ubar0(cmplx(p(1:4,1),kind=dp),hel)
         V(:)=v0(cmplx(p(1:4,2),kind=dp),-hel)
      else if(idc(1).lt.0 .and. idc(2).gt.0) then
         Ub(:)=ubar0(cmplx(p(1:4,2),kind=dp),hel)
         V(:)=v0(cmplx(p(1:4,1),kind=dp),-hel)
      endif
   endif
   ! Give ubar and v
   if(present(Ub_out)) then
      Ub_out = Ub
   endif
   if(present(V_out)) then
      V_out = V
   endif

   ! Now skip current computations if the vertex could not be identified for actual association
   ! This means there is either a W in the middle, or one of the particles is an on-shell photon, or there is really no current possible
   if(idV.eq.Not_a_particle_) then
      ! First test if the object is an on-shell photon
      if(                                                               &
         testAcoupl.eq.1 .and.                                          &
         (id(1).eq.Pho_ .or. id(2).eq.Pho_) .and.                       &
         (id(1).eq.Not_a_particle_ .or. id(2).eq.Not_a_particle_)       &
         ) then
         idV=Pho_
         if( present(isOnshellV) ) then
            isOnshellV=.true.
         endif
         if( present(prefactor_out) ) then
            prefactor_out=cone
         endif
      endif
      return
   endif
   ! 1=E,2=px,3=py,4=pz
   ! This is an expression for Ub(+/-) Gamma^\mu V(-/+)
   Vcurrent(1)=(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
   Vcurrent(2)=(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
   Vcurrent(3)=ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
   Vcurrent(4)=(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))

   if(present(useAcoupl)) then
      prefactor=CurrentPrefactor(id,hel,useA=useAcoupl)
   else
      prefactor=CurrentPrefactor(id,hel)
   endif
   if( present(prefactor_out) ) then
      prefactor_out=prefactor
   endif

   Vcurrent(:)=Vcurrent(:)*prefactor
   return
end function Vcurrent

function HFFCurrent(id,Ub,V,prefactor_s_out,prefactor_ps_out)
   integer, intent(in) :: id(1:2)
   complex(dp), intent(in) :: Ub(4)
   complex(dp), intent(in) :: V(4)
   complex(dp) :: HFFCurrent
   complex(dp), optional :: prefactor_s_out
   complex(dp), optional :: prefactor_ps_out
   complex(dp) :: prefactor_s
   complex(dp) :: prefactor_ps

   prefactor_s = ScalarFFPrefactor(id)
   prefactor_ps = PseudoScalarFFPrefactor(id)
   HFFcurrent = ((prefactor_s + prefactor_ps)*(Ub(1)*V(1) + Ub(2)*V(2))) + ((prefactor_s - prefactor_ps)*(Ub(3)*V(3) + Ub(4)*V(4)))

   if( present(prefactor_s_out) ) then
      prefactor_s_out=prefactor_s
   endif
   if( present(prefactor_ps_out) ) then
      prefactor_ps_out=prefactor_ps
   endif

   return
end function

function pol_massive(p,hel)
   implicit none
   real(dp), intent(in) :: p(1:4)
   integer, intent(in) :: hel
   real(dp) :: sincos(1:4), abs3p
   complex(dp) :: inv_mass
   complex(dp) :: pol_massive(1:4)
   real(dp) :: epsilon = 1d-14*GeV ! A small quantity ~ 1e-5 eV

   pol_massive(:)=czero
   sincos = angles(p)

   ! Transverse:
   ! pol_massive(1) ! 0
   ! pol_massive(2) ! 1/sqrt(2) * { cos(theta)*cos(phi) - lambda*i*sin(phi) }
   ! pol_massive(3) ! 1/sqrt(2) * { cos(theta)*sin(phi) + lambda*i*cos(phi) }
   ! pol_massive(4) ! -1/sqrt(2) * sin(theta)
   ! Longitudinal:
   ! pol_massive(1) ! |p|/m
   ! pol_massive(2) ! E/m*px/|p|
   ! pol_massive(3) ! E/m*py/|p|
   ! pol_massive(4) ! E/m*pz/|p|

   if(hel.ne.0) then ! lambda = +-1
      pol_massive(1)= 0d0
      pol_massive(2)= (sincos(3)*sincos(1)-ci*hel*sincos(4))/sqrt2
      pol_massive(3)= (sincos(4)*sincos(1)+ci*hel*sincos(3))/sqrt2
      pol_massive(4)= -sincos(2)/sqrt2
   else ! lambda = L
      abs3p = sqrt(abs(p(2)**2+p(3)**2+p(4)**2))
      inv_mass = (p(1)**2-abs3p**2)
      inv_mass = sqrt(inv_mass) ! make sure E**2-p**2=m**2 is satisfied by extending to the complex plane if needed since sum_h{eps^mu_h eps*^nu_h}=-(g^munu-p^mup^nu/p**2) has to be true even if the invariant mass squared is <0.
      if(abs(inv_mass).gt.epsilon) then
         pol_massive(1)= abs3p/inv_mass
         pol_massive(2)= sincos(3)*sincos(2)*p(1)/inv_mass
         pol_massive(3)= sincos(4)*sincos(2)*p(1)/inv_mass
         pol_massive(4)= sincos(1)*p(1)/inv_mass
      endif
   endif
   return
end function pol_massive

function pol_massless(p,hel)
   implicit none
   real(dp), intent(in) :: p(1:4)
   integer, intent(in) :: hel
   real(dp) :: sincos(1:4)
   complex(dp) :: pol_massless(1:4)

   pol_massless(:)=czero
   sincos = angles(p)

   ! pol_massless(1) ! 0
   ! pol_massless(2) ! 1/sqrt(2) * { cos(theta)*cos(phi) - lambda*i*sin(phi) }
   ! pol_massless(3) ! 1/sqrt(2) * { cos(theta)*sin(phi) + lambda*i*cos(phi) }
   ! pol_massless(4) ! -1/sqrt(2) * sin(theta)

   if(hel.ne.0) then ! lambda = +-1
      pol_massless(1)= 0d0
      pol_massless(2)= (sincos(3)*sincos(1)-ci*hel*sincos(4))/sqrt2
      pol_massless(3)= (sincos(4)*sincos(1)+ci*hel*sincos(3))/sqrt2
      pol_massless(4)= -sincos(2)/sqrt(2d0)
   endif
   return
end function pol_massless

function angles(vector, halfangles)
   implicit none
   real(dp), intent(in) :: vector(4)
   logical, optional :: halfangles
   logical :: usehalftheta
   real(dp) :: angles(4)
   real(dp) :: epsilon = 1d-14*GeV ! A small quantity ~ 1e-5 eV
   real(dp) :: abs3p, abspt

   ! angles(1)=cos(theta)
   ! angles(2)=sin(theta)
   ! angles(3)=cos(phi)
   ! angles(4)=sin(phi)

   usehalftheta=.false.
   if(present(halfangles)) then
      usehalftheta=halfangles
   endif

   angles(:)=0_dp

   abs3p = sqrt(abs(vector(2)**2+vector(3)**2+vector(4)**2))
   if(abs3p.lt.epsilon) then
     angles(1)=1_dp ! theta=0
   else
     angles(1)=vector(4)/abs3p
     angles(2)=sqrt(abs((1d0+angles(1))*(1d0-angles(1))))
     if(usehalftheta) then
        angles(2)=sqrt(abs((1_dp-angles(1))/2_dp))
        angles(1)=sqrt(abs((1_dp+angles(1))/2_dp))
     endif
   endif

   abspt = sqrt(abs(vector(2)**2+vector(3)**2))
   if(abspt.lt.epsilon) then
      angles(3)=1_dp
   else
      angles(3)=vector(2)/abspt
      angles(4)=vector(3)/abspt
   endif
   return
end function angles


! Basic amplitude tools
function dgmunuab_ewk(p) ! 2 g_munu g_ab - g_mua g_nub -g_mub g_nua
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4) ! Lorentz, current index
   complex(dp) :: dgmunuab_ewk

   dgmunuab_ewk = 2_dp*(p(:,1).dot.p(:,2))*(p(:,3).dot.p(:,4))-(p(:,1).dot.p(:,3))*(p(:,2).dot.p(:,4))-(p(:,1).dot.p(:,4))*(p(:,2).dot.p(:,3))
   return
end function dgmunuab_ewk

function dp_gmunu_ewk(p,current) ! (p+-p-)_mu g_ab + (p- -q)_a g_mub + (q-p+)_b g_mua
   use ModMisc
   implicit none
   complex(dp), intent(in) :: current(1:4,1:3) ! Lorentz, current Z/A W+-
   real(dp), intent(in) :: p(1:4,1:2) ! Lorentz, mom. W+-
   real(dp) :: q(1:4)
   complex(dp) :: dpm(1:4),dmq(1:4),dqp(1:4)
   complex(dp) :: dp_gmunu_ewk

   q(:) = -p(:,1)-p(:,2)
   dpm(:) = cmplx(p(:,1)-p(:,2),kind=dp)
   dmq(:) = cmplx(-q(:)+p(:,2),kind=dp)
   dqp(:) = cmplx(q(:)-p(:,1),kind=dp)

   dp_gmunu_ewk = (dpm(:).dot.current(:,1))*(current(:,2).dot.current(:,3)) &
                + (dmq(:).dot.current(:,2))*(current(:,3).dot.current(:,1)) &
                + (dqp(:).dot.current(:,3))*(current(:,1).dot.current(:,2))
   return
end function dp_gmunu_ewk

function epsilon_munuab(e1,e2,e3,e4)
   implicit none
   complex(dp), intent(in) :: e1(4), e2(4), e3(4), e4(4)
   complex(dp) :: epsilon_munuab

   epsilon_munuab =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3)  &
                     -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
                     +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
                     -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
                     +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
                     -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
                     +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
                     -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
                     +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
                     -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
                     +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
                     -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)

   return
end function epsilon_munuab

function epsilon_munuab_single(mu,nu,a,b)
   implicit none
   integer mu,nu,a,b
   real(dp) epsilon_munuab_single
   epsilon_munuab_single=real((mu-nu)*(mu-a)*(mu-b)*(nu-a)*(nu-b)*(a-b),kind=dp)/12_dp
   return
end function epsilon_munuab_single

function Id_Order(Npart,idV,idT)
   implicit none
   integer, intent(in) :: Npart
   integer, intent(in) :: idV(1:NPart),idT(1:NPart)
   integer :: Id_Order(1:Npart),target,part

   Id_Order(:)=0

   do part=1,NPart
      do target=1,NPart
         if(idV(part).eq.idT(target) .and. Id_Order(target).ne.0) then
            Id_Order(target)=part
            exit
         endif
      enddo
   enddo

   return
end function Id_Order



! WWZZ,WWZA,WWAA,WWWW
function QuarticEWKVertex(p,idV,outType) ! p(:,1:4): Currents; order: (incoming W+,incoming W-, outgoing Z or A or W, outgoing Z or A or W); outType(1:2) for Z, A, W+-
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4) ! Lorentz, current index
   integer, intent(in) :: idV(1:4) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:4)
   complex(dp) :: QuarticEWKVertex

   QuarticEWKVertex = czero
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))
   elseif( (outType(1).eq.Pho_ .and. outType(2).eq.Z0_) .or. (outType(2).eq.Pho_ .and. outType(1).eq.Z0_) ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Pho_ /))
   elseif( outType(1).eq.Pho_ .and. outType(2).eq.Pho_ ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Pho_,Pho_ /))
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Z0_,Z0_ /))
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"QuarticEWKVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   if( outType(1).eq.Wp_ .or. outType(1).eq.Wm_ ) then ! W+W-W+W-
      QuarticEWKVertex = ci*dgmunuab_ewk( (/ p(1:4,order(1)), p(1:4,order(3)), p(1:4,order(2)), p(1:4,order(4)) /) )
      QuarticEWKVertex = QuarticEWKVertex*eW**2
   else ! W+W-(ZZ/ZA/AA); symmetric in munu and ab indices
      QuarticEWKVertex = -ci*dgmunuab_ewk( (/ p(1:4,order(1)), p(1:4,order(2)), p(1:4,order(3)), p(1:4,order(4)) /) )
      if(outType(1).eq.Z0_) then
         QuarticEWKVertex = QuarticEWKVertex*eZ
      else
         QuarticEWKVertex = QuarticEWKVertex*eA
      endif
      if(outType(2).eq.Z0_) then
         QuarticEWKVertex = QuarticEWKVertex*eZ
      else
         QuarticEWKVertex = QuarticEWKVertex*eA
      endif
   endif
   return
end function QuarticEWKVertex

! (Z/A)WW
function TripleEWKVertex(p,current,idV,useAcoupl)
   implicit none
   real(dp), intent(in) :: p(1:4,1:3)
   complex(dp), intent(in) :: current(1:4,1:3) ! Lorentz, particle index; p and current are for Z/A and W+W-
   integer, intent(in) :: idV(1:3)
   logical, optional :: useAcoupl
   logical :: forceAcoupl
   integer :: order(1:3)
   complex(dp) :: TripleEWKVertex

   TripleEWKVertex = czero
   forceAcoupl=.false.
   if(present(useAcoupl)) forceAcoupl = useAcoupl
   if(forceAcoupl) then
      order(:)=Id_Order(3,idV,(/ Pho_,Wp_,Wm_ /))
   else
      order(:)=Id_Order(3,idV,(/ Z0_ ,Wp_,Wm_ /))
   endif

   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"TripleEWKVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   TripleEWKVertex = ci*dp_gmunu_ewk( (/ p(1:4,order(2)), p(1:4,order(3)) /), (/ current(1:4,order(1)), current(1:4,order(2)), current(1:4,order(3)) /) )

   if(forceAcoupl) then
      TripleEWKVertex = TripleEWKVertex*eA
   else
      TripleEWKVertex = TripleEWKVertex*eZ
   endif

   return
end function TripleEWKVertex

! (A/Z)ff,f->fA/Zf'f'
function ZAf_fZAHfpfp(diagVCurrent, diagACurrent, diagHCurrent, useAcoupl, useRootPropagator, usebsi)
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2)
   type(VACurrent),intent(in) :: diagACurrent(1:2)
   type(ScalarCurrent),intent(in) :: diagHCurrent(1:2)
   logical, optional :: useAcoupl
   logical, optional :: useRootPropagator
   integer, optional :: usebsi
   type(VACurrent) :: activeCurrent(1:2)
   complex(dp) :: ZAf_fZAHfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: testAcoupl
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: q1scprop(1:2),massPass(1:2,1:2),gFivePass(1:2,1:2),scalarcontribution(1:2),pseudoscalarcontribution(1:2),innerprefactor(1:2),compositeCurrent(1:4,1:2)
   integer :: iV, i1, i2
   logical :: hasRootProp
   integer :: bsicode

   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1 ! For the root of the diagram
   endif
   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif
   if( present(usebsi) ) then
      bsicode=usebsi
   else
      bsicode=kbsi
   endif

   ZAf_fZAHfpfp(:) = czero

   ! Calculate total "q"
   do iV=1,2
      if(diagVCurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
      else if(diagACurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagACurrent(iV)%pVertex(:)
      else
         return
      endif
   enddo


   do i1=1,2
      i2=3-i1

      if(testAcoupl.eq.0) then
         activeCurrent(1)=diagVCurrent(i1)
      else
         activeCurrent(1)=diagACurrent(i1)
      endif
      if(activeCurrent(1)%isOnshellBoson) then
         cycle
      endif
      compositeCurrent(:,:)=czero
      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      gFivePass(:,:)=czero ! gamma5^munu index, fermion/antifermion index
      q1(:,1)=q(:)-activeCurrent(1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+activeCurrent(1)%p(:,1) ! q of the anti-fermion in the middle
      do iV=1,2
         q1scprop(iV) = ScalarPropagator(activeCurrent(1)%id(iV),q1(:,iV),compmass=massPass(3-iV,iV)) ! Propagator of the fermion/anti-fermion in the middle
      enddo

      ! Z/A->f1f1, f1->f1 Z->f2f2
      ! V2=Z
      activeCurrent(2)=diagVCurrent(i2)
      if(activeCurrent(2)%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent(2)%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! Z/A->f1f1, f1->f1 A->f2f2
      ! V2=A
      activeCurrent(2)=diagACurrent(i2)
      if(activeCurrent(2)%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent(2)%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! Z/A->f1f1, f1->f1 H->f2f2
      ! V2=H
      if(diagHCurrent(i2)%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
         scalarcontribution(1) = ScalarFFPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(1) = PseudoScalarFFPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/))*diagHCurrent(i2)%propcurrent
         scalarcontribution(2) = ScalarFFPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(2) = PseudoScalarFFPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/))*diagHCurrent(i2)%propcurrent
         do iV=1,2
            massPass(iV,iV) = scalarcontribution(iV)
            gFivePass(iV,iV) = pseudoscalarcontribution(iV)
         enddo
      endif

      restmp = gamma_gamma_mu(3,(/ compositeCurrent(:,1), cmplx(q1(:,1),kind=dp) /),activeCurrent(1)%Ub,activeCurrent(1)%V, whichUndotted=3, massMatrix=massPass(:,1), gammaFiveMatrix=gFivePass(:,1) )*q1scprop(1) + &
               gamma_gamma_mu(3,(/ cmplx(q1(:,2),kind=dp), compositeCurrent(:,2) /),activeCurrent(1)%Ub,activeCurrent(1)%V, massMatrix=massPass(:,2), gammaFiveMatrix=gFivePass(:,2) )*q1scprop(2)

      if(hasRootProp) then
         restmp = VectorPropagator(activeCurrent(1)%idVertex,q,restmp) ! Add -igmunu/[V] piece
      endif
      restmp(:) = restmp(:)*activeCurrent(1)%prefactor
      ZAf_fZAHfpfp(:) = ZAf_fZAHfpfp(:) + restmp(:)
   enddo
   return
end function ZAf_fZAHfpfp

! (A/Z)ff,f->fWf'f'
! Note: Here, the Z/A received is an observed W, so there is no need to pass diagACurrent since it is not actually observed.
function ZAf_fWfpfp(diagVCurrent, useAcoupl, useRootPropagator)
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2) ! Receives ~W W as input
   logical, optional :: useAcoupl
   logical, optional :: useRootPropagator
   type(VACurrent) :: activeCurrent(1:2)
   complex(dp) :: ZAf_fWfpfp(1:4), restmp(1:4) ! Result is a 'current'
   integer :: testAcoupl
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: q1scprop(1:2),massPass(1:2,1:2),outerprefactor(1:2),innerprefactor,compositeCurrent(1:4,1:2)
   integer :: idVroot, iV, i1, i2
   logical :: hasRootProp

   testAcoupl=0
   if(present(useAcoupl)) then
      if(useAcoupl) testAcoupl=1 ! For the root of the diagram
   endif
   if(testAcoupl.eq.1) then
      idVroot=Pho_
   else
      idVroot=Z0_
   endif
   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif

   ZAf_fWfpfp(:) = czero

   ! Calculate total "q"
   do iV=1,2
      if(diagVCurrent(iV)%idVertex.ne.Wp_ .and. diagVCurrent(iV)%idVertex.ne.Wm_) then
         return
      endif
      q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
   enddo

   do i1=1,2
      i2=3-i1

      activeCurrent(1)=diagVCurrent(i1)
      innerprefactor=activeCurrent(1)%prefactor ! Coupling for the fermion in the middle to the boson branchingo out is just the coupling of the "perceived" root boson

      compositeCurrent(:,:)=czero
      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      q1(:,1)=q(:)-activeCurrent(1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+activeCurrent(1)%p(:,1) ! q of the anti-fermion in the middle
      do iV=1,2
         q1scprop(iV) = ScalarPropagator(-activeCurrent(1)%id(3-iV),q1(:,iV),compmass=massPass(3-iV,iV)) ! Propagator of the fermion/anti-fermion in the middle
      enddo

      ! Z/A->f1f1, f1->f1 W->f2f2
      ! V2=W
      activeCurrent(2)=diagVCurrent(i2)
      if(activeCurrent(2)%idVertex.ne.Not_a_particle_) then
         if(testAcoupl.eq.0) then
            outerprefactor(1) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel)
            outerprefactor(2) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel)
         else
            outerprefactor(1) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel,useA=.true.)
            outerprefactor(2) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel,useA=.true.)
         endif

         do iV=1,2
            compositeCurrent(:,iV) = activeCurrent(2)%propcurrent(:)*outerprefactor(iV)
         enddo
      endif

      restmp = gamma_gamma_mu(3,(/ compositeCurrent(:,1), cmplx(q1(:,1),kind=dp) /),activeCurrent(1)%Ub,activeCurrent(1)%V, whichUndotted=3, massMatrix=massPass(:,1))*q1scprop(1) + &
               gamma_gamma_mu(3,(/ cmplx(q1(:,2),kind=dp), compositeCurrent(:,2) /),activeCurrent(1)%Ub,activeCurrent(1)%V, massMatrix=massPass(:,2))*q1scprop(2)
      restmp(:) = restmp(:)*innerprefactor
      ZAf_fWfpfp(:) = ZAf_fWfpfp(:) + restmp(:)
   enddo
   if(hasRootProp) then
      ZAf_fWfpfp = VectorPropagator(idVroot,q,ZAf_fWfpfp) ! Add -igmunu/[V] piece
   endif
   return
end function ZAf_fWfpfp

! Wff,f->f(A/Z)f'f'
function Wf_fZAHfpfp(diagVCurrent, diagACurrent, diagHCurrent, useRootPropagator, usebsi)
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2)
   type(VACurrent),intent(in) :: diagACurrent(1:2)
   type(ScalarCurrent),intent(in) :: diagHCurrent(1:2)
   logical, optional :: useRootPropagator
   integer, optional :: usebsi
   type(VACurrent) :: activeCurrent(1:2)
   complex(dp) :: Wf_fZAHfpfp(1:4), restmp(1:4) ! Result is a 'current'
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: q1scprop(1:2),massPass(1:2,1:2),gFivePass(1:2,1:2),scalarcontribution(1:2),pseudoscalarcontribution(1:2),innerprefactor(1:2),compositeCurrent(1:4,1:2)
   integer :: iV, i1, i2
   logical :: hasRootProp
   integer :: bsicode

   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif
   if( present(usebsi) ) then
      bsicode=usebsi
   else
      bsicode=kbsi
   endif

   Wf_fZAHfpfp(:) = czero

   ! Calculate total "q"
   do iV=1,2
      if(diagVCurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
      else if(diagACurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagACurrent(iV)%pVertex(:)
      else
         return
      endif
   enddo

   do i1=1,2
      i2=3-i1

      if(diagVCurrent(i1)%idVertex.ne.Wp_ .and. diagVCurrent(i1)%idVertex.ne.Wm_) then
         cycle
      else
         activeCurrent(1)=diagVCurrent(i1)
      endif

      compositeCurrent(:,:)=czero
      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      gFivePass(:,:)=czero ! gamma5^munu index, fermion/antifermion index
      q1(:,1)=q(:)-activeCurrent(1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+activeCurrent(1)%p(:,1) ! q of the anti-fermion in the middle
      do iV=1,2
         q1scprop(iV) = ScalarPropagator(activeCurrent(1)%id(iV),q1(:,iV),compmass=massPass(3-iV,iV)) ! Propagator of the fermion/anti-fermion in the middle
      enddo

      ! W->f1f1, f1->f1 Z->f2f2
      ! V2=Z
      activeCurrent(2)=diagVCurrent(i2)
      if(activeCurrent(2)%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent(2)%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! W->f1f1, f1->f1 A->f2f2
      ! V2=A
      activeCurrent(2)=diagACurrent(i2)
      if(activeCurrent(2)%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/),activeCurrent(1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/),activeCurrent(1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent(2)%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! Z/A->f1f1, f1->f1 H->f2f2
      ! V2=H
      if(diagHCurrent(i2)%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
         scalarcontribution(1) = ScalarFFPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(1) = PseudoScalarFFPrefactor((/ activeCurrent(1)%id(1),-activeCurrent(1)%id(1)/))*diagHCurrent(i2)%propcurrent
         scalarcontribution(2) = ScalarFFPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(2) = PseudoScalarFFPrefactor((/ -activeCurrent(1)%id(2), activeCurrent(1)%id(2)/))*diagHCurrent(i2)%propcurrent
         do iV=1,2
            massPass(iV,iV) = scalarcontribution(iV)
            gFivePass(iV,iV) = pseudoscalarcontribution(iV)
         enddo
      endif

      restmp = gamma_gamma_mu(3,(/ compositeCurrent(:,1), cmplx(q1(:,1),kind=dp) /),activeCurrent(1)%Ub,activeCurrent(1)%V, whichUndotted=3, massMatrix=massPass(:,1), gammaFiveMatrix=gFivePass(:,1) )*q1scprop(1) + &
               gamma_gamma_mu(3,(/ cmplx(q1(:,2),kind=dp), compositeCurrent(:,2) /),activeCurrent(1)%Ub,activeCurrent(1)%V, massMatrix=massPass(:,2), gammaFiveMatrix=gFivePass(:,2) )*q1scprop(2)

      if(hasRootProp) then
         restmp = VectorPropagator(activeCurrent(1)%idVertex,q,restmp) ! Add -igmunu/[V] piece
      endif
      restmp(:) = restmp(:)*activeCurrent(1)%prefactor
      Wf_fZAHfpfp(:) = Wf_fZAHfpfp(:) + restmp(:)
   enddo
   return
end function Wf_fZAHfpfp

! Wff,f->fpWf'f'
! Note: Here, the root boson is observed as a ~Z!
function Wf_fWfpfp(diagVCurrent, useRootPropagator)
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2)
   logical, optional :: useRootPropagator
   type(VACurrent) :: activeCurrent(1:2)
   complex(dp) :: Wf_fWfpfp(1:4), restmp(1:4) ! Result is a 'current'
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: q1scprop,innerprefactor,outerprefactor,allprefactor(1:2),compositeCurrent(1:4,1:2),massPass(1:2,1:2)
   integer :: idRoot(1:2), i1, i2, iV
   logical :: hasRootProp
   integer :: rootPair(1:3),rp

   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif

   Wf_fWfpfp(:) = czero

   ! Calculate total "q"
   do iV=1,2
      q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
   enddo

   do i1=1,2
      i2=3-i1

      compositeCurrent(:,:)=czero
      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      rootPair(:)=Not_a_particle_
      allprefactor(:)=czero
      outerprefactor=czero
      innerprefactor=czero

      activeCurrent(1)=diagVCurrent(i1)
      activeCurrent(2)=diagVCurrent(i2)
      ! activeCurrent(1) has to be an observed ~Z
      ! activeCurrent(2) has to be an observed W
      if(activeCurrent(1)%idVertex.ne.Z0_ .and. activeCurrent(1)%idVertex_flavorviolating.ne.Z0_ .and. activeCurrent(2)%idVertex.ne.Wp_ .and. activeCurrent(2)%idVertex.ne.Wm_) then
         cycle
      endif
      q1(:,1)=q(:)-activeCurrent(1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+activeCurrent(1)%p(:,1) ! q of the anti-fermion in the middle

      ! Fermion case
      rootPair = WDaughterPair(activeCurrent(1)%id(2))
      idRoot(1) = CoupledVertex( (/ rootPair(1), activeCurrent(1)%id(2) /), activeCurrent(1)%hel )
      if(idRoot(1).eq.Wp_ .or. idRoot(1).eq.Wm_) then
         do rp=1,3
            if(rootPair(rp).ne.Not_a_particle_) then
               outerprefactor = CurrentPrefactor((/ rootPair(rp), activeCurrent(1)%id(2) /),activeCurrent(1)%hel) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/ activeCurrent(1)%id(1),-rootPair(rp) /),activeCurrent(1)%hel) ! f1 -> f1 W
               q1scprop = ScalarPropagator(rootPair(rp),q1(:,1))
               allprefactor(1) = allprefactor(1) + outerprefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      ! Anti-fermion case
      rootPair(:)=Not_a_particle_
      rootPair = WDaughterPair(activeCurrent(1)%id(1))
      idRoot(2) = CoupledVertex( (/ activeCurrent(1)%id(1), rootPair(1) /),activeCurrent(1)%hel )
      if(idRoot(2).eq.Wp_ .or. idRoot(2).eq.Wm_) then
         do rp=1,3
            if(rootPair(rp).ne.Not_a_particle_) then
               outerprefactor = CurrentPrefactor((/ activeCurrent(1)%id(1), rootPair(rp) /),activeCurrent(1)%hel) ! W->f1 fb1
               innerprefactor = CurrentPrefactor((/ -rootPair(rp), activeCurrent(1)%id(2) /),activeCurrent(1)%hel) ! fb1 -> fb1 W
               q1scprop = ScalarPropagator(rootPair(rp),q1(:,2))
               allprefactor(2) = allprefactor(2) + outerprefactor*innerprefactor*q1scprop
            endif
         enddo
      endif
      do iV=1,2
         compositeCurrent(:,iV) = activeCurrent(2)%propcurrent(:)*allprefactor(iV)
      enddo

      restmp = gamma_gamma_mu(3,(/ compositeCurrent(:,1), cmplx(q1(:,1),kind=dp) /),activeCurrent(1)%Ub,activeCurrent(1)%V, whichUndotted=3, massMatrix=massPass(:,1)) + &
               gamma_gamma_mu(3,(/ cmplx(q1(:,2),kind=dp), compositeCurrent(:,2) /),activeCurrent(1)%Ub,activeCurrent(1)%V, massMatrix=massPass(:,2))
      Wf_fWfpfp(:) = Wf_fWfpfp(:) + restmp(:)
   enddo
   if(hasRootProp) then
      Wf_fWfpfp = VectorPropagator(Wp_,q,Wf_fWfpfp) ! Add -igmunu/[V] piece
   endif
   return
end function Wf_fWfpfp


! Hff,f->fA/Zf'f'
function Hf_fZAHfpfp(diagVCurrent, diagACurrent, diagHCurrent, useRootPropagator, usebsi) ! Keep usebsi here for future use, but there is no actual use in 4V diagrams since the root itself has a Higgs
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2)
   type(VACurrent),intent(in) :: diagACurrent(1:2)
   type(ScalarCurrent),intent(in) :: diagHCurrent(1:2)
   logical, optional :: useRootPropagator
   integer, optional :: usebsi
   type(VACurrent) :: activeCurrent
   complex(dp) :: Hf_fZAHfpfp, restmp
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: curnull(1:4), q1scprop(1:2),massPass(1:3,1:2),gFivePass(1:3,1:2),scalarcontribution(1:2),pseudoscalarcontribution(1:2),innerprefactor(1:2),compositeCurrent(1:4,1:2)
   integer :: iV, i1, i2
   logical :: hasRootProp
   integer :: bsicode

   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif
   if( present(usebsi) ) then
      bsicode=usebsi
   else
      bsicode=kbsi
   endif

   Hf_fZAHfpfp = czero
   if(.not.do1HVertex) then
      return
   endif
   curnull(:)=czero ! Just a 0-current

   ! Calculate total "q"
   do iV=1,2
      if(diagVCurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
      else if(diagACurrent(iV)%idVertex.ne.Not_a_particle_) then
         q(:)=q(:)+diagACurrent(iV)%pVertex(:)
      else
         return
      endif
   enddo


   do i1=1,2
      i2=3-i1

      if(diagHCurrent(i1)%idVertex.ne.Not_a_particle_) then
         cycle
      endif
      compositeCurrent(:,:)=czero
      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      gFivePass(:,:)=czero ! gamma5^munu index, fermion/antifermion index
      q1(:,1)=q(:)-diagHCurrent(i1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+diagHCurrent(i1)%p(:,1) ! q of the anti-fermion in the middle
      do iV=1,2
         q1scprop(iV) = ScalarPropagator(diagHCurrent(i1)%id(iV),q1(:,iV),compmass=massPass(2,iV)) ! Propagator of the fermion/anti-fermion in the middle
      enddo

      ! Z/A->f1f1, f1->f1 Z->f2f2
      ! V2=Z
      activeCurrent=diagVCurrent(i2)
      if(activeCurrent%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ diagHCurrent(i1)%id(1),-diagHCurrent(i1)%id(1)/),diagHCurrent(i1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -diagHCurrent(i1)%id(2), diagHCurrent(i1)%id(2)/),diagHCurrent(i1)%hel) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! Z/A->f1f1, f1->f1 A->f2f2
      ! V2=A
      activeCurrent=diagACurrent(i2)
      if(activeCurrent%idVertex.ne.Not_a_particle_ .and. bsicode.ne.ksig) then
         innerprefactor(1) = CurrentPrefactor((/ diagHCurrent(i1)%id(1),-diagHCurrent(i1)%id(1)/),diagHCurrent(i1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         innerprefactor(2) = CurrentPrefactor((/ -diagHCurrent(i1)%id(2), diagHCurrent(i1)%id(2)/),diagHCurrent(i1)%hel, useA=.true.) ! Coupling for the fermion in the middle to the boson branchingo out
         do iV=1,2
            compositeCurrent(:,iV) = compositeCurrent(:,iV) + activeCurrent%propcurrent(:)*innerprefactor(iV)
         enddo
      endif
      ! Z/A->f1f1, f1->f1 H->f2f2
      ! V2=H
      if(diagHCurrent(i2)%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
         scalarcontribution(1) = ScalarFFPrefactor((/ diagHCurrent(i1)%id(1),-diagHCurrent(i1)%id(1)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(1) = PseudoScalarFFPrefactor((/ diagHCurrent(i1)%id(1),-diagHCurrent(i1)%id(1)/))*diagHCurrent(i2)%propcurrent
         scalarcontribution(2) = ScalarFFPrefactor((/ -diagHCurrent(i1)%id(2), diagHCurrent(i1)%id(2)/))*diagHCurrent(i2)%propcurrent
         pseudoscalarcontribution(2) = PseudoScalarFFPrefactor((/ -diagHCurrent(i1)%id(2), diagHCurrent(i1)%id(2)/))*diagHCurrent(i2)%propcurrent
         do iV=1,2
            massPass(2*iV-1,iV) = scalarcontribution(iV) ! 1,3
            gFivePass(2*iV-1,iV) = pseudoscalarcontribution(iV) ! 1,3
            massPass(5-2*iV,iV) = diagHCurrent(i1)%scprefactor ! 3,1
            gFivePass(5-2*iV,iV) = diagHCurrent(i1)%pscprefactor ! 3,1
         enddo
      endif


      restmp = gamma_dotted_contracted(3,(/ compositeCurrent(:,1), cmplx(q1(:,1),kind=dp), curnull /),diagHCurrent(i1)%Ub,diagHCurrent(i1)%V, massM=massPass(:,1), gFiveM=gFivePass(:,1) )*q1scprop(1) + &
               gamma_dotted_contracted(3,(/ curnull, cmplx(q1(:,2),kind=dp), compositeCurrent(:,2) /),diagHCurrent(i1)%Ub,diagHCurrent(i1)%V, massM=massPass(:,2), gFiveM=gFivePass(:,2) )*q1scprop(2)

      if(hasRootProp) then
         restmp = restmp*ScalarPropagator(diagHCurrent(i1)%idVertex,q) ! Add i/[H] piece
      endif
      Hf_fZAHfpfp = Hf_fZAHfpfp + restmp
   enddo
   return
end function Hf_fZAHfpfp

! Hff,f->fWf'f'
! Note: Here, the H received is an observed W, so there is no need to pass diagHCurrent since it is not actually observed.
function Hf_fWfpfp(diagVCurrent, useRootPropagator, idHChoice)
   implicit none
   type(VACurrent),intent(in) :: diagVCurrent(1:2) ! Receives ~W W as input
   logical, optional :: useRootPropagator
   integer, optional :: idHChoice
   complex(dp) :: Hf_fWfpfp, restmp
   real(dp) :: q(1:4),q1(1:4,1:2)
   complex(dp) :: curnull(1:4), q1scprop(1:2),massPass(1:3,1:2),gFivePass(1:3,1:2),scalarcontribution(1:2),pseudoscalarcontribution(1:2),innerprefactor
   integer :: idVroot, iV, i1, i2
   logical :: hasRootProp

   idVroot=Hig_
   if( present(useRootPropagator) ) then
      hasRootProp=useRootPropagator
   else
      hasRootProp=.true.
   endif
   if( present(idHChoice) ) then
      idVroot=idHChoice
   else
      idVroot=Hig_
   endif

   Hf_fWfpfp = czero
   if(.not.do1HVertex) then
      return
   endif
   curnull(:)=czero

   ! Calculate total "q"
   do iV=1,2
      if(diagVCurrent(iV)%idVertex.ne.Wp_ .and. diagVCurrent(iV)%idVertex.ne.Wm_) then
         return
      endif
      q(:)=q(:)+diagVCurrent(iV)%pVertex(:)
   enddo

   do i1=1,2
      i2=3-i1

      innerprefactor=diagVCurrent(i1)%prefactor ! Coupling for the fermion in the middle to the boson branchingo out is just the coupling of the "perceived" root boson

      massPass(:,:)=czero ! gamma^munu index, fermion/antifermion index
      q1(:,1)=q(:)-diagVCurrent(i1)%p(:,2) ! q of the fermion in the middle
      q1(:,2)=-q(:)+diagVCurrent(i1)%p(:,1) ! q of the anti-fermion in the middle
      do iV=1,2
         q1scprop(iV) = ScalarPropagator(-diagVCurrent(i1)%id(3-iV),q1(:,iV),compmass=massPass(2,iV)) ! Propagator of the fermion/anti-fermion in the middle
      enddo

      ! H->f1f1, f1->f1 W->f2f2
      ! V2=W
      ! Find out the Higgs couplings at the root
      scalarcontribution(1) = ScalarFFPrefactor((/ -diagVCurrent(i1)%id(2), diagVCurrent(i1)%id(2)/))
      scalarcontribution(2) = ScalarFFPrefactor((/ diagVCurrent(i1)%id(1),-diagVCurrent(i1)%id(1)/))
      pseudoscalarcontribution(1) = PseudoScalarFFPrefactor((/ -diagVCurrent(i1)%id(2), diagVCurrent(i1)%id(2)/))
      pseudoscalarcontribution(2) = PseudoScalarFFPrefactor((/ diagVCurrent(i1)%id(1),-diagVCurrent(i1)%id(1)/))

      do iV=1,2
         ! 1,3 for the second potential Higgs, which DNE in a W-like decay
         massPass(5-2*iV,iV) = czero ! 1,3
         gFivePass(5-2*iV-1,iV) = czero ! 1,3
         ! 3,1 for the Higgs at the root
         massPass(2*iV-1,iV) = scalarcontribution(iV) ! 3,1
         gFivePass(2*iV-1,iV) = pseudoscalarcontribution(iV) ! 3,1
      enddo
      restmp = gamma_dotted_contracted(3,(/ (diagVCurrent(i2)%propcurrent(:)*innerprefactor), cmplx(q1(:,1),kind=dp), curnull /),diagVCurrent(i1)%Ub,diagVCurrent(i1)%V, massM=massPass(:,1), gFiveM=gFivePass(:,1) )*q1scprop(1) + &
               gamma_dotted_contracted(3,(/ curnull, cmplx(q1(:,2),kind=dp), (diagVCurrent(i2)%propcurrent(:)*innerprefactor) /),diagVCurrent(i1)%Ub,diagVCurrent(i1)%V, massM=massPass(:,2), gFiveM=gFivePass(:,2) )*q1scprop(2)
      Hf_fWfpfp = Hf_fWfpfp + restmp
   enddo
   if(hasRootProp) then
      Hf_fWfpfp = Hf_fWfpfp*ScalarPropagator(idVroot,q) ! Add i/[H] piece
   endif

   return
end function Hf_fWfpfp

! (A/Z)ff,f->fA/Zf'f'
function f_NV_fp(diagVCurrent, diagACurrent, diagHCurrent, usebsi)
   implicit none
   integer, parameter :: N=4,maxbit=2**(N-1)-1,maxpaircomb=3**(N-2)-1
   type(VACurrent),intent(in) :: diagVCurrent(1:N)
   type(VACurrent),intent(in) :: diagACurrent(1:N)
   type(ScalarCurrent),intent(in) :: diagHCurrent(1:N)
   integer, optional :: usebsi
   type(VACurrent) :: activeCurrent(1:N)
   complex(dp) :: f_NV_fp, restmp
   integer :: rootPair(1:3,1:N-2),rp(1:N-2)
   integer :: idRoot(2:N), i1, iV, i, order(1:N),code,codeprime
   logical :: jarr(2:N)
   real(dp) :: q(1:4),q1(1:4,1:N-2)
   complex(dp) :: q1scprop(1:N-2)
   complex(dp) :: gammaDottingArray(1:4,1:2*(N-2)+1), gammaDottingArray_inst(1:4,1:2*(N-2)+1)
   complex(dp) :: gammaMassArray(1:2*(N-2)+1),gammaFiveArray(1:2*(N-2)+1) ! No need to have *_inst versions since these change every combination
   logical :: mtc ! These are input to get the group of permutations
   logical :: even ! These are input to get the group of permutations
   logical :: isValid, doSkip, doSkipRoot
   integer :: bsicode

   if(N.lt.3) then
      return
   endif

   if( present(usebsi) ) then
      bsicode=usebsi
   else
      bsicode=kbsi
   endif

   f_NV_fp = czero

   do iV=1,N ! Identity
      order(iV)=iV
   enddo
   mtc = .false. ! Start with this statement
60 continue
   call NEXPER(N, order(:), mtc, even)

   ! Now having the order, evaluate based on base current 1 coupled to 2..N Vs
   activeCurrent(1)=diagVCurrent(order(1))
   if( & ! Need a current as the base
      activeCurrent(1)%idVertex.ne.Not_a_particle_ .or. activeCurrent(1)%idVertex_flavorviolating.ne.Not_a_particle_ &
   ) then
      isValid=.true.

      do iV=2,N ! Loop over the remaining currents
         if(diagVCurrent(order(iV))%idVertex_flavorviolating.ne.Not_a_particle_) then ! The rest have to be simple currents
            isValid=.false.
            exit
         elseif(iV.lt.N) then
            ! Determine the propagators
            if(iV.eq.2) then
               if(diagVCurrent(order(iV))%idVertex.ne.Not_a_particle_) then ! Maybe it is a Z/W
                  q1(:,iV-1) = diagVCurrent(order(iV))%pVertex(:)
               else if(diagACurrent(order(iV))%idVertex.ne.Not_a_particle_) then ! Maybe it is an onshell A
                  q1(:,iV-1) = diagACurrent(order(iV))%pVertex(:)
               else ! If it is neither Z/W nor onshell A, skip this permutation
                  isValid=.false.
                  exit
               endif ! Either way, q1 is the same
               q1(:,iV-1) = q1(:,iV-1) + activeCurrent(1)%p(:,1)
            else
               if(diagVCurrent(order(iV))%idVertex.ne.Not_a_particle_) then ! Maybe it is a Z/W
                  q1(:,iV-1) = diagVCurrent(order(iV))%pVertex(:)
               else if(diagACurrent(order(iV))%idVertex.ne.Not_a_particle_) then ! Maybe it is an onshell A
                  q1(:,iV-1) = diagACurrent(order(iV))%pVertex(:)
               else ! If it is neither Z/W nor onshell A, skip this permutation
                  isValid=.false.
                  exit
               endif ! Either way, q1 is the same
               q1(:,iV-1) = q1(:,iV-1) + q1(:,iV-2)
            endif
            gammaDottingArray(:,2+2*(iV-2))=cmplx(q1(:,iV-1),kind=dp) ! This does not change as long as iV combination is fixed.
         endif
      enddo
      if(isValid) then
         do i=0,maxbit ! =2**(N-1)-1
            doSkip=.false.

            ! Determine possible vertices and currents
            do iV=2,N
               jarr(iV)=btest(i,-iV+2) ! Determine if the bit at location 2-iV is 0 (Z) or 1 (A)

               if(diagVCurrent(order(iV))%idVertex.ne.Not_a_particle_ .and. .not.(jarr(iV))) then
                  activeCurrent(iV)=diagVCurrent(order(iV))
               else if(diagACurrent(order(iV))%idVertex.ne.Not_a_particle_ .and. jarr(iV)) then
                  activeCurrent(iV)=diagACurrent(order(iV))
               else
                  doSkip=.true.
                  exit
               endif

               ! Find the ids of fermion propagators
               idRoot(iV)=activeCurrent(iV)%idVertex
               if(bsicode.ne.ksig) then
                  gammaDottingArray(:,1+2*(iV-2))=activeCurrent(iV)%propcurrent(:) ! This does not change as long as Z/A/W combination is fixed.
               else
                  gammaDottingArray(:,1+2*(iV-2))=czero                            ! This does not change as long as Z/A/W combination is fixed.
               endif
               if(iV.lt.N) then
                  if(idRoot(iV).eq.Wp_ .or. idRoot(iV).eq.Wm_) then
                     if(iV.eq.2) then
                        rootPair(:,iV-1) = WDaughterPair(activeCurrent(iV)%id(1),considerTops=.true.)
                     else
                        rootPair(:,iV-1) = WDaughterPair(-rootPair(1,iV-2),considerTops=.true.) ! No need to loop over the previous rootPair since a quark remains a quark etc.
                     endif
                  else
                     if(iV.eq.2) then
                        rootPair(1,iV-1) = -activeCurrent(iV)%id(1)
                     else
                        rootPair(1,iV-1) = -rootPair(1,iV-2)
                     endif
                     rootPair(2:3,iV-1)=Not_a_particle_
                  endif
               endif
            enddo

            if(doSkip) cycle

            do code=0,maxpaircomb ! =3**(N-2)-1
               codeprime=code
               doSkipRoot=.false.

               gammaDottingArray_inst=gammaDottingArray ! current2, qf1, ... ,qf(N-2),currentN
               gammaMassArray(:)=czero
               gammaFiveArray(:)=czero
               q1scprop(:)=czero
               do iV=1,N-2 ! Choose the indices in the combination
                  rp(iV) = mod(codeprime,3) ! 0, 1 or 2
                  codeprime = (codeprime-rp(iV))/3
                  rp(iV) = rp(iV)+1 ! 1, 2 or 3
                  if(rootPair(rp(iV),iV).eq.Not_a_particle_) then
                     doSkipRoot=.true.
                     exit
                  endif
                  q1scprop(iV) = ScalarPropagator(-rootPair(rp(iV),iV),q1(:,iV),compmass=gammaMassArray(2+2*(iV-1)))
                  gammaDottingArray_inst(:,2+2*(iV-1)) = gammaDottingArray_inst(:,2+2*(iV-1))*q1scprop(iV) ! current2, qf1/Pf1, ... ,qf(N-2)/Pf(N-2),currentN
                  gammaMassArray(2+2*(iV-1)) = gammaMassArray(2+2*(iV-1))*q1scprop(iV) ! 0, mf1/Pf1, ... ,mf(N-2)/Pf(N-2),0
               enddo

               if(doSkipRoot) cycle

               do iV=2,N
                  if(iV.eq.2) then
                     gammaDottingArray_inst(:,1+2*(iV-2)) = gammaDottingArray_inst(:,1+2*(iV-2))*CurrentPrefactor( (/ activeCurrent(1)%id(1), rootPair(rp(iV-1),iV-1) /), activeCurrent(1)%hel, useA=jarr(iV) )
                     if(activeCurrent(iV)%idVertex.eq.Z0_ .and. diagHCurrent(order(iV))%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
                        gammaMassArray(1+2*(iV-2)) = ScalarFFPrefactor( (/ activeCurrent(1)%id(1), rootPair(rp(iV-1),iV-1) /) )*diagHCurrent(order(iV))%propcurrent
                        gammaFiveArray(1+2*(iV-2)) = PseudoScalarFFPrefactor( (/ activeCurrent(1)%id(1), rootPair(rp(iV-1),iV-1) /) )*diagHCurrent(order(iV))%propcurrent
                     endif
                  else if(iV.eq.N) then
                     gammaDottingArray_inst(:,1+2*(iV-2)) = gammaDottingArray_inst(:,1+2*(iV-2))*CurrentPrefactor( (/ -rootPair(rp(iV-2),iV-2), activeCurrent(1)%id(2) /), activeCurrent(1)%hel, useA=jarr(iV) )
                     if(activeCurrent(iV)%idVertex.eq.Z0_ .and. diagHCurrent(order(iV))%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
                        gammaMassArray(1+2*(iV-2)) = ScalarFFPrefactor( (/ -rootPair(rp(iV-2),iV-2), activeCurrent(1)%id(2) /) )*diagHCurrent(order(iV))%propcurrent
                        gammaFiveArray(1+2*(iV-2)) = PseudoScalarFFPrefactor( (/ -rootPair(rp(iV-2),iV-2), activeCurrent(1)%id(2) /) )*diagHCurrent(order(iV))%propcurrent
                     endif
                  else
                     gammaDottingArray_inst(:,1+2*(iV-2)) = gammaDottingArray_inst(:,1+2*(iV-2))*CurrentPrefactor( (/ -rootPair(rp(iV-2),iV-2), rootPair(rp(iV-1),iV-1) /), activeCurrent(1)%hel, useA=jarr(iV) )
                     if(activeCurrent(iV)%idVertex.eq.Z0_ .and. diagHCurrent(order(iV))%idVertex.ne.Not_a_particle_ .and. do1HVertex .and. bsicode.ne.kbkg) then
                        gammaMassArray(1+2*(iV-2)) = ScalarFFPrefactor( (/ -rootPair(rp(iV-2),iV-2), rootPair(rp(iV-1),iV-1) /) )*diagHCurrent(order(iV))%propcurrent
                        gammaFiveArray(1+2*(iV-2)) = PseudoScalarFFPrefactor( (/ -rootPair(rp(iV-2),iV-2), rootPair(rp(iV-1),iV-1) /) )*diagHCurrent(order(iV))%propcurrent
                     endif
                  endif
               enddo

               restmp = gamma_dotted_contracted(2*(N-2)+1, gammaDottingArray_inst, activeCurrent(1)%Ub, activeCurrent(1)%V, massM=gammaMassArray, gFiveM=gammaFiveArray)
               f_NV_fp = f_NV_fp + restmp
            enddo
         enddo
      endif
   endif
   if(mtc) GOTO 60

   return
end function f_NV_fp


! Signal-related amplitudes
! 4H
function QuarticHVertex(p,sccurrent,idV)
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4) ! Lorentz, H index
   complex(dp), intent(in) :: sccurrent(1:4) ! H index; scalar current
   integer, intent(in) :: idV(1:4) ! current ids
   integer :: order(1:4)
   complex(dp) :: QuarticHVertex

   QuarticHVertex = czero
   if(.not.do4HVertex) then
      return
   endif
   order(:)=0
   order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Hig_,Hig_ /))
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"QuarticHVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   QuarticHVertex = (-ci*3_dp*M_Reso**2/vev)*(sccurrent(1)*sccurrent(2)*sccurrent(3)*sccurrent(4))
   return
end function QuarticHVertex

! 3H
function TripleHVertex(p,sccurrent,idV)
   implicit none
   complex(dp), intent(in) :: p(1:4,1:3) ! Lorentz, H index
   complex(dp), intent(in) :: sccurrent(1:3) ! H index; scalar current
   integer, intent(in) :: idV(1:3) ! current ids
   integer :: order(1:3)
   complex(dp) :: TripleHVertex

   TripleHVertex = czero
   if(.not.do3HVertex) then
      return
   endif
   order(:)=0
   order(:)=Id_Order(3,idV,(/ Hig_,Hig_,Hig_ /))
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"TripleHVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   TripleHVertex = (-ci*3_dp*(M_Reso/vev)**2)*(sccurrent(1)*sccurrent(2)*sccurrent(3))
   return
end function TripleHVertex

! VVHH
function DoubleVHVertex(p,sccurrent,current,idV,outType)
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:4),sccurrent(1:4),current(1:4,1:4) ! Lorentz, current index
   integer, intent(in) :: idV(1:4) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:4)
   complex(dp) :: DoubleVHVertex

   DoubleVHVertex = czero
   if(.not.do2HVertex) then
      return
   endif
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Wp_,Wm_ /))
      M_V=M_W
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(4,idV,(/ Hig_,Hig_,Z0_,Z0_ /))
      M_V=M_Z
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 .or. order(4).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"DoubleVHVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   DoubleVHVertex = (ci*2_dp*(M_V/vev)**2)*(sccurrent(order(1))*sccurrent(order(2)))* &
                    ( current(:,order(3)).dot.current(:,order(4)) )
   return
end function DoubleVHVertex

! VVH
function VVHVertex(p,current,idV,outType)
   use ModMisc
   implicit none
   complex(dp), intent(in) :: p(1:4,1:3),current(1:4,1:3) ! Lorentz, current index
   integer, intent(in) :: idV(1:3) ! current ids
   integer, intent(in) :: outType(1:2)
   integer :: order(1:3)
   complex(dp) :: VVHVertex
   complex(dp) :: aa(1:3)
   complex(dp) :: ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
   complex(dp) :: e3_e4,q3_q4,e3_q4,e4_q3
   real(dp) :: q_q,q3_q3,q4_q4

   VVHVertex = czero
   if(.not.do1HVertex) then
      return
   endif
   aa(:)=czero
   if( (outType(1).eq.Wp_ .and. outType(2).eq.Wm_) .or. (outType(2).eq.Wp_ .and. outType(1).eq.Wm_) ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Wp_,Wm_ /))
   elseif( outType(1).eq.Z0_ .and. outType(2).eq.Z0_ ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Z0_,Z0_ /))
   elseif( (outType(1).eq.Z0_ .and. outType(2).eq.Pho_) .or. (outType(1).eq.Pho_ .and. outType(2).eq.Z0_) ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Z0_,Pho_ /))
   elseif( outType(1).eq.Pho_ .and. outType(2).eq.Pho_ ) then
      order(:)=Id_Order(3,idV,(/ Hig_,Pho_,Pho_ /))
   else
      return
   endif
   if( order(1).eq.0 .or. order(2).eq.0 .or. order(3).eq.0 ) then
#if VVHOffshell_debug==1
      print *,"VVHVertex::Target mismatch (ids, order)",idV,order
#endif
      return
   endif

   q_q   = p(:,order(1)).dot.p(:,order(1))
   q3_q3 = p(:,order(2)).dot.p(:,order(2))
   q4_q4 = p(:,order(3)).dot.p(:,order(3))
   q3_q4 = p(:,order(2)).dot.p(:,order(3))
   e3_e4 = current(:,order(2)).dot.current(:,order(3))
   e3_q4 = current(:,order(2)).dot.p(:,order(3))
   e4_q3 = p(:,order(2)).dot.current(:,order(3))

   if( outType(1).ne.Pho_ .and. outType(2).ne.Pho_ ) then
      ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q_q)
      ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q_q)
      if(outType(1).eq.Z0_) then
         M_V=M_Z
      else
         M_V=M_W
      endif
      aa(1) = ghz1_dyn*M_V**2                             &
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   elseif( outType(1).ne.Pho_ .or. outType(2).ne.Pho_ ) then ! Parameterization goes as Zgs in ModHiggs
      ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q_q)
      ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q_q)

      M_V=M_Z
      aa(1) = ghz1_dyn*M_V**2                             &
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   else
      ghz1_dyn = czero
      ghz2_dyn = HVVSpinZeroDynamicCoupling(9, 0d0,0d0,q_q)
      ghz3_dyn = HVVSpinZeroDynamicCoupling(10,0d0,0d0,q_q)
      ghz4_dyn = HVVSpinZeroDynamicCoupling(11,0d0,0d0,q_q)

      aa(1) =                                             &
            + ghz2_dyn*two*q3_q4                          &
            + ghz3_dyn*(q3_q4/Lambda)**2
      aa(2) = -2d0*ghz2_dyn-ghz3_dyn/Lambda**2*q3_q4
      aa(3) = -2d0*ghz4_dyn
   endif

   aa(1) = aa(1)*e3_e4
   aa(2) = aa(2)*e3_q4*e4_q3
   aa(3) = aa(3)*epsilon_munuab(current(:,order(2)),current(:,order(3)),p(:,order(2)),p(:,order(3)))

   VVHVertex = ci*( aa(1) + aa(2) + aa(3) )/vev
   return
end function VVHVertex



! Subroutines to combine the amplitudes

! Signal consists of VV--HH + V--(H)V--VH, VV--V--VH (W+W-ZZH/W+ZWW-H/W-ZWW+H), (VV/HH)--H--(VV/HH), and HHHH

! Signal (VV/HH)--H--(VV/HH), HHHH
subroutine amp_VVHVV(diag)
   implicit none
   class(Single4VDiagram) :: diag
   complex(dp) :: ime

   real(dp) :: qV(1:4,1:4)
   integer :: idH(1:4),idV(1:4),idA(1:4),idSort(1:4)
   complex(dp) :: scalarH(1:4) ! Scalar current
   complex(dp) :: currentV(1:4,1:4) ! Current
   complex(dp) :: currentA(1:4,1:4) ! Current
   complex(dp) :: ime_tmp(1:2),ime_acc
   complex(dp) :: currentH(1:4),qBack(1:4,1:3),qFront(1:4,1:3),cBack(1:4,1:3),cFront(1:4,1:3),hBack(1:3),hFront(1:3),hprop
   integer :: outTypeBack(1:2),outTypeFront(1:2),order(1:4),idBack(1:3),idFront(1:3),ic,lineorder(1:4),iperm,nAonshells,iV,fstype(1:4)
   integer :: nzgscomb,izgscomb,fstype_zgs(1:4)
   real(dp) :: qH(1:4)

   if(.not.doSignal) then
      return
   endif

   ime=czero
   ime_acc=czero
   ime_tmp(:)=czero
   qH(:)=0_dp
   currentH(:)=czero ! This is just a dummy current
   nAonshells=0 ! Keep track of how many gammas there are.

   do ic=1,4
      idH(ic) = fullH(ic)%idVertex
      idV(ic) = fullV(ic)%idVertex
      idA(ic) = fullA(ic)%idVertex

      if(idV(iV).ne.Not_a_particle_) then
         qV(:,iV) = fullV(iV)%pVertex(:)
      elseif(idA(iV).ne.Not_a_particle_) then
         qV(:,iV) = fullA(iV)%pVertex(:)
      else ! Failure!
         return
      endif

      scalarH(ic) = fullH(ic)%propcurrent
      currentV(:,ic) = fullV(ic)%propcurrent(:)
      currentA(:,ic) = fullA(ic)%propcurrent(:)
      if(idA(ic).eq.Pho_ .and. fullA(ic)%isOnshellBoson) then
         nAonshells = nAonshells+1
         if(nAonshells.gt.3) return ! 4A states are not allowed since there have to be at least two fermions generating these diagrams
         idSort(ic)=idA(ic)
         !idH(ic)=Not_a_particle_
      else
         idSort(ic)=idV(ic)
      endif
   enddo

   ! Try (W+W-)-H-(W+W-)
   order(:)=0
   fstype(:)=(/ Wp_,Wm_,Wp_,Wm_ /)
   order(:)=Id_Order(4,idSort,fstype)
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0 .and. do1HVertex) then
   ! As long as Higgs is not charged, an observed (W+W-)x2 can only proceed through a single H in the middle.
       do iperm=1,2
         if(iperm.eq.1) then
            lineorder(:)=(/ 1,2,3,4 /)
         else
            lineorder(:)=(/ 1,4,3,2 /)
         endif
         outTypeBack(:)=(/ fstype(lineorder(1)),fstype(lineorder(2)) /);outTypeFront(:)=(/ fstype(lineorder(3)),fstype(lineorder(4)) /)

         qH(:)=qV(:,order(lineorder(3)))+qV(:,order(lineorder(4)))
         idBack(:)=(/Hig_,idV(order(lineorder(1))),idV(order(lineorder(2)))/)
         qBack(:,1)=qH(:)
         qBack(:,2)=qV(:,order(lineorder(1)))
         qBack(:,3)=qV(:,order(lineorder(2)))
         cBack(:,1)=currentH(:)
         cBack(:,2)=currentV(:,order(lineorder(1)))
         cBack(:,3)=currentV(:,order(lineorder(2)))
         idFront(:)=(/Hig_,idV(order(lineorder(3))),idV(order(lineorder(4)))/)
         qFront(:,1)=-qH(:)
         qFront(:,2)=qV(:,order(lineorder(3)))
         qFront(:,3)=qV(:,order(lineorder(4)))
         cFront(:,1)=currentH(:)
         cFront(:,2)=currentV(:,order(lineorder(3)))
         cFront(:,3)=currentV(:,order(lineorder(4)))
         ime = ime + VVHVertex(qBack,cBack,idBack,outTypeBack)*VVHVertex(qFront,cFront,idFront,outTypeFront)*ScalarPropagator(Hig_,qH(:))
      enddo
   endif

   ! Try (W+W-)-H-(ZZ/HH)
   order(:)=0
   if(nAonshells.eq.0) then
      fstype(:)=(/ Wp_,Wm_,Z0_,Z0_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else if(nAonshells.eq.1) then
      fstype(:)=(/ Wp_,Wm_,Z0_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else if(nAonshells.eq.2) then
      fstype(:)=(/ Wp_,Wm_,Pho_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      lineorder(:)=(/1,2,3,4/)
      outTypeBack(:)=(/ fstype(lineorder(1)),fstype(lineorder(2)) /)

      ime_acc=czero
      idBack(:)=(/Hig_,outTypeBack(1),outTypeBack(2)/)
      qH(:)=qV(:,order(lineorder(3)))+qV(:,order(lineorder(4)))
      qBack(:,1)=qH(:)
      qFront(:,1)=-qH(:)
      hBack(1)=cone; hBack(2:3)=czero
      cBack(:,1)=currentH(:)
      cFront(:,1)=currentH(:)
      do iV=1,2
         qBack(:,iV+1)=qV(:,order(lineorder(iV)))
         qFront(:,iV+1)=qV(:,order(lineorder(iV+2)))
         cBack(:,iV+1)=currentV(:,order(lineorder(iV)))
      enddo

      if(do1HVertex) then
         ime_tmp(1)=VVHVertex(qBack,cBack,idBack,outTypeBack)*ScalarPropagator(Hig_,qH(:))
      else
         ime_tmp(1)=czero
      endif
      ime_tmp(2)=czero

      nzgscomb = 2**(2-nAonshells)
      do izgscomb=0,(nzgscomb-1)
         if(.not.includeGammaStar .and. izgscomb.gt.0) then
            exit
         endif
         fstype_zgs=fstype
         if(nAonshells.lt.2) then
            do iV=3,(4-nAonshells)
               if(btest(izgscomb,-iV+3)) then
                  fstype_zgs(iV)=Pho_
               endif
            enddo
         endif

         outTypeFront(:)=(/ fstype_zgs(lineorder(3)),fstype_zgs(lineorder(4)) /)
         idFront(:)=(/Hig_,outTypeFront(1),outTypeFront(2)/)
         do iV=3,4
            if(outTypeFront(iV-2).eq.Z0_) then
               cFront(:,iV-1)=currentV(:,order(lineorder(iV)))
            else
               cFront(:,iV-1)=currentA(:,order(lineorder(iV)))
            endif
         enddo

         if(do1HVertex) then
            ime_tmp(2) = ime_tmp(2) + VVHVertex(qFront,cFront,idFront,outTypeFront) ! H->V3V4
         endif
         if( & ! H->H3H4->ffb f'fb'
            outTypeFront(1).eq.Z0_ .and. outTypeFront(2).eq.Z0_ .and. &
            idH(order(lineorder(3))).ne.Not_a_particle_ .and. idH(order(lineorder(4))).ne.Not_a_particle_ &
            ) then

            hFront(1)=cone
            do iV=3,4
               hFront(iV-1)=scalarH(order(lineorder(iV)))
            enddo

            if(do3HVertex) then
               ime_tmp(2) = ime_tmp(2) + TripleHVertex(qFront, hFront, (/ Hig_,Hig_,Hig_ /) )
            endif
         endif
      enddo
      ime = ime + ime_tmp(1)*ime_tmp(2)
      ime_tmp(:)=czero
   endif

   ! Try (ZZ)-H-(ZZ)
   order(:)=0
   ! It is important that Zs come first! The izgscomb algorithm depends on it.
   if(nAonshells.eq.0) then
      fstype(:)=(/ Z0_,Z0_,Z0_,Z0_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else if(nAonshells.eq.1) then
      fstype(:)=(/ Z0_,Z0_,Z0_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else if(nAonshells.eq.2) then
      fstype(:)=(/ Z0_,Z0_,Pho_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else
      fstype(:)=(/ Z0_,Pho_,Pho_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      do iperm=1,3
         if(iperm.eq.1) then
            lineorder(:)=(/1,2,3,4/)
         else if(iperm.eq.2) then
            lineorder(:)=(/1,3,2,4/)
         else
            lineorder(:)=(/1,4,2,3/)
         endif

         ime_acc=czero
         qH(:)=qV(:,order(lineorder(3)))+qV(:,order(lineorder(4)))
         hprop = ScalarPropagator(Hig_,qH(:))
         qBack(:,1)=qH(:)
         qFront(:,1)=-qH(:)
         do iV=1,2
            qBack(:,iV+1)=qV(:,order(lineorder(iV)))
            qFront(:,iV+1)=qV(:,order(lineorder(iV+2)))
         enddo
         cBack(:,1)=currentH(:)
         cFront(:,1)=currentH(:)

         nzgscomb = 2**(4-nAonshells)
         do izgscomb=0,(nzgscomb-1)
            if(.not.includeGammaStar .and. izgscomb.gt.0) then
               exit
            endif
            fstype_zgs=fstype
            do iV=1,(4-nAonshells)
               if(btest(izgscomb,-iV+1)) then
                  fstype_zgs(iV)=Pho_
               endif
            enddo
            outTypeBack(:)=(/ fstype_zgs(lineorder(1)),fstype_zgs(lineorder(2)) /)
            outTypeFront(:)=(/ fstype_zgs(lineorder(3)),fstype_zgs(lineorder(4)) /)

            idBack(:)=(/Hig_,outTypeBack(1),outTypeBack(2)/)
            idFront(:)=(/Hig_,outTypeFront(1),outTypeFront(2)/)
            do iV=1,2
               if(outTypeBack(iV).eq.Z0_) then
                  cBack(:,iV+1)=currentV(:,order(lineorder(iV)))
               else
                  cBack(:,iV+1)=currentA(:,order(lineorder(iV)))
               endif
            enddo
            do iV=3,4
               if(outTypeFront(iV-2).eq.Z0_) then
                  cFront(:,iV-1)=currentV(:,order(lineorder(iV)))
               else
                  cFront(:,iV-1)=currentA(:,order(lineorder(iV)))
               endif
            enddo

            ime_tmp(:)=czero
            if(do1HVertex) then
               ime_tmp(1) = ime_tmp(1) + VVHVertex(qBack,cBack,idBack,outTypeBack)
               ime_tmp(2) = ime_tmp(2) + VVHVertex(qFront,cFront,idFront,outTypeFront)
            endif
            if( & ! H->H1H2->ffb f'fb'
               outTypeBack(1).eq.Z0_ .and. outTypeBack(2).eq.Z0_ .and. &
               idH(order(lineorder(1))).ne.Not_a_particle_ .and. idH(order(lineorder(2))).ne.Not_a_particle_ &
               ) then
               hBack(1)=cone
               do iV=1,2
                  hBack(iV+1)=scalarH(order(lineorder(iV)))
               enddo
               if(do3HVertex) then
                  ime_tmp(1) = ime_tmp(1) + TripleHVertex(qBack, hBack, (/ Hig_,Hig_,Hig_ /) )
               endif
            endif
            if( & ! H->H3H4->ffb f'fb'
               outTypeFront(1).eq.Z0_ .and. outTypeFront(2).eq.Z0_ .and. &
               idH(order(lineorder(3))).ne.Not_a_particle_ .and. idH(order(lineorder(4))).ne.Not_a_particle_ &
               ) then
               hFront(1)=cone
               do iV=3,4
                  hFront(iV-1)=scalarH(order(lineorder(iV)))
               enddo
               if(do3HVertex) then
                  ime_tmp(2) = ime_tmp(2) + TripleHVertex(qFront, hFront, (/ Hig_,Hig_,Hig_ /) )
               endif
            endif
            ime_acc = ime_acc + ime_tmp(1)*ime_tmp(2)*hprop

            if( & ! H1H2H3H4
               iperm.eq.1 .and. & ! The other combinations are also the same, so do not pursue them!
               outTypeBack(1).eq.Z0_ .and. outTypeBack(2).eq.Z0_ .and. &
               outTypeFront(1).eq.Z0_ .and. outTypeFront(2).eq.Z0_ .and. &
               idH(order(lineorder(1))).ne.Not_a_particle_ .and. idH(order(lineorder(2))).ne.Not_a_particle_ .and. &
               idH(order(lineorder(3))).ne.Not_a_particle_ .and. idH(order(lineorder(4))).ne.Not_a_particle_       &
               ) then
               if(do4HVertex) then
                  ime_acc = ime_acc + QuarticHVertex(cmplx(qV(:,:),kind=dp),scalarH,idH) ! Does not matter if you order them by order(lineorder(:)) or not
               endif
            endif
         enddo
         ime = ime + ime_acc
      enddo
      ime_acc=czero
   endif

   diag%Amplitude = diag%Amplitude + diag%permutation_factor*ime
   return
end subroutine amp_VVHVV

! Quartic, triple-triple, and triple-HVV vertices EWK vertex diagrams, signal or background
subroutine amp_WWZA(diag)
   implicit none
   class(Single4VDiagram) :: diag
   integer :: idH(1:4),idV(1:4),idA(1:4),idSort(1:4),idSorted(1:4)
   real(dp) :: qV(1:4,1:4),qWV(1:4,1:2),qSorted(1:4,1:4)
   complex(dp) :: currentV(1:4,1:4),currentA(1:4,1:4),currentH(1:4),scalarH(1:4),currentSorted(1:4,1:4),scalarHsorted(1:4),cWV(1:4,1:2),propSc(1:2)
   integer :: iV,ic,iz,h,mu,nAonshells,order(1:4),fstype(1:4),fstype_zgs(1:4),izgscomb,nzgscomb
   complex(dp) :: ime,ime_tmp(1:2)

   ime=czero
   nAonshells=0
   currentH(:)=czero ! This is just a dummy current
   do iV=1,4
      idH(iV) = fullH(iV)%idVertex
      idV(iV) = fullV(iV)%idVertex
      idA(iV) = fullA(iV)%idVertex
      scalarH(iV) = fullH(iV)%propcurrent
      currentV(:,iV)=fullV(iV)%propcurrent(:)
      currentA(:,iV)=fullA(iV)%propcurrent(:)

      if(idV(iV).ne.Not_a_particle_) then
         qV(:,iV) = fullV(iV)%pVertex(:)
      elseif(idA(iV).ne.Not_a_particle_) then
         qV(:,iV) = fullA(iV)%pVertex(:)
      else ! Failure!
         return
      endif

      idSort(iV)=idV(iV)
      if(idA(iV).eq.Pho_ .and. fullA(iV)%isOnshellBoson) then
         nAonshells = nAonshells+1
         if(nAonshells.gt.2) return
         idSort(iV)=idA(iV)
      endif
   enddo

   order(:)=0
   if(nAonshells.eq.0) then
      fstype(:)=(/ Wp_,Wm_,Z0_,Z0_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else if(nAonshells.eq.1) then
      fstype(:)=(/ Wp_,Wm_,Z0_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   else
      fstype(:)=(/ Wp_,Wm_,Pho_,Pho_ /)
      order(:)=Id_Order(4,idSort,fstype)
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      nzgscomb = 2**(2-nAonshells)
      do izgscomb=0,(nzgscomb-1)
         if(.not.includeGammaStar .and. izgscomb.gt.0) then
            exit
         endif
         fstype_zgs=fstype
         if(nAonshells.lt.2) then
            do iV=3,(4-nAonshells)
               if(btest(izgscomb,-iV+3)) then
                  fstype_zgs(iV)=Pho_
               endif
            enddo
         endif
         do iV=1,4
            qSorted(:,iV)=qV(:,order(iV))
            scalarHsorted(iV)=czero
            if(fstype_zgs(iV).eq.Pho_) then
               currentSorted(:,iV)=currentA(:,order(iV))
               idSorted(iV)=idA(order(iV))
            else
               currentSorted(:,iV)=currentV(:,order(iV))
               idSorted(iV)=idV(order(iV))
               if(fstype_zgs(iV).eq.Z0_) then
                  scalarHsorted(iV)=scalarH(order(iV))
               endif
            endif
         enddo

         if(doBkg) then ! WW.VV'
            ime = ime + QuarticEWKVertex(currentSorted, idSorted, fstype_zgs(3:4))
         endif
         do ic=3,4
            iz=7-ic
            qWV(:,1) = -qSorted(:,1)-qSorted(:,ic)
            qWV(:,2) = -qSorted(:,2)-qSorted(:,iz)
            propSc(1) = ScalarPropagator(-idSorted(1),qWV(:,1))
            do h=-1,1
               ime_tmp(:)=czero
               cWV(:,1) = pol_massive(qWV(:,1),h)
               do mu=1,4
                  cWV(mu,2) = -conjg(cWV(mu,1))
               enddo

               if(doBkg) then ! WV--W--WV'
                  ime_tmp(1) = ime_tmp(1) + TripleEWKVertex( (/ qWV(:,1), qSorted(:,1), qSorted(:,ic) /), (/ cWV(:,1), currentSorted(:,1), currentSorted(:,ic) /), (/ -idSorted(1), idSorted(1), idSorted(ic) /), useAcoupl=(idSorted(ic).eq.Pho_) )
                  ime_tmp(2) = ime_tmp(2) + TripleEWKVertex( (/ qWV(:,2), qSorted(:,2), qSorted(:,iz) /), (/ cWV(:,2), currentSorted(:,2), currentSorted(:,iz) /), (/ -idSorted(2), idSorted(2), idSorted(iz) /), useAcoupl=(idSorted(iz).eq.Pho_) )
               endif
               if(do1HVertex) then ! WH--W--WH
                  if(idSorted(ic).eq.Z0_) then
                     ime_tmp(1) = ime_tmp(1) + VVHVertex( (/ cmplx(qWV(:,1),kind=dp), cmplx(qSorted(:,1),kind=dp), cmplx(qSorted(:,ic),kind=dp) /), (/ cWV(:,1), currentSorted(:,1), currentH(:) /), (/ -idSorted(1), idSorted(1), idH(order(ic)) /), (/ -fstype_zgs(1),fstype_zgs(1) /) )
                  endif
                  if(idSorted(iz).eq.Z0_) then
                     ime_tmp(2) = ime_tmp(2) + VVHVertex( (/ cmplx(qWV(:,2),kind=dp), cmplx(qSorted(:,2),kind=dp), cmplx(qSorted(:,iz),kind=dp) /), (/ cWV(:,2), currentSorted(:,2), currentH(:) /), (/ -idSorted(2), idSorted(2), idH(order(iz)) /), (/ -fstype_zgs(2),fstype_zgs(2) /) )
                  endif
               endif
               ime = ime + propSc(1)*ime_tmp(1)*ime_tmp(2)
            enddo
         enddo

         if(do1HVertex .and. .not.(idSorted(3).eq.Pho_ .and. idSorted(4).eq.Pho_)) then ! WW--Z/G--HZ/G; H has to derive from a Z, to avoid double-counting at the least.
            do ic=3,4
               iz=7-ic
               qWV(:,1) = -qSorted(:,1)-qSorted(:,2)
               qWV(:,2) = -qSorted(:,ic)-qSorted(:,iz)
               propSc(1)=ScalarPropagator(Z0_,qWV(:,1))
               propSc(2)=ScalarPropagator(Pho_,qWV(:,1))

               do h=-1,1
                  cWV(:,1) = pol_massive(qWV(:,1),h)
                  do mu=1,4
                     cWV(mu,2) = -conjg(cWV(mu,1))
                  enddo

                  ! Z in the middle
                  ime_tmp(:)=czero
                  ime_tmp(1) = ime_tmp(1) + TripleEWKVertex( (/ qWV(:,1), qSorted(:,1), qSorted(:,2) /), (/ cWV(:,1), currentSorted(:,1), currentSorted(:,2) /), (/ Z0_, idSorted(1), idSorted(2) /) )
                  if(idSorted(ic).eq.Z0_) then
                     ime_tmp(2) = ime_tmp(2) + VVHVertex( (/ cmplx(qWV(:,2),kind=dp), cmplx(qSorted(:,iz),kind=dp), cmplx(qSorted(:,ic),kind=dp) /), (/ cWV(:,2), currentSorted(:,iz), currentH(:) /), (/ Z0_, idSorted(iz), idH(order(ic)) /), (/ Z0_,fstype_zgs(iz) /) )
                  endif
                  if(idSorted(iz).eq.Z0_) then
                     ime_tmp(2) = ime_tmp(2) + VVHVertex( (/ cmplx(qWV(:,2),kind=dp), cmplx(qSorted(:,ic),kind=dp), cmplx(qSorted(:,iz),kind=dp) /), (/ cWV(:,2), currentSorted(:,ic), currentH(:) /), (/ Z0_, idSorted(ic), idH(order(iz)) /), (/ Z0_,fstype_zgs(ic) /) )
                  endif
                  ime = ime + propSc(1)*ime_tmp(1)*ime_tmp(2)
                  ! A in the middle
                  ime_tmp(:)=czero
                  ime_tmp(1) = ime_tmp(1) + TripleEWKVertex( (/ qWV(:,1), qSorted(:,1), qSorted(:,2) /), (/ cWV(:,1), currentSorted(:,1), currentSorted(:,2) /), (/ Pho_, idSorted(1), idSorted(2) /), useAcoupl=.true. )
                  if(idSorted(ic).eq.Z0_) then
                     ime_tmp(2) = ime_tmp(2) + VVHVertex( (/ cmplx(qWV(:,2),kind=dp), cmplx(qSorted(:,iz),kind=dp), cmplx(qSorted(:,ic),kind=dp) /), (/ cWV(:,2), currentSorted(:,iz), currentH(:) /), (/ Pho_, idSorted(iz), idH(order(ic)) /), (/ Pho_,fstype_zgs(iz) /) )
                  endif
                  if(idSorted(iz).eq.Z0_) then
                     ime_tmp(2) = ime_tmp(2) + VVHVertex( (/ cmplx(qWV(:,2),kind=dp), cmplx(qSorted(:,ic),kind=dp), cmplx(qSorted(:,iz),kind=dp) /), (/ cWV(:,2), currentSorted(:,ic), currentH(:) /), (/ Pho_, idSorted(ic), idH(order(iz)) /), (/ Pho_,fstype_zgs(ic) /) )
                  endif
                  ime = ime + propSc(2)*ime_tmp(1)*ime_tmp(2)
               enddo
            enddo
         endif
         if( & ! H->H3H4->ffb f'fb'
            do2HVertex .and. &
            fstype_zgs(3).eq.Z0_ .and. fstype_zgs(4).eq.Z0_ .and. &
            idH(order(3)).ne.Not_a_particle_ .and. idH(order(4)).ne.Not_a_particle_ &
            ) then
            ime = ime + DoubleVHVertex( &
                                       cmplx(qSorted(:,:),kind=dp), & ! W+W-HH
                                       scalarHsorted, &
                                       currentSorted, &
                                       (/ idSorted(1),idSorted(2),Hig_,Hig_ /), &
                                       fstype_zgs(1:2) )
         endif

      enddo
   endif

   diag%Amplitude = diag%Amplitude + diag%permutation_factor*ime
   return
end subroutine amp_WWZA

subroutine amp_WWWW(diag)
   implicit none
   class(Single4VDiagram) :: diag
   integer :: idV(1:4)
   real(dp) :: qV(1:4,1:4),qZA(1:4,1:2)
   complex(dp) :: currentV(1:4,1:4),cZA(1:4,1:2) ! Current
   complex(dp) :: ime,ime_WW_Z_WW,ime_WW_A_WW
   integer :: outType(1:2),order(1:4),h,mu,ic

   if(.not.doBkg) then
      return
   endif

   qZA(:,:) = 0_dp
   cZA(:,:) = czero
   ime = czero
   ime_WW_Z_WW = czero
   ime_WW_A_WW = czero

   do ic=1,4
      idV(ic) = fullV(ic)%idVertex
      currentV(:,ic) = fullV(ic)%propcurrent(:)
      if(idV(ic).ne.Not_a_particle_) then
         qV(:,ic) = fullV(ic)%pVertex(:)
      else ! Failure!
         return
      endif
   enddo

   outType(:) = (/ Wp_, Wm_ /)
   order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      ime = QuarticEWKVertex(currentV,idV,outType)

      qZA(:,1) = -qV(:,order(1)) - qV(:,order(2))
      qZA(:,2) = -qV(:,order(3)) - qV(:,order(4))
      do h=-1,1
         cZA(:,1) = pol_massive(qZA(:,1),h)
         do mu=1,4
            cZA(mu,2) = -conjg(cZA(mu,1)) ! cZA(:,2) = -pol_massive(qZA(:,2),h)
         enddo

         ime_WW_Z_WW = ime_WW_Z_WW + &
                       TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(2)) /), (/ cZA(:,1), currentV(:,order(1)), currentV(:,order(2)) /), (/ Z0_, idV(order(1)), idV(order(2)) /) ) * &
                       TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(4)) /), (/ cZA(:,2), currentV(:,order(3)), currentV(:,order(4)) /), (/ Z0_, idV(order(3)), idV(order(4)) /) )

         if(h.ne.0) then
            ! These lines have the same effect
            !cZA(:,1) = pol_massless(qZA(:,1),h)
            !do mu=1,4
            !   cZA(mu,2) = -conjg(cZA(mu,1))
            !enddo
            ime_WW_A_WW = ime_WW_A_WW + &
                          TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(2)) /), (/ cZA(:,1), currentV(:,order(1)), currentV(:,order(2)) /), (/ Pho_, idV(order(1)), idV(order(2)) /), useAcoupl=.true. ) * &
                          TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(4)) /), (/ cZA(:,2), currentV(:,order(3)), currentV(:,order(4)) /), (/ Pho_, idV(order(3)), idV(order(4)) /), useAcoupl=.true. )
         endif
      enddo
      ime_WW_Z_WW = ime_WW_Z_WW * ScalarPropagator(Z0_,qZA(:,1))
      ime_WW_A_WW = ime_WW_A_WW * ScalarPropagator(Pho_,qZA(:,1))

      qZA(:,1) = -qV(:,order(1)) - qV(:,order(4))
      qZA(:,2) = -qV(:,order(3)) - qV(:,order(2))
      do h=-1,1
         cZA(:,1) = pol_massive(qZA(:,1),h)
         do mu=1,4
            cZA(mu,2) = -conjg(cZA(mu,1))
         enddo

         ime_WW_Z_WW = ime_WW_Z_WW + &
                       TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(4)) /), (/ cZA(:,1), currentV(:,order(1)), currentV(:,order(4)) /), (/ Z0_, idV(order(1)), idV(order(4)) /) ) * &
                       TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(2)) /), (/ cZA(:,2), currentV(:,order(3)), currentV(:,order(2)) /), (/ Z0_, idV(order(3)), idV(order(2)) /) )

         if(h.ne.0) then
            cZA(:,1) = pol_massless(qZA(:,1),h)
            do mu=1,4
               cZA(mu,2) = -conjg(cZA(mu,1))
            enddo
            ime_WW_A_WW = ime_WW_A_WW + &
                          TripleEWKVertex( (/ qZA(:,1), qV(:,order(1)), qV(:,order(4)) /), (/ cZA(:,1), currentV(:,order(1)), currentV(:,order(4)) /), (/ Pho_, idV(order(1)), idV(order(4)) /), useAcoupl=.true. ) * &
                          TripleEWKVertex( (/ qZA(:,2), qV(:,order(3)), qV(:,order(2)) /), (/ cZA(:,2), currentV(:,order(3)), currentV(:,order(2)) /), (/ Pho_, idV(order(3)), idV(order(2)) /), useAcoupl=.true. )
         endif
      enddo
      ime_WW_Z_WW = ime_WW_Z_WW * ScalarPropagator(Z0_,qZA(:,1))
      ime_WW_A_WW = ime_WW_A_WW * ScalarPropagator(Pho_,qZA(:,1))
      ime = ime + ime_WW_Z_WW + ime_WW_A_WW

      diag%Amplitude = diag%Amplitude + diag%permutation_factor*ime
   endif
   return
end subroutine amp_WWWW

! Triple EWK vertex and HVV diagrams with fermions in t channel (V->4f), signal or background
subroutine amp_VV_VtVp(diag)
   implicit none
   class(Single4VDiagram) :: diag
   integer :: idV(1:4),idA(1:4),idSort(1:4),idSort_loose(1:4)
   real(dp) :: qtV(1:4),qV(1:4,1:4)
   complex(dp) :: currentV(1:4,1:4),currentA(1:4,1:4)
   integer :: ic,iz,nAonshells,order(1:4),lineorder(1:4),iperm,nzgscomb,izgscomb,idzgs(1:2),nLoose
   complex(dp) :: composite(1:4,1:2),nullcur(1:4),curzgs(1:4,1:3),scalarcomposite,ime
   integer :: ibsi

   ime=czero
   nullcur(:)=czero
   nAonshells=0
   do ic=1,4
      currentV(:,ic)=fullV(ic)%propcurrent(:)
      currentA(:,ic)=fullA(ic)%propcurrent(:)

      idV(ic) = fullV(ic)%idVertex
      if(idV(ic).eq.Not_a_particle_) then
         idV(ic) = fullV(ic)%idVertex_flavorviolating
      endif
      idA(ic) = fullA(ic)%idVertex
      if(idA(ic).eq.Not_a_particle_) then
         idA(ic) = fullA(ic)%idVertex_flavorviolating
      endif

      if(idV(ic).ne.Not_a_particle_) then
         qV(:,ic) = fullV(ic)%pVertex(:)
      elseif(idA(ic).ne.Not_a_particle_) then
         qV(:,ic) = fullA(ic)%pVertex(:)
      endif

      idSort(ic)=idV(ic)
      if(idA(ic).eq.Pho_ .and. fullA(ic)%isOnshellBoson) then
         nAonshells = nAonshells+1
         if(nAonshells.gt.3) return ! The case of 3 on-shell As is for the diagram GG--H->ff->ffA
         idSort(ic)=idA(ic)
      endif
      idSort_loose(ic)=idSort(ic)
      if(fullV(ic)%idVertex_flavorviolating.ne.Not_a_particle_) then
         nLoose = nLoose+1
         idSort_loose(ic)=fullV(ic)%idVertex_flavorviolating
      endif
   enddo

   order(:)=0
   if(nAonshells.eq.0) then
      order(:)=Id_Order(4,idSort,(/ Z0_,Z0_,Z0_,Z0_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
   else if(nAonshells.eq.1) then
      order(:)=Id_Order(4,idSort,(/ Z0_,Z0_,Z0_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
   else if(nAonshells.eq.2) then
      order(:)=Id_Order(4,idSort,(/ Z0_,Z0_,Pho_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
   else
      order(:)=Id_Order(4,idSort,(/ Z0_,Pho_,Pho_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0 .and. do1HVertex) then ! HVV or Hff always enters here
      do iperm=1,3
         if(iperm.eq.1) then
            lineorder(:)=(/1,2,3,4/) ! 1+2--H--3+4 / 3+4--H--1+2
         else if(iperm.eq.2) then
            lineorder(:)=(/1,3,2,4/) ! 1+3--H--2+4 / 2+4--H--1+3
         else
            lineorder(:)=(/1,4,2,3/) ! 1+4--H--2+3 / 2+3--H--1+4
         endif

         do ic=1,3,2
            iz=ic+1

            if(fullH(order(lineorder(4-ic)))%idVertex.eq.Not_a_particle_ .and. fullH(order(lineorder(6-iz)))%idVertex.eq.Not_a_particle_) then
               cycle
            endif

            qtV(:) = -qV(:,order(lineorder(ic)))-qV(:,order(lineorder(iz)))
            scalarcomposite = Hf_fZAHfpfp( & ! No need for bsi index since what follows always uses a vertex with a Higgs involved
                                          (/ fullV(order(lineorder(4-ic))),fullV(order(lineorder(6-iz))) /), &
                                          (/ fullA(order(lineorder(4-ic))),fullA(order(lineorder(6-iz))) /), &
                                          (/ fullH(order(lineorder(4-ic))),fullH(order(lineorder(6-iz))) /) &
                                         ) ! H at the main vertex
            ! Do (Z/A)(Z/A)H->ff(H/Z/A)
            nzgscomb=2
            if(idSort(order(lineorder(ic))).eq.Pho_) then
               nzgscomb=nzgscomb-1
            endif
            if(idSort(order(lineorder(iz))).eq.Pho_) then
               nzgscomb=nzgscomb-1
            endif
            nzgscomb=2**nzgscomb
            do izgscomb=0,nzgscomb-1
               if(izgscomb.gt.0 .and. .not.includeGammastar) then
                  exit
               endif
               curzgs(:,1)=nullcur(:)
               if(btest(izgscomb,0)) then
                  idzgs(1)=Pho_
               else
                  idzgs(1)=idSort(order(lineorder(ic))) ! Pass idSort on purpose so that id=(Z/A,A) are also supported.
               endif
               if(btest(izgscomb,1)) then
                  idzgs(2)=Pho_
               else
                  idzgs(2)=idSort(order(lineorder(iz))) ! Pass idSort on purpose so that id=(Z/A,A) are also supported.
               endif
               if(idzgs(1).eq.Pho_) then
                  curzgs(:,2)=currentA(:,order(lineorder(ic)))
               else
                  curzgs(:,2)=currentV(:,order(lineorder(ic)))
               endif
               if(idzgs(2).eq.Pho_) then
                  curzgs(:,3)=currentA(:,order(lineorder(iz)))
               else
                  curzgs(:,3)=currentV(:,order(lineorder(iz)))
               endif
               ime = ime + scalarcomposite * VVHVertex( &
                                                         (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(lineorder(ic))),kind=dp), cmplx(qV(:,order(lineorder(iz))),kind=dp) /), &
                                                         curzgs, &
                                                         (/ Hig_, idzgs(1:2) /), &
                                                         idzgs &
                                                      )
            enddo
            ! HHH->ff(Z/A/H)
            if(do3HVertex .and. idSort(order(lineorder(ic))).eq.Z0_ .and. idSort(order(lineorder(iz))).eq.Z0_) then
               ime = ime + TripleHVertex( &
                                          (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(lineorder(ic))),kind=dp), cmplx(qV(:,order(lineorder(iz))),kind=dp) /), &
                                          (/ scalarcomposite, fullH(order(lineorder(ic)))%propcurrent, fullH(order(lineorder(iz)))%propcurrent /), &
                                          (/ Hig_, fullH(order(lineorder(ic)))%idVertex, fullH(order(lineorder(iz)))%idVertex /) &
                                        )
            endif

         enddo
      enddo
   endif

   ! From this point, nAonshell>1 has no applicability.
   if(nAonshells.gt.1) then
      return
   endif
   order(:)=0
   if(nAonshells.eq.0) then
      order(:)=Id_Order(4,idSort,(/ Wp_,Wm_,Z0_,Z0_ /)) ! Z0_ is dummy here, just to ensure identical particles are present.
   else
      order(:)=Id_Order(4,idSort,(/ Wp_,Wm_,Z0_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical particles are present.
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0 .and. nLoose.lt.(2-nAonshells)) then
      ! WWZ/A/H.2fZ2f
      qtV(:) = -qV(:,order(1))-qV(:,order(2))

      ! W(ic)Z(iz)W.ffpWffp + W(ic)H(iz)W.ffpWffp
      do ic=1,2
         do iz=3,4
            if(nAonshells.eq.1 .and. iz.eq.3) then ! Avoid passing on-shell photon into the Wf_fWfpfp
               cycle
            endif
            qtV(:) = -qV(:,order(ic))-qV(:,order(iz))
            composite(:,1) = Wf_fWfpfp( (/ fullV(order(7-iz)),fullV(order(3-ic)) /) )
            if(doBkg) then
               ime = ime + TripleEWKVertex( (/ qtV(:), qV(:,order(ic)), qV(:,order(iz)) /), (/ composite(:,1), currentV(:,order(ic)), currentV(:,order(iz)) /), (/ idSort(order(3-ic)), idSort(order(ic)), idSort(order(iz)) /) )
            endif
            if(do1HVertex) then
               ime = ime + VVHVertex( &
                                       (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(ic)),kind=dp), cmplx(qV(:,order(iz)),kind=dp) /), &
                                       (/ composite(:,1), currentV(:,order(ic)), currentV(:,order(iz)) /), &
                                       (/ idSort(order(3-ic)), idSort(order(ic)), fullH(order(iz))%idVertex /), &
                                       (/ idSort(order(3-ic)), idSort(order(ic)) /) &
                                    )
            endif
         enddo
      enddo

      if(nLoose.eq.0) then
         if(.not.doBkg) then
            ibsi=ksig ! ZAf_fZAHfpfp has its separate protection against .not.do1HVertex, but this setting would be the only one to protect against a bkg. diagram.
         else if(.not.do1HVertex) then
            ibsi=kbkg ! Not necessary, but anyway...
         else
            ibsi=kbsi
         endif
         composite(:,1) = ZAf_fZAHfpfp( &
                                          (/ fullV(order(3)),fullV(order(4)) /), &
                                          (/ fullA(order(3)),fullA(order(4)) /), &
                                          (/ fullH(order(3)),fullH(order(4)) /), &
                                          usebsi=ibsi                            &
                                       ) ! Z at the main vertex
         composite(:,2) = ZAf_fZAHfpfp( &
                                          (/ fullV(order(3)),fullV(order(4)) /), &
                                          (/ fullA(order(3)),fullA(order(4)) /), &
                                          (/ fullH(order(3)),fullH(order(4)) /), &
                                          useAcoupl=.true.,                      &
                                          usebsi=ibsi                            &
                                       ) ! A at the main vertex
         ime = ime + &
               TripleEWKVertex( (/ qtV(:), qV(:,order(1)), qV(:,order(2)) /), (/ composite(:,1), currentV(:,order(1)), currentV(:,order(2)) /), (/ Z0_ , idSort(order(1)), idSort(order(2)) /) ) &
             + TripleEWKVertex( (/ qtV(:), qV(:,order(1)), qV(:,order(2)) /), (/ composite(:,2), currentV(:,order(1)), currentV(:,order(2)) /), (/ Pho_, idSort(order(1)), idSort(order(2)) /) )

         if(do1HVertex) then
            ! WWH.ffp(Z/A/H)ffp
            scalarcomposite = Hf_fZAHfpfp( &
                                       (/ fullV(order(3)),fullV(order(4)) /), &
                                       (/ fullA(order(3)),fullA(order(4)) /), &
                                       (/ fullH(order(3)),fullH(order(4)) /)  &
                                       ) ! H at the main vertex
            ime = ime + scalarcomposite*VVHVertex( &
                                                   (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(1)),kind=dp), cmplx(qV(:,order(2)),kind=dp) /), &
                                                   (/ nullcur(:), currentV(:,order(1)), currentV(:,order(2)) /), &
                                                   (/ Hig_, idSort(order(1)), idSort(order(2)) /), &
                                                   (/ idSort(order(1)), idSort(order(2)) /) &
                                                 )
            ! (Z/A)(Z/A)H.ffpWffp
            scalarcomposite = Hf_fWfpfp( (/ fullV(order(1)),fullV(order(2)) /) ) ! H at the main vertex
            nzgscomb=2**(2-nAonshells)
            do izgscomb=0,nzgscomb-1
               if(izgscomb.gt.0 .and. .not.includeGammastar) then
                  exit
               endif
               curzgs(:,1)=nullcur(:)
               if(btest(izgscomb,0)) then
                  idzgs(1)=Pho_
               else
                  idzgs(1)=idSort(order(3)) ! Pass idSort on purpose so that id=(Z/A,A) are also supported.
               endif
               if(btest(izgscomb,1)) then
                  idzgs(2)=Pho_
               else
                  idzgs(2)=idSort(order(4)) ! Pass idSort on purpose so that id=(Z/A,A) are also supported.
               endif
               if(idzgs(1).eq.Pho_) then
                  curzgs(:,2)=currentA(:,order(3))
               else
                  curzgs(:,2)=currentV(:,order(3))
               endif
               if(idzgs(2).eq.Pho_) then
                  curzgs(:,3)=currentA(:,order(4))
               else
                  curzgs(:,3)=currentV(:,order(4))
               endif
               ime = ime + scalarcomposite * VVHVertex( &
                                                         (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(3)),kind=dp), cmplx(qV(:,order(4)),kind=dp) /), &
                                                         curzgs, &
                                                         (/ Hig_, idzgs(1:2) /), &
                                                         idzgs &
                                                      )
            enddo
            ! HHH.ffpWffp
            if(do3HVertex .and. idSort(order(3)).eq.Z0_ .and. idSort(order(4)).eq.Z0_) then
               ime = ime + TripleHVertex( &
                                          (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(3)),kind=dp), cmplx(qV(:,order(4)),kind=dp) /), &
                                          (/ scalarcomposite, fullH(order(3))%propcurrent, fullH(order(4))%propcurrent /), &
                                          (/ Hig_, fullH(order(3))%idVertex, fullH(order(4))%idVertex /) &
                                        )
            endif

         endif
      endif
   endif

   order(:)=0
   if(nAonshells.eq.0) then
      order(:)=Id_Order(4,idSort,(/ Wp_,Wm_,Wp_,Wm_ /)) ! Some of W+ or W- here could actually be Zs observed as Ws!
   endif
   if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
      ! W(ic)W(iz)(Z/A).ffp(~W(4-ic))W(6-iz)ffp + W(ic)W(iz)H.ffp(~W(4-ic))W(6-iz)ffp
      do ic=1,3,2
         do iz=2,4,2
            qtV(:) = -qV(:,order(ic))-qV(:,order(iz))
            if(doBkg) then
               composite(:,1) = ZAf_fWfpfp( (/ fullV(order(4-ic)),fullV(order(6-iz)) /) ) ! Z at the main vertex
               composite(:,2) = ZAf_fWfpfp( (/ fullV(order(4-ic)),fullV(order(6-iz)) /), useAcoupl=.true. ) ! A at the main vertex
               ime = ime &
                         + TripleEWKVertex( (/ qtV(:), qV(:,order(ic)), qV(:,order(iz)) /), (/ composite(:,1), currentV(:,order(ic)), currentV(:,order(iz)) /), (/ Z0_ , idSort(order(ic)), idSort(order(iz)) /) ) &
                         + TripleEWKVertex( (/ qtV(:), qV(:,order(ic)), qV(:,order(iz)) /), (/ composite(:,2), currentV(:,order(ic)), currentV(:,order(iz)) /), (/ Pho_, idSort(order(ic)), idSort(order(iz)) /) )
            endif
            if(do1HVertex) then
               scalarcomposite = Hf_fWfpfp( (/ fullV(order(4-ic)),fullV(order(6-iz)) /) ) ! H at the main vertex
               ime = ime + scalarcomposite*VVHVertex( &
                                                      (/ cmplx(qtV(:),kind=dp), cmplx(qV(:,order(ic)),kind=dp), cmplx(qV(:,order(iz)),kind=dp) /), &
                                                      (/ nullcur(:), currentV(:,order(ic)), currentV(:,order(iz)) /), &
                                                      (/ Hig_, idSort(order(ic)), idSort(order(iz)) /), &
                                                      (/ idSort(order(ic)), idSort(order(iz)) /) &
                                                    )
            endif
         enddo
      enddo
   endif

   diag%Amplitude = diag%Amplitude + diag%permutation_factor*ime
   return
end subroutine amp_VV_VtVp

! Zs, As or Ws in the t-channel
subroutine amp_tVVV(diag)
   use ModMisc
   implicit none
   class(Single4VDiagram) :: diag
   integer :: idV(1:4),idA(1:4),idSort(1:4),idSort_loose(1:4)
   real(dp) :: qV(1:4,1:4),qtchannel(1:4)
   integer :: ic,iperm,nLoose,nAonshells,order(1:4),lineorder(1:4)
   complex(dp) :: composite(1:4,1:8),scalarcomposite(1:2),ime
   integer :: ibsi

   ime=czero
   nAonshells=0
   composite(:,:)=czero
   scalarcomposite(:)=czero
   do ic=1,4
      idV(ic) = fullV(ic)%idVertex
      if(idV(ic).eq.Not_a_particle_) then
         idV(ic) = fullV(ic)%idVertex_flavorviolating
      endif
      idA(ic) = fullA(ic)%idVertex
      if(idA(ic).eq.Not_a_particle_) then
         idA(ic) = fullA(ic)%idVertex_flavorviolating
      endif

      if(idV(ic).ne.Not_a_particle_) then
         qV(:,ic) = fullV(ic)%pVertex(:)
      elseif(idA(ic).ne.Not_a_particle_) then
         qV(:,ic) = fullA(ic)%pVertex(:)
      else ! Failure!
         return
      endif

      if(idA(ic).eq.Pho_ .and. fullA(ic)%isOnshellBoson) then
         nAonshells = nAonshells+1
         if(nAonshells.gt.3) return ! 4A states are not allowed since there have to be at least two fermions generating these diagrams
         idSort(ic)=idA(ic)
      else
         idSort(ic)=idV(ic)
      endif
      idSort_loose(ic)=idSort(ic)
      if(fullV(ic)%idVertex_flavorviolating.ne.Not_a_particle_) then
         nLoose = nLoose+1
         idSort_loose(ic)=fullV(ic)%idVertex_flavorviolating
      endif
   enddo

   if(nAonshells.eq.3) then
   ! Case 0: Z.A+A+A
      order(:)=Id_Order(4,idSort,(/ Z0_,Pho_,Pho_,Pho_ /)) ! Z0_ ensures identical fermions are present.
      if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0 .and. nLoose.eq.0) then ! nLoose=0 is to ensure that flavor violation does not ever happen at the base current
         ! Add the 3V.base current diagrams
         if(.not.doBkg) then
            ibsi=ksig
         else if(.not.do1HVertex) then
            ibsi=kbkg ! Not necessary, but anyway...
         else
            ibsi=kbsi
         endif
         ime = ime + f_NV_fp( &
                              (/ fullV(order(1)),fullV(order(2)),fullV(order(3)),fullV(order(4)) /), &
                              (/ fullA(order(1)),fullA(order(2)),fullA(order(3)),fullA(order(4)) /), &
                              (/ fullH(order(1)),fullH(order(2)),fullH(order(3)),fullH(order(4)) /), &
                              usebsi=ibsi                                                            &
                            )
      endif
   ! End Case 0 Z.A+A+A
   else ! From here on, disallow 3A states
   ! Case 1: 12.Z.34.Z.56.Z.78
      if(nAonshells.eq.0) then
         order(:)=Id_Order(4,idSort,(/ Z0_,Z0_,Z0_,Z0_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
      else if(nAonshells.eq.1) then
         order(:)=Id_Order(4,idSort,(/ Z0_,Z0_,Z0_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
      else
         order(:)=Id_Order(4,idSort,(/ Z0_,Pho_,Z0_,Pho_ /)) ! Z0_ is dummy here, just to ensure identical fermions are present.
      endif
      if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
         do iperm=1,3
            if(iperm.eq.1) then
               lineorder(:)=(/1,2,3,4/) ! 1+2--Z--3+4
            else if(iperm.eq.2) then
               lineorder(:)=(/1,3,2,4/) ! 1+3--Z--2+4
               if(nAonshells.eq.2) then
                  cycle
               endif
            else
               lineorder(:)=(/1,4,2,3/) ! 1+4--Z--2+3
            endif

            qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
            composite(:,:)=czero
            scalarcomposite(:)=czero
            do ic=1,2
               ! NOTE: composite(:,1:4) is always calculated as BSI if do1HVertex, and as Bkg if !do1HVertex.
               !       composite(:,5:8) is always calculated as Bkg.
               !       Therefore, (1.2+3.4)-(5.6+7.8) corresponds to (do1HVertex.Bkg)+(Bkg.do1HVertex)+(do1HVertex.do1HVertex).
               composite(:,ic)   = ZAf_fZAHfpfp( &
                                                (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                useRootPropagator=(ic.eq.1)                                          &
                                               )                                                ! Z at the main vertex
               composite(:,ic+2) = ZAf_fZAHfpfp( &
                                                (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                useAcoupl=.true., useRootPropagator=(ic.eq.1)                        &
                                               )                                                ! A at the main vertex
               if(.not.doBkg) then ! Subtract the bkg.bkg contributions
                  composite(:,ic+4)   = ZAf_fZAHfpfp( &
                                                   (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                   (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                   (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                   useRootPropagator=(ic.eq.1),                                         &
                                                   usebsi=kbkg                                                          &
                                                  )                                                ! Z at the main vertex
                  composite(:,ic+6) = ZAf_fZAHfpfp( &
                                                   (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                   (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                   (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                   useAcoupl=.true., useRootPropagator=(ic.eq.1),                       &
                                                   usebsi=kbkg                                                          &
                                                  )                                                ! A at the main vertex
               endif
               if(do1HVertex) then
                  scalarcomposite(ic) = Hf_fZAHfpfp( &
                                                   (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                   (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                   (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                   useRootPropagator=(ic.eq.1)                                          &
                                                   )                                               ! H at the main vertex
               endif
            enddo
            ime = ime + (composite(:,1).dot.composite(:,2)) + (composite(:,3).dot.composite(:,4)) + scalarcomposite(1)*scalarcomposite(2)
            if(.not.doBkg) then ! Subtract the bkg.bkg contributions if !doBkg
               ime = ime - (composite(:,5).dot.composite(:,6)) - (composite(:,7).dot.composite(:,8))
            endif
         enddo

      ! Add the 3V.base current diagrams
      if(.not.doBkg) then
         ibsi=ksig
      else if(.not.do1HVertex) then
         ibsi=kbkg ! Not necessary, but anyway...
      else
         ibsi=kbsi
      endif
      ime = ime + f_NV_fp( &
                           (/ fullV(order(1)),fullV(order(2)),fullV(order(3)),fullV(order(4)) /), &
                           (/ fullA(order(1)),fullA(order(2)),fullA(order(3)),fullA(order(4)) /), &
                           (/ fullH(order(1)),fullH(order(2)),fullH(order(3)),fullH(order(4)) /), &
                           usebsi=ibsi                                                            &
                         )
      endif
      ! End Case 1 obs(ZZZZ)

      ! Case 2: 12.Z.34.W.56.Z.78 / 12.W.34.W(~Z).56.W.78 / 12.W.34.Z.56.Z.78 == obs(WWZZ): W(Z) boson in the second(third) diagram looks like a ~Z(~W) due to a W branching out.
      if(nAonshells.eq.0) then
         order(:)=Id_Order(4,idSort_loose,(/ Wp_,Wm_,Z0_,Z0_ /))
      else if(nAonshells.eq.1) then
         order(:)=Id_Order(4,idSort_loose,(/ Wp_,Wm_,Z0_,Pho_ /))
      else
         order(:)=Id_Order(4,idSort_loose,(/ Wp_,Wm_,Pho_,Pho_ /))
      endif
      if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
         do iperm=1,3
            if(iperm.eq.1) then
               lineorder(:)=(/1,3,2,4/) ! W1Z3----W2Z4
            else if(iperm.eq.2) then
               lineorder(:)=(/1,4,2,3/) ! W1Z4----W2Z3
            else
               lineorder(:)=(/1,2,3,4/) ! W1W2----Z3Z4
               if(nAonshells.eq.2) then
                  cycle
               endif
            endif

            qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
            composite(:,:)=czero
            scalarcomposite(:)=czero
            do ic=1,2
               if(iperm.lt.3) then
                  ! NOTE: composite(:,1:4) is always calculated as BSI if do1HVertex, and as Bkg if !do1HVertex.
                  !       composite(:,5:8) is always calculated as Bkg.
                  !       Therefore, (1.2+3.4)-(5.6+7.8) corresponds to (do1HVertex.Bkg)+(Bkg.do1HVertex)+(do1HVertex.do1HVertex).
                  composite(:,ic)  = Wf_fZAHfpfp( &
                                                 (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                 (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                 (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                 useRootPropagator=(ic.eq.1)                                          &
                                                ) ! W at the branch
                  if( .not.(fullV(order(lineorder(2*ic)))%isOnshellBoson .or. fullV(order(lineorder(2*ic-1)))%isOnshellBoson) ) then
                     composite(:,ic+2) = Wf_fWfpfp( (/ fullV(order(lineorder(2*ic))),fullV(order(lineorder(2*ic-1))) /), useRootPropagator=(ic.eq.1) ) ! Swaps lineorder to pass ~Z first, W at the branch
                  endif
                  composite(:,ic)=composite(:,ic)+composite(:,ic+2) ! Can add bc. the W at the root is the same.
                  if(.not.doBkg) then ! Subtract the bkg.bkg contributions
                     composite(:,ic+4)  = Wf_fZAHfpfp( &
                                                       (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                       (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                       (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                       useRootPropagator=(ic.eq.1),                                         &
                                                       usebsi=kbkg                                                          &
                                                     ) ! W at the branch
                     composite(:,ic+4)=composite(:,ic+4)+composite(:,ic+2) ! Can add bc. the W at the root is the same.
                  endif
                  composite(:,ic+2)=czero
               else
                  if(ic.eq.1) then
                     composite(:,ic)   = ZAf_fWfpfp ( (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /) )                   ! Z at the main vertex
                     composite(:,ic+2) = ZAf_fWfpfp ( (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), useAcoupl=.true. ) ! A at the main vertex
                     if(do1HVertex) then
                        scalarcomposite(ic)   = Hf_fWfpfp ( (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /) )             ! H at the main vertex
                     endif
                  else
                     if(.not.doBkg) then
                        ibsi=ksig
                     else if(.not.do1HVertex) then
                        ibsi=kbkg ! Not necessary, but anyway...
                     else
                        ibsi=kbsi
                     endif
                     composite(:,ic)   = ZAf_fZAHfpfp( &
                                                      (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                      (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                      (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                      useRootPropagator=.false.,                                           &
                                                      usebsi=ibsi                                                          &
                                                      )                   ! Z at the main vertex
                     composite(:,ic+2) = ZAf_fZAHfpfp( &
                                                      (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                      (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                      (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                      useAcoupl=.true., useRootPropagator=.false.,                         &
                                                      usebsi=ibsi                                                          &
                                                      ) ! A at the main vertex
                     if(do1HVertex) then
                        scalarcomposite(ic)   = Hf_fZAHfpfp( &
                                                         (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                         (/ fullA(order(lineorder(2*ic-1))),fullA(order(lineorder(2*ic))) /), &
                                                         (/ fullH(order(lineorder(2*ic-1))),fullH(order(lineorder(2*ic))) /), &
                                                         useRootPropagator=.false. )                ! H at the main vertex
                     endif
                  endif
               endif
            enddo
            if(iperm.lt.3) then
               ime = ime + (composite(:,1).dot.composite(:,2))
               if(.not.doBkg) then ! Subtract the bkg.bkg contributions
                  ime = ime - (composite(:,5).dot.composite(:,6))
               endif
            else
               ime = ime + (composite(:,1).dot.composite(:,2)) + (composite(:,3).dot.composite(:,4)) + scalarcomposite(1)*scalarcomposite(2)
            endif
         enddo

         ! Add the 3V.base current diagrams
         if(.not.doBkg) then
            ibsi=ksig
         else if(.not.do1HVertex) then
            ibsi=kbkg ! Not necessary, but anyway...
         else
            ibsi=kbsi
         endif
         ime = ime + f_NV_fp( &
                              (/ fullV(order(1)),fullV(order(2)),fullV(order(3)),fullV(order(4)) /), &
                              (/ fullA(order(1)),fullA(order(2)),fullA(order(3)),fullA(order(4)) /), &
                              (/ fullH(order(1)),fullH(order(2)),fullH(order(3)),fullH(order(4)) /), &
                              usebsi=ibsi                                                            &
                            )
      endif
   ! End Case 2 obs(WWZZ)

   ! Case 3: 12.W+.34.Z.56.W-.78 == obs(WWWW)
      order(:)=Id_Order(4,idV,(/ Wp_,Wm_,Wp_,Wm_ /))
      if(order(1).ne.0 .and. order(2).ne.0 .and. order(3).ne.0 .and. order(4).ne.0) then
         do iperm=1,2
            if(iperm.eq.1) then
               lineorder(:)=(/1,2,3,4/) ! 1+2--Z--3+4
            else
               lineorder(:)=(/1,4,3,2/) ! 1+4--Z--3+2
            endif

            qtchannel(:) = -qV(:,order(lineorder(1)))-qV(:,lineorder(2))
            composite(:,:)=czero
            scalarcomposite(:)=czero
            do ic=1,2
               if(doBkg) then
                  composite(:,ic)   = ZAf_fWfpfp( (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), useRootPropagator=(ic.eq.1) )                   ! Z at the main vertex
                  composite(:,ic+2) = ZAf_fWfpfp( (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), useAcoupl=.true., useRootPropagator=(ic.eq.1) ) ! A at the main vertex
               endif
               if(do1HVertex) then
                  scalarcomposite(ic) = Hf_fWfpfp( &
                                                   (/ fullV(order(lineorder(2*ic-1))),fullV(order(lineorder(2*ic))) /), &
                                                   useRootPropagator=(ic.eq.1) )                ! H at the main vertex
               endif
            enddo
            ime = ime + (composite(:,1).dot.composite(:,2)) + (composite(:,3).dot.composite(:,4)) + scalarcomposite(1)*scalarcomposite(2)
         enddo

         ! Add the 3V.base current diagrams
         if(.not.doBkg) then
            ibsi=ksig
         else if(.not.do1HVertex) then
            ibsi=kbkg ! Not necessary, but anyway...
         else
            ibsi=kbsi
         endif
         ime = ime + f_NV_fp( &
                              (/ fullV(order(1)),fullV(order(2)),fullV(order(3)),fullV(order(4)) /), &
                              (/ fullA(order(1)),fullA(order(2)),fullA(order(3)),fullA(order(4)) /), &
                              (/ fullH(order(1)),fullH(order(2)),fullH(order(3)),fullH(order(4)) /), &
                              usebsi=ibsi                                                            &
                            )
      endif
   ! End Case 3 obs(WWWW)
   endif

   diag%Amplitude = diag%Amplitude + diag%permutation_factor*ime
   return
end subroutine amp_tVVV


end module ModVVHOffshell

