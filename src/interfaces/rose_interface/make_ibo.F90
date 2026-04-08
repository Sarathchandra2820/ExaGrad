module make_ibo

   ! references : 
   ! [0] : Pipek and Mezey, J. Chem. Phys. 90, 4916 (1989)
   ! [1] : Knizia, J. Chem. Theory. Comput. 9, 4834 (2013)
   ! [2] : Dubillard et al. J. Chem. Phys. 124, 154307 (2006)
   ! [3] : G. Berghold et al., Phys Rev B. 61 (2000) 10040–10048. doi:10.1103/PhysRevB.61.10040.

   use rose_global

contains

   subroutine Localization(CMO_IAO,LMO_IAO,nmo_frag,Expnt)
   ! make intrinsic bond orbitals (IBOs) (localized occupied orbitals expressed in terms of IAOs)
   ! It consists in effectively minimizing the number of atoms an orbital is centered on.
   ! For this, unitary operations (rotations) are perform amongst the occupied (or valence virtual) MOs
   ! (leaving the Slater Determinant invariant) until a function (that contains information
   ! about the # of electrons located on the IAO of each atoms) is maximised.

        use linear_algebra
        use properties

        implicit none

        integer, intent(in)        :: nmo_frag(:)    ! Size of the fragments
        real(8), intent(in)        :: CMO_IAO(:,:,:) ! Canonical MOs in terms of IAOs
        real(8), intent(inout)     :: LMO_IAO(:,:,:) ! Localized MO in terms of IAOs
        integer, intent(in)        :: Expnt
        logical                    :: conv
        integer                    :: i,j,k,ij, imo, nmoi, it
        integer                    :: n_fragments, nact, nmo2, rcq
        real(8)                    :: phi, func, fGrad, Aij, absBij, Qii, Qjj, max_condition, reQij
        real(8)                    :: damping
        real(8), allocatable       :: Bij(:,:)
        real(8), allocatable       :: Bij_previous(:)
        real(8), allocatable       :: Ai(:,:), Aj(:,:), Aj0(:,:)
        real(8), allocatable       :: Pij(:),Qij(:),res(:)

        nmo2 = size(CMO_IAO,1)
        nact = size(CMO_IAO,2)
        rcq  = size(CMO_IAO,3)
        n_fragments = size(nmo_frag)

        allocate(Pij(rcq),Qij(rcq),res(rcq))
        allocate(Bij(nact*(nact-1)/2,rcq))
        if (rcq .eq. 1) then
           allocate(Bij_previous(nact*(nact-1)/2))
           Bij_previous = 0.D0
           damping = 1.D0
        endif

        LMO_IAO = CMO_IAO

        write(luout,*)
        write(luout,*)"----- Localization procedure: 2x2 Rotations -----"
        write(luout,'(A,I1,A)')" (Exponent = ",Expnt,")"
        write(luout,*)

        conv = .false.

        write(luout,'(1A5,A16,A12)') "Iter","Function","Gradient"
        func = 0.D0
        do it=1,500
          ! Compute the function in Eq.(4) of [1]
          ! Note that it is the same as in Eq.(3) of [2] 
          ! because the IAOs are orthonormalized
          ! The goal is to maximize this function. 
          call get_localization_func(LMO_IAO,nmo_frag,Expnt,func)
  
          ! Start the orbital localization procedure:
          ! Sum over all occupied orbital pairs i,j ; j < i
          fGrad = 0.D0
          ij = 0
          do i=2,nact
             do j=1,i-1
                ij = ij + 1
                ! Aij and Bij are given in Appendix D of [1]
                ! Bij is the gradient (see g_ij in Eq.(7) of [2])
                Bij(ij,:) = 0.D0
                imo = 0
                ! Increment Bij for each fragment:
                do k=1,n_fragments
                   nmoi = nmo_frag(k)
                   allocate(Ai(nmoi,rcq))
                   allocate(Aj(nmoi,rcq))
                   Ai = LMO_IAO(imo+1:imo+nmoi,i,:)
                   Aj = LMO_IAO(imo+1:imo+nmoi,j,:)
                   ! Compute the charge matrix elements of fragment k
                   call dot_prod(Ai,Ai,res)
                   Qii = res(1)
                   call dot_prod(Aj,Aj,res)
                   Qjj = res(1)
                   call dot_prod(Ai,Aj,Qij)
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                   ! Pipek-Mezey localization for complex orbitals
                     Bij(ij,:) = Bij(ij,:) + 4.D0*Qij(:)*(Qii - Qjj)
                   else if (Expnt.eq.4) then
                   ! Appendix D of Knizia's paper adapted for complex orbitals
                     Bij(ij,:) = Bij(ij,:) + 4.D0*Qij(:)*(Qii*Qii*Qii - Qjj*Qjj*Qjj)
                   else
                     write(luout,*)"STOP: exponent should be 2, 3 or 4."
                     stop
                   end if
                   deallocate(Ai,Aj)
                   imo = imo + nmoi
                enddo

                ! Use of a damping parameter to achieve convergence (but can also slow it down).
                ! This works only for real algebra. To be extended.
                if (rcq .eq. 1) then
                   Bij(ij,1) = damping*Bij(ij,1) + (1.D0 - damping)*Bij_previous(ij)
                   Bij_previous(ij) = Bij(ij,1)
                endif

                ! We first calculate the (real,complex or quaternion) phase of Bij
                ! In this way we can always apply the equation for real orbitals,
                ! by simply scaling  orbital |j> with the conjugate of this phase factor or,
                ! equivalently, orbital |i> with the phase factor.
                absBij = euc_norm(Bij(ij,:))
                if (absBij > 1.0D-15) then
                   Pij = Bij(ij,:) / absBij
                else
                   Pij    = 0.D0
                   Pij(1) = 1.D0
                endif

                ! Aij is the second derivative at phi=0 (see Eq.(8) of [2]).
                Aij = 0.D0
                imo = 0
                ! Increment Aij and Bij for each fragment:
                do k=1,n_fragments
                   nmoi = nmo_frag(k)
                   allocate(Ai(nmoi,rcq))
                   allocate(Aj(nmoi,rcq))
                   Ai = LMO_IAO(imo+1:imo+nmoi,i,:)
                   Aj = LMO_IAO(imo+1:imo+nmoi,j,:)
                   ! Compute the charge matrix elements of fragment k
                   call dot_prod(Ai,Ai,res)
                   Qii = res(1)
                   call dot_prod(Aj,Aj,res)
                   Qjj = res(1)
                   call dot_prod(Ai,Aj,Qij)
                   ! compute Re(Qij * Pij^*)
                   if (rcq .eq. 1) then
                     reQij = Qij(1) * Pij(1)
                   else if (rcq .eq. 2) then
                     reQij = Qij(1) * Pij(1) + Qij(2) * Pij(2)
                   else if (rcq .eq. 4) then
                     reQij = Qij(1) * Pij(1) + Qij(2) * Pij(2) + Qij(3) * Pij(3) + Qij(4) * Pij(4)
                   end if
                   if ((Expnt.eq.2) .or. (Expnt.eq.3)) then
                   ! Pipek-Mezey localization for complex orbitals
                     Aij = Aij + 4.D0*reQij**2.D0 - (Qii-Qjj)*(Qii-Qjj)
                   else if (Expnt.eq.4) then
                   ! Appendix D of Knizia's paper adapted for complex orbitals
                     Aij = Aij - Qii*Qii*Qii*Qii - Qjj*Qjj*Qjj*Qjj + &
                      & 6.D0*(Qii*Qii + Qjj*Qjj)*reQij**2.D0 + Qjj*Qii*Qii*Qii + &
                      & Qii*Qjj*Qjj*Qjj
                   else
                     write(luout,*)"STOP: exponent should be 2, 3 or 4."
                     stop
                   end if
                   deallocate(Ai,Aj)
                   imo = imo + nmoi
                enddo

                !Bij is the actual gradient. Aij is effectively
                !the second derivative at phi=0.

                ! Now perform the rotation:
                if ((Aij*Aij + absBij**2.D0) > 1.0D-15) then

                  fGrad = fGrad + absBij**2.D0
                  phi = 0.25D0*datan2(absBij,-Aij)

                  allocate(Ai(nmo2,rcq))
                  allocate(Aj(nmo2,rcq))
                  allocate(Aj0(nmo2,rcq))
                  ! we need a copy of the original LMO_IAO(:,j,:) as this is overwritten in the first update
                  Aj0 = LMO_IAO(:,j,:) 

                  ! in the update of |j> we first need to multiply |i> with Pij.
                  ! This makes Bij real and positive as was assumed in the calculation of phi with absBij.
                  Aj = LMO_IAO(:,j,:) 
                  if (rcq .eq. 1) then
                     Ai = LMO_IAO(:,i,:)*Pij(1)
                  else if (rcq .eq. 2) then
                     Ai(:,1) = LMO_IAO(:,i,1)*Pij(1) - LMO_IAO(:,i,2)*Pij(2)
                     Ai(:,2) = LMO_IAO(:,i,2)*Pij(1) + LMO_IAO(:,i,1)*Pij(2)
                  else if (rcq .eq. 4) then
                     Ai(:,1) = LMO_IAO(:,i,1)*Pij(1) - LMO_IAO(:,i,2)*Pij(2) - LMO_IAO(:,i,3)*Pij(3) - LMO_IAO(:,i,4)*Pij(4)
                     Ai(:,2) = LMO_IAO(:,i,2)*Pij(1) + LMO_IAO(:,i,1)*Pij(2) + LMO_IAO(:,i,3)*Pij(4) - LMO_IAO(:,i,4)*Pij(3)
                     Ai(:,3) = LMO_IAO(:,i,3)*Pij(1) + LMO_IAO(:,i,1)*Pij(3) + LMO_IAO(:,i,4)*Pij(2) - LMO_IAO(:,i,2)*Pij(4)
                     Ai(:,4) = LMO_IAO(:,i,4)*Pij(1) + LMO_IAO(:,i,1)*Pij(4) + LMO_IAO(:,i,2)*Pij(3) - LMO_IAO(:,i,3)*Pij(2)
                  end if
                  LMO_IAO(:,j,:) =-dsin(phi)*Ai + dcos(phi)*Aj
 
                  ! in the update of |i> we first multiply |j> with Pij^*.
                  ! This again makes Bij real and positive as was assumed in the calculation of phi with absBij.
                  Ai = LMO_IAO(:,i,:)
                  if (rcq .eq. 1) then
                     Aj = Aj0*Pij(1)
                  else if (rcq .eq. 2) then
                     Aj(:,1) = Aj0(:,1)*Pij(1) + Aj0(:,2)*Pij(2)
                     Aj(:,2) = Aj0(:,2)*Pij(1) - Aj0(:,1)*Pij(2)
                  else if (rcq .eq. 4) then
                     Aj(:,1) = Aj0(:,1)*Pij(1) + Aj0(:,2)*Pij(2) + Aj0(:,3)*Pij(3) + Aj0(:,4)*Pij(4)
                     Aj(:,2) = Aj0(:,2)*Pij(1) - Aj0(:,1)*Pij(2) - Aj0(:,3)*Pij(4) + Aj0(:,4)*Pij(3)
                     Aj(:,3) = Aj0(:,3)*Pij(1) - Aj0(:,1)*Pij(3) - Aj0(:,4)*Pij(2) + Aj0(:,2)*Pij(4)
                     Aj(:,4) = Aj0(:,4)*Pij(1) - Aj0(:,1)*Pij(4) - Aj0(:,2)*Pij(3) + Aj0(:,3)*Pij(2)
                  end if
                  LMO_IAO(:,i,:) = dcos(phi)*Ai + dsin(phi)*Aj
 
                  deallocate(Ai,Aj,Aj0)
                  ! condition for the maximum in Ref. [3].
                  if (Expnt .eq. 2) then
                     max_condition = 16.D0*Aij*dcos(4.D0*phi) - 16.D0*absBij*dsin(4.D0*phi)
                     if (max_condition > 0.D0) then
                     write(luout,*)"Warning: The condition for the maximum is not fulfilled &
                              & (i.e. max_condition is not < 0) for orbital pair i,j = ",i,j
                     write(luout,'(A20,ES12.4)') "  max_condition = ", max_condition
                     endif
                  endif
                endif
             enddo
          enddo

          fGrad = (fGrad**0.5D0)/(nact*(nact-1)/2)
          write(luout,"(I5,F16.6,ES12.2)") it,func,fGrad
          if (abs(fGrad) .lt. 1.0D-15) then
             conv = .true.
             exit
          endif
        enddo
        deallocate(Pij,Qij,res)

        write(luout,*)
        write(luout,*)"******************************************"
        if (.not. conv) then
           write(luout,*) "No convergence reached, value of fGrad = ", fGrad
        else
           write(luout,*) "Final converged values (it, func, fGrad):"
           write(luout,"(I5,F16.8,ES12.4)") it,func,fGrad
        endif
        write(luout,*)"******************************************"
        write(luout,*)

   end subroutine Localization

   subroutine get_localization_func(LMO_IAO,nmo_frag,Expnt,func)
        ! Eq.(3) of [2]
        ! Actually, because the IAOs are orthonormalized, it reduces to
        ! Eq.(4) of [1].

        use linear_algebra

        implicit none
        integer                    :: i,j,imo,nmoi,n_fragments,nmo2,nact,rcq
        integer, intent(in)        :: Expnt
        integer, intent(in)        :: nmo_frag(:)    ! Size of the fragments
        real(8), intent(in)        :: LMO_IAO(:,:,:) ! The localized occupied MO in terms of IAOs
        real(8), intent(out)       :: func
        real(8), allocatable       :: vec(:,:)
        real(8), allocatable       :: population(:)

        rcq  = size(LMO_IAO,3)
        nmo2 = size(LMO_IAO,1)
        nact = size(LMO_IAO,2)
        n_fragments = size(nmo_frag)
        func = 0.D0
        imo = 0
        allocate(population(rcq))
        do i=1,n_fragments
           nmoi = nmo_frag(i)
           allocate(vec(nmoi,rcq))
           do j=1,nact
              vec = LMO_IAO(imo+1:imo+nmoi,j,:)
              call dot_prod(vec,vec,population)
              func = func + population(1)**dble(Expnt)
           enddo
           deallocate(vec)
           imo = imo + nmoi
        enddo
        deallocate(population)
   end subroutine get_localization_func

   subroutine sort_on_fragments (LMO_IAO,nmo_frag,nlo_frag,frag_bias)
   ! First assign and then sort the localized orbitals in terms of fragments 

        use linear_algebra

        implicit none

        integer, intent(in)        :: nmo_frag(:)       ! Size of the fragments (#IAOs on each)
        real(8), intent(in)        :: frag_bias(:)      ! Bias fragments when assigning LMOs (for non-polar bonds with nearly equal occupations)
        integer, intent(out)       :: nlo_frag(:)       ! Counts how many LMOs were assigned to each fragment
        real(8), intent(inout)     :: LMO_IAO(:,:,:)    ! Localized MO in terms of IAOs, will get sorted
        real(8), allocatable       :: LMO_sorted(:,:,:) ! Scratch array to sort LMO_IAO
        integer                    :: nmo2, nact, rcq, n_fragments, i_frag
        integer                    :: i, j, nmoi, imo
        real(8), allocatable       :: population(:), pop(:)
        integer, allocatable       :: map(:)

        nmo2 = size(LMO_IAO,1)
        nact = size(LMO_IAO,2)
        rcq  = size(LMO_IAO,3)
        n_fragments = size(nmo_frag)
        allocate(population(rcq))
        allocate(pop(n_fragments))
        allocate(map(nact))
        map = 0
 
        nlo_frag = 0
        do j=1,nact
           imo = 0
           do i=1,n_fragments
              nmoi = nmo_frag(i)
              call dot_prod(LMO_IAO(imo+1:imo+nmoi,j,:),LMO_IAO(imo+1:imo+nmoi,j,:),population)
              pop(i) = population(1)
              imo = imo + nmoi
           enddo
           i_frag = maxloc(pop+frag_bias,1) ! fragment to which orbital j is assigned
           map(j) = i_frag
           nlo_frag(i_frag) = nlo_frag(i_frag) + 1
        enddo

        allocate(LMO_sorted(nmo2,nact,rcq))
        i = 0
        do i_frag = 1,n_fragments
           do j=1,nact
              if (map(j) == i_frag) then
                 i = i + 1
                 LMO_sorted(:,i,:) = LMO_IAO(:,j,:)
              end if
           enddo
        enddo

        LMO_IAO = LMO_sorted

        deallocate(pop,map,population,LMO_sorted)

   end subroutine sort_on_fragments

   subroutine recanonize(LMO,CMO_energy,LMO_energy,nlo_frag)
   ! compute and recanonize diagonal blocks of the Fock matrix in LMO basis
   ! original energies are unmodified, but LMO will have recanonized orbitals on return

        use linear_algebra
        use rose_utils
        use properties
        use check

        implicit none

        integer, intent(in)        :: nlo_frag(:)    ! Number of local orbitals per fragment
        real(8), intent(in)        :: CMO_energy(:)  ! CMO energies
        real(8), intent(inout)     :: LMO_energy(:)  ! LMO energies (obtained by recanonization)
        real(8), intent(inout)     :: LMO(:,:,:)     ! Localized MOs in terms of original MOs.
        integer                    :: n_fragments, nmo1, nmo2, rcq
        integer                    :: i, j, i_frag 
        real(8), allocatable       :: LMO_bar(:,:,:), LMO_Fock(:,:,:)
        real(8), allocatable       :: HF_Fock(:,:,:), LMO_can(:,:,:)

        nmo1 = size(LMO,1)
        nmo2 = size(LMO,2)
        rcq  = size(LMO,3)
        n_fragments = size(nlo_frag)

        ! Compute the Fock matrix in LMO basis
        allocate(HF_Fock(nmo1,nmo1,rcq))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = CMO_energy(i)
        allocate(LMO_Fock(nmo2,nmo2,rcq))
        allocate(LMO_bar(nmo2,nmo1,rcq))
        call transpose_conjg(LMO,LMO_bar)
        call matrix_mul(LMO_bar,HF_Fock,LMO,LMO_Fock)

        allocate(LMO_can(nmo1,nmo2,rcq))
        i = 1
        j = 0
        do i_frag = 1, n_fragments
           j = j + nlo_frag(i_frag)
           if (j >= i) then
               ! eigenvectors will be returned in LMO_Fock and define the canonization
               call diag_matrix(LMO_Fock(i:j,i:j,:),LMO_energy(i:j))
               ! add zero's outside the fragment space
               LMO_Fock(1:i-1,i:j,:) = 0.D0
               LMO_Fock(j+1:nmo2,i:j,:) = 0.D0
               ! transform the LMOs to this (semi)canonical basis
               call matrix_mul(LMO,LMO_Fock,LMO_can)
           end if 
           i = i + nlo_frag(i_frag)
        end do
        
        LMO = LMO_can
        deallocate(HF_Fock,LMO_bar,LMO_Fock,LMO_can)

   end subroutine recanonize

   subroutine sort_LMO_by_energy(LMO,CMO_energy,LMO_energy,restricted)
   ! sort the localized molecular orbitals

        use linear_algebra
        use rose_utils
        use properties
        use check

        implicit none

        real(8), intent(inout)     :: LMO(:,:,:)
        real(8), intent(in)        :: CMO_energy(:)  ! CMO energies
        real(8), intent(inout)     :: LMO_energy(:)  ! LMO energies
        logical, intent(in)        :: restricted
        integer                    :: i
        integer                    :: nmo1, nmo2, rcq
        real(8), allocatable       :: HF_Fock(:,:,:), LMO_Fock(:,:,:), LMO_bar(:,:,:)

        nmo1 = size(LMO,1)
        nmo2 = size(LMO,2)
        rcq  = size(LMO,3)

        ! Compute the Fock matrix in LMO basis
        allocate(HF_Fock(nmo1,nmo1,rcq))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = CMO_energy(i)
        allocate(LMO_Fock(nmo2,nmo2,rcq))
        allocate(LMO_bar(nmo2,nmo1,rcq))
        call transpose_conjg(LMO,LMO_bar)
        call matrix_mul(LMO_bar,HF_Fock,LMO,LMO_Fock)

        forall (i=1:nmo2) LMO_energy(i) = LMO_Fock(i,i,1)

        ! Compute the energies (diagonal elements of Fock matrix in LMO basis) and sort
        call sort_by_ascending_eigenvalues(LMO,LMO_energy,restricted)

        deallocate(HF_Fock,LMO_Fock,LMO_bar)

   end subroutine sort_LMO_by_energy

   subroutine make_LMO_Fock(IAO,LMO,LMO_IAO,CMO_energy,nmo_frag,occupied)
   ! compute and write the Fock matrix in LMO basis

        use linear_algebra
        use rose_utils
        use properties
        use check

        implicit none

        integer, intent(in)        :: nmo_frag(:)    ! Size of the fragments
        real(8), intent(in)        :: IAO(:,:,:)
        real(8), intent(in)        :: CMO_energy(:)  ! CMO energies
        real(8), intent(inout)     :: LMO_IAO(:,:,:) ! Localized MOs in terms of IAOs
        real(8), intent(inout)     :: LMO(:,:,:)     ! Localized MOs in terms of original MOs.
        logical                    :: occupied
        integer                    :: i, j, file1, file2
        integer                    :: n_fragments, nact, nmo1, nmo2, rcq
        real(8), allocatable       :: LMO_bar(:,:,:), LMO_Fock(:,:,:), canonical_energies(:)
        real(8), allocatable       :: HF_Fock(:,:,:)
        character(4)               :: extension
        logical                    :: tobe

        nmo1 = size(IAO,1)
        nmo2 = size(IAO,2)
        rcq  = size(IAO,3)
        nact = size(LMO_IAO,2)
        n_fragments = size(nmo_frag)

        ! Compute the Fock matrix in LMO basis
        allocate(HF_Fock(nmo1,nmo1,rcq))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = CMO_energy(i)
        allocate(LMO_Fock(nact,nact,rcq))
        allocate(LMO_bar(nact,nmo1,rcq))
        call transpose_conjg(LMO,LMO_bar)
        call matrix_mul(LMO_bar,HF_Fock,LMO,LMO_Fock)

        ! TEST: write the Fock matrix in the ILMO basis, and the localization trasnformation matrix (LMO)
        call get_free_fileunit(file1)
        call get_free_fileunit(file2)
        if (occupied) extension = "_occ"
        if (.not. occupied) extension = "_vir"
        inquire (file="ILMO_Fock"//extension//".dat",exist=tobe)
        if (tobe) then
         print*, "warning: overwriting file ILMO_Fock"//extension//".dat"
         open (file1,file="ILMO_Fock"//extension//".dat",status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
        else
         open (file1,file="ILMO_Fock"//extension//".dat",status='NEW',FORM='FORMATTED',access='SEQUENTIAL')
        end if
        do i=1,nact
         write(file1,'(1000ES10.2)') (LMO_Fock(j,i,1),j=1,nact)
        enddo
        close(file1)
        inquire (file="ILMO_transformation"//extension//".dat",exist=tobe)
        if (tobe) then
         print*, "warning: overwriting file ILMO_transformation"//extension//".dat"
         open (file2,file="ILMO_transformation"//extension//".dat",status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
        else
         open (file2,file="ILMO_transformation"//extension//".dat",status='NEW',FORM='FORMATTED',access='SEQUENTIAL')
        end if
        do i=1,nmo1
         write(file2,'(1000ES10.2)') (LMO(i,j,1),j=1,nact)
        enddo
        close(file2)

        ! Check for occupied:
        if (occupied) then
        write(luout,*)
        write(luout,*)"----- Diagonalization of the Fock matrix to check the recovering of the occupied canonical MOs -----"
        write(luout,*)
        allocate(canonical_energies(nact))
        call diag_matrix(LMO_Fock,canonical_energies)
        write(luout,'(A12,A20)')"MO","Energy"
        do i=1,nact
           write(luout,'(I12,F20.12)') i, canonical_energies(i)
        end do
        endif
        
        deallocate(HF_Fock,LMO_bar,LMO_Fock)

   end subroutine make_LMO_Fock

end module make_ibo
