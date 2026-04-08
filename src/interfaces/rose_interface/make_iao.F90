module make_iao

contains

   subroutine MakeIAO_Standard_2013(A,C,S11,S21,S22,P12)
        use linear_algebra
        use rose_utils

        integer                :: i, nmo1, nmo2, nocc, rcq
        real(8), intent(in)    :: S11(:,:,:), S21(:,:,:), S22(:,:,:)
        real(8), intent(in)    :: P12(:,:,:)
        real(8), intent(in)    :: C(:,:,:)
        real(8), intent(inout) :: A(:,:,:)
        real(8), allocatable   :: c_tilde(:,:,:),P21(:,:,:),Id(:,:,:)
        real(8), allocatable   :: term11(:,:,:),term12(:,:,:),term21(:,:,:),term22(:,:,:),term1(:,:,:),term2(:,:,:)
        real(8), allocatable   :: Cbar(:,:,:),c_tilde_bar(:,:,:),CtCtbarS11(:,:,:),term22tmp(:,:,:)

        if (size(A,3) .ne. size(C,3)) then
           print*,"make_iao: A and C of different algebra."
           stop
        else if (size(A,3) .ne. size(S11,3)) then
           print*,"make_iao: A and S11 of different algebra."
           stop
        else if (size(A,3) .ne. size(S21,3)) then
           print*,"make_iao: A and S21 of different algebra."
           stop
        else if (size(A,3) .ne. size(S22,3)) then
           print*,"make_iao: A and S22 of different algebra."
           stop
        else if (size(A,3) .ne. size(P12,3)) then
           print*,"make_iao: A and P12 of different algebra."
           stop
        end if

        rcq  = size(A,3)
        nmo1 = size(A,1)
        nmo2 = size(A,2)
        nocc = size(C,2)

        allocate (c_tilde(nmo1,nocc,rcq))
        allocate (P21(nmo2,nmo1,rcq))
        allocate (Id(nmo1,nmo1,rcq))
        call solve_lineq(S22,S21,P21)

        Id = 0.D0
        forall (i=1:nmo1) Id(i,i,1) = 1.D0
        call matrix_mul(P12,P21,C,c_tilde)
        call SymOrth(c_tilde,S11)
        allocate(term1(nmo1,nmo2,rcq))
        allocate(term2(nmo1,nmo2,rcq))
        allocate(term11(nmo1,nmo1,rcq))
        allocate(term12(nmo1,nmo2,rcq))
        allocate(term21(nmo1,nmo1,rcq))
        allocate(term22(nmo1,nmo2,rcq))
        allocate(term22tmp(nmo1,nmo1,rcq))
        allocate(Cbar(nocc,nmo1,rcq))
        allocate(c_tilde_bar(nocc,nmo1,rcq))
        allocate(CtCtbarS11(nmo1,nmo1,rcq))

        ! compute term11:
        call transpose_conjg(C,Cbar)
        call matrix_mul(C,Cbar,S11,term11)
        ! compute term12:
        call transpose_conjg(c_tilde,c_tilde_bar)
        call matrix_mul(c_tilde,c_tilde_bar,S11,CtCtbarS11)
        call matrix_mul(CtCtbarS11,P12,term12)
        ! compute term21:
        term21 = Id - term11
        ! compute term22:
        term22tmp = Id - CtCtbarS11
        call matrix_mul(term22tmp,P12,term22)
        ! compute term1:
        call matrix_mul(term11,term12,term1)
        ! compute term2:
        call matrix_mul(term21,term22,term2)
        A = term1 + term2

        deallocate(P21,c_tilde,Id,term11,term12,term21,term22)
        deallocate(Cbar,c_tilde_bar,CtCtbarS11,term22tmp)

   end subroutine MakeIAO_Standard_2013

   subroutine MakeIAO_Simple_2014(A,C,S21,S22,P12)
        use linear_algebra

        implicit none 
        integer                :: ierr
        integer                :: nmo1, nmo2, nocc, rcq
        real(8), allocatable   :: coeff2_occ(:,:,:)
        real(8), allocatable   :: c_tilde(:,:,:),c_tilde_bar(:,:,:),c_tilde2_bar(:,:,:)
        real(8), allocatable   :: s_tilde(:,:,:),t4(:,:,:)
        real(8), intent(in)    :: S21(:,:,:),S22(:,:,:)
        real(8), intent(in)    :: P12(:,:,:)
        real(8), intent(in)    :: C(:,:,:)
        real(8), intent(inout) :: A(:,:,:)
        real(8), allocatable   :: A_X(:,:,:), P21(:,:,:), S_12(:,:,:), c_tilde_X(:,:,:)
        real(8), allocatable   :: c_tilde_X_bar(:,:,:), Cbar(:,:,:), CCbarS12(:,:,:)
        real(8), allocatable   :: P12ctilde(:,:,:), t4coeff2(:,:,:), coeff2_occ_bar(:,:,:)
        real(8), allocatable   :: P12CXCXbarS22(:,:,:)
        real(8)                :: norm

        if (size(A,3) .ne. size(C,3)) then
           print*,"make_iao: A and C of different algebra."
           stop
        else if (size(A,3) .ne. size(S21,3)) then
           print*,"make_iao: A and S21 of different algebra."
           stop
        else if (size(A,3) .ne. size(S22,3)) then
           print*,"make_iao: A and S22 of different algebra."
           stop
        else if (size(A,3) .ne. size(P12,3)) then
           print*,"make_iao: A and P12 of different algebra."
           stop
        end if

        rcq  = size(A,3)
        nmo1 = size(A,1)
        nmo2 = size(A,2)
        nocc = size(C,2)

        ! Make intrinsic basis orbitals for the occupied space (Simple/2014 scheme of Knizia's reference implementation)
        ! Note that our overlap matrices refer to the MO basis of the molecule and its fragments and that the coefficient
        ! matrices are simple diagonal matrices. we transform to AO basis after this step.
        allocate (coeff2_occ(nmo2,nocc,rcq))
        allocate (coeff2_occ_bar(nocc,nmo2,rcq))
        allocate (c_tilde(nmo2,nocc,rcq))
        allocate (c_tilde_bar(nocc,nmo2,rcq))
        allocate (c_tilde2_bar(nmo2,nocc,rcq))
        allocate (s_tilde(nocc,nocc,rcq))
        allocate (t4(nmo1,nocc,rcq))
        allocate (P12ctilde(nmo1,nocc,rcq))
        allocate (t4coeff2(nmo1,nmo2,rcq))

        call matrix_mul(S21,C,coeff2_occ)
        call solve_lineq(S22,coeff2_occ,c_tilde)
        call transpose_conjg(coeff2_occ,coeff2_occ_bar)
        call matrix_mul(coeff2_occ_bar,c_tilde,s_tilde)
        call transpose_conjg(c_tilde,c_tilde_bar)
        call solve_lineq(s_tilde,c_tilde_bar,c_tilde_bar)
        call transpose_conjg(c_tilde_bar,c_tilde2_bar)
        call matrix_mul(P12,c_tilde2_bar,P12ctilde)
        t4 = C - P12ctilde
        call matrix_mul(t4,coeff2_occ_bar,t4coeff2)
        A = P12 + t4coeff2

        ! Compare with reference result to check (just to check if transformation is okay). 
        ! I should get rid of this once we are sure the transformation is correct.
        print*
        print*," *** Check if Simple/2014 transformation is okay ***"
        allocate(A_X(nmo1,nmo2,rcq)) 
        allocate(P21(nmo2,nmo1,rcq)) 
        allocate(S_12(nmo1,nmo2,rcq)) 
        allocate(c_tilde_X(nmo2,nocc,rcq))
        allocate(c_tilde_X_bar(nocc,nmo2,rcq))
        allocate(Cbar(nocc,nmo1,rcq))
        allocate(CCbarS12(nmo1,nmo2,rcq))
        allocate(P12CXCXbarS22(nmo1,nmo2,rcq))
        call transpose_conjg(S21,S_12)
        call solve_lineq(S22,S21,P21)
        call matrix_mul(P21,C,c_tilde_X)
        call SymOrth(c_tilde_X,S22)
        ! Compute A_X:
        call transpose_conjg(C,Cbar)
        call matrix_mul(C,Cbar,S_12,CCbarS12)
        call transpose_conjg(c_tilde_X,c_tilde_X_bar)
        call matrix_mul(P12,c_tilde_X,c_tilde_X_bar,S22,P12CXCXbarS22)
        A_X = P12 + CCbarS12 - P12CXCXbarS22

        ! check:
        call diff_norm2_rcq(A,A_X,norm)
        if (norm > 1.D-12) then
          print*," *** norm2(A_X - A) = ",norm,&
                 ", Transformation NOT OKAY (or the test has failed... !) ***"
        else
          print*," *** norm2(A_X - A) = ",norm,&
                 ", Transformation OKAY ***"
        endif

        deallocate (coeff2_occ,c_tilde,c_tilde_bar,c_tilde2_bar,s_tilde,t4,A_X,P21,S_12,c_tilde_X)
        deallocate (coeff2_occ_bar,c_tilde_X_bar,t4coeff2,P12ctilde,Cbar,CCbarS12,P12CXCXbarS22)

   end subroutine MakeIAO_Simple_2014

   subroutine MakeIAO_Simple_2013(A,C,S11,S12,S21,S22,P12)
        use linear_algebra

        implicit none
        integer                :: ierr
        integer                :: nmo1, nmo2, nocc, rcq
        real(8), intent(inout) :: A(:,:,:)
        real(8), intent(in)    :: S11(:,:,:),S12(:,:,:),S21(:,:,:),S22(:,:,:)
        real(8), intent(in)    :: P12(:,:,:)
        real(8), intent(in)    :: C(:,:,:)
        real(8), allocatable   :: c_tilde(:,:,:),P21(:,:,:)
        real(8), allocatable   :: Cbar(:,:,:), c_tilde_bar(:,:,:), CCbarS12(:,:,:)
        real(8), allocatable   :: CtCtbarS12(:,:,:)

        if (size(A,3) .ne. size(C,3)) then
           print*,"make_iao: A and C of different algebra."
           stop
        else if (size(A,3) .ne. size(S11,3)) then
           print*,"make_iao: A and S11 of different algebra."
           stop
        else if (size(A,3) .ne. size(S12,3)) then
           print*,"make_iao: A and S12 of different algebra."
           stop
        else if (size(A,3) .ne. size(S21,3)) then
           print*,"make_iao: A and S21 of different algebra."
           stop
        else if (size(A,3) .ne. size(S22,3)) then
           print*,"make_iao: A and S22 of different algebra."
           stop
        else if (size(A,3) .ne. size(P12,3)) then
           print*,"make_iao: A and P12 of different algebra."
           stop
        end if

        rcq  = size(A,3)
        nmo1 = size(A,1)
        nmo2 = size(A,2)
        nocc = size(C,2)

        allocate (P21(nmo2,nmo1,rcq))
        allocate (c_tilde(nmo1,nocc,rcq))
        allocate (Cbar(nocc,nmo1,rcq))
        allocate (c_tilde_bar(nocc,nmo1,rcq))
        allocate (CCbarS12(nmo1,nmo2,rcq))
        allocate (CtCtbarS12(nmo1,nmo2,rcq))
        call solve_lineq(S22,S21,P21)
        call matrix_mul(P12,P21,C,c_tilde)
        call SymOrth(c_tilde,S11)
        call transpose_conjg(C,Cbar)
        call transpose_conjg(c_tilde,c_tilde_bar)
        call matrix_mul(C,Cbar,S12,CCbarS12)
        call matrix_mul(c_tilde,c_tilde_bar,S12,CtCtbarS12)
        A = P12 + CCbarS12 - CtCtbarS12

        deallocate (c_tilde,P21,Cbar,c_tilde_bar,CCbarS12,CtCtbarS12)

   end subroutine MakeIAO_Simple_2013

end module make_iao
