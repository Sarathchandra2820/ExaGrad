module rose_davidson
    use iso_c_binding, only: c_double
    use molecule_t, only: fragment_scf_info

    implicit none

contains

    ! Replaces frag_scf_info(k)%C_mo with the fragment-local LMO block
    ! [c_lmo_occ_frag_k | c_lmo_vir_frag_k], shaped (nao, nocc_k + nvir_k).
    ! After this call, frag_scf_info(k)%C_mo(:,1:nocc_k) are the occupied LMOs
    ! and frag_scf_info(k)%C_mo(:,nocc_k+1:) are the virtual LMOs for fragment k.
    subroutine partition_rose_clmo(frag_scf_info, nfrags, nlo_frag_occ, nlo_frag_vir, c_lmo_occ, c_lmo_vir)
        implicit none
        type(fragment_scf_info), intent(inout) :: frag_scf_info(0:)
        integer,        intent(in) :: nfrags
        integer,        intent(in) :: nlo_frag_occ(:), nlo_frag_vir(:)
        real(c_double), intent(in) :: c_lmo_occ(:,:), c_lmo_vir(:,:)

        integer :: k, occ_off, vir_off, nao, nocc_k, nvir_k

        nao     = size(c_lmo_occ, 1)
        occ_off = 0
        vir_off = 0

        do k = 1, nfrags
            nocc_k = nlo_frag_occ(k)
            nvir_k = nlo_frag_vir(k)

            if (allocated(frag_scf_info(k)%C_mo)) deallocate(frag_scf_info(k)%C_mo)
            allocate(frag_scf_info(k)%C_mo(nao, nocc_k + nvir_k))

            frag_scf_info(k)%C_mo(:, 1:nocc_k) = &
                c_lmo_occ(:, occ_off+1 : occ_off+nocc_k)
            frag_scf_info(k)%C_mo(:, nocc_k+1 : nocc_k+nvir_k) = &
                c_lmo_vir(:, vir_off+1 : vir_off+nvir_k)

            frag_scf_info(k)%nocc = nocc_k
            frag_scf_info(k)%nvir = nvir_k

            occ_off = occ_off + nocc_k
            vir_off = vir_off + nvir_k
        end do

    end subroutine partition_rose_clmo

end module rose_davidson
