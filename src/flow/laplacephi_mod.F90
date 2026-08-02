
MODULE laplacephi_mod
    USE core_mod

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE laplacephi_backend(res, phi, gsaw, gsae, gsas, gsan, &
                gsab, gsat, gsap, bp, mygrids_,  kkk_, jjj_, iii_, &
                ip3d_, ip1dx_, ip1dy_, ip1dz_) &
                BIND(C, name="laplacephi_c")
            USE precision_mod, ONLY: realk, intk
            REAL(realk), INTENT(inout) :: res(:)
            REAL(realk), INTENT(in) :: phi(:)
            REAL(realk), INTENT(in) :: gsaw(:)
            REAL(realk), INTENT(in) :: gsae(:)
            REAL(realk), INTENT(in) :: gsas(:)
            REAL(realk), INTENT(in) :: gsan(:)
            REAL(realk), INTENT(in) :: gsab(:)
            REAL(realk), INTENT(in) :: gsat(:)
            REAL(realk), INTENT(in) :: gsap(:)
            REAL(realk), INTENT(in) :: bp(:)
            INTEGER(intk), INTENT(in) :: mygrids_(:)
            INTEGER(intk), INTENT(in) :: kkk_(:)
            INTEGER(intk), INTENT(in) :: jjj_(:)
            INTEGER(intk), INTENT(in) :: iii_(:)
            INTEGER(intk), INTENT(in) :: ip3d_(:)
            INTEGER(intk), INTENT(in) :: ip1dx_(:)
            INTEGER(intk), INTENT(in) :: ip1dy_(:)
            INTEGER(intk), INTENT(in) :: ip1dz_(:)
        END SUBROUTINE laplacephi_backend
    END INTERFACE

    PUBLIC :: laplacephi, laplacephi_level

CONTAINS

    SUBROUTINE laplacephi(res_f, phi_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f

        ! Local variables
        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")
        CALL get_field(bp_f, "BP")

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("laplacephi")
#endif

#ifdef _MGLET_USE_BACKEND_
        CALL laplacephi_backend(res_f%arr, phi_f%arr, gsaw%arr, gsae%arr, &
            gsas%arr, gsan%arr, gsab%arr, gsat%arr, gsap%arr, bp_f%arr, &
            mygrids, kkk, jjj, iii, ip3d, ip1dx, ip1dy, ip1dz)
        !CALL laplacephi_impl(res_f%arr, phi_f%arr, gsaw%arr, gsae%arr, &
        !    gsas%arr, gsan%arr, gsab%arr, gsat%arr, gsap%arr, bp_f%arr)

#else
        CALL laplacephi_impl(res_f%arr, phi_f%arr, gsaw%arr, gsae%arr, &
            gsas%arr, gsan%arr, gsab%arr, gsat%arr, gsap%arr, bp_f%arr)
#endif

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
    END SUBROUTINE laplacephi


    SUBROUTINE laplacephi_impl(res, phi, aw, ae, as, an, ab, at, ap, bp)
        ! Subroutine arguments
        REAL(realk), CONTIGUOUS, INTENT(inout) :: res(:)
        REAL(realk), CONTIGUOUS, DIMENSION(:), INTENT(in) :: phi, aw, ae, as, &
            an, ab, at, ap, bp
        
        ! Local variables
        INTEGER(intk) :: i, igrid, kk, jj, ii, ip3, ipx, ipy, ipz

        !!$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ipx, &
        !!$omp& ipy, ipz)
        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            CALL laplacephi_grid(kk, jj, ii, res(ip3), phi(ip3), &
                aw(ipx), ae(ipx), an(ipy), as(ipy), at(ipz), ab(ipz), &
                ap(ip3), bp(ip3))
        END DO
        !!$omp end target teams distribute
    END SUBROUTINE laplacephi_impl


    SUBROUTINE laplacephi_level(ilevel, res_f, phi_f)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in):: ilevel
        TYPE(field_t), INTENT(inout) :: res_f
        TYPE(field_t), INTENT(in) :: phi_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, ip3, ipx, ipy, ipz

        TYPE(field_t), POINTER :: gsaw, gsae, gsas, gsan, gsab, gsat, gsap, bp_f

        CALL get_field(gsaw, "GSAW")
        CALL get_field(gsae, "GSAE")
        CALL get_field(gsas, "GSAS")
        CALL get_field(gsan, "GSAN")
        CALL get_field(gsab, "GSAB")
        CALL get_field(gsat, "GSAT")
        CALL get_field(gsap, "GSAP")
        ! BP for noib is 1.0. Take extra multiplications instead of branching
        ! in the kernel or duplicating code for now.
        CALL get_field(bp_f, "BP")

        ASSOCIATE ( &
            phi => phi_f%arr, &
            res => res_f%arr, &
            ap  => gsap%arr, &
            aw  => gsaw%arr, &
            ae  => gsae%arr, &
            an  => gsan%arr, &
            as  => gsas%arr, &
            at  => gsat%arr, &
            ab  => gsab%arr, &
            bp  => bp_f%arr)

#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_push("laplacephi_level")
#endif
        ! !$omp target teams distribute private(i, igrid, kk, jj, ii, ip3, ipx, &
        ! !$omp& ipy, ipz)
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_ip3(ip3, igrid)
            CALL get_ip1x(ipx, igrid)
            CALL get_ip1y(ipy, igrid)
            CALL get_ip1z(ipz, igrid)

            CALL laplacephi_grid(kk, jj, ii, res(ip3), phi(ip3), &
                aw(ipx), ae(ipx), an(ipy), as(ipy), at(ipz), ab(ipz), &
                ap(ip3), bp(ip3))
        END DO
        ! !$omp end target teams distribute
#ifdef _MGLET_PROFILE_ANNOTATIONS_
        CALL profile_range_pop()
#endif
        END ASSOCIATE
    END SUBROUTINE laplacephi_level


    PURE SUBROUTINE laplacephi_grid(kk, jj, ii, res, phi, aw, ae, an, as, &
            at, ab, ap, bp)
        !$omp declare target
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: res(kk, jj, ii)
        REAL(realk), INTENT(in) :: phi(kk, jj, ii)
        REAL(realk), INTENT(in) :: aw(ii), ae(ii), an(jj), as(jj), &
            at(kk), ab(kk)
        REAL(realk), INTENT(in) :: ap(kk, jj, ii), bp(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        !!$omp parallel do collapse(3) private(i, j, k)
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    res(k, j, i) = &
                        - aw(i)*phi(k, j, i-1)*bp(k, j, i-1)*bp(k, j, i) &
                        - ae(i)*phi(k, j, i+1)*bp(k, j, i)*bp(k, j, i+1) &
                        - as(j)*phi(k, j-1, i)*bp(k, j-1, i)*bp(k, j, i) &
                        - an(j)*phi(k, j+1, i)*bp(k, j, i)*bp(k, j+1, i) &
                        - ab(k)*phi(k-1, j, i)*bp(k-1, j, i)*bp(k, j, i) &
                        - at(k)*phi(k+1, j, i)*bp(k, j, i)*bp(k+1, j, i) &
                        - ap(k, j, i)*phi(k, j, i)
                END DO
            END DO
        END DO
        !!$omp end parallel do
    END SUBROUTINE laplacephi_grid
END MODULE laplacephi_mod
