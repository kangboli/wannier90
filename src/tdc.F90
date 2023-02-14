module tdc
  !use w90_constants, only: dp
  implicit none
  integer, parameter   :: dp = kind(1.0d0)
  integer, allocatable :: bk_perm(:,:)
  real(dp)             :: eps = 1.0e-5_dp
  public
contains

  subroutine  setup(nntot,num_kpts,eps_)
    implicit none
    integer, intent(in) :: nntot, num_kpts
    real(dp), optional  :: eps_
    if (present(eps_)) eps = eps_
    if (.not.allocated(bk_perm)) allocate(bk_perm(nntot,num_kpts))
  end subroutine setup

  subroutine  compute_permutation_of_bk(nntot,num_kpts,bk)
    implicit none
    integer,  intent(in) :: nntot, num_kpts
    real(dp), intent(in) :: bk(3,nntot,num_kpts)
    integer  :: nkp, inn, jnn
    real(dp) :: dbk
    bk_perm = 0
    !$omp parallel do
    do inn = 1, nntot
      bk_perm(inn,1) = inn
    end do
    !$omp end parallel do
    !
    !$omp parallel do collapse(2)
    do nkp = 2, num_kpts
      do inn = 1, nntot
        do jnn = 1, nntot
          if (norm2(bk(:, jnn, nkp) - bk(:, inn, 1)) > eps) then
          else
            bk_perm(inn,nkp) = jnn
          end if
        end do ! jnn
      end do ! inn
    enddo
    !$omp end parallel do
  end subroutine compute_permutation_of_bk

  subroutine  debug_print_bk_perm(nntot,num_kpts,bk)
    implicit none
    integer,  intent(in) :: nntot, num_kpts
    real(dp), intent(in) :: bk(3,nntot,num_kpts)
    integer  :: nkp, inn, unt
    open(newunit=unt,file='tdc.debug.print_bk_perm_test',action='write')
    write(unt,*)"#nkp, inn, bk(1:3,inn,1) bk(1:3,bk_perm(inn,nkp),nkp)"
    do nkp = 2, num_kpts
      do inn = 1, nntot
        write(unt,'(2I5, 3F12.5," | ",3F12.5)') nkp,inn, bk(1:3,inn,1), bk(1:3,bk_perm(inn,nkp),nkp)
      end do ! inn
    enddo
    close(unt)
  end subroutine debug_print_bk_perm

  subroutine  cleanup()
    implicit none
    if (allocated(bk_perm)) deallocate(bk_perm)
  end subroutine cleanup
end module tdc
