subroutine normalizedSpatialHet(SH, array)
    implicit none
    integer :: y1, y2, x1, x2
    integer :: y1i, y2i, x1i, x2i
    integer :: k, j, count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    real(kind=8), dimension(:, :), intent(in) :: array
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G
    real(kind=8) :: N
    integer, dimension(2) :: dim
    integer :: d1, d2

    integer :: total_permutes
    total_permutes = 0
    
    dim = shape(array)
    d1 = dim(1)
    d2 = dim(2)
    allocate(hstack_array(d1, 2*d2))
    allocate(v_hstack_array(2*d1, 2*d2))
    
    array_mean = sum(array)/max(1, size(array))
    G = 0
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [d1, 2*d2])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*d1, 2*d2]))
    N = array_mean*(1.5*d1*d2*(d1-1)*(d2-1) + d1*(d1-1) + d2*(d2-1))/((d1*(d1-1)+1)*(d2*(d2-1)+1))

    do y1 = 1, d2
        do y2 = 1, d2 
            y1i = y1
            y2i = y2

            if (y1 == y2) then
                cycle
            end if
            if (y2 < y1) then 
                y2i = y2i + d2
            end if

            subarray_mean = 0.0
            count = 0
            do k = y1i, y2i-1
                do j = 1, d1
                    ! note that loop structure is y then x so access array 
                    ! indices as j, k
                    subarray_mean = subarray_mean + v_hstack_array(j, k) 
                    count = count + 1
                end do
            end do
            subarray_mean = subarray_mean / count
            G = G + abs(subarray_mean - array_mean)
            !print *, 'lp1.1', abs(subarray_mean - array_mean), y1i-1, y2i-1

            do x1 = 1, d1
                do x2 = 1, d1 
                    x1i = x1
                    x2i = x2

                    if (x1 == x2) then
                        total_permutes = total_permutes + 1
                        cycle
                    end if
                    if (x2 < x1) then 
                        x2i = x2i + d1
                    end if

                    subarray_mean = 0.0
                    count = 0
                    do k = y1i, y2i-1
                        do j = x1i, x2i-1
                            subarray_mean = subarray_mean + v_hstack_array(k, j)
                            count = count + 1
                        end do
                    end do
                    subarray_mean = subarray_mean / count
                    G = G + abs(subarray_mean - array_mean)

                    total_permutes = total_permutes + 1


                    !print *, 'mean = ', subarray_mean, 'G = ', G
                    !print *, 'lp1.2', abs(subarray_mean - array_mean)
                end do
            end do
        end do
    end do
    !print *, 'permutes within first loop: ', total_permutes
    do x1 = 1, d1
        do x2 = 1, d1
            x1i = x1
            x2i = x2

            if (x1 == x2) then
                total_permutes = total_permutes + 1
                cycle
            end if
            if (x2 < x1) then 
                x2i = x2i + d1
            end if

            subarray_mean = 0.0
            count = 0
            do k = x1i, x2i-1
                do j = 1, d2
                    ! note that loop structure is x then y so access array 
                    ! indices as k, j
                    subarray_mean = subarray_mean + v_hstack_array(k, j)
                    count = count + 1
                end do
            end do
            subarray_mean = subarray_mean / count
            G = G + abs(subarray_mean - array_mean)
            !print *, 'lp2', abs(subarray_mean - array_mean)
            total_permutes = total_permutes + 1
        end do
    end do
    !print *, 'Total permutes: ', total_permutes
    !print *, 'N: ', N
    SH = G / (N*(d1*(d1-1)+1)*(d2*(d2-1)+1))
    deallocate(hstack_array)
    deallocate(v_hstack_array)
end subroutine normalizedSpatialHet
