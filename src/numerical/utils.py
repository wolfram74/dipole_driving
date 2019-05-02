import numpy

def read_column(rank_2_tensor, col_num):
    return [row[col_num] for row in rank_2_tensor]

def rk45_steps(gradient_function, t, state, step_size):
    # http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    kernel1 = step_size*gradient_function(t, state + 0)
    kernel2 = step_size*gradient_function(
        t+step_size/4.,
        state + kernel1/4.
        )
    kernel3 = step_size*gradient_function(
        t+3.*step_size/8.,
        state + 3.*kernel1/32. + 9.*kernel2/32.
        )
    kernel4 = step_size*gradient_function(
        t+12.*step_size/13.,
        state + 1932.*kernel1/2197. - 7200.*kernel2/2197. + 7296.*kernel3/2197.
        )
    kernel5 = step_size*gradient_function(
        t+step_size/1.,
        state + 439.*kernel1/216. - 8.*kernel2 + 3680.*kernel3/513. + 845.*kernel4/4104.
        )
    kernel6 = step_size*gradient_function(
        t+step_size/2.,
        state - 8.*kernel1/27. + 2.*kernel2 - 3544.*kernel3/2565.
        + 1859.*kernel4/4104. - 11.*kernel5/40.
        )
    delta1 = (
        25.*kernel1/216.
        +1408.*kernel3/2565.
        +2197.*kernel4/4101.
        -1.*kernel5/5.
        )
    delta2 = (
        16.*kernel1/135.
        +6656.*kernel3/12825.
        +28561.*kernel4/56430.
        -9.*kernel5/50.
        +2.*kernel6/55.
        )
    return delta1, delta2

def RK45_path(
        gradient_function, state, start_time, end_time,
        precision=10.**-4, step_size=10.**-3, max_steps=10**4
    ):
    pass

def step_scaler(precision, step_size, delta1, delta2):
    rescale = 10.**6
    for ind in range(len(delta1)):
        proposed = (precision*step_size/(
            2*abs(delta2[ind]-delta1[ind])
            ))**.25
        print(proposed)
        if proposed < rescale:
            rescale = proposed
    return rescale
