import numpy

def read_column(rank_2_tensor, col_num):
    return [row[col_num] for row in rank_2_tensor]

def rk45_steps(gradient_function, t, state, step_size):
    # http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    # cross reference
    # https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
    # https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method#Butcher_tableau_for_Fehlberg%27s_4(5)_method
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
        state + 439.*kernel1/216. - 8.*kernel2 + 3680.*kernel3/513. - 845.*kernel4/4104. #error found
        )
    kernel6 = step_size*gradient_function(
        t+step_size/2.,
        state - 8.*kernel1/27. + 2.*kernel2 - 3544.*kernel3/2565.
        + 1859.*kernel4/4104. - 11.*kernel5/40.
        )
    delta1 = (
        25.*kernel1/216.
        +1408.*kernel3/2565.
        +2197.*kernel4/4104. #error found
        -1.*kernel5/5.
        )
    delta2 = (
        16.*kernel1/135.
        +6656.*kernel3/12825.
        +28561.*kernel4/56430.
        -9.*kernel5/50.
        +2.*kernel6/55.
        ) #O
    return delta1, delta2

def RK45_simple_path(
        gradient_function, state, start_time, end_time,
        precision=10.**-4, step_size=10.**-2, max_steps=10**4
    ):
    path = [(start_time, state)]
    running = True
    last_loop = False
    time_left = end_time - start_time
    while running:
        last_time = path[-1][0]
        last_state = path[-1][1]
        if last_loop:
            running = False
        if len(path) > max_steps:
            running = False
            print('Did non finish, remaining time:: %f' % time_left)
            return path

        deltO4, deltO5 = rk45_steps(
            gradient_function, last_time, last_state, step_size
            )
        step_adjust = step_scaler(precision, step_size, deltO4, deltO5)
        if step_adjust < .9:
            step_size *= step_adjust
            continue
        new_state = last_state+deltO5
        new_time = last_time+step_size
        path.append((new_time, new_state))

        if step_adjust > 1.25:
            step_adjust = 1.25
        step_size*= step_adjust

        time_left = end_time-new_time
        if time_left < step_size:
            print('finishing early, steps= %d' % (len(path)+1))
            step_size = time_left
            last_loop = True
    return path


def step_scaler(precision, step_size, delta1, delta2):
    rescale = 10.**6
    numer = precision*step_size
    for ind in range(len(delta1)):
        diff = abs(delta2[ind]-delta1[ind])
        if diff ==0:
            proposed = 1
            continue
        proposed = (
            numer/(2.*diff)
        )**.25
        # print(proposed)
        if proposed < rescale:
            rescale = proposed
    return rescale

