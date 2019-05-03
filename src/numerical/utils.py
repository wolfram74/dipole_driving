import numpy
sin = numpy.sin
cos = numpy.cos

def read_column(rank_2_tensor, col_num):
    return [row[col_num] for row in rank_2_tensor]

def rk45_steps(gradient_function, t, state, step_size):
    # http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    # errors in k5 and delta1 coeffs
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
    # print(kernel1)
    # print(kernel2)
    # print(kernel3)
    # print(kernel4)
    # print(kernel5)
    # print(kernel6)
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
            print('Did not finish, remaining time:: %f' % time_left)
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

def RK45_bouncing_path(
        gradient_function, state, start_time, end_time,
        precision=10.**-4, step_size=10.**-2, max_steps=10**4
    ):
    '''
    if
        radial momentum is negative
        and
        radial position less than the boundary+margin after an update
        then bounce
    if
        radial momentum is negative
        and
        radial position greater than the boundary
        and
        the linear interpolated position puts it past the boundary
        set the time step so it would put it at the margin
    '''
    path = [(start_time, state)]
    running = True
    last_loop = False
    time_left = end_time - start_time
    boundary = 1.
    margin = 10**-4
    while running:
        last_time = path[-1][0]
        last_state = path[-1][1]
        if last_loop:
            running = False
        if len(path) > max_steps:
            running = False
            print('Did not finish, remaining time:: %f' % time_left)
            return path

        deltO4, deltO5 = rk45_steps(
            gradient_function, last_time, last_state, step_size
            )
        step_adjust = step_scaler(precision, step_size, deltO4, deltO5)
        if step_adjust < .9:
            step_size *= step_adjust
            continue
        #decided steps are acceptably precise

        #normal force prevents acceleration inward
        # print(last_state)
        # if last_state[3]==1.:
        #     print(last_state[3] <= 1.0 , deltO5[7]<0. , last_state[7]==0.)
        if last_state[3] <= 1.0 and deltO5[7]<0. and last_state[7]==0.:
            # print('on shell hit')
            deltO5[3] = 0.
            deltO5[7] = 0.
        new_state = last_state+deltO5
        new_time = last_time+step_size
        #bounce
        if new_state[3] < (boundary+margin) and new_state[7] < 0:
            new_state[7] *= -1.
        path.append((new_time, new_state))

        if step_adjust > 1.25:
            step_adjust = 1.25
        step_size*= step_adjust

        if new_state[7] < 0:
            t_guess = (boundary+margin-new_state[3])/(2.*new_state[7])
            # print(t_guess, step_size)
            # print(new_time, new_state[3], new_state[7], deltO5[7])
            # print(t_guess, step_size, '\n')
            #t_guess will always be larger than needed
            if t_guess<step_size:
                step_size = t_guess

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
        print(diff)
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

def reduced_dipole_equations(t, state):
    '''
    state is:
        positions:
            0)phi_d, 1)phi_t, 2)theta, 3)r
        momenta:
            4)pd, 5)pt, 6)ptht, 7)pr
    '''
    deltas = numpy.zeros(len(state))
    #velocities
    deltas[0] = 20.* state[4]
    deltas[1] = 20.* state[5]
    deltas[2] = 2. * state[6]*(state[3]**(-2.))
    deltas[3] = 2. * state[7]
    #forces
    deltas[4] = -sin(state[0])/(12.*state[3]**3)
    deltas[5] = -sin(state[1] - 2.*state[2])/(4.*state[3]**3)
    deltas[6] = 2.*sin(state[1] - 2.*state[2])/(4.*state[3]**3)
    deltas[7] = (
        2.*state[6]**2/state[3]**3
        -(
            cos(state[0]) + 3.*cos(state[1] - 2.*state[2])
        )/(4.*state[3]**4)
    )
    #forbidden region
    # if state[3] <= 1. and deltas[7] <= 0.:
    #     deltas[7]=0.
    return deltas

def total_energy(state):
    return PE(state)+KE(state)

def total_L(state):
    return 2.*state[5]+state[6]

def PE(state):
    return -(
        cos(state[0]) + 3*cos(state[1] - 2.*state[2])
        )/(12.*state[3]**3)

def KE(state):
    return (
        2.*state[7]**2
        +2.*state[6]**2/state[3]**2
        +20.*(state[4]**2+state[5]**2)
        )/2.
