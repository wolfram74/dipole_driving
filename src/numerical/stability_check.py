import utils
import numpy
from matplotlib import pyplot

def perturbed_quantity():
    gam1 = numpy.array([
        0.,
        -9.80118581243531,
        -5.54970354689117,
        6.95112240577761,
        0.,
        -0.141001484944198,
        0.282002969888397,
        1.0
    ])
    gam1 = gam1/numpy.linalg.norm(gam1)
    # print(gam1)
    gam1/=10.**6.
    # print(gam1)
    return gam1

def perturbed_quantity2():
    gam1 = numpy.array([
        0,
        -43.532147520898,
        2.8830368802245,
        0,
        0,
        0,
        0,
        1.0
    ])
    gam1 = gam1/numpy.linalg.norm(gam1)
    # print(gam1)
    gam1/=10.**3.
    # print(gam1)
    return gam1



def perturbed_quantity3():
    gam1 = numpy.array([
        0.,
        0.,
        0.,
        1.,
        0.,
        0.,
        0.,
        0.,
    ])
    gam1 = gam1/numpy.linalg.norm(gam1)
    # print(gam1)
    gam1/=10.**2.
    # print(gam1)
    return gam1



def circular_quantity(r0):
    ptht0 = 2.**(-.5)*r0**(-.5)
    pphit0 = 2.**(.5)*r0**(-2.5)/10.
    gam0 = numpy.array([
        0.,0.,0.,r0,
        0., pphit0, ptht0, 0.
        ])
    return gam0

def plot_all_variables(t, states):
    figure, subplots = pyplot.subplots(2,4)
    exp1 = 0.287723317652649
    fits = []
    figure.suptitle('residues of $\\phi_d$, $\\phi_t$, $\\theta$ and r. associated momenta plotted below. Green curve denotes expected curve from exponential growth.')
    for ti in t:
        row = []
        for ind in range(len(states[0])):
            row.append([states[0][ind]*numpy.exp(exp1*ti)])
        fits.append(row)
    # print(len(t), len(fits), utils.read_column(fits, 0))
    for i in range(4):
        subplots[0][i].plot(t, utils.read_column(states,i))
        # subplots[0][i].plot(t, utils.read_column(fits,i))
        subplots[1][i].plot(t, utils.read_column(states,i+4))
    pyplot.show()

def calc_residues(path_out,gam0):
    output = []
    Omeg0 = 20.*gam0[5]
    for step in path_out:
        state = step[1]
        t= step[0]
        circle = numpy.array([
            0.,Omeg0*t,Omeg0*t/2,gam0[3],
            0., gam0[5], gam0[6], 0.
        ])
        output.append(state-circle)
    return output

def main():
    gam1 = perturbed_quantity2()
    gam0 = circular_quantity(2.)
    print(gam1)
    print(gam0)
    state0 = gam0
    state0 = gam0+gam1
    print(state0)
    path_out = utils.RK45_simple_path(
        utils.reduced_dipole_equations, state0, 0., 20.,
        max_steps=1000, precision=10**-9)
    t = [el[0] for el in path_out]
    states = [el[1] for el in path_out]
    residues = calc_residues(path_out, gam0)
    # print(states[-1])
    # print(residues[0])
    # print(residues[-1])
    # plot_all_variables(t, states)
    plot_all_variables(t, residues)

if __name__=='__main__':
    main()

