# test for returning to position of SHO
# test for amplitude insensitivity of SHO

import unittest
import numpy
import random
import utils

class AdaptiveRK45Test(unittest.TestCase):
    def test_meta(self):
        self.assertTrue(True)

    def test_2_steps_generated(self):
        stateI = numpy.array([0.,1.])
        step = 10**-4
        prec = 10**-6
        stepO4, stepO5 = utils.rk45_steps(base_SHO, 0, stateI, step_size=10**-4)
        self.assertEqual(len(stepO4), len(stepO5))
        adjust1 = utils.step_scaler(prec, step, stepO4, stepO5)
        step *= adjust1
        stepO4, stepO5 = utils.rk45_steps(base_SHO, 0, stateI, step_size=step)
        adjust2 = utils.step_scaler(prec, step, stepO4, stepO5)
        print(adjust1, adjust2)
        self.assertTrue(abs(numpy.log(adjust2))< abs(numpy.log(adjust1)))

    def test_period_check(self):
        stateI = numpy.array([0.,1.])
        print('pre starting')
        path_out = utils.RK45_simple_path(
            base_SHO, stateI, 0., 2*numpy.pi,
            max_steps=500, precision=10**-6)
        print(path_out[-1][0]-2*numpy.pi, path_out[-1][1])
        print(path_out[-2])
        print(len(path_out))
        self.assertTrue(abs(path_out[-1][0]-2*numpy.pi) < 10**-5)
        self.assertTrue(abs(path_out[-1][1][0]) < 10**-5)

def base_SHO(t, state):
    deltas = numpy.zeros(len(state))
    deltas[0] = state[1]
    deltas[1] = -state[0]
    return deltas

