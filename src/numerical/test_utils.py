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

class BouncingTests(unittest.TestCase):
    def test_spin_stays_at_r1(self):
        stateI = numpy.array([0.,0.,0.,1.,.5,0.,0.,0.,])
        # print('spinning #####')
        path_out = utils.RK45_bouncing_path(
            utils.reduced_dipole_equations, stateI, 0., 20,
            max_steps=500, precision=10**-7)
        # print(path_out[0])
        # print(path_out[1])
        self.assertTrue(path_out[-1][1][3]==1.)

    def test_orbital_stays_at_r1(self):
        stateI = numpy.array([0.,0.,0.,1.,0.,.1,-.05,0.,])
        path_out = utils.RK45_bouncing_path(
            utils.reduced_dipole_equations, stateI, 0., 20,
            max_steps=1000, precision=10**-7)
        self.assertTrue(path_out[-1][1][3]==1.)

    @unittest.skip('guessing bounce time is hard')
    def test_bounce_stays_above_r1(self):
        stateI = numpy.array([0.,0.,0.,1.0,0.,0.,0.,0.2,])
        path_out = utils.RK45_bouncing_path(
            utils.reduced_dipole_equations, stateI, 0., 20,
            max_steps=500, precision=10**-7)
        any_inwards = False
        # for state in path_out[:5]:
        #     print state
        for ind in range(len(path_out)):
            state = path_out[ind]
            if state[1][3]<1.:
                print('failed state')
                pre = path_out[ind-1]
                print(pre)
                print(state)
                print(state[0]-pre[0])

            self.assertTrue(state[1][3]>=1.)
            any_inwards = any_inwards or state[1][7]<0.
        self.assertTrue(any_inwards)

    def test_conservations(self):
        stateI = numpy.array([
            0.,0.,0.,1.0,
            0.1, -0.1 , 0.25,0.2
            ])
        path_out = utils.RK45_bouncing_path(
            utils.reduced_dipole_equations, stateI, 0., 20,
            max_steps=1000, precision=10**-7)
        E_init = utils.total_energy(path_out[0][1])
        L_init = utils.total_L(path_out[0][1])
        E_fin = utils.total_energy(path_out[-1][1])
        L_fin = utils.total_L(path_out[-1][1])
        del_E = abs(E_fin-E_init)
        del_L = abs(L_fin-L_init)
        f_state = path_out[-1][1]
        print(del_E, E_init, E_fin)
        print(path_out[-1][1])
        print(utils.PE(f_state))
        print(utils.KE(f_state))
        # print(del_L)
        self.assertTrue(del_E < 10**-5)
        self.assertTrue(del_L < 10**-5)

def base_SHO(t, state):
    deltas = numpy.zeros(len(state))
    deltas[0] = state[1]
    deltas[1] = -state[0]
    return deltas

