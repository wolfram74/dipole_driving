from expressions import *
from ham_work import circular_orbit
import sympy
from sympy.solvers.solveset import linsolve

import sys
import pdb

def K_mat_gen(H, q_vars, equi_vals):
    # subs = []
    # for index in range(len(q_vars)):
    #     subs.append((
    #         q_vars[index],
    #         equi_vals[index]
    #         ))
    return sympy.Matrix(
        [[
            -H.diff(col_ind, row_ind).subs(equi_vals) for col_ind in q_vars
        ] for row_ind in q_vars]
        )

def M_mat_gen(H, p_vars):
    subs = []
    for pi in p_vars:
        subs.append((pi, 1))
    def Mij(i,j):
        if i==j:
            return 1/H.diff(p_vars[i]).subs(subs)
        else:
            return 0
    return sympy.Matrix(4,4, Mij)

def Hvelocs():
    return [new_Ham.diff(pi) for pi in new_momenta]

def Hforces():
    return [-new_Ham.diff(qi) for qi in new_positions]


def perturbed_quants(terms, order):
    ep = sympy.symbols('epsilon', real=True)
    replacements = []
    new_vars = []
    for term in terms:
        raw = str(term)
        expanded = [ raw+'%d'%expand for expand in range(order+1)]
        symbs = sympy.symbols(' '.join(expanded), real=True)
        total = 0
        for ind in range(order+1):
            total += (symbs[ind]*ep**ind)
        replacements.append((term, total))
        new_vars.append(symbs)
    return replacements, new_vars


def first_round():
    w = sympy.symbols('omega')
    gamma = []
    gamma += new_positions
    gamma += new_momenta
    equi_momen = circular_orbit(True)
    gamma0_subs = [
        (gamma[0],0),
        (gamma[1],0),
        (gamma[2],0),
        (gamma[3],gamma[3]),
        (gamma[4],0),
        (gamma[5],equi_momen[0]),
        (gamma[6],equi_momen[1]),
        (gamma[7],0)
    ]
    sympy.pprint(gamma0_subs)
    K_mat = K_mat_gen(new_Ham, new_positions, gamma0_subs)
    sympy.pprint(K_mat)
    M_mat = M_mat_gen(new_Ham, new_momenta)
    sympy.pprint(M_mat)
    mode_mat = K_mat+w**2*M_mat
    sympy.pprint(mode_mat)

    modes = mode_mat.eigenvects()
    for mode in modes:
        print('new mode')
        sympy.pprint(mode[0])
        freqs = sympy.solve(mode[0],w**2)
        freq_sqr = freqs[0]
        sympy.pprint(freqs)
        sympy.pprint(freq_sqr)
        eig_vecs = [[vi.subs(w**2, freq_sqr)] for vi in mode[2][0]]
        sympy.pprint(eig_vecs)

def print_subs(subs):
    for sub in subs:
        mathjaxify(sub[0])
        mathjaxify(sub[1])

def angle_subs(new_vars, time):
    output=[
        (new_vars[0][0],0),
        (new_vars[1][0],2*time),
        (new_vars[2][0],time),
    ]
    return output
def round_two():
    gamma = []
    gamma += new_positions
    gamma += new_momenta
    ep, Omeg, w, t = sympy.symbols('epsilon Omega omega t', real=True)
    replacements, new_vars = perturbed_quants(gamma, 1)
    zth_angles = angle_subs(new_vars, Omeg*t)
    sympy.pprint(zth_angles)
    # print_subs(replacements)
    # sympy.pprint(replacements)
    # sympy.pprint(new_vars)
    velocs = Hvelocs()
    forces = Hforces()
    # sympy.pprint(velocs)
    expanded_velocs =[]
    expanded_forces =[]
    for veloc in velocs:
        expanded_velocs.append(
            veloc.subs(replacements).series(ep, n=2).removeO()
            )
    for force in forces:
        expanded_forces.append(
            force.subs(replacements).series(ep, n=2).removeO()
            )
    zeroth_ord = []
    first_ord = []
    for term in expanded_velocs+expanded_forces:
        # mathjaxify(term)
        split = sympy.Poly(term.subs(zth_angles), ep).all_coeffs()
        zeroth_ord.append(split[1])
        first_ord.append(split[0])
    sympy.pprint(zeroth_ord[7])
    sympy.pprint((zeroth_ord[7]))
    ptht0 = sympy.solve(zeroth_ord[7],new_vars[6][0])[1]
    ptht0_sub = (new_vars[6][0], ptht0)
    sympy.pprint(ptht0_sub)
    H_first_mat = []
    for ind in range(len(first_ord)):
        row = []
        expr = first_ord[ind]
        expr = expr.subs([
            (Omeg, ptht0),
            ptht0_sub
            ])
        # sympy.pprint((gamma[ind], expr))
        for ind2 in range(len(first_ord)):
            vari = new_vars[ind2][1]
            coeffs = sympy.Poly(expr, vari).all_coeffs()
            if len(coeffs)>1:
                row.append(coeffs[0])
            else:
                row.append(0.)
        # sympy.pprint(row)
        H_first_mat.append(row)
    H_first_mat = sympy.Matrix(H_first_mat)
    sympy.pprint(H_first_mat)
    # charpoly = H_first_mat.charpoly(w).simplify().as_expr()
    # t1 = H_first_mat.charpoly(w).simplify()
    # sympy.pprint(charpoly)
    # sympy.pprint(sympy.collect(charpoly,w**2))
    # sympy.pprint(charpoly.simplify().collect(w))
    # sympy.pprint(sympy.Poly(charpoly, w).all_coeffs())
    # test = w**2-1
    # soltns = sympy.solve(test, w)
    # sympy.pprint(soltns)
    # soltns = sympy.solvers.solve(charpoly.simplify(), w**2)
    # sympy.pprint(soltns)
    # sympy.pprint(
    #     H_first_mat.eigenvals()
    #     )
# 13.04
    # system = H_first_mat.col_insert(-1, sympy.Matrix([0,0,0,0,0,0,0,0]))
    # linsolve(system)
    # soltn = linsolve(system, gamma)
    # soltn = sympy.solve(system)
    # soltn = sympy.solvers.solve_linear_system(system)
    # sympy.pprint(soltn)
    # pdb.set_trace()
    modes = H_first_mat.eigenvects(error_when_incomplete=False)
    for mode in modes:
        print('new mode')
        sympy.pprint(mode[0])
        freqs = mode[0]
        freq_sqr = (freqs**2).expand().simplify()
        # sympy.pprint(freqs)
        # sympy.pprint(mode[2])
        sympy.pprint(freq_sqr)
        # pdb.set_trace()
        # eig_vecs = [[vi.subs(freq_sqr,w**2).expand().simplify()] for vi in mode[2][0]]
        # sympy.pprint(eig_vecs)
def round_three():
    pairs = []
    ep, Omeg, w, t = sympy.symbols('epsilon Omega omega t', real=True)

    for ind in range(len(new_positions)):
        force = new_Lag.diff(new_positions[ind])
        momen = new_Lag.diff(new_velocities[ind]).subs(acc_sub)
        if ind ==2:
            momen = momen.subs([(atht, (atht-2*vtht*vr/r))]).expand()
        pairs.append([force, momen])
        print(ind)
        sympy.pprint(pairs[-1])
    lag_eqs = [pair[0]-pair[1] for pair in pairs]
    full_coords = new_positions+new_velocities+new_accs
    pert_subs, new_coords = perturbed_quants(full_coords, 1)
    zth_angles = angle_subs(new_coords, Omeg*t)
    sympy.pprint(lag_eqs)
    sympy.pprint(pert_subs)

    subbed_eqs = [eqn.subs(pert_subs) for eqn in lag_eqs]
    zeroth_ord = []
    first_ord = []
    for eqtn in subbed_eqs:
        # mathjaxify(term)
        eqtn = eqtn.subs(zth_angles)
        sympy.pprint(eqtn)
        taylored = eqtn.series(ep, n=2).removeO()
        sympy.pprint(taylored)
        # split = sympy.Poly(term.subs(zth_angles), ep).all_coeffs()
        # zeroth_ord.append(split[1])
        # first_ord.append(split[0])



op_codes=[
    first_round,
    round_two,
    round_three,
]

if __name__ == '__main__':
    if len(sys.argv)==1:
        for i in range(len(op_codes)):
            print(i, op_codes[i].__name__)
    else:
        op_code = int(sys.argv[1])
        print(op_code)
        op_codes[op_code]()
