import sys
import sympy
from expressions import *

''' includes
positions, velocities, momenta
Tp,Tl,U,
Ham, Lag
'''

def ham_equations():
    for qi in positions:
        sympy.pprint(qi)
        sympy.pprint(-Ham.diff(qi))
        # print(sympy.latex(-Ham.diff(qi).subs(positions[3], 1)))
        print(sympy.latex(-Ham.diff(qi)))
    for pi in momenta:
        sympy.pprint(pi)
        sympy.pprint(Ham.diff(pi))
        # print(sympy.latex(Ham.diff(pi).subs(positions[3], 1)))
        print(sympy.latex(Ham.diff(pi)))


def equi_from_jk(j, k):
    return [
        (j+k)*sympy.pi/2,
        (j-k)*sympy.pi/2,
        0
    ]

def K_mat_gen(H, q_vars, equi_vals):
    subs = []
    for index in range(len(q_vars)):
        subs.append((
            q_vars[index],
            equi_vals[index]
            ))
    return sympy.Matrix(
        [[
            -H.diff(col_ind, row_ind).subs(subs) for col_ind in q_vars
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
    return sympy.Matrix(3,3, Mij)


def mode_analysis():
    Gam00 = equi_from_jk(0,1)
    sympy.pprint(Ham)
    sympy.pprint(U)
    angles = positions[:3]
    p_angles = momenta[:3]
    w = sympy.symbols('omega')
    contact_ham = Ham.subs(positions[-1], 1)
    sympy.pprint(contact_ham.diff(positions[0]))
    sympy.pprint(contact_ham.diff(positions[0], positions[0]))
    k_mat = K_mat_gen(contact_ham, angles , Gam00)
    print('inspecting modes at')
    sympy.pprint(Gam00)
    sympy.pprint(k_mat)
    M_mat = M_mat_gen(contact_ham, p_angles)
    sympy.pprint(M_mat)
    mode_mat = k_mat+w**2*M_mat
    sympy.pprint(mode_mat)
    print( sympy.latex(mode_mat).replace('\\\\', '\\\\\n'))

    modes = mode_mat.eigenvects()
    for mode in modes:
        sympy.pprint(mode[0])
        freqs = sympy.solve(mode[0],w**2)
        freq_sqr = freqs[0]
        sympy.pprint(freqs)
        sympy.pprint(freq_sqr)
        eig_vecs = [[vi.subs(w**2, freq_sqr)] for vi in mode[2][0]]
        sympy.pprint(eig_vecs)

def new_ham_equations():
    for qi in new_positions:
        sympy.pprint(qi)
        sympy.pprint(-new_Ham.diff(qi))
        print(sympy.latex(-new_Ham.diff(qi).subs(positions[3], 1)))
        # print(sympy.latex(-new_Ham.diff(qi)))
    for pi in new_momenta:
        sympy.pprint(pi)
        sympy.pprint(new_Ham.diff(pi))
        print(sympy.latex(new_Ham.diff(pi).subs(positions[3], 1)))
        # print(sympy.latex(new_Ham.diff(pi)))

op_codes=[
    ham_equations, mode_analysis,
    new_ham_equations
]

if __name__ == '__main__':
    if len(sys.argv)==1:
        for i in range(len(op_codes)):
            print(i, op_codes[i].__name__)
    else:
        op_code = int(sys.argv[1])
        print(op_code)
        op_codes[op_code]()
