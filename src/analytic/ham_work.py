from expressions import *
import sympy
import sys


def new_ham_equations():
    mathjaxify(new_Ham)
    for qi in new_positions:
        print('#### force')

        sympy.pprint(qi)
        force = -new_Ham.diff(qi)
        sympy.pprint(force)
        print(sympy.latex(force))
        mathjaxify(force)

        # print(sympy.latex(-new_Ham.diff(qi)))
    for pi in new_momenta:
        print('#### velocity')
        sympy.pprint(pi)
        velo = new_Ham.diff(pi)
        sympy.pprint(velo)
        print(sympy.latex(new_Ham.diff(pi)))
        mathjaxify(velo)

        # print(sympy.latex(new_Ham.diff(pi)))

def circular_orbit():
    #dot(phit-2tht)=0
    dot_phit = new_Ham.diff(new_momenta[1])
    dot_tht = new_Ham.diff(new_momenta[2])
    stability1 = dot_phit-2*dot_tht
    sympy.pprint(stability1)
    p_phit_soltn = sympy.solve(stability1, new_momenta[1])[0]
    sympy.pprint(p_phit_soltn)
    #Lt = ptht+2pt
    Lt = new_momenta[2]+2*p_phit_soltn
    sympy.pprint(Lt)
    #dot(pr)=0
    r_force = -new_Ham.diff(new_positions[3])
    r_force = r_force.subs([
        (new_positions[0],0),
        (new_positions[1]-2*new_positions[2],0)
        ])
    sympy.pprint(r_force)
    p_tht_soltn = sympy.solve(r_force, new_momenta[2])[1].simplify()
    sympy.pprint(p_tht_soltn)
    sympy.pprint(Lt.subs(new_momenta[2], p_tht_soltn))
    p_phit_soltn_inr = p_phit_soltn.subs(new_momenta[2], p_tht_soltn)

    circ_energy = new_Ham.subs([
        (new_positions[0],0),
        (new_positions[1],0),
        (new_positions[2],0),
        (new_momenta[0],0),
        (new_momenta[3],0),
        (new_momenta[1],p_phit_soltn_inr),
        (new_momenta[2],p_tht_soltn),
        ])
    r = new_positions[3]
    sympy.pprint(circ_energy)
    sympy.pprint(p_tht_soltn.subs(r, 1.5).evalf())
    sympy.pprint(p_phit_soltn_inr.subs(r, 1.5).evalf())

op_codes=[
    new_ham_equations,
    circular_orbit
]

if __name__ == '__main__':
    if len(sys.argv)==1:
        for i in range(len(op_codes)):
            print(i, op_codes[i].__name__)
    else:
        op_code = int(sys.argv[1])
        print(op_code)
        op_codes[op_code]()
