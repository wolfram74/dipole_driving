from expressions import *
import sympy
import sys


def new_ham_equations():
    for qi in new_positions:
        print('#### force')

        sympy.pprint(qi)
        sympy.pprint(-new_Ham.diff(qi))
        print(sympy.latex(-new_Ham.diff(qi).subs(positions[3], 1)))
        # print(sympy.latex(-new_Ham.diff(qi)))
    for pi in new_momenta:
        print('#### velocity')
        sympy.pprint(pi)
        sympy.pprint(new_Ham.diff(pi))
        print(sympy.latex(new_Ham.diff(pi).subs(positions[3], 1)))
        # print(sympy.latex(new_Ham.diff(pi)))

op_codes=[
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
