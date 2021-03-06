import sympy
#new_positions [phi_d, phi_t, tht, r]
#new_momenta [p_phi_d, p_phi_t, p_tht, p_r]


def mathjaxify(expr):
    print(sympy.latex(expr).replace('\\', '\\\\'))


pphi1, pphi2, ptht, pr = sympy.symbols('p_phi1 p_phi2 p_theta p_r')
vphi1, vphi2, vtht, vr = sympy.symbols('v_phi1 v_phi2 v_theta v_r')
vphid, vphit, vtht, vr = sympy.symbols('v_phid v_phit v_theta v_r')
aphid, aphit, atht, ar = sympy.symbols('a_phid a_phit a_theta a_r')

phi1, phi2, tht, r = sympy.symbols('phi_1 phi_2 theta r')

phid, phit, pd, pt = sympy.symbols('phi_d phi_t p_d p_t', real=True)

positions = [phi1, phi2, tht, r]
velocities = [vphi1, vphi2, vtht, vr]
momenta = [pphi1, pphi2, ptht, pr]

new_positions = [phid, phit, tht, r]
new_velocities = [vphid, vphit, vtht, vr]
new_accs = [aphid, aphit, atht, ar]
new_momenta  = [pd, pt, ptht, pr]
acc_sub = [
    (new_velocities[ind], new_accs[ind])
    for ind in range(len(new_positions))
    ]


Tp = (2*pr**2+2*ptht**2*r**(-2) + 10*pphi1**2 + 10*pphi2**2)/2
Tl = (vr**2/2+vtht**2*r**(-2)/2 + vphi1**2/10 + vphi2**2/10)/2
U = -r**(-3)*(
    sympy.cos(phi1-phi2)
    +3*sympy.cos(phi1+phi2-2*tht)
    )/12
new_U = -r**(-3)*(
    sympy.cos(phid)
    +3*sympy.cos(phit-2*tht)
    )/12
new_Tl = (vr**2/2+vtht**2*r**(-2)/2 + (vphit**2 + vphid**2)/20)/2
new_Tp = (2*pr**2+2*ptht**2*r**(-2) + 20*(pt**2 + pd**2))/2
Ham = Tp+U
Lag = Tl-U
new_Ham = new_Tp + new_U
new_Lag = new_Tl - new_U
tot_Ln = ptht + 2*pt

if __name__ == '__main__':
    sympy.pprint(Tp)
    sympy.pprint(U)
    sympy.pprint(Ham)

