from sympy import symbols, expand, Add, S, Poly
from noether_engine import Bracket

def evaluate_on_hesse_pencil(bracket_expr, alphabet, lambda_sym):
    def expand_bracket(b):
        a, b_sym, c = b.args
        a1, a2, a3 = symbols(f'{a}1 {a}2 {a}3')
        b1, b2, b3 = symbols(f'{b_sym}1 {b_sym}2 {b_sym}3')
        c1, c2, c3 = symbols(f'{c}1 {c}2 {c}3')
        return (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)

    #Replace all bracket objects with standard polynomial variables
    expr_expanded = bracket_expr
    
    for b in bracket_expr.find(Bracket):
        expr_expanded = expr_expanded.subs(b, expand_bracket(b))
            
    # Multiply it out into a standard polynomial
    expr_expanded = expand(expr_expanded)
    
    # Collect all the variables 
    all_vars = set()
    for dummy in alphabet:
        all_vars.update(symbols(f'{dummy}1 {dummy}2 {dummy}3'))
        
    # Convert to a formal SymPy Poly object for degree checking
    poly_expr = Poly(expr_expanded, *all_vars)
    
    result = S.Zero
    
    for powers, coeff in poly_expr.terms():
        
        term_value = coeff
        valid_term = True
        var_list = list(poly_expr.gens)
        
        for dummy in alphabet:
            v1, v2, v3 = symbols(f'{dummy}1 {dummy}2 {dummy}3')
            
            deg1 = powers[var_list.index(v1)] if v1 in var_list else 0
            deg2 = powers[var_list.index(v2)] if v2 in var_list else 0
            deg3 = powers[var_list.index(v3)] if v3 in var_list else 0
            
            total_deg = deg1 + deg2 + deg3
            if total_deg == 0:
                continue
                
            # Apply Hesse identification rules
            if deg1 == 3 or deg2 == 3 or deg3 == 3:
                term_value *= 1
            elif deg1 == 1 and deg2 == 1 and deg3 == 1:
                term_value *= lambda_sym
            else:
                valid_term = False
                break
                
        if valid_term:
            result += term_value
            
    return result