from collections import deque
from sympy import symbols, S
from noether_engine import (
    Linear, Bracket, canonicalize_expression, 
    apply_plucker_recursive, faltung_ternary
)
from groebner_checker import groebner_independence_check

def is_new_generator(new_form, basis_set, all_vars):
    """
    Two-step verification for new invariants.
    1. Fast symbolic core check (catches scalar multiples instantly).
    2. Rigorous Gröbner syzygy check.
    """
    #Fast Core Check
    new_coeff, new_core = new_form.as_coeff_Mul()
    for existing_form in basis_set:
        exist_coeff, exist_core = existing_form.as_coeff_Mul()
        if new_core == exist_core:
            return False
            
    #Rigorous Gröbner Check
    print("    [Running Gröbner Independence Check...]")
    return groebner_independence_check(new_form, basis_set, all_vars)

def main():
    print("The Noether Code has Started \n")
    
    #Define alphabet and variables
    a, b, c, d = symbols('a b c d')
    alphabet = [a, b, c, d]
    x, u = symbols('x u')
    all_vars = [a, b, c, d, x, u]

    #Initialize the basis set with the foundational quartic form
    f = Linear(a, x)**4
    basis_set = [f]
    
    #Setup the BFS Queue
    queue = deque()
    
    # Fold f with itself
    f_a = Linear(a, x)**4
    f_b = Linear(b, x)**4
    queue.append((f_a, f_b))
    
    print(f"Starting with basis: {basis_set}\n")
    print("Beginning Transvection Loop...\n")

    operations = 0
    max_operations = 4  # Running k=1, 2, 3, 4 for the quartic

    while queue and operations < max_operations:
        form_F, form_G = queue.popleft()
        
        for k in [1, 2, 3, 4]:
            operations += 1
            
            # Transvect
            raw_form = faltung_ternary(form_F, form_G, k, dummy_F=a, dummy_G=b, dual_var=u)
            
            if raw_form == S.Zero:
                continue
                
            # Plücker & Symmetry
            straightened = apply_plucker_recursive(raw_form)
            canonical = canonicalize_expression(straightened, alphabet)
            
            # Syzygy Detection
            if canonical != S.Zero:
                print(f"[Op {operations}] Found Candidate (k={k}): {canonical}")
                if is_new_generator(canonical, basis_set, all_vars):
                    print(f"  -> VALIDATED! Added to Basis Set.\n")
                    basis_set.append(canonical)
                else:
                    print(f"  -> REJECTED! Syzygy detected.\n")

    print("End of Loop")
    print(f"Total independent generators found: {len(basis_set)}")
    for i, gen in enumerate(basis_set):
        print(f"{i+1}. {gen}")

if __name__ == "__main__":
    main()
