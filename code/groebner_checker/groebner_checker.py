from sympy import symbols, groebner, reduced

def groebner_independence_check(new_invariant, basis_set, variables):
    """
    Checks if a newly discovered invariant is mathematically independent.
    If it is a syzygy, it can be written entirely as a polynomial of the basis_set.
    """
    if not basis_set:
        return True
        
    # Create abstract variables for the invariants: I_0, I_1, I_2...
    I_vars = symbols(f'I_0:{len(basis_set)}')
    
    # The Ideal is I_k - base_form_k = 0
    ideal_equations = []
    for i, base_form in enumerate(basis_set):
        ideal_equations.append(I_vars[i] - base_form)
        
    try:
        # Map our custom Bracket and Linear objects to standard SymPy symbols (Z_0, Z_1...)
        # This prevents the polynomial engine from getting confused by class attributes.
        custom_objects = set()
        for eq in ideal_equations + [new_invariant]:
            for atom in eq.atoms():
                if type(atom).__name__ in ['Bracket', 'Linear']:
                    custom_objects.add(atom)
            try:
                from noether_engine import Bracket, Linear
                custom_objects.update(eq.find(Bracket))
                custom_objects.update(eq.find(Linear))
            except:
                pass
                
        sub_map = {}
        for i, obj in enumerate(custom_objects):
            sub_map[obj] = symbols(f'Z_{i}')
            
        mapped_ideal = [eq.subs(sub_map) for eq in ideal_equations]
        mapped_new = new_invariant.subs(sub_map)
        
        # Z_vars MUST come first so the algorithm eliminates them
        Z_vars = list(sub_map.values())
        mapped_vars = Z_vars + list(I_vars)
        
        # Compute Gröbner basis
        G = groebner(mapped_ideal, *mapped_vars, order='lex')
        
        # Reduce the new invariant modulo the ideal
        remainder = reduced(mapped_new, G, *mapped_vars)[1]
        
        # Check if the reduction successfully eliminated all geometric Z_vars.
        # If the remainder only contains I_vars, it is a syzygy!
        remainder_vars = remainder.free_symbols
        
        if any(z in remainder_vars for z in Z_vars):
            return True   # It contains geometric variables: it's a NEW generator!
        else:
            return False  # It is a syzygy: purely a polynomial of I_vars
            
    except Exception as e:
        print(f"  [Gröbner Warning] Could not compute basis: {e}")
        return True