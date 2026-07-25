import itertools
from sympy import Expr, Symbol, S, sympify, Add, Mul, default_sort_key

def permutation_parity(lst):
    inversions = 0
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if str(lst[i]) > str(lst[j]):
                inversions += 1
    return -1 if inversions % 2 != 0 else 1

class Linear(Expr):
    """Represents the linear product a_x"""
    is_commutative = True
    def __new__(cls, symbol, variable):
        return Expr.__new__(cls, sympify(symbol), sympify(variable))
    
    def _sympystr(self, printer):
        return f"{self.args[0]}_{self.args[1]}"

class Bracket(Expr):
    """Represents the determinant bracket (a b c)"""
    is_commutative = True
    def __new__(cls, *args):
        args = [sympify(arg) for arg in args]
        
        # Repeated symbols in a determinant evaluate to 0
        if len(set(args)) < len(args): 
            return S.Zero
            
        # Sort alphabetically and calculate parity sign
        sorted_args = sorted(args, key=str)
        sign = permutation_parity(args)
        
        obj = Expr.__new__(cls, *sorted_args)
        return obj if sign == 1 else -obj
        
    def _sympystr(self, printer):
        return f"({' '.join(printer._print(arg) for arg in self.args)})"

def canonicalize_monomial(term, alphabet):
    coeff, core_term = term.as_coeff_Mul()
    if core_term == S.One: 
        return term
        
    used_dummies = [sym for sym in core_term.free_symbols if sym in alphabet]
    k = len(used_dummies)
    if k == 0: 
        return term
        
    target_dummies = alphabet[:k]
    mapped_forms = set()
    
    for p in itertools.permutations(target_dummies):
        sub_dict = dict(zip(used_dummies, p))
        new_term = core_term.subs(sub_dict, simultaneous=True)
        mapped_forms.add(new_term)
        
    # If swapping dummy variables negates the form, it is 0
    for form in mapped_forms:
        if -form in mapped_forms: 
            return S.Zero
            
    canonical_core = sorted(list(mapped_forms), key=default_sort_key)[0]
    return coeff * canonical_core

def canonicalize_expression(expr, alphabet):
    expr = expr.expand()
    terms = Add.make_args(expr)
    canonical_terms = [canonicalize_monomial(t, alphabet) for t in terms]
    return Add(*canonical_terms)

def apply_plucker_to_monomial(term):
    coeff, core = term.as_coeff_Mul()
    if core == S.One or not isinstance(core, Mul): 
        return term
        
    factors = core.args
    brackets = [f for f in factors if isinstance(f, Bracket)]
    linears = [f for f in factors if isinstance(f, Linear)]
    
    if not brackets or not linears: 
        return term
    
    for bracket in brackets:
        B_args = list(bracket.args)
        max_B_sym = B_args[2]  # The largest symbol in the sorted bracket
        
        for linear in linears:
            l_sym = linear.args[0]
            var = linear.args[1]
            
            # If a linear symbol is smaller than the largest bracket symbol, rewrite!
            if str(l_sym) < str(max_B_sym):
                remaining_factors = list(factors)
                remaining_factors.remove(bracket)
                remaining_factors.remove(linear)
                remainder = Mul(*remaining_factors)
                
                p, q, r = B_args[0], B_args[1], B_args[2]
                s = l_sym
                
                # Identity: (p q r)s_x = (q r s)p_x - (p r s)q_x + (p q s)r_x
                term1 =  Bracket(q, r, s) * Linear(p, var)
                term2 = -Bracket(p, r, s) * Linear(q, var)
                term3 =  Bracket(p, q, s) * Linear(r, var)
                
                new_expression = coeff * remainder * (term1 + term2 + term3)
                
                # Recursively apply in case new terms need straightening
                return apply_plucker_recursive(new_expression)
                
    return term

def apply_plucker_recursive(expr):
    expr = expr.expand()
    if isinstance(expr, Add):
        return Add(*[apply_plucker_to_monomial(t) for t in expr.args])
    else:
        return apply_plucker_to_monomial(expr)


def faltung_ternary(form_F, form_G, k, dummy_F, dummy_G, dual_var):
    """
    Computes the k-th ternary transvection of form_F and form_G.
    Replaces k instances of (dummy_F)_x and (dummy_G)_x with the bracket.
    """
    def extract_linears(term, symbol):
        # 1. If it's a single item (not multiplied by anything)
        if not isinstance(term, Mul):
            if isinstance(term, Linear) and term.args[0] == symbol:
                return 1, S.One
            if term.is_Pow and isinstance(term.base, Linear) and term.base.args[0] == symbol:
                return term.exp, S.One
            return 0, term
            
        linears_count = 0
        other_factors = []
        for factor in term.args:
            if isinstance(factor, Linear) and factor.args[0] == symbol:
                linears_count += 1
            elif factor.is_Pow and isinstance(factor.base, Linear) and factor.base.args[0] == symbol:
                linears_count += factor.exp
            else:
                other_factors.append(factor)
                
        return linears_count, Mul(*other_factors)

    count_F, remainder_F = extract_linears(form_F, dummy_F)
    count_G, remainder_G = extract_linears(form_G, dummy_G)

    # Cannot transvect if we don't have enough linear factors
    if count_F < k or count_G < k: 
        return S.Zero

    new_bracket = Bracket(dummy_F, dummy_G, dual_var)
    rem_linears_F = Linear(dummy_F, 'x')**(count_F - k) if count_F - k > 0 else S.One
    rem_linears_G = Linear(dummy_G, 'x')**(count_G - k) if count_G - k > 0 else S.One

    return remainder_F * remainder_G * rem_linears_F * rem_linears_G * (new_bracket**k)