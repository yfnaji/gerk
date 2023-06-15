from enum import Enum

class GerkArrayErrorEnum(Enum):
    must_be_list = 0
    bad_item_data_type = 1
    unequal_dim_bc = 2
    condition_fail_b = 3
    condition_fail_bc = 4
    bad_low_tri_dim_A = 5
    bad_square_dim_A = 6
    bad_shape_A = 7
    condition_fail_Ac = 8
    no_real_values = 9
    condition_fail_b_star = 10
    condition_fail_b_star_c  = 11
    solve_not_run = 12


class GerkFunctionErrorEnum(Enum):
    not_callable = 0
    unknown_variables = 1


class GerkAdaptiveErrorEnum(Enum):
    bad_tolerance = 0
    cannot_plot_efficiency_graph = 1


class GerkArrayError(Exception):
    __module__ = "builtins"
    def __init__(self, err, *args):

        if err.value == 0:
            msg = "All arrays must be in the form of a list!"
        elif err.value == 1:
            msg = "Arrays must contain int, float or numeric strings!"
        elif err.value == 2:
            msg = "Arrays `b` and `c` must have the same length\n"
            msg += f"len(b) = {args[0]}, len(c) = {args[1]}"
        elif err.value == 3:
            msg = "condition_b set to True"
            msg += "Values in `b` must sum to 1"
        elif err.value == 4:
            msg = "condition_bc set to True"
            msg += "b * c^T must sum to 1/2"
        elif err.value == 5:
            msg = "Lower triangular matrix must have dimension n-1 "
            msg += "if b and c have dimension n\n"
            msg += f"len(b), len(c) = {args[1]}, len(A) = {args[0]}\n"
            msg += f"len(A) must be {args[0] - 1}"
        elif err.value == 6:
            msg = "Square matrix must have dimension equal to that of b and c\n"
            msg += f"len(b) len(c) = {args[1]}, len(A) = {args[0]}\n"
            msg += f"len(A) must be {args[1]}"
        elif err.value == 7:
            msg = "Matrix A does not confirm to a either square or lower triangular format"
        elif err.value == 8:
            msg = "condition_Ac set to True\n"
            msg += "The sum of row A[i] must be equal to c[i]"
        elif err.value == 9:
            msg = "real_values array has not been set"
        elif err.value == 10:
            msg = "condition_b_star set to True"
            msg += "values in b_star must sum to 1" 
        elif err.value == 11:
            msg = "condition_b_star_c set to True\n"
            msg += "b_star * c^T does not sum to 1/2"
        elif err.value == 12:
            msg = "solve() has not been run\n"
            msg += "y approximation values have not been created yet!"
            

        super().__init__(msg)


class GerkFunctionError(Exception):
    __module__ = "builtins"

    def __init__(self, err, *args):
        if err.value == 0:
            msg = f"Not a valid type for func: {type(*args)}\n"
            msg += "Must be callable (use def or lambda)"
        elif err.value == 1:
            msg = "func must only use arguments `x` and `y`"

        super().__init__(msg)


class GerkAdaptiveError(Exception):
    __module__ = "builtins"

    def __init__(self, err, *args):
        if err.value == 0:
            msg = "tolerance must be greater than 0!"
        elif err.value == 1:
            msg = "Cannot create Efficiency Graph for "
            msg += "Adaptive Runge-Kutta methods"

        super().__init__(msg)