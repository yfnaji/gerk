import numpy as np
import matplotlib.pyplot as plt
import math
import decimal
from _gerk_error import (
    GerkArrayError,
    GerkFunctionError, 
    GerkAdaptiveError,
    GerkArrayErrorEnum,
    GerkFunctionErrorEnum,
    GerkAdaptiveErrorEnum,
)

# decimal.getcontext().prec = 16
decimal.getcontext().rounding = decimal.ROUND_HALF_EVEN


class Gerk:
    """
    The Generalized Explicit Runge-Kutta class (GERK)

    Attributes
        A (list): A matrix of the Butcher tableau
        b (list): b array of the Butcher tableau
        c (list): c array of the Butcher tableau
        time_steps (int or float): Number of discretizations for Ordinary Runge-Kutta method
                                   Initial step size for the adaptive Runge-Kutta method
        initial_conditions (tuple): initial condition coordinate
        final (float or int): Final value before the Runge-Kutta method terminates
        func (callable): The function defined in the system dy/dx = f(x,y)
        real_values (callable or list): The function that solves the system explicitly 
                                        or list containing the value of the real y values
        b_star (list): The array of the additional b values for the adaptive Runge-Kutta method
        tolerance (float): The minimum error threshold for the adaptive Runge-Kutta method

    Additional Paramters
        condition_b (bool): Activate condition sum(b) == 1
        condition_bc (bool): Activate condition sum(b*c) == 1/2
        condition_Ac (bool): Activate condition sum(A[i]) == c[i]
        condition_b_star (bool): Activate condition sum(b_star) == 1
        condition_b_star_c (bool): Activate condition sum(b_star*c) == 1/2
    """
    def __init__(
        self, 
        A, 
        b, 
        c, 
        time_steps, 
        initial_conditions, 
        final, 
        func,
        real_values = None,
        b_star = None,
        tolerance = None,
        condition_b = False,
        condition_bc = False,
        condition_Ac = False,
        condition_b_star = False,
        condition_b_star_c = False,
    ):
    
        def _process_array(arr):
    
            if type(arr) is not list:
                raise GerkArrayError(GerkArrayErrorEnum.must_be_list)
            
            try:
                if type(arr[0]) is not list:
                    arr = [arr]

                for row in range(len(arr)):
                    for i in range(len(arr[row])):
                        x = str(arr[row][i])
                        if "/" in x:
                            num, den = x.split("/")
                            arr[row][i] = decimal.Decimal(num)/decimal.Decimal(den)
                        else:
                            arr[row][i] = decimal.Decimal(x)
            
            except (NameError, TypeError):
                raise GerkArrayError(GerkArrayErrorEnum.bad_item_data_type)
            
            if len(arr) == 1:
                return np.array(arr[0])
            else:
                return np.array(arr)

        def _is_lower_triangular(mat):
            for i in range(len(mat)):
                if len(mat[i]) != i+1:
                    return False
            return True
    
        def _is_square(mat):
            n = len(mat)
            for i in range(len(mat)):
                if len(A[i]) != n:
                    return False
            return True

        n = len(b)
        if n != len(c):
            raise GerkArrayError(GerkArrayErrorEnum.unequal_dim_bc, len(b), len(c))
        if condition_b and sum(b) != 1:
            raise GerkArrayError(GerkArrayErrorEnum.condition_fail_b)
        if condition_bc and sum([_b*_c for _b, _c in zip(b,c)]) != 1/2:
            raise GerkArrayError(GerkArrayErrorEnum.condition_fail_bc)
        if b_star:
            if condition_b_star and sum(b_star) != 1:
                raise GerkArrayError(GerkArrayErrorEnum.condition_fail_b_star)
            if condition_b_star_c and sum([_b_star*_c for _b_star, _c in zip(b_star,c)]) != 1/2:
                raise GerkArrayError(GerkArrayErrorEnum.condition_fail_b_star_c)

        if _is_lower_triangular(A):
            if len(A) + 1 != n:
                raise GerkArrayError(GerkArrayErrorEnum.bad_low_tri_dim_A)
            A = [[]] + A
        elif _is_square(A):
            if len(A) != n:
                raise GerkArrayError(GerkArrayErrorEnum.bad_square_dim_A, len(A), len(b))                
        else:
            raise GerkArrayError(GerkArrayErrorEnum.bad_shape_A)
        
        self.A = _process_array(A)
        self.c = _process_array(c)

        if condition_Ac:
            for i,row in enumerate(self.A):
                if sum(row) != self.c[i]:
                    raise GerkArrayError(GerkArrayErrorEnum.condition_fail_Ac)

        if not callable(func):
            raise GerkFunctionError(GerkFunctionErrorEnum.not_callable)
        
        self.func = func
        
        if real_values:
            if not callable(real_values):
                raise GerkFunctionError(GerkFunctionErrorEnum.not_callable)
            else:
                self.real_values = real_values
            
        self.b = _process_array(b)
        self.b_star = _process_array(b_star) if b_star else None

        if b_star is not None and not tolerance:
            self.tolerance = decimal.Decimal(0.001)
        elif tolerance and tolerance < 0:
                raise GerkAdaptiveError(GerkAdaptiveErrorEnum.bad_tolerance)
        elif tolerance:
            self.tolerance = decimal.Decimal(tolerance)

        self.initial_condition = initial_conditions
        self.final = final
        self.time_steps = time_steps

    def solve(self):
        """
        Executes the Runge-Kutta method on the system defined in the object 
        """
        k = [decimal.Decimal(0) for i in range(len(self.b))]
        if self.b_star is None:
            h = decimal.Decimal((self.final-self.initial_condition[0]))/decimal.Decimal((self.time_steps))
        else:
            h =decimal.Decimal(self.time_steps)
        
        if self.b_star is not None:
            order = 1/decimal.Decimal(min(len([_b for _b in self.b if _b]), len([_b for _b in self.b_star if _b])))
        
        x_n = decimal.Decimal(self.initial_condition[0])
        y_n = decimal.Decimal(self.initial_condition[1])
        self.X = [x_n]
        self.Y = [y_n]
        
        if x_n < self.final:
            _forward = True
            _backward = False
        else:
            _forward = False
            _backward = True

        def _k_eval(x, y, k_values):
            for i in range(len(self.A)):
                x += self.c[i]*h
                y += h*sum([self.A[i][j]*k_values[j] for j in range(len(self.A[i]))])
                
                try:
                    k_values[i] = self.func(x,y)
                except TypeError:
                    raise GerkFunctionError(GerkFunctionErrorEnum.unknown_variables)

        while (_forward and x_n < self.final) or (_backward and x_n > self.final):

            _k_eval(x_n, y_n, k)
            if self.b_star is not None:
                
                y_auth = h*sum(self.b*k)
                y_hat = h*sum(self.b_star*k)
                
                e = abs(y_auth-y_hat)

                if e > self.tolerance:
                    h *= decimal.Decimal(0.9)*(self.tolerance/e)**order
                    continue

            x_n += h
            y_n += h*sum(self.b*k)
            self.X.append(x_n)
            self.Y.append(y_n)


    def plot(
            self, 
            with_real=False,
            title="Runge-Kutta Approximation",
            x_label="x",
            y_label="y"
        ):
        """
        Plots the approximated curve created by the Runge-Kutta method
        Note: Must run solve() first

        Parameters
            with_real (bool): Plot real curve with approximated curve
                              Note: Requires real_values attribute to be set
            title (str): Title of the plot
            x_label (str): x-axis label
            y_label (str): y-axis label
        """

        
        if with_real:
            _y = [self.real_values(x) for x in self.X]
            plt.plot(self.X, _y, color="b", label="Exact Soln")
        
        plt.plot(self.X, self.Y, color="r", label="RK")
        
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        plt.show()


    @property
    def get_approximations(self):
        """
        Returns a list of the approximated values at every time step

        Returns
            list
        """
        if not self.Y:
            raise GerkArrayError(GerkArrayErrorEnum.no_real_values)
        return self.Y

    @property
    def get_errors(self):
        """
        Returns a list of the errors at every time step

        Returns
            list
        """
        if not self.real_values:
            raise GerkArrayError(GerkArrayErrorEnum.no_real_values)
        
        if callable(self.real_values):
            errors = []
            for i,x in enumerate(self.X):
                errors.append(abs(self.Y[i]-decimal.Decimal(self.real_values(x))))
            return errors
        else:
            if not self.Y:
                raise GerkArrayError(GerkArrayErrorEnum.solve_not_run)
            return abs(self.real_values-self.Y)


    def efficiency_graph(self):
        """
        Creates the efficiency graph for the system defined in the object
        """

        if self.b_star is not None:
            raise GerkAdaptiveError(GerkAdaptiveErrorEnum.cannot_plot_efficiency_graph)

        discretizations = [
            100,
            1_000,
            10_000,
            100_000,
            1_000_000
        ]        

        max_e = []
        for ts in discretizations:
            self.time_steps = ts
            self.solve()
            max_e.append(math.log10(max(self.get_errors)))
 
        t_ = [math.log10(t) for t in discretizations]
        plt.title("Efficiency Graph")
        plt.scatter(t_, max_e)
        plt.xlabel("$log_{10}(discretizations)$")
        plt.ylabel("$log_{10}(\max(error))$")
        plt.plot(t_, max_e)
        plt.show()

A = [ 
        ["1/3"], # Remember, fractions must be passed as strings
        ["-1/3", 1],
        [1, -1, 1]
]

b = ["1/8", "3/8", "3/8", "1/8"]

c = [0, "1/3", "2/3", 1]

rk_obj = Gerk(
        A=A, 
        b=b, 
        c=c, 
        initial_conditions=(0, 1), 
        time_steps=1000, 
        final=5,
        func=lambda x, y: y, # Remember, you must define both x and y as variables, even if one is not used
        real_values = lambda y: math.exp(y)
    )
rk_obj.solve()
rk_obj.plot(with_real=True)