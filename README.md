# Generalized Explicit Runge-Kutta (GERK)

A package for the curious mathematicians and engineers who want to experiment the Runge-Kutta method with their own coefficients.

# Contents

1. [Runge-Kutta Overview](#rko)
2. [What is Gerk?](#wig)
3. [How to use Gerk](#h2ug)
    * [Attributes](#attr)
    * [Condition Arguments](#cond-arg)
    * [Methods](#methods)
3. [Example](#example)
    * [Creating the Gerk Object](#cgo)
    * [plot example](#plot-eg)
    * [get_approximations example](#get-approx-eg)
    * [get_errors example](#get-err-eg)
    * [efficiency_graph example](#eff-graph-eg)
5. [Adaptive Runge-Kutta Methods](#ark)
6. [Adaptive Runge-Kutta Example](#ark-eg)

<h1 id="rko">Runge-Kutta Overview</h1>

Runge-Kutta methods aim to numerically solve ordinary differential equations of the form:

$$
\frac{\text{d}y}{\text{d}x}=f(x,y)
$$

with a given initial condition $(x_0, y_0)$. 

We define a member of the Runge-Kutta family with a *butcher tableau*:

$$
  \begin{array}{c| c c c c}
    0\\
    c_2 & a_{21}\\
    c_3 & a_{31} & a_{32}\\
    \vdots & \vdots &  & \ddots\\
    c_r & a_{r1} & a_{r2} & \cdots & a_{rr-1}\\
    \hline
      & b_1 & b_2 & \cdots & b_{r-1} & b_r
  \end{array}\\ \\
$$


The above tableau is often abbreviated to

$$
  \begin{array}{c| c}
    c & A\\
    \hline
      & b
  \end{array}\\ \\
$$

The Butcher tableau presents the category of coefficients that will be used to define our Runge-Kutta method.

The $n^\text{th}$ evaluation of the solution will be denoted as $(x_n, y_n)$. We also define $h$ as the time step i.e. the step size from the previous approximation to the next, and therefore

$$
x_{n+1} = x_{n} + h
$$

Before defining $y_n$, we must familiarize ourselves with another array $k$. We define the $i^{\text{th}}$ row of $k$ for the $n^\text{th}$ as:

$$
k_i(x_n, y_n) = f\left(x_n + c_i h, y_n + \sum^{i-1}_{j=1} a_{ij}k_j(x_n, y_n)\right)
$$

where $f$ is the function defined in the differential equation above. Note the recursion in the second argument of $f$ where we sum rows preceding $k_{i}$ and apply a scale factor of $a_{ij}$.

We now have everything we need to calculate $y_n$:

$$
y_{n+1} = y_{n} + h \sum^{r}_{i=1}b_i k_{i}(x_n, y_n)
$$

The idea is that we try to approximate the slope at point $y_n$ to ascertain the ascent (or descent) to attain the next predicted value of $y$.

<h1 id="wig">What is Gerk?</h1>

Most packages for the Runge-Kutta method usually have the coefficients $a_{ij}$, $b_i$ and $c_i$ determined beforehand for known methods such of the *Forward-Euler method*, the *1/4 rule*, *the 3/8 rule* etc, but do not allow one to create their own Runge-Kutta.

Gerk is an easy interface to allow the user to determine their own coefficient values for the Runge-Kutta method. The package can return the approximations in the form of an array, a plot (using matplotlib) and also compare with real values (if available) to produce error values and even an efficiency graph!

<h1 id="h2ug">How to use Gerk</h1>

We can simply import `Gerk` in the following way:

```
from gerk import Gerk
```

provided you have placed the file in the appropriate folder.

<h2 id="attr">Attributes</h2>

As seen in mathematics above, there is quite a bit information required to execute the Runge-Kutta method. This has been broken down into arguments to be passed into the `Gerk` class:

- `A` The $A$ matrix in the Butcher tableau. This **must** be a lower triangular matrix that is formatted as a list of lists that contain floats, decimals or integers **Note**: Fractions must be in a string i.e. "1/3"
- `b` The $b$ array. Must be a list of floats, decimals or integers
- `c` The $c$ array. Must be a list of floats, decimals or integers
- `initial_conditions` A `tuple` that acts as the coordinate of the initial condition values
- `final` The value of $x$ for which we terminate the Runge-Kutta method
- `time_steps` The number of discretizations you want to apply on the from the initial $x_0$ to `final`. Must be an integer
- `func` The function expression to be numerically integrated
- `real_values` (Optional) A `callable` which is the explicit solution to the initial value problem _**or**_ a list of real values **Note** the list must be of size `time_steps + 1`


<h2 id="cond-arg">Condition Arguments</h2>

There is no consensus to what conditions must hold regarding the coefficients you choose for your method, however, some known Runge-Kutta methods do consistently conform to some of these conditions. 

In this package, we have given you the option to overlook these conditions to provide more flexibility with experimentation.

$$\sum^{r}_{i=1}b_i=1 \ \ \ \ \sum^{r}_{i=1}b_ic_i = 1/2 \ \ \ \ \sum^{r}_{j=1}a_{ij} = c_i$$

`condition_b`, `condition_bc` and `condition_Ac`

These are all boolean arguments that are set to `False` by default to allow the user to freely experiment.

<h2 id="methods">Methods</h2>

* `solve` calculates the value of $y_n$ for each time step. It will return all the approximated values in a list

Once `solve` has been run, the following methods will be available to run:

* `plot` Creates a graph of the approximated curve. If you have set the (optional) `real_values` argument, you can plot both the approximated and exact curves in one plot by setting `with_real=True`. You also have the option to amend `x_label` and `y_label` which are arugments of the method
* `get_approximation` returns the approximated values of $y$
* `get_errors` returns the errors (only applicable if `real_values` argument was used) 

<h1 id="example">Example</h1>

<h2 id="cgo">Creating the Gerk Object</h2>

In the example below, we will use the 3/8-th rule to approximate the following initial value problem:

$$
\frac{\text{d}y}{\text{d}x} = y
$$

with the initial condition $(0, 1)$ and will be making approximations up to $x=5$. The domain we are working on $0\leq x\leq5$ will be discretized into 1000 steps.

The 3/8-th rule has the following Butcher tableau:

$$
  \begin{array}{c| c c c}
    0\\
    1/3 & 1/3\\
    2/3 & -1/3 & 1\\
    1 & 1 & -1 & 1\\
    \hline
      & 1/8 & 3/8 & 3/8 & 1/8
  \end{array}\\ \\
$$


The $A$ lower triangular matrix in the Butcher tableau above can be implemented in the following way:

```
A = [ 
        ["1/3"], # Remember, fractions must be passed as strings
        ["-1/3", 1],
        [1, -1, 1]
]
```

Note that if you would rather use "square" matrices, you can also define `A` as:

```
A = [
        [0,0,0,0], 
        ["1/3", 0, 0, 0],
        ["-1/3", 1, 0, 0],
        [1, -1, 1, 0]
]
```

Either version of `A` is acceptable in `Gerk`.

Now for the `b` and `c` vectors:

```
b = ["1/8", "3/8", "3/8", "1/8"]

c = [0, "1/3", "2/3", 1]
```

We now have sufficient information to utilise the `Gerk` class. However, we know the solution to the initial value problem is $y=e^x$. So in this case we can utilise the `real_values` parameter by setting it to the lambda function:

```
lambda x,y: math.exp(x)
```


Now we are ready to create the `Gerk` object:

```
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
```

To finally get the approximated values, we run `rk_obj.solve()`.

After running `solve`, 3 methods are now available to us:

<h2 id="plot-eg">`plot` example</h2>

We can now plot our approximated curve. Note that because we have defined `real_values`, we also have the option to plot the exact solution by setting the `with_real` flag to `True` (which is the default).

```
rk_obj.plot(with_real=True)
```

![rk_example](https://user-images.githubusercontent.com/59436765/210466156-f86bfb1e-f088-4791-a16a-01220298e6bd.png)


<h2 id="get-approx-eg">`get_approximations` example</h2>

The approximated values $y_n$ for each timestep can be obtained by running `get_approximations` which returns a list of the values:

```
rk_obj.get_approximations[2:5]
```

> [Decimal('1.010062750717580079753282335'), Decimal('1.015132034738508612169525261'), Decimal('1.020226760387164317916560525')]

_Note_ that this is a _property_ and not a method, so no need for parentheses.

<h2 id="get-err-eg">'get_errors' example</h2>

The errors for each time step can also be obtained if the `real_value` parameter has been set.

In our example, we set `real_values=lambda x: math.exp(x)` which will be used to calculate the error at each step. 

By simply calling `get_errors`, we obtain all the errors in the form of a list (again, no need for parentheses):

```
rk_obj.get_errors[2:5]
```

> [Decimal('0.00001258363341213062577786389261'), Decimal('0.00001897012278968249653135224731'), Decimal('0.00002542036040854158088425155254')]

<h2 id="eff-graph-eg">Efficiency Graphs</h2>

We can create an efficiency graph for the Runge-Kutta method provided that `real_values` has been set to a callable function. This can be done by simply running:

```
rk_obj.efficiency_graph()
```

![efficiency_graph](https://user-images.githubusercontent.com/59436765/210466181-fbc3a2ea-0fac-44b9-a067-bf01bf72ee8f.png)


_Note:_ You do not need to run `solve()` beforehand before producing the efficiency graph.

<h1 id="ark">Adaptive Runge-Kutta Methods</h1>

There is an alternate way to utilise the Runge-Kutta method by employing an additional distinct $b$ array. The Butcher tableau for such methods take the form:


<img width="242" alt="adaptive_butcher" src="https://user-images.githubusercontent.com/59436765/210662922-f5fbf612-a56b-4679-af2c-5eca4956771d.png">


where $b^*_i$ is the additional $b$ array.

This method is not too disimilar to the original Runge-Kutta method, but there is an extra step. We calculate the $k$'s in the same way as before, but there is an extra step when evaluating $y_{n+1}$. 

Although we do calculate $y_{n+1}$ in same way outlined above, we also calculate $\hat{y}_{n+1}$:

$$
y_{n+1} = y_{n} + h \sum^{r}_{i=1}b_i\cdot k_{i}(x_n, y_n) \ \ \ \ \ \ \ \ \ \hat{y}_{n+1} = y_{n} + h \sum^{r}_{i=1}b^*_i\cdot k_{i}(x_n, y_n)
$$

_Note_ that the calculation for $\hat{y}$ requires the use of $y_n$ and **not** $\hat{y}_n$.

For every time step, we calculate the following error:

$$
E := \left|y_{n+1}-\hat{y}_{n+1}\right|
$$

Now we define the tolerance, $\mathcal{E}$, which will act as the upper bound for $E$.

If $E<\mathcal{E}$, then we accept the value of $y_{n+1}$ and we increment $x_{n}$ by $h$ and start the process again for $x_{n+1}$. 

However, if $E\geq\mathcal{E}$, we will need to adjust the value of $h$ and redo the calculation with this renewed value in an attempt to statisfy the condition $E<\mathcal{E}$. 

The value of $h$ will be adjusted as follows:

$$
h \rightarrow 0.9\cdot h \cdot\sqrt[n]{\frac{h}{\mathcal{E}}}
$$

where $n=\min\left(p,q\right)$, where $p$ and $q$ are the orders $^1$ of $b_i$ and $b^*_i$ respectively.

$^1$ A Runge-Kutta method has order $p$ if

$$
\sum^{p}_{i=1}b_i = 1 \ \ \ \ 
\sum^{p}_{i=1}b_ic_i = 1/2 \ \ \ \ 
\sum^{p}_{j=1}a_{ij} = c_i  \ \ \ \ 
$$

the same conditions defined above with the parameters `condition_b`, `condition_bc` and `condition_Ac` respectively. 

However, it is difficult to calculate the order when these conditions are not met. Therefore to simplify things, $h$ will be re-calculated in the same way whether these conditions are met or not.

Note that you can impose the same conditions on $b^*_i$ with the parameters `condition_b_star` and `condition_b_star_c` for the first two conditions. Again, these are set to `False` by default.

<h1 id="ark-eg">Adaptive Runge-Kutta example</h1>

In this example, we will try to numerically integrate

$$
\frac{\text{d}y}{\text{d}x}=-2xy
$$

will initial conditions $(-5, e^{-25})$. Note that the solution is $y=e^{-x^2}$.

Here we will use the Fehlberg RK1(2) method which has the following Butcher tableau:

$$
  \begin{array}{c| c c c}
    0 \\
    1/2 & 1/2\\
    1 & 1/256 & 255/256\\
    \hline
      & 1/512 & 255/256 & 1/512\\
      & 1/512 & 255/256 & 0\\
  \end{array}\\ \\
$$

We require an additional argument for the $b^*$ vector:

```
A = [
        ["1/2"], # Remember, fractions need to be in strings
        ["1/256", "255/256"]
    ]

b = [
    "1/512", "255/256", "1/512"
]

b_star = [
    "1/256", "255/256", 0
]

c = [
    0, "1/2", 1
]
```

These conditions will be set to `False` (which is the default) in this example.

Also, the idea of discretizations is meaningless in this context since the step size is constantly changing. Therefore, the `time_steps` parameter will be used as the *starting* time step, which we will set to 0.01 in this example.

We will set the tolerance to a punitive value of 0.00000001. _Note: if `tolerance` is not set, a default value of 0.001 will be used instead.

We can now create our `Gerk` object, `solve` and `plot`:

```
adj_rk = Gerk(
    A=A, 
    b=b, 
    c=c, 
    initial_conditions=(-5.0, math.exp(-25.0)),
    time_steps=0.001, # Remeber, this needs to be the step size NOT the number of discretizations
    final=5,
    func=lambda x,y: -2*x*y,
    real_values=lambda x: math.exp(-x**2),
    b_star=b_star,
    tolerance=0.00001
)

adj_rk.solve()

adj_rk.plot(with_real=True)
```

![rk_adaptive](https://user-images.githubusercontent.com/59436765/210466230-cbee11e0-02e5-49ab-92f1-553e302a972f.png)


In this case, the Runge-Kutta approximation is so accurate that we can barely see the exact curve!
