# reaction-diffusion

This parallelised C++ code solves a reaction-diffusion equation. The particular system of PDEs that are solved is a variant of the Barkley model and is given by

$$
\begin{aligned}
\frac{\partial u}{\partial t}-\mu_1 \nabla^2 u=f_1(u, v) \\
\frac{\partial v}{\partial t}-\mu_2 \nabla^2 v=f_2(u, v)
\end{aligned}
$$

where $\mu_1, \mu_2$ are coefficients of diffusion and the reaction terms are given by

$$
\begin{aligned}
f_1(u, v) &=\epsilon u(1-u)\left(u-\frac{(v+b)}{a}\right) \\
f_2(u, v) &=u^3-v
\end{aligned}
$$

There is a Makefile which runs the code by using the command `make`. When running the code it creates a file called output.txt which contains the $u$ and $v$ velocity
at each co-ordinate in the domain.

The $\mu_1, \mu_2$ coefficients of diffusion can be changed. One example of the results of this code using $\mu_1=1, \mu_2=0.01$ is shown in graphical form below where the $u$ velocity has been plotted.

![test4, close to 1 everywhere with pockets](https://user-images.githubusercontent.com/59830842/195109738-4da837b4-4599-4633-b6cc-8e1bd987faf2.png)

