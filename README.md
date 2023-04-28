This repository contains code to solve the normalized Usadel equations

\begin{align*}
    &-\nabla\cdot\left(\hat{G}\nabla\hat{G}\right) + \left[\hat{\tau}^3E+i\hat{\Delta},\hat{G}\right] = 0,\\
    &\Delta\text{ln}(t) +t\sum_{n=-\infty}^\infty\left(\frac{\Delta}{\omega_n}-F_1(\omega_n,\mathbf{r})\right)
\end{align*}
For the full theoretical framework underlying the code we refer to (https://essay.utwente.nl/86403/).

The code can be run using the program Matlab and requires installation of the partial differential equations toolbox, the optimization toolbox and the symbolic math toolbox.

The code is designed to solve the Usadel on a layered geometry, but other geometries could be considered as well by changing the `GenerateGeometry.m` file in the general folder.

To run the code, navigate to one of the parametrization folders and run the file `main.m`.

The copyright for all work within this repository is the Creative Commons Attribution 4.0 International (CC BY 4.0) license
