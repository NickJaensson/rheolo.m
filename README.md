# rheolom


This Matlab code can be used to investigate the rheological characteristics of non-linear viscoelastic models. 

The file `main_rate.m` can be used to perform rate-controlled simulations (startup and/or steady-state)

The file `main_stress.m` can be used to perform stress-controlled simulations (startup)

Three types of flow are possible:
* simple shear flow
* planar extensional flow
* uniaxial extensional flow

All relevant data is available in the `rheodata` structure for postprocessing. Some plotting is avaiable in the file `rheoplot.m`.

**Viscoelastic models**

In all the viscoelastic models available the stress is given by

$$\boldsymbol{\tau} = G(\boldsymbol{c}-\boldsymbol{I})$$

with the evolution of $\boldsymbol{c}$ given by

$$\frac{d \boldsymbol{c}}{dt}=
(\nabla \boldsymbol{v})^T \cdot \boldsymbol{c}+\boldsymbol{c} 
\cdot (\nabla \boldsymbol{v})+\boldsymbol{f}(\boldsymbol{c})$$

with 
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}(\boldsymbol{c}-\boldsymbol{I})$ (upper convected Maxwell model)
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}(\boldsymbol{c}-\boldsymbol{I}+\alpha(\boldsymbol{c}-\boldsymbol{I})^2)$ (Giesekus model)
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}Y(\text{tr}(\boldsymbol{c}))(\boldsymbol{c}-\boldsymbol{I})$ 
    * with $Y(\text{tr}(\boldsymbol{c}))=1+\epsilon(\text{tr}(\boldsymbol{c})-3)$ (linear PTT model) 
    * with $Y(\text{tr}(\boldsymbol{c}))=e^{\epsilon(\text{tr}(\boldsymbol{c})-3)}$ (exponential PTT model)

**Adapted lambda models**

In adition to the viscoelastic models given above, there is also an option to simulate "adapted lambda" models by adapting the terms with $\lambda$ in the viscoelastic models above. The options are:

* A purely elastic model is obtained by letting $\lambda$ go to $\infty$, or alternatively replacing $1/\lambda$ by
$$\frac{1}{\lambda} \quad\rightarrow\quad 0$$

* The first Saramito model (SRM1) is obtained by replacing $1/\lambda$ by
$$\frac{1}{\lambda} \quad\rightarrow\quad \frac{1}{\lambda}\text{max}(0,\frac{\tau_\text{d}-\tau_\text{y}}{\tau_\text{d}})$$
where $\tau_\text{y}$ is the yield shear stress and $\tau_\text{d}$ the von Mises equivalent shear stres given by 
$|\boldsymbol{\tau}|=\sqrt{\frac{1}{2}\boldsymbol{\tau}:\boldsymbol{\tau}}$.

* The second Saramito model (SRM2) is obtained replacing $1/\lambda$ by
$$\frac{1}{\lambda} \quad\rightarrow\quad G\\,\text{max}(0,\frac{\tau_\text{d}-\tau_\text{y}}{K \tau_\text{d}^n})^{1/n} = 0 \quad \text{if} \\, \tau_\text{d} \leq \tau_\text{y}$$
$$\frac{1}{\lambda} \quad\rightarrow\quad G\\,\text{max}(0,\frac{\tau_\text{d}-\tau_\text{y}}{K \tau_\text{d}^n})^{1/n} = \frac{G}{\tau_\text{d}}(\frac{\tau_\text{d}-\tau_\text{y}}{K})^{1/n} \quad \text{if} \\, \tau_\text{d} > \tau_\text{y}$$
which creates a "power-law" type model with power law index $n$ and consistency index $K$.

