# rheolom

**Viscoelastic models**

$$\frac{d \boldsymbol{c}}{dt}=
(\nabla \boldsymbol{v})^T \cdot \boldsymbol{c}+\boldsymbol{c} 
\cdot (\nabla \boldsymbol{v})+\boldsymbol{f}(\boldsymbol{c})$$

with 
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}(\boldsymbol{c}-\boldsymbol{I})$ (UCM)
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}(\boldsymbol{c}-\boldsymbol{I}+\alpha(\boldsymbol{c}-\boldsymbol{I})^2)$ (Giesekus)
* $\boldsymbol{f}(\boldsymbol{c})=-\frac{1}{\lambda}Y(\text{tr}(\boldsymbol{c}))(\boldsymbol{c}-\boldsymbol{I})$ 
    * with $Y(\text{tr}(\boldsymbol{c}))=1+\epsilon(\text{tr}(\boldsymbol{c})-3)$ (linear PTT) 
    * with $Y(\text{tr}(\boldsymbol{c}))=e^{\epsilon(\text{tr}(\boldsymbol{c})-3)}$ (exponential PTT)
