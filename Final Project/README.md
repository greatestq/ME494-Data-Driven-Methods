As per my interest, I have decided to study fluid mechanics as part of my final project. In specific, Dr. Luchtenburg and I are investigating what is known as "Zone-SINDY". We believe this could lead into various fluid-mechanics related concepts such as physics-informing methodology, sensor placement test, and vortex mechanism investigation. As we decided to aim for publication, this final project has been extended to 2024 Summer. This document will be updated simulataneously during the summer of 2024.

## Introduction 
    The governing equations that we are after is the vorticity equation
    $
        \frac{D\omega}{Dt} = \frac{\partial \omega}{\partial t} + (u \cdot \nabla)\omega= (\omega \cdot \nabla) u + \nu \nabla^2\omega
    $ 
    which is also known as the \textbf{vorticity transport equation} which is valid for incompressible flow with conservative body forces. $(u \cdot \nabla) \omega$ is the convection term (motion of the fluid particle as it moves from one point to another), $(\omega \cdot \nabla) u$ is the stretching or tilting of vorticity due to the flow velocity gradients, and $\nu \nabla^2 \omega$ is the viscous effects.
    
As Dr. Brunton and his collaborators have stated in their work, SINDY has been found to be particularly useful in finding the governing PDEs. However, their examples, at least for fluid mechanics, are currently limited to Re=100 flow past circular cylinder.

As part of preliminary testing, I am diving the numerical simulation grid into 24 different "zones", then running SINDY locally to see what the governing equations may look like.
![Preliminary Zone Division](Images/First%20Grid%20Test%20Schematic.png)
We expect to see the vortex convection terms near the circular cylinder. Then further down the wakes, we should be seeing the laplcian terms starting to emerge as convection should be, physically speaking, less dominant further away from the cylinder.

Below are the results from preliminary testing:

    PDE-Find results ($w_t = 0.009951w_{xx}+0.009877w_{yy}-0.990566uw_x - 0.986856 v w_y)$:
    \begin{align*}
    \text{Zone}_{11}: w_t &= -10.373580w_x + 0.005845w_{yy} + 14.588403uw_x -5.206192u^2w_x \\
    \text{Zone}_{12}: w_t &= -0.814525uw_x -0.894574vw_y \\
    \text{Zone}_{13}: w_t &= -0.814948uw_x - 0.895500vw_y \\
    \text{Zone}_{14}: w_t &= -3.642810w_x + 0.006472w_{yy}+2.693152uw_x
    \\
    \text{Zone}_{21}: w_t &= 0.009943w_{yy}-0.960680uw_x -0.920402vw_y \\
    \text{Zone}_{22}: w_t &= 0.013251w_{yy}-0.998260uw_x -0.976469vw_y \\
    \text{Zone}_{23}: w_t &= 0.013038w_{yy}-0.995063uw_x-0.975659vw_y \\
    \text{Zone}_{24}: w_t &= 0.009897w_{yy}-0.958806uw_x-0.919833vw_y \\
    \text{Zone}_{31}: w_t &= -0.166340w_x + 0.010829w_{yy}-0.636702uw_x -0.165873u^2w_x -0.937844vw_y \\
    \text{Zone}_{32}: w_t &= 0.062814 - 0.085894w_x -0.326263v +0.623483uv-0.184306u^2 + 0.139838wu \\ &-0.688405uw_x -0.238665u^2w_x + 0.027628w^2w_x -0.702141vw_y -0.040227wvw_{xx} \\ &+0.165052vw_{xy}-0.014496ww_{xy}-0.256414uvw_{xy}+0.086349uw_{yy}-0.313592v^2w_{yy}\\ &-0.074420u^2w_{yy} \\
    \text{Zone}_{33}: w_t &= 0.012044w_{xx} - 0.099330u + 0.112511u^2 - 0.882673uw_x -0.680702vw_y \\ &- 0.654039uvw_y - 0.054315vw_{xx}+0.126289uvw_{xx} +0.131935vw_{xy}\\ &-0.268053uvw_{xy}+0.048027uw_{yy}-0.083223v^2w_{yy}-0.036844u^2w_{yy}
    \\
    \text{Zone}_{34}: w_t &= -0.744346uw_x
    \\
    \text{Zone}_{41}: w_t &= -0.438828w_x + 0.062870w_y -0.084233w_{xx} + 0.069863w_{xy} + 1.202151wv \\ &-0.437126uw_x -0.692593uvw_y + 0.070448uw_{xx}-0.063853u^2w_{xy}
    \\
    \text{Zone}_{42}: w_t &= -0.209069-0.417345w_x + 0.511569u - 0.279118u^2-0.596061u^2w_x\\ &-1.199597uvw_y+0.087385vw_{xx}-0.138471uvw_{xx}+0.212975vw_{xy} \\ &-0.025732ww_{xy}-0.334708uvw_{xy} \\
    \text{Zone}_{43}: w_t &= -0.395714w_x -0.637701u^2w_x-1.290394uvw_y-0.096956vw_{xx}+0.164380uvw_{xx} \\ &+0.311923vw_{xy}-0.483532uvw_{xy} \\
    \text{Zone}_{44}: w_t &= 0.012208w_{yy}-1.001086uw_x - 1.05546uvw_{y}
    \\
    \text{Zone}_{51}: w_t &= 0.041899w_{xy}-0.957095uw_x-1.087975vw_y =0.384334vw_{xx}+0.387556uvw_{xx} \\ &-0.049041u^2w_{xy} \\
    \text{Zone}_{52}: w_t &= -0.214171w_x + 0.708708u - 0.492641u^2 - 0.743324uw_x -1.023232uvw_y \\ &+0.348350vw_{xy}-0.051548ww_{xy}-0.465857uvw_{xy}+0.103566ww_{yy} \\ &-0.141991wuw_{yy} \\
    \text{Zone}_{53}: w_t &= 0.378146 -0.203851w_x -0.946217u +0.585654u^2 + 0.056877wu - 0.621635uw_x \\ &-0.170536u^2w_x + 0.133084wvw_x -0.914885uvw_y - 0.060613v^2w_{xx}+0.021458u^2w_{xx} \\ &-0.062937wvw_{xx}+0.340736vw_{xy}-0.119187ww_{xy}-0.444142uvw_{xy}+0.129094wuw_{xy}\\&-0.44919v^2w_{yy}+0.153305wvw_{yy} 
    \\
    \text{Zone}_{54}: w_t &= 0.105893 + 0.559789w_x -0.022969w_{xy}-0.292511u+0.191289u^2 \\ &-0.008664wu -2.297540uw_x+0.720216u^2w_x-1.048583vw_y+0.268753vw_{xx} \\ &-0.242326uvw_{xx}+0.263605vw_{xy}-0.305971uvw_{xy}+0.018159u^2w_{xy} \\ &+0.010341u^2w_{yy} \\
    Zone_{61}: w_t &=-0.886684uw_x
    \end{align*}
\begin{align*}
 Zone_{62}: w_t &=-0.708143w_x + 0.413029u -0.528047u^2 +0.132803ww_x -0.240371wuw_x \\ &+0.611733vw_{xy}-0.848452uvw_{xy} \\
 Zone_{63}: w_t &= -1.030841w_x -0.375420w_y -0.386098u +0.480280u^2 + 0.281571uw_x \\
 &-0.354807ww_x + 0.513594wuw_x +0.109391vw_{xx}-0.140507uvw_{xx} \\ &+0.456267vw_{xy} -0.302760ww_{xy}-0.639749uvw_{xy}+0.299630wuw_{xy} \\
 Zone_{64}: w_t &= 0.246817 -0.879589w_x - 0.199581w_y -0.593545u + 0.360972u^2 \\ &-0.426931wu + 0.312007vw_{xx}-0.291659uvw_{xx}+0.037428wuw_{xx}+0.505514vw_{xy} \\ &-0.566960uvw_{xy}
 \end{align*}

 Now, what is not obvious are the regions inbetween the extreme ends. Here, we are seeing "non-sparsity". In other words, there are many terms that are not present in the actual vorticity transport equation. Therefore, we played around with the $\lambda$ parameter to tweak the algorithm sensitivity. However, the inner zone's "non-sparsity" are still present. Based on these findings, we are now proposing three different perspectives that we would like to take on over the summer.

 1. Using the weak formulation to compute the spatial derivatives and seeing whether results are "derivative-invariant". In other words, can we test different methods' equivalence?

 2. Dr. Philip Yecko has provided us with a MATLAB code that simulates inner ocean geophysical flow with Munk condition. This type of flow is famous for having distinct regions of different flow structures. We are hoping to take our "zone-SIDNY" idea and see if we can extract the dominant physics in different parts of the flow.

 3. There is a possibility that the "non-sparsity" in the inner-regions represent the actual non-linear mechanism that generates vortices. Since our data is 2-dimensional and there is no vortex stretching in 2D, these non-sparse terms may represent the vortex generation/maintanence mechanism that is not well known in turbulence. 

This document will continue to be updated over the summer.

Last Update: May 10, 2024