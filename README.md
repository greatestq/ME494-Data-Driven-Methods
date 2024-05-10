# ME494-Data Driven Methods
 2024 Spring Independent Study on Data Driven Methods in Mechanical Engineering

This repository contains programming assignments and results for independent study ME493: Methods of Data-Driven Control. This independent study provides introduction to data-driven methods in the scope of mechanical engineering. Data-driven methods such as Proper Orthogonal Decomposition (POD), Dynamic Mode Decomposition (DMD), Eigen-Realization Algorithm (ERA), Sparse Identification of Non-Linear Dynamics (SINDY) are covered as part of introduction lectures. Then, a final project of one's choice is investigated.  

My colleague, Ben Aziel, has compiled a document describing the course's contents. The document can be viewed here: https://azielben.quarto.pub/me493-ddc/

Below are some of the results that I have collected when trying to be more familiarized with the concepts that were introduced during the lectures.

## POD Analysis 
To practice POD analysis, we have taken a data set from  http://deepblue.lib.umich.edu/data/collections/kk91fk98z. The data set that I chose for this assignment is low Reynolds number airfoil DNS with angle of attack of 25 degrees and pitching frequency of 0.05 degrees. 

### Plot of eigenvalues vs mode index
![Plot of Eigenvalues vs Mode Index](images/singular_values.png)
![Plot of Eigenvalues vs Mode Index Zoomed in](images/singular_values_zoom.png)

A closer look at the eigenvalue plot shows that the singular values of the data matrix are coming in a pair. This is an indication of an oscillatory pattern inside the fluid flow.

![Plot of cumulative energy](images/cumsum.png)
The above plot shows that the first 6 "modes" are enough to reconstruct 95% of the entire fluid flow.

### Plot of the first 6 temporal amplitudes. What is the frequency associated with oscillations?
![temporal amplitude](images/temporal_amplitdues.png)
The first two mode temporal amplitudes show about 45 degrees offset meaning the main dynamics of this fluid flow can be expressed by a combination of sine and cosine. Although the 2nd mode has nonlinear harmonic pattern, the two modes can be approximated to be almost perfect sinusoidal. The subsequent modes resemble the "beats" shape. I am not yet sure the physical intution behind such pattern but these subsequent modes also have the tendency of being a pair with approximately 45 degrees offset.

### Plot of the first 6 spatial modes
![ux spatial modes](images/ux_spatial%20modes.png)
![uy spatial modes](images/uy_spatial%20modes.png)

Now, the temporal amplitudes give us the frequency $f$. Then, once we have details of the flow condition, the fluid speed $c$ can also be computed. Then, $\lambda = c / f$ will give the wavelength associated with each oscillatory mode.

### Strouhal Number
The Strouhal number is a dimensionless number describing oscillating flow mechanisms. Typical fluid flow past cylinder should yied St $= 0.2$. Once the frequency $f$ associated with the dominant oscillatory mode pair is computed, the Strouhal number can be computed via $\frac{fL}{U}$ where L is the airfoil thickness, and U is the average flow velocity. 

### Fluid flow Reconstruction
Due to the file size limit, the reconstructed video can be found here: https://www.youtube.com/watch?v=-7NFyvJYAgo
