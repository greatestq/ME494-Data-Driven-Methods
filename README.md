# ME494-Data Driven Methods
 2024 Spring Independent Study on Data Driven Methods in Mechanical Engineering

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