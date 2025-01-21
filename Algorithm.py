import numpy as np
U = GLPM(n+1,q)[:-1, :]

def general_construction_algorithm(X0, b, a):
    """
    General construction algorithm for one arbitrary linear constraint.
    
    Parameters:
    X0: numpy.ndarray
        The experiment domain (design space).
    b: float
        The constant in the linear constraint equation.
    a: numpy.ndarray
        Coefficients in the linear constraint equation.
    u: Values of uniformly distributed random vectors over [0,1]^{s}
        
    Returns:
    Pn: numpy.ndarray
        The final design points with the smallest CCD.
    """
    # Step 1: Compare the magnitude of b with 0
    if b == 0:
        # Proceed to Step 2
        pass
    else:
        # Divide the constraint equation by b
        a = a / b
        b = 1
    
    # Step 2: Check all a_i > 0
    if np.all(a > 0):
        return None  # X0 is empty
    elif np.all(a <= 0):
        return None  # X0 is empty
    else:
        # Reorder
        positive_indices = np.where(a > 0)[0]
        negative_indices = np.where(a <= 0)[0]
        a_sorted_indices = np.concatenate([positive_indices, negative_indices])
        a = a[a_sorted_indices]
        X0_prime = X0[:, a_sorted_indices]  # Reorder columns

        x = [0] * q 
        F1_(1)(x[0]) = u[0]
        for j in range(1, m - 2):
           F1_(j|1,...,j-1)(x[j]) = u[j]
        x[m-1] = - 1 / a[m - 1] * sum(a[i] * x[i] for i in range(m - 1))
        for k in range(m, q-1):
          x[k] = u[k - 1]
          x
    
    # Step 3
    if np.all(a > 0):
       if a[0] < 1:
                    x = [0] * q 
                    F2_(2)(x[1]) = u[0]
                    for j in range(2, m - 1):
                      F2_(j|2,...,j-1)(x[j]) = u[j-1]
                    x[0] = 1/a[0]- 1 / a[0] * sum(a[i] * x[i] for i in range(1,m-1))
                    for k in range(m, q-1):
                      x[k] = u[k - 1]
                      x
       elif a[0] >= 1:
                    x = [0] * q  
                    prod_1 = 1  
                    for j in range(m - 1):
                      prod_1 *= (1 - u[j])
                      x[0] = (1 / a[0]) * (prod_1 ** (1 / (m - 1)))
                    for k in range(2, m + 1):
                      prod_k = 1
                    for j in range(k - 2):
                     prod_k *= (1 - u[j])
                     term_1 = 1 - (1 - u[k - 2]) ** (1 / (m - k + 1))
                     x[k - 1] = (1 / a[k - 1]) * (term_1 * prod_k ** (1 / (m - k)))
                    for k in range(m + 1, q + 1):
                     x[k - 1] = u[k - 2]
                    x
       else:
            x = [0] * q 
            F4_(1)(x[0]) = u[0]
            for j in range(1, m - 2):
             F4_(j|2,...,j-1)(x[j]) = u[j]
             x[m-1] = 1 / a[m - 1]- 1 / a[m - 1] * sum(a[i] * x[i] for i in range(m - 1))
            for k in range(m, q-1):
              x[k] = u[k - 1]
            x
    
    # Compute CCD values and select the design with the smallest CCD
    ccd_values = [CCD2(x) for design in x]
    Pn = designs[np.argmin(ccd_values)]
    
    return Pn

