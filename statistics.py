import numpy as np

def deviationtau(n):
    """calculate the prefactor for error margins of n-point sample"""
    if n<=2:
        n = 2
    P = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262,
         2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101,
         2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052,
         2.048, 2.045, 2.042]
    if n > 31:
        return 2.0 / np.sqrt(n)
    else:
        return P[n - 2] / np.sqrt(n)