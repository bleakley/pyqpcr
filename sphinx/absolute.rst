Absolute quantification
***********************

Quantity calculation
====================

Computing the linear regression of standard curve, we had :math:`a_g` (slope) and :math:`b_g` (Y-intercept) for each gene :math:`g` for the straight line of equation: 

.. math::
   c_t=a_{g}\log Q + b_g

For each replicate :math:`(s,\ g)`, the absolute quantity :math:`AQ_{sg}` is calculated using the linear regression associated with the gene :math:`g`:

.. math::
   \log(AQ_{sg}) = \dfrac{{\overline{c_t}}_{sg} - b_g}{a_g} \Longrightarrow
   \boxed{AQ_{sg} = 10^{\frac{{\overline{c_t}}_{sg} - b_g}{a_g}}}

To calculate the error we propagate the errors of a law :math:`z=10^y`. 
Using logarithm, we obtain

.. math::
   \ln z = y\ln 10

Differentiating this expression, we get

.. math::
   \dfrac{\delta z}{z} = \delta y \ln 10 \quad\text{with}\quad y = \dfrac{{\overline{c_t}}_{sg}- b_g}{a_g},

that leads to the following error 

.. math::
   \text{SE}(AQ_{sg}) = \dfrac{AQ_{sg}\ln(10)}{a_g}\sqrt{{\text{SE}({c_t}_{sg})}^2 + \text{SE}(b_g)^2 +
   \left(\dfrac{\text{SE}(a_g)}{a_g}({\overline{c_t}}_{sg}-b_g)\right)^2}
                       
