Absolute quantification
***********************

Calculation of quantities
=========================

Computing the linear regression of standard curve, we had :math:`a_g` (slope) and :math:`b_g` (Y-intercept) for each gene :math:`g` :

.. math::
   c_t =a_{g}\log_{10} Q  + b_g.

For each replicate :math:`(s,\ g)`, the absolute quantity :math:`AQ_{sg}` is
calculated using the linear regression associated with the gene :math:`g`

.. math::
   \log_{10}(AQ_{sg}) = \dfrac{{\overline{c_t}}_{sg} - b_g}{a_g} \Longrightarrow
   \boxed{AQ_{sg} = 10^{\frac{{\overline{c_t}}_{sg} - b_g}{a_g}}},
   :label: qabs

where overbars denote average of :math:`{c_t}_{sg}`. 

Standard error
==============

To calculate the corresponding standard error, we propagate the errors of a law :math:`z=10^y`, that is

.. math::
   \ln z = y\ln 10.

Differentiating this expression, we get

.. math::
   \dfrac{\delta z}{z} = \delta y \ln 10 \quad\text{with}\quad y = \dfrac{{\overline{c_t}}_{sg}- b_g}{a_g},

where :math:`\delta y` can be obtained by usual error propagations of differences and ratios

.. math::
   \dfrac{\delta y}{y} = \sqrt{\dfrac{{\text{SE}({c_t}_{sg})}^2 + 
   \text{SE}(b_g)^2}{({\overline{c_t}}_{sg}-b_g )^2}+\left( \dfrac{\text{SE}(a_g)}{a_g}\right)^2},

leading to the following error

.. math::
   \boxed{\text{SE}(AQ_{sg}) = AQ_{sg}\ln AQ_{sg}\sqrt{\dfrac{{\text{SE}({c_t}_{sg})}^2 + 
   \text{SE}(b_g)^2}{({\overline{c_t}}_{sg}-b_g )^2}+\left( \dfrac{\text{SE}(a_g)}{a_g}\right)^2}}.
   :label: qabserror

Normalisation of quantities
===========================

Unfortunately, Eq. :eq:`qabserror` directly depends on the starting amounts, implying that different starting
units would lead to different relative errors. It is thus necessary to rescale
Eq. :eq:`qabs` to work with a quantity that does not depend on the starting
units. Thus we replace :math:`AQ_{sg}` by the dimensionless ratio

.. math::
   Q_{sg} = \dfrac{AQ_{sg}}{\overline{Q_g}} \quad\text{with}\quad
   \overline{Q_g} = \dfrac{1}{n}\sum_{i=1}^{n}({Q_g}_i) 
   
where :math:`\overline{Q_g}` is the average of the starting amounts. The standard error
SE(:math:`Q_{sg}`) then becomes independant on the starting quantities

.. math::
   \text{SE}(Q_{sg}) = Q_{sg}\ln Q_{sg}\sqrt{\dfrac{{\text{SE}({c_t}_{sg})}^2 + 
   \text{SE}(b_g')^2}{({\overline{c_t}}_{sg}-b_g')^2}+\left( \dfrac{\text{SE}(a_g)}{a_g}\right)^2},

where :math:`b_g` is the rescaled Y-intercept defined by

.. math::
   b_g' = b_g + a_g\log_{10} \overline{Q_g}. 

The final error on :math:`AQ_{sg}` is then expressed by

.. math::
   \boxed{\text{SE}(AQ_{sg}) = \text{SE}(Q_{sg})\overline{Q_g}}
   :label: newqabserror

that is also independant on the starting units. 
