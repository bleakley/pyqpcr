Relative quantification
***********************

Standard curve calculation
==========================

Linear regression
-----------------


.. math::
   y\equiv c_t \quad \text{et} \quad x\equiv\log Q

For a straight line :math:`y=ax+b`, variance :math:`\sigma_x^2` is

.. math::
   \sigma_x^2=\displaystyle\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\bar{x})^2
   \quad\text{avec}\quad \bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i

Similarly, covariance :math:`\sigma_{xy}` is defined by


.. math::
   \sigma_{xy} =
   \displaystyle\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\bar{x})(y_i-\bar{y})

To obtain the slope :math:`a` and Y-intercept :math:`b` of a linear regression

.. math::
   \left\lbrace
   \begin{array}{l}
    a = \dfrac{\sigma_{xy}}{\sigma_x^2}=\dfrac{\sum(x_i-\bar{x})
        (y_i-\bar{y})}{\sum(x_i-\bar{x})^2} \\ \\
    b = \bar{y}-a\bar{x}
    \end{array}
    \right.

Slope error
-----------

For a problem with :math:`n-2` degrees of freedom, the standard deviation ? of a linear regression can be written:

.. math::
   \sigma_{\epsilon}^2 = \displaystyle\frac{1}{n-2}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2

Standard error :math:`\text{SE}(a)^2` on the slope is deducted :

.. math::
   \text{SE}(a)^2 = \dfrac{\sigma_{\epsilon}^2}{(n-1)\sigma_x^2} \Longrightarrow
   \boxed{\text{SE}(a) =
   \sqrt{\dfrac{\sum(y_i-\hat{y}_i)^2}{
   (n-2)\sum (x_i-\bar{x})^2}}}

Pearsson Coefficient
--------------------

.. math::
   R = \dfrac{\sigma_{xy}}{\sigma_x\sigma_y}

:math:`\text{SE}(a)` can be express as function of :math:`R^2` :

.. math::
   \text{SE}(a) = \dfrac{\sigma_y}{\sigma_x}\sqrt{\dfrac{1-R^2}{n-2}}

Confidence Interval
-------------------

For a Student t-test on the expectation with an unknown standard deviation, 
for a given confidence level :math:`\alpha`, the error on :math:`a` is:

.. math::
   \Delta a = \text{SE}(a) \cdot t_{(1-\alpha)/2}^{n-2}

where :math:`t_{(1-\alpha)/2}^{n-2}` is the quantile of order 
:math:`\alpha/2` of the Student law with :math:`n-2` degrees of freedom. 
In term of probabilities, we have:

.. math::
   \boxed{%
   P\left[a- \text{SE}(a) \cdot t_{(1-\alpha)/2}^{n-2}\le \beta\le a +\text{SE}(a)
   \cdot
   t_{(1-\alpha)/2}^{n-2}\right] = \alpha}

Efficiency
----------

.. math::
   \text{eff} = 10^{-1/a}-1\quad \text{or in \%:}\quad \text{eff}
   =100\cdot\left(10^{-1/a}-1\right)

The error on the slope :math:`\Delta a` can be propagated to the 
calculation of the error on the efficiency :

.. math::
   \epsilon(\text{eff}) = \ln 10(\text{eff}+100) \dfrac{\Delta a}{a^2}

ou bien ???

.. math::
   \epsilon(\text{eff}) = \ln 10(\text{eff}+100) \dfrac{\text{SE}(a)}{a^2}

Relative quantification
=======================

For a replicate :math:`(s,\ g)`, we calculate :math:`c_t` mean of all wells, i.e.

.. math::
   {\bar{c_t}}_{sg} = \dfrac{1}{n_w}\sum_{w}{c_t}_{wsg},

and the associated standard error

.. math::
   \text{SE} ({c_t}_{sg}) =
   \sqrt{\dfrac{1}{n_w(n_w-1)}\sum_{\text{w}}
   ({c_t}_{wsg}-{\bar{c_t}}_{sg})^2}

The associated confidence interval can be written

.. math::
   \left[{\bar{c_t}}_{sg}-t_{(1-\alpha)/2}^{n_w-1}\text{SE}(c_t),\
   {\bar{c_t}}_{sg}+t_{(1-\alpha)/2}^{n_w-1}\text{SE}(c_t)\right]

For each gene, the coefficient :math:`{c_t}_{\text{ref}}` is calculated

.. math::
   {c_t}_{\text{ref},s} = \dfrac{1}{n_s}\sum_{s}{c_t}_{sg} 

For each replicate we define :math:`{\Delta c_t}_{sg}`


.. math::
   {\Delta c_t}_{sg} = {c_t}_{\text{ref},s} - {\bar{c_t}}_{sg},

and :math:`\text{RQ}_{sg}`

.. math::
   \text{RQ}_{sg} = \text{eff}_g^{{\Delta c_t}_{sg}}

To calculate the error we propagate the errors of a law :math:`z=x^y`. 
Using logarithm, we obtain

.. math::
   \ln z = y\ln x

Differentiating this expression, we have

.. math::
   \dfrac{\delta z}{z} = y\dfrac{\delta x}{x}+\delta y \ln x,

which conduct to the following error

.. math::
   \dfrac{\text{SE}(z)}{\bar{z}} =\sqrt{%
   \left(\bar{y}\dfrac{\text{SE}(x)}{\bar{x}}\right)^2+\left(\text{SE}(y) \ln
   \bar{x}\right)^2}

When replacing x, y and z with magnitudes of interest

.. math::
   \text{SE}(\text{RQ}_{sg}) =\text{RQ}_{sg}\sqrt{%
   \left({\Delta c_t}_{sg}\dfrac{\text{SE}(\text{eff}_g)}{\text{eff}_g}\right)^2+
   \left(\text{SE}({\Delta c_t}_{sg}) \ln \text{eff}_g\right)^2}

Reference genes
===============

We perform a geometric mean of reference genes RQ. 
For each sample we define :math:`\text{NF}_{s}` as

.. math::
   \text{NF}_s = \left( \prod_{p=1}^{n_{\text{generef}}} \text{RQ}_{ps}
   \right)^{1/n_{\text{generef}}}

With the standard error

.. math::
   \text{SE}(\text{NF}_s) = \text{NF}_s  \sqrt{\sum_{p=1}^{n_{\text{generef}}}
   \left( \dfrac{\text{SE}(\text{RQ}_{ps})}{n_{\text{generef}}\cdot 
   \text{RQ}_{ps}} \right)^2}


RQ normalization
================

We normalize in comparison to reference sample:

.. math::
   \text{NRQ}_{gs} =
   \dfrac{\text{RQ}_{g,s}}{\text{NF}_s}
   \cdot \dfrac{\text{NF}_{\text{echref}}} {\text{RQ}_{g,\text{echref}}}

The standard error is calculated

.. math::
   \text{SE}(\text{NRQ}_{g,s}) = \sqrt{%
   \left(\dfrac{\text{SE}(\text{RQ}_{g,s})}{\text{RQ}_{g,s}} \right)^2 +
   \left(\dfrac{\text{SE}(\text{NF}_{s})}{\text{NF}_{s}} \right)^2 +
   \left(\dfrac{\text{SE}(\text{NF}_{\text{echref}})}{\text{NF}_{\text{echref}}}
   \right)^2 +
   \left(\dfrac{\text{SE}(\text{RQ}_{g,\text{echref}})}
   {\text{RQ}_{g,\text{echref}}
   } \right)^2}

