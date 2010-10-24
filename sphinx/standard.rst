Standard curve calculation
**************************

Linear regression
-----------------


.. math::
   y\equiv c_t \quad \text{and} \quad x\equiv\log Q

For a straight line :math:`y=ax+b`, variance :math:`\sigma_x^2` is

.. math::
   \sigma_x^2=\displaystyle\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\overline{x})^2
   \quad\text{with}\quad \overline{x}=\frac{1}{n}\sum_{i=1}^{n}x_i

Similarly, covariance :math:`\sigma_{xy}` is defined by


.. math::
   \sigma_{xy} =
   \displaystyle\frac{1}{n-1}\sum_{i=1}^{n}(x_i-\overline{x})(y_i-\overline{y})

To obtain the slope :math:`a` and Y-intercept :math:`b` of a linear regression

.. math::
   \left\lbrace
   \begin{array}{l}
    a = \dfrac{\sigma_{xy}}{\sigma_x^2}=\dfrac{\sum(x_i-\overline{x})
        (y_i-\overline{y})}{\sum(x_i-\overline{x})^2} \\ \\
    b = \overline{y}-a\overline{x}
    \end{array}
    \right.

Slope error
-----------

For a problem with :math:`n-2` degrees of freedom, the standard deviation of a linear regression can be written:

.. math::
   \sigma_{\epsilon}^2 = \displaystyle\frac{1}{n-2}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2

Standard error :math:`\text{SE}(a)^2` on the slope is deducted :

.. math::
   \text{SE}(a)^2 = \dfrac{\sigma_{\epsilon}^2}{(n-1)\sigma_x^2} \Longrightarrow
   \boxed{\text{SE}(a) =
   \sqrt{\dfrac{\sum(y_i-\hat{y}_i)^2}{
   (n-2)\sum (x_i-\overline{x})^2}}}

Y-intercept error
-----------------
.. math::
   \text{SE}(b) = \text{SE}(a)\sqrt{\sum_{i=1}^{n}\dfrac{(x_i)^2}{n}}


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
