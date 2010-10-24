Relative quantification
***********************

RQ calculation
===============

For a replicate :math:`(s,\ g)`, we calculate :math:`c_t` mean of all wells, i.e.

.. math::
   {\overline{c_t}}_{sg} = \dfrac{1}{n_w}\sum_{w}{c_t}_{wsg},

and the associated standard error

.. math::
   \text{SE} ({c_t}_{sg}) =
   \sqrt{\dfrac{1}{n_w(n_w-1)}\sum_{w}
   ({c_t}_{wsg}-{\overline{c_t}}_{sg})^2}

The associated confidence interval can be written

.. math::
   \left[{\overline{c_t}}_{sg}-t_{(1-\alpha)/2}^{n_w-1}\text{SE}(c_t),\
   {\overline{c_t}}_{sg}+t_{(1-\alpha)/2}^{n_w-1}\text{SE}(c_t)\right]

For each gene, the coefficient :math:`{c_t}_{\text{ref}}` is calculated

.. math::
   {c_t}_{\text{ref},s} = \dfrac{1}{n_s}\sum_{s}{c_t}_{sg} 

For each replicate we define :math:`{\Delta c_t}_{sg}`


.. math::
   {\Delta c_t}_{sg} = {c_t}_{\text{ref},s} - {\overline{c_t}}_{sg},

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
   \dfrac{\text{SE}(z)}{\overline{z}} =\sqrt{%
   \left(\overline{y}\dfrac{\text{SE}(x)}{\overline{x}}\right)^2+\left(\text{SE}(y) \ln
   \overline{x}\right)^2}

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

