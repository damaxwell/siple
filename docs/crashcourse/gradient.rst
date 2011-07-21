================
Gradient methods
================

The |siple| library focusses on regularization algorithms, and in particular
gradient-based algorithms, which apply in the case that :math:`X` and
:math:`Y` are both inner product spaces. (The reader might want
to take a moment now to brush up on :ref:`inner products, adjoints, and derivatives <appendix>`.)

Suppose we wish to solve the inverse problem

.. math:: \calF(x) = y.
  :label: inv

One approach would be to consider the functional

.. math:: J(x) = \frac{1}{2}||y-\calF(x)||_Y^2 = \frac{1}{2}\ip<y-\calF(x),y-\calF(x)>_Y.

A solution of equation :eq:`inv` corresponds to a minimizer of 
:math:`J`, since :math:`J` is always non-negative and equals
zero exactly when equation :eq:`inv` holds.  So the problem 
of finding a solution of :eq:`inv` corresponds to minimizing :math:`J`.
It turns out that there are good techniques for finding approximate
minimizers of :math:`J` that regularize the inverse problem.

A natural approach to minimizing :math:`J` is to head \"downhill\".
Recall that the gradient of :math:`J` at :math:`x`
is the vector :math:`\nabla J_x` such that the derivative 
of :math:`J` at :math:`x` in the direction :math:`h`
(i.e. :math:`J_x'(x)`) is given by

.. math:: J_x'(h) = \ip< \nabla J_x, h>_X
  :label: grad1

for all :math:`h\in X`.  The vector :math:`\nabla J_x` points in the 
direction (at :math:`x`) of steepest ascent of :math:`J`.  So the 
direction of steepest descent is :math:`-\nabla J_x`.  To compute this
direction, recall

.. math::  

  J'_x(h) &= \left. \frac{d}{dt}\right|_{t=0} J(x+th)\\
  & = \frac{d}{dt}|_{t=0} \frac{1}{2} \ip<y-\calF(x+th),y-\calF(x+th)>_Y\\
  & =  2 \frac{1}{2} \ip<y-\calF(x), -\left.\frac{d}{dt}\right|_{t=0} \calF(x+th)>_Y\\
  & = -\ip< y-\calF(x), \calF'_x(h)>_Y.

Now :math:`\calF'_x` is a linear map from :math:`X` to :math:`Y` and
it has an adjoint :math:`(\calF'_x)^*` that satsfies

.. math:: -\ip< y-\calF(x), \calF'_x(h)>_Y = -\ip< (\calF'_x)^*(y-\calF(x), h)>_X

for all :math:`h\in X`.  We have just seen that

.. math:: J_x'(h) = -\ip< (\calF'_x)^*(y-\calF(x)), h>_X
  :label: grad2

Comparing equations :eq:`grad1` and :eq:`grad2` we see that

.. math:: \nabla_x J = -(\calF'_x)^*(y-\calF(x))

and the direction of steepest descent at :math:`x` is then

.. math:: -\nabla_x J = (\calF'_x)^*(y-\calF(x)).

To use a gradient method, we need to be able to compute the following:

1. Inner products on :math:`X` and :math:`Y`.
2. The forward map :math:`\calF`.
3. The linearization of :math:`\calF` at :math:`x`, i.e. :math:`\calF_x'`.
4. The adjoint of this linearization, i.e. :math:`(\calF_x')^*`.

These computations are all encapsulated in a single object called 
a :class:`NonlinearFowardProblem <siple.gradient.forward.NonlinearForwardProblem>` through the methods

1. :attr:`NonlinearFowardProblem.domainIP <siple.gradient.forward.NonlinearForwardProblem.domainIP>` and :attr:`NonlinearFowardProblem.rangeIP <siple.gradient.forward.NonlinearForwardProblem.rangeIP>`
2. :attr:`NonlinearFowardProblem.F <siple.gradient.forward.NonlinearForwardProblem.F>`
3. :attr:`NonlinearFowardProblem.T <siple.gradient.forward.NonlinearForwardProblem.T>`
4. :attr:`NonlinearFowardProblem.TStar <siple.gradient.forward.NonlinearForwardProblem.TStar>`

For problems where the forward map :math:`\calF` is in fact
a linear map, a
:class:`LinearFowardProblem <siple.gradient.forward.LinearForwardProblem>`
should be used instead.

The |siple| library contains implementation of the following 
gradient algorithms.

Linear problems:
  * Landweber iteration
  * Steepest descent
  * Conjugate gradient method
  * Conjugate gradient method applied to the normal equation

Nonlinear problems:
  * Nonlinear conjugate gradient method
  * Incomplete Gauss-Newton


Regularization: the Morozov discrepancy principle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



