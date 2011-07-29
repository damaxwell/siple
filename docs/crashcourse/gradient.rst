.. _gradient:

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

.. math:: J(x) = ||y-\calF(x)||_Y^2 = \ip<y-\calF(x),y-\calF(x)>_Y

which is sometimes called the :dfn:`misfit functional`.
A solution of equation :eq:`inv` corresponds to a minimizer of 
:math:`J`, since :math:`J` is always non-negative and equals
zero exactly when equation :eq:`inv` holds.  So the problem 
of finding a solution of :eq:`inv` corresponds to minimizing :math:`J`.
It turns out that there are good techniques for finding approximate
minimizers of :math:`J` that regularize the inverse problem.

A natural approach to minimizing :math:`J` is to head \"downhill\".
Recall that the :ref:`gradient <gradients>` of :math:`J` at :math:`x`
is the vector :math:`\nabla J_x` such that the derivative 
of :math:`J` at :math:`x` in the direction :math:`h`
(i.e. :math:`J_x'(x)`) is given by

.. math:: J_x'(h) = \ip< \nabla J_x, h>_X
  :label: grad1

for all :math:`h\in X`.  The vector :math:`\nabla J_x` points in the 
direction (at :math:`x`) of steepest ascent of :math:`J`.  So the 
direction of steepest descent is :math:`-\nabla J_x`.  To compute this
direction, :ref:`recall <appendix>`

.. math::  

  J'_x(h) &= \left. \frac{d}{dt}\right|_{t=0} J(x+th)\\
  & = \frac{d}{dt}|_{t=0}  \ip<y-\calF(x+th),y-\calF(x+th)>_Y\\
  & =  2 \ip<y-\calF(x), -\left.\frac{d}{dt}\right|_{t=0} \calF(x+th)>_Y\\
  & = -2\ip< y-\calF(x), \calF'_x(h)>_Y.

Now the :ref:`derivative <derivatives>` :math:`\calF'_x` is a linear map from :math:`X` to :math:`Y` and it has an :ref:`adjoint <adjoints>` :math:`(\calF'_x)^*` that satisfies

.. math:: -2\ip< y-\calF(x), \calF'_x(h)>_Y = -2\ip< (\calF'_x)^*(y-\calF(x), h)>_X

for all :math:`h\in X`.  We have just seen that

.. math:: J_x'(h) = -2\ip< (\calF'_x)^*(y-\calF(x)), h>_X
  :label: grad2

Comparing equations :eq:`grad1` and :eq:`grad2` we see that

.. math:: \nabla_x J = -2(\calF'_x)^*(y-\calF(x))

and the direction of steepest descent at :math:`x` is then

.. math:: -\nabla_x J = 2(\calF'_x)^*(y-\calF(x)).

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

The |siple| library contains implementations of the following 
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

Suppose that instead of knowing :math:`y` exactly,
we only know an approximation :math:`y^\delta`
satisfying :math:`||y-y^\delta||_Y<\delta`. Then

.. math:: J(x) = ||y^\delta-\calF(x)||_Y^2.
  :label: Japprox

Optimization methods for finding a minimizer of :math:`J`
typically find a sequence of iterates :math:`x_k` starting
from an initial estimate :math:`x_0`.  The Morozov discrepancy principle
states that iteration should not proceed too far (for otherwise errors
in :math:`y^\delta` will overly pollute the approximate minimizer).
Rather, minimization occurs up to the first iterate :math:`x_k^\delta`
that satisfies

.. math:: J(x_k^\delta) = ||y^\delta-\calF(x_k^\delta)||_Y^2 < \delta^2.
  :label: Morozov

For each error level :math:`\delta` and initial estimate :math:`x_0` we
obtain a regularization strategy :math:`R_\delta(\cdot;x_0)` where

.. math:: R_\delta(y^\delta;x_0) = x_k^\delta

and :math:`x_k^\delta` is the first iterate satisfying inequality
:eq:`Morozov`.  

It has been proved (CITATIONS!) that if :math:`\calF` is a linear map and  :math:`J` is minimized using the Landweber method, or steepest descent, 
or the conjugate gradient method,
then the maps :math:`R_\delta(\cdot;x_0)` are indeed a regularization strategy 
(as defined in the :ref:`overview <regularization>`). However, the formal proofs require stopping as soon as

.. math:: J(x_k^\delta) = ||y^\delta-\calF(x_k^\delta)||_Y^2 < \mu^2\delta^2.
  :label: Morozov_correct

where :math:`\mu>1` is a fixed constant.

We can interpret the Morozov discrepancy principle as follows.  If we only know :math:`y` to within error level :math:`\delta`, then there is no purpose in minimizing :math:`J` beyond the level of the Morozov discrepancy principle, at least not without further justification.  For gradient algorithms, the principle works as a regularization strategy because the downhill
gradient directions tend to make corrections using the smallest possible
changes in :math:`x` (and using changes that are most easily transmitted through :math:`\calF`).  So stopping as soon as possible tends to introduce
the smallest change (with the fewest features) from the initial estimate :math:`x_0` that is needed to find a minimizer that is correct up to the discrepancy :math:`\delta`.

One seeming disadvantage to stopping early is that the final solution
:math:`x_k^\delta` depends on the initial estimate :math:`x_0`.
This phenomenon is a useful feature of the Morozov discrepancy principle, and it allows additional *a-priori* information to be introduced into the solution.  In effect, using the Morozov discrepancy principle finds the
least-featured change to :math:`x_0` that is justified by the error in :math:`y`. As an extreme example, if :math:`x_0` is in fact the exact solution :math:`x` of  :math:`\calF(x)=y`, then 

.. math:: J(x_0) = ||y^\delta-\calF(x_0)||_Y^2 = ||y^\delta-y||_Y^2 < \delta^2

since :math:`||y^\delta-y||_Y<\delta` by hypothesis.  Iterations will stop
immediately with the true solution :math:`x=x_0,` (regardless of the specific value of :math:`\delta`) since the true solution is always consistent to within discrepancy :math:`\delta` of :math:`y^\delta`.


