================
Gradient methods
================

The |siple| library focusses on regularization algorithms, and in particular
gradient-based algorithms, which apply in the case that :math:`X` and
:math:`Y` are both inner product spaces. (The reader might want
to take a moment to brush up on :ref:`inner products, adjoints, and derivatives <appendix>`.)

Suppose we wish to solve the inverse problem

.. math:: \calF(x) = y.
  :label: inv

One approach would be to consider the functional

.. math:: J(x) = \frac{1}{2}||y-\calF(x)||_Y^2 = \frac{1}{2}\ip<y-\calF(x),y-\calF(x)>_Y.

A solution of equation :eq:`inv` corresponds to a minimizer of 
:math:`J`  (since :math:`J` is always non-negative, and equals
zero exactly when equation :eq:`inv` holds).  So the problem 
of finding a solution of :eq:`inv` corresponds to minimizing :math:`J`,
and it turns out that there are good techniques for finding approximate
minimizers of :math:`J` that regularize the inverse problem.

  A natural approach to minimizing :math:`J` is to head \"downhill\".
Recall that the gradient of :math:`J` at :math:`x`
is the vector :math:`\nabla J_x` such that the derivative
of :math:`J` at :math:`x` in the direction :math:`h`
(i.e. :math:`J_x'(x)`) is given by

.. math:: J_x'(h) = \ip< \nabla J_x, h>_X

for all :math:`h\in X`.






  The gradient of :math:`J`

at :math:`x` is the element :math:`\nabla J_x \in X` such that

.. math:: J'_x(h) = \ip<\nabla J_x, h>_X.

.. math:: 

  J'_x(h)=\ip<y-\calF(x), \calF'_x(h)>_Y = 
  \ip<(\calF'_x)^*(y-\calF(x)),h>_X.

for all :math:`h\in X`.

