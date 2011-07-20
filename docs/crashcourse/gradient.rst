================
Gradient methods
================

The |siple| library focusses on regularization algorithms, and in particular
gradient-based algorithms, which apply in the case that :math:`X` and
:math:`Y` are both inner-product spaces. (The reader might want
to take a moment to brush up on :ref:`inner products, adjoints, and derivatives <appendix>`)

Suppose we wish to solve the inverse problem

.. math:: \calF(x) = y.
  :label: inv

One approach would be to consider the functional

.. math:: J(x) = \frac{1}{2}||y-\calF(x)||_Y^2 = \frac{1}{2}\ip<y-\calF(x),y-\calF(x)>_Y.

A solution of equation :eq:`inv` corresponds to a minimizer of 
:math:`J`  (since :math:`J` is always non-negative, and equals
zero exactly when equation :eq:`inv` holds).  




  The gradient of :math:`J`

at :math:`x` is the element :math:`\nabla J_x \in X` such that

.. math:: J'_x(h) = \ip<\nabla J_x, h>_X.

.. math:: 

  J'_x(h)=\ip<y-\calF(x), \calF'_x(h)>_Y = 
  \ip<(\calF'_x)^*(y-\calF(x)),h>_X.

for all :math:`h\in X`.

