.. _appendix:

======================================================
Appendix
======================================================


Inner products
^^^^^^^^^^^^^^

An :dfn:`inner product` on a vector space :math:`X` is a kind of
generalization of the Euclidean dot product on :math:`\Reals^n`.
For example, if :math:`f` and :math:`g` are functions defined on the interval :math:`[0,1]`, their :math:`L^2` inner product is

.. math:: \ip<f,g>_{L^2} = \int_0^1 f(x)g(x)\; dx.

If :math:`x` belongs to the vector space :math:`X`, and if
:math:`X` is equipped with the inner product :math:`\ip<\cdot,\cdot>_X`, 
we define the :dfn:`norm`  of :math:`x` to be 

.. math:: ||x||_{X} = \sqrt{\ip<x,x>_X}.

If :math:`x` and :math:`z` belong to :math:`X`, the distance
between :math:`x` and :math:`z` is

.. math:: ||x-z||_X.

All of these formulas are analogues of the familiar Euclidean ones.

A vector space equipped with an inner product is called an
inner product space.

Adjoints
^^^^^^^^

Gradient-based regularization methods make use of a mathematical
object known as an adjoint. Suppose :math:`T:X\rta Y` is a linear map between two inner product spaces. Its :dfn:`adjoint`  
is the linear map :math:`T^*` from :math:`Y` to :math:`X` that
satisfies the equation

.. math:: \ip<T(x),y>_{Y} = \ip<x,T^*(y)>_{X}

for all :math:`x\in X` and :math:`y\in Y`.  

.. admonition:: Exercise

  Suppose that :math:`X=Y=\Reals^3`, and that these spaces are equipped
  with the standard Euclidean dot product (so :math:`\ip<x,z>=x\cdot z`).
  Let :math:`T` be the map 
  
  .. math:: T(x) = \begin{pmatrix} 1 & 3 &0\\ 
                     -1& 1 & 0\\
                     7 & 0 & 1 \end{pmatrix} \begin{pmatrix} x_1\\x_2\\x_3\end{pmatrix}.

  What is :math:`T^*`?  (Hint: write the dot product :math:`x\cdot y` 
  using a transpose as :math:`x^T y`)

Every linear map between inner product spaces has an adjoint, but 
computation of the adjoint of a particular linear map is a nontrivial task.
For maps involving differential equations, it often involves integration by 
parts formulas, and there are examples of these in the :ref:`tutorial`.


Derivatives
^^^^^^^^^^^^^^

Consider the nonlinear map :math:`\calF:\Reals^3\rta\Reals^2`

.. math:: \calF(x_1,x_2,x_3) = (x_1 x_2, x_1-x_3).

Its Jacobian matrix of partial derivatives :math:`\partial\calF_i/\partial x_j` is

.. math:: T_x =  \begin{pmatrix} x_2 & x_1 & 0 \\
                                 1 & 0 & -1 \end{pmatrix}.

We can consider :math:`T_x` to be a linear map from :math:`\Reals^3` to
:math:`\Reals^2`, and for small values of :math:`h\in\Reals^3`,

.. math:: \calF(x+h) \approx \calF(x) + T_x(h).

In fact, :math:`T_x` is the linear map that works best in this approximation.
We call :math:`T_x` the derivative of :math:`\calF` at :math:`x`,
and write it as :math:`\calF'_x`.

One way to compute :math:`\calF'_x(h)` is via the formula

.. math:: \calF'_x(h) = \left.\frac{d}{dt}\right|_{t=0} \calF(x+th).

The quantity :math:`\calF'_x(h)` is sometimes called the :dfn:`directional derivative` of :math:`\calF` at :math:`x` in the direction :math:`h`. 

As a non-trivial example of computing a derivative, let :math:`\calF`
be the map taking the function :math:`f` defined on :math:`[0,1]` to
the solution :math:`u` of 

.. math::
  :label: sturm
  
   -u'' + u^3 &= f\\
   u(0) &= 1\\
   u'(1) &= 5.

That is, :math:`\calF(f)=u` where :math:`u` solves :eq:`sturm`.
To compute :math:`\calF'_f(h)` we consider for each :math:`t\in \Reals`

.. math:: u_t = \calF(f+t h)

so

.. math:: \calF_f'(h) = \left.\frac{d}{dt}\right|_{t=0} \calF(f+t h) = \left.\frac{d}{dt}\right|_{t=0} u_t.

To keep the notation tidy we'll write :math:`u_0 = u` and :math:`\left.\frac{d}{dt}\right|_{t=0} u_t = w`.
For each :math:`t`, :math:`u_t` satisfies

.. math:: 

  -u_t'' + (u_t)^3 &= f+th \\
   u_t(0) &= 1 \\
   u_t'(1) &= 5.

Taking the derivatives of these equations with respect to :math:`t` 
and setting :math:`t=0` we see

.. math:: 
  :label: linearized

  -w'' + 3u^2w &= h\\
   w(0) &= 0\\
   w'(1) &= 0
  
Unraveling everything we have just done, we see that :math:`\calF'_f(h)`
is the solution :math:`w` of equations :eq:`linearized`, where :math:`u=\calF(f)`.


Gradients
^^^^^^^^^^^^^^

Suppose :math:`J` is a map from an inner product space :math:`X` to :math:`\Reals`.  Then for each :math:`x\in X`, the derivative :math:`J'_x` is a  linear map from :math:`X` to :math:`\Reals`, and

.. math:: J(x+h) \approx J(x) + J'_x(h)

for small vectors :math:`h\in X`.  A linear map from :math:`X` to 
:math:`\Reals` can usually be written as an inner product with some
vector :math:`z` in :math:`X`.  That is, there is a vector :math:`z\in X`
such that

.. math:: J'_x(h) = \ip<z,h>_X

for all :math:`h\in X`.  We call :math:`z` the gradient of :math:`J` at :math:`x` and write it as :math:`\nabla J_x`.  






