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

When a vector space is equipped with an inner-product, we'll call it
an inner-product space.

Adjoints
^^^^^^^^

Gradient-based regularization methods make use of a mathematical
object known as an adjoint. Suppose :math:`T` is a linear map between the inner-product spaces :math:`X` and :math:`Y`. Its :dfn:`adjoint`  
is the linear map :math:`T^*` from :math:`Y` to :math:`X` that
satisfies the equation

.. math:: \ip<T(x),y>_{Y} = \ip<x,T^*(y)>_{X}

for all :math:`x\in X` and :math:`y\in Y`.  Every linear map
has an adjoint.

.. admonition:: Exercise

  Suppose that :math:`X=Y=\Reals^3`, and that these spaces are equipped
  with the standard Euclidean dot product (so :math:`\ip<x,z>=x\cdot z`).
  Let :math:`T` be the map 
  
  .. math:: T(x) = \begin{pmatrix} 1 & 3 &0\\ 
                     -1& 1 & 0\\
                     7 & 0 & 1 \end{pmatrix} \begin{pmatrix} x_1\\x_2\\x_3\end{pmatrix}.

  What is :math:`T^*`?

Derivatives
^^^^^^^^^^^^^^

Suppose :math:`\calF` is a map from the inner-product space :math:`X` to the
inner-product space :math:`Y`.  The derivative of :math:`\calF` at :math:`x`
is a linear map :math:`T` such that

.. math:: \calF(x+h) \approx \calF(x) + T(h)

for small vectors :math:`h`. The approximation is so good that the ratio

.. math:: \frac{||\calF(x+h)-[\calF(x)+T(h)]||_Y}{||h||_X}

tends to zero as the norm of :math:`h` goes to zero.  We'll
write :math:`\calF'_x` to denote the derivative of :math:`\calF` at
:math:`x`.

For example, consider :math:`\calF(x)=x^2`.  In a first-year
calculus class, you would have written :math:`\calF'(x)=2x`.
This is very compact notation for the following: at each :math:`x\in\Reals`,
the derivative of :math:`\calF` at :math:`x` is the linear map :math:`T_x` where :math:`T_x(h)=(2x)(h)`.  Note that the linear map depends on :math:`x`,
so the derivative of :math:`\calF` at :math:`x=1`
is the linear map :math:`T(h)=2h`. For small values of :math:`h`,
:math:`\calF(1+h)\approx \calF(1)+2h = 1+2h`. The compact notation works in one dimension because there is a one-to-one correspondence between linear maps and numbers.  This breaks down in higher dimensions where linear maps might correspond to matrices.
