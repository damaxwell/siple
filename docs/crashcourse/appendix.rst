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

.. _adjoints:

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

.. _derivatives:

Derivatives
^^^^^^^^^^^^^^

Suppose :math:`\calF` is a (possibly nonlinear) map from the vector space :math:`X` to the vector space :math:`Y`.  Given a point :math:`x\in X` and a direction :math:`h\in X`,
the :dfn:`directional derivative` of :math:`\calF` at :math:`x` in the direction :math:`h` is

.. math:: \calF'_x(h) = \left. \frac{d}{dt} \right|_{t=0} \calF(x+th).

For example, if :math:`\calF(x_1,x_2)=\cos(x_1 x_2^2)` 
then the directional derivative of :math:`\calF` at :math:`x=(5,3)`
in the direction :math:`h=(0,1)` is 

.. math:: \left. \frac{d}{dt} \right|_{t=0} \cos( 5\cdot(3+t)^2) = \left.-\sin(5\cdot(3+t)^2)\cdot 5\cdot 2\cdot(3+t)\right|_{t=0} =  -30 \sin(45).

It would perhaps be more familiar to write this as 

.. math:: \left.\frac{\partial}{\partial x_2} \calF(x_1,x_2)\right|_{x=(2,3)}

but the more general notation is helpful when :math:`X` doesn't have a preferred basis like :math:`\Reals^2` does.  

The map :math:`\calF'_x` taking :math:`h` to the directional derivative :math:`\calF'_x(h)` is called the :dfn:`total derivative` of :math:`\calF`
at :math:`x`, and is sometimes called the :dfn:`linearization` of 
:math:`\calF` at :math:`x`.  The linearization of :math:`\calF` at :math:`x`
is, as its name suggests, a linear map from :math:`X` to :math:`Y`.  For small values of :math:`h`,

.. math:: \calF(x+h) \approx \calF(x) + \calF'_x(h)
  :label: linapprox

and :math:`\calF'_x` is the linear map that works best in this approximation.

.. admonition:: Exercise
  
  Let :math:`\calF(x)=\sin(x)`.  In a first-year calculus class you would
  have written :math:`\calF'(x)=\cos(x)`.  In terms of the notation we
  are using here, this means that the derivative of :math:`\calF` at :math:`x` is the linear map :math:`\calF'_x` that takes :math:`h` to :math:`\cos(x)h`.
  Notice that this map depends nonlinearly on :math:`x` but linearly on :math:`h`.  
  
  With all this in mind, use :math:`x=0` and :math:`h=0.2` in equation
  :eq:`linapprox` to approximate :math:`\sin(0.2)`.  How good is the approximation?

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


.. _gradients:

Gradients
^^^^^^^^^^^^^^

Let :math:`X` be an inner product space, let :math:`z\in X`, and
define

.. math:: T_z(x) = \ip<z,x>_X.

Then :math:`T_z` is a linear map from :math:`X` to :math:`\Reals`,
For reasonable inner product spaces, if :math:`T:X\rta\Reals`
is a linear map, then there exists a :math:`z\in X` such that

.. math:: T = T_z.

That is, linear maps from :math:`X` to :math:`\Reals` can be 
written in terms of inner products.

Now suppose :math:`\calF:X\rta  \Reals` is some nonlinear function.
Then the derivative of :math:`\calF` at :math:`x` is the linear
map :math:`\calF_x'` from :math:`X` to :math:`\Reals.`  Given
what we have just discussed, there is a vector :math:`z\in X`
such that

.. math:: \calF'_x(h) = \ip<z,h>_X

for all :math:`h\in X`.  We call :math:`z` the :dfn:`gradient` of
:math:`\calF` at :math:`x` and write it as :math:`\nabla \calF_x`.

It's important to keep the notation straight: :math:`\calF'_x` is a linear map from :math:`X` to :math:`\Reals`, and :math:`\nabla \calF_x` is a vector in :math:`X`.  They are related by

.. math:: \calF'_x(h) = \ip<\nabla F_x,h>_X

for all :math:`h\in X`.  Now

.. math:: \ip<\nabla F_x,h> = ||\nabla F_x|| ||h|| \cos(\theta)

where :math:`\theta` is the angle between the vectors :math:`\nabla F_x`
and :math:`h`.  If :math:`h` is a unit vector, this equation becomes

.. math:: \ip<\nabla F_x,h> = ||\nabla F_x|| \cos(\theta)

and it is maximized exactly when :math:`\theta=0` (i.e. when :math:`h` points in the direction of :math:`\nabla F_x`).  So :math:`\nabla F_x` points
in the direction of steepest increase of :math:`\calF`  at :math:`x`, and
its length encodes how fast :math:`\calF` is changing in this direction.

