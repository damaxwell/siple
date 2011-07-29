=============================================
Nonlinear example: coefficient identification
=============================================

Consider the following differential equation on the
interval :math:`[0,1]` with periodic boundary conditions.

.. math:: 
  :label: sturm

  -u''+ e^\beta u &= 1 \\
  u(0) &= u(1) \\
  u'(0) &= u'(1).

.. admonition:: Inverse Problem

  Given a function  :math:`u` defined on an interval :math:`[0,1]`, 
  determine the coefficient :math:`\beta`.
  
The corresponding forward problem is

.. _forward:
.. admonition:: Forward Problem

  Given a coefficient :math:`\beta`, find :math:`u` solving equation :eq:`sturm`.
  
Following our standard notation for forward problems, we'll write :math:`u=\calF(\beta),`  and we'll use the :math:`L^2` norm for functions
defined on :math:`[0,1]`.

To apply a gradient method to solve the inverse problem, we need to compute
the :ref:`linearization <derivatives>` :math:`\calF'_\beta` of :math:`\calF` and its :ref:`adjoint <adjoints>` :math:`(\calF'_\beta)^*`.  Recall

.. math:: \calF'_\beta(h) = \left. \frac{d}{dt}\right|_{t=0} \calF(\beta+th).

Letting :math:`u_t=\calF(\beta+th)` we see that :math:`u_t` is periodic and

.. math:: -u_t''+ e^{\beta+th} u_t &= 1.

Taking the derivative of this equation with respect to :math:`t` and setting
:math:`t=0` we have 

.. math:: -w''+ e^{\beta} w + e^{\beta} h u &= 0.

where :math:`u=u_0` and :math:`w= d/dt|_{t=0} u_t`.  Hence
:math:`\calF'_\beta(h)` is :math:`w` solving

.. math:: 
  :label: sturmp

  -w''+ e^\beta w &= -e^{\beta} u h \\
  w(0) &= w(1) \\
  w'(0) &= w'(1).

To compute the adjoint of :math:`\calF'_\beta`, it's useful
to factor :math:`\calF'_\beta = G \circ H` where
:math:`H(h)=-e^{\beta} u h` and :math:`G(f)=w` where

.. math:: 
  :label: Gdim
  
  -w''+ e^\beta w &= f \\
  w(0) &= w(1) \\
  w'(0) &= w'(1).

Now :math:`(\calF'_\beta)^*=(G\circ H)^*=H^*\circ G^*`.  Now :math:`H^*=H`
since for any functions :math:`k` and :math:`h` in :math:`L^2([0,1])`,

.. math:: \ip<H(h),k> = \int_0^1 H(h) k \;dx = \int_0^1 -e^{\beta}u h k \;dx
= \int_0^1 h(-e^{\beta}u k) \;dx = \int_0^1 h H(k) \;dx = \ip<h,H(k)>.

On the other hand, :math:`G^*=G` also.  To see this, let :math:`f,g\in L^2([0,1])` and let :math:`w=G(f)` and :math:`v=G(g)`.  Then

.. math:: 
  :label: Gadj1

  \int_0^1 v'w' + e^{\beta} wv \; dx = \int_0^1 f v\; dx

since :math:`w` solves equation :eq:`Gdef` and is periodic.  But

.. math:: 
  :label: Gadj2

  \int_0^1 v'w' + e^{\beta} wv \; dx = \int_0^1 w g\; dx

since :math:`v` solves equation :eq:`Gdef` with right-hand side :math:`g`.
From equations :eq:`Gadj1` and :eq:`Gadj2`it follows that

.. math:: \ip<G(f),g> = \int_0^1 w g\; dx = \int_0^1 f v\; dx = \ip<f,G(g)>.

Hence :math:`G^*=G`.  We conclude that 

.. math:: (\calF_\beta')(g) = H^*(G^*(g)) = -e^{\beta u} v

where :math:`v` is periodic and solves

.. math:: -v'' + e^{\beta} v = g.


