=======================================
Overview
=======================================

Inverse problems are, as the language suggests, the inverses of forward problems. For us, a forward problem will be a function :math:`\calF` 
from a space :math:`X` to a space :math:`Y`.  If :math:`x\in X`
and 

.. math:: \calF(x) = y
  :label: basicF

then we say that :math:`y` is the solution of the forward problem (for 
:math:`x`).  We will assume that the forward problem is :dfn:`well-posed`:

1. (:dfn:`existence`) For every :math:`x\in X` there is at least one :math:`y\in Y` that solves the forward problem.
2. (:dfn:`uniqueness`) For every :math:`x\in X` there is at most one :math:`y\in Y` that solves the  forward problem.
3. (:dfn:`continuity`) The function :math:`\calF` is continuous.

The first two conditions ensure that for every :math:`x\in X` there is
one and only one solution of the forward problem.  This was already
included in our assumption that :math:`\calF` is a function.  The third
condition guarantees that if :math:`x^\delta` is an approximation of
:math:`x`, then :math:`y^\delta=\calF(x^\delta)` will be an approximation
of :math:`y=\calF(x)`.  This is a more subtle requirement than 
existence or uniqueness, but it is crucial.  Errors are inherent 
to many real problems: measurement error, modeling error, numerical
error, and so forth.  If continuity fails, then the solution 
:math:`y^\delta=\calF(x^\delta)` is not guaranteed
to be close to :math:`y=F(x)` when :math:`x^\delta` is close to :math:`x`.
That is, if :math:`\calF` is not continuous, then error in :math:`x` leads to
uncontrolled error in :math:`y`.

The inverse problem for the forward problem :math:`\calF` is, naturally, given :math:`y\in Y`, find :math:`x\in X` solving equation :eq:`basicF`.  Inverse problems are frequently :dfn:`ill-posed` (i.e. not well-posed), even when the
forward problem is well-posed.

.. admonition:: Exercise

  Suppose the forward problem is the polynomial function
  :math:`\calF(x)=x^3-x` with domain and range both the set of real numbers :math:`\Reals`. Explain why the forward problem is well-posed. Show that the inverse problem is not well-posed (it satisfies existence, but not uniqueness or continuity).

Although a failure of existence or uniqueness causes difficulty,  these
deficiencies can often be treated by expanding the notion of a solution or specification of extra data.  Dealing with a failure of continuity is
more problematic and requires a process known as :dfn:`regularization`.

.. _regularization:

Regularization
^^^^^^^^^^^^^^

Suppose we wish to solve the inverse problem

.. math:: \calF(x) = y
  :label: invreg

where :math:`y` is known and :math:`x` is to be determined.  We assume
that :math:`y` is not known exactly, and only an approximation
:math:`y^\delta` is available.  We also assume that we know something
about how good an approximation :math:`y^\delta` is of :math:`y`,
namely :math:`||y-y^\delta||_Y<\delta`. 

If the inverse problem is not continuous, we do not want to solve

.. math:: \calF(x^\delta) = y^\delta
  :label: invreg_approx

exactly, since small error in :math:`y^\delta` could lead to large
error (with uncontrolled size) in :math:`x^\delta`.  

A regularization strategy for solving equation :eq:`invreg` is
a family of maps :math:`R_\delta : Y\rta X` that are approximate
inverses of :math:`\calF`, one for each error level :math:`\delta`.
Given an approximation :math:`y^\delta` of :math:`y`,

.. math:: x^\delta = R_\delta(y^\delta)

is the associated approximate solution of equation :eq:`invreg`.  It 
typically does not solve equation :eq:`invreg_approx` exactly.  However,
a regularization strategy comes with the following promise

1. For each :math:`\delta` there is a constant :math:`C_\delta` 
   (depending only on the error estimate :math:`\delta`) such that

   .. math:: ||x-R_\delta(y^\delta)||_X < C_\delta.

2. As :math:`\delta\rta 0`, the constants :math:`C_\delta\rta 0`.

Practically speaking, this means that if we use a regularization strategy
:math:`R_\delta` that matches the error in :math:`y^\delta`, then
our approximate solution :math:`x^\delta` will not be too far from :math:`x` (no further than :math:`C_\delta` away), and improvements in our approximations :math:`y^\delta` will lead to improvements in our approximate
solutions :math:`x^\delta`.

The principal regularization algorithms in the |siple| library
employ the so-called Morozov discrepancy principle, as described in
the following section on :ref:`gradient methods <gradient>`.






