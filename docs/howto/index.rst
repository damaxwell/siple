.. _user:

|siple| HOWTO
========================

Linear Algebra
______________


Parameter Mechanism
____________________________


Gradient Algorithms
____________________________

The inverse problem to solve is 

.. math:: \calF(x) = y

where :math:`y` is known and :math:`x` is to be found.

Linear case
^^^^^^^^^^^

Suppose :math:`\calF` is a linear map (which we'll call :math:`T`
here and in the |siple| code as a reminder that it is linear). 
We are to solve

.. math:: T(x) = y.

The forward problem is specified with a :class:`LinearFowardProblem <siple.gradient.forward.LinearForwardProblem>`, which implements
the inner products on :math:`X` and :math:`Y`, as well as 
the forward and adjoint maps :math:`T` and :math:`T^*`.
