================================
Inverse Problems Overview
================================

Inverse problems are, as the language suggests, the inverses of forward problems. For us, a forward problem will be a function :math:`\calF` 
from a space :math:`X` to a space :math:`Y`, and we will insist that this function be continuous.  That is, :math:`X` and :math:`Y` will be equipped with a notion of distance, and if :math:`\hat x` is brought near :math:`x`, then :math:`\hat y=\calF(\hat x)` can be made as close to :math:`y=\calF(x)` as we wish.

For example, let :math:`\rho` be a density of matter in space :math:`\Reals^3` that
is contained inside some ball of radius :math:`R` about the origin. The forward 
problem is to find the gravitational potential of :math:`\rho`.  That is, find
a function :math:`\Phi=\calF(\rho)` solving

.. math::
    -\Laplacian \Phi = 4\pi G \rho

where :math:`G` is the universal gravitational constant.  For the problem
to be well-defined, we need to add a boundary condition, and we
require the decay condition that :math:`\Phi\rta 0` at infinity.  

Showing this  particular problem is continuous requires some technical details: a good choice of distance is needed for the spaces :math:`X` and :math:`Y` of matter distributions and gravitational potentials (so-called weighted :math:`L^2` spaces would be reasonable).  But the idea that small perturbations in :math:`\rho` make correspondingly small perturbations in :math:`\Phi` is clear.  The moon doesn't care when I walk around on the Earth.

The corresponding inverse problem is as follows.  Given a potential field :math:`\Phi`, determine :math:`\rho`.  This problem has a very different character than the original problem.  If you perturb :math:`\Phi` in a way so that its second derivatives don't change much, then you won't change :math:`\rho` much.  
But measurements of physical fields frequently measure their values but not their derivatives.  A measurement :math:`\hat \Phi` of :math:`\Phi` might have values that are close (in some sense) to :math:`\Phi` but give no control on its derivatives.
For example, :math:`\hat \Phi` could be obtained from :math:`\Phi` by adding a tiny but very high-frequency oscillation.  This difficulty is a typical feature of a true inverse problem.

In general terms, the inverse problem for :math:`\calF` is as follows.  Given :math:`y\in Y`, find :math:`x\in X` such that

.. math::
    \calF(x) = y.
    :label: inverse

The theory of inverse problems focusses on problems where the inverse problem is ill-posed.  That is, at least one of the following fails:


#. Given any :math:`y\in Y`, there is a solution :math:`x\in X` of equation :eq:`inverse`.
#. Given any :math:`y\in Y`, there is no more than one solution :math:`x\in X` of equation :eq:`inverse`.
#. The solution of :eq:`inverse` depends continuously on on :math:`y\in Y`.

The difficulties posed by a failure of requirements 1) and 2) are not the focus of ``siple``, and can sometimes be dealt with by expanding the notion of a solution or specification of additional constraints.  This was done, for example, in the forward problem for the gravitational potential where the constraint that :math:`\Phi\rta 0` at infinity was added.

Requirement 3) is more subtle and shows up frequently in cases (like the density problem) where the forward problem involves solving a differential equation, and the inverse problem involves taking a derivative.  It is 
typically dealt with by so-called regularization methods.

Regularization
^^^^^^^^^^^^^^

Suppose :math:`\calF(x)=y` and :math:`y^\delta` is an approximation of :math:`y`.


The derivative problem
----------------------

We will work with the following model problem.  Given measurements
of a function :math:`f` defined on a circle, determine its derivative.

