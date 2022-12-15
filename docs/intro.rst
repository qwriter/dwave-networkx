============
Introduction
============

D-Wave NetworkX provides tools for working with :term:`Chimera` and :term:`Pegasus`
graphs and implementations of graph-theory algorithms on the D-Wave system and other
binary quadratic model :term:`sampler`\ s; for example, functions such as
:meth:`~dwave_networkx.drawing.chimera_layout.draw_chimera`
provide easy visualization for Chimera graphs; functions such
as :meth:`~dwave_networkx.algorithms.max_cut.maximum_cut` 
or :meth:`~dwave_networkx.algorithms.cover.min_vertex_cover` provide graph algorithms useful to
optimization problems that fit well with the D-Wave system.

Like the D-Wave system, all other supported samplers must have
:meth:`~dimod.Sampler.sample_qubo` and :meth:`~dimod.Sampler.sample_ising` methods 
for solving :term:`Ising` and :term:`QUBO` models 
and return an ``iterable`` of samples in order of increasing energy. You can set
a default sampler using the :meth:`~dwave_networkx.default_sampler.set_default_sampler` function.

* For an introduction to quantum processing unit (QPU) topologies such as the
  Chimera and Pegasus graphs, see :std:doc:`Topology <oceandocs:concepts/topology>`.
* For an introduction to binary quadratic models (BQMs), see
  :std:doc:`Binary Quadratic Models <oceandocs:concepts/bqm>`.
* For an introduction to samplers, see
  :std:doc:`Samplers and Composites <oceandocs:concepts/samplers>`.

Example
=======

Below you can see how to create Chimera graphs implemented in the D-Wave 2X and D-Wave 2000Q systems:

.. code:: python

  import dwave_networkx as dnx

  # D-Wave 2X
  C = dnx.chimera_graph(12, 12, 4)

  # D-Wave 2000Q
  C = dnx.chimera_graph(16, 16, 4)
