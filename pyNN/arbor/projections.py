# encoding: utf-8
"""
Arbor implementation of the PyNN API.

:copyright: Copyright 2006-2016 by the PyNN team, see AUTHORS.
:license: CeCILL, see LICENSE for details.

"""
from copy import deepcopy
import numpy
import logging
try:
    from itertools import izip
except ImportError:
    izip = zip  # Python 3 zip returns an iterator already
from itertools import repeat, chain
from collections import defaultdict
from pyNN import common, errors, core
from pyNN.random import RandomDistribution, NativeRNG
from pyNN.space import Space
from . import simulator
from .standardmodels.synapses import StaticSynapse, TsodyksMarkramSynapse

logger = logging.getLogger("PyNN")

_projections = []  # if a Projection is created but not assigned to a variable,
                   # the connections will not exist, so we store a reference here


class Projection(common.Projection):
    __doc__ = common.Projection.__doc__
    _simulator = simulator
    _static_synapse_class = StaticSynapse

    def __init__(self, presynaptic_population, postsynaptic_population,
                 connector, synapse_type=None, source=None, receptor_type=None,
                 space=Space(), label=None):
        __doc__ = common.Projection.__init__.__doc__
        common.Projection.__init__(self, presynaptic_population, postsynaptic_population,
                                   connector, synapse_type, source, receptor_type,
                                   space, label)
        self._connections = dict((index, defaultdict(list)) for index in self.post._mask_local.nonzero()[0])
        connector.connect(self)
        self._presynaptic_components = dict((index, {}) for index in 
                                            self.pre._mask_local.nonzero()[0])
        if self.synapse_type.presynaptic_type:
            self._configure_presynaptic_components()
        _projections.append(self)
        simulator.model_factory.register(self)
        logger.info("--- Projection[%s].__init__() ---" % self.label)

    @property
    def connections(self):
        for x in self._connections.values():
            for y in x.values():
                for z in y:
                    yield z

    def __getitem__(self, i):
        __doc__ = common.Projection.__getitem__.__doc__
        if isinstance(i, int):
            if i < len(self):
                return self.connections[i]
            else:
                raise IndexError("%d > %d" % (i, len(self) - 1))
        elif isinstance(i, slice):
            if i.stop < len(self):
                return [self.connections[j] for j in range(*i.indices(i.stop))]
            else:
                raise IndexError("%d > %d" % (i.stop, len(self) - 1))

    def __len__(self):
        """Return the number of connections on the local MPI node."""
        return len(list(self.connections))

    def _convergent_connect(self, presynaptic_indices, postsynaptic_index,
                            **connection_parameters):
        """
        Connect a neuron to one or more other neurons with a static connection.

        `presynaptic_cells`     -- a 1D array of pre-synaptic cell IDs
        `postsynaptic_cell`     -- the ID of the post-synaptic cell.
        `connection_parameters` -- each parameter should be either a
                                   1D array of the same length as `sources`, or
                                   a single value.
        """
        pass # TODO: implement _convergent_connect()

    def _configure_presynaptic_components(self):
        """
        For gap junctions potentially other complex synapse types the presynaptic side of the 
        connection also needs to be initiated. This is a little tricky with sources distributed on different nodes as the parameters need to be gathered to the node where the source is hosted before it can be set
        """
        raise NotImplemented

    def _set_attributes(self, parameter_space):
        raise NotImplemented

    def _set_initial_value_array(self, variable, value):
        raise NotImplemented
