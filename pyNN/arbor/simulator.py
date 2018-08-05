# encoding: utf8
"""
Implementation of the "low-level" functionality used by the common
implementation of the API, for the Arbor simulator.

Classes and attributes useable by the common implementation:

Classes:
    ID
    Connection

Attributes:
    state -- a singleton instance of the _State class.

All other functions and classes are private, and should not be used by other
modules.

:copyright: Copyright 2006-2016 by the PyNN team, see AUTHORS.
:license: CeCILL, see LICENSE for details.

"""

try:
    xrange
except NameError:
    xrange = range
from pyNN import __path__ as pyNN_path
from pyNN import common
import logging
import numpy
import os.path
from operator import itemgetter
import pyarb as arb

logger = logging.getLogger("PyNN")
name = "ARBOR"  # for use in annotating output data

# Instead of starting the projection var-GID range from 0, the first _MIN_PROJECTION_VARGID are 
# reserved for other potential uses
_MIN_PROJECTION_VARGID = 1000000


class _Recipe(arb.recipe):
    """
    Map PyNN model description to Arbor model description.
    
    Specifically, look up all connection and cell type information in
    PyNN data structures and return them in Arbor format.
    """

    def __init__(self, TODO):
        # TODO: pass in Populations, Projections or equivalent data structures
        pass

    def num_cells(self):
        pass

    def cell_description(self, gid):
        # TODO: cell description should be created by the standardcelltype class
        pass

    def num_sources(self, gid):
        pass

    def num_targets(self, gid):
        pass

    def kind(self, gid):
        # TODO: look up gid in Populations and map to kind/celltype
        pass

    def connections_on(self, gid):
        """
        Return connections on cell with given gid.

        Returns:
            list of `arb.connection` objects describing each connection

        """
        pass


class _ModelFactory(object):
    """
    Singleton class that keeps track of all PyNN populations and projections
    and creates a corresponding Arbor recipe and Arbor model.
    """

    def __init__(self):
        self.clear()


    def clear(self):
        self.population_list = []
        self.cell_list = []


    def register(self, *items):
        """
        Add items to the list of cells/populations to be initialized. Cell
        objects must have a `memb_init()` method.
        """
        for item in items:
            if isinstance(item, (common.BasePopulation, common.Assembly)):
                if item.celltype.injectable:  # don't do memb_init() on spike sources
                    self.population_list.append(item)
            else:
                if hasattr(item._cell, "memb_init"):
                    self.cell_list.append(item)


    def build_model(self):
        """
        Finalize recipe before building a runnable model/simulation.
        """
        recipe = _Recipe(TODO) # TODO: build recipe from Populations, Projections, Connections
        decomp = arb.partition_load_balance(recipe)
        model = arb.model(recipe, decomp)
        return model




class _State(common.control.BaseState):
    """
    Represent the simulator state.
    """

    def __init__(self):
        """
        Initialize the simulator.
        """
        pass


    def _pre_run(self):
        """
        Finalize the recipe.
        """
        model = model_factory.build_model()
        recorder = arb.make_spike_recorder(model)
        model.run(500, 0.025)


    def run(self, simtime):
        """
        Advance the simulation for a certain time (quantity).
        """
        self.run_until(self.tstop + simtime)


    def run_until(self, tstop):
        """
        Advance the simulation until certain time (point in time).
        """
        self._update_current_sources(tstop)
        self._pre_run()
        self.tstop = tstop
        #logger.info("Running the simulation until %g ms" % tstop)
        if self.tstop > self.t:
            self.parallel_context.psolve(self.tstop)

# --- Initialization, and module attributes ------------------------------------

model_factory = _ModelFactory()
state = _State(model_factory)
del _ModelFactory
del _State