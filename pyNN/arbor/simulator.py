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


class _ArborRecipe(arb.recipe):
    """
    Map PyNN model description to Arbor model description.
    
    Specifically, look up all connection and cell type information in
    PyNN data structures and return them in Arbor format.
    """

    def __init__(self, populations, projections):
        self.populations = populations
        self.projections = projections


    def num_cells(self):
        """
        @override   arb::recipe::num_cells
        """
        return sum((pop.size for pop in self.populations))


    def cell_description(self, gid):
        """
        @override   arb::recipe::cell_description
        """
        for pop in self.populations:
            if gid in pop.all_cells:
                return pop.celltype.model # model serves as description
        raise ValueError("Cell GID {} not found in any Population.".format(gid))


    def num_sources(self, gid):
        """
        @override   arb::recipe::num_sources
        """
        pass


    def num_targets(self, gid):
        """
        @override   arb::recipe::num_targets
        """
        pass


    def kind(self, gid):
        """
        @override   arb::recipe::kind
        """
        # TODO: look up gid in Populations and map to kind/celltype
        #   - Q: how is the gid-cell correspondence determined by Arbor?
        for pop in self.populations:
            if gid in pop.all_cells:
                return pop.celltype.arb_cell_kind
        raise ValueError("Cell GID {} not found in any Population.".format(gid))


    def connections_on(self, gid):
        """
        Return connections on cell with given gid.
        
        @override   arb::recipe::conections_on

        Returns:
            list of `arb.connection` objects describing each connection

        """
        pass


class _ArborModelFactory(object):
    """
    Singleton class that keeps track of all PyNN populations and projections
    and creates a corresponding Arbor recipe and Arbor model.
    """

    def __init__(self):
        self.clear()


    def clear(self):
        self.population_list = []
        self.cell_list = []
        self.projection_list = []


    def register(self, *items):
        """
        Add items to the list of cells/populations to be initialized. Cell
        objects must have a `memb_init()` method.

        Parameters
        ----------

        items : *list
            Population or Projection objects
        """
        for item in items:
            if isinstance(item, (common.BasePopulation, common.Assembly)):
                if item.celltype.injectable:  # don't do memb_init() on spike sources
                    self.population_list.append(item)
            elif isinstance(item, common.Projection):
                self.projection_list.append(item)
            else:
                if hasattr(item._cell, "memb_init"):
                    self.cell_list.append(item)


    def build_model(self):
        """
        Finalize recipe before building a runnable model/simulation.
        """
        recipe = _ArborRecipe(self.population_list, self.projection_list) # TODO: build recipe from Populations, Projections, Connections
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
        self.model = model_factory.build_model()
        recorder = arb.make_spike_recorder(model)
        


    def run(self, simtime):
        """
        Advance the simulation for a certain time (quantity).
        """
        self.model.run(500, 0.025)
        # self.run_until(self.tstop + simtime)


    # def run_until(self, tstop):
    #     """
    #     Advance the simulation until certain time (point in time).
    #     """
    #     self._update_current_sources(tstop)
    #     self._pre_run()
    #     self.tstop = tstop
    #     #logger.info("Running the simulation until %g ms" % tstop)
    #     if self.tstop > self.t:
    #         self.parallel_context.psolve(self.tstop)

# --- Initialization, and module attributes ------------------------------------

model_factory = _ArborModelFactory()
state = _State(model_factory)
del _ArborModelFactory
del _State