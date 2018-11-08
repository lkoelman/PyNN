# encoding: utf-8
"""
Definition of cell models for the Arbor module.

These are the cell model classes associated with a cell type.
They represent the implementation aspect of the cell type, whereas the 
celltype class represents the specification.

:copyright: Copyright 2006-2018 by the PyNN team, see AUTHORS.
:license: CeCILL, see LICENSE for details.

"""

import logging
from math import pi
import numpy
import pyarb as arb

from pyNN import errors
from pyNN.models import BaseCellType
from .recording import recordable_pattern
from .simulator import state

try:
    reduce
except NameError:
    from functools import reduce

logger = logging.getLogger("PyNN")


# TODO: inherit from the Arbor cell description classes
# - in Arbor master they are in include/arbor/[lif_cell/mc_cell/spike_source_cell].hpp
# - in Python branch they are in src/[cell/rss_cell/lif_cell_description].hpp
#   - in pyarb.cpp the three cell descriptions (classes) exported are:
#       + lif_cell_description (for LIF cells)
#       + rss_cell (for periodic spike trains)
#       + cell (for multi-compartmental cells)
#   - in pyarb.cpp the three cell kinds (enum instances) exported are:
#       + arb.cell_kind.cable1d_neuron
#       + arb.cell_kind.lif_neuron
#       + arb.cell_kind.data_spike_source
#       + arb.cell_kind.regular_spike_source


def _new_property(obj_hierarchy, attr_name):
    """
    Returns a new property, mapping attr_name to obj_hierarchy.attr_name.

    For example, suppose that an object of class A has an attribute b which
    itself has an attribute c which itself has an attribute d. Then placing
      e = _new_property('b.c', 'd')
    in the class definition of A makes A.e an alias for A.b.c.d

    Parameters
    ----------

    obj_hierarchy : str
        Object hierarcy separated by dots, e.g. 'A.b.c.d'
    """

    def set(self, value):
        obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
        setattr(obj, attr_name, value)

    def get(self):
        obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
        return getattr(obj, attr_name)
    return property(fset=set, fget=get)


################################################################################
# Biophysical Compartmental Cell Models
################################################################################

# TODO: implement multi-compartmental neurons.
#       - see PyNN MC feature branch for implementation ideas
#       - combine with existing Ephys-based codebase
#       - Morphologies in Arbor:
#           + see example https://github.com/eth-cscs/arbor/blob/master/example/miniapp/miniapp_recipes.cpp

# NOTE: in arbor master branch 'cell' has become 'mc_cell'

class BaseCompartmentalNeuron(arb.cell):
    """docstring"""

    def __init__(self, c_m, i_offset):

        # initialise Section object with 'pas' mechanism
        diam = 1000 / pi  # gives area = 1e-3 cm2
        self.seg = soma = arb.cell(diam)
        soma.add_mechanism('hh')

        # TODO: find out how to set all passive properties

        self.source_section = self.seg # TODO: source section only necessary for NEURON backend?

        # TODO: insert current source

        # TODO: setup for recording

        self.v_init = None

    def area(self):
        """Membrane area in µm²"""
        return pi * self.L * self.seg.diam

    c_m = _new_property('seg', 'cm')
    i_offset = _new_property('stim', 'amp')

    def memb_init(self):
        assert self.v_init is not None, "cell is a %s" % self.__class__.__name__
        for seg in self:
            seg.v = self.v_init
        #self.seg.v = self.v_init

    def set_parameters(self, param_dict):
        for name in self.parameter_names:
            setattr(self, name, param_dict[name])


class SingleCompartmentNeuron(arb.cell):
    pass


class MultiCompartmentNeuron(arb.cell):
    # see include/arbor/mc_cell.hpp
    pass

################################################################################
# Artificial Cell Models
################################################################################

# Integrate-and-fire type cell models and mathematical reductions

class BaseIntegrateFireNeuron(arb.lif_cell_description):
    """Single compartment with excitatory and inhibitory synapses"""

    synapse_models = {
        'current': {'exp': h.ExpISyn, 'alpha': h.AlphaISyn},
        'conductance': {'exp': h.ExpSyn, 'alpha': h.AlphaSyn},
    }

    def __init__(self, syn_type, syn_shape, c_m, i_offset,
                 tau_e, tau_i, e_e, e_i):
        BaseSingleCompartmentNeuron.__init__(self, c_m, i_offset)

        self.syn_type = syn_type
        self.syn_shape = syn_shape

        # insert synapses
        assert syn_type in ('current', 'conductance'), "syn_type must be either 'current' or 'conductance'. Actual value is %s" % syn_type
        assert syn_shape in ('alpha', 'exp'), "syn_type must be either 'alpha' or 'exp'"
        synapse_model = self.synapse_models[syn_type][syn_shape]
        self.esyn = synapse_model(0.5, sec=self)
        self.isyn = synapse_model(0.5, sec=self)

    def area(self):
        """Membrane area in µm²"""
        return pi * self.L * self.seg.diam

    c_m = _new_property('seg', 'cm')
    i_offset = _new_property('stim', 'amp')

    def memb_init(self):
        assert self.v_init is not None, "cell is a %s" % self.__class__.__name__
        for seg in self:
            seg.v = self.v_init
        #self.seg.v = self.v_init

    def set_parameters(self, param_dict):
        for name in self.parameter_names:
            setattr(self, name, param_dict[name])

    @property
    def excitatory(self):
        return self.esyn

    @property
    def inhibitory(self):
        return self.isyn

    def _get_tau_e(self):
        return self.esyn.tau

    def _set_tau_e(self, value):
        self.esyn.tau = value
    tau_e = property(fget=_get_tau_e, fset=_set_tau_e)

    def _get_tau_i(self):
        return self.isyn.tau

    def _set_tau_i(self, value):
        self.isyn.tau = value
    tau_i = property(fget=_get_tau_i, fset=_set_tau_i)

    def _get_e_e(self):
        return self.esyn.e

    def _set_e_e(self, value):
        self.esyn.e = value
    e_e = property(fget=_get_e_e, fset=_set_e_e)

    def _get_e_i(self):
        return self.isyn.e

    def _set_e_i(self, value):
        self.isyn.e = value
    e_i = property(fget=_get_e_i, fset=_set_e_i)


class LeakySingleCompartmentNeuron(SingleCompartmentNeuron):

    def __init__(self, syn_type, syn_shape, tau_m, c_m, v_rest, i_offset,
                 tau_e, tau_i, e_e, e_i):
        SingleCompartmentNeuron.__init__(self, syn_type, syn_shape, c_m, i_offset,
                                         tau_e, tau_i, e_e, e_i)
        self.insert('pas')
        self.v_init = v_rest  # default value

    def __set_tau_m(self, value):
        #print("setting tau_m to", value, "cm =", self.seg.cm))
        self.seg.pas.g = 1e-3 * self.seg.cm / value  # cm(nF)/tau_m(ms) = G(uS) = 1e-6G(S). Divide by area (1e-3) to get factor of 1e-3

    def __get_tau_m(self):
        #print("tau_m = ", 1e-3*self.seg.cm/self.seg.pas.g, "cm = ", self.seg.cm)
        return 1e-3 * self.seg.cm / self.seg.pas.g

    def __get_cm(self):
        #print("cm = ", self.seg.cm)
        return self.seg.cm

    def __set_cm(self, value):  # when we set cm, need to change g to maintain the same value of tau_m
        #print("setting cm to", value)
        tau_m = self.tau_m
        self.seg.cm = value
        self.tau_m = tau_m

    v_rest = _new_property('seg.pas', 'e')
    tau_m = property(fget=__get_tau_m, fset=__set_tau_m)
    c_m = property(fget=__get_cm, fset=__set_cm)  # if the property were called 'cm'
                                                    # it would never get accessed as the
                                                    # built-in Section.cm would always
                                                    # be used first


class StandardIF(LeakySingleCompartmentNeuron):
    """docstring"""

    def __init__(self, syn_type, syn_shape, tau_m=20, c_m=1.0, v_rest=-65,
                 v_thresh=-55, t_refrac=2, i_offset=0, v_reset=None,
                 tau_e=5, tau_i=5, e_e=0, e_i=-70):
        if v_reset is None:
            v_reset = v_rest
        LeakySingleCompartmentNeuron.__init__(self, syn_type, syn_shape, tau_m, c_m, v_rest,
                                              i_offset, tau_e, tau_i, e_e, e_i)
        # insert spike reset mechanism
        self.spike_reset = h.ResetRefrac(0.5, sec=self)
        self.spike_reset.vspike = 40  # (mV) spike height
        self.source = self.spike_reset
        self.rec = h.NetCon(self.source, None)

        # process arguments
        self.parameter_names = ['c_m', 'tau_m', 'v_rest', 'v_thresh', 't_refrac',   # 'c_m' must come before 'tau_m'
                                'i_offset', 'v_reset', 'tau_e', 'tau_i']
        if syn_type == 'conductance':
            self.parameter_names.extend(['e_e', 'e_i'])
        self.set_parameters(locals())

    v_thresh = _new_property('spike_reset', 'vthresh')
    v_reset = _new_property('spike_reset', 'vreset')
    t_refrac = _new_property('spike_reset', 'trefrac')


class Izhikevich_(BaseSingleCompartmentNeuron):
    """docstring"""

    def __init__(self, a_=0.02, b=0.2, c=-65.0, d=2.0, i_offset=0.0):
        BaseSingleCompartmentNeuron.__init__(self, 1.0, i_offset)
        self.L = 10
        self.seg.diam = 10 / pi
        self.c_m = 1.0

        # insert Izhikevich mechanism
        self.izh = h.Izhikevich(0.5, sec=self)
        self.source = self.izh
        self.rec = h.NetCon(self.seg._ref_v, None,
                            self.get_threshold(), 0.0, 0.0,
                            sec=self)
        self.excitatory = self.inhibitory = self.source

        self.parameter_names = ['a_', 'b', 'c', 'd', 'i_offset']
        self.set_parameters(locals())
        self.u_init = None

    a_ = _new_property('izh', 'a')
    b = _new_property('izh', 'b')
    c = _new_property('izh', 'c')
    d = _new_property('izh', 'd')
    ## using 'a_' because for some reason, cell.a gives the error "NameError: a, the mechanism does not exist at PySec_170bb70(0.5)"

    def get_threshold(self):
        return self.izh.vthresh

    def memb_init(self):
        assert self.v_init is not None, "cell is a %s" % self.__class__.__name__
        assert self.u_init is not None
        for seg in self:
            seg.v = self.v_init
        self.izh.u = self.u_init


################################################################################
# Artificial Spike Sources
################################################################################


class RandomSpikeSource(hclass(h.NetStimFD)):

    parameter_names = ('start', '_interval', 'duration')

    def __init__(self, start=0, _interval=1e12, duration=0):
        self.start = start
        self.interval = _interval
        self.duration = duration
        self.noise = 1
        self.spike_times = h.Vector(0)
        self.source = self
        self.rec = h.NetCon(self, None)
        self.switch = h.NetCon(None, self)
        self.source_section = None
        self.seed(state.mpi_rank + state.native_rng_baseseed)  # should allow user to set specific seeds somewhere, e.g. in setup()

    def _set_interval(self, value):
        self.switch.weight[0] = -1
        self.switch.event(h.t + 1e-12, 0)
        self.interval = value
        self.switch.weight[0] = 1
        self.switch.event(h.t + 2e-12, 1)

    def _get_interval(self):
        return self.interval
    _interval = property(fget=_get_interval, fset=_set_interval)


class RandomPoissonRefractorySpikeSource(hclass(h.PoissonStimRefractory)):

    parameter_names = ('rate', 'tau_refrac', 'start', 'duration')

    def __init__(self, rate=1, tau_refrac=0.0, start=0, duration=0):
        self.rate = rate
        self.tau_refrac = tau_refrac
        self.start = start
        self.duration = duration
        self.spike_times = h.Vector(0)
        self.source = self
        self.rec = h.NetCon(self, None)
        self.source_section = None
        self.seed(state.mpi_rank + state.native_rng_baseseed)


class RandomGammaSpikeSource(hclass(h.GammaStim)):

    parameter_names = ('alpha', 'beta', 'start', 'duration')

    def __init__(self, alpha=1, beta=0.1, start=0, duration=0):
        self.alpha = alpha
        self.beta = beta
        self.start = start
        self.duration = duration
        self.spike_times = h.Vector(0)
        self.source = self
        self.rec = h.NetCon(self, None)
        self.switch = h.NetCon(None, self)
        self.source_section = None
        self.seed(state.mpi_rank + state.native_rng_baseseed)


class VectorSpikeSource(hclass(h.VecStim)):

    parameter_names = ('spike_times',)

    def __init__(self, spike_times=[]):
        self.spike_times = spike_times
        self.source = self
        self.source_section = None
        self.rec = None

    def _set_spike_times(self, spike_times):
        # spike_times should be a Sequence object
        try:
            self._spike_times = h.Vector(spike_times.value)
        except (RuntimeError, AttributeError):
            raise errors.InvalidParameterValueError("spike_times must be an array of floats")
        if numpy.any(spike_times.value[:-1] > spike_times.value[1:]):
            raise errors.InvalidParameterValueError("Spike times given to SpikeSourceArray must be in increasing order")
        self.play(self._spike_times)

    def _get_spike_times(self):
        return self._spike_times

    spike_times = property(fget=_get_spike_times,
                           fset=_set_spike_times)

    def clear_past_spikes(self):
        """If previous recordings are cleared, need to remove spikes from before the current time."""
        end = self._spike_times.indwhere(">", h.t)
        if end > 0:
            self._spike_times.remove(0, end - 1)  # range is inclusive

    
