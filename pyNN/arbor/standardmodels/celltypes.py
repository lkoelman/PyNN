# encoding: utf-8
"""
Definition of cell types for the Arbor backend (implemenation 
of standard base_cells).

These are the cell types passed to a Population, i.e. the specification
of the cell type. Their 'model' attribute points to a class that is
the implementation of the cell type.

:copyright: Copyright 2006-2016 by the PyNN team, see AUTHORS.
:license: CeCILL, see LICENSE for details.


Developer Notes
---------------

Each cell type implements a PyNN cell type from module pyNN.standardmodels.cells.

Its 'model' attribute points to a class implementing the cell type: for the
Arbor backend they are subclasses of an Arbor cell description.
"""

from pyNN.standardmodels import cells as base_cells, build_translations
from .celltypes import (StandardIF, Izhikevich_,
                           RandomSpikeSource, VectorSpikeSource,
                           RandomGammaSpikeSource,
                           RandomPoissonRefractorySpikeSource)
import pyarb as arb
import logging

logger = logging.getLogger("PyNN")


# LIF with alpha-function current-based synapses
class IF_curr_alpha(base_cells.IF_curr_alpha):

    __doc__ = base_cells.IF_curr_alpha.__doc__

    # Translate base celltype parameters to self.model attributes
    translations = build_translations(
        ('tau_m',      'tau_m'),
        ('cm',         'c_m'),
        ('v_rest',     'v_rest'),
        ('v_thresh',   'v_thresh'),
        ('v_reset',    'v_reset'),
        ('tau_refrac', 't_refrac'),
        ('i_offset',   'i_offset'),
        ('tau_syn_E',  'tau_e'),
        ('tau_syn_I',  'tau_i'),
    )

    model = StandardIF # subclass of Arbor cell description class
    arb_cell_kind = arb.cell_kind.lif_neuron

    extra_parameters = {'syn_type': 'current',
                        'syn_shape': 'alpha'}


# LIF with exponential-function current-based synapses
class IF_curr_exp(base_cells.IF_curr_exp):

    __doc__ = base_cells.IF_curr_exp.__doc__

    translations = build_translations(
        ('tau_m',      'tau_m'),
        ('cm',         'c_m'),
        ('v_rest',     'v_rest'),
        ('v_thresh',   'v_thresh'),
        ('v_reset',    'v_reset'),
        ('tau_refrac', 't_refrac'),
        ('i_offset',   'i_offset'),
        ('tau_syn_E',  'tau_e'),
        ('tau_syn_I',  'tau_i'),
    )
    
    model = StandardIF # subclass of Arbor cell description class
    arb_cell_kind = arb.cell_kind.lif_neuron

    extra_parameters = {'syn_type': 'current',
                        'syn_shape': 'exp'}


# LIF with alpha-function conductance-based synapses
class IF_cond_alpha(base_cells.IF_cond_alpha):

    __doc__ = base_cells.IF_cond_alpha.__doc__

    translations = build_translations(
        ('tau_m',      'tau_m'),
        ('cm',         'c_m'),
        ('v_rest',     'v_rest'),
        ('v_thresh',   'v_thresh'),
        ('v_reset',    'v_reset'),
        ('tau_refrac', 't_refrac'),
        ('i_offset',   'i_offset'),
        ('tau_syn_E',  'tau_e'),
        ('tau_syn_I',  'tau_i'),
        ('e_rev_E',    'e_e'),
        ('e_rev_I',    'e_i')
    )
    
    model = StandardIF # subclass of Arbor cell description class
    arb_cell_kind = arb.cell_kind.lif_neuron

    extra_parameters = {'syn_type': 'conductance',
                        'syn_shape': 'alpha'}


# LIF with exponential-function conductance-based synapses
class IF_cond_exp(base_cells.IF_cond_exp):

    __doc__ = base_cells.IF_cond_exp.__doc__

    translations = build_translations(
        ('tau_m',      'tau_m'),
        ('cm',         'c_m'),
        ('v_rest',     'v_rest'),
        ('v_thresh',   'v_thresh'),
        ('v_reset',    'v_reset'),
        ('tau_refrac', 't_refrac'),
        ('i_offset',   'i_offset'),
        ('tau_syn_E',  'tau_e'),
        ('tau_syn_I',  'tau_i'),
        ('e_rev_E',    'e_e'),
        ('e_rev_I',    'e_i')
    )
    
    model = StandardIF # subclass of Arbor cell description class
    arb_cell_kind = arb.cell_kind.lif_neuron

    extra_parameters = {'syn_type': 'conductance',
                        'syn_shape': 'exp'}


class Izhikevich(base_cells.Izhikevich):
    __doc__ = base_cells.Izhikevich.__doc__

    translations = build_translations(
        ('a',        'a_'),
        ('b',        'b'),
        ('c',        'c'),
        ('d',        'd'),
        ('i_offset', 'i_offset')
    )
    
    model = Izhikevich_
    arb_cell_kind = arb.cell_kind.lif_neuron # TODO: find Arbor cell_kind for Izhikevich cells


################################################################################
# Artificial Spike Sources
################################################################################

# TODO: Replace Arbor cell_kind.data_spike_source by something more suitable

class SpikeSourcePoisson(base_cells.SpikeSourcePoisson):

    __doc__ = base_cells.SpikeSourcePoisson.__doc__

    translations = build_translations(
        ('start',    'start'),
        ('rate',     '_interval',  "1000.0/rate",  "1000.0/_interval"),
        ('duration', 'duration'),
    )

    model = RandomSpikeSource
    arb_cell_kind = arb.cell_kind.data_spike_source


class SpikeSourcePoissonRefractory(base_cells.SpikeSourcePoissonRefractory):

    __doc__ = base_cells.SpikeSourcePoissonRefractory.__doc__

    translations = build_translations(
        ('start',      'start'),
        ('rate',       'rate'),
        ('tau_refrac', 'tau_refrac'),
        ('duration',   'duration'),
    )

    model = RandomPoissonRefractorySpikeSource
    arb_cell_kind = arb.cell_kind.data_spike_source


class SpikeSourceGamma(base_cells.SpikeSourceGamma):
    __doc__ = base_cells.SpikeSourceGamma.__doc__

    translations = build_translations(
        ('alpha',    'alpha'),
        ('beta',     'beta',    0.001),
        ('start',    'start'),
        ('duration', 'duration'),
    )

    model = RandomGammaSpikeSource
    arb_cell_kind = arb.cell_kind.data_spike_source


class SpikeSourceArray(base_cells.SpikeSourceArray):

    __doc__ = base_cells.SpikeSourceArray.__doc__

    translations = build_translations(
        ('spike_times', 'spike_times'),
    )

    model = VectorSpikeSource
    arb_cell_kind = arb.cell_kind.data_spike_source



