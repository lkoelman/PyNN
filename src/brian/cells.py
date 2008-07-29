# ==============================================================================
# Standard cells for nest1
# $Id: cells.py 294 2008-04-04 12:07:56Z apdavison $
# ==============================================================================

from pyNN import common
import brian_no_units_no_warnings
from brian.library.synapses import *
import brian


class IF_curr_alpha(common.IF_curr_alpha):
    """Leaky integrate and fire model with fixed threshold and alpha-function-
    shaped post-synaptic current."""
    translations = common.build_translations(
        ('v_rest',     'v_rest', mV),
        ('v_reset',    'v_reset', mV),
        ('cm',         'cm'), # C is in pF, cm in nF
        ('tau_m',      'tau_m', ms),
        ('tau_refrac', 'tau_refrac',"max(get_time_step(), tau_refrac*0.001)", "tau_refrac"),
        ('tau_syn_E',  'tau_syn_E', ms),
        ('tau_syn_I',  'tau_syn_I', ms),
        ('v_thresh',   'v_thresh', mV),
        ('i_offset',   'i_offset', pA), # I0 is in pA, i_offset in nA
        ('v_init',     'v', mV),
    )
    eqs= brian.Equations('''
        dv/dt  = (ge + gi + (v_rest-v))/tau_m : volt
        dge/dt = (ye-ge)/tau_syn_E              : 1 
        dye/dt = -ye/tau_syn_E                  : 1
        dgi/dt = (yi-gi)/tau_syn_I              : 1 
        dyi/dt = -yi/tau_syn_I                  : 1
        tau_syn_E : second
        tau_syn_I : second
        tau_m     : second
        v_rest    : volt
        '''
        )


class IF_curr_exp(common.IF_curr_exp):
    """Leaky integrate and fire model with fixed threshold and
    decaying-exponential post-synaptic current. (Separate synaptic currents for
    excitatory and inhibitory synapses."""
    
    translations = common.build_translations(
        ('v_rest',     'v_rest', 0.001),
        ('v_reset',    'v_reset', 0.001),
        ('cm',         'cm'), # C is in pF, cm in nF
        ('tau_m',      'tau_m', 0.001),
        ('tau_refrac', 'tau_refrac',"max(get_time_step(), tau_refrac*0.001)", "tau_refrac"),
        ('tau_syn_E',  'tau_syn_E', 0.001),
        ('tau_syn_I',  'tau_syn_I', 0.001),
        ('v_thresh',   'v_thresh', 0.001),
        ('i_offset',   'i_offset'), # I0 is in pA, i_offset in nA
        ('v_init',     'v', 0.001),
    )
    eqs= brian.Equations('''
        dv/dt  = (ge + gi + (v_rest-v))/tau_m : volt
        dge/dt = -ge/tau_syn_E              : 1
        dgi/dt = -gi/tau_syn_I              : 1
        tau_syn_E                           : second
        tau_syn_I                           : second
        tau_m                               : second
        v_rest                              : volt
        '''
        )
    

class IF_cond_alpha(common.ModelNotAvailable):
    translations = common.build_translations(
        ('v_rest',     'v_rest',0.001)    ,
        ('v_reset',    'v_reset',0.001),
        ('cm',         'cm'), # C_m is in pF, cm in nF
        ('tau_m',      'tau_m',0.001),
        ('tau_refrac', 'tau_refrac',"max(get_time_step(), tau_refrac*0.001)", "tau_refrac"),
        ('tau_syn_E',  'tau_syn_E',0.001),
        ('tau_syn_I',  'tau_syn_I',0.001),
        ('v_thresh',   'v_thresh',0.001),
        ('i_offset',   'i_offset',0.001), # I_e is in pA, i_offset in nA
        ('e_rev_E',    'e_rev_E',0.001),
        ('e_rev_I',    'e_rev_I',0.001),
        ('v_init',     'v',0.001),
    )
    eqs= brian.Equations('''
        dv/dt  = ((v_rest-v) + ge*(e_rev_E-v) + gi*(e_rev_I-v))/tau_m : volt
        dge/dt = ye-ge/tau_syn_E              : 1 
        dye/dt = -ye/tau_syn_E                : 1
        dgi/dt = yi-gi/tau_syn_I              : 1 
        dyi/dt = -yi/tau_syn_I                : 1
        tau_syn_E              : second
        tau_syn_I              : second
        tau_m                  : second
        v_rest                 : volt
        e_rev_E                : volt
        e_rev_I                : volt
        '''
        )


class IF_cond_exp(common.IF_cond_exp):
    """Leaky integrate and fire model with fixed threshold and 
    exponentially-decaying post-synaptic conductance."""
    translations = common.build_translations(
        ('v_rest',     'v_rest',0.001)    ,
        ('v_reset',    'v_reset',0.001),
        ('cm',         'cm'), # C_m is in pF, cm in nF
        ('tau_m',      'tau_m',0.001),
        ('tau_refrac', 'tau_refrac',"max(get_time_step(), tau_refrac*0.001)", "tau_refrac"),
        ('tau_syn_E',  'tau_syn_E',0.001),
        ('tau_syn_I',  'tau_syn_I',0.001),
        ('v_thresh',   'v_thresh',0.001),
        ('i_offset',   'i_offset',0.001), # I_e is in pA, i_offset in nA
        ('e_rev_E',    'e_rev_E',0.001),
        ('e_rev_I',    'e_rev_I',0.001),
        ('v_init',     'v',0.001),
    )
    eqs= brian.Equations('''
        dv/dt  = ((v_rest-v) + ge*(e_rev_E-v) + gi*(e_rev_I-v))/tau_m : volt
        dge/dt = -ge/tau_syn_E : 1
        dgi/dt = -gi/tau_syn_I : 1
        tau_syn_E              : second
        tau_syn_I              : second
        tau_m                  : second
        v_rest                 : volt
        e_rev_E                : volt
        e_rev_I                : volt
        '''
        )

class IF_facets_hardware1(common.IF_facets_hardware1):
    """Leaky integrate and fire model with conductance-based synapses and fixed
    threshold as it is resembled by the FACETS Hardware Stage 1. For further
    details regarding the hardware model see the FACETS-internal Wiki:
    https://facets.kip.uni-heidelberg.de/private/wiki/index.php/WP7_NNM
    """
    pass


class SpikeSourcePoisson(common.SpikeSourcePoisson):
    """Spike source, generating spikes according to a Poisson process."""
    translations = common.build_translations(
        ('rate',     'rate'),
        ('start',    'start', 0.001),
        ('duration', 'duration', 0.001),
    )
    
    def __init__(self, parameters):
        common.SpikeSourcePoisson.__init__(self, parameters)
        start    = self.parameters['start']
        duration = self.parameters['duration']
        rate     = self.parameters['rate']
        self.fct = lambda t: (start <= t <= start + duration and rate)
    
class SpikeSourceArray(common.SpikeSourceArray):
    """Spike source generating spikes at the times given in the spike_times array."""


class EIF_cond_alpha_isfa_ista(common.ModelNotAvailable):
    pass

class HH_cond_exp(common.HH_cond_exp):
    
    translations = common.build_translations(
        ('gbar_Na',    'gbar_Na'),   
        ('gbar_K',     'gbar_K'),    
        ('g_leak',     'g_leak'),    
        ('cm',         'cm'),  
        ('v_offset',   'v_offset'),
        ('e_rev_Na',   'e_rev_Na'),
        ('e_rev_K',    'e_rev_K'), 
        ('e_rev_leak', 'e_rev_leak'),
        ('e_rev_E',    'e_rev_E'),
        ('e_rev_I',    'e_rev_I'),
        ('tau_syn_E',  'tau_syn_E'),
        ('tau_syn_I',  'tau_syn_I'),
        ('i_offset',   'i_offset'),
        ('v_init',     'v'),
    )
    
    eqs= brian.Equations('''
        dv/dt = (g_leak*(e_rev_leak-v)+ge*(e_rev_E-v)+gi*(e_rev_I-v)-gbar_Na*(m*m*m)*h*(v-e_rev_Na)-gbar_K*(n*n*n*n)*(v-e_rev_K))/cm : volt
        dm/dt = alpham*(1-m)-betam*m : 1
        dn/dt = alphan*(1-n)-betan*n : 1
        dh/dt = alphah*(1-h)-betah*h : 1
        dge/dt = -ge/tau_syn_E : siemens
        dgi/dt = -gi/tau_syn_I : siemens
        alpham = 0.32*(mV**-1)*(13*mV-v+VT)/(exp((13*mV-v+VT)/(4*mV))-1.)/ms : Hz
        betam = 0.28*(mV**-1)*(v-VT-40*mV)/(exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
        alphah = 0.128*exp((17*mV-v+VT)/(18*mV))/ms : Hz
        betah = 4./(1+exp((40*mV-v+VT)/(5*mV)))/ms : Hz
        alphan = 0.032*(mV**-1)*(15*mV-v+VT)/(exp((15*mV-v+VT)/(5*mV))-1.)/ms : Hz
        betan = .5*exp((10*mV-v+VT)/(40*mV))/ms : Hz
        tau_syn_E              : second
        tau_syn_I              : second
        e_rev_E                : volt
        e_rev_I                : volt
        e_rev_Na               : volt
        e_rev_K                : volt
        e_rev_leak             : volt
        gbar_Na                : nS
        gbar_K                 : nS
        gbar_leak              : nS
        v_offset               : volt
        cm                     : nF
    ''')

class SpikeSourceInhGamma(common.ModelNotAvailable):
    pass

class IF_cond_exp_gsfa_grr(common.ModelNotAvailable):
    pass