import json
from pycatkin.classes.state import *
from pycatkin.classes.reaction import *
from pycatkin.classes.system import *
from pycatkin.classes.reactor import *
from pycatkin.classes.scaling import *
from pycatkin.classes.energy import *


def read_from_input_file(input_path='input.json'):
    """Reads simulation setup including mechanism,
    reaction conditions and solver settings from
    input json file. Creates a new simulator object.

    Returns the simulator.
    """

    print('Loading input file...')

    states = None
    reactions = None
    reactor = None
    sim_system = None

    with open(input_path) as file:
        pck_system = json.load(file)

    if 'states' in pck_system.keys():
        states = dict()
        for s in pck_system['states'].keys():
            states[s] = State(name=s, **pck_system['states'][s])

    if 'scaling_relation_states' in pck_system.keys():
        if states is None:
            states = dict()
        for s in pck_system['scaling_relation_states'].keys():
            states[s] = Scaling(name=s, **pck_system['scaling_relation_states'][s])

    if 'reactions' in pck_system.keys():
        reactions = dict()
        for r in pck_system['reactions'].keys():
            reactions[r] = Reaction(name=r, **pck_system['reactions'][r])
            reactions[r].reactants = [states[s] for s in reactions[r].reactants]
            reactions[r].products = [states[s] for s in reactions[r].products]
            if reactions[r].TS is not None:
                reactions[r].TS = [states[s] for s in reactions[r].TS]
        for r in reactions.keys():
            for s in reactions[r].reactants + reactions[r].products:
                if isinstance(s, Scaling):
                    for sr in s.scaling_reactions.keys():
                        s.scaling_reactions[sr]['reaction'] = reactions[s.scaling_reactions[sr]['reaction']]
            if reactions[r].TS is not None:
                for s in reactions[r].TS:
                    if isinstance(s, Scaling):
                        for sr in s.scaling_reactions.keys():
                            s.scaling_reactions[sr]['reaction'] = reactions[s.scaling_reactions[sr]['reaction']]

    if 'reactor' in pck_system.keys():
        if pck_system['reactor'] == 'InfiniteDilutionReactor':
            reactor = InfiniteDilutionReactor()
        elif pck_system['reactor'] == 'CSTReactor':
            reactor = CSTReactor(**pck_system['reactor']['CSTReactor'])
        else:
            print('Unknown reactor option, using default: InfiniteDilutionReactor.')
            reactor = InfiniteDilutionReactor()

    if 'system' in pck_system.keys():
        sys_params = pck_system['system']
        p = sys_params['p']
        startsites = 0.0
        if 'start_state' in sys_params.keys():
            for s in sys_params['start_state'].keys():
                if states[s].state_type == 'gas':
                    sys_params['start_state'][s] = sys_params['start_state'][s] * p / bartoPa
                elif states[s].state_type == 'surface' or states[s].state_type == 'adsorbate':
                    startsites += sys_params['start_state'][s]
            if startsites == 0.0:
                raise ValueError('Initial surface coverage cannot be zero for all states!')
        if 'inflow_state' in sys_params.keys():
            for s in sys_params['inflow_state'].keys():
                if states[s].state_type == 'gas':
                    sys_params['inflow_state'][s] = sys_params['inflow_state'][s] * p / bartoPa
                else:
                    raise TypeError('Only gas states can comprise the inflow!')
        if reactor is not None and reactions is not None:
            sim_system = System(reactor=reactor, reactions=reactions)
            sim_system.set_parameters(**sys_params)

    print('Done.')

    return sim_system
