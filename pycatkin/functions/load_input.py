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

    print('Loading input file: %s.' % input_path)

    with open(input_path) as file:
        pck_system = json.load(file)

    if 'states' in pck_system.keys():
        print('Reading states:')
        states = dict()
        for s in pck_system['states'].keys():
            print('* %s' % s)
            states[s] = State(name=s, **pck_system['states'][s])
    else:
        raise RuntimeError('Input file contains no states.')

    if 'scaling_relation_states' in pck_system.keys():
        print('Reading scaling relation states:')
        if states is None:
            states = dict()
        for s in pck_system['scaling_relation_states'].keys():
            print('* %s' % s)
            states[s] = Scaling(name=s, **pck_system['scaling_relation_states'][s])

    if 'system' in pck_system.keys():
        print('Reading system:')
        sys_params = pck_system['system']
        p = sys_params['p']
        print('* Pressure: %1.0f Pa' % p)
        T = sys_params['T']
        print('* Temperature: %1.0f K' % T)
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
        sim_system = System()
        sim_system.set_parameters(**sys_params)
        for s in states.keys():
            if states[s].gasdata is not None:
                states[s].gasdata['state'] = [states[i] for i in states[s].gasdata['state']]
            sim_system.add_state(state=states[s])

    else:
        raise RuntimeError('Input file contains no system details.')

    if 'reactions' in pck_system.keys():
        print('Reading reactions:')
        reactions = dict()
        for r in pck_system['reactions'].keys():
            print('* %s' % r)
            reactions[r] = Reaction(name=r, **pck_system['reactions'][r])
            reactions[r].reactants = [sim_system.states[s] for s in reactions[r].reactants]
            reactions[r].products = [sim_system.states[s] for s in reactions[r].products]
            if reactions[r].TS is not None:
                reactions[r].TS = [sim_system.states[s] for s in reactions[r].TS]
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
            sim_system.add_reaction(reaction=reactions[r])

    if 'reactor' in pck_system.keys():
        print('Reading reactor:')
        if pck_system['reactor'] == 'InfiniteDilutionReactor':
            print('* InfiniteDilutionReactor')
            reactor = InfiniteDilutionReactor()
        elif pck_system['reactor'] == 'CSTReactor':
            print('* CSTReactor')
            reactor = CSTReactor(**pck_system['reactor']['CSTReactor'])
        else:
            raise TypeError('Unknown reactor option, please choose InfiniteDilutionReactor or CSTReactor.')
        sim_system.add_reactor(reactor=reactor)
        sim_system.names_to_indices()
    else:
        if sim_system.reactions is not None:
            raise RuntimeError('Cannot consider reactions without reactor.' +
                               'To use constant boundary conditions, please specify InfiniteDilutionReactor.')

    if 'energy landscapes' in pck_system.keys():
        print('Reading energy landscapes:')
        for pes in pck_system['energy landscapes'].keys():
            print('* %s' % pes)
            minima = pck_system['energy landscapes'][pes]["minima"]
            labels = pck_system['energy landscapes'][pes]["labels"]
            minima = [[sim_system.states[s] for s in minima[k]] for k in range(len(minima))]
            labels = labels if labels else [i[0].name for i in minima]
            energy_landscape = Energy(name=pes, minima=minima, labels=labels)
            sim_system.add_energy_landscape(energy_landscape=energy_landscape)

    print('Done.')

    return sim_system

