'''
Project: SETO Open Energy Data Initiative (OEDI)
Author: Bo Chen @ Argonne National Lab
Date: 09/02/2022
Description: Convert imported model data into a format easier to use for optimization formulation
Note:
While it is rare, sometimes in OpenDSS multiple single-phase lines can be added between two three-phase buses.
This could result in infeasible solution if not merging them into a single line.
This could be done by adding a function inside edge_class to screen out and merge them.
'''

from model_import import func_model_import
import numpy as np



# Define the class for edge
class edge_class(object):
    # Default values:  -1 means undefined if info is not provided
    def __init__(self, Name='NoName', Bus_A='-1', Bus_B='-1', BusNum_A=0, BusNum_B=0, Length=-1, Device='-1',
                 Rperlen = np.zeros((3,3)), Xperlen = np.zeros((3,3)),Cperlen = np.zeros((3,3)),NormAmps = 1000.0, EmergAmps = 1000.0,
                 NumPhases = 3, Phases=np.array([1, 1, 1]), Control=False, Status=False, Outage=False):
        self.Name = Name
        self.Bus_A = Bus_A
        self.Bus_B = Bus_B
        self.BusNum_A = BusNum_A
        self.BusNum_B = BusNum_B
        self.Length = Length
        self.Device = Device    # Can be "line", "switch", "transformer", "regulator"
        self.Rperlen = Rperlen
        self.Xperlen = Xperlen
        self.Cperlen = Cperlen
        self.NormAmps = NormAmps
        self.EmergAmps = EmergAmps
        self.NumPhases = NumPhases
        self.Phases = Phases
        self.Control = Control  # if a line is controllable (True) or uncontrollable (False)
        self.Status = Status    # if a line is closed (True) or opened (False)
        self.Outage = Outage    # if a line is damaged (True) or normal (False)

### edge includes regular lines, switches, and transformers. Note  the line variable already includes both regular lines and switches
def make_edge_data(Circuit, Bus, Line, Transformer, Regcontrol):
    n_bus = Circuit['NumBuses']
    n_line = Circuit['NumLines']
    n_tran = Circuit['NumTrans']
    n_regctrl = Circuit['NumRegctrls']

    edge_list = []
    edge_dict = {}

    for key in Line:
        if key != "AllLineNames" and key != 'LineNumbers' and key!= 'NumLines':
            name = Line[key]['Name']
            busA = Line[key]['Bus1']
            busB = Line[key]['Bus2']
            busnumA = Bus[busA]['Number']
            busnumB = Bus[busB]['Number']
            Length = Line[key]['Length']
            if Length >= 0.01:
                Device = 'line'
                Control = False
            else:
                Device = 'switch'
                Control = True
            Rperlen = impedance_matrix_format(Line[key]['RmatrixPerLengthUnit'], Line[key]['Phases'])
            Xperlen = impedance_matrix_format(Line[key]['XmatrixPerLengthUnit'], Line[key]['Phases'])
            Cperlen = impedance_matrix_format(Line[key]['CmatrixPerLengthUnit'], Line[key]['Phases'])
            NormAmps = Line[key]['NormAmps']
            EmergAmps = Line[key]['EmergAmps']
            NumPhases = Line[key]['NumPhases']
            Phases = [0,0,0]
            for i in Line[key]['Phases']:
                Phases[i-1] = 1
            Status = Line[key]['Enabled']
            Outage = Line[key]['Outage']

            new_edge = edge_class(name,busA, busB, busnumA, busnumB, Length, Device, Rperlen, Xperlen, Cperlen, \
                                  NormAmps, EmergAmps, NumPhases, Phases, Control, Status, Outage)

            edge_list.append(new_edge)
            edge_dict[key] = new_edge

    regulator_names = [Regcontrol[key]['Transformer'] for key in Regcontrol if key != "AllTransNames" and key != 'TransNumbers' and key!= 'NumTrans']
    for key in Transformer:
        if key != "AllTransNames" and key != 'TransNumbers' and key!= 'NumTrans':
            name = Transformer[key]['Name']
            busA = Transformer[key]['Bus1']
            busB = Transformer[key]['Bus2']
            busnumA = Bus[busA]['Number']
            busnumB = Bus[busB]['Number']
            Length = 1 # Default for transformer length
            if name in regulator_names:
                Device = 'regulator'
                Control = True
            else:
                Device = 'transfromer'
                Control = False
            Rperlen = impedance_matrix_format([Transformer[key]['R']],Transformer[key]['Bus1Phases']) ####### make [R]
            Xperlen = impedance_matrix_format([Transformer[key]['Xhl']],Transformer[key]['Bus1Phases']) ##### make [Xhl]
            Cperlen = np.zeros((3,3))
            NormAmps = Transformer[key]['NormalAmps']
            EmergAmps = Transformer[key]['EmergAmps']
            NumPhases = Transformer[key]['Bus1Phases']
            Phases = [0,0,0]
            for i in Transformer[key]['Bus1Phases']:
                Phases[i-1] = 1
            Status = Transformer[key]['Enabled']
            Outage = Transformer[key]['Outage']

            new_edge = edge_class(name,busA, busB, busnumA, busnumB, Length, Device, Rperlen, Xperlen, Cperlen, \
                                  NormAmps, EmergAmps, NumPhases, Phases, Control, Status, Outage)

            edge_list.append(new_edge)
            edge_dict[key] = new_edge

    for i in range(len(edge_list)):
        edge_list[i].Number = i
        edge_dict[edge_list[i].Name].Number = i


    return edge_dict, edge_list

# convert any 1-phase, 2-phase impedance matrix to 3x3 matrix according to phase information
def impedance_matrix_format(imp_raw, phases):
    imp_matrix = np.zeros((3,3))
    if len(phases) == 3 and len(imp_raw) == 9:
        imp_matrix = np.array(imp_raw).reshape((3,3))
    elif len(phases) == 3 and len(imp_raw) == 1:
        imp_matrix[0,0] = imp_raw[0]
        imp_matrix[1,1] = imp_raw[0]
        imp_matrix[2,2] = imp_raw[0]
    elif len(phases) == 2 and len(imp_raw) == 4:
        imp_matrix[phases[0]-1, phases[0]-1] = imp_raw[0]
        imp_matrix[phases[0]-1, phases[1]-1] = imp_raw[1]
        imp_matrix[phases[1]-1, phases[0]-1] = imp_raw[2]
        imp_matrix[phases[1]-1, phases[1]-1] = imp_raw[3]
    elif len(phases) == 2 and len(imp_raw) == 1:
        imp_matrix[phases[0]-1, phases[0]-1] = imp_raw[0]
        imp_matrix[phases[0]-1, phases[1]-1] = imp_raw[0]
        imp_matrix[phases[1]-1, phases[0]-1] = imp_raw[0]
        imp_matrix[phases[1]-1, phases[1]-1] = imp_raw[0]
    elif len(phases) == 1:
        imp_matrix[phases[0]-1, phases[0]-1] = imp_raw[0]

    return imp_matrix

# merge single-phase regulators connecting to the same two buses into a three-phase regulator
# note single-phase definition will allow different local regulator setting for each phase, but this is not needed
# for centralized optimization, such as centralized VVO
# this is found in IEEE123 system
def merge_transformers(Transformer):
    keys_to_pop = []
    for key1 in Transformer:
        if key1 != "AllTransNames" and key1 != 'TransNumbers' and key1 != 'NumTrans':
            if 'checked' in Transformer[key1].keys():
                continue
            else:
                Transformer[key1]['checked'] = True
                for key2 in Transformer:
                    if key2 != "AllTransNames" and key2 != 'TransNumbers' and key2 != 'NumTrans':
                        if 'checked' in Transformer[key2].keys():
                            continue
                        else:
                            if Transformer[key1]['Bus1'] == Transformer[key2]['Bus1'] and \
                                Transformer[key1]['Bus2'] == Transformer[key2]['Bus2']:
                                new_phase = Transformer[key1]['Bus1Phases']+Transformer[key2]['Bus1Phases']
                                new_phase.sort()
                                Transformer[key1]['Bus1Phases'] = new_phase
                                new_phase = Transformer[key1]['Bus2Phases']+Transformer[key2]['Bus2Phases']
                                new_phase.sort()
                                Transformer[key1]['Bus2Phases'] = new_phase

                                Transformer[key2]['checked'] = True
                                keys_to_pop.append(key2)

    for key in keys_to_pop:
        Transformer.pop(key)
        Transformer['AllTransNames'].remove(key)
        Transformer['NumTrans'] -= 1

    Transformer['TransNumbers'] = list(range(Transformer['NumTrans']))

    return Transformer

def make_bus_data(Bus):
    bus_list = []
    bus_dict = {}
    for key in Bus:
        if key != 'AllBusNames' and key != 'BusNumbers' and key != 'NumBuses':
            a=[0,0,0]
            for i in Bus[key]['Phases']:
                a[i-1] = 1
            Bus[key]['Phases'] = a

            bus_dict[key] = Bus[key]
            bus_list.append(Bus[key])
    return bus_dict, bus_list

def make_load_data(Bus, Load):
    load_list = []
    load_dict = {}
    for key in Load:
        if key != 'AllLoadNames' and key != 'LoadNumbers' and key != 'NumLoads':
            a=[0,0,0]
            for i in Load[key]['Phases']:
                a[i-1] = 1
            Load[key]['Phases'] = a
            Load[key]['BusNum'] = Bus[Load[key]['BusName']]['Number']

            load_dict[key] = Load[key]
            load_list.append(Load[key])

    return load_dict, load_list


def make_cap_data(Bus, Capacitor):
    cap_list = []
    cap_dict = {}
    for key in Capacitor:
        if key != 'AllCapacitorNames' and key != 'CapNumbers' and key != 'NumCaps':
            a=[0,0,0]
            for i in Capacitor[key]['Phases']:
                a[i-1] = 1
            Capacitor[key]['Phases'] = a
            Capacitor[key]['BusNum'] = Bus[Capacitor[key]['BusName']]['Number']

            cap_dict[key] = Capacitor[key]
            cap_list.append(Capacitor[key])

    return cap_dict, cap_list

def make_loadcap_data(load_dict, load_list, cap_dict, cap_list):
    loadcap_dict = load_dict.copy()
    for key in loadcap_dict:
        loadcap_dict[key]['Device'] = 'load'
    for key in cap_dict:
        loadcap_dict[key] = cap_dict[key]
        loadcap_dict[key]['Device'] = 'capacitor'
    loadcap_list = load_list.copy()
    for l in loadcap_list:
        l['Device'] = 'load'
    for l in cap_list:
        loadcap_list.append(l)
        l['Device'] = 'capacitor'

    return loadcap_dict, loadcap_list

def make_gen_data(Substation):
    gen_dict = {}
    gen_list = []
    for key in Substation:
        if key != 'AllSubNames' and key != 'SubNumbers' and key != 'NumSubs':
            gen_dict[key] = Substation[key]
            gen_list.append(Substation[key])
    return gen_dict, gen_list

def format_all_data(Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation):
    Transformer = merge_transformers(Transformer)  # Regcontrol is not updated in this function for now
    [edge_dict, edge_list] = make_edge_data(Circuit, Bus, Line, Transformer, Regcontrol)
    [bus_dict, bus_list] = make_bus_data(Bus)
    [load_dict, load_list] = make_load_data(Bus, Load)
    [cap_dict, cap_list] = make_cap_data(Bus, Capacitor)
    [loadcap_dict, loadcap_list] = make_loadcap_data(load_dict, load_list, cap_dict, cap_list)
    [gen_dict, gen_list] = make_gen_data(Substation)
    [block_list, link_list, A] = make_block_data(bus_dict, edge_dict)

    return edge_dict, edge_list, bus_dict, bus_list, load_dict, load_list, cap_dict, cap_list, \
           loadcap_dict, loadcap_list, gen_dict, gen_list, block_list, link_list, A

def make_block_data(bus_dict,edge_dict):
    # Group directly connected nodes and lines based on the graph determined by edge_list[i].Control
    bus_set = set(bus_dict.keys())
    block_list = find_islands(bus_dict, edge_dict)
    [link_list, A] = find_links(block_list, bus_set, edge_dict)

    return block_list, link_list, A

class graph_node(object):
    def __init__(self, name):
        self.__name = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

def connected_components(nodes):
    # List of connected components found. The order is random.
    result = []

    # Make a copy of the set, so we can modify it.
    nodes = set(nodes)

    # Iterate while we still have nodes to process.
    while nodes:

        # Get a random node and remove it from the global set.
        n = nodes.pop()

        # This set will contain the next group of nodes connected to each other.
        group = {n}

        # Build a queue with this node in it.
        queue = [n]

        # Iterate the queue.
        # When it's empty, we finished visiting a group of connected nodes.
        while queue:
            # Consume the next item from the queue.
            n = queue.pop()

            # Fetch the neighbors.
            neighbors = n.links

            # Remove the neighbors we already visited.
            neighbors.difference_update(group)

            # Remove the remaining nodes from the global set.
            nodes.difference_update(neighbors)

            # Add them to the group of connected nodes.
            group.update(neighbors)

            # Add them to the queue, so we visit them in the next iterations.
            queue.extend(neighbors)

        # Add the group to the list of groups.
        result.append(group)

    # Return the list of groups.
    return result

def find_islands(bus_dict, edge_dict):
    # Thanks to the author of the algorthm: https://breakingcode.wordpress.com/2013/04/08/finding-connected-components-in-a-graph/
    graph_node_set = set()
    graph_dict = {}  # pay extra attention here
    for i in bus_dict:
        graph_dict[i] = graph_node(bus_dict[i]['Name'])
        graph_node_set.add(graph_dict[i])
    for i in edge_dict:
        # if edge_dict[i].Device != 'switch' and edge_dict[i].Device != 'regulator':
        # if edge_dict[i].Device != 'switch' and edge_dict[i].Name != 'reg1a':
        if edge_dict[i].Device != 'switch':
            busA = edge_dict[i].Bus_A
            busB = edge_dict[i].Bus_B
            graph_dict[busA].add_link(graph_dict[busB])
    graph_islands = connected_components(graph_node_set)
    islands = []
    for i in graph_islands:
        subislands = set()
        for j in i:
            subislands.update([j.name])  # not j.name, must be [j.name] as a string
        islands.append(subislands)

    return islands


def find_links(block_list, node_set, edge_dict):
    links = []
    A = np.zeros([len(block_list),len(block_list)])
    A = A.astype(list)
    for e in edge_dict:
        if edge_dict[e].Device != 'switch':
            continue
        else:
            for b_idx,b in enumerate(block_list):
                if edge_dict[e].Bus_A in b:
                    for bb_idx,bb in enumerate(block_list):
                        if edge_dict[e].Bus_B in bb:
                            links.append([b_idx,bb_idx, edge_dict[e].Number])
                            A[b_idx, bb_idx] = 1
                            A[bb_idx, b_idx] = 1
    return links, A



if __name__ == '__main__':

    feeder_folder = '123Bus'
    dss_master_file = 'IEEE123Master_FLISR.dss'
    Substation_names = ['150']

    [Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation] = \
        func_model_import(feeder_folder, dss_master_file, Substation_names)

    [edge_dict, edge_list, bus_dict, bus_list, load_dict, load_list, cap_dict, cap_list, loadcap_dict, loadcap_list, gen_dict, gen_list, block_list, link_list, A] = \
        format_all_data(Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation)

    a=1
