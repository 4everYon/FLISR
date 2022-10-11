'''
Project: SETO Open Energy Data Initiative (OEDI)
Author: Bo Chen @ Argonne National Lab
Date: 09/18/2022
Description: Outage isolation algorithm using MILP formulation
Note:
Ref: B. Chen, Z. Ye, C. Chen and J. Wang, "Toward a MILP Modeling Framework for Distribution System Restoration,"
in IEEE Transactions on Power Systems, vol. 34, no. 3, pp. 1749-1760, May 2019
'''

from data_preparation import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyomo.environ as pyo


def impedance_matrix(n_edge, n_bus, edge_list, Zbase):
    Z = np.zeros((3 * n_bus, 3 * n_bus)) * 1j
    for i in range(n_edge):
        Nfrom = edge_list[i].BusNum_A + 1
        Nto = edge_list[i].BusNum_B + 1
        Imp = edge_list[i].Rperlen * edge_list[i].Length + edge_list[i].Xperlen * edge_list[i].Length * 1j
        Imp_pu = Imp/Zbase  # impedance per mile
        Z[(3 * Nfrom - 2 - 1):(3 * Nfrom), (3 * Nto - 2 - 1):(3 * Nto)] = Imp_pu
        Z[(3 * Nto - 2 - 1):(3 * Nto), (3 * Nfrom - 2 - 1):3 * Nfrom] = Imp_pu
    return Z

if __name__ == '__main__':

    feeder_folder = '123Bus'
    dss_master_file = 'IEEE123Master_FLISR.dss'
    Substation_names = ['150']

    [Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation] = \
        func_model_import(feeder_folder, dss_master_file, Substation_names)

    [edge_dict, edge_list, bus_dict, bus_list, load_dict, load_list, cap_dict, cap_list, loadcap_dict, loadcap_list, gen_dict, gen_list, block_list, link_list, A] = \
        format_all_data(Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation)


    # Configurate the optimization model constraints
    st_pf_balance        = True
    st_pf_voltage        = True
    st_voltage_limit     = True
    st_line_capacity     = True
    st_gen_capacity      = True

    st_regulator_control = True
    st_cap_control = True

    iso_Vsrc = 1.05  # expected voltage in per unit of the black-start DG
    M = 1000000  # value used in the big-M method.
    iso_vmin = 0.95
    iso_vmax = 1.05

    Vbase = Substation['150']['kVBase']*1000
    Sbase = 1000 # KVA

    n_edge = len(edge_list)
    n_bus = len(bus_list)
    n_loadcap = len(loadcap_list)
    n_load = len(load_list)
    n_gen = len(gen_list)
    n_block = len(block_list)
    n_link = len(link_list)

    # voltage regulator pre-defined positions according to the original files
    for i in range(n_edge):
        if edge_list[i].Device == 'regulator':
            if edge_list[i].Bus_A == '150':
                edge_list[i].Position = np.array([7*0.00625+1, 1, 1])
            if edge_list[i].Bus_A == '160':
                edge_list[i].Position = [8*0.00625+1, 1*0.00625+1, 5*0.00625+1]
            if edge_list[i].Bus_A == '25':
                edge_list[i].Position = [1, 0, -1*0.00625+1]
            if edge_list[i].Bus_A == '9':
                edge_list[i].Position = [-1*0.00625+1, 0, 0]

    #######################################################################
    # switch initial status to be derived from OpenDSS data
    swi_status = []
    swi_idx = []
    for i in range(n_edge):
        if edge_list[i].Device == 'switch':
            if edge_list[i].Status == True:
                swi_status.append(1)
                swi_idx.append(i)
            else:
                swi_status.append(0)
    n_switch = len(swi_status)

    # in IEEE123 all switches are closed. Manually set sw7 and sw8 to be opened
    swi_status = [1,1,1,1,1,1,0,0]

    ################################################################################
    ## damage information
    edge_list[50].Outage = True
    edge_list[104].Outage = True




    ## Update bus and block outage status
    for i in range(n_edge):
        if edge_list[i].Outage == True:
            bus_list[edge_list[i].BusNum_A]['Outage'] = True
            bus_list[edge_list[i].BusNum_B]['Outage'] = True
            bus_dict[edge_list[i].Bus_A]['Outage'] = True
            bus_dict[edge_list[i].Bus_B]['Outage'] = True
    for i in range(n_gen):
        if gen_list[i]['Outage'] == True:
            bus_list[gen_list[i]['BusNum']]['Outage'] = True
            bus_dict[gen_list[i]['BusName']]['Outage'] = True
    for i in range(n_loadcap):
        if loadcap_list[i]['Outage'] == True:
            bus_list[loadcap_list[i]['BusNum']]['Outage'] = True
            bus_dict[loadcap_list[i]['BusName']]['Outage'] = True
    block_status = np.ones(n_block) # default 1 if no outages inside this block
    for i in range(n_block):
        for n in block_list[i]:
            if bus_dict[n]['Outage'] == True:
                block_status[i] = 0


    # Get Z matrix, which is 3*n_bus X 3*n_bus dimension, each 3X3 block represent the impedance between two buses (and self impedance)
    Zbase = Vbase**2/Sbase
    Z = impedance_matrix(n_edge, n_bus, edge_list, Zbase)
    
    # Develop Incident Matrices
    I_gen_bus = np.zeros((n_gen, n_bus))
    for i in range(n_gen):
        for j in range(n_bus):
            if gen_list[i]['BusNum'] == bus_list[j]['Number']:
                I_gen_bus[i,j] = 1
            else:
                I_gen_bus[i,j] = 0
    
    I_loadcap_bus = np.zeros((n_loadcap, n_bus))
    for i in range(n_loadcap):
        for j in range(n_bus):
            if loadcap_list[i]['BusNum'] == bus_list[j]['Number']:
                I_loadcap_bus[i,j] = -1
            else:
                I_loadcap_bus[i,j] = 0
    
    I_edge_bus = np.zeros((n_edge, n_bus))
    for i in range(n_edge):
        for j in range(n_bus):
            if edge_list[i].BusNum_A == bus_list[j]['Number']:
                I_edge_bus[i,j] = -1  # leave "FROM" node
            elif edge_list[i].BusNum_B == bus_list[j]['Number']:
                I_edge_bus[i,j] = 1   # inject into "TO" node
            else:
                I_edge_bus[i,j] = 0



    # Declear concrete model in pyomo
    m = pyo.ConcreteModel()


    # define sets
    m.n_gen = pyo.Set(initialize = [i for i in range(n_gen)])
    m.n_loadcap = pyo.Set(initialize = [i for i in range(n_loadcap)])
    m.n_edge = pyo.Set(initialize = [i for i in range(n_edge)])
    m.n_bus = pyo.Set(initialize = [i for i in range(n_bus)])
    m.n_block = pyo.Set(initialize = [i for i in range(n_block)])
    m.n_link = pyo.Set(initialize = [i for i in range(n_link)])

    # define decision variables
    m.x_gen = pyo.Var(m.n_gen, domain=pyo.Binary)
    m.p_gen_a = pyo.Var(m.n_gen, bounds=(0, M))
    m.p_gen_b = pyo.Var(m.n_gen, bounds=(0, M))
    m.p_gen_c = pyo.Var(m.n_gen, bounds=(0, M))
    m.q_gen_a = pyo.Var(m.n_gen, bounds=(-M, M))
    m.q_gen_b = pyo.Var(m.n_gen, bounds=(-M, M))
    m.q_gen_c = pyo.Var(m.n_gen, bounds=(-M, M))

    m.x_loadcap = pyo.Var(m.n_loadcap, domain=pyo.Binary)
    m.p_load_a = pyo.Var(m.n_loadcap, bounds=(0, M))
    m.p_load_b = pyo.Var(m.n_loadcap, bounds=(0, M))
    m.p_load_c = pyo.Var( m.n_loadcap, bounds=(0, M))
    m.q_load_a = pyo.Var(m.n_loadcap, bounds=(-M, M))
    m.q_load_b = pyo.Var(m.n_loadcap, bounds=(-M, M))
    m.q_load_c = pyo.Var(m.n_loadcap, bounds=(-M, M))

    m.x_edge = pyo.Var(m.n_edge, domain=pyo.Binary)
    m.p_edge_a = pyo.Var(m.n_edge, bounds=(-M, M))
    m.p_edge_b = pyo.Var(m.n_edge, bounds=(-M, M))
    m.p_edge_c = pyo.Var(m.n_edge, bounds=(-M, M))
    m.q_edge_a = pyo.Var(m.n_edge, bounds=(-M, M))
    m.q_edge_b = pyo.Var(m.n_edge, bounds=(-M, M))
    m.q_edge_c = pyo.Var(m.n_edge, bounds=(-M, M))
    m.rg_a2 = pyo.Var(m.n_edge, bounds=(-0.81, 1.21))
    m.rg_b2 = pyo.Var(m.n_edge, bounds=(-0.81, 1.21))
    m.rg_c2 = pyo.Var(m.n_edge, bounds=(-0.81, 1.21))

    m.x_bus = pyo.Var(m.n_bus, domain=pyo.Binary)
    m.v_bus_a2 = pyo.Var(m.n_bus, bounds=(0, M))
    m.v_bus_b2 = pyo.Var(m.n_bus, bounds=(0, M))
    m.v_bus_c2 = pyo.Var(m.n_bus, bounds=(0, M))

    m.x_block = pyo.Var(m.n_block, domain = pyo.Binary)
    m.x_link = pyo.Var(m.n_link, domain = pyo.Binary)

    # Define constraints
    m.constraint = pyo.ConstraintList()

    # Formulate the isolation model
    if st_pf_balance == True:
        # Power availability
        for e in range(n_edge):
            m.constraint.add(-1 * M * m.x_edge[e] <= m.p_edge_a[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.p_edge_a[e])
            m.constraint.add(-1 * M * m.x_edge[e] <= m.p_edge_b[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.p_edge_b[e])
            m.constraint.add(-1 * M * m.x_edge[e] <= m.p_edge_c[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.p_edge_c[e])
            m.constraint.add(-1 * M * m.x_edge[e] <= m.q_edge_a[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.q_edge_a[e])
            m.constraint.add(-1 * M * m.x_edge[e] <= m.q_edge_b[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.q_edge_b[e])
            m.constraint.add(-1 * M * m.x_edge[e] <= m.q_edge_c[e])
            m.constraint.add( 1 * M * m.x_edge[e] >= m.q_edge_c[e])
            if edge_list[e].Phases[0] == 0:
                m.constraint.add(0 == m.p_edge_a[e])
                m.constraint.add(0 == m.q_edge_a[e])
            if edge_list[e].Phases[1] == 0:
                m.constraint.add(0 == m.p_edge_b[e])
                m.constraint.add(0 == m.q_edge_b[e])
            if edge_list[e].Phases[2] == 0:
                m.constraint.add(0 == m.p_edge_c[e])
                m.constraint.add(0 == m.q_edge_c[e])
        for g in range(n_gen):
            m.constraint.add(0 <= m.p_gen_a[g])
            m.constraint.add(0 <= m.p_gen_b[g])
            m.constraint.add(0 <= m.p_gen_c[g])
            m.constraint.add(M * m.x_gen[g] >= m.p_gen_a[g])
            m.constraint.add(M * m.x_gen[g] >= m.p_gen_b[g])
            m.constraint.add(M * m.x_gen[g] >= m.p_gen_c[g])
            m.constraint.add(-1 * M * m.x_gen[g] <= m.q_gen_a[g])
            m.constraint.add(-1 * M * m.x_gen[g] <= m.q_gen_b[g])
            m.constraint.add(-1 * M * m.x_gen[g] <= m.q_gen_c[g])
            m.constraint.add(M * m.x_gen[g] >= m.q_gen_a[g])
            m.constraint.add(M * m.x_gen[g] >= m.q_gen_b[g])
            m.constraint.add(M * m.x_gen[g] >= m.q_gen_c[g])

        # power balance
        PF_BAL_PA = np.dot(m.p_gen_a, I_gen_bus) + \
                    np.dot(m.p_load_a, I_loadcap_bus) + \
                    np.dot(m.p_edge_a, I_edge_bus)
        PF_BAL_PB = np.dot(m.p_gen_b, I_gen_bus) + \
                    np.dot(m.p_load_b, I_loadcap_bus) + \
                    np.dot(m.p_edge_b, I_edge_bus)
        PF_BAL_PC = np.dot(m.p_gen_c, I_gen_bus) + \
                    np.dot(m.p_load_c, I_loadcap_bus) + \
                    np.dot(m.p_edge_c, I_edge_bus)
        PF_BAL_QA = np.dot(m.q_gen_a, I_gen_bus) + \
                    np.dot(m.q_load_a, I_loadcap_bus) + \
                    np.dot(m.q_edge_a, I_edge_bus)
        PF_BAL_QB = np.dot(m.q_gen_b, I_gen_bus) + \
                    np.dot(m.q_load_b, I_loadcap_bus) + \
                    np.dot(m.q_edge_b, I_edge_bus)
        PF_BAL_QC = np.dot(m.q_gen_c, I_gen_bus) + \
                    np.dot(m.q_load_c, I_loadcap_bus) + \
                    np.dot(m.q_edge_c, I_edge_bus)

        for n in range(n_bus):
            m.constraint.add( PF_BAL_PA[n] == 0 )
            m.constraint.add( PF_BAL_PB[n] == 0 )
            m.constraint.add( PF_BAL_PC[n] == 0 )
            m.constraint.add( PF_BAL_QA[n] == 0 )
            m.constraint.add( PF_BAL_QB[n] == 0 )
            m.constraint.add( PF_BAL_QC[n] == 0 )

    # Power flow voltage constraints
    if st_pf_voltage == True:
        for e in range(n_edge):
            Nfrom = edge_list[e].BusNum_A
            Nto   = edge_list[e].BusNum_B
            Phases= edge_list[e].Phases
            if edge_list[e].Device != 'regulator':
                ZZ = np.diag(Phases).dot(Z[(3 * Nfrom):(3 * (Nfrom+1)), (3 * Nto):(3 * (Nto+1))])
                RR = np.real(ZZ)
                XX = np.imag(ZZ)

                a  = np.array([[1], [np.exp(1j*-120/180*np.pi)], [np.exp(1j*+120/180*np.pi)]])
                Requal = np.real(a.dot(a.conj().T)) * RR + np.imag(a.dot(a.conj().T)) * XX
                Xequal = np.real(a.dot(a.conj().T)) * XX + np.imag(a.dot(a.conj().T)) * RR

                Vf2 = np.array([m.v_bus_a2[Nfrom], m.v_bus_b2[Nfrom], m.v_bus_c2[Nfrom]])
                Vt2 = np.array([m.v_bus_a2[Nto],   m.v_bus_b2[Nto],   m.v_bus_c2[Nto]])

                Pft = np.array([m.p_edge_a[e], m.p_edge_b[e], m.p_edge_c[e]])
                Qft = np.array([m.q_edge_a[e], m.q_edge_b[e], m.q_edge_c[e]])

                # Put 3x1 vector on the left then + a 1x1 number. Reverse it will not work!
                Left  =  -(Vt2 - Vf2 + 2 * (np.dot(Requal, Pft) + np.dot(Xequal,Qft))) - (M * (1 - m.x_edge[e]))
                Right =   (Vt2 - Vf2 + 2 * (np.dot(Requal, Pft) + np.dot(Xequal,Qft))) - (M * (1 - m.x_edge[e]))

                for c in Left:
                    m.constraint.add(c <= 0 )
                for c in Right:
                    m.constraint.add(c <= 0 )

                m.constraint.add(m.rg_a2[e] == 0)
                m.constraint.add(m.rg_b2[e] == 0)
                m.constraint.add(m.rg_c2[e] == 0)
                
            else:
                if st_regulator_control == False:
                    m.constraint.add(m.rg_a2[e] == edge_list[e].Position[0]**2)
                    m.constraint.add(m.rg_b2[e] == edge_list[e].Position[1]**2)
                    m.constraint.add(m.rg_c2[e] == edge_list[e].Position[2]**2)
                    if edge_list[e].Phases[0] == 1:
                        m.constraint.add( -(m.v_bus_a2[Nto] - edge_list[e].Position[0]**2 * m.v_bus_a2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_a2[Nto] - edge_list[e].Position[0]**2 * m.v_bus_a2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                    if edge_list[e].Phases[1] == 1:
                        m.constraint.add( -(m.v_bus_b2[Nto] - edge_list[e].Position[1]**2 * m.v_bus_b2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_b2[Nto] - edge_list[e].Position[1]**2 * m.v_bus_b2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                    if edge_list[e].Phases[2] == 1:
                        m.constraint.add( -(m.v_bus_c2[Nto] - edge_list[e].Position[2]**2 * m.v_bus_c2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_c2[Nto] - edge_list[e].Position[2]**2 * m.v_bus_c2[Nfrom]) - (M*(1-m.x_edge[e])) <= 0 )
                else:
                    if edge_list[e].Phases[0] == 1:
                        m.constraint.add( -(m.v_bus_a2[Nto] - (m.rg_a2[e]+1.01*m.v_bus_a2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_a2[Nto] - (m.rg_a2[e]+1.01*m.v_bus_a2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(m.rg_a2[e] >= 0.81 )
                        m.constraint.add(m.rg_a2[e] <= 1.21 )
                    else:
                        m.constraint.add(m.rg_a2[e] == 0 )
                    if edge_list[e].Phases[1] == 1:
                        m.constraint.add( -(m.v_bus_b2[Nto] - (m.rg_b2[e]+1.01*m.v_bus_b2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_b2[Nto] - (m.rg_b2[e]+1.01*m.v_bus_b2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(m.rg_b2[e] >= 0.81 )
                        m.constraint.add(m.rg_b2[e] <= 1.21 )
                    else:
                        m.constraint.add(m.rg_b2[e] == 0 )
                    if edge_list[e].Phases[2] == 1:
                        m.constraint.add( -(m.v_bus_c2[Nto] - (m.rg_c2[e]+1.01*m.v_bus_c2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(  (m.v_bus_c2[Nto] - (m.rg_c2[e]+1.01*m.v_bus_c2[Nfrom]-1.01)) - (M*(1-m.x_edge[e])) <= 0 )
                        m.constraint.add(m.rg_c2[e] >= 0.81 )
                        m.constraint.add(m.rg_c2[e] <= 1.21 )
                    else:
                        m.constraint.add(m.rg_c2[e] == 0 )

    # Voltage limit constraints
    if st_voltage_limit == True:
        # substation voltage
        for g in range(n_gen):
            if gen_list[g]['GridForming'] == True:  # if the generator is controllable
                if gen_list[g]['Outage'] == False:  # if the generator is damaged.
                    if gen_list[g]['Device'] == 'substation':
                        m.constraint.add(m.v_bus_a2[gen_list[g]['BusNum']] == m.x_bus[gen_list[g]['BusNum']] * iso_Vsrc**2)
                        m.constraint.add(m.v_bus_b2[gen_list[g]['BusNum']] == m.x_bus[gen_list[g]['BusNum']] * iso_Vsrc**2)
                        m.constraint.add(m.v_bus_c2[gen_list[g]['BusNum']] == m.x_bus[gen_list[g]['BusNum']] * iso_Vsrc**2)
        # voltage limit
        for n in range(n_bus):
            m.constraint.add(m.v_bus_a2[n] >= iso_vmin * iso_vmin * m.x_bus[n] )
            m.constraint.add(m.v_bus_a2[n] <= iso_vmax * iso_vmax * m.x_bus[n] )
            m.constraint.add(m.v_bus_b2[n] >= iso_vmin * iso_vmin * m.x_bus[n] )
            m.constraint.add(m.v_bus_b2[n] <= iso_vmax * iso_vmax * m.x_bus[n] )
            m.constraint.add(m.v_bus_c2[n] >= iso_vmin * iso_vmin * m.x_bus[n] )
            m.constraint.add(m.v_bus_c2[n] <= iso_vmax * iso_vmax * m.x_bus[n] )

    # Line capacity constraints
    if st_line_capacity == True:
        for e in range(n_edge):
            plg = edge_list[e].EmergAmps * 1000 * ( Vbase / np.sqrt(3) ) / Sbase / np.sqrt(2)
            m.constraint.add(-1 * plg - m.p_edge_a[e] <= 0)
            m.constraint.add(     plg - m.p_edge_a[e] >= 0)
            m.constraint.add(-1 * plg - m.p_edge_b[e] <= 0)
            m.constraint.add(     plg - m.p_edge_b[e] >= 0)
            m.constraint.add(-1 * plg - m.p_edge_c[e] <= 0)
            m.constraint.add(     plg - m.p_edge_c[e] >= 0)
            m.constraint.add(-1 * plg - m.q_edge_a[e] <= 0)
            m.constraint.add(     plg - m.q_edge_a[e] >= 0)
            m.constraint.add(-1 * plg - m.q_edge_b[e] <= 0)
            m.constraint.add(     plg - m.q_edge_b[e] >= 0)
            m.constraint.add(-1 * plg - m.q_edge_c[e] <= 0)
            m.constraint.add(     plg - m.q_edge_c[e] >= 0)

    # Generator capacity constraints
    if st_gen_capacity == True:
        for g in range(n_gen):
            m.constraint.add( 1 * gen_list[g]['Pmin'] * m.x_gen[g] - (m.p_gen_a[g] + m.p_gen_b[g] + m.p_gen_c[g]) <= 0)
            m.constraint.add(-1 * gen_list[g]['Pmax'] * m.x_gen[g] + (m.p_gen_a[g] + m.p_gen_b[g] + m.p_gen_c[g]) <= 0)
            m.constraint.add( 1 * gen_list[g]['Qmin'] * m.x_gen[g] - (m.q_gen_a[g] + m.q_gen_b[g] + m.q_gen_c[g]) <= 0)
            m.constraint.add(-1 * gen_list[g]['Qmax'] * m.x_gen[g] + (m.q_gen_a[g] + m.q_gen_b[g] + m.q_gen_c[g]) <= 0)

            if gen_list[g]['GridForming'] == True:  # if the generator is controllable
                if gen_list[g]['Outage'] == False:  # if the generator is damaged.
                    if gen_list[g]['Device'] != 'substation':
                        m.constraint.add( m.p_gen_a[g] == m.p_gen_b[g] )
                        m.constraint.add( m.p_gen_a[g] == m.p_gen_c[g] )
                        m.constraint.add( m.q_gen_a[g] == m.q_gen_b[g] )
                        m.constraint.add( m.q_gen_a[g] == m.q_gen_c[g] )

    # # load model
    for l in range(n_loadcap):
        if loadcap_list[l]['Device'] == 'load':
            m.constraint.add(m.p_load_a[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][0] == 0)
            m.constraint.add(m.p_load_b[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][1] == 0)
            m.constraint.add(m.p_load_c[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][2] == 0)
            m.constraint.add(m.q_load_a[l] - m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][0] == 0)
            m.constraint.add(m.q_load_b[l] - m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][1] == 0)
            m.constraint.add(m.q_load_c[l] - m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][2] == 0)
        else:
            m.constraint.add(m.p_load_a[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][0] == 0)
            m.constraint.add(m.p_load_b[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][1] == 0)
            m.constraint.add(m.p_load_c[l] - m.x_loadcap[l] * loadcap_list[l]['kW']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][2] == 0)
            m.constraint.add(m.q_load_a[l] + m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][0] == 0)
            m.constraint.add(m.q_load_b[l] + m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][1] == 0)
            m.constraint.add(m.q_load_c[l] + m.x_loadcap[l] * loadcap_list[l]['kvar']/sum(loadcap_list[l]['Phases'])*loadcap_list[l]['Phases'][2] == 0)

    # Topology
    # Substation
    for g in range(n_gen):
        m.constraint.add( m.x_gen[g] <= m.x_bus[gen_list[g]['BusNum']] )
        if gen_list[g]['Outage'] == False:
            if gen_list[g]['GridForming'] == True:
                m.constraint.add( m.x_gen[g] == 1 )
        else:
            m.constraint.add( m.x_gen[g] == 0 )
            m.constraint.add( m.x_gen[g] == m.x_bus[gen_list[g]['BusNum']] )
    # load
    for l in range(n_loadcap):
        if loadcap_list[l]['Outage'] == False:
            # if loadcap_list[l]['Control'] == True:
            #     if loadcap_list[l]['Device'] == 'capacitor' and st_cap_control == False:
            #         m.constraint.add( m.x_loadcap[l] == m.x_bus[loadcap_list[l]['BusNum']] )
            #     else:
            #         m.constraint.add( m.x_loadcap[l] <= m.x_bus[loadcap_list[l]['BusNum']] )
            # else:
            #     m.constraint.add( m.x_loadcap[l] == m.x_bus[loadcap_list[l]['BusNum']] )

            m.constraint.add( m.x_loadcap[l] == m.x_bus[loadcap_list[l]['BusNum']] )
        else:
            m.constraint.add( m.x_loadcap[l] == 0 )
    # edge
    for e in range(n_edge):
        if edge_list[e].Outage == False:
            if edge_list[e].Control == True:
                m.constraint.add( m.x_edge[e] <= m.x_bus[edge_list[e].BusNum_A] )
                m.constraint.add( m.x_edge[e] <= m.x_bus[edge_list[e].BusNum_B] )
            else:
                m.constraint.add( m.x_edge[e] == m.x_bus[edge_list[e].BusNum_A] )
                m.constraint.add( m.x_edge[e] == m.x_bus[edge_list[e].BusNum_B] )
        if edge_list[e].Outage == True:
            m.constraint.add( m.x_edge[e] == 0 )
            m.constraint.add( m.x_bus[edge_list[e].BusNum_A] == 0 )
            m.constraint.add( m.x_bus[edge_list[e].BusNum_B] == 0 )


    # only allow the switches to open for isolation purpose
    # Note for isolation function, no radial topology constraints need to be included
    # as the system was originally operated in radial topology
    for s in range(n_switch):
        m.constraint.add( m.x_edge[swi_idx[s]] <= swi_status[s] )

    # objective function
    total_bus = 0
    for i in range(n_bus):
        total_bus += m.x_bus[i]
    m.cost = pyo.Objective(expr=-total_bus, sense=pyo.minimize)

    # pyo.SolverFactory('gurobi', solver_io="python").solve(m).write()   # Seems Gurobi solver cannot be identified in PyCharm
    pyo.SolverFactory('glpk').solve(m).write() # GLPK solver works just fine



    #####################################################
    ## print switch operation decisions
    print('\n')
    print('\nSwitch Positon:')
    for s in range(n_switch):
        print(edge_list[swi_idx[s]].Bus_A, edge_list[swi_idx[s]].Bus_B, m.x_edge[swi_idx[s]]())
    print('\nRegulator Positon:')
    for e in range(n_edge):
        if edge_list[e].Device == 'regulator':
            print(edge_list[e].Name, np.sqrt(m.rg_a2[e]()), np.sqrt(m.rg_b2[e]()), np.sqrt(m.rg_c2[e]()))
    print('\nCapacitor Position:')
    for c in range(n_loadcap):
        if loadcap_list[c]['Device'] == 'capacitor':
            print(loadcap_list[c]['Name'], m.q_load_a[c](), m.q_load_b[c](), m.q_load_c[c]())

    ## plot figures
    # get voltage data
    v_abc = np.zeros((n_bus,3))
    for i in range(n_bus):
        if bus_list[i]['Phases'][0] == 1:
            v_abc[i,0] = np.sqrt(np.array(m.v_bus_a2[i]()))
        if bus_list[i]['Phases'][1] == 1:
            v_abc[i,1] = np.sqrt(np.array(m.v_bus_b2[i]()))
        if bus_list[i]['Phases'][2] == 1:
            v_abc[i,2] = np.sqrt(np.array(m.v_bus_c2[i]()))

    # create two figures
    fig, (axx,ax) = plt.subplots(2)

    # plot voltage profile
    axx.axis([-5, 140, 0.8, 1.1])
    axx.grid()
    axx.plot(v_abc[:,0], 'o', color = 'r' )
    axx.plot(v_abc[:,1], 'o', color = 'b' )
    axx.plot(v_abc[:,2], 'o', color = 'g' )

    # Plot the routes on a map
    coord_x = [bus_list[i]['Coord_X'] for i in range(n_bus)]
    coord_y = [bus_list[i]['Coord_Y'] for i in range(n_bus)]
    name_list = [bus_list[i]['Name'] for i in range(n_bus)]

    ax.scatter(coord_x, coord_y, s=50, color='blue', marker='o')
    for i, txt in enumerate(name_list):
        ax.annotate(txt, (coord_x[i]-1, coord_y[i]-1))

    gen_idx = [gen_list[i]['Number'] for i in range(n_gen)]
    for i in range(n_gen):
        if gen_list[i]['Outage'] == False:
            ax.scatter(coord_x[gen_idx[i]], coord_y[gen_idx[i]], s=100, color='green', marker='o')
        if gen_list[i]['Outage'] == True:
            ax.scatter(coord_x[gen_idx[i]], coord_y[gen_idx[i]], s=100, color='red', marker='o')

    coord_ex1 = [bus_list[edge_list[i].BusNum_A]['Coord_X'] for i in range(n_edge)]
    coord_ey1 = [bus_list[edge_list[i].BusNum_A]['Coord_Y'] for i in range(n_edge)]
    coord_ex2 = [bus_list[edge_list[i].BusNum_B]['Coord_X'] for i in range(n_edge)]
    coord_ey2 = [bus_list[edge_list[i].BusNum_B]['Coord_Y'] for i in range(n_edge)]

    for i in range(n_edge):
        if m.x_edge[i]() == 1:
            if edge_list[i].Device != 'switch':
                plt.plot([coord_ex1[i], coord_ex2[i]], [coord_ey1[i], coord_ey2[i]], color = 'b')
            else:
                plt.plot([coord_ex1[i], coord_ex2[i]], [coord_ey1[i], coord_ey2[i]], color = 'g')
        else:
            if edge_list[i].Control != 'switch':
                plt.plot([coord_ex1[i], coord_ex2[i]], [coord_ey1[i], coord_ey2[i]], color = 'b', linestyle="--")
            else:
                plt.plot([coord_ex1[i], coord_ex2[i]], [coord_ey1[i], coord_ey2[i]], color = 'g', linestyle="--")
        if edge_list[i].Outage == True:
            plt.plot([coord_ex1[i], coord_ex2[i]], [coord_ey1[i], coord_ey2[i]], color = 'r', linestyle="--")
    plt.show()


    #####
    print('\n')
    print('Gen Power:', m.p_gen_a[0](),m.p_gen_b[0](),m.p_gen_c[0]())

    p_a = 0
    p_b = 0
    p_c = 0
    for i in range(n_loadcap):
        p_a += m.p_load_a[i]()
        p_b += m.p_load_b[i]()
        p_c += m.p_load_c[i]()
    print('Load Demand:', p_a, p_b, p_c)


