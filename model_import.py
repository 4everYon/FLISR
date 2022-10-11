'''
Project: SETO Open Energy Data Initiative (OEDI)
Author: Bo Chen @ Argonne National Lab
Date: 08/30/2022
Description: Use OpenDSS COM interface to import data from existing OpenDSS test systems

Note:
Use OpenDSSSdirect package from https://github.com/NREL/OpenDSSDirect.py
Based on OpenDSSDirect 0.6.0  [0.7.0 cannot be installed since dss-python 0.12.0 is not available but required by 0.7.0]
Need to install Pandas
Make sure OpenDSS test feeder is under a path without space or special characters, otherwise cannot be found by OpenDSS COM
'''

import os
import opendssdirect as dss
import math
import copy
import pickle

### import Bus data ###########################################################
def func_import_bus(dss):
    Bus = dss.Bus
    n_bus = dss.Circuit.NumBuses()

    busData = {}
    if n_bus == 0: # no bus defined, then return an empty array
        pass
    else: # if any buses, read data
        busData['AllBusNames'] = dss.Circuit.AllBusNames()
        busData['BusNumbers'] = list(range(n_bus))
        busData['NumBuses'] = n_bus

        bus_number = -1
        for bus_name in dss.Circuit.AllBusNames():
            ### note this code is to 'active' the bus you want to get the info from
            dss.Circuit.SetActiveBus(bus_name)
            ### read data
            bus_number += 1
            busData[bus_name] = {} # 2-D dict
            busData[bus_name]['Number'] = bus_number
            busData[bus_name]['Name'] = Bus.Name()
            busData[bus_name]['NumPhases'] = Bus.NumNodes()
            busData[bus_name]['Phases'] = Bus.Nodes() # here we use phases to represent the node concept in OpenDSS
            busData[bus_name]['kVBase'] = Bus.kVBase()
            busData[bus_name]['Coorddefined'] = Bus.Coorddefined()
            busData[bus_name]['Coord_X'] = Bus.X()
            busData[bus_name]['Coord_Y'] = Bus.Y()
            busData[bus_name]['Outage'] = False


    return busData

### import load data ##########################################################
def func_import_load(dss):
    Load = dss.Loads
    n_load = len(Load.AllNames())

    loadData = {}
    if n_load == 0:
        return loadData
    else:
        loadData['AllLoadNames'] = Load.AllNames()
        loadData['LoadNumbers'] = list(range(n_load))
        loadData['NumLoads'] = n_load

        load_number = -1
        for load_name in dss.Loads.AllNames():
            ## set a load active by name
            dss.Loads.Name(load_name)
            ## import data
            load_number += 1
            loadData[load_name] = {} # 2-D dict
            loadData[load_name]['Number'] = load_number
            loadData[load_name]['Name'] = Load.Name()
            loadData[load_name]['AllocationFactor'] = Load.AllocationFactor()  ## may add more properties to use it
            loadData[load_name]['kV_LL'] = Load.kV() ## if NumPhase = 1, it's LN, otherwise, it's LL.
            loadData[load_name]['kVABase'] = Load.kVABase()
            loadData[load_name]['kW'] = Load.kW()  ## total P on all available phases
            loadData[load_name]['kvar'] = Load.kvar() ## total Q on all available phases
            loadData[load_name]['Model'] = Load.Model() ##
            loadData[load_name]['IsDelta'] = Load.IsDelta()


            loadData[load_name]['Outage'] = False
            loadData[load_name]['Status'] = True
            loadData[load_name]['Control'] = False

            ## get info from CkeElement
            bus_phase = dss.CktElement.BusNames()[0].split('.') ## assume load only connects to ONE bus.
            BusName = bus_phase[0]
            if len(bus_phase) == 1: # means it's a 3-phase load
                Phases = [1,2,3]
            else:
                Phases = [int(i) for i in bus_phase[1:]]
            loadData[load_name]['BusName'] = BusName
            loadData[load_name]['NumPhases'] = len(Phases)  ## dss.CktElement.NumPhases() is always 1 ????
            loadData[load_name]['Phases'] = Phases

            if len(Phases) == 1:
                loadData[load_name]['kV_LN'] = Load.kV()
            else:
                loadData[load_name]['kV_LN'] = Load.kV()/math.sqrt(3)


    return loadData

### import line data ##########################################################
def func_import_line(dss, busData):
    Line = dss.Lines
    n_line = len(dss.Lines.AllNames())

    lineData = {}
    if n_line == 0:
        pass
    else:
        lineData['AllLineNames'] = Line.AllNames()
        lineData['LineNumbers'] = list(range(n_line))
        lineData['NumLines'] = n_line

        line_number = -1
        for line_name in Line.AllNames():
            ## set a line active by name
            Line.Name(line_name)
            ## import data
            line_number += 1
            lineData[line_name] = {}  # 2-D dict
            lineData[line_name]['Name'] = Line.Name()
            lineData[line_name]['Number'] = line_number
            lineData[line_name]['Bus1'] = Line.Bus1().split('.')[0]
            lineData[line_name]['Bus2'] = Line.Bus2().split('.')[0]
            lineData[line_name]['Enabled'] = dss.CktElement.Enabled()
            lineData[line_name]['NormAmps'] = Line.NormAmps()
            lineData[line_name]['EmergAmps'] = Line.EmergAmps()
            lineData[line_name]['Outage'] = False # Default to be false

            bus_phase = Line.Bus1().split('.') ## assume Bus1 and Bus2 phasing is the same for a Line object
            if len(bus_phase) == 1: # means it's a 3-phase load
                Phases = [1,2,3]
            else:
                Phases = [int(i) for i in bus_phase[1:]]
            lineData[line_name]['NumPhases'] = Line.Phases()  # Line.Phases() is different from Load.Phases()
            lineData[line_name]['Phases'] = Phases

            lineData[line_name]['Length'] = Line.Length()
            lineData[line_name]['Units'] = Line.Units()  # Units= {1:mi, 2:km, 3:kft, 4:m, 5:ft, 6:in, 7:cm}

            ## in MATLAB, if 'LineCode' data is available, Rmatrix is per unit, if 'Geometry' is available, Rmatrix is total value, if both available: to be checked
            if not Line.LineCode():
                pass
            else:
                dss.LineCodes.Name( Line.LineCode() )  # set the right LineCode to active
                lineData[line_name]['LineCode'] = Line.LineCode()
                lineData[line_name]['LineCodeLengthUnits'] = dss.LineCodes.Units()  # Units= {1:mi, 2:km, 3:kft, 4:m, 5:ft, 6:in, 7:cm}

                lineData[line_name]['RmatrixPerLengthUnit'] = dss.LineCodes.Rmatrix()
                lineData[line_name]['XmatrixPerLengthUnit'] = dss.LineCodes.Xmatrix()
                lineData[line_name]['CmatrixPerLengthUnit'] = dss.LineCodes.Cmatrix()

            if not Line.Geometry():
                pass
            else:
                lineData[line_name]['LineGeometry'] = Line.Geometry()
                ## in MATLAB, if 'LineCode' data is available, Rmatrix is per unit, if 'Geometry' is available, Rmatrix is ohms value, if both available: to be checked
                lineData[line_name]['Rmatrix'] = Line.RMatrix()
                lineData[line_name]['Xmatrix'] = Line.XMatrix()
                lineData[line_name]['Cmatrix'] = Line.CMatrix()

            if not Line.LineCode() or not Line.Geometry():
                lineData[line_name]['RmatrixPerLengthUnit'] = Line.RMatrix()
                lineData[line_name]['XmatrixPerLengthUnit'] = Line.XMatrix()
                lineData[line_name]['CmatrixPerLengthUnit'] = Line.CMatrix()
                # print('Found a line with neither LineCode nor LineGeometry data: ' + line_name)
            else:
                print('Found a line with both LineCode and LineGeometry data: ' + line_name)

            ## if a line is disabled, it's end buses may be not included in Bus data: Update Missing Bus
            if lineData[line_name]['Enabled'] == False:
                # check if there are any missing buses
                if lineData[line_name]['Bus1'] in busData.keys():
                    bus1_exist = True
                else:
                    bus1_exist = False
                if lineData[line_name]['Bus2'] in busData.keys():
                    bus2_exist = True
                else:
                    bus2_exist = False
                # update Bus to include missing buses
                if bus1_exist == True and bus2_exist == False:
                    # use DEEEEEEEEEEEEEEEEEp copy!!!
                    busData.update( {lineData[line_name]['Bus2'] : copy.deepcopy(busData[ lineData[line_name]['Bus1'] ])} )
                    # still need to change something..
                    # update Bus information
                    busData['AllBusNames'].append(lineData[line_name]['Bus2'])
                    n_bus = busData['NumBuses'] + 1
                    busData['NumBuses'] = n_bus
                    busData['BusNumbers'].append(n_bus - 1) # Python count a list start from 0
                    # Bus 2 has same info with Bus 1, must correct
                    busData[ lineData[line_name]['Bus2'] ]['Number'] = n_bus - 1
                    busData[ lineData[line_name]['Bus2'] ]['Name']   = lineData[line_name]['Bus2']
                    # assume all other Bus properties are the same...

                elif bus1_exist == False and bus2_exist == True:
                    # use DEEEEEEEEEEEEEEEEEp copy!!!
                    busData.update({lineData[line_name]['Bus1']: copy.deepcopy(busData[ lineData[line_name]['Bus2'] ])})

                    busData['AllBusNames'].append(lineData[line_name]['Bus1'])
                    n_bus = busData['NumBuses'] + 1
                    busData['NumBuses'] = n_bus
                    busData['BusNumbers'].append(n_bus - 1)  # Python count a list start from 0
                    # Bus 2 has same info with Bus 1, must correct
                    busData[ lineData[line_name]['Bus1'] ]['Number'] = n_bus - 1
                    busData[ lineData[line_name]['Bus1'] ]['Name']   = lineData[line_name]['Bus1']
                    # assume all other Bus properties are the same...

                elif bus1_exist == False and bus2_exist == False:
                    print('A disabled line is found to have two missing buses: ' + line_name)


    ## check if any reactors
    idx = dss.PDElements.First()
    while idx == 1:
        r_name = dss.PDElements.Name().split('.')[0]
        if r_name == 'Reactor':
            ## update lineData
            line_name = dss.PDElements.Name().split('.')[1]

            line_number += 1
            lineData[line_name] = {}  # 2-D dict
            lineData[line_name]['Name'] = line_name
            lineData[line_name]['Number'] = line_number
            lineData[line_name]['Bus1'] = dss.CktElement.BusNames()[0].split('.')[0]
            lineData[line_name]['Bus2'] = dss.CktElement.BusNames()[1].split('.')[0]
            lineData[line_name]['Enabled'] = dss.CktElement.Enabled()
            lineData[line_name]['NormAmps'] = Line.NormAmps()
            lineData[line_name]['EmergAmps'] = Line.EmergAmps()
            lineData[line_name]['Outage'] = False

            bus_phase = dss.CktElement.Name()[0].split('.') ## assume Bus1 and Bus2 phasing is the same for a Line object
            if len(bus_phase) == 1:
                Phases = [1,2,3]
            else:
                Phases = [int(i) for i in bus_phase[1:]]
            lineData[line_name]['NumPhases'] = len(Phases)
            lineData[line_name]['Phases'] = Phases

            lineData[line_name]['Length'] = 0.001 ##############################
            lineData[line_name]['Units'] = 4      ######### Units= {1:mi, 2:km, 3:kft, 4:m, 5:ft, 6:in, 7:cm}

            lineData[line_name]['Rmatrix'] = [0, 0, 0, 0, 0, 0]
            lineData[line_name]['Xmatrix'] = [0, 0, 0, 0, 0, 0]
            lineData[line_name]['Cmatrix'] = [0, 0, 0, 0, 0, 0]


        idx = dss.PDElements.Next()



    return lineData

### import transformer data ###################################################
def func_import_transformer(dss):
    Txm = dss.Transformers
    n_txm = len(Txm.AllNames())

    txmData = {}
    if n_txm == 0:
        pass
    else:
        txmData['AllTransNames'] = Txm.AllNames()
        txmData['TransNumbers'] = list(range(n_txm))
        txmData['NumTrans'] = n_txm

        txm_number = -1
        for txm_name in Txm.AllNames():
            ## set a txm active by name
            Txm.Name(txm_name)
            ## import data
            txm_number += 1
            txmData[txm_name] = {}  # 2-D dict
            txmData[txm_name]['Name'] = Txm.Name()
            txmData[txm_name]['Number'] = txm_number
            txmData[txm_name]['Enabled'] = dss.CktElement.Enabled()
            txmData[txm_name]['IsDelta'] = Txm.IsDelta()
            txmData[txm_name]['kV'] = Txm.kV()
            txmData[txm_name]['kVA'] = Txm.kVA()
            txmData[txm_name]['MinTap'] = Txm.MinTap()
            txmData[txm_name]['MaxTap'] = Txm.MaxTap()
            txmData[txm_name]['NumTaps'] = Txm.NumTaps()
            txmData[txm_name]['NumWindings'] = Txm.NumWindings()
            txmData[txm_name]['R'] = Txm.R()
            txmData[txm_name]['Xhl'] = Txm.Xhl()
            txmData[txm_name]['NormalAmps'] = dss.CktElement.NormalAmps()
            txmData[txm_name]['EmergAmps'] = dss.CktElement.EmergAmps()
            txmData[txm_name]['Outage'] = False

            #### some test feeder has 3-winding load transformer... (from IEEE8500 test feeder)
            # Definitions of the 1-phase "center tapped" 120 / 240V service transformers.
            # These transformers are defined as 3-winding transformers(as they are in reality)
            # The primary winding is connected 7200V line-to-neutral on one of the phases.
            # The secondary windings are consistently connected 1.0 and 0.2 respectively to get the polarity correct.

            # Triplex lines connecting the secondary of the service transformers(MV/LV) to 120/240V loada
            # These lines are defined(see Triplex_linecodes.dss) as 2 - phase lines
            # See transformer connection for the meaning of nodes 1 and 2 at the secondary buses.

            # Loads are defined as two 120V loads connected line to neutral
            # For phases > 1, OpenDSS Load model assumes that a 2-phase load is 2 phases of a 3-phase system
            # and requires the voltage base to be specified same as 3 - phase(Line-Line kV = Line-Neutral * sqrt(3))
            # Thus, the base voltage is 208V to get 120V line-to-neutral.
            # Alternatively, we could have defined two separate 1-phase loads rated at 120V.
            # The service transformer connection enforces the voltage polarities and phase angles.
            # The kW load is divided equally between the two "phases"

            # if Txm.NumWindings() == 3:
            #     a=1

            bus_no = 0
            for BUS in dss.CktElement.BusNames():
                bus_phase = BUS.split('.')
                if len(bus_phase) == 1: # means it's a 3-phase load
                    Phases = [1,2,3]
                else:
                    Phases = [int(i) for i in bus_phase[1:]]
                bus_no += 1
                txmData[txm_name]['Bus' + str(bus_no)] = BUS.split('.')[0]
                txmData[txm_name]['Bus' + str(bus_no) + 'Phases'] = Phases

    return txmData

### import capacitor data #####################################################
def func_import_capacitor(dss):
    Cap = dss.Capacitors
    n_cap = len(Cap.AllNames())

    capData = {}
    if n_cap == 0:
        return capData
    else:
        capData['AllCapacitorNames'] = Cap.AllNames()
        capData['CapNumbers'] = list(range(n_cap))
        capData['NumCaps'] = n_cap

        cap_number = -1
        for cap_name in Cap.AllNames():
            ## set a load active by name
            dss.Capacitors.Name(cap_name)
            ## import data
            cap_number += 1
            capData[cap_name] = {} # 2-D dict
            capData[cap_name]['Number'] = cap_number
            capData[cap_name]['Name'] = Cap.Name()
            capData[cap_name]['kV_LL'] = Cap.kV()
            capData[cap_name]['IsDelta'] = Cap.IsDelta()
            capData[cap_name]['kvar'] = Cap.kvar()
            capData[cap_name]['AvailableSteps'] = Cap.AvailableSteps()
            capData[cap_name]['NumSteps'] = Cap.NumSteps()
            capData[cap_name]['States'] = Cap.States()
            capData[cap_name]['Outage'] = False
            capData[cap_name]['Control'] = True

            capData[cap_name]['kW'] = 0

            ## get info from CkeElement
            bus_phase = dss.CktElement.BusNames()[0].split('.') ## assume load only connects to ONE bus.
            BusName = bus_phase[0]
            if len(bus_phase) == 1: # means it's a 3-phase load
                Phases = [1,2,3]
            else:
                Phases = [int(i) for i in bus_phase[1:]]
            capData[cap_name]['BusName'] = BusName
            capData[cap_name]['NumPhases'] = len(Phases)  ## dss.CktElement.NumPhases() is always 1 ????
            capData[cap_name]['Phases'] = Phases

            if len(Phases) == 1:
                capData[cap_name]['kV_LN'] = Cap.kV()
            else:
                capData[cap_name]['kV_LN'] = Cap.kV()/math.sqrt(3)

    return capData

### import regulator data #####################################################
def func_import_regulator(dss):
    Reg = dss.RegControls
    n_reg = len(Reg.AllNames())

    regData = {}
    if n_reg == 0:
        pass
    else:
        regData['AllTransNames'] = Reg.AllNames()
        regData['TransNumbers'] = list(range(n_reg))
        regData['NumTrans'] = n_reg

        reg_number = -1
        for reg_name in Reg.AllNames():
            ## set a txm active by name
            Reg.Name(reg_name)
            ## import data
            reg_number += 1
            regData[reg_name] = {}  # 2-D dict
            regData[reg_name]['Name'] = Reg.Name()
            regData[reg_name]['Number'] = reg_number
            regData[reg_name]['Enabled'] = dss.CktElement.Enabled()

            regData[reg_name]['Transformer'] = Reg.Transformer()
            regData[reg_name]['Winding'] = Reg.Winding()

            regData[reg_name]['TapWinding'] = Reg.TapWinding()
            regData[reg_name]['Delay'] = Reg.Delay()
            regData[reg_name]['TapDelay'] = Reg.TapDelay()
            regData[reg_name]['TapNumber'] = Reg.TapNumber()
            regData[reg_name]['VoltageLimit'] = Reg.VoltageLimit()

            regData[reg_name]['PTRatio'] = Reg.PTRatio()
            regData[reg_name]['CTPrimary'] = Reg.CTPrimary()
            regData[reg_name]['ForwardR'] = Reg.ForwardR()
            regData[reg_name]['ForwardX'] = Reg.ForwardX()
            regData[reg_name]['ForwardBand'] = Reg.ForwardBand()
            regData[reg_name]['ForwardVreg'] = Reg.ForwardVreg()
            regData[reg_name]['ReverseR'] = Reg.ReverseR()
            regData[reg_name]['ReverseX'] = Reg.ReverseX()
            regData[reg_name]['ReverseBand'] = Reg.ReverseBand()
            regData[reg_name]['ReverseVreg'] = Reg.ReverseVreg()

    return regData

### get substation data #####################################################
def func_get_substation(dss, Bus, Line, Transformer, Substation_names):

    Sub = Substation_names
    n_sub = len(Sub)

    kVBase = []
    for key in Sub:
        kVBase = Bus[key]['kVBase']

    ## check if this bus is connected to any existing sub buses through a LINE
    for key_line in Line:
        if key_line != 'AllLineNames' and key_line != 'LineNumbers' and key_line != 'NumLines':
            if Line[key_line]['Bus1'] in Sub and Line[key_line]['Bus2'] in Sub:
                Sub.remove(Line[key_line]['Bus2']) ## normally a switch bus

    subData = {}
    if n_sub == 0: # no bus defined, then return an empty array
        pass
    else: # if any buses, read data
        subData['AllSubNames'] = Sub
        subData['SubNumbers'] = list(range(n_sub))
        subData['NumSubs'] = n_sub

        ## Search Coordinates...................
        sub_number = -1
        for sub_name in Sub:
            ### read data
            sub_number += 1
            subData[sub_name] = {} # 2-D dict
            subData[sub_name]['Number'] = sub_number
            subData[sub_name]['BusName'] = sub_name ## !!!! Should be a bus name
            subData[sub_name]['BusNum'] = Bus[sub_name]['Number']
            subData[sub_name]['NumPhases'] = Bus[sub_name]['NumPhases']
            subData[sub_name]['kVBase'] = Bus[sub_name]['kVBase']

            subData[sub_name]['GridForming'] = True
            subData[sub_name]['Status'] = True

            subData[sub_name]['Outage'] = False
            subData[sub_name]['Status'] = True

            subData[sub_name]['Phases'] = [1,1,1]
            subData[sub_name]['Device'] = 'substation'
            subData[sub_name]['Pmin'] = 0
            subData[sub_name]['Pmax'] = 1e8
            subData[sub_name]['Qmin'] = -1e7
            subData[sub_name]['Qmax'] = 1e7

            if Bus[sub_name]['Coorddefined'] == 1:
                subData[sub_name]['Coorddefined'] = Bus[sub_name]['Coorddefined']
                subData[sub_name]['Coord_X'] = Bus[sub_name]['Coord_X']
                subData[sub_name]['Coord_Y'] = Bus[sub_name]['Coord_Y']
            else: ## if sourcebus coordinates are not available, check the other bus of the line connected to it
                candi_name = sub_name

                max_try = 10
                try_time = 0
                found = 0
                while found == 0 and try_time < max_try:
                    try_time += 1
                    if found == 0:
                        for key in Line:
                            if key != 'AllLineNames' and key != 'LineNumbers' and key!= 'NumLines':
                                if Line[key]['Bus1'] == candi_name:
                                    bus2_name = Line[key]['Bus2']
                                    if Bus[ bus2_name ]['Coorddefined'] == 1:
                                        subData[sub_name]['Coorddefined'] = 1
                                        subData[sub_name]['Coord_X'] = Bus[bus2_name]['Coord_X']
                                        subData[sub_name]['Coord_Y'] = Bus[bus2_name]['Coord_Y']
                                        found = 1
                                        break
                                    else:
                                        candi_name = bus2_name
                                elif Line[key]['Bus2'] == candi_name:
                                    bus1_name = Line[key]['Bus1']
                                    if Bus[ bus1_name ]['Coorddefined'] == 1:
                                        subData[sub_name]['Coorddefined'] = 1
                                        subData[sub_name]['Coord_X'] = Bus[bus1_name]['Coord_X']
                                        subData[sub_name]['Coord_Y'] = Bus[bus1_name]['Coord_Y']
                                        found = 1
                                        break
                                    else:
                                        candi_name = bus1_name
                    ## if the loop is not broken yet::
                    if found == 0:
                        for key in Transformer:
                            if key != 'AllTransNames' and key != 'TransNumbers' and key!= 'NumTrans':
                                if Transformer[key]['Bus1'] == candi_name:
                                    bus2_name = Transformer[key]['Bus2']
                                    if Bus[ bus2_name ]['Coorddefined'] == 1:
                                        subData[sub_name]['Coorddefined'] = 1
                                        subData[sub_name]['Coord_X'] = Bus[bus2_name]['Coord_X']
                                        subData[sub_name]['Coord_Y'] = Bus[bus2_name]['Coord_Y']
                                        found = 1
                                        break
                                    else:
                                        candi_name = bus2_name
                                elif Transformer[key]['Bus2'] == candi_name:
                                    bus1_name = Transformer[key]['Bus1']
                                    if Bus[ bus1_name ]['Coorddefined'] == 1:
                                        subData[sub_name]['Coorddefined'] = 1
                                        subData[sub_name]['Coord_X'] = Bus[bus1_name]['Coord_X']
                                        subData[sub_name]['Coord_Y'] = Bus[bus1_name]['Coord_Y']
                                        found = 1
                                        break
                                    else:
                                        candi_name = bus1_name

    return subData

### import circuit data #####################################################
def func_import_circuit(dss):
    Circuit = {}
    Circuit['Name'] = dss.Circuit.Name()
    Circuit['DataPath'] = dss.Basic.DataPath()
    Circuit['VoltageBases'] = dss.Settings.VoltageBases()
    Circuit['NumBuses'] = dss.Circuit.NumBuses()
    Circuit['NumNodes'] = dss.Circuit.NumNodes()
    Circuit['NumLines'] = len(dss.Lines.AllNames())
    Circuit['NumLoads'] = len(dss.Loads.AllNames())
    Circuit['NumTrans'] = len(dss.Transformers.AllNames())
    Circuit['NumRegctrls']  = len(dss.RegControls.AllNames())
    Circuit['NumCaps']  = len(dss.Capacitors.AllNames())

    return Circuit

### get all the model data
def func_model_import(feeder_folder, dss_master_file, Substation_names):
    cwd = os.getcwd()
    parent_folder = os.path.abspath(os.path.join(cwd, os.pardir))
    project_folder = os.path.abspath(os.path.join(parent_folder, os.pardir))
    feeder_to_import = project_folder + '\\data\\opendss\\' + feeder_folder + '\\' + dss_master_file
    Substation_names = ['150']

    dss.run_command('clear')
    dss.run_command('compile ' + feeder_to_import)
    ### import data through dss COM interface with OpenDSS
    Bus         = func_import_bus(dss)
    Load        = func_import_load(dss)
    Line        = func_import_line(dss, Bus) # update Bus_Dict in case some diabled lines cause non-indexed buses
    Transformer = func_import_transformer(dss)
    Capacitor   = func_import_capacitor(dss)
    Regcontrol   = func_import_regulator(dss)
    Circuit     = func_import_circuit(dss)
    Substation  = func_get_substation(dss, Bus, Line, Transformer, Substation_names)

    print('Feeder model: ' + feeder_folder + ' is imported from ' + feeder_to_import)

    return Bus, Load, Line, Transformer, Capacitor, Regcontrol, Circuit, Substation

if __name__ == '__main__':
    ### select test feeder to import
    feeder_folder = '123Bus'
    dss_master_file = 'IEEE123Master_FLISR.dss'
    Substation_names = ['150']

    [Bus, Load, Line, Transformer, Capacitor, Regulator, Circuit, Substation] = \
        func_model_import(feeder_folder, dss_master_file, Substation_names)


    # ## Print & check node coordinates in case we need to draw the feeder diagram
    # for i in Bus['AllBusNames']:
    #     if Bus[i]['Coorddefined']:
    #         print(i + ' : ' + str(Bus[i]['Coord_X']) + ', ' + str(Bus[i]['Coord_Y']))
    #     else:
    #         print(i + ' : ' + str('XXXX' + ', ' + 'YYYY'))

    # ## pickle system data, saved in the opendss feeder model folderr
    # obj = [Bus, Load, Line, Transformer, Capacitor, Regulator, Substation, Circuit]
    # pickle.dump(obj, open(feeder_folder + "_system_data.dat", "wb"), True)

