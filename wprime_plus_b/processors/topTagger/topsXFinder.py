import math


top_mass_pdg = 173.1
w_mass_pdg = 80.4
muon_mass_pdg = 0.10566
electron_mass_pdg = 0.511e-3


def topsXFinder(events , top_type, types_jets, results):
    boosted_top_mass_low = 100
    boosted_top_mass_up = 400 
    
    if (top_type['boosted_top']):
        boosted_top(events, types_jets['jetsTop_M'],
                    boosted_top_mass_low,
                    boosted_top_mass_up, 
                    results
                    )
                    
                    
    # allCombinationsResults are the results with all found top candidates, including those that may share the same final state objects. Duplicates are removed in a later step.
    all_combinations_results = {'top': [], 'top_AK4jets': [], 'w': [], 'top_topology_decay': [], 'chi2': []}
    chi2_range = 5
    
    
    # 1. partially-boosted top
    partially_boosted_top_sigma = 37.14
    partially_boosted_w_sigma = 20.09  # sigma comes from reco-gen level
    partially_boosted_top_mass_low = 100
    partially_boosted_top_mass_up = 300
    partially_boosted_w_mass_low = 40
    partially_boosted_w_mass_up = 200  # Top mass range and W mass range comes from reco-gen level
    
    
    if (top_type['partially_boosted_top']):
        partially_boosted_top(events, 
                              types_wjets['jetsW_M'],
                              types_jets['jetsB_L'], 
                              partially_boosted_top_sigma,
                              partially_boosted_w_sigma,
                              partially_boosted_w_mass_low,
                              partially_boosted_w_mass_up,
                              partially_boosted_top_mass_low,
                              partially_boosted_top_mass_up, chi2_range,
                              all_combinations_results
                              )
    
    # 2. resolved hadronic top
    hadronic_top_sigma = 35.02
    hadronic_w_sigma = 24.98
    hadronic_top_mass_low = 100
    hadronic_top_mass_up = 300
    hadronic_w_mass_low = 40
    hadronic_w_mass_up = 200  # Top mass range and W mass range comes from reco-gen level
    
    
    if (top_type['hadronic_top']): 
        hadronic_top(events, 
                    types_jets['jets_cleaned'], 
                    types_jets['jetsB_L'],
                    hadronic_top_sigma, 
                    hadronic_w_sigma, 
                    hadronic_w_mass_low,
                    hadronic_w_mass_up, 
                    hadronic_top_mass_low, 
                    hadronic_top_mass_up,
                    chi2_range, 
                    all_combinations_results
                    )
                    
                    
    # 3. resolved leptonic top
    leptonic_top_sigma = 120.8
    leptonic_w_sigma = 127.2
    leptonic_top_mass_low = 100
    leptonic_top_mass_up = 600
    leptonic_w_mass_low = 40
    leptonic_w_mass_up = 600  # Top mass range and W mass range comes from reco-gen level
    
    
     if (top_type['leptonic_top']): 
        single_leptonic_top(events, 
                            select_muons(events), 
                            select_electrons(events),
                            types_jets['jetsB_L'], 
                            leptonic_top_sigma, 
                            leptonic_w_sigma,
                            leptonic_w_mass_low, 
                            leptonic_w_mass_up, 
                            leptonic_top_mass_low,
                            leptonic_top_mass_up, 
                            chi2_range, 
                            all_combinations_results
                            )
 
     # now we have all the top candidates in the last 3 cases, and also their chi2
    # First we will order them in chi2 increasing order, and then if the later one have the same elements with the previous one, we will just remove it.
    remove_order(all_combinations_results, results)
    
    
    
    
    
    
# reconstruct the boosted top
def boosted_top(events: ak.Array, year: str,  top_mass_low, top_mass_up, results):

    top_jets = selected_top_jets(events, year)
    
    for key in top_jets:
        boosted_top_mass = top_jets[key].mass
        if top_mass_low < boosted_top_mass < top_mass_up:
            results['top'].append(top_jets[key])
            results['w'].append(ROOT.TLorentzVector(0, 0, 0, 0))               # Check if we can use TLorentzVector
            eesults['top_topology_decay'].append(0)
            results['chi2'].append(0.0000001 * index)
         
    

# reconstruct partially-boosted top   
def partially_boosted_top(events: ak.Array, top_sigma, w_sigma, w_mass_low, w_mass_up,
                          top_mass_low, top_mass_up, chi2_range, all_combinations_results):
                          
    fat_jets = events.FatJet
    jets = events.Jet
    bjets = selected_bjets(events, year)   
    wjets = selected_w_jets(events, year)

     # 1 boosted W matched with several b jets
    for key in wjets:
        boosted_w_b = []  # the b jets matched with boosted W
        chi2_partially = []  # chi2_partially keep all the value of chi2 with boosted W.
        boosted_w_mass = wjets.mass
        boosted_w_p4 = wjets              
        
        if w_mass_low < boosted_w_mass < w_mass_up:
        
            for key in bjets['jetsB_L']:      
            
                b_w_mass = (boosted_w_p4 + jets).mass
                
                if top_mass_low < b_w_mass < top_mass_up:
                
                    chi2a = chi2(b_w_mass, boosted_w_mass, top_sigma, w_sigma)
                    
                    if chi2a < chi2_range:
                    
                        boosted_w_b_idx.append(selected_bjets_idx[ijet])
                        chi2_partially.append(chi2a)
                        
            if len(chi2_partially) > 0:
            
                sorted_chi2_partially, sorted_partially_boosted_b = pick_1b(chi2_partially, boosted_w_b)    
                partially_boosted_top = [sorted_partially_boosted_b[0]]
                
                add_results(sorted_chi2_partially[0], 
                            partially_boosted_top_idx, boosted_w_p4,
                            jets[sorted_partially_boosted_b[0]].p4(), 
                            all_combinations_results, 1)


    
# calculate chi2
def chi2(top_mass, w_mass, top_sigma, w_sigma):
    t = (top_mass - top_mass_pdg) / top_sigma
    w = (w_mass - w_mass_pdg) / w_sigma
    chi2 = math.pow(t, 2) + math.pow(w, 2)
    
    return chi2
    
    
    
# select only 1 b jet for each W (this b should have the lowest chi2)            
def pick_1b(chi2_list, b_candidate_idx):
    # order them in the chi2-increasing order, so that we can get the lowest one as the b candidate
    sorted_chi2 = []
    sorted_b_candidate = []
    
    sorted_chi2_b_zip = sorted(zip(chi2_list, b_candidate_idx), key=lambda x: x[0])
    sorted_chi2, sorted_b_candidate = map(list, zip(*sorted_chi2_b_zip))
    return sorted_chi2, sorted_b_candida      
            
               
    

def add_results(chi2, top_AK4jets_list, w_p4, b_p4, all_combinations_results, topology):
    all_combinations_results['chi2'].append(chi2)
    all_combinations_results['top_AK4jets'].append(top_AK4jets_list)
    all_combinations_results['w'].append(w_p4)
    all_combinations_results['top'].append(w_p4 + b_p4)
    all_combinations_results['top_topology_decay'].append(topology)
    
    

# reconstruct hadronic top
def hadronic_top(events, selected_jets_idx, selected_bjets_idx, top_sigma, w_sigma, w_mass_low, w_mass_up, top_mass_low,
                 top_mass_up, chi2_range, all_combinations_results):
    
    jets = events.Jet
    
    hadronic_wA_list = []  # hadronically decaying W has 2 jets (A and B)
    hadronic_wB_list = []
    
    hadronic_w(events, selected_jets, w_mass_low, w_mass_up, hadronic_wA_list, hadronic_wB_list)
       
    for index in hadronic_wA_list:
        hadronic_w_b_idx = []
        chi2_hadronic = []

        if (len(hadronic_wA_list) > 0) and (len(hadronic_wB_list) > 0):
            hadronic_w_mass = (jets[hadronic_wA_list[index]] + jets[hadronic_wB_list[index]]).mass()
            hadronic_w_p4 = (jets[hadronic_wA_list[index]] + jets[hadronic_wB_list[index]])
            
            for ijet in xrange(0, len(selected_bjets_idx)):
            
                if ((selected_bjets[ijet] != hadronic_wA_list[index]) and (selected_bjets_idx[ijet] !=  hadronic_wB_list[index])):
                    b_w_mass = (jets[selected_bjets[ijet]] + hadronic_w_p4).mass()
                    
                    if (top_mass_low < b_w_mass < top_mass_up):
                        chi2_a = chi2(b_w_mass, hadronic_w_mass, top_sigma, w_sigma)
                        
                        if (chi2_a < chi2_range):
                            hadronic_w_b.append(selected_bjets[ijet])
                            chi2_hadronic.append(chi2_a)
                            
            if (len(chi2_hadronic) > 0):
            
                sorted_chi2_hadronic, sorted_hadronic_b = pick_1b(chi2_hadronic, hadronic_w_b
                hadronic_top = [hadronic_wA_list[index], hadronic_wB_list[index], sorted_hadronic_b[0]]
                
                add_results(sorted_chi2_hadronic[0], 
                            hadronic_top_idx, 
                            hadronic_w_p4, 
                            jets[sorted_hadronic_b[0]],
                            all_combinations_results, 
                            2)
                 
                 
# reconstruct resolved hadronic top
def hadronic_w(events, selected_jets, w_mass_low, w_mass_up, hadronic_wA_list, hadronic_wB_list):

    jets = events.Jet
    selected_jets = []  # avoid changing origin jets list
    
    for element in selected_jets:
        selected_jets.append(element)
        
    for wjet in xrange(0, len(selected_jets_idx1) / 2):
    
        min_w_mass = 9999
        w_jet1 = -1
        w_jet2 = -1
        
        for ijet in range(0, len(selected_jets)):
        
            for jjet in range(ijet + 1, len(selected_jets)):
            
                w_mass = (jets[selected_jets[ijet]] + jets[selected_jets[jjet]]).mass()
                
                if not abs(w_mass - w_mass_pdg) < min_w_mass: continue
                
                if not w_mass_low < w_mass < w_mass_up: continue
                
                min_w_mass = abs(w_mass - w_mass_pdg)
                w_jet1 = selected_jets[ijet]
                w_jet2 = selected_jets[jjet]
                
        if w_jet1 != -1 and w_jet2 != -1:
        
            hadronic_wA_list.append(w_jet1)
            hadronic_wB_list.append(w_jet2)
            selected_jets.remove(w_jet1)
            selected_jets.remove(w_jet2)




# reconstruct resolved leptonic top
def single_leptonic_top(events, selected_muon, selected_electron, selected_bjets, top_sigma, w_sigma,
                        w_mass_low, w_mass_up, top_mass_low, top_mass_up, chi2_range, all_combinations_results):
                       
    jets = events.Jet
    muons = events.Muon
    eles = events.Electron
        
    leptonic_w_mass = []
    leptonic_b = []
    chi2_lep = []
    
    
    if (len(selected_muon_idx) == 1):
    
        lepton_p4 = muons[selected_muon[0]]
        lepton_mass_pdg = muon_mass_pdg
        
        
    if (len(selected_electron_idx) == 1):
    
        lepton_p4 = eles[selected_electron[0]]
        lepton_mass_pdg = electron_mass_pdg
        
    if ((len(selected_muon) + len(selected_electron) == 1):
    
        metpz = met_pz(events, w_mass_pdg, lepton_mass_pdg, lepton_p4)
        metE = math.sqrt(events.MET.pt * event.MET.pt + metpz * metpz)
        
        neutrino = TLorentzVector(events.MET_pt * math.cos(events.MET.phi), events.MET.pt * math.sin(events.MET.phi), metpz, metE)                             
                                       
        single_leptonic_w(lepton_p4, w_mass_low, w_mass_up, neutrino, leptonic_w_mass)
        
        if (leptonic_w_mass[0] != 9999):
        
            for ijet in range(0, len(selected_bjets)):
            
                b_w_mass = (lepton_p4 + neutrino + jets[selected_bjets[ijet]]).mass()
                
                if top_mass_low < b_w_mass < top_mass_up:
                
                    chi2a = chi2(b_w_mass, leptonic_w_mass[0], top_sigma, w_sigma)
                    
                    if chi2a < chi2_range:
                    
                        leptonic_b_idx.append(selected_bjets[ijet])
                        chi2_lep.append(chi2a)
                        
            if len(chi2_lep) > 0:
            
                sorted_chi2_leptonic, sorted_leptonic_b = pick_1b(chi2_lep, leptonic_b)
                
                leptonic_top = [sorted_leptonic_b[0]]
                
                add_results(sorted_chi2_leptonic[0], 
                            leptonic_top, 
                            lepton_p4 + neutrino,
                            jets[sorted_leptonic_b[0]], 
                            all_combinations_results, 
                            3)


                    
                        

# using MET to calculate the metpz
def met_pz(event, w_mass_pdg, lepton_mass_pdg, lepton):

    met_px = event.MET.pt * math.cos(event.MET.phi)
    met_py = event.MET.pt * math.sin(event.MET.phi)
    
    met_pz = 0.0
    
    lepton_px = lepton.pt * math.cos(lepton.phi)
    lepton_py = lepton.pt * math.sin(lepton.phi)
    lepton_pz = lepton.pt()) * math.sinh(lepton.phi)
    lepton_E = math.sqrt(lepton.pt * lepton.pt + lepton_pz * lepton_pz + lepton.mass * lepton.mass)
    
    a = w_mass_pdg * w_mass_pdg - lepton_mass_pdg * lepton_mass_pdg + 2.0 * (lepton_px * met_px + lepton_py * met_py)
    
    A = 4.0 * (lepton_E * lepton_E - lepton_pz * lepton_pz)
    B = -4.0 * a * lepton_pz
    C = 4.0 * lepton_E * lepton_E * (met_px * met_px + met_py * met_py) - a * a
    
    
    tmp_root = B * B - 4.0 * A * C
    tmp_sol1 = 0.0
    tmp_sol2 = 0.0
    
    if (tmp_root < 0):
    
        met_pz = -B / (2 * A)  # take real part of complex roots
        
    else:
    
        tmp_sol1 = (-B + math.sqrt(tmp_root)) / (2.0 * A)
        tmp_sol2 = (-B - math.sqrt(tmp_root)) / (2.0 * A)
        
        if ((abs(tmp_sol2 - lepton)) < (abs(tmp_sol1 - lepton))):
        
            met_pz = tmp_sol2
            
        else:
        
            met_pz = tmp_sol1
            
    return met_pz


# reconstruct the resolved leptonic W (muon decay)
def single_leptonic_w(lepton, w_mass_low, w_mass_up, met, leptonic_w_mass):
    
    l_nu_mass = (lepton + met).mass()
    
    if (w_mass_low < l_nu_mass < w_mass_up):
    
        w_mass = l_nu_mass
        
    else:
    
        w_mass = 9999
    
    leptonic_w_mass.append(w_mass)



# remove the same element and order in chi2 increasing order
def remove_order(all_combinations_results, results):
    
    # Now we have all the top candidates, and also their chi2. First, we will order them with chi2 increasing order.
    
    if len(all_combinations_results['top_AK4jets']) > 0:
    
        pre_result = list(zip(all_combinations_results['chi2'],
                            all_combinations_results['top_AK4jets'],    
                            all_combinations_results['top'], all_combinations_results['w'],
                            all_combinations_results['top_topology_decay']))
                              
                              
        sorted_pre_result = sorted(pre_result, key=lambda x: x[0])
        
        (sorted_chi2, sorted_top_AK4_jets, sorted_top, sorted_w, sorted_top_topology_decay) = list(map(list, zip(*sorted_pre_result)))
            
        # Now we have these sorted lists. They are in the right order. We need to remove the one who have the same elements with previous.
        list_index = get_non_overlapping_lists(sorted_top_AK4_jets)  # final_top_ak4_jets: keep those jets after the removal, list_index: keep the position of those left elements after removal

    else:
    
        sorted_chi2 = []
        sorted_top_AK4_jets = []
        sorted_top = []
        sorted_w = []
        sorted_top_topology_decay = []
        list_index = []
    
    for index in list_index:
    
        results['chi2'].append(sorted_chi2[index])
        results['top'].append(sorted_top[index])
        results['w'].append(sorted_w[index])
        results['top_topology_decay'].append(sorted_top_topology_decay[index])


# Check out the same element in 2 list and remove the same element in the later list. Also keep the position of these removed elements in order to remove corresponding element of chi2, top/w candidate and top_topology_decay
def get_non_overlapping_lists(listA):
    
    listB = []
    positions = []
    
    for i, sublist in enumerate(listA):
    
        if i == 0:
        
            listB.append(sublist)
            positions.append(i)
            
        else:
        
            overlap = False
            
            for element in sublist:
            
                if element in listB[-1]:
                
                    overlap = True
                    
                    break
                    
            if not overlap:
            
                listB.append(sublist)
                positions.append(i)
                
    return positions
