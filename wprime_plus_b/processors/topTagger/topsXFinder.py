import math


top_mass_pdg = 173.1
w_mass_pdg = 80.4
muon_mass_pdg = 0.10566
electron_mass_pdg = 0.511e-3



##################################################
########## complementary functions ###############
##################################################

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
    
    return sorted_chi2, sorted_b_candidate


def add_results(chi2, top_AK4jets_list, w_p4, b_p4, all_combinations_results, topology):
    all_combinations_results['chi2'].append(chi2)
    all_combinations_results['top_AK4jets'].append(top_AK4jets_list)
    all_combinations_results['w'].append(w_p4)
    all_combinations_results['top'].append(w_p4 + b_p4)
    all_combinations_results['top_topology_decay'].append(topology)
    







##################################################
##############  Main function ####################
##################################################

def topsXFinder(events, selected_idx, top_type, results):
    # 0. pick out the boosted top
    boosted_top_mass_low = 100
    boosted_top_mass_up = 400  # Top mass range comes from reco-gen level
    
    
    if (top_type['boosted_top']):         
        boosted_top(events, 
                    selected_idx['jetsTop_M'], 
                    boosted_top_mass_low, 
                    boosted_top_mass_up, 
                    results)
        

    ################# 1. partially-boosted top #######################
    
    partially_boosted_top_sigma = 37.14
    partially_boosted_w_sigma = 20.09  # sigma comes from reco-gen level
    partially_boosted_top_mass_low = 100
    partially_boosted_top_mass_up = 300
    partially_boosted_w_mass_low = 40
    partially_boosted_w_mass_up = 200  # Top mass range and W mass range comes from reco-gen level
    
    
    chi2_range = 5
    
    # allCombinationsResults are the results with all found top candidates, including those that may share the same final state objects. Duplicates are removed in a later step.
    all_combinations_results = {'top': [], 
                                'top_AK4jets': [], 
                                'w': [], 
                                'top_topology_decay': [], 
                                'chi2': []}


    if (top_type['partially_boosted_top']):      
        partially_boosted_top(events, 
                              selected_idx['jetsW_M'],
                              selected_idx['jetsB_L'], 
                              partially_boosted_top_sigma,
                              partially_boosted_w_sigma,
                              partially_boosted_w_mass_low,
                              partially_boosted_w_mass_up,
                              partially_boosted_top_mass_low,
                              partially_boosted_top_mass_up, chi2_range,
                              all_combinations_results)

        

    ################## 2. resolved hadronic top ######################
    
    hadronic_top_sigma = 35.02
    hadronic_w_sigma = 24.98
    hadronic_top_mass_low = 100
    hadronic_top_mass_up = 300
    hadronic_w_mass_low = 40
    hadronic_w_mass_up = 200  # Top mass range and W mass range comes from reco-gen level
    
    if (top_type['hadronic_top']): 
        hadronic_top(events, 
                     selected_idx['jets_cleaned'], 
                     selected_idx['jetsB_L'],
                     hadronic_top_sigma, 
                     hadronic_w_sigma, 
                     hadronic_w_mass_low,
                     hadronic_w_mass_up, 
                     hadronic_top_mass_low, 
                     hadronic_top_mass_up,
                     chi2_range, 
                     all_combinations_results)

    
    ################## 3. resolved leptonic top ######################
    
    leptonic_top_sigma = 120.8
    leptonic_w_sigma = 127.2
    leptonic_top_mass_low = 100
    leptonic_top_mass_up = 600
    leptonic_w_mass_low = 40
    leptonic_w_mass_up = 600  # Top mass range and W mass range comes from reco-gen level
    
    if (top_type['leptonic_top']): 
        single_leptonic_top(events, 
                            selected_idx['muon'], 
                            selected_idx['electron'],
                            selected_idx['jetsB_L'], 
                            leptonic_top_sigma, 
                            leptonic_w_sigma,
                            leptonic_w_mass_low, 
                            leptonic_w_mass_up, 
                            leptonic_top_mass_low,
                            leptonic_top_mass_up, 
                            chi2_range, 
                            all_combinations_results)

    # now we have all the top candidates in the last 3 cases, and also their chi2
    # First we will order them in chi2 increasing order, and then if the later one have the same elements with the previous one, we will just remove it.
    
    remove_order(all_combinations_results, results)



##################################################
################  Boost top ######################
##################################################

# reconstruct the boosted top
def boosted_top(events, 
                selected_top_jets_idx, 
                top_mass_low, 
                top_mass_up, 
                results):
                    
    good_top_mass = (
        (selected_top_jets_idx.mass > top_mass_low)
        & (selected_top_jets_idx.mass < top_mass_up)
    )
    
    
    boosted_top = selected_top_jets_idx[good_top_mass]   
    
    results['top'].append(boosted_top)
    results['w'].append(TLorentzVector(0,0,0,0))    # append(ROOT.TLorentzVector(0, 0, 0, 0))  # https://coffeateam.github.io/coffea/api/coffea.nanoevents.methods.vector.LorentzVector.html
    results['top_topology_decay'].append(0)
    results['chi2'].append(0.0000001)




# reconstruct partially-boosted top
def partially_boosted_top(events, 
                          selected_w_jets_idx, 
                          selected_bjets_idx, 
                          top_sigma, 
                          w_sigma, 
                          w_mass_low, 
                          w_mass_up,
                          top_mass_low, 
                          top_mass_up, 
                          chi2_range, 
                          all_combinations_results):
    
    boosted_w_b_idx = []  # the b jets matched with boosted W
    chi2_partially = []  # chi2_partially keep all the value of chi2 with boosted W.


    ###### 1 boosted W matched with several b jets #####

    
    good_w_mass = (
        (selected_w_jets_idx.mass > w_mass_low)
        & (selected_w_jets_idx.mass < w_mass_up)
    )
    
    partially_boosted_top_w = selected_w_jets_idx[good_w_mass] 

    
    b_w_idx = selected_w_jets_idx + selected_bjets_idx
    
    good_b_w_mass = (
        (b_w_idx.mass > top_mass_low)
        & (b_w_idx.mass < top_mass_up)
    )
    
    partially_boosted_top_b_w = b_w[good_b_w_mass] 
    
    chi2a = chi2(partially_boosted_top_b_w.mass, partially_boosted_top_w.mass, top_sigma, w_sigma)
    
    if chi2a < chi2_range:
        boosted_w_b_idx.append(selected_bjets_idx)
        chi2_partially.append(chi2a)
        
    if len(chi2_partially) > 0:
        
        sorted_chi2_partially, sorted_partially_boosted_b = pick_1b(chi2_partially, boosted_w_b_idx)
        partially_boosted_top_idx = [sorted_partially_boosted_b[0]]
            
        add_results(sorted_chi2_partially[0], 
                    partially_boosted_top_idx, 
                    boosted_w_p4,
                    jets[sorted_partially_boosted_b[0]].p4(), 
                    all_combinations_results, 
                    1)


##################################################
#################### top types ###################
##################################################


# using MET to calculate the metpz
def met_pz(event, w_mass_pdg, lepton_mass_pdg, lepton_p4):
    
    met_px = event.MET_pt * math.cos(event.MET_phi)
    met_py = event.MET_pt * math.sin(event.MET_phi)
    met_pz = 0.0    
    
    
    lepton_px = lepton_p4.pt * math.cos(lepton_p4.phi)  # px = pT * cos(phi)
    lepton_py = lepton_p4.pt * math.sin(lepton_p4.phi)  # py = pT * sin(phi)
    lepton_pz = lepton_p4.pt * math.sinh(lepton_p4.phi)  # ¿Por qué  pz = pT * senh(phi)
    
    lepton_E = math.sqrt(np.power(lepton_p4.pt, 2) + np.power(lepton_pz, 2) + np.power(lepton_p4.mass, 2))


    a = w_mass_pdg**2 - lepton_mass_pdg**2 + 2.0 * (lepton_px * met_px + lepton_py * met_py)  
    
    A = 4.0 * (np.power(lepton_E,2) - np.power(lepton_pz,2))
    B = -4.0 * a * lepton_pz
    C = 4.0 * np.power(lepton_E,2) * (np.power(met_px,2) + np.power(met_py,2)) - a**2
    
    
    tmp_root = np.power(B,2) - 4.0 * A * C
    tmp_sol1 = 0.0
    tmp_sol2 = 0.0
    
    
    if (tmp_root < 0):
        met_pz = -B / (2 * A)  # take real part of complex roots
        
    else:
        tmp_sol1 = (-B + math.sqrt(tmp_root)) / (2.0 * A)
        tmp_sol2 = (-B - math.sqrt(tmp_root)) / (2.0 * A)
        
        if ((np.abs(tmp_sol2 - lepton_pz)) < (np.abs(tmp_sol1 - lepton_pz))):
            met_pz = tmp_sol2
         
        else:
            met_pz = tmp_sol1
            
    return met_pz


# reconstruct the resolved leptonic W (muon decay)
def single_leptonic_w(lepton_p4, w_mass_low, w_mass_up, met, leptonic_w_mass):
    
    l_nu_mass = (lepton_p4 + met).mass()
    
    if (w_mass_low < l_nu_mass < w_mass_up):
        
        w_mass = l_nu_mass
        
    else:
        
        w_mass = 9999
        
    leptonic_w_mass.append(w_mass)


# reconstruct resolved hadronic top ****************
def hadronic_w(events, selected_jets_idx, w_mass_low, w_mass_up, hadronic_wA_list, hadronic_wB_list): 

# reconstruct resolved leptonic top
def single_leptonic_top(events, 
                        selected_muon_idx,         # selected_idx['muon']
                        selected_electron_idx,     # selected_idx['electron']
                        selected_bjets_idx,        # selected_idx['jetsB_L']
                        top_sigma,                 # leptonic_top_sigma
                        w_sigma,                   # leptonic_w_sigma
                        w_mass_low, w_mass_up,     # leptonic_w_mass_low(top)
                        top_mass_low, top_mass_up, # leptonic_top_mass_low(top)
                        chi2_range,                # chi2_range
                        all_combinations_results   # all_combinations_results -> Lista de resultados
                       ): 


# reconstruct hadronic top *****************
def hadronic_top(events, 
                 selected_jets_idx, 
                 selected_bjets_idx, 
                 top_sigma, 
                 w_sigma, 
                 w_mass_low, 
                 w_mass_up, 
                 top_mass_low,
                 top_mass_up, 
                 chi2_range, 
                 all_combinations_results):




