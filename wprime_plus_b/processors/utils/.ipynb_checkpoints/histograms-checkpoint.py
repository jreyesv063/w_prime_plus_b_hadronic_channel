import hist
import numpy as np

# --------------------------
# histogram axes definition
# --------------------------
# systematics axis
syst_axis = hist.axis.StrCategory([], name="variation", growth=True)

# dataset axis
dataset_axis = hist.axis.StrCategory([], name="dataset", growth=True)

# jet axes
jet_pt_axis = hist.axis.Variable(
    edges=[20, 60, 90, 120, 150, 180, 210, 240, 300, 500],
    name="jet_pt",
)
jet_eta_axis = hist.axis.Regular(
    bins=50,
    start=-2.4,
    stop=2.4,
    name="jet_eta",
)
jet_phi_axis = hist.axis.Regular(
    bins=50,
    start=-np.pi,
    stop=np.pi,
    name="jet_phi",
)
# met axes
ttbar_met_axis = hist.axis.Variable(
    edges=[50, 75, 100, 125, 150, 175, 200, 300, 500],
    name="met",
)
qcd_met_axis = hist.axis.Variable(
    edges=[0, 50, 75, 100, 125, 150, 175, 200, 300, 500],
    name="met",
)
met_phi_axis = hist.axis.Regular(
    bins=50,
    start=-np.pi,
    stop=np.pi,
    name="met_phi",
)
# lepton axes
lepton_pt_axis = hist.axis.Variable(
    edges=[30, 60, 90, 120, 150, 180, 210, 240, 300, 500],
    name="lepton_pt",
)
lepton_eta_axis = hist.axis.Regular(
    bins=50,
    start=-2.4,
    stop=2.4,
    name="lepton_eta",
)
lepton_phi_axis = hist.axis.Regular(
    bins=50,
    start=-np.pi,
    stop=np.pi,
    name="lepton_phi",
)
lepton_reliso = hist.axis.Regular(
    bins=25,
    start=0,
    stop=1,
    name="lepton_reliso",
)
# lepton + bjet axes
lepton_bjet_dr_axis = hist.axis.Regular(
    bins=30,
    start=0,
    stop=5,
    name="lepton_bjet_dr",
)
lepton_bjet_mass_axis = hist.axis.Variable(
    edges=[40, 75, 100, 125, 150, 175, 200, 300, 500],
    name="lepton_bjet_mass",
)
# lepton + missing energy axes
lepton_met_mass_axis = hist.axis.Variable(
    edges=[40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
    name="lepton_met_mass",
)
lepton_met_delta_phi_axis = hist.axis.Regular(
    bins=30, start=0, stop=4, name="lepton_met_delta_phi"
)
# lepton + missing energy + bjet axes
lepton_met_bjet_mass_axis = hist.axis.Variable(
    edges=[40, 75, 100, 125, 150, 175, 200, 300, 500, 800],
    name="lepton_met_bjet_mass",
)


# --------------------------
# ttbar analysis histograms
# --------------------------
# jet histogram
ttbar_jet_hist = hist.Hist(
    dataset_axis,
    jet_pt_axis,
    jet_eta_axis,
    jet_phi_axis,
    syst_axis,
    hist.storage.Weight(),
)
# met histogram
ttbar_met_hist = hist.Hist(
    dataset_axis,
    ttbar_met_axis,           # met
    met_phi_axis,             # met_phi
    syst_axis,                # variation
    hist.storage.Weight(),
)
# lepton histogram
ttbar_lepton_hist = hist.Hist(
    dataset_axis,
    lepton_pt_axis,
    lepton_eta_axis,
    lepton_phi_axis,
    syst_axis,
    hist.storage.Weight(),
)
# lepton + bjet histogram
ttbar_lepton_bjet_hist = hist.Hist(
    dataset_axis,
    lepton_bjet_dr_axis,
    lepton_bjet_mass_axis,
    syst_axis,
    hist.storage.Weight(),
)
# lepton + missing energy histogram
ttbar_lepton_met_hist = hist.Hist(
    dataset_axis,
    lepton_met_mass_axis,
    lepton_met_delta_phi_axis,
    syst_axis,
    hist.storage.Weight(),
)
# lepton + missing energy + bjet histogram
ttbar_lepton_met_bjet_hist = hist.Hist(
    dataset_axis,
    lepton_met_bjet_mass_axis,
    syst_axis,
    hist.storage.Weight(),
)


# -------------------------------
# ztoll control region histogram
# -------------------------------
# dilepton mass axis
dilepton_mass_axis = hist.axis.Regular(
    bins=40, start=50, stop=200, name="dilepton_mass"
)
# dilepton mass histogram
dilepton_mass_hist = hist.Hist(dataset_axis, dilepton_mass_axis, hist.storage.Weight())


# ------------------------
# qcd analysis histograms
# ------------------------
# lepton ID region axis
lepton_id_axis = hist.axis.StrCategory([], name="lepton_id", growth=True)

# met histogram
qcd_met_hist = hist.Hist(
    dataset_axis,
    qcd_met_axis,
    lepton_id_axis,
    hist.storage.Weight(),
)
# lepton + MET mass histogram
qcd_lepton_met_hist = hist.Hist(
    dataset_axis,
    lepton_met_mass_axis,
    lepton_id_axis,
    hist.storage.Weight(),
)


# -------------------------------
# tt CR: Hadronic channel
# -------------------------------
# Tau axes
tau_pt_axis = hist.axis.Variable(
    edges=[20,40,60,80,100,120,140,160,180,200,240,260,280,300,320,340,360,380,400,420,440,460,480,500],
    name="tau_pt",
)
tau_eta_axis = hist.axis.Regular(
    bins=24,   # each 0.2
    start=-2.4,
    stop=2.4,
    name="tau_eta",
)
tau_phi_axis = hist.axis.Regular(
    bins=32,    # each 0.2
    start=-np.pi,
    stop=np.pi,
    name="tau_phi",
)

# MET axes
met_pt_axis = hist.axis.Variable(
    edges=[200, 240, 260, 280, 300, 320, 340, 360 ,380,400,420,440,460,480,500],
    name="met_pt",
)
met_phi_axis = hist.axis.Regular(
    bins=32,    # each 0.2
    start=-np.pi,
    stop=np.pi,
    name="met_phi",
)

# B axes
b_pt_axis = hist.axis.Variable(
    edges=[10,20,40,60,80,100,120,140,160,180,200,240,260,280,300,320,340,360,380,400,420,440,460,480,500],
    name="b_pt",
)
b_eta_axis = hist.axis.Regular(
    bins=24,   # each 0.2
    start=-2.4,
    stop=2.4,
    name="b_eta",
)
b_phi_axis = hist.axis.Regular(
    bins=32,    # each 0.2
    start=-np.pi,
    stop=np.pi,
    name="b_phi",
)

# top axes
top_mrec_axis = hist.axis.Variable(
    edges=[120 ,130, 140, 150, 160, 170, 180, 190, 200, 210, 220],
    name="top_mrec",
)


# tau histogram
ttbar_h_tau_hist = hist.Hist(
    dataset_axis,              # "dataset"
    tau_pt_axis,            # "tau_pt"
    tau_eta_axis,           # "tau_eta"
    tau_phi_axis,           # "tau_phi"
    syst_axis,                 # "variation"
    hist.storage.Weight(),
)

ttbar_h_met_hist = hist.Hist(
    dataset_axis,              # "dataset"
    met_pt_axis,            # "met_pt"
    met_phi_axis,           # "met_phi"
    syst_axis,                 # "variation"
    hist.storage.Weight(),
)

ttbar_h_b_hist = hist.Hist(
    dataset_axis,              # "dataset"
    b_pt_axis,            # "b_pt"
    b_eta_axis,           # "b_eta"
    b_phi_axis,           # "b_phi"
    syst_axis,                 # "variation"
    hist.storage.Weight(),
)


ttbar_h_top_hist = hist.Hist(
    dataset_axis,              # "dataset"
    top_mrec_axis,            # "top_mrec"
    syst_axis,                 # "variation"
    hist.storage.Weight(),
)

