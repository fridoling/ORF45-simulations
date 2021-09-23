import simplesbml
import re
import numpy as np

model = simplesbml.loadSBMLFile("../models/model_invitro_incell.xml")

def format_reaction(x):
    ss = x.split(" ")
    ss = [s for s in ss if '_factor' not in s]
    ss = [s for s in ss if '_on' not in s]
    if ss[0] == "*":
        ss = ss[1:]
    x_out = " ".join([format_reaction_component(s) for s in ss])
    x_out = x_out.replace("\\cdot \\cdot", ") \\cdot")
    return(x_out)

def format_reaction_component(s):
    mapping = {"*": "\cdot", "_": "\_"}
    if s not in "*-+":
        s = re.sub("(\(*)\s*(\w+)\s*(\)*)", r"\1 \\texttt{\2} \3", s)
    for k, v in mapping.items():
        s = s.replace(k, v)
    if '_a' in s:
        s = re.sub(r"(.*)\\_a", r"\\texttt{a} \\cdot \1", s)
    if '_d' in s:
        s = re.sub(r"(.*)\\_d", r"\\texttt{d} \\cdot \1", s)
    return(s)

def format_rate_law(x, remove_factors=True):
    ss = x.split(" ")
    if remove_factors:
        ss = [s for s in ss if '_factor' not in s]
        ss = [s for s in ss if '_on' not in s]
    if ss[0] == "*":
        ss = ss[1:]
    x_out = " ".join([format_rate_law_component(s) for s in ss])
    x_out = x_out.replace("\\cdot \\cdot", ") \\cdot")
    return(x_out)

def format_rate_law_component(s):
    mapping = {"*": "\cdot", "_": "\_"}
    if s not in "*-+":
        if s[0] in '(k':
            s = re.sub("(\(*)\s*(\w+)\s*(\)*)", r"\1 \\texttt{\2} \3", s)
        else:
            s = re.sub("(\(*)\s*(\w+)\s*(\)*)", r"\1 \\texttt{[\2]} \3", s)
    for k, v in mapping.items():
        s = s.replace(k, v)
    if '_a' in s:
        s = re.sub(r"(.*)\\_a", r"\\texttt{a} \\cdot \1", s)
    if '_d' in s:
        s = re.sub(r"(.*)\\_d", r"\\texttt{d} \\cdot \1", s)
    return(s)

groups = {
    "ER": re.compile("koff_ER"),
    "EpR": re.compile("koff_EpR"),
    "EO": re.compile("koff_EO"),
    "pEO": re.compile("koff_pEO"),
    "OR": re.compile("koff_OR"),
    "OpR": re.compile("koff_OpR"),
    "EP": re.compile("koff_EP|kon_EP"),
    "EK": re.compile("koff_EK|kon_EK"),
    "pE": re.compile("kp_E|kdp_E"),
    "RP2": re.compile("koff_RP2|kon_RP2"),
    "pR": re.compile("kp_R|kdp_R"),
    "pK": re.compile("kp_K|kdp_K")
}
f = {}
for group in groups.keys():
    f[group] = open('../../doc/reactions_'+group+'.tex', 'w')
    f[group].write("\\begin{tabular}{"+'ll'+"}\n")
    f[group].write("\\textbf{Reaction} & \\textbf{Rate law} \\\\\n")
    f[group].write("\\midrule\n")

r_ids = model.getListOfReactionIds()
counter = 0
for r_id in r_ids:
    reactants = []
    products = []
    for i_r in range(model.getNumReactants(r_id)):
        reactants.append(format_reaction(model.getReactant(r_id, i_r)))
    for i_p in range(model.getNumProducts(r_id)):
        products.append(format_reaction(model.getProduct(r_id, i_p)))
    rate_law = model.getRateLaw(r_id)
    if ('kp' in rate_law or 'kdp' in rate_law) and not 'kp_K' in rate_law:
        operator = " \\rightarrow"
    else:
        operator = " \\leftrightharpoons "
    group = [g for g, pat in groups.items() if pat.search(rate_law)]
    if len(group)>0:
        f[group[0]].write(r"$"+" + ".join(reactants)+operator+" + ".join(products)+" $ & $"+format_rate_law(rate_law)+" $ \\\\\n")
    else:
        counter+=1
for group in groups:
    f[group].write("\\end{tabular}")
    f[group].close()
