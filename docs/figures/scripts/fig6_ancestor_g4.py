#!/usr/bin/env python3
# fig6_ancestor_g4.py -- Figure 6: the reconstructed pre-duplication Saccharomycotina PIF1 ancestor folded
# on the parallel G4 (AlphaFold3, top-ranked model_0), with the ancestral Arg-wedge (position 89) engaging
# the 5' G-tetrad. Renders with PyMOL open-source (arm64; ~/miniforge3/envs/pymol).
# Run from repo root:  ~/miniforge3/envs/pymol/bin/pymol -cq docs/figures/scripts/fig6_ancestor_g4.py
import math
from pymol import cmd, util

CIF = "results/asr/af3/folds/g4/fold_ancestor_dupnode_g4_model_0.cif"
OUT = "docs/figures/fig6_ancestor_g4.png"
coral, teal, purple = "0xD55E00", "0x009E73", "0x9C5A8C"

cmd.reinitialize()
cmd.load(CIF, "anc")
cmd.bg_color("white")
cmd.hide("everything")

cmd.select("prot", "polymer.protein")
cmd.select("dna", "polymer.nucleic")
cmd.select("kion", "resn K")
cmd.select("wedge", "prot and resi 89")
# engaged tetrad = the 4 guanines nearest the wedge (robust to chain/numbering)
_wa = cmd.get_model("wedge").atom
_dg = sorted(set(a.resi for a in cmd.get_model("dna and resn DG").atom),
             key=lambda ri: min(math.dist(a.coord, b.coord) for a in _wa
                                 for b in cmd.get_model("dna and resn DG and resi %s" % ri).atom))
cmd.select("tetrad", "dna and resn DG and resi " + "+".join(_dg[:4]))

# protein: faint cartoon for context (thin, mostly transparent)
cmd.show("cartoon", "prot")
cmd.set("cartoon_transparency", 0.7, "prot")
cmd.set("cartoon_side_chain_helper", 1)
cmd.color("grey80", "prot")

# G4: thin backbone tube (no filled rings -> no blob); coral
cmd.set("cartoon_ring_mode", 0)
cmd.set("cartoon_tube_radius", 0.3)
cmd.show("cartoon", "dna")
cmd.color(coral, "dna")

# the engaged 5' tetrad: guanines as sticks, orange, so the tetrad face reads clearly
cmd.show("sticks", "tetrad")
cmd.set("stick_radius", 0.2, "tetrad")
cmd.color("orange", "tetrad")
util.cnc("tetrad")

# K+ channel ions
cmd.show("spheres", "kion")
cmd.set("sphere_scale", 0.4, "kion")
cmd.color(purple, "kion")

# the ancestral wedge
cmd.show("sticks", "wedge")
cmd.set("stick_radius", 0.28, "wedge")
cmd.color(teal, "wedge")
util.cnc("wedge")

# nearest wedge-atom -> tetrad-atom contact (the workflow/22 metric), dashed + labelled
_wa2 = cmd.get_model("wedge").atom
_ta = cmd.get_model("tetrad").atom
best = min(((math.dist(a.coord, b.coord), a, b) for a in _wa2 for b in _ta), key=lambda t: t[0])
cmd.distance("wcontact", "wedge and name %s" % best[1].name,
             "tetrad and resi %s and name %s" % (best[2].resi, best[2].name))
cmd.hide("labels", "wcontact")
cmd.set("dash_color", "grey20"); cmd.set("dash_width", 3.5); cmd.set("dash_gap", 0.35)
print("[fig6] nearest wedge->5'tetrad = %.2f A" % best[0])

# small label on the wedge
cmd.set("label_size", 15)
cmd.set("label_color", "black")
cmd.set("label_font_id", 7)
cmd.label("wedge and name CZ", '"R89"')

# view + render: frame the wedge/tetrad interface
cmd.orient("wedge or tetrad")
cmd.zoom("wedge or tetrad", 7)
cmd.turn("y", 12)
cmd.set("ray_shadows", 0)
cmd.set("ray_trace_mode", 1)
cmd.set("antialias", 2)
cmd.set("ray_trace_color", "grey20")
cmd.ray(2200, 1650)
cmd.png(OUT, dpi=300)
print("[fig6] wrote", OUT)
