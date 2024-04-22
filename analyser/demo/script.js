// 'HB': '3' ,# hydrogen bond - orange
// 'PS': '1', # pi stacking - red
// 'PC': '9', # oi cation - pink
// 'SB': '11',# saltbridge - purple
// 'HA': '4', # halogen bond - yellow
// 'WB': '10',# water bridge - cyan
// 'MC': '14' # metal complex - ochra


document.addEventListener("DOMContentLoaded", function () {
var motifs_color = NGL.ColormakerRegistry.addSelectionScheme(
[
["purple", "265"],
["red", ""],
["pink", ""],
["orange", "61 or 69 or 88 or 85 or 249 or 252 or 261"],
["yellow", ""],
["cyan", ""],
["ochra", ""],
],
"Motifs Color Scheme"
);

// Create NGL Stage object
var stage = new NGL.Stage("viewport", { backgroundColor: "#ffffff" });

// Handle window resizing
window.addEventListener("resize", function (event) {
stage.handleResize();
}, false);

// Code for example: showcase/rhodopsin
stage.loadFile("https://files.rcsb.org/download/7XP6.pdb").then(function (o) {
//receptor
o.addRepresentation("cartoon", {
sele: "protein and :R",
color: "#c1c1c1",
material: "flat",
});

//ligand
o.addRepresentation("licorice", {
sele: "SY9",
color: "#4193d5",
radius: 0.4,
});

//motif
o.addRepresentation("licorice", {
sele: " (61 or 69 or 88 or 85 or 249 or 252 or 261 or 265) and :R ",
color: motifs_color,
});

o.addRepresentation("label", {
sele: "((61 or 69 or 88 or 85 or 249 or 252 or 261 or 265) and :R) and .CA",
color: motifs_color,
labelType: "res",
xOffset: 1,
yOffset: 1,
zOffset: 1,
radius: 4,
});

stage.autoView();
});

// Handle window scrolling
window.addEventListener("scroll", function (event) {
stage.mouseControls.remove("scroll");
}, false);
});