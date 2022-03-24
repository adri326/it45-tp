# == Parameter declaration ==
set ELEMS;
set ALLOYS;

# The min- and max- proportion of the different elements in our final alloy
param min_content{elem in ELEMS} >= 0;
param max_content{elem in ELEMS} >= 0;

# How much of an element each input alloy contains
param content{alloy in ALLOYS, elem in ELEMS} >= 0;
# How much of that input alloy we have (this value scales with the `quantite` parameter)
param stocks{alloy in ALLOYS} >= 0;
# How much each input alloy costs
param couts{alloy in ALLOYS} >= 0;

# How much alloy we're making
param quantite >= 0;

# == Problem definition ==

# The proportion of each input alloy in the final alloy
var x{alloy in ALLOYS} >= 0;

minimize cout: sum{alloy in ALLOYS} couts[alloy] * x[alloy];

# We don't want bubbles in our alloy, so we require that the sum of the ratios of the input alloys adds up to 1
s.t. bubbles: sum{alloy in ALLOYS} x[alloy] = 1.0;

# Enforce the minimum and maximum proportion for the different elements in the final alloy
s.t. min_c{elem in ELEMS}: sum{alloy in ALLOYS} content[alloy, elem] * x[alloy] >= min_content[elem];
s.t. max_c{elem in ELEMS}: sum{alloy in ALLOYS} content[alloy, elem] * x[alloy] <= max_content[elem];

# Enforce the maximum quantity of each input alloy
s.t. stock{alloy in ALLOYS}: x[alloy] * quantite <= stocks[alloy];

solve; # == Result section ==

printf "Cout total: ";
printf cout * quantite;
printf "\n";
display{alloy in ALLOYS} (alloy, x[alloy] * quantite);
display{elem in ELEMS}: (elem, sum{alloy in ALLOYS} content[alloy, elem] * x[alloy]);

data; # == Parameter definition ==

# It looks like this is how one writes sets now; every piece of documentation I found had it differently, but eh
set ELEMS := "C" "Cu" "Mn";
set ALLOYS := "Fer 1" "Fer 2" "Fer 3" "Cuivre 1" "Cuivre 2" "Aluminium 1" "Aluminium 2";

param min_content :=
    "C" 0.02
    "Cu" 0.004
    "Mn" 0.012;
param max_content :=
    "C" 0.03
    "Cu" 0.006
    "Mn" 0.0165;

# `param VAR: i1 i2 ... in :=`
# is short for `param VAR := j1 i1 x11   j1 i2 x12 ...`
param content: "C" "Cu" "Mn" :=
    "Fer 1"        0.025 0     0.013
    "Fer 2"        0.03  0     0.008
    "Fer 3"        0     0.003 0
    "Cuivre 1"     0     0.9   0
    "Cuivre 2"     0     0.96  0.04
    "Aluminium 1"  0     0.004 0.012
    "Aluminium 2"  0     0.006 0.0;

# `param: VAR1 VAR2 := ...` is short for `param VAR1 := j1 x11  j2 x12 ...` and `param VAR2 := j1 x21  j2 x22 ...`
param: stocks couts :=
    "Fer 1"        4000 1.2
    "Fer 2"        3000 1.5
    "Fer 3"        6000 0.9
    "Cuivre 1"     5000 1.3
    "Cuivre 2"     2000 1.2
    "Aluminium 1"  3000 1.2
    "Aluminium 2"  2500 1;

param quantite := 5000;
