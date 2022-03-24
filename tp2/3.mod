param n_wagons := 10;
param n_crates := 9;
param max_vol := 100;

param vol{crate in 1..n_crates};

# Index of the wagon in which it'll end up
var x{crate in 1..n_crates, wagon in 1..n_wagons}, binary;

# An auxiliary variable to keep track of which wagon contains a crate
var populated{wagon in 1..n_wagons}, binary;
s.t. c_populated{wagon in 1..n_wagons, crate in 1..n_crates}: populated[wagon] >= x[crate, wagon]; # populated cannot be false if any crate is in it
s.t. c_populated_2{wagon in 1..n_wagons}: sum{crate in 1..n_crates} x[crate, wagon] >= populated[wagon]; # populated cannot be true if no crate is in it
s.t. c_populated_3{wagon in 2..n_wagons}: populated[wagon] <= populated[wagon - 1]; # we make sure that the first wagons are populated first

minimize z: sum{wagon in 1..n_wagons} populated[wagon];

s.t. c_vol{wagon in 1..n_wagons}: sum{crate in 1..n_crates} x[crate, wagon] * vol[crate] <= max_vol;
s.t. c_loaded{crate in 1..n_crates}: sum{wagon in 1..n_wagons} x[crate, wagon] = 1;

solve;

display z;

display{wagon in 1..n_wagons, crate in 1..n_crates: x[crate, wagon] == 1} (crate, wagon);

display{wagon in 1..n_wagons} (wagon, sum{crate in 1..n_crates} x[crate, wagon] * vol[crate]);

data;

param vol :=
    1 56
    2 23
    3 15
    4 34
    5 12
    6 37
    7 13
    8 46
    9 24;

end;
