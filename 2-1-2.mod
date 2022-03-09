var x1 >= 0, <= 1;
var x2 >= 0, <= 1;

param proteins{i in 1..3};
param lipids{i in 1..3};
param cost{i in 1..3};

minimize z: cost[1]*x1 + cost[2]*x2 + cost[3]*(1-x1-x2);

s.t. x3_exists: x1 + x2 <= 1;
s.t. sum_proteins: proteins[1] * x1 + proteins[2] * x2 + proteins[3] * (1 - x1 - x2) >= 0.22;
# s.t. sum_proteins: (proteins[1] - proteins[3]) * x1 + (proteins[2] - proteins[3]) * x2 >= 0.22 - proteins[3];
s.t. sum_lipids: lipids[1] * x1 + lipids[2] * x2 + lipids[3] * (1 - x1 - x2) >= 0.036;

solve;
data;

param proteins :=
    1 0.12
    2 0.52
    3 0.42;
param lipids :=
    1 0.02
    2 0.02
    3 0.1;
param cost :=
    1 25
    2 41
    3 39;
