param capacity;
set objects := 1..6;

param value{i in objects};
param weight{i in objects};

var x{i in objects} binary;
# var x{i in objects} integer >= 0;

maximize z: sum{i in objects} x[i] * value[i];
s.t. capacity_accepted: sum{i in objects} x[i] * weight[i] <= capacity;

solve;

display{i in objects} x[i];
display z;

data;

param capacity := 27

param value :=
    1 20
    2 16
    3 11
    4 12
    5 9
    6 1;

param weight :=
    1 9
    2 8
    3 6
    4 5
    5 4
    6 1;
