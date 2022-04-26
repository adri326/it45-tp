# The matrix of the arcs in the graph G
param n;
param coords{i in 1..52, d in 0..1};
# param c{i in 1..n, j in 1..n};
var x{i in 1..n, j in 1..n} binary;

param c{i in 1..n, j in 1..n} := if i != j then sqrt(
    (coords[i, 0] - coords[j, 0]) * (coords[i, 0] - coords[j, 0])
    + (coords[i, 1] - coords[j, 1]) * (coords[i, 1] - coords[j, 1])
) else 10000000;

minimize distance: sum{i in 1..n, j in 1..n} c[i, j] * x[i, j];

s.t. contrainte_1{i in 1..n}: sum{j in 1..n} x[i, j] = 1;
s.t. contrainte_2{j in 1..n}: sum{i in 1..n} x[i, j] = 1;

# MTZ secundary variable to guarantee that only one circuit is made
var u{i in 1..n} integer >= 0;

s.t. init: u[1] = 1;
s.t. range{i in 2..n}: 2 <= u[i] <= n;
s.t. conn{i in 1..n, j in 2..n: i != j}: u[i] - u[j] + n*x[i, j] <= n-1;

solve;

display{i in 1..n, j in 1..n: x[i, j] = 1} (i-1, j-1);
display sum{i in 1..n, j in 1..n} x[i, j]*c[i, j];

data;

param coords: 0 1 :=
    1 565.0 575.0
    2 25.0 185.0
    3 345.0 750.0
    4 945.0 685.0
    5 845.0 655.0
    6 880.0 660.0
    7 25.0 230.0
    8 525.0 1000.0
    9 580.0 1175.0
    10 650.0 1130.0
    11 1605.0 620.0
    12 1220.0 580.0
    13 1465.0 200.0
    14 1530.0 5.0
    15 845.0 680.0
    16 725.0 370.0
    17 145.0 665.0
    18 415.0 635.0
    19 510.0 875.0
    20 560.0 365.0
    21 300.0 465.0
    22 520.0 585.0
    23 480.0 415.0
    24 835.0 625.0
    25 975.0 580.0
    26 1215.0 245.0
    27 1320.0 315.0
    28 1250.0 400.0
    29 660.0 180.0
    30 410.0 250.0
    31 420.0 555.0
    32 575.0 665.0
    33 1150.0 1160.0
    34 700.0 580.0
    35 685.0 595.0
    36 685.0 610.0
    37 770.0 610.0
    38 795.0 645.0
    39 720.0 635.0
    40 760.0 650.0
    41 475.0 960.0
    42 95.0 260.0
    43 875.0 920.0
    44 700.0 500.0
    45 555.0 815.0
    46 830.0 485.0
    47 1170.0 65.0
    48 830.0 610.0
    49 605.0 625.0
    50 595.0 360.0
    51 1340.0 725.0
    52 1740.0 245.0;

param n := 10;
