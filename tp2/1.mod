param n integer >= 0; # Nombre de produits
param m integer >= 0; # Nombre de composants

param c{j in 1..n}; # Gain par produit
param b{i in 1..m}; # Nombre de composants
param a{i in 1..m, j in 1..n}; # Nombre de composants i nécessaires pour fabriquer j

var x{j in 1..n} >= 0;

maximize gain: sum{j in 1..n} c[j] * x[j];

s.t. stock{i in 1..m}: sum{j in 1..n} a[i, j] * x[j] <= b[i];

solve;

display gain;
display{j in 1..n} x[j];
display{i in 1..m} (i, sum{j in 1..n} a[i, j] * x[j]);

data;

param n := 2;
param m := 3;
param c :=
    1 3
    2 4;
param b :=
    1 180
    2 120
    3 150;

# Yes, it's backwards and I hate it
param a: 1 2 :=
    1  2 3
    2  2 1
    3  1 3;

end;
