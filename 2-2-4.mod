set INGREDIENT := 1..6;

param cout{i in INGREDIENT};
param mg{i in INGREDIENT};
param sucre{i in INGREDIENT};
param oeuf{i in INGREDIENT};
param eau{i in INGREDIENT};

var x{i in INGREDIENT} >= 0;

minimize z: sum{i in INGREDIENT} x[i] * cout[i];

s.t. mg_sum: sum{i in INGREDIENT} x[i]*mg[i] = 21.5;
s.t. sucre_sum: sum{i in INGREDIENT} x[i]*sucre[i] = 21.0;
s.t. oeuf_sum: sum{i in INGREDIENT} x[i]*oeuf[i] = 1.2;
s.t. eau_sum: sum{i in INGREDIENT} x[i]*eau[i] = 56.3;
s.t. conservation_matiere: sum{i in INGREDIENT} x[i] = 100.0;

solve;

display sum{i in INGREDIENT} x[i] * cout[i];
display{i in INGREDIENT} x[i];

data;

# 1 = Crème
# 2 = Jaune d'oeuf frais
# 3 = Lait entier en poudre
# 4 = Jaune d'oeuf surgelé, sucré
# 5 = Sirop de sucre de canne
# 6 = Eau

param cout :=
    1 0.45
    2 0.6
    3 0.15
    4 0.3
    5 0.12
    6 0;

# For the last exercise:
# param cout :=
#     1 0.6
#     2 1
#     3 0.15
#     4 0.3
#     5 0.12
#     6 0;

param mg :=
    1 0.4
    2 0.5
    3 0.12
    4 0.3
    5 0.0
    6 0.0;
param sucre :=
    1 0.0
    2 0.0
    3 0.0
    4 0.14
    5 0.7
    6 0.0;
param oeuf :=
    1 0.0
    2 0.4
    3 0.0
    4 0.14
    5 0.0
    6 0.0;
param eau :=
    1 0.6
    2 0.1
    3 0.88
    4 0.16
    5 0.3
    6 1.0;
