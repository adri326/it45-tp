# === Problem section ===

param n integer;
param m integer;

set SITES := 1..n;

param capacity{i in SITES};
param distance{i in SITES, j in SITES};
param max_distance;

var x{i in SITES} binary;

# We want to maximize the capacity of our owned places
maximize z: sum{i in SITES} x[i]*capacity[i];

# There can only be m owned places
s.t. x_m: sum{i in SITES} x[i] = m;

# Compute δᵢⱼ, which is equal to 1 if distance[i, j] > max_distance
param delta{i in SITES, j in SITES} := if distance[i, j] > max_distance then 1 else 0;

# Verify that the distance between two different owned places is always less than 30km away from other owned places
# NOTE: This *could* be written as:
# s.t. dist_max{i in SITES, j in SITES: i != j}: x[i] * x[j] * delta[i, j] == 0;
# But GMPL won't allow it because it isn't in linear form.
# NOTE: The i != j isn't necessary here, because `distance[i, i] <= 30`, which wouldn't trigger the condition.
# I left it here in case someone puts a big value in the distance matrix.
s.t. dist_max{i in SITES, j in SITES: i != j}: x[i] + x[j] <= 2 - delta[i, j];

solve; # === Result section ===

printf "Capacity = ";
printf sum{i in SITES} x[i] * capacity[i];
printf "\n";

# Print which places were chosen
display{i in SITES: x[i]} i;

data; # === Data section ===

param m := 3;
param n := 6;
param capacity :=
    1 134
    2 167
    3 189
    4 182
    5 136
    6 192;

param distance: 1 2 3 4 5 6 :=
    1  0 55 52 20 24 24
    2 55  0 43 40 43 36
    3 52 43  0 53 34 27
    4 20 40 53  0 51 25
    5 24 43 34 51  0 17
    6 24 36 27 25 17  0;
param max_distance := 30;

end;
