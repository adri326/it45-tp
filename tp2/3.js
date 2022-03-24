"use strict";

const N_WAGONS = 10;
const MAX_VOLUME = 100;

const VOLUMES = [56, 23, 15, 34, 12, 37, 13, 46, 24];

function minimize_wagons() {
    for (let n_wagons = 1; n_wagons <= N_WAGONS; n_wagons++) {
        let current = new Array(VOLUMES.length).fill(0);
        let sum = new Array(n_wagons).fill(0);

        w: while (true) {
            // Count in base `n_wagons`
            current[0]++;
            for (let n = 0; n < current.length && current[n] >= n_wagons; n++) {
                if (n == current.length - 1) break w;
                current[n] = 0;
                current[n + 1]++;
            }

            // Assert the volume condition
            sum.fill(0);
            for (let n = 0; n < current.length; n++) {
                sum[current[n]] += VOLUMES[n];
            }
            let valid = true;
            for (let n = 0; n < sum.length; n++) {
                valid = valid && sum[n] <= MAX_VOLUME;
            }

            if (valid) return [n_wagons, sum, current];
        }
    }
}

function minimize_difference() {
    function evaluate(sum) {
        let count = 0;
        let average = 0;
        for (let n = 0; n < sum.length; n++) {
            if (sum[n] > 0) {
                count++;
                average += sum[n];
            }
        }
        average /= count;

        let res = 0;
        for (let n = 0; n < sum.length; n++) {
            if (sum[n] > 0) {
                res += (sum[n] - average) * (sum[n] - average);
            }
        }

        return res;
    }

    function list_actions(state, sum, n) {
        let res = [];
        let o = 0;
        for (; o < sum.length && state[o] >= 0; o++) {
            if (sum[o] + VOLUMES[n] <= 100) res.push(o);
        }
        if (o < sum.length) {
            res.push(o);
        }
        return res;
    }

    function rec(state, sum, n) {
        if (n == VOLUMES.length) return [state, sum, evaluate(sum)];

        let best = [null, null, Infinity];

        for (let action of list_actions(state, sum, n)) {
            let state2 = state.slice();
            state2.push(action);
            let sum2 = sum.slice();
            sum2[action] += VOLUMES[n];

            let [best_state, best_sum, score] = rec(state2, sum2, n + 1);


            if (score < best[2]) {
                best = [best_state, best_sum, score];
            }
        }

        return best;
    }

    let [best_state, best_sum, best_score] = rec([], new Array(N_WAGONS).fill(0), 0);

    return [best_state.reduce((acc, act) => Math.max(acc, act)) + 1, best_sum, best_state];
}

let p = 0;
for (let [n_wagons, sum, solution] of [minimize_wagons(), minimize_difference()]) {
    console.log(`== ${["Minimize wagons", "Minimize difference"][p]} ==`);
    p += 1;
    console.log(`Found solution for ${n_wagons} wagons:`);
    for (let n = 0; n < n_wagons; n++) {
        console.log(`Wagon ${n+1}: ${sum[n]}m³`);
        for (let c = 0; c < VOLUMES.length; c++) {
            if (solution[c] == n) {
                console.log(`| ${c} (${VOLUMES[c]}m³)`);
            }
        }
    }
    console.log();
}
