from pathlib import Path
import json
from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

@dataclass
class AlgorithmData:
     atoms: list[int] = field(default_factory=list)
     search: list[float] = field(default_factory = list)
     build: list[float] = field(default_factory = list)
     name: str = field(default_factory = str)

input_dir = Path(__file__).resolve().parent.parent / "data" / "benchmark_statistics"

#passes string
def extract_algorithm (algorithm):

    result = AlgorithmData()

    result.name = algorithm

    for file in input_dir.glob("*.json"):
    
            with open(file) as f:
                data = json.load(f)
    
                entry = data["algorithms"].get(algorithm)

                #brute force may not exist for some files
                if entry is None: 
                    continue

                result.atoms.append(data["atoms"])
                result.search.append(entry["search"])

                build_time = entry.get("build")

                #brute force has no build time
                if build_time is None:
                    continue

                result.build.append(entry["build"])

    return result

def fit(atoms, times, degree = 1):

    coefficients = np.polyfit(atoms, times, degree)

    return coefficients

def predict(coefficients, atoms):

    prediction = np.polyval(coefficients, atoms)

    return prediction

def plot(predicted, atoms, search):
     plt.scatter(atoms, search)
     plt.plot(atoms, predicted)
     plt.show()

def compute_r2_score(predicted, search):
     r2 = r2_score(search, predicted)
     return r2

def compute_rmse(predicted, search):
     rmse = np.sqrt(np.mean((np.array(search) - predicted) **2))
     return rmse

def get_times(algorithm, metric):
    if metric == "search":
        times = algorithm.search

    elif metric == "build":
        times = algorithm.build

    else:
        raise ValueError("Unknown metric")

    #brute force has no build
    if len(times) == 0:
        return None

    return times

def main():

    #type AlgorithmData
    brute_force = extract_algorithm("brute_force")
    kd_tree = extract_algorithm("kd_tree")
    cell_list = extract_algorithm("cell_list")
    OMP_cell_list = extract_algorithm("OMP_cell_list")

    algorithms = [brute_force, kd_tree, cell_list, OMP_cell_list]

    metrics = ["search", "build"]

    best_models = []

    for metric in metrics:
        print (f"Fitting the best model for the {metric} function")

        for algorithm in algorithms:

            atoms = algorithm.atoms
            time = get_times(algorithm, metric)

            if time == None:
                continue
            nlogn = np.array(atoms) * np.log(atoms)

            models = [("Linear", atoms, fit(atoms, time, 1)),
            ("Quadratic", atoms, fit(atoms, time, 2)),
            ("N log N", nlogn, fit(nlogn, time, 1))]

            best_r2 = -1
            print("_________________________")
            print(f"Results for {algorithm.name}_{metric}: ")
            print("_________________________")

            for name, x, coeffs in models:

                pred = predict(coeffs, x)

                r2 = compute_r2_score(pred, time)
                rmse = compute_rmse(pred, time)

                penalties = {
                    "Linear": 0.0,
                    "Quadratic": 0.002,
                    "N log N": 0.001
                }

                if name == "Linear":

                    a, b = coeffs

                    print(f"Linear: f(x) = {a}x, + {b}")
                    print(f"r2_score: {r2}")
                    print(f"rmse: {rmse}")
                    print()

                # to penalize extraneous complexity
                if name == "Quadratic":

                    a, b, c = coeffs

                    print(f"Quadratic: f(x) = {a}x^2, + {b}x + {c}")
                    print(f"r2_score: {r2}")
                    print(f"rmse: {rmse}")
                    print()

                if name == "N log N":

                    a, b = coeffs
            
                    print(f"N log N: f(x) = {a} xlog x, + {b}")
                    print(f"r2_score: {r2}")
                    print(f"rmse: {rmse}")
                    print()

                penalty = penalties[name]

                if r2 - penalty > best_r2:
                    best_prediction = pred
                    best_coeffs = coeffs
                    best_name = name
                    best_r2 = r2
                    rmse = compute_rmse(pred, time)

            best_model = {
                        "algorithm": algorithm.name,
                        "metric": metric,
                        "name": best_name,
                        "coeffs": best_coeffs,
                        "prediction": best_prediction,
                        "r2": best_r2,
                        "rmse": rmse,
                        "atoms": atoms,
                        "time": time
                        }
            
            best_models.append(best_model)

        for best_model in best_models:
            order = np.argsort(best_model["atoms"])

            x = np.array(best_model["atoms"])[order]
            y = np.array(best_model["time"])[order]
            prediction = np.array(best_model["prediction"][order])

            plt.scatter(x, y)
            plt.plot(x, prediction)

            plt.title = (f"{best_model["algorithm"]}: {best_model["metric"]}")
            plt.xlabel("Atoms")
            plt.ylabel("Milliseconds")

            filename = f"{best_model["algorithm"]}_{best_model["metric"]}.png"
            plt.savefig(filename, dpi = 300, bbox_inches="tight")
            plt.close()

            print(f"{best_model["algorithm"], best_model["metric"], best_model["name"]}")

if __name__ == "__main__":
    main()
    