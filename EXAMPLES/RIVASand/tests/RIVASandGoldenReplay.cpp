/* Replay the five frozen V8 material-point histories against the native
 * OpenSees kernel copy. This test intentionally has no OpenSees dependency. */
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using Row = std::unordered_map<std::string, double>;

static std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    return values;
}

static bool readRows(const std::string &path, std::vector<Row> &rows)
{
    std::ifstream input(path);
    std::string line;
    if (!input || !std::getline(input, line)) return false;
    const std::vector<std::string> header = split(line);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split(line);
        if (fields.size() != header.size()) return false;
        Row row;
        for (std::size_t i = 0; i < header.size(); ++i) {
            if (header[i] != "case_id") row[header[i]] = std::stod(fields[i]);
        }
        rows.push_back(row);
    }
    return !rows.empty();
}

static bool close(double actual, double expected)
{
    /* The oracle was generated with NumPy/Python math. A 1e-6 relative
     * tolerance accommodates libm operation ordering while remaining much
     * tighter than any constitutive calibration tolerance. */
    return std::abs(actual-expected) <= 1.0e-10+1.0e-6*std::abs(expected);
}

static bool compare(const std::string &caseName, int step,
                    const char *field, double actual, double expected)
{
    if (close(actual, expected)) return true;
    std::cerr << caseName << " step " << step << " " << field
              << ": actual=" << std::setprecision(17) << actual
              << " expected=" << expected
              << " error=" << std::abs(actual-expected) << '\n';
    return false;
}

static bool replay(const std::string &directory, const std::string &file)
{
    std::vector<Row> rows;
    if (!readRows(directory + "/" + file, rows)) {
        std::cerr << "cannot read " << directory << "/" << file << '\n';
        return false;
    }

    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters);
    const Row &first = rows.front();
    riva_tensor_t stress = {
        first.at("stress_xx"), first.at("stress_yy"),
        first.at("stress_zz"), first.at("stress_xy"),
        first.at("stress_yz"), first.at("stress_xz")};
    riva_state_t state = {};
    if (!riva_initialize_material(&parameters, &material, stress,
                                 first.at("void_ratio"), &state))
        return false;
    riva_tensor_t reference = stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_begin_dynamic_phase(&parameters, reference, &state)) return false;

    const std::string caseName = file.substr(0, file.size()-4);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const Row &row = rows[i];
        riva_update_info_t info = {};
        if (i > 0) {
            const riva_tensor_t increment = {
                row.at("deps_xx"), row.at("deps_yy"), row.at("deps_zz"),
                row.at("deps_xy"), row.at("deps_yz"), row.at("deps_xz")};
            const int substeps = static_cast<int>(row.at("nsub"));
            if (!riva_update_material(&parameters, &material, increment,
                                     substeps, &state, &stress, &info)) {
                std::cerr << caseName << " failed at step " << i << '\n';
                return false;
            }
            if (info.accepted_substeps !=
                    static_cast<int>(row.at("accepted_substeps")) ||
                info.reversal_registered !=
                    static_cast<int>(row.at("reversal_registered"))) {
                std::cerr << caseName << " event mismatch at step " << i
                          << '\n';
                return false;
            }
        } else {
            stress = state.stress;
        }

        bool ok = true;
        ok &= compare(caseName, static_cast<int>(i), "stress_xx",
                      stress.xx, row.at("stress_xx"));
        ok &= compare(caseName, static_cast<int>(i), "stress_yy",
                      stress.yy, row.at("stress_yy"));
        ok &= compare(caseName, static_cast<int>(i), "stress_zz",
                      stress.zz, row.at("stress_zz"));
        ok &= compare(caseName, static_cast<int>(i), "stress_xy",
                      stress.xy, row.at("stress_xy"));
        ok &= compare(caseName, static_cast<int>(i), "stress_yz",
                      stress.yz, row.at("stress_yz"));
        ok &= compare(caseName, static_cast<int>(i), "stress_xz",
                      stress.xz, row.at("stress_xz"));
        ok &= compare(caseName, static_cast<int>(i), "mean_effective_pressure",
                      riva_pressure(stress), row.at("mean_effective_pressure"));
        ok &= compare(caseName, static_cast<int>(i), "q",
                      riva_q(stress), row.at("q"));
        ok &= compare(caseName, static_cast<int>(i), "compatibility_residual",
                      riva_compatibility_residual(&state),
                      row.at("compatibility_residual"));
        ok &= compare(caseName, static_cast<int>(i), "void_ratio",
                      state.void_ratio, row.at("void_ratio"));
        ok &= compare(caseName, static_cast<int>(i), "D_ir",
                      state.D_ir, row.at("D_ir"));
        ok &= compare(caseName, static_cast<int>(i), "D_re",
                      state.D_re, row.at("D_re"));
        ok &= compare(caseName, static_cast<int>(i), "cyclic_amplitude",
                      state.cyclic_amplitude, row.at("cyclic_amplitude"));
        ok &= compare(caseName, static_cast<int>(i), "bias_ratchet_strain",
                      state.bias_ratchet_strain,
                      row.at("bias_ratchet_strain"));
        if (state.reversals != static_cast<int64_t>(row.at("reversals")) ||
            state.amplitude_reversals !=
                static_cast<int64_t>(row.at("amplitude_reversals"))) {
            std::cerr << caseName << " counter mismatch at step " << i
                      << '\n';
            ok = false;
        }
        if (!ok) return false;
    }
    std::cout << "PASS: " << caseName << " (" << rows.size()-1
              << " increments)\n";
    return true;
}

int main(int argc, char **argv)
{
    std::string directory;
    if (argc > 1) {
        directory = argv[1];
    } else {
        const std::string source = __FILE__;
        const std::size_t separator = source.find_last_of("/\\");
        directory = source.substr(0, separator) + "/../golden_data";
    }
    const std::vector<std::string> files = {
        "cyclic_zero_bias_reference.csv",
        "cyclic_bias025_reference.csv",
        "cyclic_bias0375_dense.csv",
        "cyclic_bias025_four_substeps.csv",
        "mixed_3d_bias015.csv"};
    bool ok = true;
    for (const std::string &file : files) ok &= replay(directory, file);
    if (!ok) return 1;
    std::cout << "PASS: all frozen RIVA-Sand golden histories\n";
    return 0;
}
