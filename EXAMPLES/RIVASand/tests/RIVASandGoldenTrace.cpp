/* Golden-oracle (re)generation helper: step the CURRENT kernel through the
 * strain schedule of an existing golden CSV and print, per step, the event
 * info + compatibility residual + the full state vector (schema order).
 * Used by regen_goldens.py to assemble V9 reference CSVs. */
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

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
        for (std::size_t i = 0; i < header.size(); ++i)
            if (header[i] != "case_id") row[header[i]] = std::stod(fields[i]);
        rows.push_back(row);
    }
    return !rows.empty();
}

int main(int argc, char **argv)
{
    if (argc < 2) { std::cerr << "usage: GoldenTrace <case.csv>\n"; return 1; }
    std::vector<Row> rows;
    if (!readRows(argv[1], rows)) { std::cerr << "cannot read\n"; return 1; }

    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters);
    const Row &first = rows.front();
    riva_tensor_t stress = {
        first.at("stress_xx"), first.at("stress_yy"), first.at("stress_zz"),
        first.at("stress_xy"), first.at("stress_yz"), first.at("stress_xz")};
    riva_state_t state = {};
    if (!riva_initialize_material(&parameters, &material, stress,
                                 first.at("void_ratio"), &state)) return 1;
    riva_tensor_t reference = stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_begin_dynamic_phase(&parameters, reference, &state)) return 1;

    std::cout << std::setprecision(17);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        riva_update_info_t info = {};
        if (i > 0) {
            const Row &row = rows[i];
            const riva_tensor_t increment = {
                row.at("deps_xx"), row.at("deps_yy"), row.at("deps_zz"),
                row.at("deps_xy"), row.at("deps_yz"), row.at("deps_xz")};
            const int substeps = static_cast<int>(row.at("nsub"));
            if (!riva_update_material(&parameters, &material, increment,
                                     substeps, &state, &stress, &info)) {
                std::cerr << "update failed at step " << i << '\n';
                return 1;
            }
        } else {
            stress = state.stress;
        }
        double values[RIVA_STATE_VALUE_COUNT];
        riva_state_values(&state, values);
        std::cout << info.accepted_substeps << ','
                  << info.reversal_registered << ','
                  << riva_compatibility_residual(&state);
        for (int v = 0; v < RIVA_STATE_VALUE_COUNT; ++v)
            std::cout << ',' << values[v];
        std::cout << '\n';
    }
    return 0;
}
