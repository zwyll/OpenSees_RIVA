/* Reproducible native runtime baseline for the production RIVA-Sand kernel.
 *
 * This deliberately benchmarks the solver-independent constitutive kernel,
 * not OpenSees assembly, the DSS stress servo, file I/O, or Python.  The
 * research loose-biased successor has no native implementation and therefore
 * must not be represented by this production-kernel timing.
 */
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

constexpr int kPointsPerCycle = 64;
constexpr int kCyclesPerHistory = 12;
constexpr int kUpdatesPerHistory = kPointsPerCycle * kCyclesPerHistory;
constexpr int kMeasuredRepeats = 7;
/* Keeps the complete four-resolution benchmark below roughly 20 seconds on
 * an Apple-silicon laptop while retaining hundreds of thousands of local
 * updates in every median. */
constexpr int kHistoriesPerRepeat = 250;
constexpr double kPi = 3.1415926535897932384626433832795;

using Clock = std::chrono::steady_clock;

riva_state_t initial_state(const riva_parameters_t &parameters,
                           const riva_material_parameters_t &material)
{
    /* Representative shallow biased DSS state in kPa. */
    const riva_tensor_t reference = {-19.4, -19.4, -40.0, 0, 0, 0};
    const riva_tensor_t biased = {-19.4, -19.4, -40.0, 6.0, 0, 0};
    riva_state_t state = {};
    if (!riva_initialize_material(
            &parameters, &material, biased, 0.657, &state) ||
        !riva_begin_dynamic_phase(&parameters, reference, &state)) {
        std::cerr << "native benchmark initialization failed\n";
        std::exit(2);
    }
    return state;
}

std::array<riva_tensor_t, kUpdatesPerHistory> strain_history()
{
    std::array<riva_tensor_t, kUpdatesPerHistory> history = {};
    double previous_gamma = 0.0;
    for (int step = 0; step < kUpdatesPerHistory; ++step) {
        const double cycle = static_cast<double>(step + 1) / kPointsPerCycle;
        const double gamma = 5.0e-4 * std::sin(2.0 * kPi * cycle);
        const double increment = gamma - previous_gamma;
        previous_gamma = gamma;
        /* Tensor shear is half engineering shear. */
        history[step] = {0, 0, 0, 0.5 * increment, 0, 0};
    }
    return history;
}

double run_repeat(const riva_parameters_t &parameters,
                  const riva_material_parameters_t &material,
                  const std::array<riva_tensor_t, kUpdatesPerHistory> &history,
                  int substeps, volatile double &checksum)
{
    const auto started = Clock::now();
    for (int repeat = 0; repeat < kHistoriesPerRepeat; ++repeat) {
        riva_state_t state = initial_state(parameters, material);
        for (const riva_tensor_t increment : history) {
            riva_tensor_t stress = {};
            riva_update_info_t info = {};
            if (!riva_update_material(
                    &parameters, &material, increment, substeps,
                    &state, &stress, &info) ||
                info.accepted_substeps != substeps ||
                !riva_finite_tensor(stress)) {
                std::cerr << "native benchmark update failed\n";
                std::exit(3);
            }
        }
        checksum += state.stress.xy + 1.0e-3 * state.lambda_total;
    }
    return std::chrono::duration<double>(Clock::now() - started).count();
}

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

} // namespace

int main()
{
    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters);
    const auto history = strain_history();
    volatile double checksum = 0.0;

    std::cout << "kernel,substeps,host_updates,constitutive_substeps,"
                 "median_seconds,ns_per_host_update,ns_per_constitutive_substep\n";
    for (const int substeps : {1, 2, 4, 8}) {
        /* Warm-up is outside every measured sample. */
        (void)run_repeat(parameters, material, history, substeps, checksum);
        std::vector<double> samples;
        samples.reserve(kMeasuredRepeats);
        for (int sample = 0; sample < kMeasuredRepeats; ++sample)
            samples.push_back(run_repeat(
                parameters, material, history, substeps, checksum));
        const double seconds = median(samples);
        const double host_updates =
            static_cast<double>(kUpdatesPerHistory) * kHistoriesPerRepeat;
        std::cout << "production," << substeps << ','
                  << static_cast<long long>(host_updates) << ','
                  << static_cast<long long>(host_updates * substeps) << ','
                  << std::setprecision(9) << seconds << ','
                  << 1.0e9 * seconds / host_updates << ','
                  << 1.0e9 * seconds / (host_updates * substeps) << '\n';
    }
    std::cerr << "checksum=" << std::setprecision(17) << checksum << '\n';
    return 0;
}
