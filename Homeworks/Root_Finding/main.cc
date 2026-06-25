// Nikolaos Antonios Iakynthos Nemet
// This project was created using a combination of Claude Code and GitHub Copilot.

// "main.cc" implementation file.
// Orchestrates all three parts of the Root Finding homework.

#include "hydrogen.h"
#include "newton.h"
#include "reporting.h"
#include "test_problems.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace pp {

namespace {

// Prints a one-line summary of a Newton solve result including the solution
// point, residual norm, status, and all diagnostic counters.
void print_newton_summary(const char* label, const NewtonResult& result) {
	std::printf("  %-28s x=%s  ||f||=%.3e  status=%s  iter=%d f_calls=%d J_calls=%d ls_steps=%d J_alloc=%d\n",
			label,
			result.x.to_string().c_str(),
			result.fx.norm(),
			newton_status_cstr(result.status),
			result.stats.iterations,
			result.stats.f_evaluations,
			result.stats.jacobian_evaluations,
			result.stats.line_search_iterations,
			result.stats.jacobian_allocations);
}

// Part A: demonstrate the Newton solver on scalar and 2D debug problems,
// then locate the minima of the Rosenbrock and Himmelblau functions by
// finding zeros of their analytic gradients.
void run_part_a() {
	std::puts("=== Part A (6 points): Newton with numerical Jacobian and line-search ===");
	std::puts("Debug systems + extrema from analytic gradients\n");

	NewtonOptions options;
	options.acc = 1e-10;
	options.alpha_min = 1e-3;
	options.max_iter = 100;
	options.line_search = LineSearchType::backtracking;

	// 1D sanity check.
	const NewtonResult scalar_res = newton_solve(simple_scalar_equation(), vector{5.0}, options);
	print_newton_summary("1D: x^2-2=0", scalar_res);
	std::printf("    expected root ~ %.15g\n", std::sqrt(2.0));

	// 2D sanity check.
	const NewtonResult twod_res = newton_solve(simple_two_dimensional_system(), vector{0.3, 2.8}, options);
	print_newton_summary("2D: linear+circle system", twod_res);
	std::printf("    expected root near {1,2}\n");

	// Rosenbrock extremum via roots of analytic gradient.
	const NewtonResult rosen_res = newton_solve(rosenbrock_gradient(), vector{1.2, 1.2}, options);
	print_newton_summary("Rosenbrock grad=0", rosen_res);
	std::printf("    Rosenbrock value at root = %.12g\n", rosenbrock_value(rosen_res.x));

	// Himmelblau minima from multiple starts to expose different basins.
	const VectorFunction himmel_grad = himmelblau_gradient();
	const std::vector<vector> starts{
			vector{3.0, 3.0},
			vector{-3.0, 3.0},
			vector{-3.0, -3.0},
			vector{3.0, -2.0},
	};

	std::puts("  Himmelblau minima from multiple initial guesses:");
	for (std::size_t i = 0; i < starts.size(); i++) {
		const NewtonResult r = newton_solve(himmel_grad, starts[i], options);
		std::printf("    start%-2zu -> x=%s  H(x)=%.3e  ||grad||=%.3e  status=%s\n",
				i + 1,
				r.x.to_string().c_str(),
				himmelblau_value(r.x),
				r.fx.norm(),
				newton_status_cstr(r.status));
	}

	std::puts("");
}

// Runs the shooting solver for the given config and bundles the result into
// a ConvergenceRow tagged with the varied parameter value for table output.
ConvergenceRow solve_for_config(const ShootingConfig& cfg, double parameter_value) {
	const ShootingResult root = find_ground_state_energy(cfg);
	const double exact = hydrogen_exact_energy_ground();

	ConvergenceRow row;
	row.parameter = parameter_value;
	row.energy = root.energy;
	row.abs_error = root.converged ? std::abs(root.energy - exact) : std::numeric_limits<double>::infinity();
	row.converged = root.converged;
	return row;
}

// Prints a parameter-convergence table showing the energy error as a
// function of the varied parameter (rmin, rmax, acc, or eps).
void print_convergence_table(const char* title, const std::vector<ConvergenceRow>& rows) {
	std::puts(title);
	std::puts("  parameter        E0                 |E0+1/2|         converged");
	for (const ConvergenceRow& row : rows) {
		std::printf("  %-12.5g  %-17.12g  %-13.3e   %s\n",
				row.parameter,
				row.energy,
				row.abs_error,
				row.converged ? "yes" : "no");
	}
	std::puts("");
}

// Part B: find the hydrogen ground-state energy via the shooting method,
// plot the wavefunction against the exact solution, and study convergence
// as rmin, rmax, and ODE tolerances are varied.
void run_part_b() {
	std::puts("=== Part B (3 points): Hydrogen bound state via shooting method ===");

	ShootingConfig base;
	base.rmin = 1e-3;
	base.rmax = 8.0;
	base.ode_options.acc = 1e-6;
	base.ode_options.eps = 1e-6;
	base.ode_options.h_initial = 0.05;
	base.ode_options.h_min = 1e-8;
	base.ode_options.max_steps = 160000;
	base.e_min = -0.8;
	base.e_max = -0.2;
	base.scan_points = 120;
	base.root_acc = 1e-10;
	base.root_max_iter = 120;

	const ShootingResult ground = find_ground_state_energy(base);
	const double exact = hydrogen_exact_energy_ground();

	std::printf("  baseline parameters: rmin=%.1e rmax=%.2f acc=%.1e eps=%.1e\n",
			base.rmin,
			base.rmax,
			base.ode_options.acc,
			base.ode_options.eps);
	std::printf("  computed E0 = %.16f  exact = %.16f  |error|=%.3e  mismatch=%.3e  status=%s\n\n",
			ground.energy,
			exact,
			std::abs(ground.energy - exact),
			ground.mismatch,
			ground.converged ? "converged" : "not_converged");

	const std::vector<WavefunctionPoint> wf = compute_wavefunction(base, ground.energy);
	write_wavefunction_data("hydrogen_wavefunction.data", wf);

	// Convergence with respect to rmax.
	const std::vector<double> rmax_values{4.0, 6.0, 8.0, 10.0, 12.0};
	std::vector<ConvergenceRow> rmax_rows;
	rmax_rows.reserve(rmax_values.size());
	for (double rmax : rmax_values) {
		ShootingConfig cfg = base;
		cfg.rmax = rmax;
		rmax_rows.push_back(solve_for_config(cfg, rmax));
	}
	write_convergence_data("hydrogen_convergence_rmax.data", "rmax", rmax_rows);
	print_convergence_table("  Convergence vs rmax:", rmax_rows);

	// Convergence with respect to rmin.
	const std::vector<double> rmin_values{1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4};
	std::vector<ConvergenceRow> rmin_rows;
	rmin_rows.reserve(rmin_values.size());
	for (double rmin : rmin_values) {
		ShootingConfig cfg = base;
		cfg.rmin = rmin;
		rmin_rows.push_back(solve_for_config(cfg, rmin));
	}
	write_convergence_data("hydrogen_convergence_rmin.data", "rmin", rmin_rows);
	print_convergence_table("  Convergence vs rmin:", rmin_rows);

	// Convergence with respect to ODE absolute tolerance.
	const std::vector<double> acc_values{1e-3, 3e-4, 1e-4, 3e-5, 1e-5};
	std::vector<ConvergenceRow> acc_rows;
	acc_rows.reserve(acc_values.size());
	for (double acc : acc_values) {
		ShootingConfig cfg = base;
		cfg.ode_options.acc = acc;
		acc_rows.push_back(solve_for_config(cfg, acc));
	}
	write_convergence_data("hydrogen_convergence_acc.data", "acc", acc_rows);
	print_convergence_table("  Convergence vs ODE acc:", acc_rows);

	// Convergence with respect to ODE relative tolerance.
	const std::vector<double> eps_values{1e-3, 3e-4, 1e-4, 3e-5, 1e-5};
	std::vector<ConvergenceRow> eps_rows;
	eps_rows.reserve(eps_values.size());
	for (double eps : eps_values) {
		ShootingConfig cfg = base;
		cfg.ode_options.eps = eps;
		eps_rows.push_back(solve_for_config(cfg, eps));
	}
	write_convergence_data("hydrogen_convergence_eps.data", "eps", eps_rows);
	print_convergence_table("  Convergence vs ODE eps:", eps_rows);
}

// Packs Newton solve statistics into a LineSearchRow for table output.
LineSearchRow make_line_search_row(const std::string& problem, const std::string& method, const NewtonResult& r) {
	return LineSearchRow{
			problem,
			method,
			r.stats.iterations,
			r.stats.f_evaluations,
			r.stats.jacobian_evaluations,
			r.stats.line_search_iterations,
			r.fx.norm(),
	};
}

// Part C: compare backtracking vs quadratic-interpolation line search on
// Rosenbrock and Himmelblau, measuring the total number of f-evaluations
// and line-search steps.  Also verifies the Jacobian-reuse policy.
void run_part_c() {
	std::puts("=== Part C (1 point): Quadratic interpolation line-search + Jacobian reuse ===");

	NewtonOptions bt;
	bt.acc = 1e-10;
	bt.alpha_min = 1e-3;
	bt.max_iter = 120;
	bt.line_search = LineSearchType::backtracking;

	NewtonOptions qi = bt;
	qi.line_search = LineSearchType::quadratic_interpolation;

	std::vector<LineSearchRow> rows;

	const vector rosen_start{1.2, 1.2};
	const NewtonResult rosen_bt = newton_solve(rosenbrock_gradient(), rosen_start, bt);
	const NewtonResult rosen_qi = newton_solve(rosenbrock_gradient(), rosen_start, qi);

	rows.push_back(make_line_search_row("rosenbrock", "backtracking", rosen_bt));
	rows.push_back(make_line_search_row("rosenbrock", "quadratic", rosen_qi));

	const vector himmel_start{0.0, 0.0};
	const NewtonResult himmel_bt = newton_solve(himmelblau_gradient(), himmel_start, bt);
	const NewtonResult himmel_qi = newton_solve(himmelblau_gradient(), himmel_start, qi);

	rows.push_back(make_line_search_row("himmelblau", "backtracking", himmel_bt));
	rows.push_back(make_line_search_row("himmelblau", "quadratic", himmel_qi));

	write_line_search_data("line_search_comparison.data", rows);

	std::puts("  Method comparison summary:");
	for (const LineSearchRow& row : rows) {
		std::printf("  %-10s %-12s iter=%-4d f_calls=%-5d J_calls=%-4d ls_steps=%-4d ||f||=%.3e\n",
				row.problem.c_str(),
				row.method.c_str(),
				row.iterations,
				row.f_evaluations,
				row.jacobian_evaluations,
				row.line_search_iterations,
				row.residual_norm);
	}

	std::puts("\n  Jacobian allocation policy: exactly one matrix allocation per Newton solve (updated in-place each iteration).");
	std::puts("");
}

} // namespace

int run() {
	run_part_a();
	run_part_b();
	run_part_c();

	std::puts("Data files written:");
	std::puts("  hydrogen_wavefunction.data");
	std::puts("  hydrogen_convergence_rmax.data");
	std::puts("  hydrogen_convergence_rmin.data");
	std::puts("  hydrogen_convergence_acc.data");
	std::puts("  hydrogen_convergence_eps.data");
	std::puts("  line_search_comparison.data");

	return EXIT_SUCCESS;
}

} // namespace pp

int main() {
	return pp::run();
}
