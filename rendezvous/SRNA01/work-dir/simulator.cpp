/* AUTO-GENERATED SIMULATOR BY C2E2 */

# include <iostream>
# include <vector>
# include <boost/numeric/odeint.hpp>
# include <math.h>

using namespace std;
using namespace boost::numeric::odeint;

typedef vector<double> state_t;

//INTEGRATOR OBSERVER
class IntObs {
	private:
		vector<state_t> &io_states;
		vector<double> &io_times;

	public:
		IntObs(vector<state_t> &states, vector<double> &times)
			: io_states(states), io_times(times) { }

		void operator()(const state_t &x, double t) {
			io_states.push_back(x);
			io_times.push_back(t);
		}
};

//ODE FUNCTIONS
void ProxA(const state_t &x, state_t &dxdt, const double t) {
			dxdt[0]=x[2];
			dxdt[1]=x[3];
			dxdt[2]=-0.0575997658817729*x[0] + 0.000200959896519766*x[1] - 2.89995083970656*x[2] + 0.00877200894463775*x[3];
			dxdt[3]=-0.000174031357370456*x[0] - 0.0665123984901026*x[1] - 0.00875351105536225*x[2] - 2.90300269286856*x[3];
			dxdt[4]=1;
			dxdt[5]=0;
}

void ProxB(const state_t &x, state_t &dxdt, const double t) {
			dxdt[0]=x[2];
			dxdt[1]=x[3];
			dxdt[2]=-0.575999943070835*x[0] + 0.000262486079431672*x[1] - 19.2299795908647*x[2] + 0.00876275931760007*x[3];
			dxdt[3]=-0.000262486080737868*x[0] - 0.575999940191886*x[1] - 0.00876276068239993*x[2] - 19.2299765959399*x[3];
			dxdt[4]=1;
			dxdt[5]=0;
}

//ODE FUNCTION POINTER
void (*rhs[2])(const state_t &x, state_t &dxdt, const double t) =
	{ProxA, ProxB};

int main() {
	//VARIABLES
	double ts, dt, te;
	double abs_err, rel_err;
	int cur_mode;
	state_t x(6);
	vector<double> times;
	vector<state_t> trace;

	//PARSING CONFIG
	cin >> ts;
	for (int i = 0; i < 6; i++) {
		cin >> x[i];
	}
	cin >> abs_err >> rel_err >> dt >> te >> cur_mode;
	cur_mode--;

	//INTEGRATING
	runge_kutta4<state_t> stepper;
	size_t steps = integrate_const(stepper, rhs[cur_mode], x, ts, te, dt,	IntObs(trace, times));

	//PRINTING STEPS
	for (size_t i = 0; i <= steps; i++) {
		cout << fixed;
		cout << setprecision(9) << times[i];
		for (int j = 0; j < 6; j++) {
			cout << setprecision(10) << ' ' << trace[i][j];
		}
		cout << endl;

		if (i != 0 && i != steps) {
			cout << fixed;
			cout << setprecision(9) << times[i];
			for (int j = 0; j < 6; j++) {
				cout << setprecision(10) << ' ' << trace[i][j];
			}
			cout << endl;
		}
	}
}