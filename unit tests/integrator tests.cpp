#include "catch.hpp"
#include "integrator.h"

using namespace Waves;

// Some basic sample unit testing using Catch https://github.com/philsquared/Catch

TEST_CASE("Initialise() constructs an appropriate number of cells, position and adjacency information.", "[Integrator::Initialise]") {
	g_waves.Initialise(2, 2, 1e-2, 1, 1);

	REQUIRE(g_waves.cells.size() == g_waves.position_information.size());
	REQUIRE(g_waves.cells.size() == g_waves.adjacency_information.size());
}
