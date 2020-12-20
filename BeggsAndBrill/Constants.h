#pragma once
namespace PhysicConstants
{
	const double g = 9.80665;
}

namespace Flow
{
	enum class FlowDirection {
		DownUp,
		UpDown,
	};
}

namespace TubeStreamParaeters
{
	enum class StreamType {
		Laminar,
		Turbulent
	};

	enum class TubeRoughType {
		Smooth,
		Rough
	};
}