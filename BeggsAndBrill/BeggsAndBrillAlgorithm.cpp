#include "BeggsAndBrillAlgorithm.h"


void BeggsAndBrillAlgorithm::Execute(const IReadableBeggsAndBrill& tube)
{
	FlowBorders flow_borders;		//границы, расчитанные режимом потока	
	double H_L_0;					//объёмное содержание жидкости в трубе
	double N_R_e;					//число Рейнольдса

	PipeParameters pipe_parameters;
	ReadPipeParameters(tube,
		pipe_parameters.diameter, pipe_parameters.roughness,
		pipe_parameters.temperature, pipe_parameters.pressure);

	FlowInfo oilFlow{ "oil" };
	FlowInfo waterFlow{ "water" };
	FlowInfo gasFlow{ "gas" };

	ReadFlowNormalDebits(tube,
		oilFlow.volume_debit_normal,
		waterFlow.volume_debit_normal,
		gasFlow.volume_debit_normal);

}

void BeggsAndBrillAlgorithm::Update()
{

}

void BeggsAndBrillAlgorithm::ReadPipeParameters(const IReadableBeggsAndBrill & tube,
	double & out_diameter, double & out_roughness, double & out_temperature, double & out_pressure)
{
	out_diameter = tube.GetDiameter();
	out_roughness = tube.GetRoughness();
	out_temperature = tube.GetTemperature();
	out_pressure = tube.GetPressure();
}

void BeggsAndBrillAlgorithm::ReadFlowNormalDebits(const IReadableBeggsAndBrill & tube, 
	double & out_oil_normal_debit, double & out_water_normal_debit, double & out_gas_normal_debit)
{
	out_oil_normal_debit = tube.GetNormalOilDebit();
	out_water_normal_debit = tube.GetNormalWaterDebit();
	out_gas_normal_debit = tube.GetNormalGasDebit();
}

void BeggsAndBrillAlgorithm::VolumeCoefficientCount(const IPVTBeggsAndBrill & PVT,
	const PipeParameters& pipe_parameters, 
	double & out_volume_coefficient_oil, double & out_volume_coefficient_water, double & out_volume_coefficient_gas)
{
	const double& temperature = pipe_parameters.temperature;
	const double& pressure = pipe_parameters.pressure;
	out_volume_coefficient_oil = PVT.OilVolumeCoefficientCount(temperature,pressure);
	out_volume_coefficient_water = PVT.WaterVolumeCoefficientCount(temperature, pressure);
	out_volume_coefficient_gas = PVT.GasVolumeCoefficientCount(temperature, pressure);
}

void BeggsAndBrillAlgorithm::SolubilityCount(const IPVTBeggsAndBrill & PVT, 
	const PipeParameters& pipe_parameters,
	double & out_gas_sollubility_oil, double & out_gas_sollubility_water)
{
	const double& temperature = pipe_parameters.temperature;
	const double& pressure = pipe_parameters.pressure;
	out_gas_sollubility_oil = PVT.OilSolubilityCount(temperature, pressure);
	out_gas_sollubility_water = PVT.WaterSolubilityCount(temperature, pressure);
}

void BeggsAndBrillAlgorithm::VolumeDebitCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
	double& out_volume_debit_oil, double& out_volume_debit_water, double& out_volume_debit_gas)
{
	out_volume_debit_oil = oilFlow.volume_debit_normal*oilFlow.volume_coefficient / 86400;
	out_volume_debit_water = waterFlow.volume_debit_normal*waterFlow.volume_coefficient / 86400;
	out_volume_debit_gas = (gasFlow.volume_debit_normal - oilFlow.volume_debit_normal*oilFlow.gas_solubility - waterFlow.volume_debit_normal*waterFlow.gas_solubility)*gasFlow.volume_coefficient / 86400;
}

double BeggsAndBrillAlgorithm::ApCount(const double& diameter)
{
	return M_PI * pow(diameter, 2)*1.0e-5 / 4;
}

void BeggsAndBrillAlgorithm::ReducedSpeedCount(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
	const PipeParameters& pipe_parameters,
	double& out_reduced_speed_liquid, double& out_reduced_speed_gas)
{
	double q_l = (oilFlow.volume_debit + waterFlow.volume_debit);
	double A_p = ApCount(pipe_parameters.diameter);
	out_reduced_speed_liquid = q_l / A_p;
	out_reduced_speed_gas = (gasFlow.volume_debit / A_p);
}

void BeggsAndBrillAlgorithm::MixtureSpeedCount(const FlowInfo & waterFlow, const FlowInfo & gasFlow, double & out_mixture_speed)
{
	out_mixture_speed = waterFlow.reduced_speed + gasFlow.reduced_speed;
}

void BeggsAndBrillAlgorithm::VolumeContainityCount(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, double & out_lambda)
{
	double q_l = oilFlow.volume_debit + waterFlow.volume_debit;
	out_lambda = q_l / (q_l + gasFlow.volume_debit);
}

void BeggsAndBrillAlgorithm::ViscosityCount(const IPVTBeggsAndBrill& PVT,
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
	double& out_viscosity_oil, double& out_viscosity_water, double& out_viscosity_gas)
{
	out_viscosity_oil = PVT.OilViscosityCount(42);
	out_viscosity_water = PVT.WaterViscosityCount(42);
	out_viscosity_gas = PVT.GasViscosityCount(42);
}

void BeggsAndBrillAlgorithm::ToughnessCount(const IPVTBeggsAndBrill & PVT, 
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	double & out_toughness_oil, double & out_toughness_water, double & out_toughness_gas)
{
	out_toughness_oil = PVT.OilToughnessCount(42);
	out_toughness_water = PVT.WaterToughnessCount(42);
	out_toughness_gas = PVT.GasToughnessCount(42);
}

void BeggsAndBrillAlgorithm::SurfaceTensionCount(const IPVTBeggsAndBrill & PVT,
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	double & out_surface_tension_oil, double & out_surface_tension_water)
{
	out_surface_tension_oil = PVT.OilSurfaceTensionCount(42);
	out_surface_tension_water = PVT.WaterSurfaceTensionCount(42);
}

void BeggsAndBrillAlgorithm::MixtureViscosity_Toughness_SurfaceTension_Count(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	const MixtureInfo& mixture,
	double & out_mixture_viscosity_secondphase, double & out_mixture_toughness_secondphase,
	double & out_mixture_viscosity_liquid, double & out_mixture_toughness_liquid, double & out_mixture_surfacetansion_liquid)
{
	double oil_proportion = oilFlow.volume_debit / (oilFlow.volume_debit + waterFlow.volume_debit);	//f_o - доля нефти в жидкой фазе
	out_mixture_viscosity_liquid = 1.0e-2 * (oilFlow.viscosity*oil_proportion + waterFlow.viscosity*(1 - oil_proportion));
	out_mixture_viscosity_secondphase = out_mixture_viscosity_liquid * mixture.lambda + 1.0e-2 * gasFlow.viscosity*(1 - mixture.lambda);
	out_mixture_toughness_liquid = oilFlow.toughness*oil_proportion + waterFlow.toughness*(1 - oil_proportion);
	out_mixture_toughness_secondphase = out_mixture_toughness_liquid * mixture.lambda + gasFlow.toughness*(1 - mixture.lambda);
	out_mixture_surfacetansion_liquid = oilFlow.surface_tension*oil_proportion + waterFlow.surface_tension*(1 - oil_proportion);
}

void BeggsAndBrillAlgorithm::FroudNumberCount(const MixtureInfo & mixture, const PipeParameters & pipe_parameters, double & out_froud_number)
{
	out_froud_number = pow(mixture.speed, 2) / PhysicConstants::g*pipe_parameters.diameter * 1.0e-2;
}

void BeggsAndBrillAlgorithm::LiquidIndexCount(const MixtureInfo & mixture, double & out_liquid_speed_index)
{
	out_liquid_speed_index = mixture.reduced_speed*pow((mixture.toughness_liquid / PhysicConstants::g*mixture.surface_tension_liquid), 1 / 4.f);
}

void BeggsAndBrillAlgorithm::FlowBordersCount(const MixtureInfo & mixture, double & out_L1, double & out_L2, double & out_L3, double & out_L4)
{
	out_L1 = 316 * pow(mixture.lambda, 0.302);
	out_L2 = 0.000925 * pow(mixture.lambda, -2.468);
	out_L3 = 0.1 * pow(mixture.lambda, -1.452);
	out_L4 = 0.5 * pow(mixture.lambda, -6.738);
}

void BeggsAndBrillAlgorithm::VolumeCorrectionLiquidCount(const IReadableBeggsAndBrill & tube, const MixtureInfo & mixture, const FlowBorders& borders,
	double& out_volume_correction_liquid)
{
	const MixtureInfo * ptr_mixture = &mixture;
	std::vector<FlowState*> flowMods;
	flowMods.reserve(3);
	flowMods.push_back(new DividedFlowState(ptr_mixture));
	flowMods.push_back(new IntermittentFlowState(ptr_mixture));
	flowMods.push_back(new DistributedFlowState(ptr_mixture));

	//ищем среди них подохдящий режим
	FlowState* correctMode = nullptr;
	for (auto ptr_flow : flowMods) {
		if (ptr_flow->ModeChecking(borders)) {
			correctMode = ptr_flow;
			break;
		}
	}
	if (correctMode == nullptr) {
		std::cout << "Не обнаржуено подохдящих модов, выбран самый первый BeggsAndBrillAlgorithm::VolumeCorrectionLiquid()"<<std::endl;
		correctMode = flowMods[0];
	}

	out_volume_correction_liquid = correctMode->GetVolumetric_content_liquid_angled(tube);
}

void BeggsAndBrillAlgorithm::PayneVolumeCorrectionCount(const IReadableBeggsAndBrill & tube, const double & volume_correction_liquid, double & out_payne_corrected_vlume_correction_liquid)
{
	double teta = tube.GetRealInclineAngle();
	if (teta > 0) {
		out_payne_corrected_vlume_correction_liquid = 0.924*volume_correction_liquid;
	}
	else if(teta<0)
	{
		out_payne_corrected_vlume_correction_liquid = 0.685*volume_correction_liquid;
	}
	else {
		std::cout << "Не предусмотренный алгоритм вариант! Ошибка либо в алгоритме, либо в коде! BeggsAndBrillAlgorithm::PayneVolumeCorrection()" << std::endl;
	}
}

void BeggsAndBrillAlgorithm::RaynoldsNumberCount(const MixtureInfo & mixture, const PipeParameters & pipe_parameters, 
	double& out_raynolds_number)
{
	out_raynolds_number = mixture.toughness_secondphase*mixture.speed*pipe_parameters.diameter*1.0e-2 / mixture.viscosity_secondphase;

}

void BeggsAndBrillAlgorithm::NormalizingFrictionCoefficientCount(const IReadableBeggsAndBrill & tube, const PipeParameters & pipe_parameters, 
	const double & raynolds_number, 
	double& out_normalizing_friction_coefficient)
{
	using namespace TubeStreamParaeters;
	StreamType stream_type = tube.GetStramType();
	TubeRoughType rough_type = tube.GetTubeRoughType();
	
	switch (stream_type) {
	case StreamType::Laminar:
		out_normalizing_friction_coefficient = 64 * 1 / raynolds_number;
		break;
	case StreamType::Turbulent:
		switch (rough_type) {
		case TubeRoughType::Smooth:
			if (raynolds_number < 1.0e+6) {
				out_normalizing_friction_coefficient = 0.316*pow(raynolds_number, -0.25);
			}
			else {
				out_normalizing_friction_coefficient = 0.0056 + 0.5*pow(raynolds_number, -0.32);
			}
			break;
		case TubeRoughType::Rough:
			double f = 0.08;
			const int UNDECLARED_ITERATION_COUNT = 10;//количество итераций неизвестно (я не смог найти)
			for (int i = 0; i < UNDECLARED_ITERATION_COUNT; ++i) {
				f = 1 / (4 * pow(log(pipe_parameters.roughness / (3.7*pipe_parameters.diameter) + 2.51 / (raynolds_number*sqrt(f))),2));
			}
			out_normalizing_friction_coefficient = f;
			break;
		}
		break;
	}
}

double BeggsAndBrillAlgorithm::GetY(const double & lambda, const double & volume_correction_liquid)
{
	return lambda / pow(volume_correction_liquid, 2);
}

double BeggsAndBrillAlgorithm::GetS(const double & y)
{
	double log_y = log(y);

	return log_y / (-0.0523 + 3.182*log_y - 0.8725*pow(log_y, 2) + 0.01853*(pow(log_y, 4)));
}

void BeggsAndBrillAlgorithm::FrictionCoefficientCount(const MixtureInfo & mixture, const double & volume_correction_liquid,
	const double & normilizing_friction_coefficient, 
	double & out_friction_coefficient)
{
	double y = GetY(mixture.lambda, volume_correction_liquid);
	double s = GetS(y);
	out_friction_coefficient = normilizing_friction_coefficient * pow(std::exp(1.0), s);
}

void BeggsAndBrillAlgorithm::PressureGradientCount(const IReadableBeggsAndBrill & tube,
	const PipeParameters & pip_parameters, const MixtureInfo & mixture, const FlowInfo& gasFlow,
	const double & friction_coefficient, const double & volume_correction_liquid, 
	double & out_pressure_gradient)
{
	double E_k = (mixture.speed*gasFlow.reduced_speed*mixture.toughness_secondphase) /
		(1.0e+7*pip_parameters.pressure);
	double ro_s = mixture.toughness_liquid*volume_correction_liquid + gasFlow.toughness*(1 - volume_correction_liquid);
	double temp = (friction_coefficient*mixture.toughness_secondphase*pow(mixture.speed, 2)) / (2 * pip_parameters.diameter*1.0e-2);
	out_pressure_gradient = 1.0e-5*(temp + ro_s * PhysicConstants::g*sin(tube.GetRealInclineAngle())) / (1 - E_k);
}





//FloState
BeggsAndBrillAlgorithm::FlowState::FlowState(const MixtureInfo * ptr_mixture_info)
{
	ptr_mixture = ptr_mixture_info;

	//все режимы в напр. сверху вниз
	e_up_down = 4.7;
	f_up_down = -0.3692;
	g_up_down = 0.1244;
	h_up_down = -0.5056;
}

BeggsAndBrillAlgorithm::FlowState::~FlowState()
{
	ptr_mixture = nullptr;
}

double BeggsAndBrillAlgorithm::FlowState::GetVolumetric_content_liquid_angled(const IReadableBeggsAndBrill & tube)
{
	double H_L_0 = GetVolumetric_content_liquid_not_angled();
	double psi = GetAngleCorrectionCoefficient(tube);
	double result = H_L_0 * psi;
	return result;
}

double BeggsAndBrillAlgorithm::FlowState::GetVolumetric_content_liquid_not_angled()
{
	 double result = a * pow(ptr_mixture->lambda, b) / pow(ptr_mixture->froud_number, c);
	 //assert(result >= ptr_mixture->lambda && "Не выполняется ограничение объёмного содержание, видимо что-то не так в расчётах(коде)...");
	 return result;
}

double BeggsAndBrillAlgorithm::FlowState::GetC_Coefficient(const Flow::FlowDirection& flow_direction)
{
	double e, f, g, h;						//параметры для расчёта C
	switch (flow_direction) {
	case Flow::FlowDirection::DownUp:
		e = e_down_up;
		f = f_down_up;
		g = g_down_up;
		h = h_down_up;
		break;
	case Flow::FlowDirection::UpDown:
		e = e_up_down;
		f = f_up_down;
		g = g_up_down;
		h = h_up_down;
	}
	

	double result = (1 - ptr_mixture->lambda)*log(e*
		pow(ptr_mixture->lambda, f)*
		pow(ptr_mixture->liquid_speed_index, g)*
		pow(ptr_mixture->froud_number, h));

	return result;
}

double BeggsAndBrillAlgorithm::FlowState::GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube)
{
	Flow::FlowDirection flow_direction = tube.GetFlowDirection();
	double teta = tube.GetRealInclineAngle();
	double sin_repl = sin(1.8*teta);
	double result = 1.0 + GetC_Coefficient(flow_direction)*(sin_repl - 0.333*pow(sin_repl, 3));

	return result;
}

//DividedFlowState
BeggsAndBrillAlgorithm::DividedFlowState::DividedFlowState(const MixtureInfo * ptr_mixture_info) : FlowState(ptr_mixture_info)
{
	a = 0.980;
	b = 0.4846;
	c = 0.0868;
	e_down_up = 0.011;
	f_down_up = -3.7680;
	g_down_up = 3.5390;
	h_down_up = -1.6140;
}

BeggsAndBrillAlgorithm::DividedFlowState::~DividedFlowState()
{
	//ничего на текущий момент
}

bool BeggsAndBrillAlgorithm::DividedFlowState::ModeChecking(const FlowBorders & borders)
{
	const double& lambda = ptr_mixture->lambda;
	const double& froud_number = ptr_mixture->froud_number;
	return lambda < 0.01 && froud_number < borders.L1 ||
		lambda >= 0.01 && froud_number <= borders.L2;
}

//IntermittentFlowState
BeggsAndBrillAlgorithm::IntermittentFlowState::IntermittentFlowState(const MixtureInfo * ptr_mixture_info) : FlowState(ptr_mixture_info)
{
	a = 0.845;
	b = 0.5351;
	c = 0.0173;
	e_down_up = 2.960;
	f_down_up = 0.3050;
	g_down_up = -0.4473;
	h_down_up = 0.0978;
}

BeggsAndBrillAlgorithm::IntermittentFlowState::~IntermittentFlowState()
{
	//ничего на текущий момент
}

bool BeggsAndBrillAlgorithm::IntermittentFlowState::ModeChecking(const FlowBorders & borders)
{
	const double& lambda = ptr_mixture->lambda;
	const double& froud_number = ptr_mixture->froud_number;
	return ((0.01 <= lambda && lambda < 0.4) &&
		(borders.L3 < froud_number && froud_number <= borders.L1)) ||
		((lambda >= 0.4) &&
		(borders.L3 < froud_number && froud_number <= borders.L4));
}

//DistributedFlowState
BeggsAndBrillAlgorithm::DistributedFlowState::DistributedFlowState(const MixtureInfo * ptr_mixture_info) : FlowState(ptr_mixture_info)
{
	a = 0.845;
	b = 0.5351;
	c = 0.0173;
;
}

BeggsAndBrillAlgorithm::DistributedFlowState::~DistributedFlowState()
{
	//ничего на текущий момент
}

bool BeggsAndBrillAlgorithm::DistributedFlowState::ModeChecking(const FlowBorders & borders)
{
	const double& lambda = ptr_mixture->lambda;
	const double& froud_number = ptr_mixture->froud_number;
	return lambda<0.4 && froud_number >= borders.L1 ||
		lambda >= 0.4 && froud_number>borders.L4;
}

double BeggsAndBrillAlgorithm::DistributedFlowState::GetC_Coefficient(const Flow::FlowDirection & flow_direction)
{
	return 0.0;
}

double BeggsAndBrillAlgorithm::DistributedFlowState::GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube)
{
	return 1.0;
}
