#include "TestBeggsAndBrillReadableTube.h"

double TestBeggsAndBrillReadableTube::GetDiameter() const
{
	return 10.0;
}

double TestBeggsAndBrillReadableTube::GetRoughness() const
{
	return 1.0;
}

double TestBeggsAndBrillReadableTube::GetNormalOilDebit() const
{
	return 1.0;
}

double TestBeggsAndBrillReadableTube::GetNormalWaterDebit() const
{
	return 2.0;
}

double TestBeggsAndBrillReadableTube::GetNormalGasDebit() const
{
	return 3.0;
}

double TestBeggsAndBrillReadableTube::GetTemperature(const int & dot_index, const int & dot_count) const
{
	return 12.0;
}

double TestBeggsAndBrillReadableTube::GetPressure(const int & dot_index, const int & dot_count) const
{
	return 120.0;
}

double TestBeggsAndBrillReadableTube::GetRealInclineAngle() const
{
	return 10.0;
}

Flow::FlowDirection TestBeggsAndBrillReadableTube::GetFlowDirection() const
{
	return Flow::FlowDirection();
}

TubeStreamParaeters::StreamType TestBeggsAndBrillReadableTube::GetStramType() const
{
	return TubeStreamParaeters::StreamType();
}

TubeStreamParaeters::TubeRoughType TestBeggsAndBrillReadableTube::GetTubeRoughType() const
{
	return TubeStreamParaeters::TubeRoughType();
}

double TestBeggsAndBrillReadableTube::GetLength() const
{
	return 100.0;
}

double TestBeggsAndBrillReadableTube::GetPwh() const
{
	return 10.0;
}

double TestBeggsAndBrillReadableTube::GetHdyn() const
{
	return 20.0;
}
