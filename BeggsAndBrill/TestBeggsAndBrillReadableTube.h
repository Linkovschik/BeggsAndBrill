#pragma once
#include "IReadableBeggsAndBrill.h"
class TestBeggsAndBrillReadableTube :
	public IReadableBeggsAndBrill
{
	// Унаследовано через IReadableBeggsAndBrill
	virtual double GetDiameter() const override;
	virtual double GetRoughness() const override;
	virtual double GetNormalOilDebit() const override;
	virtual double GetNormalWaterDebit() const override;
	virtual double GetNormalGasDebit() const override;
	virtual double GetTemperature(const int & dot_index, const int & dot_count) const override;
	virtual double GetPressure(const int & dot_index, const int & dot_count) const override;
	virtual double GetRealInclineAngle() const override;
	virtual Flow::FlowDirection GetFlowDirection() const override;
	virtual TubeStreamParaeters::StreamType GetStramType() const override;
	virtual TubeStreamParaeters::TubeRoughType GetTubeRoughType() const override;
	virtual double GetLength() const override;
	virtual double GetPwh() const override;
	virtual double GetHdyn() const override;
};

