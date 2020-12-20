#pragma once
#include "Constants.h"
class IHydraulic�omputationStrategy;
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;									// d - �������
	virtual double GetRoughness() const = 0;								// e - �������������
	virtual double GetNormalOilDebit() const = 0;							//q_o - �������� ����� �����
	virtual double GetNormalWaterDebit() const = 0;							//q_w - �������� ����� ����
	virtual double GetNormalGasDebit() const = 0;							//q_g - �������� ����� ����
	virtual double GetTemperature() const = 0;								//t_initial - �����������
	virtual double GetPressure() const = 0;									//p_initial - ��������
	virtual double GetRealInclineAngle() const = 0;							//teta - ����������� ���� ������� �����
	virtual Flow::FlowDirection GetFlowDirection() const = 0;				//���������� ������
	virtual TubeStreamParaeters::StreamType GetStramType() const = 0;		//��� ������
	virtual TubeStreamParaeters::TubeRoughType GetTubeRoughType() const = 0;//��� ������������� �����
	virtual double GetLength() const = 0;									//����� �����

	virtual void SetHydraulicComputationStrategy(IHydraulic�omputationStrategy* straegy) = 0; //��������� ���������-��������� ��������������� ������� ������ � �����
private:
	virtual void NotifyObservingHydraulicComputationStrategy() = 0;		//����������� ��������� �� ���������� ��������� �����
};

