#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm
{
public:
	double Execute();
protected:
	std::string file_name;

	BeggsAndBrillAlgorithm();
	virtual ~BeggsAndBrillAlgorithm();
private:
	//��������� �����
	struct PipeParameters {

		double d;			//�������
		double e;			//�������������

	} pipe_parameters;

	//���������� � ������
	struct FlowInfo {
		std::string name;	//�������� (�����, ����, ���)
		double q;			//�������� �����
		double q_normal;	//�������� ����� ��� �.�.
		double v_reduced;	//���������� �������� (����������� � ����� = 0.0) 
							//(v_sl � v_sg � ������� � ���� ��������������)

	};
	FlowInfo oilFlow{ "oil" };
	FlowInfo waterFlow{ "water" };
	FlowInfo gasFlow{ "gas" };


	double t_initial;		// ����������� � �������� �����
	double p_initial;		//�������� � �������� �����

	double volume_coefficient_oil;				//�������� ����������� ����� (B_o)
	double volume_coefficient_water;			//�������� ����������� ���� (B_w)
	double volume_coefficient_gas;				//�������� ����������� ���� (B_w)

	double solubility_oil;			//������������� ���� � ����� (R_s)
	double solubility_water;		//������������� ���� � ���� (R_sw)

	double mixture_speed;			//�������� ����� (v_m)

	//���������� ��������� �������� (���������� ����� � �������� ������� ������)
	//(�� ����������������)
	virtual void ReadInitialParameters(const IReadableBeggsAndBrill& tube) final;

	//������ ��������� ������������ ����� � ����
	virtual void VolumeCoefficientCount(const IPVTBeggsAndBrill& PVT);

	//������ ������������� ���� � ����� � � ����
	virtual void SolubilityCount(const IPVTBeggsAndBrill& PVT);

	//������ �������� ������� �����, ���� � ���� � �������� �����
	virtual void VolumeCount();

	//���������� ����������� �������� �������� � ����
	//���������� �������� �����
	//���������� ��������� ���������� ��������
	void ReducedSpeedCount();

	
};

