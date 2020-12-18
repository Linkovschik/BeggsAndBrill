#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include "Constants.h"
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm
{
public:
	//Execute - ��������� �����, �� ������ ����������������
	double Execute(const IReadableBeggsAndBrill& tube);
protected:
	std::string file_name;

	//�����������
	BeggsAndBrillAlgorithm();

	//����������
	virtual ~BeggsAndBrillAlgorithm();

	//��������� ����� (���������� �����)
	struct PipeParameters {

		double diameter;						//������� (d)
		double roughness;						//������������� (e)
		double temperature;						//����������� � ��������� ����� (T_initial)
		double pressure;						//�������� � ��������� ����� (P_initial)
	};

	//���������� � ������ (���������� ����� � ����������� � ���������)
	struct FlowInfo {

		std::string name;						//�������� (�����, ����, ���)
		double volume_debit;					//�������� ����� (q)
		double volume_debit_normal;				//�������� ����� ��� �.�. (q_normal)
		double volume_coefficient;				//�������� ����������� (B_o � B_w � ����� � ���� �����.)
		double gas_solubility;					//������������� ���� (R_s, R_sw � ����� � ���� �����.)
		double speed_reduced;					//���������� �������� (v_sl � v_sg � ������� � ���� �����.)
		double viscosity;						//�������� (mu)
		double toughness;						//��������� (ro)
		double surface_tension;					//������������� ��������� (delt)
	};

	//���������� � ������ ������, ����������� � ���������
	struct MixtureInfo {

		double speed;								//�������� ����� (v_m)
		double lambda;								//�������� ���������� �������� (lambda)
		double viscosity_liquid;					//�������� ����� (mu)
		double toughness_liquid;					//��������� ����� (ro)
		double surface_tension_liquid;				//������������� ��������� ����� (delt)
		double viscosity_secondphase;				//�������� �����  (��������.) (mu)
		double toughness_secondphase;				//��������� ����� (��������) (ro) 
		double froud_number;						//����� ����� (N_fr)
		double liquid_speed_index;					//���������� �������� �������� (N_lv)
	};
private:

	//���������� ��������� �������� (���������� �����)
	virtual void ReadPipeParameters(const IReadableBeggsAndBrill& tube, 
		double& out_diameter, double& out_roughness, double& out_temperature, double& out_pressure);

	//���������� ���������� ������
	virtual void ReadFlowNormalDebits(const IReadableBeggsAndBrill& tube,
		double& out_oil_normal_debit, double& out_water_normal_debit, double& out_gas_normal_debit);

	//������ ��������� ������������ ����� � ����
	virtual void VolumeCoefficientCount(const IPVTBeggsAndBrill& PVT,
		const PipeParameters& pipe_parameters,
		double& out_volume_coefficient_oil, double& out_volume_coefficient_water, double& out_volume_coefficient_gas);

	//������ ������������� ���� � ����� � � ����
	virtual void SolubilityCount(const IPVTBeggsAndBrill& PVT,
		const PipeParameters& pipe_parameters,
		double& out_gas_sollubility_oil, double& out_gas_sollubility_water);

	//������ �������� ������� �����, ���� � ���� � �������� �����
	virtual void VolumeDebitCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		double& out_volume_debit_oil, double& out_volume_debit_water, double& out_volume_debit_gas);

	double ApCount(const double& diameter);

	//���������� ����������� �������� �������� � ����
	virtual void ReducedSpeedCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		const PipeParameters& pipe_parameters,
		double& out_reduced_speed_water, double& out_reduced_speed_gas);

	//���������� �������� �����
	virtual void MixtureSpeedCount(const FlowInfo& waterFlow, const FlowInfo& gasFlow, 
		double& out_mixture_speed);

	//���������� ��������� ���������� ��������
	virtual void VolumeContainityCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		double& out_lambda);

	//������ �������� �����, ����, ����
	virtual void ViscosityCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_viscosity_oil, double& out_viscosity_water, double& out_viscosity_gas);

	//������ ��������� �����, ����, ����
	virtual void ToughnessCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_toughness_oil, double& out_toughness_water, double& out_toughness_gas);

	//������ �������������� ��������� ����� � ����
	virtual void SurfaceTensionCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_surface_tension_oil, double& out_surface_tension_water);

	//��������� ��������, ��������� � ������������� ��������� �����
	virtual void MixtureViscosity_Toughness_SurfaceTension_Count(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		const MixtureInfo& mixture,
		double& out_mixture_viscosity_secondphase, double& out_mixture_toughness_secondphase,
		double& out_mixture_viscosity_liquid, double& out_mixture_toughness_liquid, double& out_mixture_surfacetansion_liquid);

	//������ ����� ����� �� ������ �������� ����� � �������� �����
	virtual void FroudNumberCount(const MixtureInfo& mixture, const PipeParameters& pipe_parameters,
		double& out_froud_number);

	//������ ���������� �������� ��������
	virtual
};

