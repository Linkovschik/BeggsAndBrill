#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "IHydraulic�omputationStrategy.h"
#include "Constants.h"
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm : public IHydraulic�omputationStrategy
{
public:
	//Execute - ��������� �����, �� ������ ����������������
	virtual void Execute(const IReadableBeggsAndBrill& tube) override;
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
		double reduced_speed;					//���������� �������� (v_sg � ����)
		double viscosity;						//�������� (mu)
		double toughness;						//��������� (ro)
		double surface_tension;					//������������� ��������� (delt)
	};

	//���������� � ������ ������, ����������� � ���������
	struct MixtureInfo {

		double reduced_speed;					//���������� �������� (v_sl � ��������)
		double speed;							//�������� ����� (v_m)
		double lambda;							//�������� ���������� �������� (lambda)
		double viscosity_liquid;				//�������� ����� (mu_l)
		double toughness_liquid;				//��������� ����� (ro_l)
		double surface_tension_liquid;			//������������� ��������� ����� (delt)
		double viscosity_secondphase;			//�������� �����  (��������.) (mu_n)
		double toughness_secondphase;			//��������� ����� (��������) (ro_n) 
		double froud_number;					//����� ����� (N_fr)
		double liquid_speed_index;				//���������� �������� �������� (N_lv)
	};

	//������ ������
	struct FlowBorders {
		double L1, L2, L3, L4;
	};

	//��������� ������
	class FlowState {
	public:
		//�����������
		FlowState(const MixtureInfo* ptr_mixture_info);
		//����������
		virtual ~FlowState();

		//���������� �������� ���������� �������� � ��������� �� ���� �������
		double GetVolumetric_content_liquid_angled(const IReadableBeggsAndBrill& tube);

		//���������� true, ���� �������� ������ �������� ��� ����� ������
		virtual bool ModeChecking(const FlowBorders& borders) = 0;
	protected:
		const MixtureInfo* ptr_mixture;							//��������� �� ��������� "�����"
		double a, b, c;											//��������� ��� ������� ��������� ���������� ��������
		double e_down_up, f_down_up, g_down_up, h_down_up;		//��������� ��� ������� ������� ������������ C (�����_�����)
		double e_up_down, f_up_down, g_up_down, h_up_down;		//��������� ��� ������� ������� ������������ C (������_����)
	private:
		//������ ��������� ���������� ������� ��� �������� �� ����
		virtual double GetVolumetric_content_liquid_not_angled();
		//���������� �������� ������������ C
		virtual double GetC_Coefficient(const Flow::FlowDirection& flow_direction);
		//������ ������������ ������������ �� ���� ������� ����� (���)
		virtual double GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube);
	};

	//���������� �����
	class DividedFlowState : public FlowState {
	public:
		//�����������
		DividedFlowState(const MixtureInfo* ptr_mixture_info);
		//����������
		virtual ~DividedFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;
	};

	//����������� �����
	class IntermittentFlowState : public FlowState {
	public:
		//�����������
		IntermittentFlowState(const MixtureInfo* ptr_mixture_info);
		//����������
		virtual ~IntermittentFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;
	};

	//������������� �����
	class DistributedFlowState : public FlowState {
	public:
		//�����������
		DistributedFlowState(const MixtureInfo* ptr_mixture_info);
		//����������
		virtual ~DistributedFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;

		virtual double GetC_Coefficient(const Flow::FlowDirection& flow_direction) override;
		virtual double GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube) override;
	};
private:
	virtual void Update() override;

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
		double& out_reduced_speed_liquid, double& out_reduced_speed_gas);

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
	virtual void LiquidIndexCount(const MixtureInfo& mixture,
		double& out_liquid_speed_index);

	//������ ������ ������� ������
	virtual void FlowBordersCount(const MixtureInfo& mixture,
		double& out_L1, double& out_L2, double& out_L3, double& out_L4);

	//�������� ���������� �������� � ��������� �� ���� �������
	virtual void VolumeCorrectionLiquidCount(const IReadableBeggsAndBrill& tube, const MixtureInfo& mixture, const FlowBorders& borders,
		double& out_volume_correction_liquid);

	//������������ �������� ��� ��������� ���������� �������� (Payne)
	virtual void PayneVolumeCorrectionCount(const IReadableBeggsAndBrill& tube, const double& volume_correction_liquid,
		double &out_payne_corrected_vlume_correction_liquid);

	//���������� ����� ����������
	virtual void RaynoldsNumberCount(const MixtureInfo& mixture, const PipeParameters& pipe_parameters,
		double& out_raynolds_number);

	//�� ������ ������������� � ����� ���������� ���������� �������� ������������ ������������ ������
	virtual void NormalizingFrictionCoefficientCount(const IReadableBeggsAndBrill& tube, const PipeParameters& pipe_parameters, 
		const double& raynolds_number,
		double& out_normalizing_friction_coefficient);

	//��������� �������� y
	double GetY(const double& lambda, const double& volume_correction_liquid);

	//��������� �������� s
	double GetS(const double& y);

	//��������� ����������� ������ ����������� ������
	virtual void FrictionCoefficientCount(const MixtureInfo& mixture, const double& volume_correction_liquid,
		const double& normilizing_friction_coefficient,
		double& out_friction_coefficient);

	//��������� �������� ��������
	virtual void PressureGradientCount(const IReadableBeggsAndBrill & tube, 
		const PipeParameters& pip_parameters, const MixtureInfo& mixture, const FlowInfo& gasFlow,
		const double& friction_coefficient, const double& volume_correction_liquid,
		double& out_pressure_gradient);
}; 
