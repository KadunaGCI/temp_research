#pragma once
/*****************************************************************

dam_break_test�ɂ�����σp�����[�^ �����œn������

*****************************************************************/

//#define WAVE_HEIGHT 0.5
//#define CENTER_CUBE_X 0.5
//#define DNS_RIGID0 700	//���̖��x	

/*****************************************************************

	mk_particle.cpp�ŗ��p

*****************************************************************/

//���́i�ƍ��́j�͈̔�
#define mk_MIN_X  0.0
#define mk_MIN_Y  0.0
#define mk_MIN_Z  0.0
#define mk_MAX_X  2.8
#define mk_MAX_Y  0.3
#define mk_MAX_Z  0.2

#define GHOST -1
#define DUMMY 0
#define WALL  1
#define FRONTWALL 2
#define SMWALL 3
#define FLUID 4
#define RIGID0 5
#define NUM_TYP  7		//���q�̎�ސ�

#define RE  3  //3 4 5		PARTICLE_DISTANCE*(RE+0.1)���e�����a�ł��邱�Ƃɒ���

//#define WAVE_HEIGHT 0.5
//#define WAVE_WIDTH 0.25
// ���s����mk_MAX_Y

#define PARTICLE_DISTANCE 0.01 

#define LENGTH_CUBE_X 0.02
#define LENGTH_CUBE_Y 0.12
#define LENGTH_CUBE_Z 0.02
#define CENTER_CUBE_X 0.5
#define CENTER_CUBE_Y mk_MAX_Y/2
#define CENTER_CUBE_Z LENGTH_CUBE_Z/2


/*****************************************************************

	emps.cpp�ŗ��p

*****************************************************************/

//#define PARTICLE_DISTANCE 0.002					//���ϗ��q�ԋ���

#define MIN_X  (mk_MIN_X - PARTICLE_DISTANCE*4)	//��͗̈��x�����̍ŏ��l
#define MIN_Y  (mk_MIN_Y - PARTICLE_DISTANCE*4 - PARTICLE_DISTANCE/2)	//��͗̈��y�����̍ŏ��l and for periodic boundary condition
#define MIN_Z  (mk_MIN_Z - PARTICLE_DISTANCE*4)	//��͗̈��z�����̍ŏ��l
#define MAX_X  (mk_MAX_X + PARTICLE_DISTANCE*4)	//��͗̈��x�����̍ő�l
#define MAX_Y  (mk_MAX_Y + PARTICLE_DISTANCE*4 + PARTICLE_DISTANCE/2)	//��͗̈��y�����̍ő�l
#define MAX_Z  (mk_MAX_Z + PARTICLE_DISTANCE*15)	//��͗̈��z�����̍ő�l

/*
#define GHOST -1
#define FLUID 0
#define WALL  1
#define RIGID0 2
*/

//#define NUM_TYP  4		//���q�̎�ސ�

#define DNS_FLUID 1000		//���̗��q�̖��x
#define DNS_WALL 1000		//�Ǘ��q�̖��x
//#define DNS_RIGID0 500	//���̖��x	

#define DT 0.0004			//���ԍ��ݕ�
#define dt_inv   double (1/DT)	
#define FIN_TIM 4.0  //���Ԃ̏��
#define SND 22			//����
#define OPT_FQC 100		//�o�͊Ԋu�����߂锽����
#define KNM_VSC_FRUID 0.000001	//���S���W��
//#define KNM_VSC_WALL 0.000001	//���S���W��
//#define KNM_VSC_RIGID0 0.000001	//���S���W��
#define DIM 3				//������
#define CRT_NUM 0.1		//�N�[����������
#define COL_RAT 0.2		//�ڋ߂������q�̔�����
#define DST_LMT_RAT 0.9	//����ȏ�̗��q�Ԃ̐ڋ߂������Ȃ������̌W��
#define G_X 0.0			//�d�͉����x��x����
#define G_Y 0.0			//�d�͉����x��y����
#define G_Z -9.8			//�d�͉����x��z����
#define WEI(dist, re) ((re/dist) - 1.0)	//�d�݊֐�

#define CORRECTION MAX_Y - MIN_Y

