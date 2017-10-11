#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <io.h>
#include <time.h>
#include <malloc.h>
#include "../params.h"

//#define IN_FILE "../data/init_position.prof"
char IN_FILE[100];
#define IN_DIR_VTU "../result/10102017_52008"
#define OUT_DIR_VTU "../result/10102017_52008"

FILE* fp;
char filename[256];
double *Prs_Series;
// vtuへの変換で利用
int NumberOfParticle;
int ProfFileNumber;
double *Position;
double *Velocity;
double *Pressure;
double *pressave;
int *ParticleType;

int fordebug = 0;

void read_data(int iFile) {
	sprintf_s(filename, "%s/output%05d.prof", IN_DIR_VTU, iFile);
	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d", &NumberOfParticle);
	//printf("NumberOfParticle: %d\n", NumberOfParticle);

	Position = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Velocity = (double*)malloc(sizeof(double)*NumberOfParticle * 3);
	Pressure = (double*)malloc(sizeof(double)*NumberOfParticle);
	pressave = (double*)malloc(sizeof(double)*NumberOfParticle);
	ParticleType = (int*)malloc(sizeof(int)*NumberOfParticle);

	for (int i = 0; i < NumberOfParticle; i++) {
		int a[2];
		double b[8];
		fscanf_s(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]);
		ParticleType[i] = a[1];
		Position[i * 3] = b[0];	Position[i * 3 + 1] = b[1];	Position[i * 3 + 2] = b[2];
		Velocity[i * 3] = b[3];	Velocity[i * 3 + 1] = b[4];	Velocity[i * 3 + 2] = b[5];
		Pressure[i] = b[6];
		pressave[i] = b[7];
	}
	fclose(fp);
}

void setProfFileNum() {

	ProfFileNumber = 0;

	struct _finddata_t c_file;
	long   hFile;

	// 最初のファイルを探索
	sprintf_s(filename, "%s/output*.prof", IN_DIR_VTU);
	if ((hFile = _findfirst(filename, &c_file)) == -1L)
		printf("ファイルは存在しません。\n");
	else {
		ProfFileNumber++;	//printf("%s\n", c_file.name);

							// 残りのファイルを探索
		while (_findnext(hFile, &c_file) == 0) {
			ProfFileNumber++;	//printf("%s\n", c_file.name);
		}
		_findclose(hFile);
	}
	printf("fileNum:%i\n", ProfFileNumber);
}

int main() {

	setProfFileNum();
	Prs_Series = (double*)malloc(sizeof(double)*ProfFileNumber);

	for (int iFile = 0; iFile < ProfFileNumber; iFile++) {

		read_data(iFile);
		Prs_Series[iFile] = 0;

		for (int j = 0; j < NumberOfParticle; j++) {
			if (ParticleType[j] == SURFACEWALL) {
				Prs_Series[iFile] += pressave[j] * PARTICLE_DISTANCE*PARTICLE_DISTANCE / 3;
				//printf("%f\n", Prs_Series[iFile]);
			}
		}

		free(Position);
		free(Velocity);
		free(Pressure);
		free(pressave);
		free(ParticleType);
	}
	sprintf_s(filename, "%s/Series_Pressure.csv", IN_DIR_VTU);
	fopen_s(&fp, filename, "w");
	for (int iFile = 0; iFile < ProfFileNumber; iFile++) {
		fprintf(fp, "%lf\n", Prs_Series[iFile]);
	}

	fclose(fp);

	return 0;
}