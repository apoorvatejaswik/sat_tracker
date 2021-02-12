#pragma once



typedef struct { double x, y, z; } v3;
typedef struct { double v[3][3]; } v33;
typedef struct { double az, el; } ae;

v3 operator * (const v33 &M, const v3& V);//matrix vector product
v3 operator - (const v3& V1, const v3& V2);
v33 operator * (const v33 &L, const v33 &R);//left matrix, right matrix product

typedef struct
{
	int NORAD;
	char intldesg[20];
	double MJD;
	double ndot;
	double nddot;
	int nexp;
	double bstar;
	int ibexp;
	double inclo;
	double nodeo;
	double ecco;
	double argpo;
	double mo;
	double no;
} str_tle;

