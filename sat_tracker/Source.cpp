////////////////////////////////////////////////
// Satellite Tracker
// this code serves educational purposes only
//
// Project developers and contributors:
// Daniel Kucharski, since Feb 12 2021, allsky@utexas.edu
//
// 
////////////////////////////////////////////////



#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include <limits>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <limits.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <fstream>
#include <fcntl.h>
#include <sys/types.h>
#include <io.h>
#include <istream>
#include <ostream>
#include <vector>

#include <filesystem>
#include <experimental/filesystem>
#include <iterator>

#define MAXBUF 200
#define MAX_PATH 260
namespace fs = std::experimental::filesystem;


#include "Header.h"



////////////////////////////////////////////////
//Astronomical Algorithms, Jean Meeus
//source: http://www.naughter.com/aa.html
//example functions in AATest.cpp
#include "AA+.h"

#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(x) ((void)(x))
#endif //#ifndef UNREFERENCED_PARAMETER
////////////////////////////////////////////////

v3 operator * (const v33 & M, const v3& V)//matrix vector product
{
	v3 out;
	out.x = M.v[0][0] * V.x + M.v[0][1] * V.y + M.v[0][2] * V.z;
	out.y = M.v[1][0] * V.x + M.v[1][1] * V.y + M.v[1][2] * V.z;
	out.z = M.v[2][0] * V.x + M.v[2][1] * V.y + M.v[2][2] * V.z;
	return out;
}

v3 operator - (const v3& V1, const v3& V2)
{
	v3 Out;
	Out.x = V1.x - V2.x;
	Out.y = V1.y - V2.y;
	Out.z = V1.z - V2.z;
	return Out;
}

v33 operator * (const v33 &L, const v33 &R)//left matrix, right matrix product
{
	v33 M;//out matrix

	M.v[0][0] = L.v[0][0] * R.v[0][0] + L.v[0][1] * R.v[1][0] + L.v[0][2] * R.v[2][0];
	M.v[0][1] = L.v[0][0] * R.v[0][1] + L.v[0][1] * R.v[1][1] + L.v[0][2] * R.v[2][1];
	M.v[0][2] = L.v[0][0] * R.v[0][2] + L.v[0][1] * R.v[1][2] + L.v[0][2] * R.v[2][2];

	M.v[1][0] = L.v[1][0] * R.v[0][0] + L.v[1][1] * R.v[1][0] + L.v[1][2] * R.v[2][0];
	M.v[1][1] = L.v[1][0] * R.v[0][1] + L.v[1][1] * R.v[1][1] + L.v[1][2] * R.v[2][1];
	M.v[1][2] = L.v[1][0] * R.v[0][2] + L.v[1][1] * R.v[1][2] + L.v[1][2] * R.v[2][2];

	M.v[2][0] = L.v[2][0] * R.v[0][0] + L.v[2][1] * R.v[1][0] + L.v[2][2] * R.v[2][0];
	M.v[2][1] = L.v[2][0] * R.v[0][1] + L.v[2][1] * R.v[1][1] + L.v[2][2] * R.v[2][1];
	M.v[2][2] = L.v[2][0] * R.v[0][2] + L.v[2][1] * R.v[1][2] + L.v[2][2] * R.v[2][2];

	return M;
}

//SPG4
typedef struct
{
	int error;
	bool deep_space;
	bool init;


	//Near Earth
	int isimp;
	double aycof, con41, cc1, cc4, cc5, d2, d3, d4,
		delmo, eta, argpdot, omgcof, sinmao, t2cof, t3cof,
		t4cof, t5cof, x1mth2, x7thm1, mdot, nodedot, xlcof, xmcof,
		nodecf;

	//Deep Space
	int irez;
	double d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232,
		d5421, d5433, dedt, del1, del2, del3, didt, dmdt,
		dnodt, domdt, e3, ee2, peo, pgho, pho, pinco,
		plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3,
		si2, si3, sl2, sl3, sl4, gsto, xfact, xgh2,
		xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3,
		xl4, xlamo, zmol, zmos, atime, xli, xni;

	double a, altp, alta, nddot, ndot,
		bstar, rcse, inclo, nodeo, ecco, argpo, mo,
		no;

	double MJD;
	int NORAD;

} elsetrec;

void read_TLE_catalogue(
	char *TLE_filepath,
	std::vector <str_tle> *tle
);

void read_TLE(
	char *TLE_line1,//in
	char *TLE_line2,//in
	str_tle *tle//out
);

double MJD_full_DoY1_d(
	int year,//2006
	double doy1//day of year starting from 1, Jan 1=1doy
);

void date_from_year_and_doy1(//doy1=day of year starting from 1, Jan 1=1doy
	int year,
	int doy1,
	int *month,
	int *day
);

long double MJD_full(
	int year,//2006
	int month,
	int day,
	int hour,
	int minute,
	double second
);

int sgp4_satrec_init(
	str_tle *tle,//in
	elsetrec *satrec//out
);

void getgravconst
(
	double& tumin,
	double& radiusearthkm,
	double& xke,
	double& j2,
	double& j3,
	double& j4,
	double& j3oj2
);

int sgp4init
(
	const double epoch,
	const double xbstar, const double xecco, const double xargpo,
	const double xinclo, const double xmo, const double xno,
	const double xnodeo, elsetrec& satrec
);

void initl
(
	double ecco, double epoch, double inclo, double& no,
	bool &deep_space,
	double& ainv, double& ao, double& con41, double& con42, double& cosio,
	double& cosio2, double& eccsq, double& omeosq, double& posq,
	double& rp, double& rteosq, double& sinio, double& gsto
);

void dscom(
	double epoch, double ep, double argpp, double tc, double inclp,
	double nodep, double np,
	double& snodm, double& cnodm, double& sinim, double& cosim, double& sinomm,
	double& cosomm, double& day, double& e3, double& ee2, double& em,
	double& emsq, double& gam, double& peo, double& pgho, double& pho,
	double& pinco, double& plo, double& rtemsq, double& se2, double& se3,
	double& sgh2, double& sgh3, double& sgh4, double& sh2, double& sh3,
	double& si2, double& si3, double& sl2, double& sl3, double& sl4,
	double& s1, double& s2, double& s3, double& s4, double& s5,
	double& s6, double& s7, double& ss1, double& ss2, double& ss3,
	double& ss4, double& ss5, double& ss6, double& ss7, double& sz1,
	double& sz2, double& sz3, double& sz11, double& sz12, double& sz13,
	double& sz21, double& sz22, double& sz23, double& sz31, double& sz32,
	double& sz33, double& xgh2, double& xgh3, double& xgh4, double& xh2,
	double& xh3, double& xi2, double& xi3, double& xl2, double& xl3,
	double& xl4, double& nm, double& z1, double& z2, double& z3,
	double& z11, double& z12, double& z13, double& z21, double& z22,
	double& z23, double& z31, double& z32, double& z33, double& zmol,
	double& zmos
);

void dpper(
	double e3, double ee2, double peo, double pgho, double pho,
	double pinco, double plo, double se2, double se3, double sgh2,
	double sgh3, double sgh4, double sh2, double sh3, double si2,
	double si3, double sl2, double sl3, double sl4, double t,
	double xgh2, double xgh3, double xgh4, double xh2, double xh3,
	double xi2, double xi3, double xl2, double xl3, double xl4,
	double zmol, double zmos, double inclo,
	bool init,
	double& ep, double& inclp, double& nodep, double& argpp, double& mp
);

void dsinit(
	double cosim, double emsq, double argpo, double s1, double s2,
	double s3, double s4, double s5, double sinim, double ss1,
	double ss2, double ss3, double ss4, double ss5, double sz1,
	double sz3, double sz11, double sz13, double sz21, double sz23,
	double sz31, double sz33, double t, double tc, double gsto,
	double mo, double mdot, double no, double nodeo, double nodedot,
	double xpidot, double z1, double z3, double z11, double z13,
	double z21, double z23, double z31, double z33, double ecco,
	double eccsq, double& em, double& argpm, double& inclm, double& mm,
	double& nm, double& nodem,
	int& irez,
	double& atime, double& d2201, double& d2211, double& d3210, double& d3222,
	double& d4410, double& d4422, double& d5220, double& d5232, double& d5421,
	double& d5433, double& dedt, double& didt, double& dmdt, double& dndt,
	double& dnodt, double& domdt, double& del1, double& del2, double& del3,
	double& xfact, double& xlamo, double& xli, double& xni
);

int sgp4(elsetrec& satrec, double MJD, double r[3], double v[3]);
double gstime(double jdut1);
double atan2_check(double y, double x);

void dspace
(
	int irez,
	double d2201, double d2211, double d3210, double d3222, double d4410,
	double d4422, double d5220, double d5232, double d5421, double d5433,
	double dedt, double del1, double del2, double del3, double didt,
	double dmdt, double dnodt, double domdt, double argpo, double argpdot,
	double t, double tc, double gsto, double xfact, double xlamo,
	double no,
	double& atime, double& em, double& argpm, double& inclm, double& xli,
	double& mm, double& xni, double& nodem, double& dndt, double& nm
);

double GAST_rad(double MJD);

v33 rz(double Angle_rad);
v33 ry(double Angle_rad);
v33 inv(v33 M);//inverse of 3x3 matrix
double det(v33 M);

ae topo_AzEl_rad(//topocentric
	v3 StaTCSm,
	ae StaTCSrad,
	v3 PointTCSm
);

ae get_ae(v3 v);
double IArad(v3 A, v3 B);
v3 normalize(v3 v);
double norm(v3 v);
double dot(v3 left, v3 right);

v3 get_v3(ae Vrad);//get unit vector at this orientation
v3 get_v3(ae Vrad, double Length);

ae Sun_ICS_rad(double MJD);
double Distance_to_the_SUN_AU(double MJD);
double shadow_function(double MJD, v3 Sat_ICS_m, v3 Sun_ICS_m);

#define pi 3.14159265358979323846
const double deg2rad = pi / 180.0;
const double rad2deg = 180.0 / pi;
double pi2 = 2.*pi;
const double AU_to_m = 149597870700.;//1 AU [m]

int main(int argc, const char *argv[])
{
	printf("\nSatellite Tracker\n");


	std::vector <str_tle> tle;
	tle.resize(0);

	char TLE_filepath[MAX_PATH];
	sprintf_s(TLE_filepath, "%s", "");
	sprintf_s(TLE_filepath, "%s", "C:\\Predictions\\3le.txt");//this is TLE file path 

	read_TLE_catalogue(TLE_filepath, &tle);

	//time span to test
	double MJD_start_UTC = MJD_full(2021, 2, 12, 0, 0, 0);
	double MJD_stop_UTC = MJD_full(2021, 2, 13, 0, 0, 0);
	double search_step_s = 1.;
	int search_steps = int((MJD_stop_UTC - MJD_start_UTC)*86400. / search_step_s) + 1;

	//UAEUni
	v3 Station_TCS_m;//WGS84
	ae Station_TCS_rad;//East-North	
	Station_TCS_m.x = 3282240.3;
	Station_TCS_m.y = 4807257.6;
	Station_TCS_m.z = 2598574.1;
	Station_TCS_rad.az = 55.6760*deg2rad;
	Station_TCS_rad.el = 24.2006*deg2rad;

	FILE *FFout;//output file pointer
	errno_t err;
	char OutFilename[MAX_PATH];
	char txt_line[MAX_PATH];

	sprintf_s(OutFilename, "%s", "");//output file path
	sprintf_s(OutFilename, "%s_out.txt", TLE_filepath);
	err = fopen_s(&FFout, OutFilename, "w");

	//add header to output file
	sprintf_s(txt_line, "%s", "");
	sprintf_s(txt_line, "NORAD\tTime_MJD_UTC\tECI_X_m\tECI_Y_m\tECI_Z_m\tRA_deg\tDec_deg\tTopo_Az_deg\tTopo_El_deg\tSlant_range_m");
	fprintf(FFout, "%s\n", txt_line);
	fclose(FFout);

	bool sunlit_only = false;

	for (int iSat = 0; iSat < tle.size(); iSat++)//Satellite loop
	{
		elsetrec tle_satrec;
		int ReturnError = sgp4_satrec_init(&tle[iSat], &tle_satrec);
		//0 - no error
		//1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
		//2 - mean motion less than 0.0
		//3 - pert elements, ecc < 0.0  or  ecc > 1.0
		//4 - semi-latus rectum < 0.0
		//5 - epoch elements are sub-orbital
		//6 - satellite has decayed

		double pass_start_MJD = -1.;
		double pass_stop_MJD = -1.;

		err = fopen_s(&FFout, OutFilename, "a");//open output file to append new records

		for (int iStep = 0; iStep < search_steps; iStep++)
		{
			double MJD = MJD_start_UTC + double(iStep)*search_step_s / 86400.;
			bool step_accepted = false;

			double sat_ICS_km[3], sat_ICS_kmps[3];//Inertial system, ECI J2000
			int err = sgp4(tle_satrec, MJD, sat_ICS_km, sat_ICS_kmps);
			v3 sat_ICS_m;//Inertial system, ECI J2000
			sat_ICS_m.x = sat_ICS_km[0] * 1000.;
			sat_ICS_m.y = sat_ICS_km[1] * 1000.;
			sat_ICS_m.z = sat_ICS_km[2] * 1000.;
			ae sat_ICS_rad = get_ae(sat_ICS_m);

			double gast_rad = GAST_rad(MJD);//GAST
			v33 ICS_to_TCS = rz(gast_rad);//Inertial to Terrestrial
			v3 sat_TCS_m = ICS_to_TCS * sat_ICS_m;
			ae sat_topo_rad = topo_AzEl_rad(Station_TCS_m, Station_TCS_rad, sat_TCS_m);//topocentric coordinates			
			double slant_range_m = norm(sat_TCS_m - Station_TCS_m);

			if (sat_topo_rad.el > 0.)//satellite above horizon
			{
				step_accepted = true;

				if (sunlit_only)
				{
					//check for illumination			
					ae sun_ICS_rad = Sun_ICS_rad(MJD);
					double Sun_Earth_m = Distance_to_the_SUN_AU(MJD) * AU_to_m;
					v3 sun_ICS_m = get_v3(sun_ICS_rad, Sun_Earth_m);

					double ShadowFunction = shadow_function(MJD, sat_ICS_m, sun_ICS_m);

					if (ShadowFunction < 0.2)//satellite in shadow
						step_accepted = false;
				}
			}

			if (step_accepted)
			{
				if (pass_start_MJD < 0.)
					pass_start_MJD = MJD;

				pass_stop_MJD = MJD;

				//satellite trackable, add record				
				sprintf_s(txt_line, "%s", "");
				sprintf_s(txt_line, "%d\t%.6lf\t%.1lf\t%.1lf\t%.1lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.1lf",
					tle[iSat].NORAD,
					MJD,
					sat_ICS_m.x,
					sat_ICS_m.y,
					sat_ICS_m.z,
					sat_ICS_rad.az*rad2deg,
					sat_ICS_rad.el*rad2deg,
					sat_topo_rad.az*rad2deg,
					sat_topo_rad.el*rad2deg,
					slant_range_m);

				fprintf(FFout, "%s\n", txt_line);
			}

			if ((!step_accepted) || (iStep == (search_steps - 1)))
			{
				if ((pass_start_MJD > 0.) && (pass_stop_MJD > pass_start_MJD))
				{
					double pass_duration_s = (pass_stop_MJD - pass_start_MJD)*86400.;

					/*
					//pass finished, save info
					err = fopen_s(&FFout, OutFilename, "a");
					sprintf_s(txt_line, "%s", "");
					sprintf_s(txt_line, "%d\t%.6lf\t%.6lf\t%.1lf",
						tle[iSat].NORAD,
						pass_start_MJD,
						pass_stop_MJD,
						pass_duration_s);
					fprintf(FFout, "%s\n", txt_line);
					fclose(FFout);
					*/
				}
				pass_start_MJD = -1.;
				pass_stop_MJD = -1.;
			}

		}//iStep

		fclose(FFout);
	}


	tle.resize(0);
	return 0;
}

void read_TLE_catalogue(
	char *TLE_filepath,
	std::vector <str_tle> *tle
)
{
	(*tle).resize(0);

	if (fs::exists(TLE_filepath))
	{
		//printf("Reading TLE catalogue\n");

		FILE *FFin;
		errno_t err;
		char buf[MAXBUF];
		char Line[MAXBUF], L0[MAXBUF], L1[MAXBUF], L2[MAXBUF];//3LE

		sprintf_s(buf, "%s", "");
		sprintf_s(Line, "%s", "");
		sprintf_s(L0, "%s", "");
		sprintf_s(L1, "%s", "");
		sprintf_s(L2, "%s", "");

		err = fopen_s(&FFin, TLE_filepath, "r");

		int lines = 0;

		while (fgets(Line, sizeof(Line), FFin) != NULL)
		{
			// strip trailing '\n' if it exists
			int length = strlen(Line) - 1;
			if (Line[length] == '\n')
				Line[length] = 0;

			if (length > 1)
			{
				if (!memcmp(Line, "0 ", 2))
				{
					lines++;
					sprintf_s(L0, "%s", Line);
				}
				else if (!memcmp(Line, "1 ", 2))
				{
					lines++;
					sprintf_s(L1, "%s", Line);
				}
				else if (!memcmp(Line, "2 ", 2))
				{
					lines++;
					sprintf_s(L2, "%s", Line);

					int L1_length = strlen(L1);
					int L2_length = strlen(L2);

					if ((L1_length >= 68) && (L2_length >= 68))
					{
						str_tle tle_temp;
						read_TLE(L1, L2, &tle_temp);
						(*tle).push_back(tle_temp);
					}

					sprintf_s(buf, "%s", "");
					sprintf_s(Line, "%s", "");
					sprintf_s(L0, "%s", "");
					sprintf_s(L1, "%s", "");
					sprintf_s(L2, "%s", "");
				}
			}
		}
		fclose(FFin);
	}
}

void read_TLE(
	char *TLE_line1,//in
	char *TLE_line2,//in
	str_tle *tle//out
)
{
	/*
	Line 0
	01 Line Number of Element Data
	xx Satellite name

	Line 1
	Column Description
	01 Line Number of Element Data
	03-07 Satellite Number
	08 Classification (U=Unclassified)
	10-11 International Designator (Last two digits of launch year)
	12-14 International Designator (Launch number of the year)
	15-17 International Designator (Piece of the launch)
	19-20 Epoch Year (Last two digits of year)
	21-32 Epoch (Day of the year and fractional portion of the day)
	34-43 First Time Derivative of the Mean Motion
	45-52 Second Time Derivative of Mean Motion (decimal point assumed)
	54-61 BSTAR drag term (decimal point assumed)
	63 Ephemeris type
	65-68 Element number
	69 Checksum (Modulo 10)
	(Letters, blanks, periods, plus signs = 0; minus signs = 1)

	Line 2
	Column Description
	01 Line Number of Element Data
	03-07 Satellite Number
	09-16 Inclination [Degrees]
	18-25 Right Ascension of the Ascending Node [Degrees]
	27-33 Eccentricity (decimal point assumed)
	35-42 Argument of Perigee [Degrees]
	44-51 Mean Anomaly [Degrees]
	53-63 Mean Motion [Revs per day]
	64-68 Revolution number at epoch [Revs]
	69 Checksum (Modulo 10)
	*/

	//TLE parameters for SGP4
	(*tle).NORAD = 0;
	(*tle).MJD = 0.;
	(*tle).ndot = 0.;
	(*tle).nddot = 0.;
	(*tle).nexp = 0;
	(*tle).bstar = 0.;
	(*tle).ibexp = 0;
	(*tle).inclo = 0.;
	(*tle).nodeo = 0.;
	(*tle).ecco = 0.;
	(*tle).argpo = 0.;
	(*tle).mo = 0.;
	(*tle).no = 0.;

	//read TLE components
	char longstr1[130];
	char longstr2[130];
	sprintf_s(longstr1, "%s", "");
	sprintf_s(longstr2, "%s", "");
	sprintf_s(longstr1, "%s", TLE_line1);
	sprintf_s(longstr2, "%s", TLE_line2);

	int cardnumb, numb, j;
	long revnum = 0, elnum = 0;
	int year = 0;
	int mon, day, hr, minute, nexp, ibexp;
	int epochyr;
	long int satnum;
	double epochdays, ndot, nddot, bstar, inclo, nodeo, ecco, argpo, mo, no;

	int cospar_year = 0;
	char cospar_launch[4];
	char cospar_piece[4];
	sprintf_s(cospar_launch, "%s", "");
	sprintf_s(cospar_piece, "%s", "");
	sscanf_s(longstr1, "%*9c%2d%3c%3c", &cospar_year, cospar_launch, cospar_piece);
	cospar_launch[3] = '\0';
	cospar_piece[3] = '\0';
	for (int i = 0; i < 4; i++)
	{
		if (cospar_piece[i] == ' ')
			cospar_piece[i] = '\0';
	}

	if (cospar_year < 50)
		cospar_year += 2000;
	else
		cospar_year += 1900;

	sprintf_s((*tle).intldesg, "%s", "");
	sprintf_s((*tle).intldesg, "%4d-%s%s", cospar_year, cospar_launch, cospar_piece);

	// set the implied decimal points since doing a formated read
	// fixes for bad input data values (missing, ...)
	for (j = 10; j <= 15; j++)
		if (longstr1[j] == ' ')
			longstr1[j] = '_';

	if (longstr1[44] != ' ')
		longstr1[43] = longstr1[44];
	longstr1[44] = '.';
	if (longstr1[7] == ' ')
		longstr1[7] = 'U';
	if (longstr1[9] == ' ')
		longstr1[9] = '.';
	for (j = 45; j <= 49; j++)
		if (longstr1[j] == ' ')
			longstr1[j] = '0';
	if (longstr1[51] == ' ')
		longstr1[51] = '0';
	if (longstr1[53] != ' ')
		longstr1[52] = longstr1[53];
	longstr1[53] = '.';
	longstr2[25] = '.';
	for (j = 26; j <= 32; j++)
		if (longstr2[j] == ' ')
			longstr2[j] = '0';
	if (longstr1[62] == ' ')
		longstr1[62] = '0';
	if (longstr1[68] == ' ')
		longstr1[68] = '0';

	//longstr1[69] = '\0';

	sscanf_s(longstr1, "%2d %5ld %*11c %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld",
		&cardnumb, &satnum, &epochyr,
		&epochdays, &ndot, &nddot, &nexp, &bstar,
		&ibexp, &numb, &elnum);

	sscanf_s(longstr2, "%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld \n",
		&cardnumb, &satnum, &inclo,
		&nodeo, &ecco, &argpo, &mo, &no,
		&revnum);


	if (epochyr < 50)
		year = epochyr + 2000;
	else
		year = epochyr + 1900;

	(*tle).NORAD = satnum;
	(*tle).MJD = MJD_full_DoY1_d(year, epochdays);
	(*tle).ndot = ndot;
	(*tle).nddot = nddot;
	(*tle).nexp = nexp;
	(*tle).bstar = bstar;
	(*tle).ibexp = ibexp;
	(*tle).inclo = inclo;
	(*tle).nodeo = nodeo;
	(*tle).ecco = ecco;
	(*tle).argpo = argpo;
	(*tle).mo = mo;
	(*tle).no = no;

}

double MJD_full_DoY1_d(
	int year,//2006
	double doy1//day of year starting from 1, Jan 1=1doy
)
{
	int month, day, Hour, Minute;
	date_from_year_and_doy1(year, int(doy1), &month, &day);
	double SOD = (doy1 - double(int(doy1)))*86400.;
	Hour = int(SOD) / 3600;
	Minute = int(SOD) - Hour * 3600;
	double Second = SOD - double(Hour * 3600) - double(Minute * 60);
	double Full_MJD = MJD_full(year, month, day, Hour, Minute, Second);
	return Full_MJD;
}

void date_from_year_and_doy1(//doy1=day of year starting from 1, Jan 1=1doy
	int year,
	int doy1,
	int *month,
	int *day
)
{
	int d[12] = { 31,59,90,120,151,181,212,243,273,304,334,365 };
	//is leap year? XXI century
	if (year % 4 == 0)
	{
		for (unsigned int i = 1; i < 12; i++)
			d[i] = d[i] + 1;
	}

	*month = 1;
	for (unsigned int i = 0; i < 12; i++)
	{
		if (doy1 > d[i])
			*month = i + 2;
	}
	if (*month > 1)
		*day = doy1 - d[*month - 2];
	else
		*day = doy1;
}

long double MJD_full(
	int year,//2006
	int month,
	int day,
	int hour,
	int minute,
	double second
)
{
	long double result = 0.;
	static int k[] = { 0,31,59,90,120,151,181,212,243,273,304,334 };
	if (year % 4 == 0 && month > 2) ++day;
	result =
		double((k[--month] + day + (year - 1972) * 365 + (year - 1969) / 4) + (long)41316.) +
		(double(hour * 3600 + minute * 60) + second) / 86400.;
	return result;
}


int sgp4_satrec_init(
	str_tle *tle,//in
	elsetrec *satrec//out
)
{
	int ReturnError = 0;
	//0 - no error
	//1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
	//2 - mean motion less than 0.0
	//3 - pert elements, ecc < 0.0  or  ecc > 1.0
	//4 - semi-latus rectum < 0.0
	//5 - epoch elements are sub-orbital
	//6 - satellite has decayed

	(*satrec).NORAD = (*tle).NORAD;
	(*satrec).MJD = (*tle).MJD;
	(*satrec).ndot = (*tle).ndot;
	(*satrec).nddot = (*tle).nddot;
	int nexp = (*tle).nexp;
	(*satrec).bstar = (*tle).bstar;
	int ibexp = (*tle).ibexp;
	(*satrec).inclo = (*tle).inclo;
	(*satrec).nodeo = (*tle).nodeo;
	(*satrec).ecco = (*tle).ecco;
	(*satrec).argpo = (*tle).argpo;
	(*satrec).mo = (*tle).mo;
	(*satrec).no = (*tle).no;

	const double rad = 180.0 / pi;       //  57.29577951308230
	const double xpdotp = 1440.0 / pi2;  // 229.1831180523293
	double tumin, radiusearthkm, xke, j2, j3, j4, j3oj2;
	getgravconst(tumin, radiusearthkm, xke, j2, j3, j4, j3oj2);

	(*satrec).error = 0;

	// ---- find no, ndot, nddot ----
	(*satrec).no = (*satrec).no / xpdotp; //* rad/min
	(*satrec).nddot = (*satrec).nddot * std::pow(10.0, nexp);
	(*satrec).bstar = (*satrec).bstar * std::pow(10.0, ibexp);
	// ---- convert to sgp4 units ----
	(*satrec).a = pow((*satrec).no*tumin, (-2.0 / 3.0));
	(*satrec).ndot = (*satrec).ndot / (xpdotp*1440.0);  //* ? * minperday
	(*satrec).nddot = (*satrec).nddot / (xpdotp*1440.0*1440.0);
	// ---- find standard orbital elements ----
	(*satrec).inclo = (*satrec).inclo / rad;
	(*satrec).nodeo = (*satrec).nodeo / rad;
	(*satrec).argpo = (*satrec).argpo / rad;
	(*satrec).mo = (*satrec).mo / rad;
	(*satrec).alta = (*satrec).a*(1.0 + (*satrec).ecco*(*satrec).ecco) - 1.0;
	(*satrec).altp = (*satrec).a*(1.0 - (*satrec).ecco*(*satrec).ecco) - 1.0;

	// ---------------- initialize the orbit at sgp4epoch -------------------
	double JD = (*satrec).MJD + 2400000.5;
	sgp4init(JD - 2433281.5, (*satrec).bstar, (*satrec).ecco, (*satrec).argpo, (*satrec).inclo, (*satrec).mo, (*satrec).no, (*satrec).nodeo, (*satrec));

	return ReturnError;
}
//------------------------------------------------------------------------------

/* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  outputs       :
*    tumin       - minutes in one time unit
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*    mu          - earth gravitational parameter
*
*  coupling      :
*    none
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */

void getgravconst
(
	double& tumin,
	double& radiusearthkm,
	double& xke,
	double& j2,
	double& j3,
	double& j4,
	double& j3oj2
)
{
	// ------------ wgs-84 constants ------------
	double mu = 398600.5;            // in km3 / s2
	radiusearthkm = 6378.137;     // km
	double sqrt_argument = radiusearthkm * radiusearthkm*radiusearthkm / mu;
	if (sqrt_argument >= 0.)
		xke = 60.0 / sqrt(sqrt_argument);
	else
		xke = 0.;
	tumin = 1.0 / xke;
	j2 = 0.00108262998905;
	j3 = -0.00000253215306;
	j4 = -0.00000161098761;
	j3oj2 = j3 / j2;
}   // end getgravconst

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4init
*
*  this procedure initializes variables for sgp4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*
*  outputs       :
*    satrec      - common values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
*    cc1sq  , cc2    , cc3
*    coef   , coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    nodem       -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    getgravconst-
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*    sgp4        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

int sgp4init
(
	const double epoch,
	const double xbstar, const double xecco, const double xargpo,
	const double xinclo, const double xmo, const double xno,
	const double xnodeo, elsetrec& satrec
)
{
	/* --------------------- local variables ------------------------ */
	double ao, ainv, con42, cosio, sinio, cosio2, eccsq,
		omeosq, posq, rp, rteosq,
		cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq,
		cc2, cc3, coef, coef1, cosio4, day, dndt,
		em, emsq, eeta, etasq, gam, argpm, nodem,
		inclm, mm, nm, perige, pinvsq, psisq, qzms24,
		rtemsq, s1, s2, s3, s4, s5, s6,
		s7, sfour, ss1, ss2, ss3, ss4, ss5,
		ss6, ss7, sz1, sz2, sz3, sz11, sz12,
		sz13, sz21, sz22, sz23, sz31, sz32, sz33,
		tc, temp, temp1, temp2, temp3, tsi, xpidot,
		xhdot1, z1, z2, z3, z11, z12, z13,
		z21, z22, z23, z31, z32, z33,
		qzms2t, ss, j2, j3oj2, j4, x2o3, r[3], v[3],
		tumin, radiusearthkm, xke, j3;

	/* ------------------------ initialization --------------------- */
	// sgp4fix divisor for divide by zero check on inclination
	const double temp4 = 1.0 + cos(pi - 1.0e-9);

	/* ----------- set all near earth variables to zero ------------ */
	satrec.isimp = 0;   satrec.deep_space = false; satrec.aycof = 0.0;
	satrec.con41 = 0.0; satrec.cc1 = 0.0; satrec.cc4 = 0.0;
	satrec.cc5 = 0.0; satrec.d2 = 0.0; satrec.d3 = 0.0;
	satrec.d4 = 0.0; satrec.delmo = 0.0; satrec.eta = 0.0;
	satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao = 0.0;
	satrec.t2cof = 0.0; satrec.t3cof = 0.0;
	satrec.t4cof = 0.0; satrec.t5cof = 0.0; satrec.x1mth2 = 0.0;
	satrec.x7thm1 = 0.0; satrec.mdot = 0.0; satrec.nodedot = 0.0;
	satrec.xlcof = 0.0; satrec.xmcof = 0.0; satrec.nodecf = 0.0;

	/* ----------- set all deep space variables to zero ------------ */
	satrec.irez = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
	satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
	satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
	satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt = 0.0;
	satrec.del1 = 0.0; satrec.del2 = 0.0; satrec.del3 = 0.0;
	satrec.didt = 0.0; satrec.dmdt = 0.0; satrec.dnodt = 0.0;
	satrec.domdt = 0.0; satrec.e3 = 0.0; satrec.ee2 = 0.0;
	satrec.peo = 0.0; satrec.pgho = 0.0; satrec.pho = 0.0;
	satrec.pinco = 0.0; satrec.plo = 0.0; satrec.se2 = 0.0;
	satrec.se3 = 0.0; satrec.sgh2 = 0.0; satrec.sgh3 = 0.0;
	satrec.sgh4 = 0.0; satrec.sh2 = 0.0; satrec.sh3 = 0.0;
	satrec.si2 = 0.0; satrec.si3 = 0.0; satrec.sl2 = 0.0;
	satrec.sl3 = 0.0; satrec.sl4 = 0.0; satrec.gsto = 0.0;
	satrec.xfact = 0.0; satrec.xgh2 = 0.0; satrec.xgh3 = 0.0;
	satrec.xgh4 = 0.0; satrec.xh2 = 0.0; satrec.xh3 = 0.0;
	satrec.xi2 = 0.0; satrec.xi3 = 0.0; satrec.xl2 = 0.0;
	satrec.xl3 = 0.0; satrec.xl4 = 0.0; satrec.xlamo = 0.0;
	satrec.zmol = 0.0; satrec.zmos = 0.0; satrec.atime = 0.0;
	satrec.xli = 0.0; satrec.xni = 0.0;

	// sgp4fix - note the following variables are also passed directly via satrec.
	// it is possible to streamline the sgp4init call by deleting the "x"
	// variables, but the user would need to set the satrec.* values first. we
	// include the additional assignments in case twoline2rv is not used.
	satrec.bstar = xbstar;
	satrec.ecco = xecco;
	satrec.argpo = xargpo;
	satrec.inclo = xinclo;
	satrec.mo = xmo;
	satrec.no = xno;
	satrec.nodeo = xnodeo;

	double time_since = 0.;

	/* ------------------------ earth constants ----------------------- */
	// sgp4fix identify constants and allow alternate values
	getgravconst(tumin, radiusearthkm, xke, j2, j3, j4, j3oj2);
	ss = 78.0 / radiusearthkm + 1.0;
	qzms2t = std::pow(((120.0 - 78.0) / radiusearthkm), 4);
	x2o3 = 2.0 / 3.0;

	satrec.init = true;

	initl
	(
		satrec.ecco, epoch, satrec.inclo, satrec.no, satrec.deep_space,
		ainv, ao, satrec.con41, con42, cosio, cosio2, eccsq, omeosq,
		posq, rp, rteosq, sinio, satrec.gsto
	);
	satrec.error = 0;

	if (rp < 1.0)
	{
		//         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
		satrec.error = 5;
	}

	if ((omeosq >= 0.0) || (satrec.no >= 0.0))
	{
		satrec.isimp = 0;
		if (rp < (220.0 / radiusearthkm + 1.0))
			satrec.isimp = 1;
		sfour = ss;
		qzms24 = qzms2t;
		perige = (rp - 1.0) * radiusearthkm;

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if (perige < 156.0)
		{
			sfour = perige - 78.0;
			if (perige < 98.0)
				sfour = 20.0;
			qzms24 = pow(((120.0 - sfour) / radiusearthkm), 4.0);
			sfour = sfour / radiusearthkm + 1.0;
		}
		pinvsq = 1.0 / posq;

		tsi = 1.0 / (ao - sfour);
		satrec.eta = ao * satrec.ecco * tsi;
		etasq = satrec.eta * satrec.eta;
		eeta = satrec.ecco * satrec.eta;
		psisq = fabs(1.0 - etasq);
		coef = qzms24 * pow(tsi, 4.0);
		coef1 = coef / pow(psisq, 3.5);
		cc2 = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *
			(4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satrec.con41 *
			(8.0 + 3.0 * etasq * (8.0 + etasq)));
		satrec.cc1 = satrec.bstar * cc2;
		cc3 = 0.0;
		if (satrec.ecco > 1.0e-4)
			cc3 = -2.0 * coef * tsi * j3oj2 * satrec.no * sinio / satrec.ecco;
		satrec.x1mth2 = 1.0 - cosio2;
		satrec.cc4 = 2.0* satrec.no * coef1 * ao * omeosq *
			(satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *
			(0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
				(-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
				(1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
					(2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
		satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
			(etasq + eeta) + eeta * etasq);
		cosio4 = cosio2 * cosio2;
		temp1 = 1.5 * j2 * pinvsq * satrec.no;
		temp2 = 0.5 * temp1 * j2 * pinvsq;
		temp3 = -0.46875 * j4 * pinvsq * pinvsq * satrec.no;
		satrec.mdot = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 *
			temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
		satrec.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 *
			(7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
			temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
		xhdot1 = -temp1 * cosio;
		satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
			2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
		xpidot = satrec.argpdot + satrec.nodedot;
		satrec.omgcof = satrec.bstar * cc3 * cos(satrec.argpo);
		satrec.xmcof = 0.0;
		if (satrec.ecco > 1.0e-4)
			satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
		satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
		satrec.t2cof = 1.5 * satrec.cc1;
		// sgp4fix for divide by zero with xinco = 180 deg
		if (fabs(cosio + 1.0) > 1.5e-12)
			satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
		else
			satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
		satrec.aycof = -0.5 * j3oj2 * sinio;
		satrec.delmo = std::pow((1.0 + satrec.eta * cos(satrec.mo)), 3);
		satrec.sinmao = sin(satrec.mo);
		satrec.x7thm1 = 7.0 * cosio2 - 1.0;

		/* --------------- deep space initialization ------------- */
		if ((pi2 / satrec.no) >= 225.0)
		{
			satrec.deep_space = true;
			satrec.isimp = 1;
			tc = 0.0;
			inclm = satrec.inclo;

			dscom
			(
				epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo,
				satrec.no, snodm, cnodm, sinim, cosim, sinomm, cosomm,
				day, satrec.e3, satrec.ee2, em, emsq, gam,
				satrec.peo, satrec.pgho, satrec.pho, satrec.pinco,
				satrec.plo, rtemsq, satrec.se2, satrec.se3,
				satrec.sgh2, satrec.sgh3, satrec.sgh4,
				satrec.sh2, satrec.sh3, satrec.si2, satrec.si3,
				satrec.sl2, satrec.sl3, satrec.sl4, s1, s2, s3, s4, s5,
				s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3,
				sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33,
				satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
				satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
				satrec.xl3, satrec.xl4, nm, z1, z2, z3, z11,
				z12, z13, z21, z22, z23, z31, z32, z33,
				satrec.zmol, satrec.zmos
			);
			dpper
			(
				satrec.e3, satrec.ee2, satrec.peo, satrec.pgho,
				satrec.pho, satrec.pinco, satrec.plo, satrec.se2,
				satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
				satrec.sh2, satrec.sh3, satrec.si2, satrec.si3,
				satrec.sl2, satrec.sl3, satrec.sl4, time_since,
				satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
				satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
				satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos, inclm, satrec.init,
				satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo
			);

			argpm = 0.0;
			nodem = 0.0;
			mm = 0.0;

			dsinit
			(
				cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
				ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, time_since, tc,
				satrec.gsto, satrec.mo, satrec.mdot, satrec.no, satrec.nodeo,
				satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
				satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
				satrec.irez, satrec.atime,
				satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222,
				satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232,
				satrec.d5421, satrec.d5433, satrec.dedt, satrec.didt,
				satrec.dmdt, dndt, satrec.dnodt, satrec.domdt,
				satrec.del1, satrec.del2, satrec.del3, satrec.xfact,
				satrec.xlamo, satrec.xli, satrec.xni
			);
		}

		/* ----------- set variables if not deep space ----------- */
		if (satrec.isimp != 1)
		{
			cc1sq = satrec.cc1 * satrec.cc1;
			satrec.d2 = 4.0 * ao * tsi * cc1sq;
			temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
			satrec.d3 = (17.0 * ao + sfour) * temp;
			satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
				satrec.cc1;
			satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
			satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
				(12.0 * satrec.d2 + 10.0 * cc1sq));
			satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
				12.0 * satrec.cc1 * satrec.d3 +
				6.0 * satrec.d2 * satrec.d2 +
				15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
		}
	} // if omeosq = 0 ...

	/* finally propogate to zero epoch to initialise all others. */
	if (satrec.error == 0)
		sgp4(satrec, 0.0, r, v);

	satrec.init = false;

	//#include "debug6.cpp"
	return satrec.error;
}  // end sgp4init


/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*    satn        - satellite number
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

  //static
void initl
(
	double ecco, double epoch, double inclo, double& no,
	bool &deep_space,
	double& ainv, double& ao, double& con41, double& con42, double& cosio,
	double& cosio2, double& eccsq, double& omeosq, double& posq,
	double& rp, double& rteosq, double& sinio, double& gsto
)
{
	/* --------------------- local variables ------------------------ */
	double ak, d1, del, adel, po, x2o3, j2, xke,
		tumin, radiusearthkm, j3, j4, j3oj2;

	/* ----------------------- earth constants ---------------------- */
	// sgp4fix identify constants and allow alternate values
	getgravconst(tumin, radiusearthkm, xke, j2, j3, j4, j3oj2);
	x2o3 = 2.0 / 3.0;

	/* ------------- calculate auxillary epoch quantities ---------- */
	eccsq = ecco * ecco;
	omeosq = 1.0 - eccsq;

	if (omeosq >= 0.)
		rteosq = sqrt(omeosq);
	else
		rteosq = 0.;

	cosio = cos(inclo);
	cosio2 = cosio * cosio;

	/* ------------------ un-kozai the mean motion ----------------- */
	ak = pow(xke / no, x2o3);
	d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
	del = d1 / (ak * ak);
	adel = ak * (1.0 - del * del - del *
		(1.0 / 3.0 + 134.0 * del * del / 81.0));
	del = d1 / (adel * adel);
	no = no / (1.0 + del);

	ao = pow(xke / no, x2o3);
	sinio = sin(inclo);
	po = ao * omeosq;
	con42 = 1.0 - 5.0 * cosio2;
	con41 = -con42 - cosio2 - cosio2;
	ainv = 1.0 / ao;
	posq = po * po;
	rp = ao * (1.0 - ecco);
	deep_space = false;

	gsto = gstime(epoch + 2433281.5);

	//#include "debug5.cpp"
}  // end initl


/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
  //static
void dscom(
	double epoch, double ep, double argpp, double tc, double inclp,
	double nodep, double np,
	double& snodm, double& cnodm, double& sinim, double& cosim, double& sinomm,
	double& cosomm, double& day, double& e3, double& ee2, double& em,
	double& emsq, double& gam, double& peo, double& pgho, double& pho,
	double& pinco, double& plo, double& rtemsq, double& se2, double& se3,
	double& sgh2, double& sgh3, double& sgh4, double& sh2, double& sh3,
	double& si2, double& si3, double& sl2, double& sl3, double& sl4,
	double& s1, double& s2, double& s3, double& s4, double& s5,
	double& s6, double& s7, double& ss1, double& ss2, double& ss3,
	double& ss4, double& ss5, double& ss6, double& ss7, double& sz1,
	double& sz2, double& sz3, double& sz11, double& sz12, double& sz13,
	double& sz21, double& sz22, double& sz23, double& sz31, double& sz32,
	double& sz33, double& xgh2, double& xgh3, double& xgh4, double& xh2,
	double& xh3, double& xi2, double& xi3, double& xl2, double& xl3,
	double& xl4, double& nm, double& z1, double& z2, double& z3,
	double& z11, double& z12, double& z13, double& z21, double& z22,
	double& z23, double& z31, double& z32, double& z33, double& zmol,
	double& zmos
)
{
	//constants
	const double zes = 0.01675;
	const double zel = 0.05490;
	const double c1ss = 2.9864797e-6;
	const double c1l = 4.7968065e-7;
	const double zsinis = 0.39785416;
	const double zcosis = 0.91744867;
	const double zcosgs = 0.1945905;
	const double zsings = -0.98088458;
	const double twopi = pi2;
	//local variables
	int lsflg;
	double a1, a2, a3, a4, a5, a6, a7,
		a8, a9, a10, betasq, cc, ctem, stem,
		x1, x2, x3, x4, x5, x6, x7,
		x8, xnodce, xnoi, zcosg, zcosgl, zcosh, zcoshl,
		zcosi, zcosil, zsing, zsingl, zsinh, zsinhl, zsini,
		zsinil, zx, zy;
	nm = np;
	em = ep;
	snodm = sin(nodep);
	cnodm = cos(nodep);
	sinomm = sin(argpp);
	cosomm = cos(argpp);
	sinim = sin(inclp);
	cosim = cos(inclp);
	emsq = em * em;
	betasq = 1.0 - emsq;

	if (betasq >= 0.)
		rtemsq = sqrt(betasq);
	else
		rtemsq = 0.;

	//initialize lunar solar terms
	peo = 0.0;
	pinco = 0.0;
	plo = 0.0;
	pgho = 0.0;
	pho = 0.0;
	day = epoch + 18261.5 + tc / 1440.0;
	xnodce = fmod(4.5236020 - 9.2422029e-4 * day, twopi);
	stem = sin(xnodce);
	ctem = cos(xnodce);
	zcosil = 0.91375164 - 0.03568096 * ctem;

	double sqrt_argument = 1.0 - zcosil * zcosil;
	if (sqrt_argument >= 0.)
		zsinil = sqrt(sqrt_argument);
	else
		zsinil = 0.;

	zsinhl = 0.089683511 * stem / zsinil;


	sqrt_argument = 1.0 - zsinhl * zsinhl;
	if (sqrt_argument >= 0.)
		zcoshl = sqrt(sqrt_argument);
	else
		zcoshl = 0.;

	gam = 5.8351514 + 0.0019443680 * day;
	zx = 0.39785416 * stem / zsinil;
	zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
	zx = atan2_check(zx, zy);
	zx = gam + zx - xnodce;
	zcosgl = cos(zx);
	zsingl = sin(zx);
	//do solar terms
	zcosg = zcosgs;
	zsing = zsings;
	zcosi = zcosis;
	zsini = zsinis;
	zcosh = cnodm;
	zsinh = snodm;
	cc = c1ss;
	xnoi = 1.0 / nm;
	for (lsflg = 1; lsflg <= 2; lsflg++)
	{
		a1 = zcosg * zcosh + zsing * zcosi * zsinh;
		a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
		a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
		a8 = zsing * zsini;
		a9 = zsing * zsinh + zcosg * zcosi * zcosh;
		a10 = zcosg * zsini;
		a2 = cosim * a7 + sinim * a8;
		a4 = cosim * a9 + sinim * a10;
		a5 = -sinim * a7 + cosim * a8;
		a6 = -sinim * a9 + cosim * a10;
		x1 = a1 * cosomm + a2 * sinomm;
		x2 = a3 * cosomm + a4 * sinomm;
		x3 = -a1 * sinomm + a2 * cosomm;
		x4 = -a3 * sinomm + a4 * cosomm;
		x5 = a5 * sinomm;
		x6 = a6 * sinomm;
		x7 = a5 * cosomm;
		x8 = a6 * cosomm;
		z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
		z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
		z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
		z1 = 3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
		z2 = 6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
		z3 = 3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
		z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
		z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
		z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
		z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
		z22 = 6.0 *  (a4 * a5 + a2 * a6) + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
		z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
		z1 = z1 + z1 + betasq * z31;
		z2 = z2 + z2 + betasq * z32;
		z3 = z3 + z3 + betasq * z33;
		s3 = cc * xnoi;
		s2 = -0.5 * s3 / rtemsq;
		s4 = s3 * rtemsq;
		s1 = -15.0 * em * s4;
		s5 = x1 * x3 + x2 * x4;
		s6 = x2 * x3 + x1 * x4;
		s7 = x2 * x4 - x1 * x3;
		//do lunar terms
		if (lsflg == 1)
		{
			ss1 = s1;
			ss2 = s2;
			ss3 = s3;
			ss4 = s4;
			ss5 = s5;
			ss6 = s6;
			ss7 = s7;
			sz1 = z1;
			sz2 = z2;
			sz3 = z3;
			sz11 = z11;
			sz12 = z12;
			sz13 = z13;
			sz21 = z21;
			sz22 = z22;
			sz23 = z23;
			sz31 = z31;
			sz32 = z32;
			sz33 = z33;
			zcosg = zcosgl;
			zsing = zsingl;
			zcosi = zcosil;
			zsini = zsinil;
			zcosh = zcoshl * cnodm + zsinhl * snodm;
			zsinh = snodm * zcoshl - cnodm * zsinhl;
			cc = c1l;
		}
	}

	zmol = fmod(4.7199672 + 0.22997150  * day - gam, twopi);
	zmos = fmod(6.2565837 + 0.017201977 * day, twopi);

	//do solar terms
	se2 = 2.0 * ss1 * ss6;
	se3 = 2.0 * ss1 * ss7;
	si2 = 2.0 * ss2 * sz12;
	si3 = 2.0 * ss2 * (sz13 - sz11);
	sl2 = -2.0 * ss3 * sz2;
	sl3 = -2.0 * ss3 * (sz3 - sz1);
	sl4 = -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
	sgh2 = 2.0 * ss4 * sz32;
	sgh3 = 2.0 * ss4 * (sz33 - sz31);
	sgh4 = -18.0 * ss4 * zes;
	sh2 = -2.0 * ss2 * sz22;
	sh3 = -2.0 * ss2 * (sz23 - sz21);

	//do lunar terms
	ee2 = 2.0 * s1 * s6;
	e3 = 2.0 * s1 * s7;
	xi2 = 2.0 * s2 * z12;
	xi3 = 2.0 * s2 * (z13 - z11);
	xl2 = -2.0 * s3 * z2;
	xl3 = -2.0 * s3 * (z3 - z1);
	xl4 = -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
	xgh2 = 2.0 * s4 * z32;
	xgh3 = 2.0 * s4 * (z33 - z31);
	xgh4 = -18.0 * s4 * zel;
	xh2 = -2.0 * s2 * z22;
	xh3 = -2.0 * s2 * (z23 - z21);
}  // end dscom

/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
  //static
void dpper(
	double e3, double ee2, double peo, double pgho, double pho,
	double pinco, double plo, double se2, double se3, double sgh2,
	double sgh3, double sgh4, double sh2, double sh3, double si2,
	double si3, double sl2, double sl3, double sl4, double t,
	double xgh2, double xgh3, double xgh4, double xh2, double xh3,
	double xi2, double xi3, double xl2, double xl3, double xl4,
	double zmol, double zmos, double inclo,
	bool init,
	double& ep, double& inclp, double& nodep, double& argpp, double& mp
)
{
	//local variables
	const double twopi = pi2;
	// char ildm;
	double alfdp, betdp, cosip, cosop, dalf, dbet, dls,
		f2, f3, pe, pgh, ph, pinc, pl,
		sel, ses, sghl, sghs, shll, shs, sil,
		sinip, sinop, sinzf, sis, sll, sls, xls,
		xnoh, zf, zm, zel, zes, znl, zns;
	//constants
	zns = 1.19459e-5;
	zes = 0.01675;
	znl = 1.5835218e-4;
	zel = 0.05490;
	//calculate time varying periodics
	zm = zmos + zns * t;
	//be sure that the initial call has time set to zero
	if (init)
		zm = zmos;

	zf = zm + 2.0 * zes * sin(zm);
	sinzf = sin(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * cos(zf);
	ses = se2 * f2 + se3 * f3;
	sis = si2 * f2 + si3 * f3;
	sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
	sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
	shs = sh2 * f2 + sh3 * f3;
	zm = zmol + znl * t;
	if (init)
		zm = zmol;

	zf = zm + 2.0 * zel * sin(zm);
	sinzf = sin(zf);
	f2 = 0.5 * sinzf * sinzf - 0.25;
	f3 = -0.5 * sinzf * cos(zf);
	sel = ee2 * f2 + e3 * f3;
	sil = xi2 * f2 + xi3 * f3;
	sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
	sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
	shll = xh2 * f2 + xh3 * f3;
	pe = ses + sel;
	pinc = sis + sil;
	pl = sls + sll;
	pgh = sghs + sghl;
	ph = shs + shll;
	if (!init)
	{
		//0.2 rad = 11.45916 deg
		//sgp4fix for lyddane choice
		//add next three lines to set up use of original inclination per strn3 ver
	  //  ildm = 'y';
	  //  if (inclo >= 0.2) ildm = 'n';
		pe = pe - peo;
		pinc = pinc - pinco;
		pl = pl - plo;
		pgh = pgh - pgho;
		ph = ph - pho;
		inclp = inclp + pinc;
		ep = ep + pe;
		sinip = sin(inclp);
		cosip = cos(inclp);
		//apply periodics directly
		//sgp4fix for lyddane choice
		//strn3 used original inclination - this is technically feasible
		//gsfc used perturbed inclination - also technically feasible
		//probably best to readjust the 0.2 limit value and limit discontinuity
		//use next line for original strn3 approach and original inclination
		//if (inclo >= 0.2)
		//use next line for gsfc version and perturbed inclination
		if (inclp >= 0.2)
		{
			ph = ph / sinip;
			pgh = pgh - cosip * ph;
			argpp = argpp + pgh;
			nodep = nodep + ph;
			mp = mp + pl;
		}
		else
		{
			//apply periodics with lyddane modification
			sinop = sin(nodep);
			cosop = cos(nodep);
			alfdp = sinip * sinop;
			betdp = sinip * cosop;
			dalf = ph * cosop + pinc * cosip * sinop;
			dbet = -ph * sinop + pinc * cosip * cosop;
			alfdp = alfdp + dalf;
			betdp = betdp + dbet;
			nodep = fmod(nodep, twopi);
			xls = mp + argpp + cosip * nodep;
			dls = pl + pgh - pinc * nodep * sinip;
			xls = xls + dls;
			xnoh = nodep;
			nodep = atan2_check(alfdp, betdp);
			if (fabs(xnoh - nodep) > pi)
				if (nodep < xnoh)
					nodep = nodep + twopi;
				else
					nodep = nodep - twopi;
			mp = mp + pl;
			argpp = xls - mp - cosip * nodep;
		}
	}
}  // end dpper


/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
  //static
void dsinit(
	double cosim, double emsq, double argpo, double s1, double s2,
	double s3, double s4, double s5, double sinim, double ss1,
	double ss2, double ss3, double ss4, double ss5, double sz1,
	double sz3, double sz11, double sz13, double sz21, double sz23,
	double sz31, double sz33, double t, double tc, double gsto,
	double mo, double mdot, double no, double nodeo, double nodedot,
	double xpidot, double z1, double z3, double z11, double z13,
	double z21, double z23, double z31, double z33, double ecco,
	double eccsq, double& em, double& argpm, double& inclm, double& mm,
	double& nm, double& nodem,
	int& irez,
	double& atime, double& d2201, double& d2211, double& d3210, double& d3222,
	double& d4410, double& d4422, double& d5220, double& d5232, double& d5421,
	double& d5433, double& dedt, double& didt, double& dmdt, double& dndt,
	double& dnodt, double& domdt, double& del1, double& del2, double& del3,
	double& xfact, double& xlamo, double& xli, double& xni
)
{
	/* --------------------- local variables ------------------------ */
	const double twopi = pi2;

	double ainv2, aonv = 0.0, cosisq, eoc, f220, f221, f311,
		f321, f322, f330, f441, f442, f522, f523,
		f542, f543, g200, g201, g211, g300, g310,
		g322, g410, g422, g520, g521, g532, g533,
		ses, sgs, sghl, sghs, shs, shll, sis,
		sini2, sls, temp, temp1, theta, xno2, q22,
		q31, q33, root22, root44, root54, rptim, root32,
		root52, x2o3, xke, znl, emo, zns, emsqo,
		tumin, radiusearthkm, j2, j3, j4, j3oj2;

	q22 = 1.7891679e-6;
	q31 = 2.1460748e-6;
	q33 = 2.2123015e-7;
	root22 = 1.7891679e-6;
	root44 = 7.3636953e-9;
	root54 = 2.1765803e-9;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	root32 = 3.7393792e-7;
	root52 = 1.1428639e-7;
	x2o3 = 2.0 / 3.0;
	znl = 1.5835218e-4;
	zns = 1.19459e-5;

	// sgp4fix identify constants and allow alternate values
	getgravconst(tumin, radiusearthkm, xke, j2, j3, j4, j3oj2);

	/* -------------------- deep space initialization ------------ */
	irez = 0;
	if ((nm < 0.0052359877) && (nm > 0.0034906585))
		irez = 1;
	if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
		irez = 2;

	/* ------------------------ do solar terms ------------------- */
	ses = ss1 * zns * ss5;
	sis = ss2 * zns * (sz11 + sz13);
	sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
	sghs = ss4 * zns * (sz31 + sz33 - 6.0);
	shs = -zns * ss2 * (sz21 + sz23);
	// sgp4fix for 180 deg incl
	if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
		shs = 0.0;
	if (sinim != 0.0)
		shs = shs / sinim;
	sgs = sghs - cosim * shs;

	/* ------------------------- do lunar terms ------------------ */
	dedt = ses + s1 * znl * s5;
	didt = sis + s2 * znl * (z11 + z13);
	dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
	sghl = s4 * znl * (z31 + z33 - 6.0);
	shll = -znl * s2 * (z21 + z23);
	// sgp4fix for 180 deg incl
	if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
		shll = 0.0;
	domdt = sgs + sghl;
	dnodt = shs;
	if (sinim != 0.0)
	{
		domdt = domdt - cosim / sinim * shll;
		dnodt = dnodt + shll / sinim;
	}

	/* ----------- calculate deep space resonance effects -------- */
	dndt = 0.0;
	theta = fmod(gsto + tc * rptim, twopi);
	em = em + dedt * t;
	inclm = inclm + didt * t;
	argpm = argpm + domdt * t;
	nodem = nodem + dnodt * t;
	mm = mm + dmdt * t;
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//if (inclm < 0.0)
	//  {
	//    inclm  = -inclm;
	//    argpm  = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if (irez != 0)
	{
		aonv = pow(nm / xke, x2o3);

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if (irez == 2)
		{
			cosisq = cosim * cosim;
			emo = em;
			em = ecco;
			emsqo = emsq;
			emsq = eccsq;
			eoc = em * emsq;
			g201 = -0.306 - (em - 0.64) * 0.440;

			if (em <= 0.65)
			{
				g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
				g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
				g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
				g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
				g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
				g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
			}
			else
			{
				g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
				g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
				g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
				g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
				g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
				if (em > 0.715)
					g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
				else
					g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
			}
			if (em < 0.7)
			{
				g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
				g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
				g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;
			}
			else
			{
				g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
				g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
				g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
			}

			sini2 = sinim * sinim;
			f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
			f221 = 1.5 * sini2;
			f321 = 1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
			f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
			f441 = 35.0 * sini2 * f220;
			f442 = 39.3750 * sini2 * sini2;
			f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
				0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
			f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
				10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
			f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq *
				(-12.0 + 8.0 * cosim + 10.0 * cosisq));
			f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq *
				(12.0 + 8.0 * cosim - 10.0 * cosisq));
			xno2 = nm * nm;
			ainv2 = aonv * aonv;
			temp1 = 3.0 * xno2 * ainv2;
			temp = temp1 * root22;
			d2201 = temp * f220 * g201;
			d2211 = temp * f221 * g211;
			temp1 = temp1 * aonv;
			temp = temp1 * root32;
			d3210 = temp * f321 * g310;
			d3222 = temp * f322 * g322;
			temp1 = temp1 * aonv;
			temp = 2.0 * temp1 * root44;
			d4410 = temp * f441 * g410;
			d4422 = temp * f442 * g422;
			temp1 = temp1 * aonv;
			temp = temp1 * root52;
			d5220 = temp * f522 * g520;
			d5232 = temp * f523 * g532;
			temp = 2.0 * temp1 * root54;
			d5421 = temp * f542 * g521;
			d5433 = temp * f543 * g533;
			xlamo = fmod(mo + nodeo + nodeo - theta - theta, twopi);
			xfact = mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
			em = emo;
			emsq = emsqo;
		}

		/* ---------------- synchronous resonance terms -------------- */
		if (irez == 1)
		{
			g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
			g310 = 1.0 + 2.0 * emsq;
			g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
			f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
			f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
			f330 = 1.0 + cosim;
			f330 = 1.875 * f330 * f330 * f330;
			del1 = 3.0 * nm * nm * aonv * aonv;
			del2 = 2.0 * del1 * f220 * g200 * q22;
			del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
			del1 = del1 * f311 * g310 * q31 * aonv;
			xlamo = fmod(mo + nodeo + argpo - theta, twopi);
			xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		xli = xlamo;
		xni = no;
		atime = 0.0;
		nm = no + dndt;
	}

	//#include "debug3.cpp"
}  // end dsinit



/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the nasa
*    release on the internet. there are a few fixes that are added to correct
*    known errors in the existing implementations.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time since epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    getgravconst-
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

int sgp4(elsetrec& satrec, double MJD, double r[3], double v[3])
{
	double am, axnl, aynl, betal, cosim, cnod,
		cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu,
		delm, delomg, em, emsq, ecose, el2, eo1,
		ep, esine, argpm, argpp, argpdf, pl, mrt = 0.0,
		mvt, rdotl, rl, rvdot, rvdotl, sinim,
		sin2u, sineo1, sini, sinip, sinsu, sinu,
		snod, su, t2, t3, t4, tem5, temp,
		temp1, temp2, tempa, tempe, templ, u, ux,
		uy, uz, vx, vy, vz, inclm, mm,
		nm, nodem, xinc, xincp, xl, xlm, mp,
		xmdf, xmx, xmy, nodedf, xnode, nodep, tc, dndt,
		twopi, x2o3, j2, j3, tumin, j4, xke, j3oj2, radiusearthkm,
		vkmpersec;
	int ktr;

	/* ------------------ set mathematical constants --------------- */
	// sgp4fix divisor for divide by zero check on inclination
	const double temp4 = 1.0 + cos(pi - 1.0e-9);
	twopi = pi2;
	x2o3 = 2.0 / 3.0;
	// sgp4fix identify constants and allow alternate values
	getgravconst(tumin, radiusearthkm, xke, j2, j3, j4, j3oj2);
	vkmpersec = radiusearthkm * xke / 60.0;

	/* --------------------- clear sgp4 error flag ----------------- */

	double time_since_minute = (MJD - satrec.MJD)*1440.;
	satrec.error = 0;

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf = satrec.mo + satrec.mdot * time_since_minute;
	argpdf = satrec.argpo + satrec.argpdot * time_since_minute;
	nodedf = satrec.nodeo + satrec.nodedot * time_since_minute;
	argpm = argpdf;
	mm = xmdf;
	t2 = time_since_minute * time_since_minute;
	nodem = nodedf + satrec.nodecf * t2;
	tempa = 1.0 - satrec.cc1 * time_since_minute;
	tempe = satrec.bstar * satrec.cc4 * time_since_minute;
	templ = satrec.t2cof * t2;

	if (satrec.isimp != 1)
	{
		delomg = satrec.omgcof * time_since_minute;
		delm = satrec.xmcof *
			(std::pow((1.0 + satrec.eta * cos(xmdf)), 3) -
				satrec.delmo);
		temp = delomg + delm;
		mm = xmdf + temp;
		argpm = argpdf - temp;
		t3 = t2 * time_since_minute;
		t4 = t3 * time_since_minute;
		tempa = tempa - satrec.d2 * t2 - satrec.d3 * t3 -
			satrec.d4 * t4;
		tempe = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -
			satrec.sinmao);
		templ = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
			time_since_minute * satrec.t5cof);
	}

	nm = satrec.no;
	em = satrec.ecco;
	inclm = satrec.inclo;
	if (satrec.deep_space) //deep space?
	{
		tc = time_since_minute;
		dspace
		(
			satrec.irez,
			satrec.d2201, satrec.d2211, satrec.d3210,
			satrec.d3222, satrec.d4410, satrec.d4422,
			satrec.d5220, satrec.d5232, satrec.d5421,
			satrec.d5433, satrec.dedt, satrec.del1,
			satrec.del2, satrec.del3, satrec.didt,
			satrec.dmdt, satrec.dnodt, satrec.domdt,
			satrec.argpo, satrec.argpdot, time_since_minute, tc,
			satrec.gsto, satrec.xfact, satrec.xlamo,
			satrec.no, satrec.atime,
			em, argpm, inclm, satrec.xli, mm, satrec.xni,
			nodem, dndt, nm
		);
	} // if method = d

	if (nm <= 0.0)
	{
		//         printf("# error nm %f\n", nm);
		satrec.error = 2;
	}
	am = pow((xke / nm), x2o3) * tempa * tempa;
	nm = xke / pow(am, 1.5);
	em = em - tempe;

	// fix tolerance for error recognition
	if ((em >= 1.0) || (em < -0.001) || (am < 0.95))
	{
		//         printf("# error em %f\n", em);
		satrec.error = 1;
	}
	if (em < 0.0)
		em = 1.0e-6;
	mm = mm + satrec.no * templ;
	xlm = mm + argpm + nodem;
	emsq = em * em;
	temp = 1.0 - emsq;

	nodem = fmod(nodem, twopi);
	argpm = fmod(argpm, twopi);
	xlm = fmod(xlm, twopi);
	mm = fmod(xlm - argpm - nodem, twopi);

	/* ----------------- compute extra mean quantities ------------- */
	sinim = sin(inclm);
	cosim = cos(inclm);

	/* -------------------- add lunar-solar periodics -------------- */
	ep = em;
	xincp = inclm;
	argpp = argpm;
	nodep = nodem;
	mp = mm;
	sinip = sinim;
	cosip = cosim;
	if (satrec.deep_space)
	{
		dpper
		(
			satrec.e3, satrec.ee2, satrec.peo,
			satrec.pgho, satrec.pho, satrec.pinco,
			satrec.plo, satrec.se2, satrec.se3,
			satrec.sgh2, satrec.sgh3, satrec.sgh4,
			satrec.sh2, satrec.sh3, satrec.si2,
			satrec.si3, satrec.sl2, satrec.sl3,
			satrec.sl4, time_since_minute, satrec.xgh2,
			satrec.xgh3, satrec.xgh4, satrec.xh2,
			satrec.xh3, satrec.xi2, satrec.xi3,
			satrec.xl2, satrec.xl3, satrec.xl4,
			satrec.zmol, satrec.zmos, satrec.inclo,
			'n', ep, xincp, nodep, argpp, mp
		);
		if (xincp < 0.0)
		{
			xincp = -xincp;
			nodep = nodep + pi;
			argpp = argpp - pi;
		}
		if ((ep < 0.0) || (ep > 1.0))
		{
			//            printf("# error ep %f\n", ep);
			satrec.error = 3;
		}
	} // if method = d

  /* -------------------- long period periodics ------------------ */
	if (satrec.deep_space)
	{
		sinip = sin(xincp);
		cosip = cos(xincp);
		satrec.aycof = -0.5*j3oj2*sinip;
		// sgp4fix for divide by zero for xincp = 180 deg
		if (fabs(cosip + 1.0) > 1.5e-12)
			satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
		else
			satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
	}
	axnl = ep * cos(argpp);
	temp = 1.0 / (am * (1.0 - ep * ep));
	aynl = ep * sin(argpp) + temp * satrec.aycof;
	xl = mp + argpp + nodep + temp * satrec.xlcof * axnl;

	/* --------------------- solve kepler's equation --------------- */
	u = fmod(xl - nodep, twopi);
	eo1 = u;
	tem5 = 9999.9;
	ktr = 1;
	//   sgp4fix for kepler iteration
	//   the following iteration needs better limits on corrections
	while ((fabs(tem5) >= 1.0e-12) && (ktr <= 10))
	{
		sineo1 = sin(eo1);
		coseo1 = cos(eo1);
		tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
		tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
		if (fabs(tem5) >= 0.95)
			tem5 = tem5 > 0.0 ? 0.95 : -0.95;
		eo1 = eo1 + tem5;
		ktr = ktr + 1;
	}

	/* ------------- short period preliminary quantities ----------- */
	ecose = axnl * coseo1 + aynl * sineo1;
	esine = axnl * sineo1 - aynl * coseo1;
	el2 = axnl * axnl + aynl * aynl;
	pl = am * (1.0 - el2);
	if (pl < 0.0)
	{
		//         printf("# error pl %f\n", pl);
		satrec.error = 4;
	}
	else
	{
		rl = am * (1.0 - ecose);

		if (am >= 0.)
			rdotl = sqrt(am) * esine / rl;
		else
			rdotl = 0.;

		if (pl >= 0.)
			rvdotl = sqrt(pl) / rl;
		else
			rvdotl = 0.;

		double sqrt_argument = 1.0 - el2;
		if (sqrt_argument >= 0.)
			betal = sqrt(sqrt_argument);
		else
			betal = 0.;

		temp = esine / (1.0 + betal);
		sinu = am / rl * (sineo1 - aynl - axnl * temp);
		cosu = am / rl * (coseo1 - axnl + aynl * temp);
		su = atan2_check(sinu, cosu);
		sin2u = (cosu + cosu) * sinu;
		cos2u = 1.0 - 2.0 * sinu * sinu;
		temp = 1.0 / pl;
		temp1 = 0.5 * j2 * temp;
		temp2 = temp1 * temp;

		/* -------------- update for short period periodics ------------ */
		if (satrec.deep_space)
		{
			cosisq = cosip * cosip;
			satrec.con41 = 3.0*cosisq - 1.0;
			satrec.x1mth2 = 1.0 - cosisq;
			satrec.x7thm1 = 7.0*cosisq - 1.0;
		}
		mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) +
			0.5 * temp1 * satrec.x1mth2 * cos2u;
		su = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
		xnode = nodep + 1.5 * temp2 * cosip * sin2u;
		xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
		mvt = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / xke;
		rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
			1.5 * satrec.con41) / xke;

		/* --------------------- orientation vectors ------------------- */
		sinsu = sin(su);
		cossu = cos(su);
		snod = sin(xnode);
		cnod = cos(xnode);
		sini = sin(xinc);
		cosi = cos(xinc);
		xmx = -snod * cosi;
		xmy = cnod * cosi;
		ux = xmx * sinsu + cnod * cossu;
		uy = xmy * sinsu + snod * cossu;
		uz = sini * sinsu;
		vx = xmx * cossu - cnod * sinsu;
		vy = xmy * cossu - snod * sinsu;
		vz = sini * cossu;

		/* --------- position and velocity (in km and km/sec) ---------- */
		r[0] = (mrt * ux)* radiusearthkm;
		r[1] = (mrt * uy)* radiusearthkm;
		r[2] = (mrt * uz)* radiusearthkm;
		v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
		v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
		v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
	}  // if pl > 0

  // sgp4fix for decaying satellites
	if (mrt < 1.0)
	{
		//         printf("# decay condition %11.6f \n",mrt);
		satrec.error = 6;
	}

	//#include "debug7.cpp"
	return satrec.error;
}  // end sgp4

/* -----------------------------------------------------------------------------
*
*                           function gstime
*
*  this function finds the greenwich sidereal time.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jdut1       - julian date in ut1             days from 4713 bc
*
*  outputs       :
*    gstime      - greenwich sidereal time        0 to 2pi rad
*
*  locals        :
*    temp        - temporary variable for doubles   rad
*    tut1        - julian centuries from the
*                  jan 1, 2000 12 h epoch (ut1)
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2004, 191, eq 3-45
* --------------------------------------------------------------------------- */

double gstime(double jdut1)
{
	const double twopi = pi2;
	const double deg2rad = pi / 180.0;
	double       temp, tut1;

	tut1 = (jdut1 - 2451545.0) / 36525.0;
	temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
		(876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
	temp = fmod(temp * deg2rad / 240.0, twopi); //360/86400 = 1/240, to deg, to rad

	// ------------------------ check quadrants ---------------------
	if (temp < 0.0)
		temp += twopi;

	return temp;
}  // end gstime

double atan2_check(double y, double x)
{
	double Result_rad = 0.;

	if (x != 0.)
	{
		Result_rad = atan2(y, x);
	}
	else
	{
		if (y >= 0.) Result_rad = pi / 2.;
		else
			Result_rad = -pi / 2.;
	}
	return Result_rad;
}
//------------------------------------------------------------------------------

/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/

  //static
void dspace
(
	int irez,
	double d2201, double d2211, double d3210, double d3222, double d4410,
	double d4422, double d5220, double d5232, double d5421, double d5433,
	double dedt, double del1, double del2, double del3, double didt,
	double dmdt, double dnodt, double domdt, double argpo, double argpdot,
	double t, double tc, double gsto, double xfact, double xlamo,
	double no,
	double& atime, double& em, double& argpm, double& inclm, double& xli,
	double& mm, double& xni, double& nodem, double& dndt, double& nm
)
{
	const double twopi = pi2;
	int iretn, iret;
	double delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, g22, g32,
		g44, g52, g54, fasx2, fasx4, fasx6, rptim, step2, stepn, stepp;

	ft = 0.0;
	fasx2 = 0.13130908;
	fasx4 = 2.8843198;
	fasx6 = 0.37448087;
	g22 = 5.7686396;
	g32 = 0.95240898;
	g44 = 1.8014998;
	g52 = 1.0508330;
	g54 = 4.4108898;
	rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
	stepp = 720.0;
	stepn = -720.0;
	step2 = 259200.0;

	/* ----------- calculate deep space resonance effects ----------- */
	dndt = 0.0;
	theta = fmod(gsto + tc * rptim, twopi);
	em = em + dedt * t;

	inclm = inclm + didt * t;
	argpm = argpm + domdt * t;
	nodem = nodem + dnodt * t;
	mm = mm + dmdt * t;

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//    inclm = -inclm;
	//    argpm = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted

	ft = 0.0;
	atime = 0.0;
	if (irez != 0)
	{
		if ((atime == 0.0) || ((t >= 0.0) && (atime < 0.0)) ||
			((t < 0.0) && (atime >= 0.0)))
		{
			if (t >= 0.0)
				delt = stepp;
			else
				delt = stepn;
			atime = 0.0;
			xni = no;
			xli = xlamo;
		}
		iretn = 381; // added for do loop
		iret = 0; // added for loop
		while (iretn == 381)
		{
			if ((fabs(t) < fabs(atime)) || (iret == 351))
			{
				if (t >= 0.0)
					delt = stepn;
				else
					delt = stepp;
				iret = 351;
				iretn = 381;
			}
			else
			{
				if (t > 0.0)  // error if prev if has atime:=0.0 and t:=0.0 (ge)
					delt = stepp;
				else
					delt = stepn;
				if (fabs(t - atime) >= stepp)
				{
					iret = 0;
					iretn = 381;
				}
				else
				{
					ft = t - atime;
					iretn = 0;
				}
			}

			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if (irez != 2)
			{
				xndt = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) +
					del3 * sin(3.0 * (xli - fasx6));
				xldot = xni + xfact;
				xnddt = del1 * cos(xli - fasx2) +
					2.0 * del2 * cos(2.0 * (xli - fasx4)) +
					3.0 * del3 * cos(3.0 * (xli - fasx6));
				xnddt = xnddt * xldot;
			}
			else
			{
				/* --------- near - half-day resonance terms -------- */
				xomi = argpo + argpdot * atime;
				x2omi = xomi + xomi;
				x2li = xli + xli;
				xndt = d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +
					d3210 * sin(xomi + xli - g32) + d3222 * sin(-xomi + xli - g32) +
					d4410 * sin(x2omi + x2li - g44) + d4422 * sin(x2li - g44) +
					d5220 * sin(xomi + xli - g52) + d5232 * sin(-xomi + xli - g52) +
					d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54);
				xldot = xni + xfact;
				xnddt = d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +
					d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +
					d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +
					2.0 * (d4410 * cos(x2omi + x2li - g44) +
						d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
						d5433 * cos(-xomi + x2li - g54));
				xnddt = xnddt * xldot;
			}

			/* ----------------------- integrator ------------------- */
			if (iretn == 381)
			{
				xli = xli + xldot * delt + xndt * step2;
				xni = xni + xndt * delt + xnddt * step2;
				atime = atime + delt;
			}
		}  // while iretn = 381

		nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
		xl = xli + xldot * ft + xndt * ft * ft * 0.5;
		if (irez != 1)
		{
			mm = xl - 2.0 * nodem + 2.0 * theta;
			dndt = nm - no;
		}
		else
		{
			mm = xl - nodem - argpm + theta;
			dndt = nm - no;
		}
		nm = no + dndt;
	}

	//#include "debug4.cpp"
}  // end dsspace


double GAST_rad(double MJD)
{
	//GAST/sideral time, increases over time
	double Result = 0.;
	long double Diff_time = MJD - 33282.;
	long double Parameter_1 = (100.075542 + 360.*(Diff_time - double(int(Diff_time))) + 0.985647346*Diff_time + 2.9e-13*std::pow(Diff_time, 2.)) / 360.;
	long double Parameter_2;
	Result = modfl(Parameter_1, &Parameter_2)*pi2;
	return Result;//rad
}

v33 rz(double Angle_rad)
{
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	v33 out;
	out.v[0][0] = c;  out.v[0][1] = s;  out.v[0][2] = 0.0;
	out.v[1][0] = -s;  out.v[1][1] = c;  out.v[1][2] = 0.0;
	out.v[2][0] = 0.0;  out.v[2][1] = 0.0;  out.v[2][2] = 1.0;
	return out;
}

v33 ry(double Angle_rad)
{
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	v33 out;
	out.v[0][0] = c;  out.v[0][1] = 0.0;  out.v[0][2] = -s;
	out.v[1][0] = 0.0;  out.v[1][1] = 1.0;  out.v[1][2] = 0.0;
	out.v[2][0] = s;  out.v[2][1] = 0.0;  out.v[2][2] = c;
	return out;
}

ae topo_AzEl_rad(//topocentric
	v3 StaTCSm,
	ae StaTCSrad,
	v3 PointTCSm
)
{
	//Azimuth - from North (true North) CW towards East
	//North Az=0, East Az=90, South Az=180, West Az=270
	v3 toptrt = PointTCSm - StaTCSm;
	v33 m3 = rz(StaTCSrad.az);
	v33 m2 = ry(pi / 2. - StaTCSrad.el);
	v3 out = m2 * m3*toptrt;
	out.x *= -1.;//because azimuth is counted from North vector towards East
	ae Out_rad = get_ae(out);
	return Out_rad;
}

ae get_ae(v3 v)
{
	ae Vrad;
	Vrad.az = 0.;
	Vrad.el = 0.;

	double sqrt_argument = v.x*v.x + v.y*v.y;
	double Rxy = 0.;
	if (sqrt_argument >= 0.)
		Rxy = sqrt(sqrt_argument);

	if (Rxy < 10.e-20)
		Rxy = 0.;

	if (Rxy > 0.0)
	{
		double Argument = v.x / Rxy;
		//if ((Argument < -1.) || (Argument > 1.)) Msg("acos error avi_59");
		Vrad.az = acos(Argument);
		if (v.y < 0.0)
			Vrad.az = pi2 - Vrad.az;
	}

	sqrt_argument = std::pow(v.x, 2.) + std::pow(v.y, 2.) + std::pow(v.z, 2.);
	double Rxyz = 0.;
	if (sqrt_argument >= 0.)
		Rxyz = sqrt(sqrt_argument);

	if (Rxyz > 0.)
	{
		double Argument = v.z / Rxyz;
		//if ((Argument < -1.) || (Argument > 1.)) Msg("acos error main_4671");
		Vrad.el = pi / 2. - acos(Argument);
	}
	return Vrad;
}

double IArad(v3 A, v3 B)
{
	//the angle between two vectors is given by acos of the dot product of the two (normalised) vectors
	v3 Auv = normalize(A);
	v3 Buv = normalize(B);
	return acos(dot(Auv, Buv));
}

v3 normalize(v3 v)//returns unit vector
{
	v3 Out = v;
	double Length = norm(v);
	if (Length > 0.)
	{
		Out.x /= Length;
		Out.y /= Length;
		Out.z /= Length;
	}
	return Out;
}

double norm(v3 v)//returns vector length
{
	return std::pow(v.x*v.x + v.y*v.y + v.z*v.z, 0.5);
}

double dot(v3 left, v3 right)
{
	double Sum = left.x*right.x + left.y*right.y + left.z*right.z;
	return Sum;
}

v3 get_v3(ae Vrad)//get unit vector at this orientation
{
	double Radius_XY = fabs(cos(Vrad.el));
	v3 Out;
	Out.x = Radius_XY * cos(Vrad.az);
	Out.y = Radius_XY * sin(Vrad.az);
	Out.z = sin(Vrad.el);
	return Out;
}

ae Sun_ICS_rad(double MJD)
{
	double JD = MJD + 2400000.5;
	double Epsilon = CAANutation::TrueObliquityOfEcliptic(JD);

	//sun apparent equatorial coordinates
	double Apparent_ecliptic_Lon_deg = CAASun::ApparentEclipticLongitude(JD, false);
	double Apparent_ecliptic_Lat_deg = CAASun::ApparentEclipticLatitude(JD, false);
	CAA2DCoordinate Sun1 = CAACoordinateTransformation::Ecliptic2Equatorial(Apparent_ecliptic_Lon_deg, Apparent_ecliptic_Lat_deg, Epsilon);

	ae ICSrad;
	ICSrad.az = CAACoordinateTransformation::HoursToRadians(Sun1.X);
	ICSrad.el = CAACoordinateTransformation::DegreesToRadians(Sun1.Y);

	return ICSrad;
}

double Distance_to_the_SUN_AU(double MJD)
{
	//astronomical algorithms, chapter 24, solar coordinates
	double Distance = 0.;
	double JD = MJD + 2400000.5;
	double T = (JD - 2451545.) / 36525.;
	double T2 = T * T;
	double T3 = T * T * T;
	// double L0=280.46645+36000.76983*T+0.0003032*T2;//deg
	double M_deg = 357.52910 + 35999.05030 * T - 0.0001559 * T2 - 0.00000048 * T3;//deg
	double M_rad = M_deg * deg2rad;
	double e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T2;
	double C_deg = (1.914600 - 0.004817 * T - 0.000014 * T2) * sin(M_rad) + (0.019993 - 0.000101 * T) * sin(2. * M_rad) + 0.000290 * sin(3. * M_rad);
	double C_rad = C_deg * deg2rad;
	double v_rad = M_rad + C_rad;
	double Numerator = 1.000001018 * (1. - e * e);
	double Denominator = 1 + e * cos(v_rad);
	if (Denominator != 0.)
		Distance = Numerator / Denominator;
	return Distance;
}

v3 get_v3(ae Vrad, double Length)
{
	double Radius_XY = fabs(cos(Vrad.el));
	v3 Out;
	Out.x = Length * Radius_XY*cos(Vrad.az);
	Out.y = Length * Radius_XY*sin(Vrad.az);
	Out.z = Length * sin(Vrad.el);
	return Out;
}

double shadow_function(double MJD, v3 Sat_ICS_m, v3 Sun_ICS_m)
{
	double Out_ShadowFunction = 0.;
	//1=sun, penumbra, 0=umbra/shadow

	double Earth_radius_km = 6378.137;
	double Sun_radius_km = 696000.;

	double r_km = norm(Sat_ICS_m) / 1000.;

	if (r_km > Earth_radius_km)
	{
		//Sun ICS
		double Distance_Sun_Earth_AU = Distance_to_the_SUN_AU(MJD);
		double AU_2_km = 149597870.7;//1 AU = km
		double Distance_Sun_Earth_km = Distance_Sun_Earth_AU * AU_2_km;

		double valueIArad = IArad(Sun_ICS_m, Sat_ICS_m);//always positive, result is 0..PI

		double Sc = r_km * cos(valueIArad) + std::pow(r_km*r_km - Earth_radius_km * Earth_radius_km, 0.5);//eq5

		//in the paper there is a mistake, no need for 'r' in the denominator
		double Alpha_rad = atan((Sun_radius_km - Earth_radius_km) / Distance_Sun_Earth_km);
		double Beta_rad = atan((Sun_radius_km + Earth_radius_km) / Distance_Sun_Earth_km);

		double argument = std::pow(std::pow(r_km, 2.) - std::pow(Earth_radius_km*cos(Alpha_rad), 2.), 0.5);
		double Su = r_km * cos(valueIArad) + cos(Alpha_rad)*(argument + Earth_radius_km * sin(Alpha_rad));//eq6

		argument = std::pow(std::pow(r_km, 2.) - std::pow(Earth_radius_km*cos(Beta_rad), 2.), 0.5);
		double Sp = r_km * cos(valueIArad) + cos(Beta_rad)*(argument - Earth_radius_km * sin(Beta_rad));//eq7

		double Delta_h = Su - Sp;

		//penumbra function
		double delta = 8;//was 8
	   //  double vp=0.5*(1.+tanh(delta*PI2*Earth_radius_km*Sc/Delta_h));//original eq13
		double vp2 = 0.5*(1. + tanh(delta*Sc / Delta_h));//without constant factor of 'PI2*Earth_radius_km'

		if (Su <= 0.) Out_ShadowFunction = 0.; else
			if (Sp <= 0.) Out_ShadowFunction = vp2; else
				Out_ShadowFunction = 1.;
	}

	return Out_ShadowFunction;
}

v33 inv(v33 M)//inverse of 3x3 matrix
{
	v33 A;
	A.v[0][0] = 0.;
	A.v[0][1] = 0.;
	A.v[0][2] = 0.;
	A.v[1][0] = 0.;
	A.v[1][1] = 0.;
	A.v[1][2] = 0.;
	A.v[2][0] = 0.;
	A.v[2][1] = 0.;
	A.v[2][2] = 0.;

	double D = det(M);//get determinant

	if (D != 0.)
	{
		double a11 = M.v[0][0];
		double a12 = M.v[0][1];
		double a13 = M.v[0][2];
		double a21 = M.v[1][0];
		double a22 = M.v[1][1];
		double a23 = M.v[1][2];
		double a31 = M.v[2][0];
		double a32 = M.v[2][1];
		double a33 = M.v[2][2];

		A.v[0][0] = (a22*a33 - a23 * a32) / D;
		A.v[0][1] = (a13*a32 - a12 * a33) / D;
		A.v[0][2] = (a12*a23 - a13 * a22) / D;

		A.v[1][0] = (a23*a31 - a21 * a33) / D;
		A.v[1][1] = (a11*a33 - a13 * a31) / D;
		A.v[1][2] = (a13*a21 - a11 * a23) / D;

		A.v[2][0] = (a21*a32 - a22 * a31) / D;
		A.v[2][1] = (a12*a31 - a11 * a32) / D;
		A.v[2][2] = (a11*a22 - a12 * a21) / D;
	}
	return A;
}

double det(v33 M)
{
	//calculate determinant of 3x3 matrix	
	double D =
		M.v[0][0] * M.v[1][1] * M.v[2][2] -
		M.v[0][0] * M.v[1][2] * M.v[2][1] -
		M.v[0][1] * M.v[1][0] * M.v[2][2] +
		M.v[0][1] * M.v[1][2] * M.v[2][0] +
		M.v[0][2] * M.v[1][0] * M.v[2][1] -
		M.v[0][2] * M.v[1][1] * M.v[2][0];

	return D;
}
