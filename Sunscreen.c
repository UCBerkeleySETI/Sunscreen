//Sunscreen
//=========
//v1.21 - 29 October 2018
//Writted by Brian C. Lacki for Breakthrough Listen
//A program that calculates unextincted spectra of model stellar populations with
//all stars below or above a threshold quantity screened by megastructures or some
//other way.  For example, it can calculate what a galaxy would look like if all stars
//with L* < 1,000 L_sun are screened.
//It also calculates photometric magnitudes of the screened stellar population if
//CalculateMagnitudes is set to TRUE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define CalculateMagnitudes TRUE
//Number of distinct "filters" or photometric bandpasses (besides bolometric) that the 
//program calculates absolute magnitudes in.
#define NFilt 48
#define NRedshiftOutput 2

//Number of distinct effective temperatures for which model spectra exist
#define NT 68
//Number of distinct surface gravities for which model spectra exist
//(Though not all values are available for all T_eff.)
#define NG 19
//Number of wavelength values in model spectra
//I assume that all the model spectra use the same wavelength values.
//This is true for the Lejeune spectra
#define NLam 1221

//Number of distinct population ages for the stellar ages.
//614 corresponds to log_10 t_iso = 4.00 to 10.13 inclusive, with step size of
//dlog_10 t_iso = 0.01.
//I assume this value is the same for all metallicities
#define NtIso 614

//Number of kinds of ways to divide the stellar population into screened and unscreened
//populations.
//Currently the 4 ways are: luminosity, luminosity-to-mass, mass, and stellar stage
#define NThresh 4
//Number of distinct metallicities.
//Currently, they are: log_10 (Z/Z_sun) = 0.0, -1.0, +0.301
#define NZ 3
//Number of distinct IMFs for the stellar populations
//Currently, they are: Chabrier, Salpeter, bottom-heavy, 
//and Chabrier with a M* > 1 M_sun, 2 M_sun, and 5 M_sun
#define NIMF 3
//Number of distinct SFHs for the stellar populations
//Currently they are: constant SFR, four SFHs from Behroozi et al. 2013,
//one SFH from McDiamird et al. 2015, a single-age burst, and a 
#define NSFH 8
//Number of "ages" of the stellar population, measured from the Big Bang
#define NtPop 9
//Number of temperatures of the Dyson spheres for which thermal emission
//spectra are calculated
//Currently, they are: 1000 K and 300 K.
#define ND 2

typedef struct
{
   double t;
   double M0;
   double M;
   double L;
   double Teff;
   double R;
   double logg;
   double A;
   int Stage;
   double Z;
} StellarData;

typedef struct
{
	int t0;
	int g0;
	int g1;
	double f00;
	double f01;
	double f10;
	double f11;
} sInterpConstants;

double pi = 3.14159265358979323846264;
double kB = 1.381e-16;
double c = 2.998e10;
double hPlanck = 6.626e-27;
double LSun = 3.826e33;
double sigmaSB = 5.671e-5;
double pc = 3.0857e18;
double LBol0 = 3.0128e35;
double yr = 3.1557e7;
float xPre_Planck = 6.626e-27 * 2.998e10 / 1.381e-16;
float FnuPre_RJ = 2.0 * 3.14159265358979323846264 * 1.381e-16;
float FnuPre_Planck = 2.0 * 3.14159265358979323846264 * 6.626e-27 * 2.998e10;
float Watt = 1e7;
double Omega_m_0 = 0.315;
double H_0 = 67.31 * 1e5 / 3.0857e24;
double t0_Cosmic = 1.381125e10 * 3.1557e7; //Lookback time to z = 1000 in seconds
double G_N = 6.674e-8;
double TCMB_0 = 2.73;

enum IMFList {chabrier, salpeter, botheavy, chabrier_minmass_1msun, chabrier_minmass_2msun, chabrier_minmass_5msun};
enum StageList {prems, ms, sgb, rgb, cehb1, cehb2, cehb3, eagb, tpagb, remnant};
enum ThreshTypeList {luminosity, lumtomass, mass, stellartype};
enum OutputMagnitudeType {rest, observed_redshifted};

int CountLinesInFile (FILE *f);
void ReadLamList (char *FileName, float *Spectrum);
void ReadSpectrum (char *FileName, float *Spectrum);
void ReadIsochrone (char *FileName, StellarData **Isochrone, int *NMasses);
void ReadSFH (char *FileName, int *NtSFH, double **t, double **SFH);

double SFH (double tBeginInt, double tEndInt);
double BurstSFH (double tBeginInt, double tEndInt, double tBurstStart, double tBurstEnd, double MStarBurst);
double CtsSFH (double tBeginInt, double tEndInt, double tSFStart, double SFR);
double TableSFH (double tBeginInt, double tEndInt, int NtSFH, double *tTable, double *SFRTable);

double IMF (double M, double IMFNorm, int IMFVariation);
double ChabrierIMF (double M, double IMFNorm);
double SalpeterIMF (double M, double IMFNorm);
double BotHeavyIMF (double M, double IMFNorm);
double IMFMassIntegral (double M1, double M2, double IMFNorm, int IMFVariation);
double IMFNumberIntegral (double M1, double M2, double IMFNorm, int IMFVariation);

float FnuPlanck (float lambda, float T);

int BinarySearchForFirstGTE (int N, double *Table, double x);

void ReadFilterN0 (char *FileName, int *NFLam, double **Lambda, double **N0);
double MBolObs (double *Lambda, double *Spectrum);
double InterpolatedN0 (double lambda, int NFLam, double *FilterLambda, double *FilterN0);

double RedshiftAtLookbackTime (double tL, double H0, double Omega_M, double T_CMB);
double rho_c (double H0);
double Omega_CMB (double H0, double T);

void PrintMagnitudeHeader (FILE *OutFile, char ColumnNames[NFilt][100]);

int main (void)
{
	int i_TStar, i_gStar, i_MStar, iCompact_TStar, iCompact_gStar;
	int i_IsoAge, i_ThreshL, i_lambda, i_IMF, i_DysonTemp, i_ThreshType, i_GalaxyAge, i_Z, i_SFH;
	int i_Filter, i_Redshift, i_ObsSpectrum;
	
	double t, logt;
	double tGalaxyAge[NtPop];
	double logtGalaxyAge[NtPop] = {6.00, 7.00, 8.00, 9.00, 9.20, 9.50, 9.80, 10.10, 10.13};

	double logtIsoMin = 4.00; 
	double logtIsoMax = 10.13; 
	double dlogtIso = 0.01;
	double logtIsoList[NtIso];
	StellarData *Isochrones[NtIso];
	StellarData *ThisStar;
	char StellarStageNames[10][10] = {"PreMS", "MS", "SGB", "RGB", "CEHB1", "CEHB2", "CEHB3", "EAGB", "TPAGB", "Remnant"};
	
	double logZ[NZ] = {0.00, -1.00, 0.30};
	
	int NMasses[NtIso];
	double MStar_0, AStar, LStar, TStar, MStar, loggStar;
	int StageStar;
	
	char lambdaListFileName[255] = "LejeuneSpectra/lambda-LejeuneSpectra.txt";
	
	double HotSpectrum[NLam];
	double ColdSpectrum[NLam];
	double MStarMin = 0.1;
	double MStarMax = 120.0;
	double IMFNorm[NIMF];
	double MLo, MHi;
	double tLo[NtIso], tHi[NtIso];
	double **dNdMIso[NIMF];
	double **dMIso[NSFH];
	double ThisdMIso;
	float ThisdNdMIso;
	
	int Ng;
	//int	iT;
	int igCold, igHot;
	double fTCold, fTHot;
	double fgColdLo, fgColdHi, fgHotLo, fgHotHi;
	double f00, f01, f10, f11;

	float lamList[NLam];
	double dnuList[NLam];
	double FullTGrid[NT] = {2000, 2200, 2500, 2800, 3000,\
					  3200, 3350, 3500, 3750, 4000,\
 				      4250, 4500, 4750, 5000, 5250,\
					  5500, 5750, 6000, 6250, 6500,\
					  6750, 7000, 7250, 7500, 7750,\
					  8000, 8250, 8500, 8750, 9000,\
					  9250, 9500, 9750, 10000, 10500,\
					  11000, 11500, 12000, 12500, 13000,\
					  14000, 15000, 16000, 17000, 18000,\
					  19000, 20000, 21000, 22000, 23000,\
					  24000, 25000, 26000, 27000, 28000,\
					  29000, 30000, 31000, 32000, 33000,\
					  34000, 35000, 37500, 40000, 42500,\
					  45000, 47500, 50000};
	double FullloggGrid[NG] = {-1.02, -0.70, -0.51, -0.29, +0.00,\
	                      +0.28, +0.50, +0.60, +0.87, +1.00,\
						  +1.50, +2.00, +2.50, +3.00, +3.50,\
						  +4.00, +4.50, +5.00, +5.50};
	int GridModelExists[NT][NG];
	int TExists[NT];
	int nTForZ;
	int *ngAtT;
	double *TGrid;
	double **loggGrid;
	float ***ModelSpectrum;
	float *ModelSpectrumPointer_00, *ModelSpectrumPointer_01, *ModelSpectrumPointer_10, *ModelSpectrumPointer_11;
	
	char ModelSpectrumName[255];
	char IsochroneName[255];
	FILE *InFile;

	float *ThisStellarSpectrum;
	float **IndividualStellarSpectrum[NtIso];
	float **IndividualStellarMagnitudes[NtIso];
	int *IncludeThisStar[NtIso];
	int *i_ThreshL_Cut[NThresh][NtIso];
	int i_ThisCut;

	float *ThisIsochroneSpectrum;			
	float *IsochroneUnscreenedSpectrum[NtIso], *IsochroneScreenedSpectrum_NoDim[NtIso], *IsochroneScreenedSpectrum_NoBright[NtIso];
	double IsochroneLDyson_NoDim[NtIso], IsochroneLDyson_NoBright[NtIso];
	
	float IntegratedUnscreenedSpectrum[NLam], IntegratedScreenedSpectrum_NoDim[NLam], IntegratedScreenedSpectrum_NoBright[NLam];
	double IntegratedLDyson_NoDim, IntegratedLDyson_NoBright;
	
	char ModelIDString[200];
	char DimScreenedSpectrumName[255], BrightScreenedSpectrumName[255], UnscreenedSpectrumName[255];
	FILE *DimScreenedSpectrumFile, *BrightScreenedSpectrumFile, *UnscreenedSpectrumFile;
	
	char ThreshName[NThresh][15] = {"logL", "logLtoM0", "logM0", "Stage"};
	int NL[NThresh] = {99, 69, 33, 10};
	double logThreshMin[NThresh] = {-3.5, -2.4, -1.1, 0.};
	double dlogThresh[NThresh] = {0.1, 0.1, 0.1, 1.};
	double logThresh;
		
	char SFModeName[NSFH][30] = {"ConstSF", "BurstSF", "B13-Mh+11.0", "B13-Mh+12.0", "B13-Mh+13.0", "B13-Mh+14.0", "M15-MJAM+11.0", "W14-dIrr"};
	char SFHTableFileNames[NSFH - 2][255] = {"SFRs/B13-log10Mh+11.0.txt", "SFRs/B13-log10Mh+12.0.txt", "SFRs/B13-log10Mh+13.0.txt", "SFRs/B13-log10Mh+14.0.txt", "SFRs/M15-MJAM+11.0.txt", "SFRs/W14-dIrr.txt"};	

	int NtSFH[NSFH - 2];
	double *tSFHTable[NSFH - 2], *SFHTable[NSFH - 2];
	char IMFName[NIMF][50] = {"Chabrier", "Salpeter", "BottomHeavy"};
	double tSFGyr;
	double BurstSF_MStar_0 = 1e11;
	double ConstSF_SFR = 1.;
	
	double *dMPopdMIso[NIMF];
	
	double MPop;
	FILE *MassFile;
	char MassFileName[255];
	
	double TDysonList[ND] = {1000.0, 300.0};
	double BlackbodyFlux[ND][NLam];
	double ADysonDimScreened[ND], ADysonBrightScreened[ND];
	double ThisTDyson4;
	
	char FilterName[NFilt][100] = {"FUV", "NUV", \
							   "U_JC", "B_JC", "V_JC", "R_JC", "I_JC", \
							   "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS", "z_SDSS", \
							   "g_DES", "r_DES", "i_DES", "z_DES", "Y_DES", \
							   "g_PS1", "r_PS1", "i_PS1", "z_PS1", "y_PS1", "w_PS1", \
							   "u_LSST", "g_LSST", "r_LSST", "i_LSST", "z_LSST", "y_LSST", \
							   "G_Gaia", "GBP_Gaia", "GRP_Gaia",  \
							   "J_2MASS", "H_2MASS", "Ks_2MASS", \
							   "Z_UKIRT", "Y_UKIRT", "J_UKIRT", "H_UKIRT", "K_UKIRT", \
							   "IRAC-1", "IRAC-2", "IRAC-3", "IRAC-4", \
							   "W1_WISE", "W2_WISE", "W3_WISE", "W4_WISE"};
	double DFactor_AbsMagnitude = 4.0 * pi * 100.0 * pc * pc;
	char FilterN0Name[255];
	int NFLam; 
	double *FilterLambda;
	double *FilterN0;
	
	double ObsSpectrumLambda[NtPop + 1][NLam];
	double ObsSpectrumLnu_Unscreened[NLam], ObsSpectrumLnu_NoDim[NLam], ObsSpectrumLnu_NoBright[NLam];
	double *FilterKernel[NFilt][NtPop + 1];
	double Redshift_GalaxyAge[NtPop + 1];
	double aScale;
	char zOutputName[NRedshiftOutput][20] = {"Rest", "Obs"};
	double nu_Blue, nu_Red;
	int i_FilterLambda_Blue, i_FilterLambda_Red;
	int i_FilterLambda;
	
	double *FilterKernelPointer;
	double *ObsLambdaListPointer;
	double lambda_Blue, lambda_Red;
	double dlambda;
	double Integrand_Blue, Integrand_Red;
	double ThisFilterKernel;
	
	double RatioSum_Unscreened, RatioSum_NoDim, RatioSum_NoBright;
	int i_lambda_P1;
	
	char DimScreenedMagnitudeName[255], BrightScreenedMagnitudeName[255], UnscreenedMagnitudeName[255];
	FILE *DimScreenedMagnitudeFile[NSFH][NtPop][NRedshiftOutput], *BrightScreenedMagnitudeFile[NSFH][NtPop][NRedshiftOutput];
	FILE *UnscreenedMagnitudeFile;
	FILE *MagnitudeFilePointer_NoDim, *MagnitudeFilePointer_NoBright;
	
	//Memory allocation and metallicity-independent pre-calculations
	//==============================================================
	for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
	{
		logtIsoList[i_IsoAge] = logtIsoMin + (double) i_IsoAge * dlogtIso;
		IsochroneUnscreenedSpectrum[i_IsoAge] = calloc (NLam, sizeof(float));
		IsochroneScreenedSpectrum_NoDim[i_IsoAge] = calloc (NLam, sizeof(float));		
		IsochroneScreenedSpectrum_NoBright[i_IsoAge] = calloc (NLam, sizeof(float));
		
		if (CalculateMagnitudes == FALSE) continue;
	}
	
	for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
	{
		dMIso[i_SFH] = calloc (NtPop, sizeof(double*));
		for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++) dMIso[i_SFH][i_GalaxyAge] = calloc (NtIso, sizeof(double));
	}
	
	for (i_IMF = 0; i_IMF < NIMF; i_IMF++)
	{
		dNdMIso[i_IMF] = calloc (NtIso, sizeof(double**));
		dMPopdMIso[i_IMF] = calloc (NtIso, sizeof(double*));
	}
	
	ReadLamList (lambdaListFileName, lamList);
	for (i_lambda = 0; i_lambda < NLam; i_lambda++) 
	{
		if (i_lambda > 0) dnuList[i_lambda] = c / lamList[i_lambda - 1] - c / lamList[i_lambda];
		for (i_DysonTemp = 0; i_DysonTemp < ND; i_DysonTemp++) BlackbodyFlux[i_DysonTemp][i_lambda] = FnuPlanck (lamList[i_lambda], TDysonList[i_DysonTemp]);
		ObsSpectrumLambda[0][i_lambda] = (double) lamList[i_lambda];
	}

	Redshift_GalaxyAge[0] = 0.;
	for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
	{
		tGalaxyAge[i_GalaxyAge] = pow(10., logtGalaxyAge[i_GalaxyAge]);
	
		Redshift_GalaxyAge[i_GalaxyAge + 1] = RedshiftAtLookbackTime(t0_Cosmic - tGalaxyAge[i_GalaxyAge] * yr, H_0, Omega_m_0, TCMB_0);

		//ObsSpectrumLambda is the wavelength of the redshifted spectrum from the galaxy,
		//assuming it is at cosmic age tGalaxyAge 
		for (i_lambda = 0; i_lambda < NLam; i_lambda++) ObsSpectrumLambda[i_GalaxyAge + 1][i_lambda] = (double) lamList[i_lambda] * (1. + Redshift_GalaxyAge[i_GalaxyAge + 1]); 
	}
	
	if (CalculateMagnitudes)
	{
		for (i_Filter = 0; i_Filter < NFilt; i_Filter++)
		{
			sprintf (FilterN0Name, "FilterResponses/PhotonZeroPoint-%s.txt", FilterName[i_Filter]);
			printf ("%d %s\n", i_Filter, FilterN0Name);
			
			ReadFilterN0 (FilterN0Name, &NFLam, &FilterLambda, &FilterN0);

			for (i_GalaxyAge = 0; i_GalaxyAge < NtPop + 1; i_GalaxyAge++)
			{
				FilterKernel[i_Filter][i_GalaxyAge] = calloc (NLam, sizeof(double));
				//FilterKernel is convolved with the spectrum to calculate the magnitudes.
				i_FilterLambda_Blue = i_FilterLambda_Red = 0;
				for (i_lambda = 0; i_lambda < NLam; i_lambda++)
				{
					lambda_Blue = (i_lambda == 0) ? ObsSpectrumLambda[i_GalaxyAge][0] : (0.5 * (ObsSpectrumLambda[i_GalaxyAge][i_lambda - 1] + ObsSpectrumLambda[i_GalaxyAge][i_lambda]));
					lambda_Red = (i_lambda == NLam - 1) ? ObsSpectrumLambda[i_GalaxyAge][i_lambda] : (0.5 * (ObsSpectrumLambda[i_GalaxyAge][i_lambda + 1] + ObsSpectrumLambda[i_GalaxyAge][i_lambda]));
					
					if (lambda_Red <= FilterLambda[0]) {FilterKernel[i_Filter][i_GalaxyAge][i_lambda] = 0.; continue;}
					if (lambda_Blue >= FilterLambda[NFLam - 1]) {FilterKernel[i_Filter][i_GalaxyAge][i_lambda] = 0.; continue;}

					//Index of first wavelength in FilterLambda >= lambda_Blue
					while ((lambda_Blue > FilterLambda[i_FilterLambda_Blue]) && (i_FilterLambda_Blue < NFLam - 1)) i_FilterLambda_Blue++;
					while ((lambda_Red > FilterLambda[i_FilterLambda_Red]) && (i_FilterLambda_Red < NFLam - 1)) i_FilterLambda_Red++;

					if (i_FilterLambda_Red == i_FilterLambda_Blue)
					{
						//Filter wavelengths coarser than model spectrum wavelengths
						FilterKernel[i_Filter][i_GalaxyAge][i_lambda] = InterpolatedN0 ((double) ObsSpectrumLambda[i_GalaxyAge][i_lambda], NFLam, FilterLambda, FilterN0) / (hPlanck * (double) ObsSpectrumLambda[i_GalaxyAge][i_lambda] * DFactor_AbsMagnitude);
					}
					else
					{
						//Wavelength-weighted
						Integrand_Blue = InterpolatedN0 (lambda_Blue, NFLam, FilterLambda, FilterN0) / (hPlanck * lambda_Blue * DFactor_AbsMagnitude);

						Integrand_Red = FilterN0[i_FilterLambda_Blue] / (hPlanck * FilterLambda[i_FilterLambda_Blue] * DFactor_AbsMagnitude);
						
						FilterKernel[i_Filter][i_GalaxyAge][i_lambda] = 0.5 * (Integrand_Blue + Integrand_Red) * (FilterLambda[i_FilterLambda_Blue] - lambda_Blue);
						
						Integrand_Blue = Integrand_Red;
						
						for (i_FilterLambda = i_FilterLambda_Blue; i_FilterLambda < i_FilterLambda_Red - 1; i_FilterLambda++)
						{			
							Integrand_Red = FilterN0[i_FilterLambda + 1] / (hPlanck * FilterLambda[i_FilterLambda + 1] * DFactor_AbsMagnitude);
							
							FilterKernel[i_Filter][i_GalaxyAge][i_lambda] += 0.5 * (Integrand_Blue + Integrand_Red) * (FilterLambda[i_FilterLambda + 1] - FilterLambda[i_FilterLambda]);	
							
							Integrand_Blue = Integrand_Red;
						}
						
						Integrand_Red = InterpolatedN0 (lambda_Red, NFLam, FilterLambda, FilterN0) / (hPlanck * lambda_Red * DFactor_AbsMagnitude);

						FilterKernel[i_Filter][i_GalaxyAge][i_lambda] += 0.5 * (Integrand_Blue + Integrand_Red) * (lambda_Red - FilterLambda[i_FilterLambda_Red - 1]);
						
						FilterKernel[i_Filter][i_GalaxyAge][i_lambda] /= (lambda_Red - lambda_Blue); 
					}
					
					
					//This version (1.2) of FilterKernel includes the dlambda from the integral over wavelength
					//when calculating magnitudes.
					//FilterKernel now DEPENDS ON THE WAY SUNSCREEN DOES INTEGRALS!
					//Note the integral in this code uses the trapezoidal rule.
					//Doing this means only only value of the integrand only needs to be evaluated at a time,
					//and a multiplication in the integral can be skipped.
					lambda_Blue = (i_lambda == 0) ? ObsSpectrumLambda[i_GalaxyAge][0] : ObsSpectrumLambda[i_GalaxyAge][i_lambda - 1];
					lambda_Red = (i_lambda == NLam - 1) ? ObsSpectrumLambda[i_GalaxyAge][i_lambda] : ObsSpectrumLambda[i_GalaxyAge][i_lambda + 1];
					dlambda = lambda_Red - lambda_Blue;
					FilterKernel[i_Filter][i_GalaxyAge][i_lambda] *= 0.5 * dlambda;	
					
				}
			}
			
			free (FilterLambda);
			free (FilterN0);
		}
	}
	
	//Universal lookup tables
	//=======================
	//For the Chabrier IMF, IMFNorm[0] gives the constant factor such that IMF() returns dN/dM for 1 M_sun/yr.
	//For the other IMFs, I normalize so that dN/dM (1 M_sun) equals its value for Chabrier.  This is because
	//I want the bottom-heavy tracks to coincide with the Chabrier ones for red galaxies, and red giants
	//generate most of the light in the early-type galaxies.
	IMFNorm[chabrier] = 1.0 / IMFMassIntegral (MStarMin, MStarMax, 1.0, chabrier);
	IMFNorm[salpeter] = IMFNorm[chabrier] * IMF (1.0, 1.0, chabrier) / IMF (1.0, 1.0, salpeter);
	IMFNorm[botheavy] = IMFNorm[botheavy] * IMF (1.0, 1.0, chabrier) / IMF (1.0, 1.0, botheavy);
	IMFNorm[chabrier_minmass_5msun] = IMFNorm[chabrier_minmass_2msun] = IMFNorm[chabrier_minmass_1msun] = IMFNorm[chabrier];

	for (i_SFH = 0; i_SFH < NSFH - 2; i_SFH++) ReadSFH(SFHTableFileNames[i_SFH], &(NtSFH[i_SFH]), &(tSFHTable[i_SFH]), &(SFHTable[i_SFH]));	
	

	//Main loop
	//=========	
	for (i_Z = 0; i_Z < NZ; i_Z++)
	{
		//Metallicity
		//Do metallicity loop first, since both stellar spectra and isochrones
		//depend on it, and need to be read in for each metallicity.
		
		//Isochrone read-in
		//=================
		for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
		{
			sprintf (IsochroneName, "CMD-Isochrones/Isochrone.logZ%+.2f.logt+%.2f.txt", logZ[i_Z], logtIsoList[i_IsoAge]);
			printf ("%s\n", IsochroneName);
			ReadIsochrone (IsochroneName, &(Isochrones[i_IsoAge]), &(NMasses[i_IsoAge]));
			if (i_IsoAge > 0) tLo[i_IsoAge] = tHi[i_IsoAge - 1] = 0.5 * (Isochrones[i_IsoAge][0].t + Isochrones[i_IsoAge - 1][0].t);	
		}
		tLo[0] = 0.0;
		tHi[NtIso - 1] = 1.5 * Isochrones[NtIso - 1][0].t - 0.5 * Isochrones[NtIso - 2][0].t;
		
		//Stellar population parameters
		//=============================
		for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
		{
			for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
			{
				//CANNOT simply let spectra accumulate with age, at least in the burst case.
				//There are two integrals over t -- one t is a parameter defining the SFH,
				//the other denotes stellar age within that history.
					
				t = tGalaxyAge[i_GalaxyAge];
				for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
				{
					//L_nu = Int d(L_nu)/dM*(t) * dM*(t)
					if (i_SFH == 0) dMIso[i_SFH][i_GalaxyAge][i_IsoAge] = CtsSFH (-tHi[i_IsoAge], -tLo[i_IsoAge], -t, ConstSF_SFR);
					if (i_SFH == 1) dMIso[i_SFH][i_GalaxyAge][i_IsoAge] = BurstSFH (-tHi[i_IsoAge], -tLo[i_IsoAge], -(t + 1.0), -t, BurstSF_MStar_0);
					//We want stars with ages between tLo and tHi at t.
					if (i_SFH >= 2) dMIso[i_SFH][i_GalaxyAge][i_IsoAge] = TableSFH (t - tHi[i_IsoAge], t - tLo[i_IsoAge], NtSFH[i_SFH - 2], tSFHTable[i_SFH - 2], SFHTable[i_SFH - 2]);
				}
			}
		}
		
		
		//Stellar spectra for individual stars
		//====================================
		nTForZ = 0;
		for (i_TStar = 0; i_TStar < NT; i_TStar++)
		{
			TExists[i_TStar] = FALSE;
			for (i_gStar = 0; i_gStar < NG; i_gStar++)
			{
				sprintf (ModelSpectrumName, "LejeuneSpectra/LejeuneSpectrum.logZ%+.2f.Teff+%.0f.logg%+.2f", logZ[i_Z], FullTGrid[i_TStar], FullloggGrid[i_gStar]);
				InFile = fopen (ModelSpectrumName, "r");
				if (InFile) 
				{
					GridModelExists[i_TStar][i_gStar] = TRUE; 
					fclose(InFile);
					TExists[i_TStar] = TRUE;
				} 
				else 
				{
					GridModelExists[i_TStar][i_gStar] = FALSE;
				}
			}
			nTForZ += TExists[i_TStar];
		}
		
		ngAtT = calloc (nTForZ, sizeof(int));
		TGrid = calloc (nTForZ, sizeof(double));
		loggGrid = calloc (nTForZ, sizeof(double*));
		ModelSpectrum = calloc (nTForZ, sizeof(float**));
		
		iCompact_TStar = 0;
		for (i_TStar = 0; i_TStar < NT; i_TStar++)
		{
			if (TExists[i_TStar] == FALSE) continue;
			TGrid[iCompact_TStar] = FullTGrid[i_TStar];
	
			for (i_gStar = 0; i_gStar < NG; i_gStar++) ngAtT[iCompact_TStar] += GridModelExists[i_TStar][i_gStar];

			loggGrid[iCompact_TStar] = calloc (NG, sizeof(double));
			ModelSpectrum[iCompact_TStar] = calloc (ngAtT[iCompact_TStar], sizeof(float*));

			iCompact_gStar = 0;
			for (i_gStar = 0; i_gStar < NG; i_gStar++)
			{
				if (GridModelExists[i_TStar][i_gStar] == TRUE)
				{
					ModelSpectrum[iCompact_TStar][iCompact_gStar] = calloc (NLam, sizeof(float));
					
					loggGrid[iCompact_TStar][iCompact_gStar] = FullloggGrid[i_gStar];
					sprintf (ModelSpectrumName, "LejeuneSpectra/LejeuneSpectrum.logZ%+.2f.Teff+%.0f.logg%+.2f", logZ[i_Z], TGrid[i_TStar], FullloggGrid[i_gStar]);
					printf ("%s\n", ModelSpectrumName);
					ReadSpectrum (ModelSpectrumName, ModelSpectrum[iCompact_TStar][iCompact_gStar]);
					iCompact_gStar++;
				}
				if (iCompact_gStar >= ngAtT[iCompact_TStar]) break;
			}
			iCompact_TStar++;
		}
		
	
		//Isochronic stellar population spectra
		//=====================================
		//For memory availability, I can neither precalculate the stellar spectra 
		//nor store population spectra for all threshold types simultaneously.  At the moment,
		//this recalculates the stellar spectrum for each threshold type.
		for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
		{
			//Loop over isochrone age
			printf ("%d %e  %d\n", i_IsoAge, Isochrones[i_IsoAge][0].t, NMasses[i_IsoAge]);
			
			for (i_IMF = 0; i_IMF < NIMF; i_IMF++)
			{
				dNdMIso[i_IMF][i_IsoAge] = calloc(NMasses[i_IsoAge], sizeof(double));
				dMPopdMIso[i_IMF][i_IsoAge] = 0.;
			}
				
			IndividualStellarSpectrum[i_IsoAge] = calloc (NMasses[i_IsoAge], sizeof(float*));
			IndividualStellarMagnitudes[i_IsoAge] = calloc (NMasses[i_IsoAge], sizeof(float*));
			IncludeThisStar[i_IsoAge] = calloc(NMasses[i_IsoAge], sizeof(int));
			
				
			for (i_ThreshType = 0; i_ThreshType < NThresh; i_ThreshType++) i_ThreshL_Cut[i_ThreshType][i_IsoAge] = calloc (NMasses[i_IsoAge], sizeof(int));
			
			for (i_MStar = 0; i_MStar < NMasses[i_IsoAge]; i_MStar++)
			{
				ThisStellarSpectrum = IndividualStellarSpectrum[i_IsoAge][i_MStar] = NULL;
				IncludeThisStar[i_IsoAge][i_MStar] = FALSE; //Star can be out of mass range.
				
				MStar_0 = Isochrones[i_IsoAge][i_MStar].M0;
				if (i_MStar == 0) {MLo = Isochrones[i_IsoAge][0].M0;} else {MLo = 0.5 * (Isochrones[i_IsoAge][i_MStar - 1].M0 + Isochrones[i_IsoAge][i_MStar].M0);}
				if (i_MStar == NMasses[i_IsoAge] - 1) {MHi = Isochrones[i_IsoAge][i_MStar].M0;} else {MHi = 0.5 * (Isochrones[i_IsoAge][i_MStar + 1].M0 + Isochrones[i_IsoAge][i_MStar].M0);}
				if ((MHi <= MStarMin) || (MLo >= MStarMax)) continue;
				if (MLo < MStarMin) MLo = MStarMin;
				if (MHi > MStarMax) MHi = MStarMax;
				if (MHi == MLo) continue;
				
				ThisStar = &(Isochrones[i_IsoAge][i_MStar]);
				TStar = (*ThisStar).Teff;
				AStar = (*ThisStar).A;
				LStar = (*ThisStar).L;
				MStar = (*ThisStar).M;
				StageStar = (*ThisStar).Stage;
				
				IncludeThisStar[i_IsoAge][i_MStar] = TRUE;
				ThisStellarSpectrum = IndividualStellarSpectrum[i_IsoAge][i_MStar] = calloc (NLam, sizeof(float));
				IndividualStellarMagnitudes[i_IsoAge][i_MStar] = calloc (NFilt, sizeof(float));
				
				//This is where the cut on luminosity (or whatever) happens.
				for (i_ThreshType = 0; i_ThreshType < NThresh; i_ThreshType++)
				{
					if (i_ThreshType == 0) i_ThisCut = floor((log10(LStar) - logThreshMin[i_ThreshType]) / dlogThresh[i_ThreshType]);
					if (i_ThreshType == 1) i_ThisCut = floor((log10(LStar / MStar_0) - logThreshMin[i_ThreshType]) / dlogThresh[i_ThreshType]);
					if (i_ThreshType == 2) i_ThisCut = floor((log10(MStar_0) - logThreshMin[i_ThreshType]) / dlogThresh[i_ThreshType]);
					//Sunscreen also has the ability to calculate the emission from stars in each stage
					//of their lives -- protostars, MS, subgiants, RGB, etc.
					//This is because the CMD isochrones have this type given by an integer within them.
					//Instead of dividing stars into those before and during a particular stage and those after,
					//Sunscreen outputs the luminosity from stars of that type and stars not of that type.
					//It's treated as a special case using the same variables as the threshold cuts.
					if (i_ThreshType == stellartype) i_ThisCut = StageStar;
					

					if (i_ThisCut < 0) i_ThisCut = 0;
					if (i_ThisCut > NL[i_ThreshType] - 1) i_ThisCut = NL[i_ThreshType] - 1;

					i_ThreshL_Cut[i_ThreshType][i_IsoAge][i_MStar] = i_ThisCut;
				}

				//dNIso normalized to 1 M_sun all with this age.
				//dNdMIso = d^2 N/[dM_iso dM*] * (MHi - MLo), the number of stars in the stellar
				//mass interval (MLo, MHi) per 1 M_sun of stars with this age.
				//dMPopdMIso = Int d^2 N/[dM_iso dM*] * dM*, the amount of stellar mass
				//per 1 M_sun of stars with this age.
				for (i_IMF = 0; i_IMF < NIMF; i_IMF++)
				{
					dNdMIso[i_IMF][i_IsoAge][i_MStar] = IMFNumberIntegral (MLo, MHi, IMFNorm[i_IMF], i_IMF);
					dMPopdMIso[i_IMF][i_IsoAge] += MStar * dNdMIso[i_IMF][i_IsoAge][i_MStar];
				}
				
				//Calculation of individual stellar spectrum
				//==========================================					
				//Interpolation of stellar model library
				//--------------------------------------
				if ((TStar < TGrid[0]) || (TStar > TGrid[nTForZ - 1]))
				{
					i_TStar = -1;
					igCold = igHot = -1;
					fTCold = fTHot = fgColdLo = fgColdHi = fgHotLo = fgHotHi = 0.0;
				}
				else
				{
					loggStar = (*ThisStar).logg;
					for (i_TStar = 0; i_TStar < nTForZ - 1; i_TStar++) {if (TGrid[i_TStar + 1] > TStar){break;}}
					
					fTCold = (TStar - TGrid[i_TStar]) / (TGrid[i_TStar + 1] - TGrid[i_TStar]);
					fTHot = (TGrid[i_TStar + 1] - TStar) / (TGrid[i_TStar + 1] - TGrid[i_TStar]);
					
					Ng = ngAtT[i_TStar];
					if (loggStar <= loggGrid[i_TStar][0]) 
					{
						//These need to be flipped from their "real" ordering.
						//(Because when doing the linear interpolation, f(x_Hi)
						//is multiplied by the "Lo" coefficient (x - x_Lo)/(x_Hi - x_Lo),
						//and vice-versa.)
						
						//This way handles the case where Ng = 1 for i_TStar + 1.
						igCold = -1;
						fgColdLo = 1.0;
						fgColdHi = 0.0;
					}
					else if (loggStar >= loggGrid[i_TStar][Ng - 1])
					{
						igCold = Ng - 2;
						fgColdLo = 1.0; 
						fgColdHi = 0.0;	
					}
					else
					{
						for (igCold = 0; igCold < Ng - 1; igCold++) {if (loggGrid[i_TStar][igCold + 1] > loggStar) {break;}}
						
						fgColdLo = (loggStar - loggGrid[i_TStar][igCold]) / (loggGrid[i_TStar][igCold + 1] - loggGrid[i_TStar][igCold]);
						fgColdHi = (loggGrid[i_TStar][igCold + 1] - loggStar) / (loggGrid[i_TStar][igCold + 1] - loggGrid[i_TStar][igCold]);
					}	

					Ng = ngAtT[i_TStar + 1];
					if (loggStar <= loggGrid[i_TStar + 1][0]) 
					{
						//These need to be flipped from their "real" ordering.
						//This way handles the case where Ng = 1 for i_TStar + 1.
						igHot = -1;
						fgHotLo = 1.0;
						fgHotHi = 0.0;
					}
					else if (loggStar >= loggGrid[i_TStar + 1][Ng - 1])
					{
						igHot = Ng - 2;
						fgHotLo = 1.0; 
						fgHotHi = 0.0; 	
					}
					else
					{
						//g[ig] <= gStar < g[ig + 1]
						for (igHot = 0; igHot < Ng - 1; igHot++) {if (loggGrid[i_TStar + 1][igHot + 1] > loggStar) {break;}}
						
						fgHotLo = (loggStar - loggGrid[i_TStar + 1][igHot]) / (loggGrid[i_TStar + 1][igHot + 1] - loggGrid[i_TStar + 1][igHot]);
						fgHotHi = (loggGrid[i_TStar + 1][igHot + 1] - loggStar) / (loggGrid[i_TStar + 1][igHot + 1] - loggGrid[i_TStar + 1][igHot]);	
						
					}
				}
				
				f00 = fTHot * fgColdHi;
				f01 = fTHot * fgColdLo;
				f10 = fTCold * fgHotHi;
				f11 = fTCold * fgHotLo;

				//Spectrum calculated from interpolation
				//--------------------------------------
				if (i_TStar >= 0)
				{
					//Much of Sunscreen's time is spent on the inner loop over wavelength when calculating a star's spectrum
					//This saves on pointer arithmetic
					if (igCold >= 0) ModelSpectrumPointer_00 = ModelSpectrum[i_TStar][igCold]; else ModelSpectrumPointer_00 = NULL;
					if (igCold >= -1) ModelSpectrumPointer_01 = ModelSpectrum[i_TStar][igCold + 1]; else ModelSpectrumPointer_01 = NULL;
					if (igHot >= 0) ModelSpectrumPointer_10 = ModelSpectrum[i_TStar + 1][igHot]; else ModelSpectrumPointer_10 = NULL;
					if (igHot >= -1) ModelSpectrumPointer_11 = ModelSpectrum[i_TStar + 1][igHot + 1]; else ModelSpectrumPointer_11 = NULL;
				}
				
				for (i_lambda = 0; i_lambda < NLam; i_lambda++)
				{
					if (i_TStar < 0)
					{
						ThisStellarSpectrum[i_lambda] = FnuPlanck (lamList[i_lambda], (float) TStar) * AStar;
					}
					else
					{
						//calloc on IndividualStellarSpectrum should now take care of this.
						//ThisStellarSpectrum[i_lambda] = 0.0;
						//Can't assume that any one of these will be executed, so have to pre-initialize ThisStellarSpectrum[i_lambda] to 0.
						if (igCold >= 0) ThisStellarSpectrum[i_lambda] += f00 * ModelSpectrumPointer_00[i_lambda];
						if (igCold >= -1) ThisStellarSpectrum[i_lambda] += f01 * ModelSpectrumPointer_01[i_lambda];					
						if (igHot >= 0) ThisStellarSpectrum[i_lambda] += f10 * ModelSpectrumPointer_10[i_lambda];					
						if (igHot >= -1) ThisStellarSpectrum[i_lambda] += f11 * ModelSpectrumPointer_11[i_lambda];

						ThisStellarSpectrum[i_lambda] *= AStar;
						
						if (ThisStellarSpectrum[i_lambda] > 1e99) printf ("%d %d  %d  %e %e %e %e  %e  %e\n", i_IsoAge, i_MStar, f00, f01, f10, f11, AStar, ThisStellarSpectrum[i_lambda]);
					}
				}
			}
		}


		//Inner loops
		//===========
		//Integration into single stellar population SED						
		//----------------------------------------------
		for (i_IMF = 0; i_IMF < NIMF; i_IMF++)
		{
			//Stellar population mass calculation
			//-----------------------------------
			for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
			{
				for (i_GalaxyAge = NtPop - 1; i_GalaxyAge < NtPop; i_GalaxyAge++) //0; i_GalaxyAge < NtPop; i_GalaxyAge++)
				{	
					logt = logtGalaxyAge[i_GalaxyAge];
								
					MPop = 0.0;
					for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++) MPop += dMPopdMIso[i_IMF][i_IsoAge] * dMIso[i_SFH][i_GalaxyAge][i_IsoAge];

					sprintf (MassFileName, "StellarPopulations/PopulationMass.%s.%s.logZ%+.2f.logt+%.2f.txt", IMFName[i_IMF], SFModeName[i_SFH], logZ[i_Z], logt);
					MassFile = fopen (MassFileName, "w");
					fprintf (MassFile, "%e\n", MPop);
					fclose (MassFile);
				}
			}

			//printf ("Unscreened stellar spectra.");
	
			//Unscreened single-age population calculation
			//--------------------------------------------
			for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
			{	
		
				ThisIsochroneSpectrum = IsochroneUnscreenedSpectrum[i_IsoAge];
				for (i_lambda = 0; i_lambda < NLam; i_lambda++) ThisIsochroneSpectrum[i_lambda] = 0.;
					
				for (i_MStar = 0; i_MStar < NMasses[i_IsoAge]; i_MStar++)
				{
					if (IncludeThisStar[i_IsoAge][i_MStar] == FALSE) continue;
					ThisStellarSpectrum = IndividualStellarSpectrum[i_IsoAge][i_MStar];
					ThisdNdMIso = (float) (dNdMIso[i_IMF][i_IsoAge][i_MStar]);
					for (i_lambda = 0; i_lambda < NLam; i_lambda++) ThisIsochroneSpectrum[i_lambda] += ThisdNdMIso * ThisStellarSpectrum[i_lambda];
				}
			}
			
			for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
			{
				//Loop over star formation histories.
				//for (i_GalaxyAge = NtPop - 1; i_GalaxyAge < NtPop; i_GalaxyAge++) //0; i_GalaxyAge < NtPop; i_GalaxyAge++)
				//i_GalaxyAge = NtPop - 2;
				for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
				{
					//Loop over galaxy total age.
					t = tGalaxyAge[i_GalaxyAge];
					logt = logtGalaxyAge[i_GalaxyAge];
					
					for (i_lambda = 0; i_lambda < NLam; i_lambda++) IntegratedUnscreenedSpectrum[i_lambda] = 0.;
					
					for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
					{		
						//Loop over isochrone age					
						ThisdMIso = dMIso[i_SFH][i_GalaxyAge][i_IsoAge];
							
						if (ThisdMIso == 0.0) continue;
						
						ThisIsochroneSpectrum = IsochroneUnscreenedSpectrum[i_IsoAge];
						
						for (i_lambda = 0; i_lambda < NLam; i_lambda++) IntegratedUnscreenedSpectrum[i_lambda] += ThisIsochroneSpectrum[i_lambda] * ThisdMIso;
					}
	
					sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt);			
					sprintf (UnscreenedSpectrumName, "UnscreenedResults/logZ%+.2f/%s/Unscreened_Spectrum.%s.txt", logZ[i_Z], IMFName[i_IMF], ModelIDString);
					
					printf ("\t%s\n", ModelIDString);
							
					UnscreenedSpectrumFile = fopen (UnscreenedSpectrumName, "w");
					for (i_lambda = 0; i_lambda < NLam; i_lambda++) fprintf (UnscreenedSpectrumFile, "%e %e\n", lamList[i_lambda], IntegratedUnscreenedSpectrum[i_lambda]);
					fclose (UnscreenedSpectrumFile);
					
					if (CalculateMagnitudes == FALSE) continue; 
							
					//Calculation of stellar population magnitudes
					//============================================
					for (i_Redshift = 0; i_Redshift < NRedshiftOutput; i_Redshift++)
					{
						sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f.%s", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt, zOutputName[i_Redshift]);							
						sprintf (UnscreenedMagnitudeName, "UnscreenedResults/logZ%+.2f/%s/Unscreened_Magnitudes.%s.txt", logZ[i_Z], IMFName[i_IMF], ModelIDString);
						UnscreenedMagnitudeFile = fopen (UnscreenedMagnitudeName, "w");	
						PrintMagnitudeHeader (UnscreenedMagnitudeFile, FilterName);
						
						if (i_Redshift == rest) i_ObsSpectrum = 0; else i_ObsSpectrum = i_GalaxyAge + 1;
						aScale = Redshift_GalaxyAge[i_ObsSpectrum] + 1.;
						ObsLambdaListPointer = ObsSpectrumLambda[i_ObsSpectrum];
					
						//I'm assuming D_L is used when calculating magnitude.
						for (i_lambda = 0; i_lambda < NLam; i_lambda++) ObsSpectrumLnu_Unscreened[i_lambda] = IntegratedUnscreenedSpectrum[i_lambda] * aScale;

						fprintf (UnscreenedMagnitudeFile, "%10s  %f  ", "All", MBolObs(ObsLambdaListPointer, ObsSpectrumLnu_Unscreened));
						
						for (i_Filter = 0; i_Filter < NFilt; i_Filter++)
						{
							FilterKernelPointer = FilterKernel[i_Filter][i_ObsSpectrum];

							RatioSum_Unscreened = 0.0;														
							for (i_lambda = 0; i_lambda < NLam; i_lambda++)
							{
								//Int dN/(dt dnu) * T(nu) dnu = Int dE/(dt dnu) * 1/(h nu) * T(nu) dnu
								//                            = Int dE/(dt dnu) * 1/(h lambda) * T(lambda) * dlambda
								//since dln nu = dln lambda		
								ThisFilterKernel = FilterKernelPointer[i_lambda];
								RatioSum_Unscreened += ObsSpectrumLnu_Unscreened[i_lambda] * ThisFilterKernel;
							}
													
							fprintf (UnscreenedMagnitudeFile, "%f ", -2.5 * log10(RatioSum_Unscreened));
						}
						fprintf (UnscreenedMagnitudeFile, "\n");
						fclose (UnscreenedMagnitudeFile);
					}
				
				}
			}

			for (i_ThreshType = 0; i_ThreshType < NThresh; i_ThreshType++)
			{		
				//Magnitudes for every threshold value of a given threshold type are
				//grouped into one file.
				//For example, there's a file containing magnitudes for log L = -3.50, -3.40, ...
				//But since single-age stellar population spectra are calculated for only one
				//threshold value at a time, I need to keep many magnitude files open.  Otherwise
				//I'd have to re-open and close them each time I wrote a line to them, drastically
				//slowing down the program.
				if (CalculateMagnitudes)
				{
					for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
					{
						for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
						{
							logt = logtGalaxyAge[i_GalaxyAge];
							
							for (i_Redshift = 0; i_Redshift < NRedshiftOutput; i_Redshift++)
							{
								if (i_ThreshType == stellartype)
								{
									sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f.%s", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt, zOutputName[i_Redshift]);
									sprintf (DimScreenedMagnitudeName, "Only_%s_Unscreened_Results/logZ%+.2f/%s/Only_%s_Unscreened_Magnitudes.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
									sprintf (BrightScreenedMagnitudeName, "AllExcept_%s_Unscreened_Results/logZ%+.2f/%s/AllExcept_%s_Unscreened_Magnitudes.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);								
								}
								else
								{
									sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f.%s", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt, zOutputName[i_Redshift]);							
									sprintf (DimScreenedMagnitudeName, "Small_%s_Screened_Results/logZ%+.2f/%s/Small_%s_Screened_Magnitudes.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
									sprintf (BrightScreenedMagnitudeName, "Large_%s_Screened_Results/logZ%+.2f/%s/Large_%s_Screened_Magnitudes.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
								}
								
								DimScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift] = fopen (DimScreenedMagnitudeName, "w");
								BrightScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift] = fopen (BrightScreenedMagnitudeName, "w");
								
								PrintMagnitudeHeader (DimScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift], FilterName);
								PrintMagnitudeHeader (BrightScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift], FilterName);
							}
						}
					}
				}
		
		
		
				for (i_ThreshL = 0; i_ThreshL < NL[i_ThreshType]; i_ThreshL++)
				{ 
					//Screened single-age population spectrum calculation
					//---------------------------------------------------
					//Go through isochrone stellar spectrum library
					for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
					{
						for (i_lambda = 0; i_lambda < NLam; i_lambda++) IsochroneScreenedSpectrum_NoDim[i_IsoAge][i_lambda] = IsochroneScreenedSpectrum_NoBright[i_IsoAge][i_lambda] = 0.;
						
						IsochroneLDyson_NoDim[i_IsoAge] = IsochroneLDyson_NoBright[i_IsoAge] = 0.;
											
						for (i_MStar = 0; i_MStar < NMasses[i_IsoAge]; i_MStar++)
						{
							if (IncludeThisStar[i_IsoAge][i_MStar] == FALSE) continue;
							
							//In the loop over spectrum wavelength, a lot of time would be spent on
							//pointer arithmetic; storing it in an auxiliary variable could speed things
							//up.
							ThisdNdMIso = (float) (dNdMIso[i_IMF][i_IsoAge][i_MStar]);

							//For each type of cut and each cut value, Sunscreen divides stars into those that
							//fall on one side and those that fall on the other side.  For example, if the cut
							//is luminosity, one side is those below L_cut and the other is those above L_cut.
							//ScreenedSpectrumPointer determines which spectrum (side 1 or side 2) this star's
							//luminosity is added to.

							if (i_ThreshType == stellartype)
							{
								if (i_ThreshL == i_ThreshL_Cut[i_ThreshType][i_IsoAge][i_MStar])
								{
									ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoDim[i_IsoAge];
									IsochroneLDyson_NoBright[i_IsoAge] += ThisdNdMIso * Isochrones[i_IsoAge][i_MStar].L;								
								}
								else
								{
									ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoBright[i_IsoAge];
									IsochroneLDyson_NoDim[i_IsoAge] += ThisdNdMIso * Isochrones[i_IsoAge][i_MStar].L;
								}
							}
							else
							{
								if (i_ThreshL <= i_ThreshL_Cut[i_ThreshType][i_IsoAge][i_MStar])
								{
									//This condition means that the threshold [luminosity] is LOWER than the star [luminosity]
									//So the star is [brighter] than the cut.  Which means that the star IS visible if
									//stars [dimmer] than the threshold [luminosity] are screened.
									ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoDim[i_IsoAge];
									IsochroneLDyson_NoBright[i_IsoAge] += ThisdNdMIso * Isochrones[i_IsoAge][i_MStar].L;
								}
								else
								{
									ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoBright[i_IsoAge];
									IsochroneLDyson_NoDim[i_IsoAge] += ThisdNdMIso * Isochrones[i_IsoAge][i_MStar].L;
								}	
							}
							
							ThisStellarSpectrum = IndividualStellarSpectrum[i_IsoAge][i_MStar];

							for (i_lambda = 0; i_lambda < NLam; i_lambda++) ThisIsochroneSpectrum[i_lambda] += ThisdNdMIso * ThisStellarSpectrum[i_lambda];
						}
						
					}
					
					//Now a library of SSP screened spectra for this luminosity cut are available.

					//Full stellar population spectrum integration
					//============================================
					for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
					{
						//Loop over star formation histories.
					
						//for (i_GalaxyAge = NtPop - 1; i_GalaxyAge < NtPop; i_GalaxyAge++) //0; i_GalaxyAge < NtPop; i_GalaxyAge++)
						//i_GalaxyAge = NtPop - 2; //2;
						//for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)

						for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
						{
							//Loop over galaxy total age.
						
							t = tGalaxyAge[i_GalaxyAge];
							logt = logtGalaxyAge[i_GalaxyAge];
					
							for (i_lambda = 0; i_lambda < NLam; i_lambda++) IntegratedScreenedSpectrum_NoDim[i_lambda] = IntegratedScreenedSpectrum_NoBright[i_lambda] = 0.;
							IntegratedLDyson_NoDim = IntegratedLDyson_NoBright = 0.;						
						
							for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
							{		
								//Loop over isochrone age					
								ThisdMIso = dMIso[i_SFH][i_GalaxyAge][i_IsoAge];
							
								if (ThisdMIso == 0.0) continue;
							
								ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoDim[i_IsoAge];
								for (i_lambda = 0; i_lambda < NLam; i_lambda++) IntegratedScreenedSpectrum_NoDim[i_lambda] += ThisIsochroneSpectrum[i_lambda] * ThisdMIso;
								IntegratedLDyson_NoDim += IsochroneLDyson_NoDim[i_IsoAge] * ThisdMIso;	
								
								ThisIsochroneSpectrum = IsochroneScreenedSpectrum_NoBright[i_IsoAge];						
								for (i_lambda = 0; i_lambda < NLam; i_lambda++) IntegratedScreenedSpectrum_NoBright[i_lambda] += ThisIsochroneSpectrum[i_lambda] * ThisdMIso;
								IntegratedLDyson_NoBright += IsochroneLDyson_NoBright[i_IsoAge] * ThisdMIso;
					
							}

							//Spectra output to files
							//-----------------------
							//The spectra are given in the form
							//lambda [cm]		L_nu (stars) [erg/s/Hz]		L_nu (Dyson(T1)) [erg/s/Hz]		L_nu (Dyson(T2)) [erg/s/Hz]	...
							logThresh = logThreshMin[i_ThreshType] + (double) i_ThreshL * dlogThresh[i_ThreshType];
							if (i_ThreshType == stellartype)
							{
								sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f.%s_%s", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt, ThreshName[i_ThreshType], StellarStageNames[i_ThreshL]);
								sprintf (DimScreenedSpectrumName, "Only_%s_Unscreened_Results/logZ%+.2f/%s/Only_%s_Unscreened_Spectrum.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
								sprintf (BrightScreenedSpectrumName, "AllExcept_%s_Unscreened_Results/logZ%+.2f/%s/AllExcept_%s_Unscreened_Spectrum.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);								
							}
							else
							{
								sprintf (ModelIDString, "logZ%+.2f.%s.%s.logt+%.2f.%s%+.2f", logZ[i_Z], IMFName[i_IMF], SFModeName[i_SFH], logt, ThreshName[i_ThreshType], logThresh);
								sprintf (DimScreenedSpectrumName, "Small_%s_Screened_Results/logZ%+.2f/%s/Small_%s_Screened_Spectrum.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
								sprintf (BrightScreenedSpectrumName, "Large_%s_Screened_Results/logZ%+.2f/%s/Large_%s_Screened_Spectrum.%s.txt", ThreshName[i_ThreshType], logZ[i_Z], IMFName[i_IMF], ThreshName[i_ThreshType], ModelIDString);
							}
							
							printf ("\t%s\n", ModelIDString);
							
							DimScreenedSpectrumFile = fopen (DimScreenedSpectrumName, "w");
							BrightScreenedSpectrumFile = fopen (BrightScreenedSpectrumName, "w");
							for (i_DysonTemp = 0; i_DysonTemp < ND; i_DysonTemp++) 
							{
								ThisTDyson4 = TDysonList[i_DysonTemp];
								ThisTDyson4 *= ThisTDyson4;
								ThisTDyson4 *= ThisTDyson4;
								ADysonDimScreened[i_DysonTemp] = IntegratedLDyson_NoDim * LSun / (sigmaSB * ThisTDyson4);
								ADysonBrightScreened[i_DysonTemp] = IntegratedLDyson_NoBright * LSun / (sigmaSB * ThisTDyson4);
							}

							for (i_lambda = 0; i_lambda < NLam; i_lambda++) 
							{
								fprintf (DimScreenedSpectrumFile, "%e %e", lamList[i_lambda], IntegratedScreenedSpectrum_NoDim[i_lambda]);
								for (i_DysonTemp = 0; i_DysonTemp < ND; i_DysonTemp++) fprintf (DimScreenedSpectrumFile, " %e", ADysonDimScreened[i_DysonTemp] * BlackbodyFlux[i_DysonTemp][i_lambda]);
								fprintf (DimScreenedSpectrumFile, "\n");
								
								fprintf (BrightScreenedSpectrumFile, "%e %e", lamList[i_lambda], IntegratedScreenedSpectrum_NoBright[i_lambda]);
								for (i_DysonTemp = 0; i_DysonTemp < ND; i_DysonTemp++) fprintf (BrightScreenedSpectrumFile, " %e", ADysonBrightScreened[i_DysonTemp] * BlackbodyFlux[i_DysonTemp][i_lambda]);
								fprintf (BrightScreenedSpectrumFile, "\n");
							}
							fclose (DimScreenedSpectrumFile);
							fclose (BrightScreenedSpectrumFile);
							
							
							if (CalculateMagnitudes == FALSE) continue; 
							
							//Calculation of stellar population magnitudes
							//============================================
							//Doing the convolution with the integrated spectrum rather than the 
							//individual stellar or isochrone spectra:
							//	1.) Allows for the possibility of eventually adding nebular, dust, 
							//	    megastructure thermal, etc. emission.
							//	2.) Saves about 300 MB of memory, since otherwise the code would have
							//	    to store the magnitude of each star at several redshifts.
							//Magnitude calculation at rest for each star is still useful for making
							//HR diagrams.
							for (i_Redshift = 0; i_Redshift < NRedshiftOutput; i_Redshift++)
							{
								if (i_Redshift == rest) i_ObsSpectrum = 0; else i_ObsSpectrum = i_GalaxyAge + 1;
								aScale = Redshift_GalaxyAge[i_ObsSpectrum] + 1.;
								ObsLambdaListPointer = ObsSpectrumLambda[i_ObsSpectrum];
								//printf ("%d %e\n", i_Redshift, aScale);
							
								//I'm assuming D_L is used when calculating magnitude.
								for (i_lambda = 0; i_lambda < NLam; i_lambda++)
								{
									ObsSpectrumLnu_NoDim[i_lambda] = (double) IntegratedScreenedSpectrum_NoDim[i_lambda] * aScale;
									ObsSpectrumLnu_NoBright[i_lambda] = (double) IntegratedScreenedSpectrum_NoBright[i_lambda] * aScale;
								}
								
								MagnitudeFilePointer_NoDim = DimScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift];
								MagnitudeFilePointer_NoBright = BrightScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift];
								
								if (i_ThreshType == stellartype)
								{
									fprintf (MagnitudeFilePointer_NoDim, "%10s  %f  ", StellarStageNames[i_ThreshL], MBolObs(ObsLambdaListPointer, ObsSpectrumLnu_NoDim));
									fprintf (MagnitudeFilePointer_NoBright, "%10s  %f  ", StellarStageNames[i_ThreshL], MBolObs(ObsLambdaListPointer, ObsSpectrumLnu_NoBright));
								}
								else
								{
									fprintf (MagnitudeFilePointer_NoDim, "%+f  %f  ", logThresh, MBolObs(ObsLambdaListPointer, ObsSpectrumLnu_NoDim));
									fprintf (MagnitudeFilePointer_NoBright, "%+f  %f  ", logThresh, MBolObs(ObsLambdaListPointer, ObsSpectrumLnu_NoBright));
								}
								
																
								for (i_Filter = 0; i_Filter < NFilt; i_Filter++)
								{
									FilterKernelPointer = FilterKernel[i_Filter][i_ObsSpectrum];

									RatioSum_NoDim = RatioSum_NoBright = 0.0;									
									for (i_lambda = 0; i_lambda < NLam; i_lambda++)
									{
										//Int dN/(dt dnu) * T(nu) dnu = Int dE/(dt dnu) * 1/(h nu) * T(nu) dnu
										//                            = Int dE/(dt dnu) * 1/(h lambda) * T(lambda) * dlambda
										//since dln nu = dln lambda		
										ThisFilterKernel = FilterKernelPointer[i_lambda];

										RatioSum_NoDim += ObsSpectrumLnu_NoDim[i_lambda] * ThisFilterKernel;
										RatioSum_NoBright += ObsSpectrumLnu_NoBright[i_lambda] * ThisFilterKernel;
									}
									
									fprintf (MagnitudeFilePointer_NoDim, "%f ", -2.5 * log10(RatioSum_NoDim));
									fprintf (MagnitudeFilePointer_NoBright, "%f ", -2.5 * log10(RatioSum_NoBright));
								}

								fprintf (MagnitudeFilePointer_NoDim, "\n");
								fprintf (MagnitudeFilePointer_NoBright, "\n");
								fflush (MagnitudeFilePointer_NoDim);
								fflush (MagnitudeFilePointer_NoBright);
							}	
							
						}
					}
				}
				
				if (CalculateMagnitudes)
				{
					for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
					{
						for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++)
						{
							for (i_Redshift = 0; i_Redshift < NRedshiftOutput; i_Redshift++) {fclose (DimScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift]); fclose (BrightScreenedMagnitudeFile[i_SFH][i_GalaxyAge][i_Redshift]);}
						}
					}
				}
			}
		}
				
		//Deallocation
		//============
		ThisStellarSpectrum = ThisIsochroneSpectrum = NULL;
			
		for (i_IsoAge = 0; i_IsoAge < NtIso; i_IsoAge++)
		{
			for (i_MStar = 0; i_MStar < NMasses[i_IsoAge]; i_MStar++) {free (IndividualStellarSpectrum[i_IsoAge][i_MStar]);  free (IndividualStellarMagnitudes[i_IsoAge][i_MStar]);}
			for (i_ThreshType = 0; i_ThreshType < NThresh; i_ThreshType++) free (i_ThreshL_Cut[i_ThreshType][i_IsoAge]);
			free (IncludeThisStar[i_IsoAge]);
			for (i_IMF = 0; i_IMF < NIMF; i_IMF++) free(dNdMIso[i_IMF][i_IsoAge]);
		}

		for (i_TStar = 0; i_TStar < nTForZ; i_TStar++) {for (i_gStar = 0; i_gStar < ngAtT[i_TStar]; i_gStar++) {free (ModelSpectrum[i_TStar][i_gStar]);} free (loggGrid[i_TStar]); free(ModelSpectrum[i_TStar]);}
		free (ModelSpectrum);
		free (loggGrid);
		free (TGrid);
		free (ngAtT);
	}

	if (CalculateMagnitudes)
	{
		for (i_Filter = 0; i_Filter < NFilt; i_Filter++) {for (i_GalaxyAge = 0; i_GalaxyAge < NtPop + 1; i_GalaxyAge++) free(FilterKernel[i_Filter][i_GalaxyAge]);}
	}
	
	for (i_SFH = 0; i_SFH < NSFH; i_SFH++)
	{
		for (i_GalaxyAge = 0; i_GalaxyAge < NtPop; i_GalaxyAge++) free(dMIso[i_SFH][i_GalaxyAge]);
		free (dMIso[i_SFH]);
	}

	for (i_IMF = 0; i_IMF < NIMF; i_IMF++){free (dNdMIso[i_IMF]); free (dMPopdMIso[i_IMF]);}

	return 0;
}

int CountLinesInFile (FILE *f)
{
	int n = 0;
	char c;
	
	if (f == NULL) {return 0;}
	//http://stackoverflow.com/questions/4278845/what-is-the-easiest-way-to-count-the-newlines-in-an-ascii-file
    while ((c=fgetc(f)) != EOF) {if (c == '\n') {n++;}}
    rewind (f);

	return n;
} 

void ReadIsochrone (char *FileName, StellarData **Isochrone, int *NMasses)
{
	int i;
	FILE *InFile;
	int iTemp;
	double t, M0, M, logL, logTeff, logR, logg, logZ;
	int Stage;

	InFile = fopen (FileName, "r");	
	(*NMasses) = CountLinesInFile (InFile);
	//printf ("%s %d\n", FileName, *NAges);
	
	*Isochrone = calloc (*NMasses, sizeof(StellarData));
	
	for (i = 0; i < (*NMasses); i++)
	{
		fscanf (InFile, "%le %lf  %lf %lf %lf %lf %lf %d  %lf\n", &t, &M0, &M, &logL, &logTeff, &logR, &logg, &Stage, &logZ);
		//fscanf (InFile, "%le %lf  %lf %lf %lf %lf %lf %d\n", &t, &M0, &M, &logL, &logTeff, &logR, &logg, &Stage);
		(*Isochrone)[i].t = t;
		(*Isochrone)[i].M0 = M0;
		(*Isochrone)[i].M = M;
		(*Isochrone)[i].Teff = pow(10.0, logTeff);
		(*Isochrone)[i].R = pow(10.0, logR);
		(*Isochrone)[i].logg = logg;
		(*Isochrone)[i].A = 4.0 * pi * pow(10.0, 2.0 * logR);
		//The log L = -9.999 entries actually make up ~10% of the luminosity
		//of older galaxies, because of their relatively high actual luminosities
		//and a plateau of ~1.90 M_sun stars at this log L at 1.74 Gyr (10^9.24 yr)
		(*Isochrone)[i].L = sigmaSB * pow(10.0, 4.0 * logTeff) * (*Isochrone)[i].A / LSun;
		(*Isochrone)[i].Stage = Stage;
		(*Isochrone)[i].Z = pow(10.0, logZ);
	}
	fclose (InFile);
}

void ReadLamList (char *FileName, float *Spectrum)
{
	int i;
	FILE *InFile;
	
	InFile = fopen (FileName, "r");
	for (i = 0; i < NLam; i++) {fscanf (InFile, "%e\n", &(Spectrum[i]));} // printf ("%d %e\n", i, Spectrum[i]);}
	fclose (InFile);
}

void ReadSpectrum (char *FileName, float *Spectrum)
{
	int i;
	float lam;
	FILE *InFile;
	
	InFile = fopen (FileName, "r");
	for (i = 0; i < NLam; i++) fscanf (InFile, "%e %e\n", &lam, &(Spectrum[i]));
	fclose (InFile);
}

void ReadSFH (char *FileName, int *NtSFH, double **t, double **SFH)
{
	int i;
	FILE *InFile;

	InFile = fopen (FileName, "r");	
	//printf ("%s\n", FileName);
	(*NtSFH) = CountLinesInFile (InFile);
	//printf ("%d\n", *NtSFH);
	*t = calloc (*NtSFH, sizeof(double));
	*SFH = calloc (*NtSFH, sizeof(double));
	
	for (i = 0; i < (*NtSFH); i++) fscanf (InFile, "%lf %lf\n", &((*t)[i]), &((*SFH)[i]));
	//for (i = 0; i < (*NtSFH); i++) printf ("%d %e %e\n", i, (*t)[i], (*SFH)[i]);
	fclose (InFile);
}

void ReadFilterN0 (char *FileName, int *NFLam, double **Lambda, double **N0)
{
	int i;
	FILE *InFile;

	InFile = fopen (FileName, "r");
	
	*NFLam = CountLinesInFile(InFile);
	*Lambda = calloc (*NFLam, sizeof(double));
	*N0 = calloc (*NFLam, sizeof(double));
	
	for (i = 0; i < *NFLam; i++) fscanf (InFile, "%le %le\n", &((*Lambda)[i]), &((*N0)[i])); 
	fclose (InFile);
}

//SFH returns the total mass of stars formed between tBeginInt and tEndInt

double BurstSFH (double tBeginInt, double tEndInt, double tBurstStart, double tBurstEnd, double MStarBurst)
{
	//MStarBurst is in M_sun
	if (tBurstEnd <= tBeginInt) {return 0.0;}
	if (tEndInt <= tBurstStart) {return 0.0;}
	if ((tBeginInt <= tBurstStart) && (tBurstEnd <= tEndInt)) {return MStarBurst;}
	if ((tBeginInt <= tBurstStart) && (tEndInt < tBurstEnd)) {return MStarBurst * (tEndInt - tBurstStart) / (tBurstEnd - tBurstStart);}
	if ((tBurstStart < tBeginInt) && (tBurstEnd <= tEndInt)) {return MStarBurst * (tBurstEnd - tBeginInt) / (tBurstEnd - tBurstStart);}
	if ((tBurstStart < tBeginInt) && (tEndInt < tBurstEnd)) {return MStarBurst * (tEndInt - tBeginInt) / (tBurstEnd - tBurstStart);}
	return 0.0;
}

double CtsSFH (double tBeginInt, double tEndInt, double tSFStart, double SFR)
{
	//Returns the mass of stars formed during tBeginInt <= t <= tEndInt, where t = 0 at the present.	
	//|tBeginInt| > |tEndInt|, assuming they are both less than 0.
	//SFR is in M_sun/yr
	if (tBeginInt >= tSFStart) {return SFR * (tEndInt - tBeginInt);}
	if (tEndInt <= tSFStart) {return 0.0;}
	return SFR * (tEndInt - tSFStart);
}

double TableSFH (double tBeginInt, double tEndInt, int NtSFH, double *tTable, double *SFRTable)
{
	int i;
	int iBeginGTE, iEndGTE;
	double SFRBegin, SFREnd;
	double dM = 0.0;

	//printf ("%e %e %d %e %e\n", tBeginInt, tEndInt, NtSFH, tTable[0], tTable[NtSFH - 1]);
	
	if (tEndInt < tTable[0]) return 0.0;
	if (tBeginInt < tTable[0]) tBeginInt = tTable[0];
	if (tBeginInt > tTable[NtSFH - 1]) return 0.0;
	if (tEndInt > tTable[NtSFH - 1]) tEndInt = tTable[NtSFH - 1];
	
	iBeginGTE = BinarySearchForFirstGTE (NtSFH, tTable, tBeginInt);
	iEndGTE = BinarySearchForFirstGTE (NtSFH, tTable, tEndInt);

	//printf ("%d %d\n", iBeginGTE, iEndGTE);
	
	if (iBeginGTE == iEndGTE) return 0.5 * (tEndInt - tBeginInt) * (SFRTable[iBeginGTE - 1] + SFRTable[iBeginGTE]);
	
	SFRBegin = ((tTable[iBeginGTE] - tBeginInt) * SFRTable[iBeginGTE - 1] + (tBeginInt - tTable[iBeginGTE - 1]) * SFRTable[iBeginGTE]) / (tTable[iBeginGTE] - tTable[iBeginGTE - 1]);
	SFREnd = ((tTable[iEndGTE] - tEndInt) * SFRTable[iEndGTE - 1] + (tEndInt - tTable[iEndGTE - 1]) * SFRTable[iEndGTE]) / (tTable[iEndGTE] - tTable[iEndGTE - 1]);

	dM = 0.5 * (tTable[iBeginGTE] - tBeginInt) * (SFRTable[iBeginGTE] + SFRBegin);
	for (i = iBeginGTE + 1; i < iEndGTE; i++) dM += 0.5 * (tTable[i] - tTable[i - 1]) * (SFRTable[i] + SFRTable[i + 1]);
	dM += 0.5 * (tEndInt - tTable[iEndGTE - 1]) * (SFREnd + SFRTable[iEndGTE - 1]);
	
	return dM;

}

double IMF (double M, double IMFNorm, int IMFVariation)
{
	//M is in units of M_sun
	if (IMFVariation == chabrier) return ChabrierIMF (M, IMFNorm);
	if (IMFVariation == salpeter) return SalpeterIMF (M, IMFNorm);
	if (IMFVariation == botheavy) return BotHeavyIMF (M, IMFNorm);
	if (IMFVariation == chabrier_minmass_1msun) {if (M < 1.0) return 0.0; else return ChabrierIMF (M, IMFNorm);}
	if (IMFVariation == chabrier_minmass_2msun) {if (M < 2.0) return 0.0; else return ChabrierIMF (M, IMFNorm);}
	if (IMFVariation == chabrier_minmass_5msun) {if (M < 5.0) return 0.0; else return ChabrierIMF (M, IMFNorm);}
	
	return 0.0;
}
	
double ChabrierIMF (double M, double IMFNorm)
{
	//Returns dN/dM in [M_sun^-1]
	//M is in M_sun.
	double mPrime;
	double ln10 = log(10.0);
	double lgMc = log10(0.079);
	
	//Chabrier dN/dm
	if (M <= 1.0)
	{
		mPrime = (log10(M) - lgMc) / 0.69;
		return IMFNorm * 0.158 * exp(-mPrime * mPrime / 2.0) / ln10 / M;
	}
	return IMFNorm * 0.0443 * pow(M, -1.3) / ln10 / M;
}

double SalpeterIMF (double M, double IMFNorm)
{
	//Returns dN/dM in [M_sun^-1] for a pure Salpeter power-law IMF M^-2.35
	//M is in M_sun.
	return IMFNorm * pow(M, -2.35);
}

double BotHeavyIMF (double M, double IMFNorm)
{
	//Returns dN/dM in [M_sun^-1] for a pure power-law IMF M^-3
	//M is in M_sun.
	double alpha = 3.0;
	return IMFNorm * pow(M, -alpha);
}

double IMFMassIntegral (double M1, double M2, double IMFNorm, int IMFVariation)
{
	//Returns integral of m dN/dm from M1 to M2
	double lnM1 = log(M1);
	double lnM2 = log(M2);
	double dlnM = (lnM2 - lnM1) / 100.0;
	double lnM;
	double M;
	double MSum = 0.0;
	
	for (lnM = lnM1 + dlnM; lnM < lnM2 - 0.5 * dlnM; lnM += dlnM)
	{
		//Int m dN/dm dm = Int m^2 dN/dm dlnm
		M = exp(lnM);
		MSum += M * M * IMF(M, IMFNorm, IMFVariation) * dlnM;
	}
	MSum += 0.5 * M1 * M1 * IMF(M1, IMFNorm, IMFVariation) * dlnM;
	MSum += 0.5 * M2 * M2 * IMF(M2, IMFNorm, IMFVariation) * dlnM;
	return MSum;
}

double IMFNumberIntegral (double M1, double M2, double IMFNorm, int IMFVariation)
{
	//Returns integral of dN/dm from M1 to M2
	double lnM1 = log(M1);
	double lnM2 = log(M2);
	double dlnM = (lnM2 - lnM1) / 100.0;
	double lnM;
	double M;
	double MSum = 0.0;
	
	for (lnM = lnM1 + dlnM; lnM < lnM2 - 0.5 * dlnM; lnM += dlnM)
	{
		//Int m dN/dm dm = Int m^2 dN/dm dlnm
		M = exp(lnM);
		MSum += M * IMF(M, IMFNorm, IMFVariation) * dlnM;
	}
	MSum += 0.5 * M1 * IMF(M1, IMFNorm, IMFVariation) * dlnM;
	MSum += 0.5 * M2 * IMF(M2, IMFNorm, IMFVariation) * dlnM;
	return MSum;
}

float FnuPlanck (float lambda, float T) 
{
	//double nu = c / lambda;
	//double x = hPlanck * nu / (kB * T);
	//x = hPlanck c / (kB T lambda) = [hPlanck c / kB] * 1/(lambda T)
	float x = xPre_Planck / (lambda * T);
	
	if (x > 50.0) return 0.0;
	//if (x <= 1e-6) return 2.0 * pi * kB * T / (lambda * lambda);
	if (x <= 1e-6) return FnuPre_RJ * T / (lambda * lambda);
	
	//return 2.0 * pi * hPlanck * nu / (lambda * lambda) / (exp(x) - 1.0);
	return FnuPre_Planck / (lambda * lambda * lambda) / (exp(x) - 1.0);
}

int BinarySearchForFirstGTE (int N, double *Table, double x)
{
    int iLo, iHi, iMid;

    iLo = 0;
    iHi = N - 1;
    while (iHi - iLo > 1)
    {
       //This prevents overflows when iHi + iLo > 32767.
       iMid = iLo + (iHi - iLo) / 2;
       if (x <= Table[iMid]) {iHi = iMid;} else {iLo = iMid;}
    }
    return iHi;
}

double MBolObs (double *Lambda, double *Spectrum)
{
	//Does NOT include the "missing" luminosity in the Dyson spheres.
	int i;
	double LBol = 0.0;
	double nuLo, nuHi;
	
	for (i = 1; i < NLam; i++)
	{
		nuLo = c / Lambda[i];
		nuHi = c / Lambda[i - 1];
		LBol += 0.5 * (double) (Spectrum[i] + Spectrum[i - 1]) * (nuHi - nuLo);
	}
	return (float) (-2.5 * log10(LBol / LBol0));
}

double InterpolatedN0 (double lambda, int NFLam, double *FilterLambda, double *FilterN0)
{
	int i;
	double N0;

	if (lambda < FilterLambda[0]) return 0.0;
	if (lambda > FilterLambda[NFLam - 1]) return 0.0;
	
	i = BinarySearchForFirstGTE (NFLam, FilterLambda, lambda);

	N0 = ((lambda - FilterLambda[i - 1]) * FilterN0[i] + (FilterLambda[i] - lambda) * FilterN0[i - 1]) / (FilterLambda[i] - FilterLambda[i - 1]);
	
	return N0;
}

double RedshiftAtLookbackTime (double tL, double H0, double Omega_M, double T_CMB)
{
	//This assumes a flat Universe with Omega = 1.
	int NSteps = 10000;
	int i;
	double tH = 1.0 / H_0;
	double Omega_r, Omega_L;

	double dtL = tL / (double) NSteps;
	double z = 0.;
	double zP1_1, zP1_2, zP1_3, zP1_4;
	double Ez_1, Ez_2, Ez_3, Ez_4;
	double k_1, k_2, k_3, k_4;
		
	Omega_r = Omega_CMB (H0, T_CMB);
	Omega_L = 1.0 - Omega_M - Omega_r;
	
	for (i = 1; i <= NSteps; i++)
	{
		tL = (double) i * dtL;
		
		zP1_1 = 1. + z;
		Ez_1 = sqrt(Omega_L + Omega_M * zP1_1 * zP1_1 * zP1_1 + Omega_r * zP1_1 * zP1_1 * zP1_1 * zP1_1);
		k_1 = dtL * (Ez_1 * zP1_1 / tH);
		
		zP1_2 = 1. + (z + 0.5 * k_1);
		Ez_2 = sqrt(Omega_L + Omega_M * zP1_2 * zP1_2 * zP1_2 + Omega_r * zP1_2 * zP1_2 * zP1_2 * zP1_2);
		k_2 = dtL * (Ez_2 * zP1_2 / tH);
		
		zP1_3 = 1. + (z + 0.5 * k_2);
		Ez_3 = sqrt(Omega_L + Omega_M * zP1_3 * zP1_3 * zP1_3 + Omega_r * zP1_3 * zP1_3 * zP1_3 * zP1_3);
		k_3 = dtL * (Ez_3 * zP1_3 / tH);
		
		zP1_4 = 1. + (z + k_3);
		Ez_4 = sqrt(Omega_L + Omega_M * zP1_4 * zP1_4 * zP1_4 + Omega_r * zP1_4 * zP1_4 * zP1_4 * zP1_4);
		k_4 = dtL * (Ez_4 * zP1_4 / tH);

		z += (k_1 + 2. * k_2 + 2. * k_3 + k_4) / 6.;
	}
	return z;
}

double rho_c (double H0) {return 3.0 * H0 * H0 / (8.0 * pi * G_N);}

double Omega_CMB (double H0, double T)
{
	double u_CMB = (sigmaSB * T * T * T * T) * 4.0 / c;
	return u_CMB / (rho_c (H0) * c * c);
}

void PrintMagnitudeHeader (FILE *OutFile, char ColumnNames[NFilt][100])
{
	int j;
	fprintf (OutFile, "#log10(L) \tBol(Obs)");
	for (j = 0; j < NFilt; j++) fprintf (OutFile, " \t%s", ColumnNames[j]);
	fprintf (OutFile, "\n");
}
