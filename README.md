# Sunscreen
Stellar population synthesis for Kardashev III-hosting galaxies where only some stars are in Dyson spheres.

Sunscreen (v1.21)
=================

Brian Lacki

Breakthrough Listen

29 October 2018

---

# TABLE OF CONTENTS
* I. Introduction
* II. Setup
  * A. Requirements
  * B. Installing
  * C. Running
* III. Input
  * &#48;. Allowed metallicity values
  * A. Model stellar spectra
  * B. Isochrones
  * C. Star formation histories
  * D. Filter responses and photometry
  * E. Initial mass functions
  * F. Galaxy ages
* IV. Output
  * A. Output location
  * B. Spectra
    * 1. File names
    * 2. Format
  * C. Synthetic photometry
    * 1. Conventions
    * 2. Format
  * D. Stellar masses
* V. References and Credit
  * A. Stellar spectrum library
  * B. Isochrones
  * C. Star formation histories
  * D. Filter responses
  * E. Citing Sunscreen
  * F. Support and acknowledgements

# I. Introduction
Sunscreen is a stellar population synthesis program.  It is originally designed for work in SETI to calculate what certain kinds of unnatural stellar populations would look like.  Specifically, what would a galaxy look like if only stars less luminous than the Sun were cloaked in some way, such as being enclosed in a Dyson sphere?  While the program was written for a SETI project, it can also be used to gain insight into what stars contribute to each part of an SED.  

Sunscreen calculates the spectrum of a galaxy whose stars are screened (rendered invisible) according to a cut based on physical properties of the stars.  The available cuts are stellar luminosity, luminosity-to-mass ratio, mass, and evolutionary stage.  It then calculate synthetic photometry from these spectra, both in the rest frame and redshifted according to the galaxy's age.  If you want to know where a single age stellar population would appear on a color-color diagram if all of its stars brighter than 100 L<sub>sun</sub> were blotted out, or what the spectrum is of the AGB stars in a galaxy with constant star-formation, Sunscreen is your program.

Sunscreen has several limitations as described in Lacki (2018).  It assumes all the stars have the same metallicity.  It does not model binarity of stars, or stellar types due to binary interactions (e.g. extreme horizontal branch stars).  Nor does it take into account stellar rotation.  Sunscreen does not model dust extinction, neither interstellar nor circumstellar.  Finally, it only includes stellar emission, with an option for thermal emission from Dyson spheres.  There is no dust or nebular emission included in the output spectra or photometry.


# II. Setup

## A. Requirements

Sunscreen uses 120 MB of permanent storage (hard drive) when uncompressed with the default set of auxiliary files.

When calculating the default set of models, Sunscreen uses:
+ 1.1 GB memory
+ ~20 GB permanent storage (hard drive)

Running the default set of models takes of order 10 hours on my 2 GHz laptop.

## B. Installing

1. Download Sunscreen-Base-v1_2_1.zip and unzip it in the parent directory where you want Sunscreen to be placed.
2. Download Isochrones.zip and unzip it in Sunscreen/CMD-Isochrones.
3. Download StellarSpectra.zip and unzip it in Sunscreen/LejeuneSpectra.
4. Download PhotometryZeroPoints.zip and unzip it in Sunscreen/FilterResponses.
5. Compile Sunscreen.c.  The code for Sunscreen itself is written in C and requires no special libraries.  On a UNIX-like environment, it can be compiled simply:
```
gcc Sunscreen.c -o Sunscreen
```

Sunscreen requires the following output directory structure:
```
	Small_logL_Screened_Results
		logZ+0.00
			BottomHeavy
			Chabrier
			Salpeter
		logZ+0.30
			...
		logZ-1.00
			...
	Small_logLtoM0_Screened_Results
		...
	Small_logM0_Screened_Results
		...
	Large_logL_Screened_Results
		...
	Large_logLtoM0_Screened_Results
		...
	Large_logM0_Screened_Results
		...
	Only_Stage_Unscreened_Results
		...
	AllExcept_Stage_Unscreened_Results
		...
	UnscreenedResults
		...
	StellarPopulations
```
The subdirectories listed under `Small_logL_Screened_Results` and `Small_logL_Screened_Results/logZ+0.00` should be present in all these directories.  In the subdirectories of the form `logZ+0.00`, the quantity should be the base 10 logarithm of the metallicity relative to Solar ([M/H]). 

These output directories can be recreated in a UNIX-like environment by running the included CreateOutputDirectories.sh script.
In addition, input data files &ndash; namely, the star formation histories, stellar model spectra, stellar isochrones, and filter bandpasses &ndash; need to be included in the directories:
```
	SFRs
	LejeuneSpectra
	CMD-Isochrones
	FilterResponses
```

# C. Running
To run Sunscreen after install it, simply enter
```
./Sunscreen
```
or whatever file you named the executable after compiling it.  No command line parameters are used.


# III. Input

## 0. Allowed metallicity values

The three values of metallicity for which Sunscreen produces results are 0.1 Z<sub>sun</sub>, 1.0 Z<sub>sun</sub>, and 2.0 Z<sub>sun</sub>.

Changing them requires
1. Adding the subdirectories to the output directory structure (see the included `CreateOutputDirectories.sh` script).
2. Adding model stellar spectra for the desired metallicities to the `LejeuneSpectra` directory (see III-A).
3. Adding isochrones for the desired metallicities to the `CMD-Isochrones` directory (see III-B).
4. Changing the `NZ` constant to the desired number of metallicities to calculate results for, by looking for the `#define NZ` line in Sunscreen.c.
5. Changing the initial values of the `logZ` array in Sunscreen.c to the desired values.  

## A. Model stellar spectra

Sunscreen uses the BaSeL library of theoretical stellar spectra.  These are stored in the `LejeuneSpectra` directory.  The model spectra are for model stellar atmospheres of a specified metallicity, effective temperature, and surface gravity.

I have processed the data of these spectra into a form that Sunscreen can read.  Each spectrum has a name of the format:
```
LejeuneSpectrum.logZ%+.2f.Teff+%.0f.logg%+.2f
```
where metallicity (logZ) is given relative to Solar on a base-10 log scale (so Solar metallicity is +0.00, 1/10 Solar is -1.00), also denoted as [M/H]; temperature (Teff) is given without leading zeros in Kelvin; and surface gravity (logg) is given in cgs units, cm/s<sup>2</sup>.  The format of the file is a text file with two columns: wavelength in cm, and emergent specific flux from the stellar surface (F<sub>&nu;</sub>) in erg/cm<sup>2</sup>/s/Hz.  The included `ExtractLejeuneSpectra.sh` script in the `LejeuneSpectra` directory can be used to create these files from downloadable BaSeL data.

Note that **all** the spectra are assumed to use the same wavelengths in their spectra.  There are 1,221 wavelengths given in each spectrum, all identical from file to file, and they have been listed in `LejeuneSpectra/lambda-LejeuneSpectra.txt`.  Sunscreen reads this file, and uses a global constant `NLam` set to 1,221 in the `#define NLam` line.

The BaSeL models form a partial grid in temperature and surface gravity, but not one with every pair of T<sub>eff</sub> and g.  All values of surface temperature for which there are any models are stored in the `FullTGrid` array.  Likewise, all values of log<sub>10</sub> g for which any model exists is stored in `FullloggGrid`.  Sunscreen first sorts models by temperature, and then by surface gravity.  It checks each pair of T<sub>eff</sub> and g, and notes which temperatures have a model for the current metallicity.  These are stored in `TGrid`.  Then for each temperature, Sunscreen checks whether a model spectrum exists for each value of log<sub>10</sub> g, and if so that model will be stored in the `loggGrid`, an uneven 2D array.  The number of extant models (values of surface gravity) for each temperature at the current metallicity is stored in the `ngAtT array`.

When calculating the spectrum of a star of some temperature and surface gravity, Sunscreen generally uses 2D interpolation.  If its temperature falls outside the range spanned by `TGrid`, Sunscreen will use a blackbody spectrum.  If its surface gravity falls outside the range spanned by `loggGrid` at a relevant temperature, it will use the maximum or minimum surface gravity, whichever is closer to the actual value.


## B. Isochrones

Sunscreen uses the PARSEC v1.2S + COLIBRI PR16 isochrones available from the CMD 3.0 webserver.  These are stored in the CMD-Isochrones directory.  They are given in steps of log<sub>10</sub> t<sub>Iso</sub>, where t<sub>iso</sub> is the age of the stars on the isochrone.  Currently the isochrones span log<sub>10</sub> (t<sub>iso</sub> / yr) from 4.00 to 10.13.  The lower limit is somewhat arbitrarily chosen, but generally predates any star falling on the main sequence.  The upper limit is given by the maximum age available from the CMD webserver.  log<sub>10</sub> t<sub>iso</sub> values from 4.00 to 10.12 were downloaded in two steps (4.00 to 7.00, 7.00 to 10.12) using the isochrone grid option, while 10.13 was not available this way and was downloaded as a single isochrone of age 13.4896 Gyr.  The number of isochrone ages is given in the `NtIso` constant.

The bulk stellar data in these isochrones have been processed into a special text format that Sunscreen can read.  Each name has the form:
```
Isochrone.logZ%+.2f.logt+%.2f.txt
```
where metallicity (logZ) is given in the [M/H] form, a base-10 logarithm relative to Solar, and age (logt) is given in the base 10 logarithm of age in years. Within each file is a table with several columns.  The columns are, in order: 
1. age in years
2. initial stellar mass in M<sub>sun</sub> 
3. current stellar mass in M<sub>sun</sub>
4. base-10 logarithm of luminosity in L<sub>sun</sub>
5. base-10 logarithm of effective temperature of the star in Kelvin
6. base-10 logarithm of stellar radius in cm
7. base-10 logarithm of stellar surface gravity in cm/s<sup>2</sup>.
8. an integer specifying the approximate stage of stellar evolution the star is in, and
9. [M/H], given as a check that the processing has occured correctly.  
These files have been extracted from the `output*` files I downloaded from CMD using the ExtractIsochrones script in the CMD-Isochrones folder.  

Some things to be aware of:
* The ages in these files are taken directly from the `output*` files from the CMD webserver.  As such, they only have three significant digits &ndash; i.e., they are not calculated from the logt value in the file name.
* Solar metallicity is assumed to be 0.0152, the default value for the current iteration of CMD.
* Stellar radius is calculated from the star's surface gravity and current mass, under the star is spherical.
* Sunscreen actually disregards the value of log L given in these files to ensure physical consistency.  Instead, Sunscreen calculate L using the Steffan-Boltzmann law: L = 4 pi R<sup>2</sup> &sigma;<sub>SB</sub> T<sub>eff</sub><sup>4</sup>
* The CMD isochrones for large t have a terminating point that correspond to a post-AGB star.  The luminosity of these "stars" is set to log L = -9.999, but the other physical parameters indicate that they are bright and blue.  According to the CMD 3.0 FAQ, they are only included to signal the end of the isochrone to some codes.  I regard them as spurious, and they are not included when ExtractIsochrones is run.
* The CMD FAQ indicates that "stage" is only a rough indication. I've noticed that stars of mass around ~15 M<sub>sun</sub> seem to remain in the pre-MS phase for much longer than stars of much higher or lower stellar masses.

Each stage of stellar evolution is represented in these isochrones and has a short identifier given in the StellarStageNames file.  The stages are

| Value | Short name | Stage
| --- | --- | --- |
| 0 | PreMS | Pre-main sequence star (probably including some high mass MS stars) |
| 1 | MS | Main sequence star |
| 2 | SGB | Subgiant branch |
| 3 | RGB | Red giant branch |
| 4 | CEHB1 | Horizontal branch |
| 5 | CEHB2 | Horizontal branch |
| 6 | CEHB3 | Horizontal branch |
| 7 | EAGB | Early asymptotic giant branch phase |
| 8 | TPAGB | Thermally pulsing giant branch phase |
| 9 | Remnant | Post-AGB stars and compact stellar remnant |

PARSEC, and likewise, Sunscreen, does not include any post-AGB/remnant stars.

To extend (or contract) the range of isochrone ages:
1. Set `logtIsoMin` to the new minimum of log<sub>10</sub> (t<sub>iso</sub> / yr).
2. Set `logtIsoMax` to the new maximum of log<sub>10</sub> (t<sub>iso</sub> / yr).
3. Set `dlogtIso` to the new increment in log<sub>10</sub> t<sub>iso</sub>.  Be aware that the AGB boosting phase described in Girardi et al. (2013) is resolved only if this is less than ~0.05.
4. Change the `NtIso` constant to the new number of isochrones, which includes `logtIsoMin` and `logtIsoMax`, as given in the `#define NtIso` line.

Note that Sunscreen calculates the stellar spectrum of every star on every isochrone.  This accounts for the bulk of memory usage by Sunscreen, so adding many more isochrones will increase Sunscreen's memory usage.  However, Sunscreen only reads the isochrones for whichever metallicity value it is currently working on, and only retains spectra for stars of the current metallicity.

The `CalculateIsochronePhotometry.c` C program can be used to create additional files, `IsochroneMagnitudes.*`, containing synthetic photometry of each star on the isochrones.  These files are not directly used by Sunscreen, but they can be useful for populating color-magnitude diagrams.  They are text files with a header line followed by several columns of numbers.  The columns are: log<sub>10</sub> stellar age, initial stellar mass, log<sub>10</sub> stellar luminosity, absolute bolometric magnitude, and then absolute magnitude in various bands given in the header.  In order to run `CalculateIsochronePhotometry.c`, the content of FilterResponses need to be present (section III-D), as well as the `Isochrone.*` files.  


## C. Star formation histories

Sunscreen calculates the stellar population by integrating the star-formation rate (SFR) over small intervals into stellar mass.  The width of these intervals is the same as for the input isochrones, as given in the dlogtIso variable.  So if isochrones are given in increments of dlog<sub>10</sub> t = 0.01 (as in the default), Sunscreen integrates the mass within log<sub>10</sub> t intervals of width 0.01.  

When a star formation history (SFH) function is called by Sunscreen, two variables required are `tBeginInt` and `tEndInt`.  These are the start and end of the integration interval, with the present set to t = 0 and with t in chronological order.  These values are negative.  

Two star-formation histories (SFHs) are hardcoded into Sunscreen: a constant star-formation rate (SFR) and a burst.  The constant SFR is calculated by calling the `CtsSFH` function.  It has two free parameters: `tSFStart`, which gives the time the star formation started, and `SFR`, the star-formation rate since `tSFStart` in M<sub>sun</sub>/yr.  `tSFStart` is also measured with respect to the present in cosmic time, with t = 0 at the present, and takes on a negative value.  The burst SFH is calculated by calling the `BurstSFH` function.  It has three parameters: `tBurstStart`, `tBurstEnd`, and `MStarBurst`.  `tBurstStart` and `tBurstEnd` give the start and end times of the starburst (with present having t = 0, and these parameters being negative), while `MStarBurst` specifies the mass of the starburst in M<sub>sun</sub>.  Nothing actually requires the "burst" to be short, but it is currently called with a duration of 1 year to simulate an instantaneous burst.

After those two SFHs, Sunscreen uses SFHs stored in auxiliary data files, whose names are stored in the `SFHTableFileNames` array.  I've stored SFHs in the `SFRs` directory, but it is not necessary as long as the path is given in the `SFHTableFileNames`.  The format of the SFH table files is simply two columns: cosmic time since the Big Bang in years, and SFR in M<sub>sun</sub>/yr.  Note the convention on the time column is different than the ones used by Sunscreen: it measures time since the Big Bang, not the present, although in both cases t for older stars is strictly less than t for younger stars. 

The available auxiliary SFHs are:

```
B13-log10Mh+11.0
B13-log10Mh+12.0
B13-log10Mh+13.0
B13-log10Mh+14.0
```
These are average star formation histories from Behroozi et al. (2013) for present-day halo masses of 10<sup>11</sup>, 10<sup>12</sup>, 10<sup>13</sup>, and 10<sup>14</sup> M<sub>sun</sub>.  These correspond to the values plotted in the right panel of Figure 6, or those in the `sfh` directory of `release-sfh_z0_z8_052913` downloaded from Behroozi's website.  The `sfr` directory contain SFRs of galaxies where the halo mass is measured at the observed redshift, not the present day.  The difference is that a halo that has M<sub>h</sub> = 10<sup>12</sup> M<sub>sun</sub> at z = 5 is going to be a lot bigger at z = 0.  Those values (plotted in the left panel of their Figure 6) do not correspond to the SFRs over time of individual galaxies.

All of these models have nonzero star-formation at z = 0.  The first two can be considered blue star-forming galaxies, while the latter two are green valley galaxies where the quenching process has not completed.

The SFRs directory also includes a SFH, `B13-log10Mh+15.0`, for galaxies in present-day halo masses of 10<sup>15</sup> M<sub>sun</sub>, but it is not calculated by Sunscreen.

```
M15-MJAM+11.0
```
This is one of the average star formation histories from McDermid et al. (2015), Figure 16 top panel, for early type galaxies as derived from the ATLAS3D survey.  It is the SFH for log<sub>10</sub> (M<sub>JAM</sub> / M<sub>sun</sub>) between 11.0 and 11.5 (yellow-orange line), where M<sub>JAM</sub> is taken in that paper to be the stellar mass.  From what I can tell, the line in the plot is actually the rate that z = 0 stellar mass accumulates, not the total physical SFR.  I convert it to SFR using the prescription of Leitner & Kravtsov (2012) for a Chabrier IMF, which gives a present day stellar mass of 10<sup>11</sup> M<sub>sun</sub>.  
I use it as a model of a more realistic quenched/red galaxy than a single age burst.

```
W14-dIrr
```
The mean star formation history from Weisz et al. (2014) for Local Group dIrr galaxies.  It is calculated from the thick blue line shown in Figure 13, after differentiating, averaging on 1.0 Gyr timescales, and converting to physical star formation rate using Leitner & Kravtsov (2012).
The star-formation in this model is higher shortly after the Big Bang, low for the next few Gyr, and then rising again around z = 0.  It is a blue star-forming galaxy, but one with a very small stellar mass.

For these auxiliary SFHs, Sunscreen integrates the mass formed between `tBeginInt` and `tEndInt` by calling `TableSFH`, where `tBeginInt` and `tEndInt` are still measured with respect to the present being t = 0.  Three additional variables must be supplied: `NtSFH`, the number of data points in the table, a `tTable` array storing a list of the times given in the table, and a `SFRTable` array storing SFRs in M<sub>sun</sub>/yr at those times.  These will be handled automatically by the way Sunscreen reads in the files.

To add another SFH in table form:
1. Place a SFH of the correct format somewhere like the `SFRs` directory where Sunscreen can read it.
2. Change the `NSFH` variable to the correct number of SFHs, which is set in the `#define NSFH` line.
3. Add the name of the SFH to be used in files to the `SFModeName` array.
4. Add the path to the SFH file to the `SFHTableFileNames` array.

Note that Sunscreen always does the constant and burst SFRs first before moving on to the table SFHs.  So `NSFH` must be >= 3 for any of the Table SFHs to run.  To skip the constant and burst SFHs, change the initial values of `i_SFH` in the for loops to 2, but make sure `NSFH` still includes the constant and burst SFH.


## D. Filter responses and photometry

Sunscreen can calculate synthetic photometry of galaxies in many filters, given in Table 1 of Lacki (2018).

Each filter has a short name given in the `FilterName` array (e.g., "u_SDSS" for the u-band filter used by SDSS).  A kind of filter response curve is stored in files in the FilterResponses directory named:
```
PhotonZeroPoint-%s.txt
```
where `%s` represents this short name in the `FilterName` array.

The format of the `PhotonZeroPoint*` files is a text file with two columns.  The first column is wavelength in cm.  The second column parameterizes the filter response and the zero point of the magnitude system for that filter.  To be more precise, the second column is the **reciprocal** of the number flux of **monochromatic** photons at that wavelength **only** (in photons/cm<sup>2</sup>/s) from a source with zero magnitude in that band.  Because it is a reciprocal, if a filter does NOT respond to the photons of that wavelength, this value will be zero; if it is highly responsive, the value will be relatively large.  

Magnitude is calculated from a spectrum by convolving the **photon** flux spectrum of a galaxy with the second column of `PhotonZeroPoint`.

Note that Sunscreen does NOT include any additional correction for atmospheric extinction &ndash; these must be accounted for in the calculations of the zero magnitude fluxes when producing these files. 

To add a new filter for Sunscreen to calculate photometry in:
1. Place a `PhotonZeroPoint-*.txt` file corresponding to the filter response in the `FilterResponses` directory.
2. Add the name of the filter, the same as used in the name of the `PhotonZeroPoint` file, to the `FilterName` array.
3. Update the `NFilt` constant to the correct number of filters.  It can be found in the `#define NFilt` line.


## E. Initial mass functions

IMFs are hardcoded into Sunscreen.

The IMF is calculated by calling the IMF wrapper function.  It returns a quantity of the form dN/dM&ast;, the number of stars per infitesimal stellar mass range per unit solar mass of stellar population.  All these masses are in units of M<sub>sun</sub>.  Most of the time, Sunscreen doesn't call it directly, but calls IMFMassIntegral or IMFNumberIntegral to return the mass or number of stars between `MLo` and `MHi` per unit solar mass of stellar population.

Three IMF functions exists, but there are six IMFs that Sunscreen is set up to calculate.  The three IMF functional forms are a Chabrier IMF, a Salpeter IMF, and a bottom-heavy IMF of the form dN/dM<sub>&ast;</sub> ~ M<sub>&ast;</sub><sup>-3</sup>.  At present, these are the three IMFs that Sunscreen calculates.  The other three IMFs are Chabrier IMFs, but with all star formation below 1 M<sub>sun</sub>, 2 M<sub>sun</sub>, and 5 M<sub>sun</sub> suppressed.  

For the Chabrier IMF, the normalization is set such that it gives the SFH currently being calculated.  This normalization is calculated first.  The Salpeter and bottom-heavy IMFs are then normalized so that the value dN/dM<sub>&ast;</sub> (M<sub>&ast;</sub> = 1 M<sub>sun</sub>) is equal to that quantity's value for the Chabrier IMF using the precalculated normalization.  This normalization is chosen so that an old (~10 Gyr) stellar population will be roughly as bright for the Salpeter and bottom-heavy IMFs as for the Chabrier IMF, since most of the luminosity will be dominated by red giants of mass ~1 M<sub>sun</sub>.  This means the true SFR assumed for these IMFs is NOT the value given by the SFHs.  The Chabrier IMFs with mass cutoffs have the same normalization as the Chabrier IMF.

To add a new IMF:
1. Write a function that calculates the IMF (dN/dM<sub>&ast;</sub> in units of M<sub>sun</sub>, where N is the number of stars per unit solar mass of stellar population) and put it in Sunscreen.
2. Put a new variable that will index this IMF in the enum `IMFList`, in the order that you want the IMF to be calculated.  For example, change this:
```
enum IMFList {chabrier, salpeter, botheavy, chabrier_minmass_1msun, chabrier_minmass_2msun, chabrier_minmass_5msun};
```
to this:
```
enum IMFList {chabrier, salpeter, botheavy, new_imf};
```
3. In the `IMF (double, double, int)` function, write a line of the form:
```
if (IMFVariation == new_imf) return NewIMF (M, IMFNorm);
```
The "`new_imf`" should be the same name you put in `enum IMFList`.

4. Under the "Universal lookup tables" heading in Sunscreen's code, add a line that sets the new IMF's normalization and assigns it to its place in the `IMFNorm` array.  For example, the code for the bottom-heavy IMF is:
```
IMFNorm[botheavy] = IMFNorm[botheavy] * IMF (1.0, 1.0, chabrier) / IMF (1.0, 1.0, botheavy);
```
5. Add the name of the IMF to the `IMFName` array, using the same order as `IMFList`.
6. Update the `NIMF` constant to the number of IMFs that you want calculated, in the `#define NIMF` line.

Sunscreen assumes the stellar population lies strictly between 0.1 and 120 M<sub>sun</sub>.  These values are stored in `MStarMin` and `MStarMax`, respectively.  They are currently the same for all IMFs.  


## F. Output ages

Sunscreen calculates the spectrum of each model galaxy at several different ages: 1 Myr, 10 Myr, 100 Myr, 1 Gyr, 3.16 Gyr, 6.31 Gyr, 12.6 Gyr, and 13.5 Gyr.  The base 10 logarithm of these ages in years is stored in the `logtGalaxyAge` array.

For the constant SFH, the age is the time elapsed since star formation started.  For the burst SFH, the age is the time elapsed since the starburst ended.  For each of the other SFHs stored in data files, the age instead represents the age of the Universe, with the SFH assumed to start at the Big Bang.

To change the galaxy ages for which spectra are output:
1. Change the initial values of the `logtGalaxyAge` array to include the base 10 logarithms of the desired ages in years.  
2. Update `NtPop` to the number of ages in this array, in the `#define NtPop` line.

Remember that Sunscreen integrates stellar populations by stepping through each of the isochrones (section II-B).  Ages younger than 10<sup>`logtIsoMin`</sup> or older than 10<sup>`logtIsoMax`</sup> will not be calculated properly.  Sunscreen does NOT check this.


# IV. Output
Sunscreen calculates full grids of models: one spectrum and one set of magnitudes calculated for every combination of parameters.

## A. Output location
Models are stored in directories, sorted first according to the type of model, then metallicity, then IMF.  Results stored in each of the highest level directories are:
1. `Small_logL_Screened_Results`
  Results for stellar populations where all stars with (log) luminosities below a threshold log L are screened.
2. `Small_logLtoM0_Screened_Results`
  Results for stellar populations where all stars with (log) luminosity-to-initial-mass ratios below a threshold log (L/M<sub>0</sub>) are screened.
3. `Small_logM0_Screened_Results`
  Results for stellar populations where all stars with (log) initial mass below a threshold log M<sub>0</sub> are screened.
4. `Large_logL_Screened_Results`
  Results for stellar populations where all stars with (log) luminosities greater above a threshold log L are screened.
5. `Large_logLtoM0_Screened_Results`
  Results for stellar populations where all stars with (log) luminosity-to-initial-mass ratios above a threshold log (L/M<sub>0</sub>) are screened.
6. `Large_logM0_Screened_Results`
  Results for stellar populations where all stars with (log) initial mass above a threshold log M<sub>0</sub> are screened.
7. `Only_Stage_Unscreened_Results`
  Results for stellar populations where the only unscreened stars are those in a certain stage of their lives.  For example, here is where spectra and magnitudes for stellar populations with only the main sequence stars visible.
8. `AllExcept_Stage_Unscreened_Results`
  Results for stellar populations where all stars are unscreened, except those in a particular stage of their lives.  Here's where you would find spectra and magnitudes for stellar populations that exclude all main sequence stars, but include stars in all other stages of their lives.
9. `UnscreenedResults`
  Results for stellar populations with no stars screened at all (i.e., a natural galaxy without extinction).
10. `StellarPopulations`
  Results for bulk parameters for the stellar population are stored here.  At present, the only results stored here are the populations' stellar massses (both screened and unscreened).

Within each of these directories (except `StellarPopulations`) are subdirectories, named as:
```
logZ%+.2f
```
which is the base 10 logarithm of the metallicity of the stellar population relative to the Sun ([M/H]).

Lastly, within each of these subdirectories are sub-subdirectories, using the names of each IMF, as given in the `IMFName` array.  For the default behavior, these subdirectories are:
```
Chabrier
Salpeter
BottomHeavy
```


## B. Spectra

### 1. File names
Sunscreen calculates the specific luminosity (L<sub>&nu;</sub>) of the model galaxy (and its putative Dyson spheres) and writes them to files in the subdirectory specified by type of result, IMF, and then metallicity (see IV-A).  The name of the spectrum files has the format:
```
%s_%s_Screened_Spectrum.logZ%+.2f.%s.%s.logt+%.2f.%s_%s.txt
```
The first two substrings ("`%s_%s`") are the same used in the name of the base directory storing the results.  For example, a spectrum of a galaxy where stars with luminosities below some threshold are screened is in the directory `Small_logL_Screened_Results`, and its spectrum name begins with `Small_logL_Screened_Spectrum`.  The quantity after the logZ is the base 10 logarithm of the metallicity relative to Solar ([M/H]).  The next substring is the IMF name.  After that comes the name of the SFH.  Next is the base 10 logarithm of the galaxy's age in years.  Next is generally the quantity that serves as a threshold, and finally the base 10 logarithm of the threshold value.  

Consider this file name:
```
Large_logLtoM0_Screened_Spectrum.logZ-1.00.Salpeter.BurstSF.logt+9.00.logLtoM0+1.50.txt
```
This is the spectrum of a galaxy where
* Stars with a log (L/M<sub>0</sub>) > 10<sup>1.5</sup> (L<sub>sun</sub>/M<sub>sun</sub>) are screened
* The galaxy has a Salpeter IMF.
* The galaxy's SFH is a single burst.
* The galaxy has an age of 1 Gyr.
* The galaxy has a metallicity of [M/H] = -1.0, or Z = 0.1 Z<sub>sun</sub>.

There are two exceptions to this rule.  When the screening is of individual stellar stages rather than by quantity, the name instead has the format:
```
%s_Stage_Unscreened_Spectrum.logZ%+.2f.%s.%s.logt+%.2f.Stage_%s.txt
```
The first substring is either "`Only`" or "`AllExcept`", the same as in the base output directory, and specifies whether stars of this stage are the only ones visible ("`Only_Stage_Unscreened`") or are the only ones that are invisible ("`AllExcept_Stage_Unscreened`").  Next comes metallicity in the form [M/H], IMF name, SFH name, and base 10 logarithm of the galaxy age in years.  The final substring is the short name of the stellar stage, as listed in Section III-B.

The other exception is that Sunscreen writes the unscreened/natural spectrum of the galaxy to a file in the UnscreenedResults directory, within the subdirectories corresponding to its IMF and then metallicity.  Their names have the form:
```
Unscreened_Spectrum.logZ%+.2f.%s.%s.logt+%.2f.txt
```
Again, the parameters in that name are [M/H], IMF name, SFH name, and base 10 logarithm of galaxy age.


### 2. Format
The format of the spectra are relatively simple text files.  The spectrum is output at the same wavelength used by the input BaSeL stellar model spectra (section III-A), which are listed in `LejeuneSpectra/lambda-LejeuneSpectra.txt`.  The first column of the file are these wavelengths.

The next column is the specific luminosity (L<sub>&nu;</sub>) of the galaxy, after taking into account any screening of its stars.  It has units of erg/s/Hz.  

After that are possibly more columns, which are the specific luminosity of the galaxy's Dyson spheres, which are assumed to capture and reradiate all the light of the screened stars at a particular temperature.  By default, Sunscreen gives the thermal emission from T = 1000 K and T = 300 K Dyson spheres.  These are also supposed to radiate with a blackbody spectrum.  

It is possible to change the number of megastructure temperatures for which this thermal emission is calculated and the values of those temperatures.  This can be done by:
1. Changing the `ND` constant to the number of desired output temperatures using the `#define ND` line.
2. Changing the temperatures listed in the `TDysonList` array's initial values.


## C. Synthetic photometry

### 1. Conventions

If the `CalculateMagnitudes` constant is defined to be TRUE, as it is by default, then Sunscreen calculates the magnitude of the galaxy from its spectrum.  

This simulated photometry is stored in a file in the same directory as the model's spectrum.  The file name is of the form:
```
%s_%s_Screened_Magnitudes.logZ%+.2f.%s.%s.logt+%.2f.%s.txt
```
The first two substrings determine the type of screening (e.g. "`Small_logL`" for galaxies where dim stars are screened).  Then comes metallicity in the form of [M/H], IMF name, SFH name, and base 10 logarithm of galaxy age.  The final substring is either "`Rest`" or "`Obs`".  "`Rest`" means that the photometry is calculated without any redshifting, by assuming that the galaxy is viewed in its rest frame.  "`Obs`" means the photometry is calculated *with* redshifting, as it would be observed from Earth if its age was the cosmic age, the same as the Universe's age.  So, for example, "`logt+9.50.Obs`" indicates that its redshift is the redshift where the Universe's age was 3.16 Gyr.

In converting from age to redshift, Sunscreen assumes a LambdaCDM cosmology with &Omega;<sub>m</sub> = 0.315, H<sub>0</sub> = 67.31 km/s, and the radiation energy density of the Universe to be entirely contained in the CMB with temperature 2.73 K.  These values can be altered by changing the `Omega_m_0`, `H_0`, and `TCMB_0` global variables, as well as the `t0_Cosmic` global variable giving the age of the Universe in seconds.  

Unlike the spectra, each value of the threshold quantity, like log L, is given only one line in the photometry file.  So the photometry for models of mass thresholds of 0.01 M<sub>sun</sub>, 1 M<sub>sun</sub>, and 1000 M<sub>sun</sub> are all given in the same file, if all other parameters (metallicity, IMF, SFH, age) are the same.  Essentially, the photometry files contain "tracks" of the galaxy in magnitude space as the threshold quantity is varied.

Unscreened magnitudes are also given in the `UnscreenedResults` directory, in subdirectories organized by metallicity and IMF.  The names of the files are similar but start with "`Unscreened_Magnitudes`".  There is only one set of magnitudes in these files.


### 2. Format

Each photometry file is a text file.  The first line is a header line listing each column's contents.  It starts with a "`#`" so it can be read as a comment by plotting programs.

After that are several columns.  The first column is the base 10 logarithm of whatever quantity determines screening.  For example, if all stars with luminosity above 100 L<sub>sun</sub> are screened, the value in the first column will be `+2.000000`.  When stellar stage is the cut, the short name of the stage being considered (e.g. "`RGB`") is listed instead.  In the unscreened spectra, this first column is just listed as "`All`".

The next column is the absolute bolometric magnitude of the visible/unscreened stars, as defined by IAU resolution.  

Finally are the absolute magnitudes of the visible/unscreened stellar population in each of the filter-bands.  These bands are listed in the `FilterName` array in Sunscreeen, and again at the header line of the file.  The magnitudes are separated by spaces.  Note the magnitudes do NOT include Dyson sphere thermal emission, even if the thermal emission is hot enough to appear in the optical.  

In some cases the model galaxy will appear to have zero unscreened luminosity, either in some bands or even bolometrically.  In these cases the magnitude is infinite and listed as "`inf`".


## D. Stellar masses
Within the StellarPopulations directory are files listing the stellar mass (both screened and unscreened) of model galaxies.  These are NOT sorted into further subdirectories of metallicity and IMF.  The file name is
```
PopulationMass.%s.%s.logZ%+.2f.logt+%.2f.txt
```
where the first substring is IMF name, the second is SFH name, the next quantity is metallicity [M/H], and the final quantity is base 10 logarithm of galaxy age in years.

These files contain a single number, which is the mass of the stellar population in Solar masses.


# V. References and Credit

Sunscreen includes data and models from the following sources.

## A. Stellar spectrum library

Sunscreen uses the BaSeL library of spectra from model stellar atmospheres.  These are described in the following papers:

* Lejeune, T., Cuisinier, F., & Buser, R., "Standard stellar library for evolutionary synthesis. I. Calibration of theoretical spectra", 1997, A&AS, 125, 229 [arXiv:astro-ph/9701019]
* Lejeune, T., Cuisinier, F., & Buser, R., "A standard stellar library for evolutionary synthesis. II. The M dwarf extension", 1998, A&AS, 130, 65 [arXiv:astro-ph/9710350] 

Data available online by following the "Online data" links from:
http://adsabs.harvard.edu/abs/1998yCat..41300065L


## B. Isochrones

Sunscreen uses PARSEC v1.2S + COLIBRI PR16 isochrones.  They are described in the following papers:

* Bressan, A., Marigo, P., Girardi, L., et al., "PARSEC: stellar tracks and isochrones with the PAdova and TRieste Stellar Evolution Code", 2012, MNRAS, 427, 127 [arXiv:1208.4498]
* Chen, Y., Girardi, L., Bressan, A., et al., "Improving PARSEC models for very low mass stars", 2014, MNRAS, 444, 2525 [arXiv:1409.0322]
* Tang, J., Bressan, A., Rosenfield, P., et al., "New PARSEC evolutionary tracks of massive stars at low metallicity: testing canonical stellar evolution in nearby star-forming dwarf galaxies",  2014, MNRAS, 445, 4287 [arXiv:1410.1745]
* Chen, Y., Bressan, A., Girardi, L., et al., "PARSEC evolutionary tracks of massive stars up to 350 M⊙ at metallicities 0.0001 ≤ Z ≤ 0.04", 2015, MNRAS, 452, 1068 [arXiv:1506.01681]
* Marigo, P., Girardi, L., Bressan, A., et al., "A New Generation of PARSEC-COLIBRI Stellar Isochrones Including the TP-AGB Phase", 2017, ApJ, 835, 77 [arXiv:1701.08510]

The "AGB boosting" effect mentioned in section III-B is described in:

* Girardi, L., Marigo, P., Bressan, A., & Rosenfield, P., "The Insidious Boosting of Thermally Pulsing Asymptotic Giant Branch Stars in Intermediate-age Magellanic Cloud Clusters", 2013, ApJ, 777, 142 [arXiv:1308.6088]

I downloaded isochrones from the CMD 3.0 input form.  The most recent version of CMD is available at:
http://stev.oapd.inaf.it/cgi-bin/cmd


## C. Star formation histories

Sunscreen uses SFHs compiled from several sources.

The "B13" SFHs are described in:
* Behroozi, P. S., Wechsler, R. H., & Conroy, C., "The Average Star Formation Histories of Galaxies in Dark Matter Halos from z = 0-8", 2013, ApJ, 770, 57 [arXiv:1207.6105]

They can be downloaded at:
https://www.peterbehroozi.com/data.html

The "M15" SFH comes from:

* McDermid, R. M., Alatalo, K., Blitz, L., et al., "The ATLAS3D Project - XXX. Star formation histories and stellar population scaling relations of early-type galaxies", 2015, MNRAS, 448, 3484 [arXiv:1501.03723]

  I used the yellow-orange line in the top panel of Figure 16.  As far as I can tell, the "SFR" is not the physical star-ormation rate, but the rate at which the present day stellar mass accumulated over time.  I have converted these to physical star-formation rates using:
* Leitner, S. N., & Kravtsov, A. V., "Fuel Efficient Galaxies: Sustaining Star Formation with Stellar Mass Loss", 2011, ApJ, 734, 48 [arXiv:1011.1252]

The "W14" SFH comes from:
* Weisz, D. R., Dolphin, A. E., Skillman, E. D., et al., "The Star Formation Histories of Local Group Dwarf Galaxies. I. Hubble Space Telescope/Wide Field Planetary Camera 2 Observations", 2014, ApJ, 789, 147 [arXiv:1404.7144]

  I used the thick blue line in Figure 13.  I took its derivative and averaged it over time.


## D. Filter responses
Further details on the filters can be found in Table 1 of Lacki (2018).

##### FUSE
* Morrissey, P., Schiminovich, D., Barlow, T. A., et al., "The On-Orbit Performance of the Galaxy Evolution Explorer", 2005, ApJL, 619, L7 [arXiv:astro-ph/0411310]

##### Johnson-Cousins (UBVRI)
* Bessell, M. S., "UBVRI passbands", 1990, PASP, 102, 1181 

##### SDSS (ugriz)
* Doi, M., Tanaka, M., Fukugita, M., et al., "Photometric Response Functions of the Sloan Digital Sky Survey Imager", 2010, AJ, 139, 1628 [arXiv:1002.3701]

##### DES (grizY)
* Burke, D. L., Rykoff, E. S., Allam, S., et al., "Forward Global Photometric Calibration of the Dark Energy Survey", 2018, AJ, 155, 41 [arXiv:1706.01542]

##### PanSTARRS-1 (grizyw)
* Tonry, J. L., Stubbs, C. W., Lykke, K. R., et al., "The Pan-STARRS1 Photometric System", 2012, ApJ, 750, 99 [arXiv:1203.0297]

##### LSST  (ugrizy)
* LSST Science Collaboration, Abell, P. A., Allison, J., et al., "LSST Science Book, Version 2.0" , 2009, arXiv:0912.0201
* https://github.com/lsst/throughputs/tree/master/baseline

##### Gaia (G, G<sub>BP</sub>, G<sub>RP</sub>)
* Carrasco, J. M., Evans, D. W., Montegriffo, P., et al., "Gaia Data Release 1. Principles of the photometric calibration of the G band", 2016, A&A, 595, A7 [arXiv:1611.02036]
* Evans, D. W., Riello, M., De Angeli, F., et al., "Gaia Data Release 2. Photometric content and validation", 2018, A&A, 616, A4 [arXiv:1804.09368]

##### 2MASS (J, H, K<sub>s</sub>)
* Cohen, M., Wheaton, W. A., & Megeath, S. T., "Spectral Irradiance Calibration in the Infrared. XIV. The Absolute Calibration of 2MASS", 2003, AJ, 126, 1090 [arXiv:astro-ph/0304350]
* Skrutskie, M. F., Cutri, R. M., Stiening, R., et al., "The Two Micron All Sky Survey (2MASS)", 2006, AJ, 131, 1163

##### UKIRT (ZYJHK)
* Hewett, P. C., Warren, S. J., Leggett, S. K., & Hodgkin, S. T., "The UKIRT Infrared Deep Sky Survey ZY JHK photometric system: passbands and synthetic colours", 2006, MNRAS, 367, 454 [arXiv:astro-ph/0601592]

##### Spitzer IRAC (IRAC-1, IRAC-2, IRAC-3, IRAC-4)
* Fazio, G. G., Hora, J. L., Allen, L. E., et al., "The Infrared Array Camera (IRAC) for the Spitzer Space Telescope", 2004, ApJS, 154, 10 [arXiv:astro-ph/0405616]

##### WISE (W1, W2, W3, W4)
Wright, E. L., Eisenhardt, P. R. M., Mainzer, A. K., et al., "The Wide-field Infrared Survey Explorer (WISE): Mission Description and Initial On-orbit Performance", 2010, AJ, 140, 1868 [arXiv:1008.0031]

##### JWST
Although they are not calculated by default when Sunscreen is run, Sunscreen does include properly formatted PhotonZeroPoint files for JWST filters.  These come from:
* https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters
  I used the mean for modules A & B.  See Section III-D for instructions on how to add them to the list of filters Sunscreen calculates.


## E. Citing Sunscreen

Cite this paper:
* Lacki, B. C., "Sunscreen: Photometric Signatures of Galaxies Partially Cloaked in Dyson Sphere", 2018, arXiv:1807.00077 


## F. Support and acknowledgements

I was supported by the Breakthrough Listen program.  Funding for Breakthrough Listen research is sponsored by the Breakthrough Prize Foundation:
https://breakthroughprize.org/

See the following reference for Breakthrough Listen:
* Worden, S. P., Drew, J., Siemion, A., et al., "Breakthrough Listen - A new search for life in the universe", 2017, Acta Astronautica, 139, 98 
