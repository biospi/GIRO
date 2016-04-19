# GIRO
Groupwise Image Registration and nOrmalisation: 

In shotgun LCMS proteomical experiment, it is well known that the retention time of different 
chromatography runs may vary. Moreover the quantitation of the peptide may also vary. GIRO provides
two functions: retention time alignment and abundance normalisation.  

GIRO is a novel signal processing based method for LCMS retention time adjustment. It takes every
LCMS samples as an image, and by applying groupwise image registration, it aligns the retention time
of all LCMS images into the same coordinate system so that the RT drift can be recovered in Silico. 
Moreover, the abundance normalisation is performed after RT alignment, which is also essential towards
a more accurate quantitation outcome.

In this version 2.0 the original class @GIRO is factorised into three lines of classes to gain reuseability. 
@Data is the base class providing interfaces for generic data application.
@DataSeaMass <- @Data is the subclass specially used to get SeaMass coefficients from smo files.

@Register is the base class providing interfaces for generic registration algorithms to deform the data.
@MultiresRegister <- @Register is the subclass using multiresolution scheme to do registration.

@Normalise is the base class providing interfaces for generic normalisation algorithms to normalise the data.

and finally 
@GIRO <- @MultiresRegister & @Normalise is the class derived from both registration and normalisation to 
analyze the LCMS data.
